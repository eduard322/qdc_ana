import ROOT,os,subprocess
import atexit
import time
import ctypes
from array import array
import numpy as np
import ana_kit as ana
from multiprocessing import Pool
# for fixing a root bug,  will be solved in the forthcoming 6.26 release.
ROOT.gInterpreter.Declare("""
#include "MuFilterHit.h"
#include "AbsMeasurement.h"
#include "TrackPoint.h"

void fixRoot(MuFilterHit& aHit,std::vector<int>& key,std::vector<float>& value, bool mask) {
   std::map<int,float> m = aHit.GetAllSignals(false);
   std::map<int, float>::iterator it = m.begin();
   while (it != m.end())
    {
        key.push_back(it->first);
        value.push_back(it->second);
        it++;
    }
}
void fixRootT(MuFilterHit& aHit,std::vector<int>& key,std::vector<float>& value, bool mask) {
   std::map<int,float> m = aHit.GetAllTimes(false);
   std::map<int, float>::iterator it = m.begin();
   while (it != m.end())
    {
        key.push_back(it->first);
        value.push_back(it->second);
        it++;
    }
}
void fixRoot(MuFilterHit& aHit, std::vector<TString>& key,std::vector<float>& value, bool mask) {
   std::map<TString, float> m = aHit.SumOfSignals();
   std::map<TString, float>::iterator it = m.begin();
   while (it != m.end())
    {
        key.push_back(it->first);
        value.push_back(it->second);
        it++;
    }
}

void fixRoot(std::vector<genfit::TrackPoint*>& points, std::vector<int>& d,std::vector<int>& k, bool mask) {
      for(std::size_t i = 0; i < points.size(); ++i) {
        genfit::AbsMeasurement*  m = points[i]->getRawMeasurement();
        d.push_back( m->getDetId() );
        k.push_back( int(m->getHitId()/1000) );
    }
}
""")
def pyExit():
       print("Make suicide until solution found for freezing")
       os.system('kill '+str(os.getpid()))
atexit.register(pyExit)

Tkey  = ROOT.std.vector('TString')()
Ikey   = ROOT.std.vector('int')()
Value = ROOT.std.vector('float')()

def map2Dict(aHit,T='GetAllSignals',mask=True):
     if T=="SumOfSignals":
        key = Tkey
     elif T=="GetAllSignals" or T=="GetAllTimes":
        key = Ikey
     else: 
           print('use case not known',T)
           1/0
     key.clear()
     Value.clear()
     if T=="GetAllTimes": ROOT.fixRootT(aHit,key,Value,mask)
     else:                         ROOT.fixRoot(aHit,key,Value,mask)
     theDict = {}
     for k in range(key.size()):
         if T=="SumOfSignals": theDict[key[k].Data()] = Value[k]
         else: theDict[key[k]] = Value[k]
     return theDict

import rootUtils as ut
import shipunit as u
h={}
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-r", "--runNumber", dest="runNumber", help="run number", type=int,required=True)
parser.add_argument("-p", "--path", dest="path", help="run number",required=False,default="")
parser.add_argument("-f", "--inputFile", dest="fname", help="file name for MC", type=str,default=None,required=False)
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geofile", required=True)
parser.add_argument("-b", "--heartBeat", dest="heartBeat", help="heart beat", default=10000,type=int)
parser.add_argument("-c", "--command", dest="command", help="command", default="")
parser.add_argument("-n", "--nEvents", dest="nEvents", help="number of events", default=-1,type=int)
parser.add_argument("-t", "--trackType", dest="trackType", help="DS or Scifi", default="DS")
parser.add_argument("--remakeScifiClusters", dest="remakeScifiClusters", help="remake ScifiClusters", default=False)
options = parser.parse_args()

runNr   = str(options.runNumber).zfill(6)
path     = options.path
partitions = 0
if options.runNumber > 0: 
 if path.find('eos')>0:
    path     = os.environ['EOSSHIP']+options.path
    dirlist  = str( subprocess.check_output("xrdfs "+os.environ['EOSSHIP']+" ls "+options.path,shell=True) )
# single file, before Feb'22
    data = "sndsw_raw_"+runNr+".root"
    if  dirlist.find(data)<0:
# check for partitions
       dirlist  = str( subprocess.check_output("xrdfs "+os.environ['EOSSHIP']+" ls "+options.path+"run_"+runNr,shell=True) )
       while 1>0:
        data = "raw-"+ str(partitions).zfill(4)
        if dirlist.find(data)>0:
            partitions+=1
        else: break
 else:
# check for partitions
       data = "sndsw_raw_"+runNr+".root"
       dirlist = os.listdir(options.path)
       if  not data in dirlist:
          dirlist  = os.listdir(options.path+"run_"+runNr)
          for x in dirlist:
            data = "raw-"+ str(partitions).zfill(4)
            if x.find(data)>0:
               partitions+=1

import SndlhcGeo
if (options.geoFile).find('../')<0: geo = SndlhcGeo.GeoInterface(path+options.geoFile)
else:                                         geo = SndlhcGeo.GeoInterface(options.geoFile[3:])
MuFilter = geo.modules['MuFilter']
Scifi       = geo.modules['Scifi']
nav = ROOT.gGeoManager.GetCurrentNavigator()

A,B,locA,locB,globA,globB    = ROOT.TVector3(),ROOT.TVector3(),ROOT.TVector3(),ROOT.TVector3(),ROOT.TVector3(),ROOT.TVector3()
latex = ROOT.TLatex()


# MuFilter mapping of planes and bars 
systemAndPlanes  = {1:2,2:5,3:7}
systemAndBars     = {1:7,2:10,3:60}
def systemAndOrientation(s,plane):
   if s==1 or s==2: return "horizontal"
   if plane%2==1 or plane == 6: return "vertical"
   return "horizontal"

systemAndChannels     = {1:[8,0],2:[6,2],3:[1,0]}
sdict                     = {1:'Veto',2:'US',3:'DS'}

freq      = 160.316E6
TDC2ns = 1E9/freq

# some helper functions

def getAverageZpositions():
   zPos={'MuFilter':{},'Scifi':{}}
   for s in systemAndPlanes:
       for plane in range(systemAndPlanes[s]):
          bar = 4
          p = plane
          if s==3 and (plane%2==0 or plane==7): 
             bar = 90
             p = plane//2
          elif s==3 and plane%2==1:
             bar = 30
             p = plane//2
          MuFilter.GetPosition(s*10000+p*1000+bar,A,B)
          zPos['MuFilter'][s*10+plane] = (A.Z()+B.Z())/2.
   for s in range(1,6):
      mat   = 2
      sipm = 1
      channel = 64
      for o in range(2):
          Scifi.GetPosition(channel+1000*sipm+10000*mat+100000*o+1000000*s,A,B)
          zPos['Scifi'][s*10+o] = (A.Z()+B.Z())/2.
   return zPos

def smallSiPMchannel(i):
    if i==2 or i==5 or i==10 or i==13: return True
    else: return False

def fit_langau(hist,o,bmin,bmax,opt=''):
   if opt == 2:
     params = {0:'Width(scale)',1:'mostProbable',2:'norm',3:'sigma',4:'N2'}
     F = ROOT.TF1('langau',langaufun,0,200,len(params))
   else:
     params = {0:'Width(scale)',1:'mostProbable',2:'norm',3:'sigma'}
     F = ROOT.TF1('langau',twoLangaufun,0,200,len(params))
   for p in params: F.SetParName(p,params[p])
   rc = hist.Fit('landau','S'+o,'',bmin,bmax)
   res = rc.Get()
   if not res: return res
   F.SetParameter(2,res.Parameter(0))
   F.SetParameter(1,res.Parameter(1))
   F.SetParameter(0,res.Parameter(2))
   F.SetParameter(3,res.Parameter(2))
   F.SetParameter(4,0)
   F.SetParLimits(0,0,100)
   F.SetParLimits(1,0,100)
   F.SetParLimits(3,0,10)
   rc = hist.Fit(F,'S'+o,'',bmin,bmax)
   res = rc.Get()
   return res

def twoLangaufun(x,par):
   N1 = langaufun(x,par)
   par2 = [par[0],par[1]*2,par[4],par[3]]
   N2 = langaufun(x,par2)
   return N1+abs(N2)

def  langaufun(x,par):
   #Fit parameters:
   #par[0]=Width (scale) parameter of Landau density
   #par[1]=Most Probable (MP, location) parameter of Landau density
   #par[2]=Total area (integral -inf to inf, normalization constant)
   #par[3]=Width (sigma) of convoluted Gaussian function
   #
   #In the Landau distribution (represented by the CERNLIB approximation),
   #the maximum is located at x=-0.22278298 with the location parameter=0.
   #This shift is corrected within this function, so that the actual
   #maximum is identical to the MP parameter.
#
      # Numeric constants
      invsq2pi = 0.3989422804014   # (2 pi)^(-1/2)
      mpshift  = -0.22278298       # Landau maximum location
#
      # Control constants
      np = 100.0      # number of convolution steps
      sc =   5.0      # convolution extends to +-sc Gaussian sigmas
#
      # Variables
      summe = 0.0
#
      # MP shift correction
      mpc = par[1] - mpshift * par[0]
#
      # Range of convolution integral
      xlow = max(0,x[0] - sc * par[3])
      xupp = x[0] + sc * par[3]
#
      step = (xupp-xlow) / np
#
      # Convolution integral of Landau and Gaussian by sum
      i=1.0
      if par[0]==0 or par[3]==0: return 9999
      while i<=np/2:
         i+=1
         xx = xlow + (i-.5) * step
         fland = ROOT.TMath.Landau(xx,mpc,par[0]) / par[0]
         summe += fland * ROOT.TMath.Gaus(x[0],xx,par[3])
#
         xx = xupp - (i-.5) * step
         fland = ROOT.TMath.Landau(xx,mpc,par[0]) / par[0]
         summe += fland * ROOT.TMath.Gaus(x[0],xx,par[3])
#
      return (par[2] * step * summe * invsq2pi / par[3])

def myPrint(tc,name,withRootFile=True):
     tc.Update()
     tc.Print(name+'-run'+str(options.runNumber)+'.png')
     tc.Print(name+'-run'+str(options.runNumber)+'.pdf')
     if withRootFile: tc.Print(name+'-run'+str(options.runNumber)+'.root')

def makeAnimation(histname,j0=1,j1=2,animated=True, findMinMax=True, lim = 50):
    ut.bookCanvas(h,'ani','',900,800,1,1)
    tc = h['ani'].cd()
    jmin,jmax = j0,j1
    if findMinMax:
       jmin,jmax = 0,0
       for j in range(j0,j1):
            tmp = histname.replace('$$$',str(j))
            if tmp in h:
                  if h[tmp].GetEntries()>lim:
                       jmin = j
                       print(j,tmp,h[tmp].GetEntries())
                       break
       for j in range(j1,j0,-1):
            tmp = histname.replace('$$$',str(j))
            if tmp in h:
                  if h[tmp].GetEntries()>lim:
                       jmax = j
                       break
    for j in range(jmin,jmax):
            tmp = histname.replace('$$$',str(j))
            h[tmp].Draw()
            tc.Update()
            stats = h[tmp].FindObject('stats')
            stats.SetOptFit(1111111)
            h[tmp].Draw()
            if animated: 
               h['ani'].Print('picAni'+str(j)+'.png')
            else:
               rc = input("hit return for next event or q for quit: ")
               if rc=='q': break
    if animated and jmax > jmin: 
           os.system("convert -delay 60 -loop 0 picAni*.png "+histname+".gif")
           os.system('rm picAni*.png')



# initialize 

zPos = getAverageZpositions()

if options.runNumber>0:
              eventChain = ROOT.TChain('rawConv')
              if partitions==0:
                   eventChain.Add(path+'sndsw_raw_'+str(options.runNumber).zfill(6)+'.root')
              else:
                   for p in range(partitions):
                       eventChain.Add(path+'run_'+runNr+'/sndsw_raw-'+str(p).zfill(4)+'.root')

else:
# for MC data and other files
              f=ROOT.TFile.Open(options.fname)
              if f.Get('rawConv'):   eventChain = f.rawConv
              else:                        eventChain = f.cbmsim
if options.remakeScifiClusters: eventChain.SetBranchStatus("Cluster_Scifi*",0)
rc = eventChain.GetEvent(0)
run      = ROOT.FairRunAna()
ioman = ROOT.FairRootManager.Instance()
ioman.SetTreeName(eventChain.GetName())
outFile = ROOT.TMemFile('dummy','CREATE')
source = ROOT.FairFileSource(eventChain.GetCurrentFile())

if partitions>0:
    for p in range(1,partitions):
                       source.AddFile(path+'run_'+runNr+'/sndsw_raw-'+str(p).zfill(4)+'.root')
run.SetSource(source)
sink = ROOT.FairRootFileSink(outFile)
run.SetSink(sink)

houghTransform = False # under construction, not yet tested
if houghTransform:
   muon_reco_task = SndlhcMuonReco.MuonReco()
   muon_reco_task.SetName("houghTransform")
   run.AddTask(muon_reco_task)
else:
   import SndlhcTracking
   trackTask = SndlhcTracking.Tracking()
   trackTask.SetName('simpleTracking')
   run.AddTask(trackTask)

#avoiding some error messages
xrdb = ROOT.FairRuntimeDb.instance()
xrdb.getContainer("FairBaseParSet").setStatic()
xrdb.getContainer("FairGeoParSet").setStatic()

run.Init()
if partitions>0:  eventTree = ioman.GetInChain()
else:                 eventTree = ioman.GetInTree()

OT = ioman.GetSink().GetOutTree()
Reco_MuonTracks = trackTask.fittedTracks
Cluster_Scifi         = trackTask.clusScifi
# wait for user action 

def help():
    print(" following methods exist")
    print("     make and plot  hitmaps, signal distributions for MuFIlter and Scifi:")
    print("              Scifi_hitMaps(Nev) and Mufi_hitMaps(Nev)     if no number of events is given, loop over all events in file.")
    print(" ")
    print("     plot time between events and time since first event")
    print("              eventTime(Nev=-1)")
    print(" ")
    print("     MuFilter residuals, efficiency estimates with DS or Scifi tracks")
    print("              Mufi_Efficiency(Nev=-1,optionTrack='DS' or 'Scifi")
    print(" ")
    print("     analyze and plot historgams made withMufi_Efficiency ")
    print("              analyze_EfficiencyAndResiduals(readHists=False), with option True, histograms are read from file")
    print(" ")
    print("     Scifi unbiased residuals for an optional set of alignment constants")
    print("              Scifi_residuals(Nev=-1,NbinsRes=100,xmin=-500.,alignPar=False), xmin = low edge of residual histos in microns")
    print(" ")
    print("     Minimization of Scifi alignment constants")
    print("              minimizeAlignScifi(first=True,level=1,minuit=False)")
    print(" ")
    print("     first attempt to look at hadronic showers ")
    print("              USshower(Nev=-1)")
    print(" ")
    print("     make histograms with QDC distributions and fit all distributions with Langau ")
    print("              mips()")
    print("     plot the above distributions directly or reading from file ")
    print("              plotMips(readhisto=True)")
    print(" ")
    print("     beam illumination ")
    print("             scifi_beamSpot() and beamSpot() for MuFilter")
    print(" ")
    print("     rough estimate of Scifi resolution by comparing slopes  build with two stations, x and y projection")
    print("             Scifi_slopes(Nev=-1)")


def Mufi_Efficiency(Nev=options.nEvents,optionTrack=options.trackType,withReco='True',NbinsRes=100,X=10.):
 
 projs = {1:'Y',0:'X'}
 for s in range(1,6):
     for p in projs:
         ut.bookHist(h,'dtScifivsX_'+str(s)+projs[p],'dt vs x track '+projs[p]+";X [cm]; dt [ns]",100,-10.,40.,260,-8.,5.)
         ut.bookHist(h,'clN_'+str(s)+projs[p],'cluster size '+projs[p],10,-0.5,9.5)
     ut.bookHist(h,'dtScifivsdL_'+str(s),'dt vs dL '+str(s)+";X [cm]; dt [ns]",100,-40.,0.,200,-5.,5.)

 for s in systemAndPlanes:
    for l in range(systemAndPlanes[s]):
      ut.bookHist(h,'dtLRvsX_'+sdict[s]+str(s*10+l),'dt vs x track '+str(s*10+l)+";X [cm]; dt [ns]",100,-70.,30.,260,-8.,5.)
      ut.bookHist(h,'atLRvsX_'+sdict[s]+str(s*10+l),'mean time - T0track vs x '+str(s*10+l)+";X [cm]; dt [ns]",20,-70.,30.,250,-10.,15.0)
      ut.bookHist(h,'VetoatLRvsX_'+sdict[s]+str(s*10+l),'mean time - T0Veto vs x '+str(s*10+l)+";X [cm]; dt [ns]",20,-70.,30.,250,-10.,15.0)

      scale = 1.
      if s==3: scale = 0.4
      for side in ['','L','R','S']:
        for proj in ['X','Y']:
          xmin = -X*NbinsRes/100. * scale
          xmax = -xmin
          ut.bookHist(h,'res'+proj+'_'+sdict[s]+side+str(s*10+l),'residual  '+proj+str(s*10+l),NbinsRes,xmin,xmax,100,-100.,100.)
          ut.bookHist(h,'gres'+proj+'_'+sdict[s]+side+str(s*10+l),'residual  '+proj+str(s*10+l),NbinsRes,xmin,xmax,100,-100.,100.)
          if side=='S': continue
          if side=='':
             if s==1: ut.bookHist(h,'resBar'+proj+'_'+sdict[s]+str(s*10+l),'residual '+proj+str(s*10+l),NbinsRes,xmin,xmax,7,-0.5,6.5)
        if side=='':
             ut.bookHist(h,'track_'+sdict[s]+str(s*10+l),'track x/y '+str(s*10+l)+';X [cm];Y [cm]',100,-90.,10.,100,-20.,80.)
             ut.bookHist(h,'locBarPos_'+sdict[s]+str(s*10+l),'bar sizes;X [cm];Y [cm]',100,-100,100,100,-100,100)
             ut.bookHist(h,'locEx_'+sdict[s]+str(s*10+l),'loc track pos;X [cm];Y [cm]',100,-100,100,100,-100,100)
        for bar in range(systemAndBars[s]):
             key = sdict[s]+str(s*10+l)+'_'+str(bar)
             if side=="":
                ut.bookHist(h,'dtLRvsX_'+key,'dt vs x track '+str(s*10+l)+";X [cm]; dt [ns]",100,-70.,30.,260,-8.,5.)
                ut.bookHist(h,'dtF1LRvsX_'+key,'dt vs x track '+str(s*10+l)+";X [cm]; dt [ns]",100,-70.,30.,260,-8.,5.)
                ut.bookHist(h,'dtfastLRvsX_'+key,'dt vs x track '+str(s*10+l)+";X [cm]; dt [ns]",100,-70.,30.,260,-8.,5.)
                ut.bookHist(h,'atLRvsX_'+key,'dt vs x track '+str(s*10+l)+";X [cm]; dt [ns]",100,-70.,30.,260,-8.,5.)
             else:
                ut.bookHist(h,'nSiPMs'+side+'_'+key,'#sipms',16,-0.5,15.5,20,0.,100.)
                ut.bookHist(h,'tvsX'+side+'_'+key,"t-t0 vs x track;X [cm]; dt [ns]",100,-70.,30.,200,-12.,12.)
                ut.bookHist(h,'tFastvsX'+side+'_'+key,"t-t0 vs x track;X [cm]; dt [ns]",100,-70.,30.,200,-12.,12.)
                for i in range(systemAndChannels[s][1]+systemAndChannels[s][0]):
                        if s==2 and smallSiPMchannel(i):
                            ut.bookHist(h,'signalS'+side+'_'+key+'-c'+str(i),'signal',200,0.,100.,20,0.,100.)
                        else:
                            ut.bookHist(h,'signal'+side+'_'+key+'-c'+str(i),'signal',200,0.,100.,20,0.,100.)
                        if s==3: continue
                        ut.bookHist(h,'sigmaTDC'+side+'_'+key+'-c'+str(i),'rms TDC ;dt [ns]',200,-10.0,10.0)
                        ut.bookHist(h,'TDCcalib'+side+'_'+key+'-c'+str(i),'rms TDC ;dt [ns]',200,-10.0,10.0)
                        ut.bookHist(h,'sigmaQDC'+side+'_'+key+'-c'+str(i),'rms QDC ; QDC ',200,-50.0,50.)
                        ut.bookHist(h,'tvsX'+side+'_'+key+'-c'+str(i),"t-t0 vs x track;X [cm]; dt [ns]",100,-70.,30.,200,-12.,12.)

                ut.bookHist(h,'signalT'+side+'_'+key,'signal',400,0.,400.,20,0.,100.)
                ut.bookHist(h,'signalTS'+side+'_'+key,'signal',400,0.,400.,20,0.,100.)
                ut.bookHist(h,'signal'+side+'_'+key,'signal',200,0.,100.,20,0.,100.)

 ut.bookHist(h,'resVETOY_1','channel vs residual  1',NbinsRes,xmin,xmax,112,-0.5,111.5)
 ut.bookHist(h,'resVETOY_2','channel vs residual  2',NbinsRes,xmin,xmax,112,-0.5,111.5)

 ut.bookHist(h,'trackslxy','track direction',200,-0.1,0.1,200,-0.1,0.1)
 ut.bookHist(h,'trackslxy_badChi2','track direction',200,-0.1,0.1,200,-0.1,0.1)
 ut.bookHist(h,'tracksChi2Ndof','chi2/ndof',100,0.0,100.,10,-0.5,9.5)
 ut.bookHist(h,'NdofvsNMeas','ndof Nmeas',20,-0.5,19.5,20,-0.5,19.5)

 v = 15.* u.cm/u.ns # signal propagation in fibre, assumption
 if Nev < 0 : Nev = eventTree.GetEntries()
 N=0
 for event in eventTree:
    rc = eventTree.GetEvent(N)
    N+=1
    if N>Nev: break
    if withReco:
       for aTrack in Reco_MuonTracks: aTrack.Delete()
       Reco_MuonTracks.Clear()
       if optionTrack=='DS': rc = trackTask.ExecuteTask("DS")
       else                         : rc = trackTask.ExecuteTask("Scifi")
    if not Reco_MuonTracks.GetEntries()==1: continue
    theTrack = Reco_MuonTracks[0]
    if not theTrack.getFitStatus().isFitConverged() and optionTrack!='DS':   # for H8 where only two planes / proj were avaiable
           continue
# only take horizontal tracks
    state = theTrack.getFittedState(0)
    pos   = state.getPos()
    mom = state.getMom()
    slopeX= mom.X()/mom.Z()
    slopeY= mom.Y()/mom.Z()
    if abs(slopeX)>0.25: continue   # 4cm distance, 250mrad = 1cm
    if abs(slopeY)>0.1: continue

# now extrapolate to US and check for hits.
    fitStatus = theTrack.getFitStatus()
    chi2Ndof = fitStatus.getChi2()/(fitStatus.getNdf()+1E-10)
    rc = h['tracksChi2Ndof'].Fill(chi2Ndof,fitStatus.getNdf())
    rc = h['NdofvsNMeas'].Fill(fitStatus.getNdf(),theTrack.getNumPointsWithMeasurement())
# map clusters to hit keys
    DetID2Key={}
    if hasattr(event,"Cluster_Scifi"):
               clusters = event.Cluster_Scifi
    else:
               clusters = Cluster_Scifi

    for aCluster in clusters:
        for nHit in range(event.Digi_ScifiHits.GetEntries()):
            if event.Digi_ScifiHits[nHit].GetDetectorID()==aCluster.GetFirst():
               DetID2Key[aCluster.GetFirst()] = nHit
    for aCluster in clusters:
         detID = aCluster.GetFirst()
         s = int(detID/1000000)
         p= int(detID/100000)%10
         rc = h['clN_'+str(s)+projs[p]].Fill(aCluster.GetN())

    if chi2Ndof> 9 and optionTrack!='DS': 
       rc = h['trackslxy_badChi2'].Fill(mom.x()/mom.Mag(),mom.y()/mom.Mag())
       continue
    rc = h['trackslxy'].Fill(mom.x()/mom.Mag(),mom.y()/mom.Mag())
# get T0 from Track
    if optionTrack=='DS':
    # define T0 using mean TDC L/R of the horizontal planes
         T0track = 0
         Z0track = -1
         T0s = []
         for nM in range(theTrack.getNumPointsWithMeasurement()):
            M = theTrack.getPointWithMeasurement(nM)
            W = M.getRawMeasurement()
            detID = W.getDetId()
            hkey  = W.getHitId()
            aHit = event.Digi_MuFilterHits[ hkey ]
            if aHit.isVertical(): continue
            allTDCs = map2Dict(aHit,'GetAllTimes')
            if (0 in allTDCs) and (1 in allTDCs) : T0s.append( (allTDCs[0]+allTDCs[1])*TDC2ns/2. )
         if len(T0s)==2:
            lam = (zPos['MuFilter'][32]-pos.z())/mom.z()
            xEx,yEx = lam*mom.x(),pos.y()+lam*mom.y()
            # need track length
            L = ROOT.TMath.Sqrt( (lam*mom.x())**2 + (lam*mom.y())**2 + (zPos['MuFilter'][32]-pos.z())**2)
            ToF = L / u.speedOfLight
            T0track = (T0s[0]+T0s[1]-ToF)/2.
            Z0track = pos[2]
            deltaZ02 = (zPos['MuFilter'][32]-zPos['MuFilter'][30])
            # print('debug 1',L,ToF,T0track,Z0track,deltaZ02)
    else:
         M = theTrack.getPointWithMeasurement(0)
         W = M.getRawMeasurement()
         detID = W.getDetId()
         aHit = event.Digi_ScifiHits[ DetID2Key[detID] ]
         geo.modules['Scifi'].GetSiPMPosition(detID,A,B)
         X = B-pos
         L0 = X.Mag()/v
         # need to correct for signal propagation along fibre
         clkey  = W.getHitId()
         aCl = event.Cluster_Scifi[clkey]
         T0track = aCl.GetTime() - L0
         TZero    = aCl.GetTime()
         Z0track = pos[2]
         times = {}
         for nM in range(theTrack.getNumPointsWithMeasurement()):
            state   = theTrack.getFittedState(nM)
            posM   = state.getPos()
            M = theTrack.getPointWithMeasurement(nM)
            W = M.getRawMeasurement()
            detID = W.getDetId()
            clkey = W.getHitId()
            aCl = event.Cluster_Scifi[clkey]
            aHit = event.Digi_ScifiHits[ DetID2Key[detID] ]
            geo.modules['Scifi'].GetSiPMPosition(detID,A,B)
            if aHit.isVertical(): X = B-posM
            else: X = A-posM
            L = X.Mag()/v
         # need to correct for signal propagation along fibre
            dT = aCl.GetTime() - L - T0track - (posM[2] -Z0track)/u.speedOfLight
            ss = str(aHit.GetStation())
            prj = 'X'
            l = posM[0]
            if aHit.isVertical(): 
                 prj='Y'
                 l = posM[1]
            rc = h['dtScifivsX_'+ss+prj].Fill(X.Mag(),dT)
            times[ss+prj]=[aCl.GetTime(),L*v,detID,l]
         for s in range(1,6):
            if str(s)+'X' in times and  str(s)+'Y' in times:
               deltaT = times[str(s)+'X'][0] - times[str(s)+'Y'][0]
               deltaL = times[str(s)+'X'][1] - times[str(s)+'Y'][1]
               rc = h['dtScifivsdL_'+str(s)].Fill(deltaL,deltaT)
         deltaZ02 = 40. # to be fixed
         ToF = 1.
            #print(detID,aHit.GetDetectorID(),aHit.GetTime()*TDC2ns-TZero,dT,L,aHit.GetTime()*TDC2ns - L,T0)

    muHits = {}
    for s in systemAndPlanes:
       for p in range(systemAndPlanes[s]):
          if s == 3:
             muHits[s*10+2*p]=[]
             muHits[s*10+2*p+1]=[]
          else:
             muHits[s*10+p]=[]
       for p in range(systemAndPlanes[s]): muHits[s*10+p]=[]
    for aHit in event.Digi_MuFilterHits:
         if not aHit.isValid(): continue
         s = aHit.GetDetectorID()//10000
         p = (aHit.GetDetectorID()//1000)%10
         bar = (aHit.GetDetectorID()%1000)%60
         plane = s*10+p
         if s==3:
           if aHit.isVertical(): plane = s*10+2*p+1
           else:                     plane = s*10+2*p
         muHits[plane].append(aHit)

# get T0 from VETO
    s = 1
    Z0Veto = zPos['MuFilter'][1*10+0]
    dZ = zPos['MuFilter'][1*10+1] - zPos['MuFilter'][1*10+0]
    avT = {}
    for p in range(systemAndPlanes[s]): 
         plane = s*10+p
         if len(muHits[plane])!=1: continue
         aHit = muHits[plane][0]
# check if hit within track extrapolation
         zEx = zPos['MuFilter'][s*10+plane]
         lam = (zEx-pos.z())/mom.z()
         xEx,yEx = pos.x()+lam*mom.x(),pos.y()+lam*mom.y()
         detID = aHit.GetDetectorID()
         MuFilter.GetPosition(detID,A,B)
         D = (A[1]+B[1])/2. - yEx
         if abs(D)>5: continue
         avT[plane] = aHit.GetImpactT()
    T0Veto = -999
    if len(avT)==2:
         T0Veto = (avT[10]+(avT[11]-dZ/u.speedOfLight))/2.

    vetoHits = {0:[],1:[]}
    for s in sdict:
     name = str(s)
     for plane in range(systemAndPlanes[s]):
         zEx = zPos['MuFilter'][s*10+plane]
         lam = (zEx-pos.z())/mom.z()
         xEx,yEx = pos.x()+lam*mom.x(),pos.y()+lam*mom.y()
         # tag with station close by
         if plane ==0: tag = 1
         else: tag = plane -1
         tagged = False
         for aHit in muHits[s*10+tag]:
              detID = aHit.GetDetectorID()
              MuFilter.GetPosition(detID,A,B)
              if aHit.isVertical() : D = (A[0]+B[0])/2. - xEx
              else:                      D = (A[1]+B[1])/2. - yEx
              if abs(D)<5: tagged = True
         #if not tagged: continue
         rc = h['track_'+sdict[s]+str(s*10+plane)].Fill(xEx,yEx)
         for aHit in muHits[s*10+plane]:
              detID = aHit.GetDetectorID()
              bar = (detID%1000)%60
              nSiPMs = aHit.GetnSiPMs()
              nSides  = aHit.GetnSides()
              MuFilter.GetPosition(detID,globA,globB)
              MuFilter.GetLocalPosition(detID,locA,locB)
              globEx = array('d',[xEx,yEx,zEx])
              locEx   = array('d',[0,0,0])
              nav.MasterToLocal(globEx,locEx)
              locPos   = 0.5*  (locA+locB)
              globPos = 0.5 * (globA+globB)
              dy = locPos[1] - locEx[1]
              dx = locPos[0] - locEx[0]
              gdy = globPos[1] - globEx[1]
              gdx = globPos[0] - globEx[0]
              rc = h['locBarPos_'+sdict[s]+str(s*10+plane)].Fill( locPos[0],locPos[1])
              rc = h['locEx_'+sdict[s]+str(s*10+plane)].Fill( locEx[0],locEx[1])
              rc = h['resY_'+sdict[s]+str(s*10+plane)].Fill(dy,locEx[0])
              rc = h['resX_'+sdict[s]+str(s*10+plane)].Fill(dx,locEx[1])
              rc = h['gresY_'+sdict[s]+str(s*10+plane)].Fill(gdy,globEx[0])
              rc = h['gresX_'+sdict[s]+str(s*10+plane)].Fill(gdx,globEx[1])
              S = map2Dict(aHit,'GetAllSignals')
              # check for signal in left / right or small sipm
              left,right,smallL,smallR,Sleft,Sright,SleftS,SrightS = 0,0,0,0,0,0,0,0
              if  s==1:        
                    vetoHits[plane].append( [gdy,bar] )
                    rc = h['resBarY_'+sdict[s]+str(s*10+plane)].Fill(gdy,bar)
              for x in S:
                  if  s==1:
                      nc = x + 2*nSiPMs*bar
                      h['resVETOY_'+str(plane+1)].Fill(dy,nc)
                  if x<nSiPMs: 
                       if s==2 and smallSiPMchannel(x):  smallL+=1
                       else:    left+=1
                  else:           
                       if s==2 and smallSiPMchannel(x):  smallR+=1
                       else:   right+=1
              if left>0:
                    rc = h['resY_'+sdict[s]+'L'+str(s*10+plane)].Fill(dy,locEx[1])
                    rc = h['resX_'+sdict[s]+'L'+str(s*10+plane)].Fill(dx,locEx[0])
                    rc = h['gresY_'+sdict[s]+'L'+str(s*10+plane)].Fill(gdy,globEx[1])
                    rc = h['gresX_'+sdict[s]+'L'+str(s*10+plane)].Fill(gdx,globEx[0])
              if right>0:
                    rc = h['resY_'+sdict[s]+'R'+str(s*10+plane)].Fill(dy,locEx[1])
                    rc = h['resX_'+sdict[s]+'R'+str(s*10+plane)].Fill(dx,locEx[0])
                    rc = h['gresY_'+sdict[s]+'R'+str(s*10+plane)].Fill(gdy,globEx[1])
                    rc = h['gresX_'+sdict[s]+'R'+str(s*10+plane)].Fill(gdx,globEx[0])
              if s==2 and (smallL>0 or smallR>0): 
                     rc = h['resY_'+sdict[s]+'S'+str(s*10+plane)].Fill(dy,locEx[1])
                     rc = h['resX_'+sdict[s]+'S'+str(s*10+plane)].Fill(dx,locEx[0])
                     rc = h['gresY_'+sdict[s]+'S'+str(s*10+plane)].Fill(gdy,globEx[1])
                     rc = h['gresX_'+sdict[s]+'S'+str(s*10+plane)].Fill(gdx,globEx[0])
              dist = abs(dy)
              if aHit.isVertical() : dist = abs(dx)
              if dist<3.0:   # check channels
                  if aHit.isVertical():
                     dL  = locA[1]- locEx[1]
                     dR = locEx[1] - locB[1]
                  else:
                     dR = locA[0] - locEx[0]
                     dL =  locEx[0] - locB[0]
                  barName = sdict[s]+str(s*10+plane)+'_'+str(bar)
                  rc = h['nSiPMsL_'+barName].Fill(left,dL)
                  rc = h['nSiPMsR_'+barName].Fill(right,dR)
                  for x in S:
                      qcd = S[x]
                      if x<nSiPMs:
                         if s==2 and smallSiPMchannel(x): 
                               rc = h['signalSL_'+barName+'-c'+str(x)].Fill(qcd,dL)
                               SleftS+=qcd
                         else:
                               rc = h['signalL_'+barName+'-c'+str(x)].Fill(qcd,dL)
                               Sleft+=qcd
                      else: 
                         if s==2 and smallSiPMchannel(x): 
                               rc = h['signalSR_'+barName+'-c'+str(x-nSiPMs)].Fill(qcd,dR)
                               SrightS+=qcd
                         else:
                               rc = h['signalR_'+barName+'-c'+str(x-nSiPMs)].Fill(qcd,dR)
                               Sright+=qcd
                  rc = h['signalTL_'+barName].Fill(Sleft,dL)
                  rc = h['signalTR_'+barName].Fill(Sright,dR)
                  rc = h['signalTSL_'+barName].Fill(SleftS,dL)
                  rc = h['signalTSR_'+barName].Fill(SrightS,dR)

#   look at delta time vs track X, works only for horizontal planes.
                  allTDCs = map2Dict(aHit,'GetAllTimes')
                  if not aHit.isVertical():
                     dt    = aHit.GetDeltaT()
                     dtF  = aHit.GetFastDeltaT()
                     mtVeto  = aHit.GetImpactT() - T0Veto - (globPos[2] - Z0Veto)/u.speedOfLight
                     h['dtLRvsX_'+sdict[s]+str(s*10+plane)].Fill(xEx,dt*TDC2ns)
                     h['dtLRvsX_'+barName].Fill(xEx,dt*TDC2ns)
                     if (1 in allTDCs) and (9 in allTDCs):  h['dtF1LRvsX_'+barName].Fill(xEx,(allTDCs[1]-allTDCs[9])*TDC2ns)
                     h['dtfastLRvsX_'+barName].Fill(xEx,dtF*TDC2ns)
                     if Z0track>0:
                        tcor = (globPos[2] - Z0track)/deltaZ02 * ToF
                        mtTrack = aHit.GetImpactT() - T0track - tcor
                        h['atLRvsX_'+sdict[s]+str(s*10+plane)].Fill(xEx,mtTrack)
                        h['atLRvsX_'+barName].Fill(xEx,mtTrack)
                        if s <3:
                           for i in allTDCs:
                              dt = allTDCs[i]*TDC2ns-T0track - tcor
                              #print('debug 2',tcor,dt)
                              if i<nSiPMs: h['tvsXL_'+barName+'-c'+str(i)].Fill(xEx,dt)
                              else:
                                                 h['tvsXR_'+barName+'-c'+str(i-nSiPMs)].Fill(xEx,dt)

                     h['VetoatLRvsX_'+sdict[s]+str(s*10+plane)].Fill(xEx,mtVeto)
# QDC/TDC channel variations
                  if s==3: continue
                  meanL,meanR,nL,nR=0,0,0,0
                  t0Left   = 999
                  t0Right = 999
                  for i in allTDCs:
                      if i==4: t0Left   = allTDCs[i]
                      if i==12: t0Right = allTDCs[i]

                  for i in allTDCs:
                      if s==2 and smallSiPMchannel(i): continue
                      if  i < nSiPMs:  # left side
                          nL+=1
                          meanL+=allTDCs[i]
                      else:
                          nR+=1
                          meanR+=allTDCs[i]
                  for i in allTDCs:
                     if s==2 and smallSiPMchannel(i): continue
                     if i<nSiPMs and nL>0:
                          key =  sdict[s]+str(s*10+plane)+'_'+str(bar)+'-c'+str(i)
                          rc =                     h['sigmaTDCL_'+key].Fill( (allTDCs[i]-meanL/nL)*TDC2ns )
                          if t0Left<900: rc = h['TDCcalibL_'+key].Fill( (allTDCs[i]-t0Left)*TDC2ns )
                     elif not(i<nSiPMs) and nR>0:
                          key =  sdict[s]+str(s*10+plane)+'_'+str(bar)+'-c'+str(i-nSiPMs)
                          rc = h['sigmaTDCR_'+key].Fill((allTDCs[i]-meanR/nR)*TDC2ns )
                          if t0Right<900: rc = h['TDCcalibR_'+key].Fill( (allTDCs[i]-t0Right)*TDC2ns )

                  meanL,meanR,nL,nR=0,0,0,0
                  for i in S:
                      if s==2 and smallSiPMchannel(i): continue
                      if  i < nSiPMs:  # left side
                          nL+=1
                          meanL+=S[i]
                      else:
                          nR+=1
                          meanR+=S[i]
                  for i in S:
                     if s==2 and smallSiPMchannel(i): continue
                     if i<nSiPMs and nL>0:
                          key =  sdict[s]+str(s*10+plane)+'_'+str(bar)+'-c'+str(i)
                          rc = h['sigmaQDCL_'+key].Fill( (S[i]-meanL/nL) )
                     elif not(i<nSiPMs) and nR>0:
                          key =  sdict[s]+str(s*10+plane)+'_'+str(bar)+'-c'+str(i-nSiPMs)
                          rc = h['sigmaQDCR_'+key].Fill((S[i]-meanR/nR) )

 for s in [1,2]:
    for l in range(systemAndPlanes[s]):
        for side in ['L','R']:
           for bar in range(systemAndBars[s]):
              key = 'sigmaTDC'+side+'_'+sdict[s]+str(s*10+l)+'_'+str(bar)
              h[key]=h[key+'-c0'].Clone(key)
              h[key].Reset()
              for i in range(systemAndChannels[s][1]+systemAndChannels[s][0]):
                    h[key].Add(h[key+'-c'+str(i)])

 ut.writeHists(h,'MuFilterEff_run'+str(options.runNumber)+'.root')




def qdc_dist(Nev = options.nEvents, fit = "langau", title = "US_QDC_distributions"):
 # veto system 2 layers with 7 bars and 8 sipm channels on both ends
 # US system 5 layers with 10 bars and 8 sipm channels on both ends
 # DS system horizontal(3) planes, 60 bars, readout on both sides, single channel
 #                         vertical(4) planes, 60 bar, readout on top, single channel
 bin_min = 0.
 bin_max = 50.
 hist_list = {}
 hist_list_lr = {}
 for l in range(20, 25):
    #ut.bookHist(h,'sig_'+str(l),'signal / plane '+str(l),200,0.0,50.)
    h =  ROOT.TH1I("plane" + f"_{l}", "plane" + f"_{l}", 200, bin_min, bin_max)
    h_lr =  ROOT.TH2I("plane" + f"_{l}_lr", "plane" + f"_{l}_lr", 100,-0.5,200.,100,-0.5,200.)
    hist_list[l] = h
    hist_list_lr[l] = h_lr

 for l in range(30, 37):
    h_ds =  ROOT.TH1I("plane" + f"_{l}", "plane" + f"_{l}", 200, bin_min, bin_max)
    hist_list[l] = h_ds


 N=-1
 if Nev < 0 : Nev = eventTree.GetEntries()
 eventTree.GetEvent(0)

 for event in eventTree:
    N+=1
    if N%options.heartBeat == 0: print('event ',N,' ',time.ctime())
    if N>Nev: break
    for aHit in event.Digi_MuFilterHits:
        if not aHit.isValid(): continue
        detID = aHit.GetDetectorID()

        s = detID//10000 # det number (s = 1 -- veto, s = 2 -- US, s = 3 -- DS)
        l  = (detID%10000)//1000  # plane number
        bar = (detID%1000)

        if s != 3 : continue
        nSiPMs = aHit.GetnSiPMs()
        nSides  = aHit.GetnSides()
        qdc_value = ana.av_qdc(aHit)
        if qdc_value == -1:
            continue
        print(s*10 + l)
        hist_list[s*10 + l].Fill(qdc_value)
        q_l, q_r = ana.qdc_left_right(aHit)
        hist_list_lr[s*10 + l].Fill(q_l, q_r)

 # langau fit
#  c  = ROOT.TCanvas("US QDC distribution","US QDC distribution",0,0,1000,800)
#  c.Divide(2,3)
#  for i, plane in enumerate(hist_list.keys()):
#     #print(i)
#     c.cd(i+1)
#    #  if fit == "langau":
#    #    ana.fit_langau(hist_list[plane], str(plane), bin_min, bin_max)
#    #  else:
#    #    hist_list[plane].Fit("pol5")
#     hist_list_lr[plane].Draw("COLZ")
#  c.SaveAs(title + "_" + fit + ".root")

 File = ROOT.TFile.Open(f"{options.runNumber}_run.root", "RECREATE")
 def write_dict_to_file(dict_obj, File):
    for key in dict_obj.keys():
      File.WriteObject(dict_obj[key], dict_obj[key].GetTitle())
 
 hists_all = [hist_list]
 for H in hists_all:
   write_dict_to_file(H, File)
 File.Close()

def qdc_dist_per_bar(Nev = options.nEvents, fit = "langau", title = "US_QDC_distributions"):
 # veto system 2 layers with 7 bars and 8 sipm channels on both ends
 # US system 5 layers with 10 bars and 8 sipm channels on both ends
 # DS system horizontal(3) planes, 60 bars, readout on both sides, single channel
 #                         vertical(4) planes, 60 bar, readout on top, single channel
 bin_min = 0.
 bin_max = 50.
 hist_list = {}
 hist_list_lr = {}
 for l in range(20, 25):
    #ut.bookHist(h,'sig_'+str(l),'signal / plane '+str(l),200,0.0,50.)
    h =  ROOT.TH1I("plane" + f"_{l}", "plane" + f"_{l}", 200, bin_min, bin_max)
    hist_list[l] = h
    hist_list_lr[l] = {}
    for bar in range(10):
      h_lr =  ROOT.TH2I("plane" + f"_{l}_{bar}_lr", "plane" + f"_{l}_{bar}_lr", 100,-0.5,200.,100,-0.5,200.)
      hist_list_lr[l][bar] = h_lr
 N=-1
 if Nev < 0 : Nev = eventTree.GetEntries()
 eventTree.GetEvent(0)

 for event in eventTree:
    N+=1
    if N%options.heartBeat == 0: print('event ',N,' ',time.ctime())
    if N>Nev: break
    for aHit in event.Digi_MuFilterHits:
        if not aHit.isValid(): continue
        detID = aHit.GetDetectorID()

        s = detID//10000 # det number (s = 1 -- veto, s = 2 -- US, s = 3 -- DS)
        l  = (detID%10000)//1000  # plane number
        bar = (detID%1000)

        if s != 2: continue
        nSiPMs = aHit.GetnSiPMs()
        nSides  = aHit.GetnSides()
        qdc_value = ana.av_qdc(aHit)
        if qdc_value == -1:
            continue
        hist_list[s*10 + l].Fill(qdc_value)
        q_l, q_r = ana.qdc_left_right(aHit)
        hist_list_lr[s*10 + l][bar].Fill(q_l, q_r)

#  langau fit
#  c  = ROOT.TCanvas("US QDC distribution","US QDC distribution",0,0,1000,800)
#  c.Divide(2,3)
#  for i, plane in enumerate(hist_list.keys()):
#     #print(i)
#     c.cd(i+1)
#    #  if fit == "langau":
#    #    ana.fit_langau(hist_list[plane], str(plane), bin_min, bin_max)
#    #  else:
#    #    hist_list[plane].Fit("pol5")
#     hist_list_lr[plane].Draw("COLZ")
#  c.SaveAs(title + "_" + fit + ".root")
# langau fit
 c  = ROOT.TCanvas("US QDC distribution","US QDC distribution",0,0,1200,1200)
 c.Divide(3,4)
 for i, plane in enumerate(hist_list_lr[23].keys()):
    #print(i)
    c.cd(i+1)
   #  if fit == "langau":
   #    ana.fit_langau(hist_list[plane], str(plane), bin_min, bin_max)
   #  else:
   #    hist_list[plane].Fit("pol5")
    hist_list_lr[23][i].Draw("COLZ")
 c.SaveAs(title + "_" + fit + ".root")
 c.SaveAs(title + "_" + fit + ".pdf")

 File = ROOT.TFile.Open(f"{options.runNumber}_run.root", "RECREATE")
 def write_dict_to_file(dict_obj, File):
    for key in dict_obj.keys():
      File.WriteObject(dict_obj[key], dict_obj[key].GetTitle())
 
 #hists_all = [hist_list, hist_list_lr]
 hists_all = [hist_list_lr[key] for key in hist_list_lr.keys()]
 for H in hists_all:
   write_dict_to_file(H, File)
 File.Close()


def beamSpot():
   trackTask.ExecuteTask()
   Xbar = -10
   Ybar = -10
   for  aTrack in Reco_MuonTracks:
         state = aTrack.getFittedState()
         pos    = state.getPos()
         rc = h['bs'].Fill(pos.x(),pos.y())
         points = aTrack.getPoints()
         keys     = ROOT.std.vector('int')()
         detIDs = ROOT.std.vector('int')()
         ROOT.fixRoot(points, detIDs,keys,True)
         for k in range(keys.size()):
             #                                     m = p.getRawMeasurement()
             detID =detIDs[k] # m.getDetId()
             key = keys[k]          # m.getHitId()//1000 # for mufi
             aHit = eventTree.Digi_MuFilterHits[key]
             if aHit.GetDetectorID() != detID: continue # not a Mufi hit
             s = detID//10000
             l  = (detID%10000)//1000  # plane number
             bar = (detID%1000)
             if s>2: 
               l=2*l
               if bar>59:
                    bar=bar-60
                    if l<6: l+=1
             if s==3 and l%2==0: Ybar=bar
             if s==3 and l%2==1: Xbar=bar
             nSiPMs = aHit.GetnSiPMs()
             nSides  = aHit.GetnSides()
             for p in range(nSides):
                 c=bar*nSiPMs*nSides + p*nSiPMs
                 for i in range(nSiPMs):
                      signal = aHit.GetSignal(i+p*nSiPMs)
                      if signal > 0:
                           rc  = h['Tsig_'+str(s)+str(l)].Fill(signal)
         mom = state.getMom()
         slopeY= mom.X()/mom.Z()
         slopeX= mom.Y()/mom.Z()
         h['slopes'].Fill(slopeX,slopeY)
         if not Ybar<0 and not Xbar<0 and abs(slopeY)<0.01: rc = h['bsDS'].Fill(Xbar,Ybar)

def Mufi_hitMaps(Nev = options.nEvents):
 # veto system 2 layers with 7 bars and 8 sipm channels on both ends
 # US system 5 layers with 10 bars and 8 sipm channels on both ends
 # DS system horizontal(3) planes, 60 bars, readout on both sides, single channel
 #                         vertical(4) planes, 60 bar, readout on top, single channel
 for s in systemAndPlanes:
    for l in range(systemAndPlanes[s]):
       ut.bookHist(h,'hitmult_'+str(s*10+l),'hit mult / plane '+str(s*10+l),61,-0.5,60.5)
       ut.bookHist(h,'hit_'+str(s*10+l),'channel map / plane '+str(s*10+l),160,-0.5,159.5)
       if s==3:  ut.bookHist(h,'bar_'+str(s*10+l),'hit map / plane '+str(s*10+l),60,-0.5,59.5)
       else:       ut.bookHist(h,'bar_'+str(s*10+l),'hit map / plane '+str(s*10+l),10,-0.5,9.5)
       ut.bookHist(h,'sig_'+str(s*10+l),'signal / plane '+str(s*10+l),200,0.0,200.)
       if s==2:    ut.bookHist(h,'sigS_'+str(s*10+l),'signal / plane '+str(s*10+l),200,0.0,200.)
       ut.bookHist(h,'sigL_'+str(s*10+l),'signal / plane '+str(s*10+l),200,0.0,200.)
       ut.bookHist(h,'sigR_'+str(s*10+l),'signal / plane '+str(s*10+l),200,0.0,200.)
       ut.bookHist(h,'Tsig_'+str(s*10+l),'signal / plane '+str(s*10+l),200,0.0,200.)
       ut.bookHist(h,'occ_'+str(s*10+l),'channel occupancy '+str(s*10+l),100,0.0,200.)
       ut.bookHist(h,'occTag_'+str(s*10+l),'channel occupancy '+str(s*10+l),100,0.0,200.)

 ut.bookHist(h,'leftvsright_1','Veto hits in left / right',10,-0.5,9.5,10,-0.5,9.5)
 ut.bookHist(h,'leftvsright_2','US hits in left / right',10,-0.5,9.5,10,-0.5,9.5)
 ut.bookHist(h,'leftvsright_3','DS hits in left / right',2,-0.5,1.5,2,-0.5,1.5)
 ut.bookHist(h,'leftvsright_signal_1','Veto signal in left / right',100,-0.5,200.,100,-0.5,200.)
 ut.bookHist(h,'leftvsright_signal_2','US signal in left / right',100,-0.5,200.,100,-0.5,200.)
 ut.bookHist(h,'leftvsright_signal_3','DS signal in left / right',100,-0.5,200.,100,-0.5,200.)

 ut.bookHist(h,'dtime','delta event time; dt [ns]',100,0.0,1000.)
 ut.bookHist(h,'dtimeu','delta event time; dt [us]',100,0.0,1000.)
 ut.bookHist(h,'dtimem','delta event time; dt [ms]',100,0.0,1000.)

 ut.bookHist(h,'bs','beam spot',100,-100.,10.,100,0.,80.)
 ut.bookHist(h,'bsDS','beam spot',60,-0.5,59.5,60,-0.5,59.5)
 ut.bookHist(h,'slopes','track slopes',100,-0.1,0.1,100,-0.1,0.1)

 for s in range(1,4):
    ut.bookHist(h,sdict[s]+'Mult','QDCs vs nr hits',200,0.,800.,200,0.,300.)

 N=-1
 Tprev = 0
 if Nev < 0 : Nev = eventTree.GetEntries()
 eventTree.GetEvent(0)
 t0 =  eventTree.EventHeader.GetEventTime()/freq
 listOfHits = {1:[],2:[],3:[]}
 mult = {}
 for event in eventTree:
    N+=1
    if N%options.heartBeat == 0: print('event ',N,' ',time.ctime())
    if N>Nev: break
    withX = False
    planes = {}
    listOfHits[1].clear()
    listOfHits[2].clear()
    listOfHits[3].clear()
    for s in systemAndPlanes:
        for l in range(systemAndPlanes[s]):   mult[s*10+l]=0
    for aHit in event.Digi_MuFilterHits:
        if not aHit.isValid(): continue
        detID = aHit.GetDetectorID()
        if aHit.isVertical():     withX = True
        s = detID//10000
        l  = (detID%10000)//1000  # plane number
        bar = (detID%1000)
        if s>2:
               l=2*l
               if bar>59:
                    bar=bar-60
                    if l<6: l+=1
        mult[s*10+l]+=1
        key = s*100+l
        if not key in planes: planes[key] = {}
        sumSignal = map2Dict(aHit,'SumOfSignals')
        planes[key][bar] = [sumSignal['SumL'],sumSignal['SumR']]
        nSiPMs = aHit.GetnSiPMs()
        nSides  = aHit.GetnSides()
# check left/right
        allChannels = map2Dict(aHit,'GetAllSignals')
        for c in allChannels:
           listOfHits[s].append(allChannels[c])
        if nSides==2:
           Nleft    = 0
           Nright = 0
           Sleft    = 0
           Sright = 0
           for c in allChannels:
              if  nSiPMs > c:  # left side
                    Nleft+=1
                    Sleft+=allChannels[c]
              else:
                    Nright+=1
                    Sright+=allChannels[c]
           rc = h['leftvsright_'+str(s)].Fill(Nleft,Nright)
           rc = h['leftvsright_signal_'+str(s)].Fill(Sleft,Sright)
#
        for c in allChannels:
            channel = bar*nSiPMs*nSides + c
            rc = h['hit_'+str(s)+str(l)].Fill( int(channel))
            rc = h['bar_'+str(s)+str(l)].Fill(bar)
            if s==2 and smallSiPMchannel(c) : rc  = h['sigS_'+str(s)+str(l)].Fill(allChannels[c])
            elif c<nSiPMs: rc  = h['sigL_'+str(s)+str(l)].Fill(allChannels[c])
            else             :             rc  = h['sigR_'+str(s)+str(l)].Fill(allChannels[c])
            rc  = h['sig_'+str(s)+str(l)].Fill(allChannels[c])
        allChannels.clear()
#
    for s in listOfHits:
         nhits = len(listOfHits[s])
         qcdsum = 0
         for i in range(nhits):
             rc = h[sdict[s]+'Mult'].Fill(nhits, listOfHits[s][i])
    for s in systemAndPlanes:
        for l in range(systemAndPlanes[s]):   
           rc = h['hitmult_'+str(s*10+l)].Fill(mult[s*10+l])

    maxOneBar = True
    for key in planes:
        if len(planes[key]) > 2: maxOneBar = False
    if withX and maxOneBar:  beamSpot()
    
 S = {1:[1800,800,2,1],2:[1800,1500,2,3],3:[1800,1800,2,4]}
 for s in S:
   ut.bookCanvas(h,'hitmaps' +str(s),' ',S[s][0],S[s][1],S[s][2],S[s][3])
   ut.bookCanvas(h,'barmaps'+str(s),' ',S[s][0],S[s][1],S[s][2],S[s][3])
   ut.bookCanvas(h,'signal'    +str(s),' ',S[s][0],S[s][1],S[s][2],S[s][3])
   ut.bookCanvas(h,'Tsignal'   +str(s),' ',S[s][0],S[s][1],S[s][2],S[s][3])

   for l in range(systemAndPlanes[s]):
      n = l+1
      if s==3 and n==7: n=8
      tc = h['hitmaps'+str(s)].cd(n)
      tag = str(s)+str(l)
      h['hit_'+tag].Draw()
      tc = h['barmaps'+str(s)].cd(n)
      h['bar_'+tag].Draw()
      tc = h['signal'+str(s)].cd(n)
      h['sig_'+tag].Draw()
      tc = h['Tsignal'+str(s)].cd(n)
      h['Tsig_'+tag].Draw()

 ut.bookCanvas(h,'hitmult','hit multiplicities per plane',2000,1600,4,3)
 k=1
 for s in systemAndPlanes:
     for l in range(systemAndPlanes[s]):
           tc = h['hitmult'].cd(k)
           tc.SetLogy(1)
           k+=1
           rc = h['hitmult_'+str(s*10+l)].Draw()

 ut.bookCanvas(h,'VETO',' ',1200,1800,1,2)
 for l in range(2):
    tc = h['VETO'].cd(l+1)
    hname = 'hit_'+str(1)+str(l)
    h[hname].SetStats(0)
    h[hname].Draw()
    for n in range(7):
       x = (n+1)*16-0.5
       y = h['hit_'+str(1)+str(l)].GetMaximum()
       lname = 'L'+str(n)+hname
       h[lname] = ROOT.TLine(x,0,x,y)
       h[lname].SetLineColor(ROOT.kRed)
       h[lname].SetLineStyle(9)
       h[lname].Draw('same')

 ut.bookCanvas(h,'USBars',' ',1200,900,1,1)
 colours = {0:ROOT.kOrange,1:ROOT.kRed,2:ROOT.kGreen,3:ROOT.kBlue,4:ROOT.kMagenta,5:ROOT.kCyan,
                  6:ROOT.kAzure,7:ROOT.kPink,8:ROOT.kSpring}
 for i in range(5): 
       h['bar_2'+str(i)].SetLineColor(colours[i])
       h['bar_2'+str(i)].SetLineWidth(2)
       h['bar_2'+str(i)].SetStats(0)
 h['bar_20'].Draw()
 h['bar_21'].Draw('same')
 h['bar_22'].Draw('same')
 h['bar_23'].Draw('same')
 h['bar_24'].Draw('same')
 h['lbar2']=ROOT.TLegend(0.6,0.6,0.99,0.99)
 for i in range(5): 
    h['lbar2'].AddEntry(h['bar_2'+str(i)],'plane '+str(i+1),"f")
 h['lbar2'].Draw()
 for i in range(7): 
       h['hit_3'+str(i)].SetLineColor(colours[i])
       h['hit_3'+str(i)].SetLineWidth(2)
       h['hit_3'+str(i)].SetStats(0)
 h['hit_30'].Draw()
 for i in range(1,7):
     h['hit_3'+str(i)].Draw('same')
 h['lbar3']=ROOT.TLegend(0.6,0.6,0.99,0.99)
 for i in range(7): 
    h['lbar3'].AddEntry(h['hit_3'+str(i)],'plane '+str(i+1),"f")
 h['lbar3'].Draw()

 ut.bookCanvas(h,'LR',' ',1800,900,3,2)
 h['LR'].cd(1)
 h['leftvsright_'+str(1)].Draw('textBox')
 h['LR'].cd(2)
 h['leftvsright_'+str(2)].Draw('textBox')
 h['LR'].cd(3)
 h['leftvsright_'+str(3)].Draw('textBox')
 h['LR'].cd(4)
 h['leftvsright_signal_1'].SetMaximum(h['leftvsright_signal_1'].GetBinContent(10,10))
 h['leftvsright_signal_2'].SetMaximum(h['leftvsright_signal_2'].GetBinContent(10,10))
 h['leftvsright_signal_3'].SetMaximum(h['leftvsright_signal_3'].GetBinContent(10,10))
 h['leftvsright_signal_'+str(1)].Draw('colz')
 h['LR'].cd(5)
 h['leftvsright_signal_'+str(2)].Draw('colz')
 h['LR'].cd(6)
 h['leftvsright_signal_'+str(3)].Draw('colz')

 ut.bookCanvas(h,'LRinEff',' ',1800,450,3,1)
 for s in range(1,4):
   h['lLRinEff'+str(s)]=ROOT.TLegend(0.6,0.54,0.99,0.93)
   name = 'leftvsright_signal_'+str(s)
   h[name+'0Y'] = h[name].ProjectionY(name+'0Y',1,1)
   h[name+'0X'] = h[name].ProjectionX(name+'0X',1,1)
   h[name+'1X'] = h[name].ProjectionY(name+'1Y')
   h[name+'1Y'] = h[name].ProjectionX(name+'1X')
   tc = h['LRinEff'].cd(s)
   tc.SetLogy()
   h[name+'0X'].SetStats(0)
   h[name+'0Y'].SetStats(0)
   h[name+'1X'].SetStats(0)
   h[name+'1Y'].SetStats(0)
   h[name+'0X'].SetLineColor(ROOT.kRed)
   h[name+'0Y'].SetLineColor(ROOT.kGreen)
   h[name+'1X'].SetLineColor(ROOT.kMagenta)
   h[name+'1Y'].SetLineColor(ROOT.kCyan)
   h[name+'0X'].SetMaximum(max(h[name+'1X'].GetMaximum(),h[name+'1Y'].GetMaximum()))
   h[name+'0X'].Draw()
   h[name+'0Y'].Draw('same')
   h[name+'1X'].Draw('same')
   h[name+'1Y'].Draw('same')
   # Fill(Sleft,Sright)
   h['lLRinEff'+str(s)].AddEntry(h[name+'0X'],'left with no signal right',"f")
   h['lLRinEff'+str(s)].AddEntry(h[name+'0Y'],'right with no signal left',"f")
   h['lLRinEff'+str(s)].AddEntry(h[name+'1X'],'left all',"f")
   h['lLRinEff'+str(s)].AddEntry(h[name+'1Y'],'right all',"f")
   h['lLRinEff'+str(s)].Draw()

 ut.bookCanvas(h,'signalUSVeto',' ',1200,1600,3,7)
 s = 1
 l = 1
 for plane in range(2):
     for side in ['L','R','S']:
         tc = h['signalUSVeto'].cd(l)
         l+=1
         if side=='S': continue
         h['sig'+side+'_'+str( s*10+plane)].Draw()
 s=2
 for plane in range(5):
     for side in ['L','R','S']:
         tc = h['signalUSVeto'].cd(l)
         l+=1
         h['sig'+side+'_'+str( s*10+plane)].Draw()
 ut.bookCanvas(h,'signalDS',' ',900,1600,2,7)
 s = 3
 l = 1
 for plane in range(7):
     for side in ['L','R']:
         tc = h['signalDS'].cd(l)
         l+=1
         h['sig'+side+'_'+str( s*10+plane)].Draw()


 for canvas in ['signalUSVeto','LR','USBars']:
      h[canvas].Update()
      myPrint(h[canvas],canvas)
 for canvas in ['hitmaps','barmaps','signal','Tsignal']:
      for s in range(1,4):
         h[canvas+str(s)].Update()
         myPrint(h[canvas+str(s)],canvas+sdict[s])

def pid_ana_mc(Nev = options.nEvents, fit = "langau", title = "US_QDC_distributions"):
 # veto system 2 layers with 7 bars and 8 sipm channels on both ends
 # US system 5 layers with 10 bars and 8 sipm channels on both ends
 # DS system horizontal(3) planes, 60 bars, readout on both sides, single channel
 #                         vertical(4) planes, 60 bar, readout on top, single channel

 bin_min = 0.
 bin_max = 100.
 hist_list = {}
 hist_list_mc = {}
 hist_list_mc_full = {}
 list_of_selected_pid = [11, 2212, 2112, 22, 211]
 for l in range(20, 25):
    #ut.bookHist(h,'sig_'+str(l),'signal / plane '+str(l),200,0.0,50.)
    h =  ROOT.TH1I("plane" + f"_{l}", "plane" + f"_{l}", 200, bin_min, bin_max)
    h_mc_full =  ROOT.TH1I("plane" + f"_{l}_mc_full", "plane_mc_full" + f"_{l}", 200, 0, 1e-2)
    hist_list_mc_full[l] = h_mc_full
    hist_list[l] = h
    for pid in list_of_selected_pid:
      h_mc =  ROOT.TH1I("plane" + f"_{l}_mc_{pid}", "plane_mc" + f"_{l}_{pid}", 200, 0, 1e-2)
      hist_list_mc[f"_{l}_mc_{pid}"] = h_mc

 N=-1
 if Nev < 0 : Nev = eventTree.GetEntries()
 eventTree.GetEvent(0)
 
 list_of_pid = {}
 for event in eventTree:
    L = event.Digi_MuFilterHits2MCPoints[0]
    w = event.MCTrack[0].GetWeight() # assume these are muon background events.
    list_of_detectors = []
    N+=1
    if N%2000 == 0: print('event ',N,' ',time.ctime())
    if N>Nev: break
    for i, aHit in enumerate(event.Digi_MuFilterHits):
        if not aHit.isValid(): continue
        detID = aHit.GetDetectorID()

        s = detID//10000 # det number (s = 1 -- veto, s = 2 -- US, s = 3 -- DS)
        l  = (detID%10000)//1000  # plane number
        bar = (detID%1000)
        if s != 2: continue
        nSiPMs = aHit.GetnSiPMs()
        nSides  = aHit.GetnSides()

        points = L.wList(detID)
        for p in points:
          pid = event.MuFilterPoint[p[0]].PdgCode()
          hist_list_mc_full[s*10 + l].Fill(event.MuFilterPoint[p[0]].GetEnergyLoss())
          if np.abs(pid) in list_of_selected_pid:
            hist_list_mc[f"_{s*10 + l}_mc_{np.abs(pid)}"].Fill(event.MuFilterPoint[p[0]].GetEnergyLoss())
          if pid in list_of_pid:
            list_of_pid[pid] += 1 
            continue
          list_of_pid[pid] = 0 

        qdc_value = ana.av_qdc(aHit)
        if qdc_value == -1:
            continue
        
        hist_list[s*10 + l].Fill(qdc_value)

        
 print(list_of_pid)
 # langau fit
#  c  = ROOT.TCanvas("US QDC distribution","US QDC distribution",0,0,1000,800)
#  c.Divide(2,3)
#  for i, plane in enumerate(hist_list.keys()):
#     #print(i)
#     c.cd(i+1)
#     if fit == "langau":
#       ana.fit_langau(hist_list[plane], str(plane), bin_min, bin_max)
#     else:
#       hist_list[plane].Fit("pol5")
#     hist_list[plane].Draw()
#  c.SaveAs(title + "_" + fit + ".root")

#  c1  = ROOT.TCanvas("US QDC distribution MC","US QDC distribution MC",0,0,1000,800)
#  c1.Divide(2,3)
#  for i, plane in enumerate(hist_list_mc.keys()):
#     #print(i)
#     c1.cd(i+1)
#     if fit == "langau":
#       ana.fit_langau(hist_list_mc[plane], str(plane), bin_min, bin_max)
#     else:
#       hist_list_mc[plane].Fit("pol5")
#     hist_list_mc[plane].Draw()
#  c1.SaveAs(title + "_" + fit + "_mc_enloss.root")


#  c2  = ROOT.TCanvas("US QDC distribution MC","US QDC distribution MC",0,0,1000,800)
#  c2.Divide(2,3)
#  for i, plane in enumerate(hist_list_mc_full.keys()):
#     #print(i)
#     c2.cd(i+1)
#     if fit == "langau":
#       ana.fit_langau(hist_list_mc_full[plane], str(plane), bin_min, bin_max)
#     else:
#       hist_list_mc_full[plane].Fit("pol5")
#     hist_list_mc_full[plane].Draw()
#  c2.SaveAs(title + "_" + fit + "_mc_enloss_full.root")


 File = ROOT.TFile.Open(f"{options.runNumber}_run.root", "RECREATE")

 def write_dict_to_file(dict_obj, File):
    for key in dict_obj.keys():
      File.WriteObject(dict_obj[key], dict_obj[key].GetTitle())
 
 hists_all = [hist_list, hist_list_mc, hist_list_mc_full]
 for H in hists_all:
   write_dict_to_file(H, File)
 File.Close()

def pid_create_output(Nev = options.nEvents, filename = "default", eos = True):
 # veto system 2 layers with 7 bars and 8 sipm channels on both ends
 # US system 5 layers with 10 bars and 8 sipm channels on both ends
 # DS system horizontal(3) planes, 60 bars, readout on both sides, single channel
 #                         vertical(4) planes, 60 bar, readout on top, single channel

 list_of_selected_pid = [11, 2212, 2112, 22, 211]


 N=-1
 if Nev < 0 : Nev = eventTree.GetEntries()
 eventTree.GetEvent(0)
 
 list_of_pid = {}
 output = {"pid":[], "ev_num":[], "en_dep":[], "en_dep_norm":[], "hit_qdc":[], "event_ID":[], "det_id":[]}
 for event in eventTree:
    if not ana.OneHitUS1(eventTree.Digi_MuFilterHits): continue  
    L = event.Digi_MuFilterHits2MCPoints[0]
    w = event.MCTrack[0].GetWeight() # assume these are muon background events.
    list_of_detectors = []
    N+=1
    if N%2000 == 0: print('event ',N,' ',time.ctime())
    if N>Nev: break
    for i, aHit in enumerate(event.Digi_MuFilterHits):
        if not aHit.isValid(): continue
        detID = aHit.GetDetectorID()

        s = detID//10000 # det number (s = 1 -- veto, s = 2 -- US, s = 3 -- DS)
        l  = (detID%10000)//1000  # plane number
        bar = (detID%1000)
        if s != 2: continue
        nSiPMs = aHit.GetnSiPMs()
        nSides  = aHit.GetnSides()

        qdc_value = ana.av_qdc(aHit)
        if qdc_value == -1:
            continue
        points = L.wList(detID)
        for p in points:
          pid = event.MuFilterPoint[p[0]].PdgCode()
          if (np.abs(pid) in list_of_selected_pid) or True:
            output["pid"].append(np.abs(pid))
            output["en_dep"].append(event.MuFilterPoint[p[0]].GetEnergyLoss())
            output["en_dep_norm"].append(p[1])
            output["hit_qdc"].append(qdc_value)
            output["ev_num"].append(N)
            output["event_ID"].append((N+1)*100000 + i)
            output["det_id"].append(s*100 + l*10 + bar)
 types = [np.int32, np.int32, np.float64,np.float64,np.float64, np.int64, np.int64]
 df = ROOT.RDF.MakeNumpyDataFrame({key: np.array(output[key], dtype=types[k]) for k, key in enumerate(output.keys())})
        # ... or print the content
        #df.Display().Print()
        # ... or save the data as a ROOT file
 if eos:
   output = "root://eosuser.cern.ch//eos/user/u/ursovsnd/private/SND_Data/qdc_study/preprocessed_data/"
 else:
   output = ""
 df.Snapshot('tree', f'{output}{filename}.root')


# def exec_for_event(num):
#     #print(num, "!!!")
#     output = {"pid":[], "ev_num":[], "en_dep":[], "en_dep_norm":[], "hit_qdc":[], "event_ID":[]}
#     event = eventTree
#     event.GetEvent(num)
#     L = event.Digi_MuFilterHits2MCPoints[0]
#     w = event.MCTrack[0].GetWeight() # assume these are muon background events.
#     list_of_detectors = []
#     for i, aHit in enumerate(event.Digi_MuFilterHits):
#         if not aHit.isValid(): continue
#         detID = aHit.GetDetectorID()

#         s = detID//10000 # det number (s = 1 -- veto, s = 2 -- US, s = 3 -- DS)
#         l  = (detID%10000)//1000  # plane number
#         bar = (detID%1000)
#         if s != 2: continue
#         nSiPMs = aHit.GetnSiPMs()
#         nSides  = aHit.GetnSides()

#         qdc_value = ana.av_qdc(aHit)
#         if qdc_value == -1:
#             continue
#         points = L.wList(detID)
#         for p in points:
#           pid = event.MuFilterPoint[p[0]].PdgCode()
#           if True:
#             output["pid"].append(np.abs(pid))
#             output["en_dep"].append(event.MuFilterPoint[p[0]].GetEnergyLoss())
#             output["en_dep_norm"].append(p[1])
#             output["hit_qdc"].append(qdc_value)
#             output["ev_num"].append(num)
#             output["event_ID"].append((num)*100000 + i)
#     return output
         
# def concat_dicts(dict1, dict2):
#    for key in dict1.keys():
#       dict1[key] += dict2[key]
#    return dict1

# def test_multiprocessing(Nev = options.nEvents):
#    output = {"pid":[], "ev_num":[], "en_dep":[], "en_dep_norm":[], "hit_qdc":[], "event_ID":[]}
#    LIST = [k for k in range(Nev+1)]
#    #print(LIST)
#    #event_slice = [event for k, event in enumerate(eventTree) if k < Nev]
#    with Pool(processes=2) as p:
#       #p.starmap(exec_for_event, [(k, output) for k in range(10)])
#       #concat_dicts(output, dict(p.map(exec_for_event, [k for k in range(Nev)])))
#       # out = p.map(exec_for_event, [k for k in range(Nev)])[0]
#       # print(type(out))
#       # print(out)
#       #out = p.map(exec_for_event, LIST)[0]
#       for x in p.map(exec_for_event, LIST):
#          concat_dicts(output, x)
#       #print(out["ev_num"][0])
#       #concat_dicts(output, out)
#    types = [np.int32, np.int32, np.float64,np.float64,np.float64, np.int64]
#    df = ROOT.RDF.MakeNumpyDataFrame({key: np.array(output[key], dtype=types[k]) for k, key in enumerate(output.keys())})
#         # ... or print the content
#         #df.Display().Print()
#         # ... or save the data as a ROOT file
#    #output = "root://eosuser.cern.ch//eos/user/u/ursovsnd/private/SND_Data/qdc_study/preprocessed_data/"
#    output = ""
#    df.Snapshot('tree', f'{output}' + 'test_multi_1.root')


#test_multiprocessing(1000)

pid_create_output(-1, "output_300GeV_exp_full_det_id_onehitUS1_1", eos = True)

#pid_ana_mc(5000, "langau", "US_QDC_distributions_MC")


#qdc_dist(1000, "langau", "DS_QDC_run_90")
#Mufi_hitMaps(1000)
#Mufi_Efficiency(1000)
#qdc_dist_per_bar(1000000, "langau", "US_QDC_distributions_run_54")
