#python -i MufiScifi.py -r 46 -p /eos/experiment/sndlhc/testbeam/MuFilter/TB_data_commissioning/sndsw/ -g geofile_sndlhc_H6.root

import ROOT,os,subprocess
import atexit
import time
import ctypes
from array import array
import numpy as np
import AnalysisFunctions as muAna
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

def Scifi_hitMaps(Nev=options.nEvents):
 scifi = geo.modules['Scifi']
 A=ROOT.TVector3()
 B=ROOT.TVector3()
 
 for s in range(10):
    ut.bookHist(h,'posX_'+str(s),'x',2000,-100.,100.)
    ut.bookHist(h,'posY_'+str(s),'y',2000,-100.,100.)
    if s%2==1: ut.bookHist(h,'mult_'+str(s),'mult vertical station '+str(s//2+1),100,-0.5,99.5)
    else: ut.bookHist(h,'mult_'+str(s),'mult horizontal station '+str(s//2+1),100,-0.5,99.5)
 for mat in range(30):
    ut.bookHist(h,'mat_'+str(mat),'hit map / mat',512,-0.5,511.5)
    ut.bookHist(h,'sig_'+str(mat),'signal / mat',200,-50.0,150.)
    ut.bookHist(h,'tdc_'+str(mat),'tdc / mat',200,-1.,4.)
 N=-1
 if Nev < 0 : Nev = eventTree.GetEntries()
 for event in eventTree:
    N+=1
    if N%options.heartBeat == 0: print('event ',N,' ',time.ctime())
    if N>Nev: break
    mult = [0]*10
    for aHit in event.Digi_ScifiHits:
        if not aHit.isValid(): continue
        X =  Scifi_xPos(aHit.GetDetectorID())
        rc = h['mat_'+str(X[0]*3+X[1])].Fill(X[2])
        rc  = h['sig_'+str(X[0]*3+X[1])].Fill(aHit.GetSignal(0))
        rc  = h['tdc_'+str(X[0]*3+X[1])].Fill(aHit.GetTime(0))
        scifi.GetSiPMPosition(aHit.GetDetectorID(),A,B)
        if aHit.isVertical(): rc = h['posX_'+str(X[0])].Fill(A[0])
        else:                     rc = h['posY_'+str(X[0])].Fill(A[1])
        mult[X[0]]+=1
    for s in range(10):
       rc = h['mult_'+str(s)].Fill(mult[s])

 ut.bookCanvas(h,'hitmaps',' ',1024,768,6,5)
 ut.bookCanvas(h,'signal',' ',1024,768,6,5)
 ut.bookCanvas(h,'tdc',' ',1024,768,6,5)
 for mat in range(30):
    tc = h['hitmaps'].cd(mat+1)
    A = h['mat_'+str(mat)].GetSumOfWeights()/512.
    if h['mat_'+str(mat)].GetMaximum()>10*A: h['mat_'+str(mat)].SetMaximum(10*A)
    h['mat_'+str(mat)].Draw()
    tc = h['signal'].cd(mat+1)
    h['sig_'+str(mat)].Draw()
    tc = h['tdc'].cd(mat+1)
    h['tdc_'+str(mat)].Draw()

 ut.bookCanvas(h,'positions',' ',2048,768,5,2)
 ut.bookCanvas(h,'mult',' ',2048,768,5,2)
 for s in range(5):
    tc = h['positions'].cd(s+1)
    h['posY_'+str(2*s)].Draw()
    tc = h['positions'].cd(s+6)
    h['posX_'+str(2*s+1)].Draw()
    tc = h['mult'].cd(s+1)
    tc.SetLogy(1)
    h['mult_'+str(2*s)].Draw()
    tc = h['mult'].cd(s+6)
    tc.SetLogy(1)
    h['mult_'+str(2*s+1)].Draw()

 for canvas in ['hitmaps','signal','mult']:
      h[canvas].Update()
      myPrint(h[canvas],"Scifi-"+canvas)

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
       ut.bookHist(h,'sig_'+str(s*10+l),'signal / plane '+str(s*10+l),200,0.0,50.)
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

 # mult left and right
 ut.bookHist(h,'mult_signal','mult signal',100,0,1000)

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
           rc = h['mult_signal'].Fill(Sleft,Sright)
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


 ### mult signal
 ut.bookCanvas(h,'mult_signal1', ' ', 1200,900,1,1)
 h['mult_signal'].Draw("C")
 h['mult_signal1'].Update()
 myPrint(h['mult_signal1'],'mult_signal')



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
def smallVsLargeSiPMs(Nev=-1):
 S = 2
 for l in range(systemAndPlanes[S]):
    ut.bookHist(h,'SVSl_'+str(l),'QDC large vs small sum',200,0.,200.,200,0.,200.)
    ut.bookHist(h,'sVSl_'+str(l),'QDC large vs small average',200,0.,200.,200,0.,200.)
    for side in ['L','R']:
          for i1 in range(7):
             for i2 in range(i1+1,8):
               tag=''
               if S==2 and smallSiPMchannel(i1): tag = 's'+str(i1)
               else:                              tag = 'l'+str(i1)
               if S==2 and smallSiPMchannel(i2): tag += 's'+str(i2)
               else:                              tag += 'l'+str(i2)
               ut.bookHist(h,'cor'+tag+'_'+side+str(l),'QDC channel i vs channel j',200,0.,200.,200,0.,200.)
               for bar in range(systemAndBars[S]):
                     ut.bookHist(h,'cor'+tag+'_'+side+str(l)+str(bar),'QDC channel i vs channel j',200,0.,200.,200,0.,200.)

 N=-1
 if Nev < 0 : Nev = eventTree.GetEntries()
 for event in eventTree:
    N+=1
    if N%options.heartBeat == 0: print('event ',N,' ',time.ctime())
    if N>Nev: break
    for aHit in event.Digi_MuFilterHits:
        if not aHit.isValid(): continue
        detID = aHit.GetDetectorID()
        s = detID//10000
        bar = (detID%1000)
        if s!=2: continue
        l  = (detID%10000)//1000  # plane number
        sumL,sumS,SumL,SumS = 0,0,0,0
        allChannels = map2Dict(aHit,'GetAllSignals',mask=False)
        nS = 0
        nL = 0
        for c in allChannels:
            if s==2 and smallSiPMchannel(c) : 
                sumS+= allChannels[c]
                nS += 1
            else:                                              
                sumL+= allChannels[c]
                nL+=1
        if nL>0: SumL=sumL/nL
        if nS>0: SumS=sumS/nS
        rc = h['sVSl_'+str(l)].Fill(SumS,SumL)
        rc = h['SVSl_'+str(l)].Fill(sumS/4.,sumL/12.)
#
        for side in ['L','R']:
          offset = 0
          if side=='R': offset = 8
          for i1 in range(offset,offset+7):
             if not i1 in allChannels: continue
             qdc1 = allChannels[i1]
             for i2 in range(i1+1,offset+8):
               if not (i2) in allChannels: continue
               if s==2 and smallSiPMchannel(i1): tag = 's'+str(i1-offset)
               else: tag = 'l'+str(i1-offset)
               if s==2 and smallSiPMchannel(i2): tag += 's'+str(i2-offset)
               else: tag += 'l'+str(i2-offset)
               qdc2 = allChannels[i2]
               rc = h['cor'+tag+'_'+side+str(l)].Fill(qdc1,qdc2)
               rc = h['cor'+tag+'_'+side+str(l)+str(bar)].Fill(qdc1,qdc2)
        allChannels.clear()

 ut.bookCanvas(h,'TSL','',1800,1400,3,2)
 ut.bookCanvas(h,'STSL','',1800,1400,3,2)
 for l in range(systemAndPlanes[S]):
     tc = h['TSL'].cd(l+1)
     tc.SetLogz(1)
     aHist = h['sVSl_'+str(l)]
     aHist.SetTitle(';small SiPM QCD av:large SiPM QCD av')
     nmax = aHist.GetBinContent(aHist.GetMaximumBin())
     aHist.SetMaximum( 0.1*nmax )
     tc = h['sVSl_'+str(l)].Draw('colz')
 myPrint(h['TSL'],"largeSiPMvsSmallSiPM")
 for l in range(systemAndPlanes[S]):
     tc = h['STSL'].cd(l+1)
     tc.SetLogz(1)
     aHist = h['SVSl_'+str(l)]
     aHist.SetTitle(';small SiPM QCD sum/2:large SiPM QCD sum/6')
     nmax = aHist.GetBinContent(aHist.GetMaximumBin())
     aHist.SetMaximum( 0.1*nmax )
     tc = h['SVSl_'+str(l)].Draw('colz')
 myPrint(h['STSL'],"SumlargeSiPMvsSmallSiPM")

 for l in range(systemAndPlanes[S]):
    for side in ['L','R']:
      ut.bookCanvas(h,'cor'+side+str(l),'',1800,1400,7,4)
      k=1
      for i1 in range(7):
             for i2 in range(i1+1,8):
               tag=''
               if S==2 and smallSiPMchannel(i1): tag = 's'+str(i1)
               else:                              tag = 'l'+str(i1)
               if S==2 and smallSiPMchannel(i2): tag += 's'+str(i2)
               else:                              tag += 'l'+str(i2)
               tc = h['cor'+side+str(l)].cd(k)
               for bar in range(systemAndBars[S]):
                    if bar == 0: h['cor'+tag+'_'+side+str(l)+str(bar)].Draw('colz')
                    else: h['cor'+tag+'_'+side+str(l)+str(bar)].Draw('colzsame')
               k+=1
      myPrint(h['cor'+side+str(l)],'QDCcor'+side+str(l))


def makeIndividualPlots(run=options.runNumber):
   ut.bookCanvas(h,'dummy','',900,800,1,1)
   if not "run"+str(run) in os.listdir(): os.system("mkdir run"+str(run))
   for l in range(5):
       for side in ['L','R']:
           f=ROOT.TFile('QDCcor'+side+str(l)+'-run'+str(run)+'.root')
           tcanv = f.Get('cor'+side+str(l))
           for pad in tcanv.GetListOfPrimitives():
              if not hasattr(pad,"GetListOfPrimitives"): continue
              for aHist in pad.GetListOfPrimitives():
                 if not aHist.ClassName() == 'TH2D': continue
                 hname = aHist.GetName()
                 tmp = hname.split('_')
                 bar = tmp[1][2]
                 pname = 'corUS'+str(l)+'-'+str(bar)+side+'_'+tmp[0][3:]
                 aHist.SetDirectory(ROOT.gROOT)
                 ROOT.gROOT.cd()
                 tc=h['dummy'].cd()
                 tc.SetLogz(1)
                 aHist.Draw('colz')
                 tc.Update()
                 stats = aHist.FindObject('stats')
                 stats.SetOptStat(11)
                 stats.SetX1NDC(0.15)
                 stats.SetY1NDC(0.75)
                 stats.SetX2NDC(0.35)
                 stats.SetY2NDC(0.88)
                 h['dummy'].Update()
                 tc.Print('run'+str(run)+'/'+pname+'.png')
                 tc.Print('run'+str(run)+'/'+pname+'.root')
   #os.system("convert -delay 120 -loop 0 run"+str(run)+"/corUS*.png corUS-"+str(run)+".gif")

def makeQDCcorHTML(run=options.runNumber):
   F = ROOT.TFile.Open('QDCcorrelations-run'+str(run)+'.root','recreate')
   for l in range(5):
       for side in ['L','R']:
           key = side+str(l)
           f=ROOT.TFile('QDCcor'+key+'-run'+str(run)+'.root')
           tcanv = f.Get('cor'+key).Clone()
           F.mkdir(key)
           F.cd(key)
           tc.Write()
def makeLogVersion(run):
   for l in range(5):
      for side in ['L','R']:
         fname = 'QDCcor'+side+str(l)+'-run'+str(run)
         f=ROOT.TFile('QDCcor'+side+str(l)+'-run'+str(run)+'.root')
         c = 'cor'+side+str(l)
         h['X'] = f.Get(c).Clone(c)
         for pad in h['X'].GetListOfPrimitives():
             pad.SetLogz(1)
         h['X'].Draw()
         h['X'].Print(fname+'.pdf')


def eventTime(Nev=options.nEvents):
 if Nev < 0 : Nev = eventTree.GetEntries()
 ut.bookHist(h,'Etime','delta event time; dt [s]',100,0.0,1.)
 ut.bookHist(h,'EtimeZ','delta event time; dt [ns]',1000,0.0,10000.)
 ut.bookCanvas(h,'T',' ',1024,3*768,1,3)
 
# need to make extra gymnastiques since absolute time is missing
 Ntinter = []
 N = 0
 for f in eventTree.GetListOfFiles():
    dN =  f.GetEntries()
    rc = eventTree.GetEvent(N)
    t0 = eventTree.EventHeader.GetEventTime()/freq
    rc = eventTree.GetEvent(N+dN-1)
    tmax = eventTree.EventHeader.GetEventTime()/freq
    Ntinter.append([t0,tmax])
    N+=dN

 Tduration = 0
 for x in Ntinter:
    Tduration += (x[1]-x[0])
 tsep = 3600.
 t0 =  Ntinter[0][0]
 tmax = Tduration+(tsep*(len(Ntinter)-1)) 
 #nbins = int(tmax-t0)
 nbins = 500
 yunit = "events per %5.0F s"%( (tmax-t0)/nbins)
 if 'time' in h: h.pop('time').Delete()
 #ut.bookHist(h,'time','elapsed time; t [s];'+yunit,nbins,0,tmax-t0)
 ut.bookHist(h,'time','elapsed time; t [s];'+yunit,nbins,0,tmax-t0 - 3000)
 N=-1
 Tprev  = 0
 Toffset = 0
 for event in eventTree:
    N+=1
    if N>Nev: break
    T   = event.EventHeader.GetEventTime()
    dT = T-Tprev
    if N>0 and dT >0:
           rc = h['Etime'].Fill( dT/freq )
           rc = h['EtimeZ'].Fill( dT*1E9/freq )
           rc = h['time'].Fill( (T+Toffset)/freq-t0 )
    elif dT<0: 
           Toffset+=tsep*freq+Tprev
           rc = h['time'].Fill( (T+Toffset)/freq-t0 )
    else: rc = h['time'].Fill( T/freq-t0 ) # very first event
    Tprev = T

 tc = h['T'].cd(1)

 bin_max = 0
 for i in range(100):
   if h['time'].GetBinContent(i) > bin_max:
      bin_max = h['time'].GetBinContent(i)
 time_selected = []
 for i in range(nbins):
   event_num = h['time'].GetBinContent(i)
   if event_num >= 0.0*bin_max and event_num < 1.5*bin_max:
      time_selected.append(h['time'].GetBinCenter(i))
   else:
      h['time'].SetBinContent(i, 0)
      #pass


 
 h['time'].SetStats(0)
 h['time'].Draw()
 tend = 0
 for x in Ntinter:
    tend += x[1]+tsep/2.
    m = str(int(tend))
    h['line'+m]=ROOT.TLine(tend,0,tend,h['time'].GetMaximum())
    h['line'+m].SetLineColor(ROOT.kRed)
    h['line'+m].Draw()
    tend += tsep/2.
 tc = h['T'].cd(2)
 tc.SetLogy(1)
 h['EtimeZ'].Draw()
 rc = h['EtimeZ'].Fit('expo','S','',0.,250.)
 h['T'].Update()
 stats = h['EtimeZ'].FindObject('stats')
 stats.SetOptFit(1111111)
 tc = h['T'].cd(3)
 tc.SetLogy(1)
 h['Etime'].Draw()
 rc = h['Etime'].Fit('expo','S')
 h['T'].Update()
 stats = h['Etime'].FindObject('stats')
 stats.SetOptFit(1111111)
 h['T'].Update()
 myPrint(h['T'],'time')

def TimeStudy(Nev=options.nEvents,withDisplay=False):
 if Nev < 0 : Nev = eventTree.GetEntries()
 ut.bookHist(h,'Vetotime','time',1000,0.,50.)
 ut.bookHist(h,'UStime','time',1000,0.,50.)
 ut.bookHist(h,'DStime','time',1000,0.,50.)
 ut.bookHist(h,'Stime','time',1000,0.,50.)
 ut.bookHist(h,'SvsDStime','; mean Scifi time [ns];mean Mufi time [ns]',100,0.,50.,100,0.,50.)
 ut.bookHist(h,'VEvsUStime','; mean US time [ns];mean VE time [ns]',100,0.,50.,100,0.,50.)
 ut.bookCanvas(h,'T','',900,1200,1,2)
 tc = h['T'].cd(1)
 h['Vetotime'].SetLineColor(ROOT.kOrange)
 h['UStime'].SetLineColor(ROOT.kGreen)
 h['DStime'].SetLineColor(ROOT.kRed)
 N=-1
 for event in eventTree:
   N+=1
   if N>Nev: break
   h['UStime'].Reset()
   h['DStime'].Reset()
   h['Vetotime'].Reset()
   h['Stime'].Reset()
   for aHit in eventTree.Digi_MuFilterHits:
     T = aHit.GetAllTimes()
     s = aHit.GetDetectorID()//10000
     for x in T:
       t = x.second*TDC2ns
       if t>0: 
           if s==1: rc = h['Vetotime'].Fill(t)
           if s==2: rc = h['UStime'].Fill(t)
           if s==3: rc = h['DStime'].Fill(t)
   stations = {}
   for aHit in eventTree.Digi_ScifiHits:
      t = aHit.GetTime()*TDC2ns
      rc = h['Stime'].Fill(t)
      stations[aHit.GetDetectorID()//1000000] = 1
   if len(stations)>3:
       rc = h['SvsDStime'].Fill(h['Stime'].GetMean(),h['DStime'].GetMean())
       rc = h['VEvsUStime'].Fill(h['UStime'].GetMean(),h['Vetotime'].GetMean())
   if withDisplay:
     tc = h['T'].cd(1)
     h['UStime'].Draw()
     h['DStime'].Draw('same')
     h['Vetotime'].Draw('same')
     tc = h['T'].cd(2)
     h['Stime'].Draw()
     rc = input("hit return for next event or q for quit: ")
     if rc=='q': break
 tc = h['T'].cd(1)
 h['SvsDStime'].Draw('colz')
 tc = h['T'].cd(2)
 h['SvsDStime_mufi'] = h['SvsDStime'].ProjectionY('SvsDStime_mufi')
 h['SvsDStime_scifi'] = h['SvsDStime'].ProjectionX('SvsDStime_scifi')
 h['Vetime'] = h['VEvsUStime'].ProjectionY('Vetime')
 h['UStime'] = h['VEvsUStime'].ProjectionX('UStime')
 h['SvsDStime_mufi'].SetLineColor(ROOT.kRed)
 h['SvsDStime_scifi'].SetLineColor(ROOT.kGreen)
 h['UStime'].SetLineColor(ROOT.kBlue)
 h['Vetime'].SetLineColor(ROOT.kOrange)
 h['UStime'].SetStats(0)
 h['Vetime'].SetStats(0)
 h['SvsDStime_mufi'].SetStats(0)
 h['SvsDStime_scifi'].SetStats(0)
 h['SvsDStime_mufi'].Draw()
 h['SvsDStime_scifi'].Draw('same')
 h['UStime'].Draw('same')
 h['Vetime'].Draw('same')

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

def USshower(Nev=options.nEvents, firstevent=0):
    zUS0 = zPos['MuFilter'][20] -  10
    zUS4 = zPos['MuFilter'][24] + 10
    for x in ['','-small']:
       ut.bookHist(h,'shower'+x,'energy vs z',200,0.,10000.,20,zUS0,zUS4)
       ut.bookHist(h,'showerX'+x,'energy vs z',200,0.,10000.,20,zUS0,zUS4)
       ut.bookHist(h,'wshower'+x,'z weighted energy ',100,-300.,0.)
       ut.bookHist(h,'zyshower'+x,'y vs z weighted energy ',20,zUS0,zUS4,11,-0.5,10.5)
    for p in range(systemAndPlanes[2]):
       ut.bookHist(h,'SvsL'+str(p),'small vs large Sipms plane' +str(p)+';large   [QCD]; small   [QCD] ',100,0.1,250.,100,0.1,100.)
    if Nev < 0 : Nev = eventTree.GetEntries()
    N=0
    ###################
    OUTPUT = open("qdc_dist_" + str(options.runNumber), 'w')
    SUM_SIGNAL = []
    ###################
    for event in eventTree:
       SUM_SIGNAL.append(0)
       N+=1
       if (N) % 1000 == 0:
          print(N)
       if N < firstevent:
          continue
       if N>Nev: break
       UShits = {}
       UShitsBar = {}
       for aHit in eventTree.Digi_MuFilterHits:
           if not aHit.isValid(): continue
           detID = aHit.GetDetectorID()
           s = aHit.GetDetectorID()//10000
           if s!=2: continue
           p = (aHit.GetDetectorID()//1000)%10
           S = map2Dict(aHit,'SumOfSignals')
           rc = h['SvsL'+str(p)].Fill(S['SumL'],S['SumS'])
           plane = (aHit.GetDetectorID()//1000)%10
           bar = aHit.GetDetectorID()%100
           if not plane in UShits: 
               UShits[plane]=0
               UShitsBar[plane]={}
               UShits[100+plane]=0
               UShitsBar[100+plane]={}
           if not bar in UShitsBar[plane]: 
                 UShitsBar[plane][bar]=0
                 UShitsBar[100+plane][bar]=0
           UShits[plane]+=S['Sum']
           UShitsBar[plane][bar]+=S['Sum']
           UShits[100+plane]+=S['SumS']
           UShitsBar[100+plane][bar]+=S['SumS']
       s = 2
       for plane in UShits:
           z = zPos['MuFilter'][s*10+plane%100]
           x = ''
           if plane > 99: x='-small'
           rc = h ['shower'+x].Fill(UShits[plane],z)
           if 0 in UShits:
               if UShits[0]>750: rc = h['showerX'+x].Fill(UShits[plane],z)
           rc = h ['wshower'+x].Fill(z,UShits[plane])
           for bar in UShitsBar[plane]:
                rc = h ['zyshower'+x].Fill(z,bar,UShitsBar[plane][bar])
                SUM_SIGNAL[N-1] += UShitsBar[plane][bar]
    for signal in SUM_SIGNAL:
      OUTPUT.write(str(signal) + "\n")
    OUTPUT.close()
    ut.bookCanvas(h,'lego','',900,1600,1,2)
    energy = {46:180,49:180,56:140,58:140,72:240,73:240,74:240,89:300,90:300,91:300}
    gain       = {46:2.5,49:3.65,52:3.65,54:2.5,56:2.5,58:3.65,72:1.0,73:2.5,74:3.65,86:3.65,87:2.5,88:1.0,89:3.65,90:2.5,91:1.0}
    tc = h['lego'].cd(1)
    tc.SetPhi(-20.5)
    tc.SetTheta(21.1)

    text = ""
    if options.runNumber in energy: text +="E=%5.1F"%(energy[options.runNumber])
    if options.runNumber in gain: text +="  with gain=%5.2F"%(gain[options.runNumber])
    h ['zyshower'].SetTitle(text+';z [cm]; y [bar number];QDC')
    h ['zyshower'].SetStats(0)
    h ['zyshower'].Draw('lego2')
    tc = h['lego'].cd(2)
    tc.SetPhi(-20.5)
    tc.SetTheta(21.1)
    h ['zyshower-small'].SetTitle('small sipms;z [cm]; y [bar number];QDC')
    h ['zyshower-small'].SetStats(0)
    h ['zyshower-small'].Draw('lego2')
    myPrint(h['lego'],'shower',withRootFile=True)
    
    ut.bookCanvas(h,'CorLS','',900,1200,1,5)
    h['SvsL'] =  h['SvsL0'].Clone('SvsL')
    h['SvsL'].SetTitle('small vs large Sipms all planes')
    for p in range(1,systemAndPlanes[2]):
        h['SvsL'].Add(h['SvsL'+str(p)])
    h['SvsL'].SetStats(0)
    for p in range(systemAndPlanes[2]):
        tc = h['CorLS'].cd(p+1)
        h['SvsL'+str(p)].SetStats(0)
        h['SvsL'+str(p)].Draw('colz')
    myPrint(h['CorLS'],'QCDsmallCSQCDlarge')

if options.command:
    tmp = options.command.split(';')

    if tmp[0].find('.C')>0:    # execute a C macro
        ROOT.gROOT.LoadMacro(tmp[0])
        ROOT.Execute()
        exit()
    command = tmp[0]+"("
    for i in range(1,len(tmp)):
         command+=tmp[i]
         if i<len(tmp)-1: command+=","
    command+=")"
    print('executing ' + command + "for run ",options.runNumber)
    eval(command)
    print('finished ' + command + "for run ",options.runNumber)
else:
    print ('waiting for command. Enter help() for more infomation')

def custom_ana_1(Nev=options.nEvents, firstevent=0):
   for i in range(5):
      for j in range(10):
         ut.bookHist(h,'plane_'+ str(i) + "_bar_" + str(j),'Num of fired sipm distribution, plane ' + str(i+1) + ", bar " + str(j+1), 16,0,16)
   ut.bookHist(h,'fired_bar','Fired bar distribution', 51,0,50)
   ut.bookHist(h,'fired_bar_1','Fired bar distribution', 51,0,50)
    
   if Nev < 0 : Nev = eventTree.GetEntries() 
   N=0 
   channel_info = {i:{j:[] for j in range(10)} for i in range(5)}

   for event in eventTree: 
      N+=1
      if (N) % 100000 == 0:
         print(N)
      if N < firstevent:
         continue
      if N>Nev: break 
      #if (N-1000000) % 100000 == 0:
      #   print(N-1000000)
      UShits = {} 
      UShitsBar = {} 
      M = 1
      listOfHits = []
      bar_num = 0
      bar_dist = {i:{} for i in range(5)}
      for aHit in eventTree.Digi_MuFilterHits: 
         if not aHit.isValid(): continue 
         detID = aHit.GetDetectorID() 
         s = aHit.GetDetectorID()//10000 
         if s!=2: continue 
         p = (aHit.GetDetectorID()//1000)%10 
         allChannels = map2Dict(aHit,'GetAllSignals') 
         #print("event", N, "hit", M, allChannels)
         M += 1
         plane = (aHit.GetDetectorID()//1000)%10 
         bar = aHit.GetDetectorID()%100 
         chan = 0 
         for c in allChannels: 
            channel_info[plane][bar].append(allChannels[c])
            if allChannels[c] != 0:
               chan += 1
         #print(N, plane, bar)
         h['plane_'+ str(plane) + "_bar_" + str(bar)].Fill(len(allChannels.keys()))
         #if bar not in bar_dist[plane] and chan >= 10:
         if bar not in bar_dist[plane]:
             bar_dist[plane][bar] = 0
             bar_num += 1
      h['fired_bar'].Fill(bar_num)
      if bar_num > 0:
         h['fired_bar_1'].Fill(bar_num)

   for i in range(5):
      ut.bookCanvas(h,'Fired_sipms_' + str(i),'',2000,1000,4,3)
   for i in range(5):
      for j in range(10):
        tc = h['Fired_sipms_' + str(i)].cd(j+1)
        #h['SvsL'+str(p)].SetStats(0)
        h['plane_'+ str(i) + "_bar_" + str(j)].Draw()

   for i in range(5):
      myPrint(h['Fired_sipms_' + str(i)],'Fired_sipms_' + str(i))


   ut.bookCanvas(h,'Fired_bar','',1000,800)
   h['fired_bar'].Draw()
   myPrint(h['Fired_bar'],'Fired_bars')
   ut.bookCanvas(h,'Fired_bar_morethan0','',1000,800)
   h['fired_bar_1'].Draw()
   myPrint(h['Fired_bar_morethan0'],'Fired_bar_morethan0')
   #print(channel_info)


def isShort(i):
   if i%8==3 or i%8==6: 
      return ROOT.kTRUE
   else: 
      return ROOT.kFALSE

def USshower_custom(Nev=options.nEvents):
    zUS0 = zPos['MuFilter'][20] -  10
    zUS4 = zPos['MuFilter'][24] + 10
    for x in ['','-small']:
       ut.bookHist(h,'shower'+x,'energy vs z',200,0.,10000.,20,zUS0,zUS4)
       ut.bookHist(h,'showerX'+x,'energy vs z',200,0.,10000.,20,zUS0,zUS4)
       ut.bookHist(h,'wshower'+x,'z weighted energy ',100,-300.,0.)
       ut.bookHist(h,'zyshower'+x,'y vs z weighted energy ',20,zUS0,zUS4,11,-0.5,10.5)
    for p in range(systemAndPlanes[2]):
       ut.bookHist(h,'SvsL'+str(p),'small vs large Sipms plane' +str(p)+';large   [QCD]; small   [QCD] ',100,0.1,250.,100,0.1,100.)
    if Nev < 0 : Nev = eventTree.GetEntries()
    N=0
    for event in eventTree:
       N+=1
       if (N) % 100000 == 0:
          print(N)
       #if N < 1000000:
       #   continue
       if N>Nev: break
       UShits = {}
       UShitsBar = {}
       for aHit in eventTree.Digi_MuFilterHits:
           if not aHit.isValid(): continue
           detID = aHit.GetDetectorID()
           s = aHit.GetDetectorID()//10000
           if s!=2: continue
           p = (aHit.GetDetectorID()//1000)%10
           S = map2Dict(aHit,'SumOfSignals')
           allChannels = map2Dict(aHit,'GetAllSignals')
           rc = h['SvsL'+str(p)].Fill(S['SumL'],S['SumS'])
           plane = (aHit.GetDetectorID()//1000)%10
           bar = aHit.GetDetectorID()%100
           if not plane in UShits: 
               UShits[plane]=0
               UShitsBar[plane]={}
           if not bar in UShitsBar[plane]: 
                 UShitsBar[plane][bar]=0

           chan = 0 
           channel_info = {i:{j:[] for j in range(10)} for i in range(5)}
           for c in allChannels:              
                 if allChannels[c] != 0 and not isShort(c):
                     channel_info[plane][bar].append(allChannels[c])
                     chan += 1
           if chan > 1:
                 #UShits[plane]+= sum(channel_info[plane][bar])
                 UShitsBar[plane][bar]+= sum(channel_info[plane][bar])
       s = 2
       for plane in UShits:
           z = zPos['MuFilter'][s*10+plane%100]
           x = ''
           if plane > 99: x='-small'
           #rc = h ['shower'+x].Fill(UShits[plane],z)
           if 0 in UShits:
               if UShits[0]>750: rc = h['showerX'+x].Fill(UShits[plane],z)
           #rc = h ['wshower'+x].Fill(z,UShits[plane])
           for bar in UShitsBar[plane]:
                rc = h ['zyshower'+x].Fill(z,bar,UShitsBar[plane][bar])
    ut.bookCanvas(h,'lego','',900,1600,1,2)
    energy = {46:180,49:180,56:140,58:140,72:240,73:240,74:240,89:300,90:300,91:300}
    gain       = {46:2.5,49:3.65,52:3.65,54:2.5,56:2.5,58:3.65,72:1.0,73:2.5,74:3.65,86:3.65,87:2.5,88:1.0,89:3.65,90:2.5,91:1.0}
    tc = h['lego'].cd(1)
    tc.SetPhi(-20.5)
    tc.SetTheta(21.1)

    text = ""
    if options.runNumber in energy: text +="E=%5.1F"%(energy[options.runNumber])
    if options.runNumber in gain: text +="  with gain=%5.2F"%(gain[options.runNumber])
    h ['zyshower'].SetTitle(text+';z [cm]; y [bar number];QDC')
    h ['zyshower'].SetStats(0)
    h ['zyshower'].Draw('lego2')
    tc = h['lego'].cd(2)
    tc.SetPhi(-20.5)
    tc.SetTheta(21.1)
    h ['zyshower-small'].SetTitle('small sipms;z [cm]; y [bar number];QDC')
    h ['zyshower-small'].SetStats(0)
    h ['zyshower-small'].Draw('lego2')
    myPrint(h['lego'],'shower_filtered',withRootFile=True)
    
    ut.bookCanvas(h,'CorLS','',900,1200,1,5)
    h['SvsL'] =  h['SvsL0'].Clone('SvsL')
    h['SvsL'].SetTitle('small vs large Sipms all planes')
    for p in range(1,systemAndPlanes[2]):
        h['SvsL'].Add(h['SvsL'+str(p)])
    h['SvsL'].SetStats(0)
    for p in range(systemAndPlanes[2]):
        tc = h['CorLS'].cd(p+1)
        h['SvsL'+str(p)].SetStats(0)
        h['SvsL'+str(p)].Draw('colz')
    myPrint(h['CorLS'],'QCDsmallCSQCDlarge')


def eventTime_for_filter(Nev=options.nEvents):
 if Nev < 0 : Nev = eventTree.GetEntries()
 ut.bookHist(h,'Etime','delta event time; dt [s]',100,0.0,1.)
 ut.bookHist(h,'EtimeZ','delta event time; dt [ns]',1000,0.0,10000.)
 ut.bookCanvas(h,'T',' ',1024,3*768,1,3)
 
# need to make extra gymnastiques since absolute time is missing
 Ntinter = []
 N = 0
 for f in eventTree.GetListOfFiles():
    dN =  f.GetEntries()
    rc = eventTree.GetEvent(N)
    t0 = eventTree.EventHeader.GetEventTime()/freq
    rc = eventTree.GetEvent(N+dN-1)
    tmax = eventTree.EventHeader.GetEventTime()/freq
    Ntinter.append([t0,tmax])
    N+=dN

 Tduration = 0
 for x in Ntinter:
    Tduration += (x[1]-x[0])
 tsep = 3600.
 t0 =  Ntinter[0][0]
 tmax = Tduration+(tsep*(len(Ntinter)-1)) 
 #nbins = int(tmax-t0)
 nbins = 1000
 yunit = "events per %5.0F s"%( (tmax-t0)/nbins)
 if 'time' in h: h.pop('time').Delete()
 #ut.bookHist(h,'time','elapsed time; t [s];'+yunit,nbins,0,tmax-t0)
 ut.bookHist(h,'time','elapsed time; t [s];'+yunit,nbins,0,tmax-t0)
 N=-1
 Tprev  = 0
 Toffset = 0
 for event in eventTree:
    N+=1
    if N>Nev: break
    T   = event.EventHeader.GetEventTime()
    dT = T-Tprev
    if N>0 and dT >0:
           rc = h['Etime'].Fill( dT/freq )
           rc = h['EtimeZ'].Fill( dT*1E9/freq )
           rc = h['time'].Fill( (T+Toffset)/freq-t0 )
    elif dT<0: 
           Toffset+=tsep*freq+Tprev
           rc = h['time'].Fill( (T+Toffset)/freq-t0 )
    else: rc = h['time'].Fill( T/freq-t0 ) # very first event
    Tprev = T

 tc = h['T'].cd(1)
 
 bin_max = 0
 for i in range(100):
   if h['time'].GetBinContent(i) > bin_max:
      bin_max = h['time'].GetBinContent(i)
 time_selected = []
 for i in range(nbins):
   event_num = h['time'].GetBinContent(i)
   if event_num >= 0.6*bin_max and event_num < 1.5*bin_max:
      time_selected.append(h['time'].GetBinCenter(i))
   else:
      h['time'].SetBinContent(i, 0)


 
 h['time'].SetStats(0)
 h['time'].Draw()
 


 tend = 0
 for x in Ntinter:
    tend += x[1]+tsep/2.
    m = str(int(tend))
    h['line'+m]=ROOT.TLine(tend,0,tend,h['time'].GetMaximum())
    h['line'+m].SetLineColor(ROOT.kRed)
    h['line'+m].Draw()
    tend += tsep/2.
 tc = h['T'].cd(2)
 tc.SetLogy(1)
 h['EtimeZ'].Draw()
 rc = h['EtimeZ'].Fit('expo','S','',0.,250.)
 h['T'].Update()
 stats = h['EtimeZ'].FindObject('stats')
 stats.SetOptFit(1111111)
 tc = h['T'].cd(3)
 tc.SetLogy(1)
 h['Etime'].Draw()
 rc = h['Etime'].Fit('expo','S')
 h['T'].Update()
 stats = h['Etime'].FindObject('stats')
 stats.SetOptFit(1111111)
 h['T'].Update()
 myPrint(h['T'],'time')
 return time_selected, h['time'].GetXaxis().GetBinWidth(0), t0



def USshower_for_filter(time_selected, binwidth, t0, Nev=options.nEvents, firstevent=0):
    zUS0 = zPos['MuFilter'][20] -  10
    zUS4 = zPos['MuFilter'][24] + 10
    for x in ['','-small']:
       ut.bookHist(h,'shower'+x,'energy vs z',200,0.,10000.,20,zUS0,zUS4)
       ut.bookHist(h,'showerX'+x,'energy vs z',200,0.,10000.,20,zUS0,zUS4)
       ut.bookHist(h,'wshower'+x,'z weighted energy ',100,-300.,0.)
       ut.bookHist(h,'zyshower'+x,'y vs z weighted energy ',20,zUS0,zUS4,11,-0.5,10.5)
    for p in range(systemAndPlanes[2]):
       ut.bookHist(h,'SvsL'+str(p),'small vs large Sipms plane' +str(p)+';large   [QCD]; small   [QCD] ',100,0.1,250.,100,0.1,100.)
    if Nev < 0 : Nev = eventTree.GetEntries()
    N=0
    SUM_SIGNAL = []
    N_filtered = 0
    ###########################
    Tprev  = 0
    Toffset = 0
    ###########################
    print("Before loop...............")
    for event in eventTree:
       SUM_SIGNAL.append(0)
       N+=1
       if (N) % 100000 == 0:
          print(N)
       #print("N = ", N, "Nev = ", Nev, "firstevent = ", firstevent)
       if N < firstevent:
          continue
       if N>Nev: break
       #print("In loop...............")
      #########################
       T   = event.EventHeader.GetEventTime()
       dT = T-Tprev
       if dT >0:
           TIME = (T+Toffset)/freq-t0
       elif dT<0: 
           Toffset+=tsep*freq+Tprev
           TIME = (T+Toffset)/freq-t0
       else:
           TIME = T/freq-t0
       Tprev = T
       check_key = 0
       #print(time_selected, binwidth)
       for bincenter in time_selected:
           if TIME >= bincenter - binwidth/2. and TIME < bincenter + binwidth/2:
               check_key = 1
               break 
       if check_key == 0:
           continue
       N_filtered += 1
      #########################
       UShits = {}
       UShitsBar = {}
       for aHit in eventTree.Digi_MuFilterHits:
           if not aHit.isValid(): continue
           detID = aHit.GetDetectorID()
           s = aHit.GetDetectorID()//10000
           if s!=2: continue
           p = (aHit.GetDetectorID()//1000)%10
           S = map2Dict(aHit,'SumOfSignals')
           rc = h['SvsL'+str(p)].Fill(S['SumL'],S['SumS'])
           plane = (aHit.GetDetectorID()//1000)%10
           bar = aHit.GetDetectorID()%100
           if not plane in UShits: 
               UShits[plane]=0
               UShitsBar[plane]={}
               UShits[100+plane]=0
               UShitsBar[100+plane]={}
           if not bar in UShitsBar[plane]: 
                 UShitsBar[plane][bar]=0
                 UShitsBar[100+plane][bar]=0
           UShits[plane]+=S['Sum']
           UShitsBar[plane][bar]+=S['Sum']
           UShits[100+plane]+=S['SumS']
           UShitsBar[100+plane][bar]+=S['SumS']
       s = 2
       for plane in UShits:
           z = zPos['MuFilter'][s*10+plane%100]
           x = ''
           if plane > 99: x='-small'
           rc = h ['shower'+x].Fill(UShits[plane],z)
           if 0 in UShits:
               if UShits[0]>750: rc = h['showerX'+x].Fill(UShits[plane],z)
           rc = h ['wshower'+x].Fill(z,UShits[plane])
           for bar in UShitsBar[plane]:
                rc = h ['zyshower'+x].Fill(z,bar,UShitsBar[plane][bar])
                SUM_SIGNAL[N-1] += UShitsBar[plane][bar]
    print("QDC =",np.array(SUM_SIGNAL).mean())
    print("Number of filtered events: ", N_filtered)
    ut.bookCanvas(h,'lego','',900,1600,1,2)
    energy = {46:180,49:180,56:140,58:140,72:240,73:240,74:240,89:300,90:300,91:300}
    gain       = {46:2.5,49:3.65,52:3.65,54:2.5,56:2.5,58:3.65,72:1.0,73:2.5,74:3.65,86:3.65,87:2.5,88:1.0,89:3.65,90:2.5,91:1.0}
    tc = h['lego'].cd(1)
    tc.SetPhi(-20.5)
    tc.SetTheta(21.1)

    text = ""
    if options.runNumber in energy: text +="E=%5.1F"%(energy[options.runNumber])
    if options.runNumber in gain: text +="  with gain=%5.2F"%(gain[options.runNumber])
    h ['zyshower'].SetTitle(text+';z [cm]; y [bar number];QDC')
    h ['zyshower'].SetStats(0)
    h ['zyshower'].Draw('lego2')
    tc = h['lego'].cd(2)
    tc.SetPhi(-20.5)
    tc.SetTheta(21.1)
    h ['zyshower-small'].SetTitle('small sipms;z [cm]; y [bar number];QDC')
    h ['zyshower-small'].SetStats(0)
    h ['zyshower-small'].Draw('lego2')
    myPrint(h['lego'],'shower',withRootFile=True)
    
    ut.bookCanvas(h,'CorLS','',900,1200,1,5)
    h['SvsL'] =  h['SvsL0'].Clone('SvsL')
    h['SvsL'].SetTitle('small vs large Sipms all planes')
    for p in range(1,systemAndPlanes[2]):
        h['SvsL'].Add(h['SvsL'+str(p)])
    h['SvsL'].SetStats(0)
    for p in range(systemAndPlanes[2]):
        tc = h['CorLS'].cd(p+1)
        h['SvsL'+str(p)].SetStats(0)
        h['SvsL'+str(p)].Draw('colz')
    myPrint(h['CorLS'],'QCDsmallCSQCDlarge')
    
if options.command:
    tmp = options.command.split(';')

    if tmp[0].find('.C')>0:    # execute a C macro
        ROOT.gROOT.LoadMacro(tmp[0])
        ROOT.Execute()
        exit()
    command = tmp[0]+"("
    for i in range(1,len(tmp)):
         command+=tmp[i]
         if i<len(tmp)-1: command+=","
    command+=")"
    print('executing ' + command + "for run ",options.runNumber)
    eval(command)
    print('finished ' + command + "for run ",options.runNumber)
else:
    print ('waiting for command. Enter help() for more infomation')


def execute_usshower_filtered(Nev=options.nEvents):
   time_selected, binwidth, t0 = eventTime_for_filter(Nev)
   USshower_for_filter(time_selected, binwidth, t0, Nev)



def Mufi_hitMaps_qdc(Nev = options.nEvents):
 # veto system 2 layers with 7 bars and 8 sipm channels on both ends
 # US system 5 layers with 10 bars and 8 sipm channels on both ends
 # DS system horizontal(3) planes, 60 bars, readout on both sides, single channel
 #                         vertical(4) planes, 60 bar, readout on top, single channel
 for s in systemAndPlanes:
    for l in range(systemAndPlanes[s]):
       ut.bookHist(h,'sig_'+str(s*10+l),'signal / plane '+str(s*10+l),200,0.0,50.)
 
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
        nSiPMs = aHit.GetnSiPMs()
        nSides  = aHit.GetnSides()
# check left/right
        allChannels = map2Dict(aHit,'GetAllSignals')
#
        for c in allChannels:
            channel = bar*nSiPMs*nSides + c
            rc  = h['sig_'+str(s)+str(l)].Fill(allChannels[c])
        allChannels.clear()
#
    
 S = {1:[1800,800,2,1],2:[1800,1500,2,3],3:[1800,1800,2,4]}
 for s in S:
   ut.bookCanvas(h,'signal'    +str(s),' ',S[s][0],S[s][1],S[s][2],S[s][3])
   h['signal'    +str(s)].SetLogy()
   for l in range(systemAndPlanes[s]):
      n = l+1
      if s==3 and n==7: n=8
      tag = str(s)+str(l)
      tc = h['signal'+str(s)].cd(n)
      h['sig_'+tag].Draw()

 for s in S:
   h['signal'    +str(s)].Update()
   myPrint(h['signal'    +str(s)],'signal_extr_'+sdict[s])



def Mufi_hitMaps_qdc_filtered(time_selected, binwidth, t0, Nev = options.nEvents):
 # veto system 2 layers with 7 bars and 8 sipm channels on both ends
 # US system 5 layers with 10 bars and 8 sipm channels on both ends
 # DS system horizontal(3) planes, 60 bars, readout on both sides, single channel
 #                         vertical(4) planes, 60 bar, readout on top, single channel
 for s in systemAndPlanes:
    for l in range(systemAndPlanes[s]):
       ut.bookHist(h,'sig_'+str(s*10+l),'signal / plane '+str(s*10+l),200,0.0,50.)
 
 N=-1
 N_filtered = 0
 if Nev < 0 : Nev = eventTree.GetEntries()
 ###########################
 Tprev  = 0
 Toffset = 0
 ###########################
 print("Before loop...............")
 eventTree.GetEvent(0)
 #t0 =  eventTree.EventHeader.GetEventTime()/freq
 listOfHits = {1:[],2:[],3:[]}
 mult = {}
 for event in eventTree:
    N+=1
    if N%options.heartBeat == 0: print('event ',N,' ',time.ctime())
    if N>Nev: break
#########################
    T   = event.EventHeader.GetEventTime()
    dT = T-Tprev
    if dT >0:
      TIME = (T+Toffset)/freq-t0
    elif dT<0: 
      Toffset+=tsep*freq+Tprev
      TIME = (T+Toffset)/freq-t0
    else:
      TIME = T/freq-t0
    Tprev = T
    check_key = 0
   #print(time_selected, binwidth)
    for bincenter in time_selected:
      if TIME >= bincenter - binwidth/2. and TIME < bincenter + binwidth/2:
         check_key = 1
         break 
    if check_key == 0:
      continue
    N_filtered += 1
#########################
    

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
        nSiPMs = aHit.GetnSiPMs()
        nSides  = aHit.GetnSides()
# check left/right
        allChannels = map2Dict(aHit,'GetAllSignals')
#
        for c in allChannels:
            channel = bar*nSiPMs*nSides + c
            rc  = h['sig_'+str(s)+str(l)].Fill(allChannels[c])
        allChannels.clear()
#
 print("Number of filtered events: ", N_filtered)   
 S = {1:[1800,800,2,1],2:[1800,1500,2,3],3:[1800,1800,2,4]}
 for s in S:
   ut.bookCanvas(h,'signal'    +str(s),' ',S[s][0],S[s][1],S[s][2],S[s][3])

   for l in range(systemAndPlanes[s]):
      n = l+1
      if s==3 and n==7: n=8
      tag = str(s)+str(l)
      tc = h['signal'+str(s)].cd(n)
      h['sig_'+tag].Draw()
      
 for s in S:
   h['signal'    +str(s)].Update()
   myPrint(h['signal'    +str(s)],'signal'+sdict[s])

def execute_mufi_hit_maps_filtered(Nev=options.nEvents):
   time_selected, binwidth, t0 = eventTime_for_filter(Nev)
   Mufi_hitMaps_qdc_filtered(time_selected, binwidth, t0, Nev)


def eventTime_for_filter_trash(Nev=options.nEvents):
 if Nev < 0 : Nev = eventTree.GetEntries()
 ut.bookHist(h,'Etime','delta event time; dt [s]',100,0.0,1.)
 ut.bookHist(h,'EtimeZ','delta event time; dt [ns]',1000,0.0,10000.)
 ut.bookCanvas(h,'T',' ',1024,3*768,1,3)
 
# need to make extra gymnastiques since absolute time is missing
 Ntinter = []
 N = 0
 for f in eventTree.GetListOfFiles():
    dN =  f.GetEntries()
    rc = eventTree.GetEvent(N)
    t0 = eventTree.EventHeader.GetEventTime()/freq
    rc = eventTree.GetEvent(N+dN-1)
    tmax = eventTree.EventHeader.GetEventTime()/freq
    Ntinter.append([t0,tmax])
    N+=dN

 Tduration = 0
 for x in Ntinter:
    Tduration += (x[1]-x[0])
 tsep = 3600.
 t0 =  Ntinter[0][0]
 tmax = Tduration+(tsep*(len(Ntinter)-1)) 
 #nbins = int(tmax-t0)
 nbins = 500
 yunit = "events per %5.0F s"%( (tmax-t0)/nbins)
 if 'time' in h: h.pop('time').Delete()
 #ut.bookHist(h,'time','elapsed time; t [s];'+yunit,nbins,0,tmax-t0)
 ut.bookHist(h,'time','elapsed time; t [s];'+yunit,nbins,0,tmax-t0 - 3000)
 N=-1
 Tprev  = 0
 Toffset = 0
 for event in eventTree:
    N+=1
    if N>Nev: break
    T   = event.EventHeader.GetEventTime()
    dT = T-Tprev
    if N>0 and dT >0:
           rc = h['Etime'].Fill( dT/freq )
           rc = h['EtimeZ'].Fill( dT*1E9/freq )
           rc = h['time'].Fill( (T+Toffset)/freq-t0 )
    elif dT<0: 
           Toffset+=tsep*freq+Tprev
           rc = h['time'].Fill( (T+Toffset)/freq-t0 )
    else: rc = h['time'].Fill( T/freq-t0 ) # very first event
    Tprev = T

 tc = h['T'].cd(1)
 
 bin_max = 0
 for i in range(100):
   if h['time'].GetBinContent(i) > bin_max:
      bin_max = h['time'].GetBinContent(i)
 time_selected = []
 for i in range(nbins):
   event_num = h['time'].GetBinContent(i)
   if event_num < 0.4*bin_max:
      time_selected.append(h['time'].GetBinCenter(i))
   else:
      h['time'].SetBinContent(i, 0)


 
 h['time'].SetStats(0)
 h['time'].Draw()
 


 tend = 0
 for x in Ntinter:
    tend += x[1]+tsep/2.
    m = str(int(tend))
    h['line'+m]=ROOT.TLine(tend,0,tend,h['time'].GetMaximum())
    h['line'+m].SetLineColor(ROOT.kRed)
    h['line'+m].Draw()
    tend += tsep/2.
 tc = h['T'].cd(2)
 tc.SetLogy(1)
 h['EtimeZ'].Draw()
 rc = h['EtimeZ'].Fit('expo','S','',0.,250.)
 h['T'].Update()
 stats = h['EtimeZ'].FindObject('stats')
 stats.SetOptFit(1111111)
 tc = h['T'].cd(3)
 tc.SetLogy(1)
 h['Etime'].Draw()
 rc = h['Etime'].Fit('expo','S')
 h['T'].Update()
 stats = h['Etime'].FindObject('stats')
 stats.SetOptFit(1111111)
 h['T'].Update()
 myPrint(h['T'],'time_trash')
 return time_selected, h['time'].GetXaxis().GetBinWidth(0), t0


def Mufi_hitMaps_qdc_sipm_filter(Nev = options.nEvents, sipm_num = 10):
 # veto system 2 layers with 7 bars and 8 sipm channels on both ends
 # US system 5 layers with 10 bars and 8 sipm channels on both ends
 # DS system horizontal(3) planes, 60 bars, readout on both sides, single channel
 #                         vertical(4) planes, 60 bar, readout on top, single channel
 for s in systemAndPlanes:
    for l in range(systemAndPlanes[s]):
       ut.bookHist(h,'sig_'+str(s*10+l),'signal / plane '+str(s*10+l),200,0.0,50.)
 
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
        nSiPMs = aHit.GetnSiPMs()
        nSides  = aHit.GetnSides()
# check left/right
        allChannels = map2Dict(aHit,'GetAllSignals')
#
        ch = 0 
        for c in allChannels:
            if allChannels[c] != 0:
               ch += 1
        if ch >= sipm_num:
            for c in allChannels:
               channel = bar*nSiPMs*nSides + c
               rc  = h['sig_'+str(s)+str(l)].Fill(allChannels[c])
        allChannels.clear()
#
    
 S = {1:[1800,800,2,1],2:[1800,1500,2,3],3:[1800,1800,2,4]}
 for s in S:
   ut.bookCanvas(h,'signal'    +str(s),' ',S[s][0],S[s][1],S[s][2],S[s][3])
   for l in range(systemAndPlanes[s]):
      n = l+1
      if s==3 and n==7: n=8
      tag = str(s)+str(l)
      tc = h['signal'+str(s)].cd(n)
      h['sig_'+tag].Draw()

 for s in S:
   h['signal'    +str(s)].Update()
   myPrint(h['signal'    +str(s)],'signal_sipm_filtered'+sdict[s])


def Mufi_hitMaps_qdc_5planes_hits(Nev = options.nEvents):
 # veto system 2 layers with 7 bars and 8 sipm channels on both ends
 # US system 5 layers with 10 bars and 8 sipm channels on both ends
 # DS system horizontal(3) planes, 60 bars, readout on both sides, single channel
 #                         vertical(4) planes, 60 bar, readout on top, single channel
 for s in systemAndPlanes:
    for l in range(systemAndPlanes[s]):
       ut.bookHist(h,'sig_'+str(s*10+l),'signal / plane '+str(s*10+l),200,0.0,50.)
 
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

    ##########
    list_of_fired_planes = []
    for aHit in event.Digi_MuFilterHits:
        if not aHit.isValid(): continue
        detID = aHit.GetDetectorID()
        if aHit.isVertical():     withX = True
        s = detID//10000
        l  = (detID%10000)//1000  # plane number
        list_of_fired_planes.append(str(s) + str(l))
    if list_of_fired_planes != [str(i) for i in range(20,25)]:
      continue

   ###########
    for aHit in event.Digi_MuFilterHits:
        if not aHit.isValid(): continue
        detID = aHit.GetDetectorID()
        if aHit.isVertical():     withX = True
        s = detID//10000
        l  = (detID%10000)//1000  # plane number
        list_of_fired_planes.append(str(s) + str(l))
        bar = (detID%1000)
        if s>2:
               l=2*l
               if bar>59:
                    bar=bar-60
                    if l<6: l+=1
        mult[s*10+l]+=1
        key = s*100+l
        nSiPMs = aHit.GetnSiPMs()
        nSides  = aHit.GetnSides()
# check left/right
        allChannels = map2Dict(aHit,'GetAllSignals')
        for c in allChannels:
         channel = bar*nSiPMs*nSides + c
         rc  = h['sig_'+str(s)+str(l)].Fill(allChannels[c])
        allChannels.clear()

#
    
 S = {1:[1800,800,2,1],2:[1800,1500,2,3],3:[1800,1800,2,4]}
 for s in S:
   ut.bookCanvas(h,'signal'    +str(s),' ',S[s][0],S[s][1],S[s][2],S[s][3])
   for l in range(systemAndPlanes[s]):
      n = l+1
      if s==3 and n==7: n=8
      tag = str(s)+str(l)
      tc = h['signal'+str(s)].cd(n)
      h['sig_'+tag].Draw()

 for s in S:
   h['signal'    +str(s)].Update()
   myPrint(h['signal'    +str(s)],'signal_5planes_filter_'+sdict[s])



def Mufi_hitMaps_mean_qdc(Nev = options.nEvents, key_filter = 0):
 # veto system 2 layers with 7 bars and 8 sipm channels on both ends
 # US system 5 layers with 10 bars and 8 sipm channels on both ends
 # DS system horizontal(3) planes, 60 bars, readout on both sides, single channel
 #                         vertical(4) planes, 60 bar, readout on top, single channel
 for s in systemAndPlanes:
    for l in range(systemAndPlanes[s]):
       ut.bookHist(h,'sig_'+str(s*10+l),'signal / plane '+str(s*10+l),200,0.0,50.)
       ut.bookHist(h,'left_sig_'+str(s*10+l),'left_signal / plane '+str(s*10+l),200,0.0,50.)
       ut.bookHist(h,'right_sig_'+str(s*10+l),'right_signal / plane '+str(s*10+l),200,0.0,50.)
 
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
        nSiPMs = aHit.GetnSiPMs()
        nSides  = aHit.GetnSides()
# check left/right
        allChannels = map2Dict(aHit,'GetAllSignals')

        # check the last US plane:
        #if s == 2 and l == 4:
         



        if key_filter == 0:
         channels = [allChannels[c] for c in allChannels]
         rc  = h['sig_'+str(s)+str(l)].Fill(np.array(channels).mean())
         if nSides==2:
          Nleft    = 0
          Nright = 0
          Sleft    = []
          Sright = []
          for c in allChannels:
            if  nSiPMs > c:  # left side
                  Nleft+=1
                  Sleft.append(allChannels[c])
            else:
                  Nright+=1
                  Sright.append(allChannels[c])
         rc  = h['left_sig_'+str(s)+str(l)].Fill(np.array(Sleft).mean())
         rc  = h['right_sig_'+str(s)+str(l)].Fill(np.array(Sright).mean())
        else:
         ch = 0 
         for c in allChannels:
               if allChannels[c] != 0:
                  ch += 1
         if ch >= 10:
               channels = [allChannels[c] for c in allChannels]
               rc  = h['sig_'+str(s)+str(l)].Fill(np.array(channels).mean())
               rc  = h['left_sig_'+str(s)+str(l)].Fill(np.array(Sleft).mean())
               rc  = h['right_sig_'+str(s)+str(l)].Fill(np.array(Sright).mean())
        
        allChannels.clear()
#
    
 S = {1:[1800,800,2,1],2:[1800,1500,2,3],3:[1800,1800,2,4]}
 type_of_hist = ["", "left_", "right_"]
 for TYPE in type_of_hist:
   for s in S:
      ut.bookCanvas(h,TYPE + 'signal'    +str(s),' ',S[s][0],S[s][1],S[s][2],S[s][3])
      h[TYPE + 'signal'    +str(s)].SetLogy()
      for l in range(systemAndPlanes[s]):
         n = l+1
         if s==3 and n==7: n=8
         tag = str(s)+str(l)
         tc = h[TYPE + 'signal'+str(s)].cd(n)
         h[TYPE + 'sig_'+tag].Draw()

   for s in S:
      h[TYPE + 'signal'    +str(s)].Update()
      myPrint(h[TYPE + 'signal'    +str(s)],TYPE + 'signal_mean_key_' + str(key_filter) + "_" + sdict[s])




def MIP_study(Nev = options.nEvents, key_filter = 0):
 # veto system 2 layers with 7 bars and 8 sipm channels on both ends
 # US system 5 layers with 10 bars and 8 sipm channels on both ends
 # DS system horizontal(3) planes, 60 bars, readout on both sides, single channel
 #                         vertical(4) planes, 60 bar, readout on top, single channel
 for s in systemAndPlanes:
    for l in range(systemAndPlanes[s]):
       ut.bookHist(h,'sig_'+str(s*10+l),'signal / plane '+str(s*10+l),200,0.0,50.)
       ut.bookHist(h,'left_sig_'+str(s*10+l),'left_signal / plane '+str(s*10+l),200,0.0,50.)
       ut.bookHist(h,'right_sig_'+str(s*10+l),'right_signal / plane '+str(s*10+l),200,0.0,50.)
       ut.bookHist(h,'sqrt_sig_'+str(s*10+l),'sqrt_signal / plane '+str(s*10+l),200,0.0,50.)
 
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
    good_event_US = False
    k_good_event_US = 0
    bar_good_event_US_id = []
    good_event_DS = False
    k_good_event_DS = 0
    bar_good_event_DS_id = []
    hits_list = []
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
        nSiPMs = aHit.GetnSiPMs()
        nSides  = aHit.GetnSides()
# check left/right
        allChannels = map2Dict(aHit,'GetAllSignals')

        # check the last US plane:
        #if s == 2 and l == 4:
        

        if not good_event_US:
         if s != 2 or l != 4:
            continue
         else:
            k_good_event_US += 1
            bar_good_event_US_id.append(detID%1000)
            if k_good_event_US == 2:
               if bar_good_event_US_id[0] == bar_good_event_US_id[1]:
                  good_event_US = True
        

        if not good_event_DS:
         if s != 3 or l != 0 or aHit.isVertical():
            continue
         else:
            k_good_event_DS += 1
            bar_good_event_DS_id.append(detID%1000)
            if k_good_event_DS == 2:
               if bar_good_event_DS_id[0] != bar_good_event_DS_id[1]:
                  good_event_DS = True
    if good_event_US:
      print(N, good_event_US, good_event_DS, bar_good_event_US_id, bar_good_event_DS_id)   
    if good_event_US and good_event_DS:
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
         nSiPMs = aHit.GetnSiPMs()
         nSides  = aHit.GetnSides()
   # check left/right
         allChannels = map2Dict(aHit,'GetAllSignals')
         if key_filter == 0:
            channels = [allChannels[c] for c in allChannels]
            rc  = h['sig_'+str(s)+str(l)].Fill(np.array(channels).mean())
            if nSides==2:
               Sleft    = []
               Sright = []
               for c in allChannels:
                  if  nSiPMs > c:  # left side
                        Sleft.append(allChannels[c])
                  else:
                        Sright.append(allChannels[c])
               rc  = h['left_sig_'+str(s)+str(l)].Fill(np.array(Sleft).mean())
               rc  = h['right_sig_'+str(s)+str(l)].Fill(np.array(Sright).mean())
               rc  = h['sqrt_sig_'+str(s)+str(l)].Fill(np.sqrt(np.array(Sleft).mean()*np.array(Sright).mean()))

#
    
 S = {1:[1800,800,2,1],2:[1800,1500,2,3],3:[1800,1800,2,4]}
 type_of_hist = ["", "left_", "right_", "sqrt_"]
 for TYPE in type_of_hist:
   for s in S:
      ut.bookCanvas(h,TYPE + 'signal'    +str(s),' ',S[s][0],S[s][1],S[s][2],S[s][3])
      h[TYPE + 'signal'    +str(s)].SetLogy()
      for l in range(systemAndPlanes[s]):
         n = l+1
         if s==3 and n==7: n=8
         tag = str(s)+str(l)
         tc = h[TYPE + 'signal'+str(s)].cd(n)
         h[TYPE + 'sig_'+tag].Draw()

   for s in S:
      h[TYPE + 'signal'    +str(s)].Update()
      myPrint(h[TYPE + 'signal'    +str(s)],TYPE + 'mip_signal_mean_key_' + str(key_filter) + "_" + sdict[s])



def MIP_study_1(Nev = options.nEvents, key_filter = 0):
 # veto system 2 layers with 7 bars and 8 sipm channels on both ends
 # US system 5 layers with 10 bars and 8 sipm channels on both ends
 # DS system horizontal(3) planes, 60 bars, readout on both sides, single channel
 #                         vertical(4) planes, 60 bar, readout on top, single channel
 for s in systemAndPlanes:
    for l in range(systemAndPlanes[s]):
       ut.bookHist(h,'sig_'+str(s*10+l),'signal / plane '+str(s*10+l),200,0.0,50.)
       ut.bookHist(h,'left_sig_'+str(s*10+l),'left_signal / plane '+str(s*10+l),200,0.0,50.)
       ut.bookHist(h,'right_sig_'+str(s*10+l),'right_signal / plane '+str(s*10+l),200,0.0,50.)
       ut.bookHist(h,'sqrt_sig_'+str(s*10+l),'sqrt_signal / plane '+str(s*10+l),200,0.0,50.)
 
 N=-1
 Tprev = 0
 if Nev < 0 : Nev = eventTree.GetEntries()
 eventTree.GetEvent(0)
 t0 =  eventTree.EventHeader.GetEventTime()/freq
 listOfHits = {1:[],2:[],3:[]}
 mult = {}
 US_data = []
 DS_data = []
 for event in eventTree:
    N+=1
    if N%options.heartBeat == 0: print('event ',N,' ',time.ctime())
    if N>Nev: break
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
        nSiPMs = aHit.GetnSiPMs()
        nSides  = aHit.GetnSides()
# check left/right
        allChannels = map2Dict(aHit,'GetAllSignals')

        # check the last US plane:
        #if s == 2 and l == 4:
        

        if not good_event_US:
         if s != 2 or l != 4:
            continue
         else:
            k_good_event_US += 1
            bar_good_event_US_id.append(detID%1000)
            if k_good_event_US == 2:
               if bar_good_event_US_id[0] == bar_good_event_US_id[1]:
                  good_event_US = True
        

        if not good_event_DS:
         if s != 3 or l != 0 or aHit.isVertical():
            continue
         else:
            k_good_event_DS += 1
            bar_good_event_DS_id.append(detID%1000)
            if k_good_event_DS == 2:
               if bar_good_event_DS_id[0] != bar_good_event_DS_id[1]:
                  good_event_DS = True
    if good_event_US:
      print(N, good_event_US, good_event_DS, bar_good_event_US_id, bar_good_event_DS_id)   
    if good_event_US and good_event_DS:
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
         nSiPMs = aHit.GetnSiPMs()
         nSides  = aHit.GetnSides()
   # check left/right
         allChannels = map2Dict(aHit,'GetAllSignals')
         if key_filter == 0:
            channels = [allChannels[c] for c in allChannels]
            rc  = h['sig_'+str(s)+str(l)].Fill(np.array(channels).mean())
            if nSides==2:
               Sleft    = []
               Sright = []
               for c in allChannels:
                  if  nSiPMs > c:  # left side
                        Sleft.append(allChannels[c])
                  else:
                        Sright.append(allChannels[c])
               rc  = h['left_sig_'+str(s)+str(l)].Fill(np.array(Sleft).mean())
               rc  = h['right_sig_'+str(s)+str(l)].Fill(np.array(Sright).mean())
               rc  = h['sqrt_sig_'+str(s)+str(l)].Fill(np.sqrt(np.array(Sleft).mean()*np.array(Sright).mean()))

#
    
 S = {1:[1800,800,2,1],2:[1800,1500,2,3],3:[1800,1800,2,4]}
 type_of_hist = ["", "left_", "right_", "sqrt_"]
 for TYPE in type_of_hist:
   for s in S:
      ut.bookCanvas(h,TYPE + 'signal'    +str(s),' ',S[s][0],S[s][1],S[s][2],S[s][3])
      h[TYPE + 'signal'    +str(s)].SetLogy()
      for l in range(systemAndPlanes[s]):
         n = l+1
         if s==3 and n==7: n=8
         tag = str(s)+str(l)
         tc = h[TYPE + 'signal'+str(s)].cd(n)
         h[TYPE + 'sig_'+tag].Draw()

   for s in S:
      h[TYPE + 'signal'    +str(s)].Update()
      myPrint(h[TYPE + 'signal'    +str(s)],TYPE + 'mip_signal_mean_key_' + str(key_filter) + "_" + sdict[s])


def DS_track(HITS):
	# check for low occupancy and enough hits in DS stations
   stations = {}
   # for s in systemAndPlanes:
   for plane in range(systemAndPlanes[3]): 

      stations[30+plane] = {}
    # k=-1
   # for i, aHit in enumerate(eventTree.Digi_MuFilterHit):
   for i, aHit in enumerate(HITS):
      # k+=1
      if not aHit.isValid(): continue
      detID=aHit.GetDetectorID()
      subsystem, plane, bar = muAna.parseDetID(detID)
      if subsystem != 3: continue
      key=subsystem*10+plane
      #print(aHit)
      stations[key][i]=aHit
   if not len(stations[30])*len(stations[31])*len(stations[32])*len(stations[33]) == 1: return (-1,-1) # If not 1 hit in each DS plane
	#	build trackCandidate
   hitlist = {}
   for p in range(30,34):
      k = list(stations[p].keys())[0]
      hitlist[k] = stations[p][k]
   theTrack = trackTask.fitTrack(hitlist)
   if theTrack.getFitStatus().isFitConverged():
      print("FIT IS CONVERGED")
   return theTrack	

def MIP_study_2(Nev = options.nEvents, file_name = "", oneUShitperPlane = True, key_filter = 0):
 # with DS tracking 
 # veto system 2 layers with 7 bars and 8 sipm channels on both ends
 # US system 5 layers with 10 bars and 8 sipm channels on both ends
 # DS system horizontal(3) planes, 60 bars, readout on both sides, single channel
 #                         vertical(4) planes, 60 bar, readout on top, single channel
 for s in systemAndPlanes:
    for l in range(systemAndPlanes[s]):
       ut.bookHist(h,'sig_'+str(s*10+l),'signal / plane '+str(s*10+l),200,0.0,50.)
       ut.bookHist(h,'left_sig_'+str(s*10+l),'left_signal / plane '+str(s*10+l),200,0.0,50.)
       ut.bookHist(h,'right_sig_'+str(s*10+l),'right_signal / plane '+str(s*10+l),200,0.0,50.)
       ut.bookHist(h,'sqrt_sig_'+str(s*10+l),'sqrt_signal / plane '+str(s*10+l),200,0.0,50.)
 
 N=-1
 Tprev = 0
 if Nev < 0 : Nev = eventTree.GetEntries()
 eventTree.GetEvent(0)
 US_data = []
 DS_data = []
 for i, event in enumerate(eventTree):
    N+=1
    if N%options.heartBeat == 0: print('event ',N,' ',time.ctime())
    if N>Nev: break
    if oneUShitperPlane: 
      if not muAna.OneHitPerUS(eventTree.Digi_MuFilterHits): 
         continue
   #  #print("!!!!!!!!!!!!!!!!!!")
   #  theTrack = DS_track(eventTree.Digi_MuFilterHits)
   #  #print(theTrack, stations)
   #  if theTrack == (-1,-1):
   #    continue
   #  #print("STRANGE TRACK")
   #  if not hasattr(theTrack, "getFittedState"):
   #    continue
   #  #print("ALMOST GOOD TRACK")
   #  if not theTrack.getFitStatus().isFitConverged():
   #    theTrack.Delete()
   #    continue			
	# 	# if not hasattr(theTrack, "getFittedState"): continue
   #  state=theTrack.getFittedState()
   #  pos=state.getPos()
   #  mom=state.getMom()
   #  print("GOOD TRACK")


   
    for j, aHit in enumerate(eventTree.Digi_MuFilterHits):
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
      nSiPMs = aHit.GetnSiPMs()
      nSides  = aHit.GetnSides()
# check left/right
      allChannels = map2Dict(aHit,'GetAllSignals')
      if key_filter == 0:
         channels = [allChannels[c] for c in allChannels]
         ch = 0 
         for c in allChannels:
               if allChannels[c] != 0:
                  ch += 1
         if ch < 10: continue
         rc  = h['sig_'+str(s)+str(l)].Fill(np.array(channels).mean())
         if nSides==2:
            Sleft    = []
            Sright = []
            for c in allChannels:
               if  nSiPMs > c:  # left side
                     Sleft.append(allChannels[c])
               else:
                     Sright.append(allChannels[c])
            rc  = h['left_sig_'+str(s)+str(l)].Fill(np.array(Sleft).mean())
            rc  = h['right_sig_'+str(s)+str(l)].Fill(np.array(Sright).mean())
            rc  = h['sqrt_sig_'+str(s)+str(l)].Fill(np.sqrt(np.array(Sleft).mean()*np.array(Sright).mean()))

#
    
 S = {1:[1800,800,2,1],2:[1800,1500,2,3],3:[1800,1800,2,4]}
 type_of_hist = ["", "left_", "right_", "sqrt_"]
 for TYPE in type_of_hist:
   for s in S:
      ut.bookCanvas(h,TYPE + 'signal'    +str(s),' ',S[s][0],S[s][1],S[s][2],S[s][3])
      h[TYPE + 'signal'    +str(s)].SetLogy()
      for l in range(systemAndPlanes[s]):
         n = l+1
         if s==3 and n==7: n=8
         tag = str(s)+str(l)
         tc = h[TYPE + 'signal'+str(s)].cd(n)
         h[TYPE + 'sig_'+tag].Draw()

   for s in S:
      h[TYPE + 'signal'    +str(s)].Update()
      myPrint(h[TYPE + 'signal'    +str(s)],TYPE + 'mip_signal_mean_key_' + str(key_filter) + "_" + sdict[s])




#eventTime_for_filter_trash(5000000)
#eventTime(5000000)

#custom_ana_1(1000000)
#USshower(10000)
#USshower_custom(1000000)
#Mufi_hitMaps(1000000)
#eventTime(-1)

#execute_usshower_filtered(500000*2) # 3d bar distribution filtered
#USshower(1000) # ordinary USshower + output files with qdc/event distribution
#Mufi_hitMaps(1000)
#Mufi_hitMaps_qdc(1000000)
#Mufi_hitMaps_qdc_sipm_filter(1000000) # qdc distribution sipm filter
#Mufi_hitMaps_qdc_5planes_hits(1000000)
#Mufi_hitMaps_mean_qdc(1000000, 1)
MIP_study_2(-1, oneUShitperPlane = False)
#execute_mufi_hit_maps_filtered(1000000) # qdc distribution filtered
   
   
   #print("plane", plane, "bar", bar) 
