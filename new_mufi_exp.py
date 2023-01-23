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



# def OneHitPerUS(DigiHits):
#    USPlanes={}
# 	# for i, aHit in enumerate(eventTree.Digi_MuFilterHit):
#    for i, aHit in enumerate(DigiHits):
#       if not aHit.isValid(): continue
#       detID=aHit.GetDetectorID()
#       subsystem, plane, bar = parseDetID(detID)
#       if subsystem != 2: continue
#       if plane in USPlanes: 
#         continue
#       else:
#         USPlanes[plane] = 0
#       # !!!!!!!!!!!!
#       USPlanes[plane] += 1
#       # !!!!!!!!!!!1 :)
#    for plane in USPlanes:
#       if USPlanes[plane] != 1:
#          return False
#    return True


def MIP_study(Nev = options.nEvents, oneUShitperPlane = True, withReco = False, DSTrack = True, optionTrack = "DS", title = "US_QDC_distributions"):

 # veto system 2 layers with 7 bars and 8 sipm channels on both ends
 # US system 5 layers with 10 bars and 8 sipm channels on both ends
 # DS system horizontal(3) planes, 60 bars, readout on both sides, single channel
 #                         vertical(4) planes, 60 bar, readout on top, single channel
 
 bin_min = 0.
 bin_max = 50.
 hist_list = {}
 hist_list_lr = {}
 QDC_list = []
 h_slope =  ROOT.TH2I("slope","slope", 100,-0.5,1.,100,-0.5,1.)
 for l in range(20, 25):
    #ut.bookHist(h,'sig_'+str(l),'signal / plane '+str(l),200,0.0,50.)
    h =  ROOT.TH1I("plane" + f"_{l}", "plane" + f"_{l}; E_dep [MIP];", 200, bin_min, bin_max)
    hist_list_lr[l] = {}
    for bar in range(10):
      h_lr =  ROOT.TH2I("plane" + f"_{l}_{bar}_lr", "plane" + f"_{l}_{bar}_lr", 100,-0.5,200.,100,-0.5,200.)
      hist_list_lr[l][bar] = h_lr
    #h_lr =  ROOT.TH2I("plane" + f"_{l}_lr", "plane" + f"_{l}_lr", 100,-0.5,200.,100,-0.5,200.)
    hist_list[l] = h
    #hist_list_lr[l] = h_lr

#  for l in range(30, 37):
#     h_ds =  ROOT.TH1I("plane" + f"_{l}", "plane" + f"_{l}", 200, bin_min, bin_max)
#     hist_list[l] = h_ds


 N=-1
 if Nev < 0 : Nev = eventTree.GetEntries()
 eventTree.GetEvent(0)
 #import pdb; pdb.set_trace()
 for event in eventTree:
    
    N+=1
    if N%options.heartBeat == 0: print('event ',N,' ',time.ctime())
    if N>Nev: break
    #print("!!!!1")
    #if oneUShitperPlane: 
    if not ana.OneHitUS1(eventTree.Digi_MuFilterHits): continue  
    #print("!!!!2")
    if DSTrack:
      if withReco:
         for aTrack in Reco_MuonTracks: aTrack.Delete()
         Reco_MuonTracks.Clear()
         if optionTrack=='DS': rc = trackTask.ExecuteTask("DS")
         else                         : rc = trackTask.ExecuteTask("Scifi")
      #print("!!!!3")
      if not Reco_MuonTracks.GetEntries()==1:
         #print(f"{N}: {trackTask.fittedTracks.GetEntries()}") 
         continue
      theTrack = Reco_MuonTracks[0]
      if not theTrack.getFitStatus().isFitConverged() and optionTrack!='DS':   # for H8 where only two planes / proj were avaiable
            continue
   # only take horizontal tracks
      state = theTrack.getFittedState(0)
      pos   = state.getPos()
      mom = state.getMom()
      slopeX= mom.X()/mom.Z()
      slopeY= mom.Y()/mom.Z()
      h_slope.Fill(slopeX, slopeY)
      if abs(slopeX)>0.25: continue   # 4cm distance, 250mrad = 1cm
      if abs(slopeY)>0.1: continue



    for aHit in event.Digi_MuFilterHits:
        if not aHit.isValid(): continue
        detID = aHit.GetDetectorID()

        s = detID//10000 # det number (s = 1 -- veto, s = 2 -- US, s = 3 -- DS)
        l  = (detID%10000)//1000  # plane number
        bar = (detID%1000)

        if s != 2 : continue
        nSiPMs = aHit.GetnSiPMs()
        nSides  = aHit.GetnSides()
        qdc_value = ana.av_qdc(aHit)
        if qdc_value == -1:
            continue
        QDC_list.append(qdc_value)
        #print(s*10 + l)
        hist_list[s*10 + l].Fill(qdc_value)
        q_l, q_r = ana.qdc_left_right(aHit)
        hist_list_lr[s*10 + l][bar].Fill(q_l, q_r)


 median = np.median(np.array(QDC_list))
 #langau fit
 c  = ROOT.TCanvas("US QDC distribution","US QDC distribution",0,0,1000,1000)
 c.Divide(2,3)
 for i, plane in enumerate(hist_list.keys()):
    for j in range(200):
      bin_old = hist_list[plane].GetBinContent(j)
      hist_list[plane].SetBinContent(j, bin_old/median)
    #print(i)
    c.cd(i+1)
   #  if fit == "langau":
   #    ana.fit_langau(hist_list[plane], str(plane), bin_min, bin_max)
   #  else:
       #hist_list[plane].Fit("pol5")
    hist_list[plane].Draw()
 c.SaveAs(title + "_qdc" + ".root")
 c.SaveAs(title + "_qdc" + ".png")

 c  = ROOT.TCanvas("Q_L vs. Q_R distribution","Q_L vs. Q_R distribution",0,0,1200,1200)
 c.Divide(3,4)
 for i, plane in enumerate(hist_list_lr[23].keys()):
    #print(i)
    c.cd(i+1)
   #  if fit == "langau":
   #    ana.fit_langau(hist_list[plane], str(plane), bin_min, bin_max)
   #  else:
   #    hist_list[plane].Fit("pol5")
    hist_list_lr[23][i].Draw("COLZ")
 c.SaveAs(title + "_qlqr" + ".root")
 c.SaveAs(title + "_qlqr" + ".pdf")

 File = ROOT.TFile.Open(f"{options.runNumber}_run_1.root", "RECREATE")
 File.WriteObject(h_slope, h_slope.GetTitle())
 def write_dict_to_file(dict_obj, File):
    for key in dict_obj.keys():
      File.WriteObject(dict_obj[key], dict_obj[key].GetTitle())
 
 hists_all = [hist_list] + [hist_list_lr[key] for key in hist_list_lr.keys()]
 #write_dict_to_file(hist_list, File)
 for H in hists_all:
   write_dict_to_file(H, File)
 File.Close()











#test_multiprocessing(1000)
#pid_create_output(500)

#pid_ana_mc(5000, "langau", "US_QDC_distributions_MC")


#qdc_dist(1000, "langau", "DS_QDC_run_90")
#Mufi_hitMaps(1000)
#qdc_dist_per_bar(1000000, "langau", "US_QDC_distributions_run_54")

#MIP_study(Nev = 10000, oneUShitperPlane = True, withReco = True, DSTrack = True, optionTrack = "DS", title = "100_gev_muon_ds_on_us_on_sipm_on_reco_on")

