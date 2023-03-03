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


def MIP_study(Nev = options.nEvents, oneUShitperPlane = True, withReco = False, DSTrack = True, res_cut = True, optionTrack = "DS", title = "US_QDC_distributions", label = "QDC", sipm_cut = "all", cut = 11, convert_sipm = False, fit = False, pion_mc = False):

 # veto system 2 layers with 7 bars and 8 sipm channels on both ends
 # US system 5 layers with 10 bars and 8 sipm channels on both ends
 # DS system horizontal(3) planes, 60 bars, readout on both sides, single channel
 #                         vertical(4) planes, 60 bar, readout on top, single channel
 
#  ut.bookHist(h,'bs','beam spot',100,-100.,10.,100,0.,80.)
#  ut.bookHist(h,'bsDS','beam spot',60,-0.5,59.5,60,-0.5,59.5)
#  for s in systemAndPlanes:
#     for l in range(systemAndPlanes[s]):
#       ut.bookHist(h,'Tsig_'+str(s*10+l),'signal / plane '+str(s*10+l),200,0.0,200.)
#  ut.bookHist(h,'slopes','track slopes',100,-0.1,0.1,100,-0.1,0.1)


 
 #MPVs_raw = [8.43, 8.14, 8.97, 7.32, 9.01]
 #MPVs_raw = [8.38, 8.16, 8.96, 7.31, 9.01]
 if pion_mc:
   MPVs_raw = [5 for i in range(5)]
 else:
   MPVs_raw = [7.97705, 7.90099, 7.96975, 8.14279, 7.98463] 
 MPVs = {plane: MPVs_raw[i] for i, plane in enumerate([20,21,22,23,24])} 
 MPVs_raw_l = [8.47, 8.93, 7.24, 10.03, 7.9][::-1]
 MPVs_raw_l = MPVs_raw
 #MPVs_raw_l = [np.array(MPVs_raw_l).mean() for _ in range(5)]
 MPVs_l = {plane: MPVs_raw_l[i] for i, plane in enumerate([20,21,22,23,24])}
 MPVs_raw_r = [9.50, 6.27, 10.12, 6,18, 9.46]
 MPVs_raw_r = MPVs_raw
 #MPVs_raw_r = [np.array(MPVs_raw_r).mean() for _ in range(5)]
 MPVs_r = {plane: MPVs_raw_r[i] for i, plane in enumerate([20,21,22,23,24])}
 
 if convert_sipm:
   MPVs_sipm = []
   with open("sipm_mpvs_full", "r") as f:
      x = 1
      while x:
            x = f.readline().split()
            if len(x) > 4 or len(x) == 0:
                  continue
            MPVs_sipm.append(list(map(int, x[:-1])) + [float(x[-1])])
   #print(np.array(MPVs_sipm))
   MPV_median = np.median(np.array(MPVs_sipm)[:,-1])
   MPVs_sipm = np.array([x[:-1] + [-x[-1] + MPV_median] for x in MPVs_sipm])
   # print(MPV_median)
   # print("*********")
   # print(MPVs_sipm[MPVs_sipm[:,0] == 200])
   # exit(0)


 bin_min = 0.
 if label == "MIP":
   bin_max = 250.
   bin_max_ql = 5.
 else:
   bin_max = 50.
   bin_max_ql = 200.
 hist_list = {}
 hist_list_l = {}
 hist_list_r = {}
 hist_list_lr = {}
 hist_list_bar = {}
 hist_list_sipm = {}

 QDC_list = []
 h_qdc_hit = ROOT.TH2I("qdc_vs_hit","QDC vs. Number of fired bars;Number of fired bars;QDC [MIP]", 100,0,50.,100,0,200)
 h_qdc_hit_norm = ROOT.TH2I("qdc_vs_hit_norm","QDC/Ntot vs. Number of fired bars;Number of fired bars;Etot/Ntot [MIP]", 100,0,50.,100,0,25)
 h_slope =  ROOT.TH2I("slope","slope;slope_x;slope_y", 100,-0.5,1.,100,-0.5,1.)
 h_xy_track =  ROOT.TH2I("track_xy","track_xy;x;y", 100,-100,100.,100,-100,100.)
 h_xy_track_res =  ROOT.TH2I("track_xy_res","track_xy_res;x;y", 100,-100,100.,100,-100,100.)
 h_xy_track_slope_abs =  ROOT.TH2I("track_xy_slope_abs","track_xy_slope_abs;x;y", 100,-1,1.,100,-1,1.)
 for l in range(20, 25):
    #ut.bookHist(h,'sig_'+str(l),'signal / plane '+str(l),200,0.0,50.)
    h =  ROOT.TH1I("plane" + f"_{l}", "plane" + f"_{l}; {label};", 200, bin_min, bin_max)
    h_l =  ROOT.TH1I("plane" + f"_{l}", "plane" + f"_{l}; {label};", 200, bin_min, bin_max)
    h_r =  ROOT.TH1I("plane" + f"_{l}", "plane" + f"_{l}; {label};", 200, bin_min, bin_max)
    hist_list_lr[l] = {}
    hist_list_bar[l] = {}
    hist_list_sipm[l] = {}
    for bar in range(10):
      h_lr =  ROOT.TH2I("plane" + f"_{l}_{bar}_lr", "plane" + f"_{l}_{bar}_lr; Q_L [{label}]; Q_R [{label}]", 100,-0.5,bin_max_ql,100,-0.5,bin_max_ql)
      h_bar =  ROOT.TH1I("plane" + f"_{l}_{bar}", "plane" + f"_{l}_{bar}; {label};", 200, bin_min, bin_max)
      hist_list_lr[l][bar] = h_lr
      hist_list_bar[l][bar] = h_bar
      hist_list_sipm[l][bar] = {}
      for sipm in range(16):
         h_sipm =  ROOT.TH1I("plane" + f"_{l}_{bar}_{sipm}", "plane" + f"_{l}_{bar}_{sipm}; {label};", 200, bin_min, bin_max)
         hist_list_sipm[l][bar][sipm] = h_sipm
    #h_lr =  ROOT.TH2I("plane" + f"_{l}_lr", "plane" + f"_{l}_lr", 100,-0.5,200.,100,-0.5,200.)
    hist_list[l] = h
    hist_list_l[l] = h_l
    hist_list_r[l] = h_r
    #hist_list_lr[l] = h_lr





 N=-1
 if Nev < 0 : Nev = eventTree.GetEntries()
 eventTree.GetEvent(0)
 #import pdb; pdb.set_trace()
 for event in eventTree:
    
    N+=1
    if N%1000 == 0: print('event ',N,' ',time.ctime())
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
      h_xy_track_slope_abs.Fill(mom.x()/mom.Mag(),mom.y()/mom.Mag())

    number_of_hits = 0
    qdc_per_event = 0
    for aHit in event.Digi_MuFilterHits:
        if not aHit.isValid(): continue
        detID = aHit.GetDetectorID()

        s = detID//10000 # det number (s = 1 -- veto, s = 2 -- US, s = 3 -- DS)
        l  = (detID%10000)//1000  # plane number
        bar = (detID%1000)

        if s != 2 : continue
        if res_cut:
         resx, resy = ana.residual(theTrack, detID, MuFilter, h_xy_track)
         h_xy_track_res.Fill(resx, resy)
         #res = np.sqrt(resx**2 + resy**2)
         res = np.abs(resy)
         if res >= 6.:
            continue
        nSiPMs = aHit.GetnSiPMs()
        nSides  = aHit.GetnSides()
      ####################
        if convert_sipm:
         MPVs_sipm_slice = MPVs_sipm[MPVs_sipm[:,0] == (s*10 + l)*10 + bar]
         #print((s*10 + l)*10 + bar)
         sipm_ids = MPVs_sipm_slice[:, 1]
         shift_coef = MPVs_sipm_slice[:, -1]
         MPVs_sipm_slice = {int(key): shift_coef[i] for i, key in enumerate(sipm_ids)}
        else:
         MPVs_sipm_slice = {i:0 for i in range(16)}
        #print(MPVs_sipm_slice)
      ####################
        qdc_value = ana.av_qdc(aHit, sipm_cut, cut, MPVs_sipm_slice)
      #   if pion_mc:
      #    if qdc_value < 2.:
      #       continue
        if qdc_value == -1:
            continue

        allChannels = map2Dict(aHit,'GetAllSignals')
        for Si in allChannels:
               if smallSiPMchannel(Si):
                  qdc_1 = allChannels[Si] + np.array(list(MPVs_sipm_slice.values())).mean()
               else:
                  qdc_1 = allChannels[Si] + MPVs_sipm_slice[Si]
               if qdc_1 == -1:
                  continue
               hist_list_sipm[s*10 + l][bar][Si].Fill(qdc_1)
        q_l, q_r = ana.qdc_left_right(aHit, sipm_cut, cut, MPVs_sipm_slice)
        if label == "MIP":
            qdc_value /= MPVs[s*10 + l]
            q_l /= MPVs_l[s*10 + l]
            q_r /= MPVs_r[s*10 + l]
        #QDC_list.append(qdc_value)
        #print(s*10 + l)
        hist_list[s*10 + l].Fill(qdc_value)
        hist_list_l[s*10 + l].Fill(q_l)
        hist_list_r[s*10 + l].Fill(q_r)
        hist_list_lr[s*10 + l][bar].Fill(q_l, q_r)
        hist_list_bar[s*10 + l][bar].Fill(qdc_value)
        qdc_per_event += qdc_value
        number_of_hits += 1 
    if number_of_hits > 0:
      h_qdc_hit.Fill(number_of_hits, qdc_per_event)
      h_qdc_hit_norm.Fill(number_of_hits, qdc_per_event/(number_of_hits))  


 File = ROOT.TFile.Open(f"{title}_{options.runNumber}_run_2_test.root", "RECREATE")

 c  = ROOT.TCanvas("QDC vs. NOH","QDC vs. NOH",0,0,1000,1000)
 ROOT.gPad.SetLogz()
 h_qdc_hit.Draw('colz')
 c.Write()
 c  = ROOT.TCanvas("QDC/NOH vs. NOH","QDC/NOH vs. NOH",0,0,1000,1000)
 ROOT.gPad.SetLogz()
 h_qdc_hit_norm.Draw('colz')
 c.Write()

 c  = ROOT.TCanvas("US QDC distribution","US QDC distribution",0,0,1000,1000)
 c.Divide(2,3)
 for i, plane in enumerate(hist_list.keys()):

    c.cd(i+1)
    ana.fit_langau(hist_list[plane], str(plane), bin_min, bin_max)
    hist_list[plane].Draw()
#  c.SaveAs(title + "_qdc" + ".root")
#  c.SaveAs(title + "_qdc" + ".pdf")
 c.Write()


 for pl in hist_list_lr.keys():
   c  = ROOT.TCanvas(f"Q_L vs. Q_R distribution. Plane {pl}","Q_L vs. Q_R distribution",0,0,1200,1200)
   c.Divide(3,4)
   for i, br in enumerate(hist_list_lr[pl].keys()):
      c.cd(i+1)
      hist_list_lr[pl][br].Draw("COLZ")
   c.Write()
#  c.SaveAs(title + "_qlqr" + ".root")
#  c.SaveAs(title + "_qlqr" + ".pdf")

 c  = ROOT.TCanvas("US QDC_l distribution","US QDC_l distribution",0,0,1000,1000)
 c.Divide(2,3)
 for i, plane in enumerate(hist_list_l.keys()):
    c.cd(i+1)
    #ana.fit_langau(hist_list_l[plane], str(plane), bin_min, bin_max)
    hist_list_l[plane].Draw()
 c.Write()
#  c.SaveAs(title + "_qdc_l" + ".root")
#  c.SaveAs(title + "_qdc_l" + ".pdf")


 c  = ROOT.TCanvas("US QDC_r distribution","US QDC_r distribution",0,0,1000,1000)
 c.Divide(2,3)
 for i, plane in enumerate(hist_list_r.keys()):
    c.cd(i+1)
    #ana.fit_langau(hist_list_r[plane], str(plane), bin_min, bin_max)
    hist_list_r[plane].Draw()
 c.Write()
#  c.SaveAs(title + "_qdc_r" + ".root")
#  c.SaveAs(title + "_qdc_r" + ".pdf")

 mpv_bars = []
 barID = []
 with open("bar_mpvs", "w") as f:
   for pl in hist_list_bar.keys():
      c  = ROOT.TCanvas(f"US QDC distribution. Plane {pl}", f"US QDC distribution. Bar {pl}",0,0,1000,1000)
      c.Divide(3,4)
      for i, br in enumerate(hist_list_bar[pl].keys()):
         c.cd(i+1)
         if fit:
            res = ana.fit_langau(hist_list_bar[pl][br], str(pl*10 + br), bin_min, bin_max)
            mpv_bars.append(res)
            barID.append(pl*10 + br)
            print(pl*10 + br, res, file = f)
         hist_list_bar[pl][br].Draw("COLZ")
      c.Write()
 
 if fit:
   c = ROOT.TCanvas("Bar mpvs", "Bar mpvs", 1000, 1000)
   gr = ROOT.TGraph(len(barID), array("d", barID), array("d", mpv_bars))
   gr.SetMarkerStyle(23)
   gr.SetMarkerColor(ROOT.kBlue)
   gr.GetXaxis().SetTitle("barID")
   gr.GetYaxis().SetTitle("MPVs [QDC]")
   gr.SetTitle("Bar MPVs")
   gr.Draw("AP")
   c.Write()



 mpv_sipms = []
 sipmID = []
 k_sipm = 0
 with open("sipm_mpvs", "w") as f:
   for pl in hist_list_sipm.keys():
      for br in hist_list_sipm[pl].keys():
         c  = ROOT.TCanvas(f"US QDC SiPM distribution. Bar {pl*10 + br}", f"US QDC distribution. Bar {pl*10 + br}",0,0,1000,1000)
         c.Divide(4,4)
         #with open("output_par_dsoff_" + str(bar_id), "w") as f:
         for i, Si in enumerate(hist_list_sipm[pl][br].keys()):
            c.cd(i+1)
            if not smallSiPMchannel(Si):
               if fit:
                  res = ana.fit_langau(hist_list_sipm[pl][br][Si], str(Si), bin_min, bin_max)
                  sipmID.append(k_sipm)
                  mpv_sipms.append(res)
                  print(pl*10 + br, Si, k_sipm, res, file = f)
            hist_list_sipm[pl][br][Si].Draw()           
            k_sipm += 1
         c.Write()
         
         
 if fit:
   c = ROOT.TCanvas("SiPM mpvs", "SiPM mpvs", 1000, 1000)
   gr = ROOT.TGraph(len(sipmID), array("d", sipmID), array("d", mpv_sipms))
   gr.SetMarkerStyle(23)
   gr.SetMarkerColor(ROOT.kBlue)
   gr.GetXaxis().SetTitle("SiPMID")
   gr.GetYaxis().SetTitle("MPVs [QDC]")
   gr.SetTitle("SiPM MPVs")
   gr.Draw("AP")
   #c.Write()



 File.WriteObject(h_slope, h_slope.GetTitle())
 File.WriteObject(h_xy_track, h_xy_track.GetTitle())
 File.WriteObject(h_xy_track_res, h_xy_track_res.GetTitle())
 File.WriteObject(h_xy_track_slope_abs, h_xy_track_slope_abs.GetTitle())
 File.WriteObject(h_slope, h_slope.GetTitle())
#  def write_dict_to_file(dict_obj, File):
#     for key in dict_obj.keys():
#       File.WriteObject(dict_obj[key], dict_obj[key].GetTitle())
 
#  hists_all = [hist_list] + [hist_list_lr[key] for key in hist_list_lr.keys()]
#  #write_dict_to_file(hist_list, File)
#  for H in hists_all:
#    write_dict_to_file(H, File)
 File.Close()


def SiPM_study(Nev = options.nEvents, oneUShitperPlane = True, withReco = False, DSTrack = True, optionTrack = "DS", title = "US_QDC_distributions", label = "QDC"):

 # veto system 2 layers with 7 bars and 8 sipm channels on both ends
 # US system 5 layers with 10 bars and 8 sipm channels on both ends
 # DS system horizontal(3) planes, 60 bars, readout on both sides, single channel
 #                         vertical(4) planes, 60 bar, readout on top, single channel
 
   bin_min = 0.
   bin_max = 50.
   hist_list = {}
   QDC_list = []
   h_xy_track =  ROOT.TH2I("track_xy","track_xy;x;y", 100,-100,100.,100,-100,100.)
   h_slope =  ROOT.TH2I("slope","slope", 100,-0.5,1.,100,-0.5,1.)
   plane_id = 20
   bar_id = 209
   for Si in range(16):
      h =  ROOT.TH1I("SiPM" + f"_{Si}", "SiPM" + f"_{Si}; {label};", 200, bin_min, bin_max)
      hist_list[Si] = h


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
         if abs(slopeX)>0.1: continue   # 4cm distance, 250mrad = 1cm
         if abs(slopeY)>0.1: continue


      for aHit in event.Digi_MuFilterHits:
         if not aHit.isValid(): continue
         detID = aHit.GetDetectorID()

            # residual cut
            # residual cut
         resx, resy = ana.residual(theTrack, detID, MuFilter, h_xy_track)
         res = np.sqrt(resx**2 + resy**2)
         if res >= 6.:
            continue

         s = detID//10000 # det number (s = 1 -- veto, s = 2 -- US, s = 3 -- DS)
         l  = (detID%10000)//1000  # plane number
         bar = (detID%1000)

         if s != 2 : continue
         if l != 0 : continue
         if bar != 9 : continue
         nSiPMs = aHit.GetnSiPMs()
         nSides  = aHit.GetnSides()
         allChannels = map2Dict(aHit,'GetAllSignals')
         for Si in allChannels:
               qdc_value = allChannels[Si]
               if qdc_value == -1:
                  continue
               hist_list[Si].Fill(qdc_value)
   #langau fit
   c  = ROOT.TCanvas("US QDC distribution","US QDC distribution",0,0,1000,1000)
   c.Divide(4,4)
   with open("output_par_dsoff_" + str(bar_id), "w") as f:
      for i, Si in enumerate(hist_list.keys()):
         c.cd(i+1)
         res = ana.fit_langau(hist_list[Si], str(Si), bin_min, bin_max)
         hist_list[Si].Draw()
         print(bar_id, Si, res, "\n", file = f)
   c.SaveAs(title + "_qdc" + ".root")
   c.SaveAs(title + "_qdc" + ".pdf")



   File = ROOT.TFile.Open(f"{options.runNumber}_run_1.root", "RECREATE")
   File.WriteObject(h_slope, h_slope.GetTitle())
   File.WriteObject(h_xy_track, h_xy_track.GetTitle())
   def write_dict_to_file(dict_obj, File):
      for key in dict_obj.keys():
         File.WriteObject(dict_obj[key], dict_obj[key].GetTitle())
   
   hists_all = [hist_list]
   #write_dict_to_file(hist_list, File)
   for H in hists_all:
      write_dict_to_file(H, File)
   File.Close()






def convert_data(Nev = options.nEvents, oneUShitperPlane = True, withReco = True, DSTrack = True, optionTrack = "DS"):

 # veto system 2 layers with 7 bars and 8 sipm channels on both ends
 # US system 5 layers with 10 bars and 8 sipm channels on both ends
 # DS system horizontal(3) planes, 60 bars, readout on both sides, single channel
 #                         vertical(4) planes, 60 bar, readout on top, single channel
 

 def add_zeros(channels, LEN = 16):
   return {key: channels[key] if key in channels.keys() else 0 for key in range(LEN)}

 h_xy_track =  ROOT.TH2I("track_xy","track_xy;x;y", 100,-100,100.,100,-100,100.)
 output = "/eos/user/u/ursovsnd/private/SND_Data/qdc_study/preprocessed_data/"
 #output = ""
 N=-1
 if Nev < 0 : Nev = eventTree.GetEntries()
 eventTree.GetEvent(0)
 #import pdb; pdb.set_trace()
 with open(output + "converted_data_muon_100_GeV_H8_cut_on", "w") as f:
   for event in eventTree:
      
      N+=1
      if N%options.heartBeat == 0: print('event ',N,' ',time.ctime())
      if N>Nev: break

      if not ana.OneHitUS1(eventTree.Digi_MuFilterHits): continue  
      if DSTrack:
         if withReco:
            for aTrack in Reco_MuonTracks: aTrack.Delete()
            Reco_MuonTracks.Clear()
            rc = trackTask.ExecuteTask("DS")
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
         if abs(slopeX)>0.25: continue   # 4cm distance, 250mrad = 1cm
         if abs(slopeY)>0.1: continue


      for hit_number, aHit in enumerate(event.Digi_MuFilterHits):
         if not aHit.isValid(): continue
         detID = aHit.GetDetectorID()

         s = detID//10000 # det number (s = 1 -- veto, s = 2 -- US, s = 3 -- DS)
         l  = (detID%10000)//1000  # plane number
         bar = (detID%1000)

         if s != 2 : continue
         resx, resy = ana.residual(theTrack, detID, MuFilter, h_xy_track)
         #res = np.sqrt(resx**2 + resy**2)
         res = np.abs(resy)
         if res >= 6.:
            continue


         allChannels = add_zeros(ana.map2Dict(aHit,'GetAllSignals'))
         allTimes = add_zeros(ana.map2Dict(aHit,'GetAllTimes'))
         #print(allChannels)
         channels = [allChannels[c] if allChannels[c] != 0 else 0 for c in allChannels]
         times = [allTimes[c] if allTimes[c] != 0 else 0 for c in allTimes]
         nSiPMs = aHit.GetnSiPMs()
         #print(nSiPMs)
         # if allTimes[c] != 0 else 0
         meta_data = [N, hit_number, detID]
         string_output = ""
         string_output += "\t".join(list(map(str, meta_data)))
         string_output += "\t"
         string_output += "\t".join(list(map(str, channels)))
         string_output += "\t"
         string_output += "\t".join(list(map(str, times)))
         string_output += "\n"
         
         f.write(string_output)





#test_multiprocessing(1000)
#pid_create_output(500)

#pid_ana_mc(5000, "langau", "US_QDC_distributions_MC")


#qdc_dist(1000, "langau", "DS_QDC_run_90")
#Mufi_hitMaps(1000)
#qdc_dist_per_bar(1000000, "langau", "US_QDC_distributions_run_54")

#title = "100_gev_muon_ds_on_us_on_sipm_on_reco_on_largesipm_cut_11_res_cut"
# title = "300_gev_pion_ds_off_us_on_sipm_on_reco_on_largesipm_cut_11_test"
# title = "pion_mc_pe_conversion"
# title = "muon_100_gev_H8_no_cut"
title = "pion_mc_no_cut_pe_saturation_2_low_cut"
# title = "pion_300_H8_large_scale"
cut_key = False
MIP_study(Nev = 10000, 
   oneUShitperPlane = cut_key, 
   withReco = cut_key, 
   DSTrack = cut_key,
   res_cut = cut_key, 
   optionTrack = "DS", 
   title = title,
   label = "MIP",
   sipm_cut = "large111", 
   cut = 11,
   convert_sipm = False,
   pion_mc = True)

# SiPM_study(Nev = -1, oneUShitperPlane = True, withReco = True, DSTrack = True, optionTrack = "DS", title = "SiPM_100_gev_muon_pl_0_bar_9_ds_on_us_on_res_on_all")

#convert_data(-1)
