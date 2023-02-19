import ROOT,os,csv
import numpy as np


Tkey  = ROOT.std.vector('TString')()
Ikey   = ROOT.std.vector('int')()
Value = ROOT.std.vector('float')()
A,B,locA,locB,globA,globB    = ROOT.TVector3(),ROOT.TVector3(),ROOT.TVector3(),ROOT.TVector3(),ROOT.TVector3(),ROOT.TVector3()
verticalBarDict={0:1, 1:3, 2:5, 3:6}



def smallSiPMchannel(i):
    if i==2 or i==5 or i==10 or i==13: return True
    else: return False

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

def parseDetID(detID):
   subsystem=detID//10000
   if subsystem ==1 or subsystem==2:
      plane=detID%10000//1000
      bar=detID%1000
      return subsystem, plane, bar
   if subsystem == 3:
      bar=detID%1000
      if bar>59:
         plane=verticalBarDict[detID%10000//1000]
      elif bar<60:
         plane=2*detID%10000//1000
      return subsystem, plane, bar


def OneHitPerUS(DigiHits):
   #US_planes_all = {i for i in range(5)}
   USPlanes={}
	# for i, aHit in enumerate(eventTree.Digi_MuFilterHit):
   for i, aHit in enumerate(DigiHits):
      if not aHit.isValid(): continue
      detID=aHit.GetDetectorID()
      subsystem, plane, bar = parseDetID(detID)
      if subsystem != 2: continue
      if plane in USPlanes: 
        continue
      else:
        USPlanes[plane] = 1
      # !!!!!!!!!!!!
      USPlanes[plane] += 1
      # !!!!!!!!!!!1 :)
   #print(len(USPlanes.keys()))
   if len(USPlanes.keys()) != 5: return False
   for plane in USPlanes:
      if USPlanes[plane] != 1:
         return False
   print("!!!!!")
   return True

def OneHitUS1(DigiHits):
   #US_planes_all = {i for i in range(5)}
   USPlanes={}
	# for i, aHit in enumerate(eventTree.Digi_MuFilterHit):
   for i, aHit in enumerate(DigiHits):
      if not aHit.isValid(): continue
      detID=aHit.GetDetectorID()
      subsystem, plane, bar = parseDetID(detID)
      if subsystem != 2: continue
      if plane in USPlanes: 
        continue
      else:
        USPlanes[plane] = 1
      # !!!!!!!!!!!!
      USPlanes[plane] += 1
      # !!!!!!!!!!!1 :)
   #print(len(USPlanes.keys()))
   if 0 in USPlanes.keys():
      return True
   else:
      return False


def av_qdc(aHit, sipm_cut = "all", cut = 11):
    nSiPMs = aHit.GetnSiPMs()
    nSides  = aHit.GetnSides()
    allChannels = map2Dict(aHit,'GetAllSignals')
    channels = [allChannels[c] for c in allChannels]
    ch = 0 
    if sipm_cut == "all":
        for c in allChannels:
            if allChannels[c] != 0:
                ch += 1
        if ch < cut: return -1
    else:
        for c in allChannels:
            if allChannels[c] != 0 and not smallSiPMchannel(c):
                ch += 1
        if ch < cut: return -1
    if nSides==2:
        Sleft    = []
        Sright = []
    for c in allChannels:
        if allChannels[c] == 0: continue
        if  nSiPMs > c:  # left side
                Sleft.append(allChannels[c])
        else:
                Sright.append(allChannels[c])
    return np.sqrt(np.array(Sleft).mean()*np.array(Sright).mean())




def extrapolate(theTrack,z_mid):
    state = theTrack.getFittedState()
    pos = state.getPos()
    mom = state.getMom()
    slope_x = mom.x()/mom.z()
    slope_y = mom.y()/mom.z()
    x=pos.x()+slope_x*(z_mid-pos.z())
    y=pos.y()+slope_y*(z_mid-pos.z())
    return x,y

def residual(theTrack, detID, MuFilter, h_xy_track):
    MuFilter.GetPosition(detID,A,B)
    Z = 0.5*(A[2]+B[2])
    Y = 0.5*(A[1]+B[1])
    X = 0.5*(A[0]+B[0])
    xEx, yEx = extrapolate(theTrack,Z)
    h_xy_track.Fill(xEx, yEx)
    resY = yEx-Y
    resX = xEx-X
    return resX, resY



def qdc_left_right(aHit, sipm_cut = "all", cut = 11):
    nSiPMs = aHit.GetnSiPMs()
    nSides  = aHit.GetnSides()
    allChannels = map2Dict(aHit,'GetAllSignals')
    channels = [allChannels[c] for c in allChannels]
    ch = 0 
    if sipm_cut == "all":
        for c in allChannels:
            if allChannels[c] != 0:
                ch += 1
        if ch < cut: return -1
    else:
        for c in allChannels:
            if allChannels[c] != 0 and not smallSiPMchannel(c):
                ch += 1
        if ch < cut: return -1
    if nSides==2:
        Sleft    = []
        Sright = []
    for c in allChannels:
        if allChannels[c] == 0: continue
        if  nSiPMs > c:  # left side
                Sleft.append(allChannels[c])
        else:
                Sright.append(allChannels[c])
    return (np.array(Sleft).mean(), np.array(Sright).mean())

def fit_langau(hist,o,bmin,bmax):
    params = {0:'Width(scale)',1:'mostProbable',2:'norm',3:'sigma'}
    F = ROOT.TF1('langau',langaufun,0,200,4)
    for p in params: F.SetParName(p,params[p])
    rc = hist.Fit('landau','S'+o,'',bmin,bmax)
    res = rc.Get()
    if not res: return res
    F.SetParameter(2,res.Parameter(0))
    F.SetParameter(1,res.Parameter(1))
    F.SetParameter(0,res.Parameter(2))
    F.SetParameter(3,res.Parameter(2))
    F.SetParLimits(0,0,10)
    F.SetParLimits(1,0,100)
    F.SetParLimits(3,0,10)

    rc = hist.Fit(F,'S','',bmin,bmax)
    res = rc.Get().Parameter(1)
    return res



def langaufun(x,par):
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
    xlow = x[0] - sc * par[3]
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



