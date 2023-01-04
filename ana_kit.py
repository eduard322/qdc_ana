import ROOT,os,csv
import numpy as np


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



def av_qdc(aHit):
    nSiPMs = aHit.GetnSiPMs()
    nSides  = aHit.GetnSides()
    allChannels = map2Dict(aHit,'GetAllSignals')
    channels = [allChannels[c] for c in allChannels]
    ch = 0 
    for c in allChannels:
        if allChannels[c] != 0:
            ch += 1
    if ch < 10: return -1
    if nSides==2:
        Sleft    = []
        Sright = []
    for c in allChannels:
        if  nSiPMs > c:  # left side
                Sleft.append(allChannels[c])
        else:
                Sright.append(allChannels[c])
    return np.sqrt(np.array(Sleft).mean()*np.array(Sright).mean())


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
    res = rc.Get()
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



