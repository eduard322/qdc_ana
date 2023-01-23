import ROOT,os,csv

# class mufi_analysis: # I could make several classes for generic functions like parseDetID and more specific ones for 
	# mapping SiPM channels to the whole system

A, B=ROOT.TVector3(), ROOT.TVector3()
systemAndPlanes = {1:2,2:5,3:6}
verticalBarDict={0:1, 1:3, 2:5, 3:6}

def BuildBarLengths(MuFilter):
   Vetobarlength = MuFilter.GetConfParF('MuFilter/VetoBarX')
   USbarlength = MuFilter.GetConfParF('MuFilter/UpstreamBarX')
   DSbarlength_hor = MuFilter.GetConfParF('MuFilter/DownstreamBarX')
   DSbarlength_vert = MuFilter.GetConfParF('MuFilter/UpstreamBarY_ver')
   barlengths={1:Vetobarlength, 2:USbarlength, 3:DSbarlength_hor, 4:DSbarlength_vert}
   return barlengths

def BuildzPos(MuFilter):
   zPos={}
   for s in systemAndPlanes:
      for plane in range(systemAndPlanes[s]):
         bar = 4
         p = plane
         if s==3 and plane%2==0:  
            bar = 90
            p = plane//2
         if s==3 and plane%2==1:
            bar = 30
            p = plane//2
         MuFilter.GetPosition(s*10000+p*1000+bar,A,B)
         zPos[s*10+plane] = (A.Z()+B.Z())/2.	
   return zPos

def IsSmallSiPMchannel(i):
	if i==2 or i==5 or i==10 or i==13: return True

	else: return False

def dist2BarEnd(MuFilter, detID, nSides, pred, dim):
	if nSides == 2 and dim=='x':
		MuFilter.GetPosition(detID, A, B)
		dxL=pred-B.x()
		dxR=pred-A.x()
		return dxL, dxR

	elif nSides==2 and dim=='y':
		MuFilter.GetPosition(detID, A, B)
		dyT=A.y()-pred
		return dyT		

	elif nSides == 1 and dim=='y':
		MuFilter.GetPosition(detID, A, B)
		dyT=A.y()-pred
		return dyT

	elif nSides == 1 and dim=='y':
		MuFilter.GetPosition(detID, A, B)
		dyT=A.y()-pred
		return dyT   

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

	      
def GetAvgT(mufiHit, side='both'):
   value=[0, 0]
   count=[0, 0]
   nSiPMs=mufiHit.GetnSiPMs()
   times=mufiHit.GetAllTimes()
   s, p, b = parseDetID(mufiHit.GetDetectorID())
   for element in times:
      SiPM, time = element 
      if SiPM<nSiPMs:
         value[0]+=time
         count[0]+=1
      else:
         value[1]+=time
         count[1]+=1
   if s == 2:
      if count[0] != 0 and count[1] != 0:
         if side=='both':
            average = 0.5*(value[0]/count[0]+value[1]/count[1])
            return average
         if side=='L':
            return value[0]/count[0]
         elif side == 'R':
            return value[1]/count[1]
      else: return -999.
   elif s == 3 and b < 60: 
      if side=='both':
         average = 0.5*(value[0]/count[0]+value[1]/count[1])
         return average
      elif side == 'L': return value[0]/count[0]
      elif side == 'R': return value[1]/count[1]

   # elif s == 3 and b > 59: 
   #    return

def GetChannelVal(SiPM, chs):
	for entry in chs:
	   fSiPM, val = entry
	   if fSiPM == SiPM:
	      return val
	return -999.

def TotalChVal_sides(chs):
   left, right=0., 0.
   for entry in chs:
      SiPM, val = entry
      if SiPM<8: left+=val
      elif SiPM>=8: right+=val
   return left, right

def GetOverallSiPMNumber(detID, SiPM):
   s, p, b = parseDetID(detID)
   if s==1: nSiPMs, SiPMs_plane=6, 52 # is it?
   elif s==2: nSiPMs, SiPMs_plane=16, 160
   elif s==3: nSiPMs, SiPMs_plane=1, 60 # wrong

   return SiPM+nSiPMs*b+p*SiPMs_plane   

def OneHitPerUS(DigiHits):
   USPlanes={}
	# for i, aHit in enumerate(eventTree.Digi_MuFilterHit):
   for i, aHit in enumerate(DigiHits):
      if not aHit.isValid(): continue
      detID=aHit.GetDetectorID()
      subsystem, plane, bar = parseDetID(detID)
      if subsystem != 2: continue
      if plane in USPlanes: continue
      USPlanes[plane]=1
   for plane in USPlanes:
      if USPlanes[plane] != 1:
         return False
   return True

# def DS_track(DigiHits):
#    systemAndPlanes = {1:2,2:5,3:6}
# 	# check for low occupancy and enough hits in DS stations
#    stations = {}
#    # for s in systemAndPlanes:
#    for plane in range(systemAndPlanes[3]+1): 

#       stations[30+plane] = {}
#     # k=-1
#    # for i, aHit in enumerate(eventTree.Digi_MuFilterHit):
#    for i, aHit in enumerate(DigiHits):
#       # k+=1
#       if not aHit.isValid(): continue
#       detID=aHit.GetDetectorID()
#       subsystem, plane, bar = parseDetID(detID)
#       if subsystem != 3: continue
#       key=subsystem*10+plane
#       #print(key)
#       stations[key][i]=aHit
#    if not len(stations[30])*len(stations[31])*len(stations[32])*len(stations[33]) == 1: return -1 # If not 1 hit in each DS plane
# 	#	build trackCandidate
#    hitlist = {}
#    for p in range(30,34):
#       k = list(stations[p].keys())[0]
#       hitlist[k] = stations[p][k]
#    theTrack = trackTask.fitTrack(hitlist)
#    print(theTrack)
#    return theTrack

def DS_track(DigiHits):
	# check for low occupancy and enough hits in DS stations
   stations = {}
   # for s in systemAndPlanes:
   for plane in range(systemAndPlanes[3]): 

      stations[30+plane] = {}
    # k=-1
   # for i, aHit in enumerate(eventTree.Digi_MuFilterHit):
   for i, aHit in enumerate(DigiHits):
      # k+=1
      if not aHit.isValid(): continue
      detID=aHit.GetDetectorID()
      subsystem, plane, bar = parseDetID(detID)
      if subsystem != 3: continue
      key=subsystem*10+plane
      stations[key][i]=aHit
   if not len(stations[30])*len(stations[31])*len(stations[32])*len(stations[33]) == 1: return -1 # If not 1 hit in each DS plane
	#	build trackCandidate
   hitlist = {}
   print(stations[30].keys(), len(stations[30].keys()))
   for p in range(30,34): 
      k = list(stations[p].keys())[0]
      hitlist[k] = stations[p][k]
   theTrack = trackTask.fitTrack(hitlist)
   return theTrack		

def delta_min_t(aHit):
   times = aHit.GetAllTimes()
   if len(times)==0: return -999.
   nSiPMs = aHit.GetnSiPMs()
   ts_L, ts_R = [], []
   for channel in times:
      SiPM, time = channel
      if SiPM<nSiPMs: ts_L.append(time)
      elif SiPM>=nSiPMs: ts_R.append(time)  
   return min(ts_L)-min(ts_R)

def GetnFiredSiPMs(aHit):
	nSiPMs=aHit.GetnSiPMs()

	nFiredSiPMs_left=0
	nFiredSiPMs_right=0
	channels=aHit.GetAllSignals()
	for ch in channels:
		SiPM, qdc = ch 
		if SiPM<nSiPMs: nFiredSiPMs_left+=1
		elif SiPM>=nSiPMs: nFiredSiPMs_right+=1
	return nFiredSiPMs_left, nFiredSiPMs_right 

def GetBarSlice(L, sliceL, xpred):
   nslices=int(L/sliceL)
   slice_num = int(xpred//sliceL)
   return slice_num   

def GetSiPMNumberInSystem_LandR(detID, SiPM): # 20000 SiPM 8 -> 8
   s, p, b = parseDetID(detID)
   if s==1: nSiPMs, SiPMs_plane=6, 52 # is it?
   elif s==2:
      nSiPMs, SiPMs_plane=16, 160
      return SiPM+nSiPMs*b+p*SiPMs_plane
   elif s==3: # Count left and right horizontal SiPMs consecutively
      nSiPMs, SiPMs_hor_plane, SiPMs_vert_plane=1, 120, 60
      if p not in verticalPlanes:
         total_SiPM = (p//2)*SiPMs_hor_plane+(p//2)*SiPMs_vert_plane+SiPM+2*b
         return total_SiPM
      else:
         if p==1: return SiPMs_hor_plane+b 
         elif p==3: return 2*SiPMs_hor_plane+SiPMs_vert_plane+b
         elif p==5: return 3*SiPMs_hor_plane+2*SiPMs_vert_plane+b
         elif p==6: return 3*SiPMs_hor_plane+3*SiPMs_vert_plane+b    
         
def GetSiPMNumberInPlane_LandR(detID, SiPM):
   s, p, b = parseDetID(detID)
   if s != 2: 
      print('AAAAAAAHHHHHH')
   return SiPM+b*16 

def GetSiPMNumberInPlane_LTR(detID, SiPM):
   s, p, b = parseDetID(detID)
   if s != 2: 
      print('AAAAAAAHHHHHH')
   return SiPM+b*8

def GetSiPMNumberInSystem_LTR(detID, SiPM): # 20000 SiPM 8 -> 400
   s, p, b = parseDetID(detID)
   if s==1: nSiPMs, SiPMs_plane=6, 52 # is it?
   elif s==2: nSiPMs, SiPMs_plane=8, 80
   elif s==3: nSiPMs, SiPMs_plane=1, 60 # wrong 

   if SiPM<nSiPMs:
      return SiPM+nSiPMs*b+p*SiPMs_plane
   elif SiPM>=nSiPMs:
      SiPM_r=400+SiPM%nSiPMs
      return SiPM_r+nSiPMs*b+p*SiPMs_plane

def GetDSH_average(stations):
   ts=[0., 0]
   for key in stations:
      if len(stations[key])==0: continue # By default, a key is made for all DS planes regardless of if they are present in the geometry
      if len(stations[key].values()) != 1:
         print('\nLen stations values not 1\n')
         return -999. # return if there are 2 hits asigned to 1 bar in stations
      hit_num, hit = list(stations[key].keys())[0], list(stations[key].values())[0]
      s,p,b=parseDetID(hit.GetDetectorID())
      if b>59:continue # Only considering horizontal bars
      times=hit.GetAllTimes()
      if len(times) != 2: continue # Only considering hit horizontal bars with both SiPMs firing
      for ch in times:
         SiPM, tdc = ch
         ts[0]+=tdc
         ts[1]+=1
   if ts[1]==0: return -998. # Can only happen if every entry in stations fails all selection criteria
   return ts[0]/ts[1]

def GetDeltaT(times, one_channel=None):
   # nSiPMs=aHit.GetnSiPMs()
   nSiPMs=8
   mean = [0,0]
   count = [0,0]
   # channels = aHit.GetAllTimes()
   for ch in times:
      SiPM, val = ch
      if one_channel != None:
         if not (SiPM == one_channel or SiPM == one_channel+nSiPMs): continue
      if IsSmallSiPMchannel(SiPM): continue
      if SiPM < nSiPMs:
         mean[0]+=val
         count[0]+=1
      else:
         mean[1]+=val
         count[1]+=1
   # print(count, mean)
   if count[0] != 0 and count[1] != 0:
      return (mean[0]/count[0]-mean[1]/count[1])/2.
   else: return -999.

def GetdtCalc(xpred, L, cs):
   left, right=list(filter(lambda x : x[0]<8, cs)), list(filter(lambda x : x[0]>7, cs))
   NL, NR=len(left), len(right)
   if NL == 0 or NR == 0: return -999.
   sumOfInverses=lambda x : sum( [1/i[1] for i in x] )
   return xpred/NL*sumOfInverses(left) - (L-xpred)/NR*sumOfInverses(right)

def Getcscint_i(path, fixed_ch, iteration):
   with open(path+'cscintvalues/cscint_'+fixed_ch+'.csv', 'r') as handle:
      reader=csv.reader(handle)
      alldata=[row for row in reader]
      if len(alldata)<iteration or len(alldata) == 0: return -999.
      else:
         data=alldata[iteration]
         return (float(data[1]), float(data[2]))

def Getcscint_chi2pNDF(path, fixed_ch, iteration):
   with open(path+'cscintvalues/cscint_'+fixed_ch+'.csv', 'r') as handle:
      reader=csv.reader(handle)
      alldata=[row for row in reader]
      if len(alldata)<iteration or len(alldata) == 0: return -999.
      else:
         data=alldata[iteration]
         return data[-1]

def GetAvgcscint_i(path, detID, SiPMs, iteration):
   cs=[]
   for SiPM in SiPMs:
      fixed_ch=('_').join( (str(detID), str(SiPM)) ) 
      val=Getcscint_i(path, fixed_ch, iteration)
      cs.append(val)
   
   cbar = sum([i[0] for i in cs])/len(cs)
   cbar_err = ROOT.TMath.Sqrt( sum( [i[1]**2 for i in cs] ) )
   return cbar, cbar_err

def GetLogParams_i(path, fixed_ch, iteration):
   if iteration==0:
      print('No log params calculated for iteration 0.')
      return 0
   with open(path+'Logparams/logparams_'+fixed_ch+'.csv', 'r') as handle:
      reader=csv.reader(handle)
      alldata=[row for row in reader]
      if len(alldata)<iteration-1 or len(alldata) == 0: return -999.
      else:
         data=alldata[iteration-1]
         return (float(data[1]),float(data[2]),float(data[3]),float(data[4]),float(data[5]),float(data[6]))

def GetPolyparams_i(path, fixed_ch, iteration):
   if iteration==0:
      print('No log params calculated for iteration 0.')
      return 0
   with open(path+'Polyparams/polyparams_'+fixed_ch+'.csv', 'r') as handle:
      reader=csv.reader(handle)
      alldata=[row for row in reader]
      if len(alldata)<iteration-1 or len(alldata)==0: return -999.
      else: 
         data=alldata[iteration-1][1:]
         data.pop() # We don't need to get the chi2 here.
         return tuple( [float(i) for i in data] )

def Gettimeresolution(path, fixed_ch, iteration):
   appendixdict={0:'_uncorrected', 1:'_corrected'}
   fname='timeresolution_'+fixed_ch+appendixdict[iteration]+'.csv'
   if not fname in os.listdir(path): return -999
   with open(path+fname, 'r') as f:
      reader=csv.reader(f)
      alldata=[row for row in reader]
      data=alldata[0]
   timeresolution=float(data[1]), float(data[2])
   return timeresolution

def GetMIPpeak(path, fixed_ch):
   fname='timewalk_dists_'+fixed_ch+'_poly.root'
   f=ROOT.TFile.Open(path+fname, 'READ')

   histname='dtcorrectedvqdc_'+fixed_ch+'_iteration1'
   if not histname in [k.GetName() for k in f.GetListOfKeys()]: return -999
   hist=f.Get(histname)
   xproj=hist.ProjectionX()
   #maxbin=xproj. # Currently all QDC histogram ranges are 0, 60 but this needs to change!
   res=fit_langau(xproj, 'Q', 0, 60)
   MPV=res.Parameter(1)
   MPV_err=res.ParError(1)
   return (MPV, MPV_err)

def Getchi2pNDF(path, fixed_ch):
   
   fname=path+'Polyparams/polyparams_'+fixed_ch+'.csv'
   if not os.path.exists(fname): return -999.
   with open(fname, 'r') as handle:
      reader=csv.reader(handle)
      alldata=[row for row in reader]
      data=alldata[0]
      chi2pNDF=data.pop() # We don't need to get the chi2 here.
      return float(chi2pNDF)

def Makechi2pNDFdict(path):
   path+='Polyparams/'
   res={}
   for filename in os.listdir(path):
      fixed_ch=filename[filename.find('2'):filename.find('.csv')]
      
      with open(path+filename, 'r') as handle:
         reader=csv.reader(handle)
         alldata=[row for row in reader]
         data=alldata[0]
         chi2pNDF=data.pop()
      res[fixed_ch]=float(chi2pNDF)
   sorted_tuples=sorted(res.items(), key=lambda x:x[1])
   sorted_d={k:v for k,v in sorted_tuples}
   return sorted_d

def GetAttenuationLengthData(path, fixed_ch):
   fname=path+'attenuationlengths/csvfiles/attenuationlength_'+fixed_ch+'.csv'
   if not os.path.exists(fname): return -999.
   with open(fname, 'r') as handle:
      reader=csv.reader(handle)
      alldata=[row for row in reader]
      data=alldata[0]
      return data

def GetPltThr(func, parameters):
   a, b, c =parameters
   func.SetParameters(a, b, c)
   thr = 0.
   for x in range(0, 6000):
      x_val=x/1000
      y_val = func.Eval(x_val)
      if ROOT.TMath.IsNaN(y_val) or not ROOT.TMath.Finite(y_val):
         thr=x_val
         break
   return x_val

def GetXcalculated(dt, L, cs, wanted=None):
   
   left, right=list(filter(lambda x : x[0]<8, cs)), list(filter(lambda x : x[0]>7, cs))
   NL, NR=len(left), len(right)
   if NL==0 or NR == 0: return -999.
   sumOfInverses=lambda x : sum( [1/i[1] for i in x] )
   A, B = 1/NL*sumOfInverses(left), 1/NR*sumOfInverses(right)
   xcalc = (dt+L*B)/(A+B)

   return xcalc

def MakeFixedCh(fixed):
   fixed_plane, fixed_bar, fixed_SiPM = fixed
   return '2'+str(fixed_plane)+'00'+str(fixed_bar)+'_'+str(fixed_SiPM)

def correct_ToF(fixed, clock, xpred, cs, xref):
   
   fixed_plane, fixed_bar, fixed_SiPM = fixed
   time=clock*6.25

   c_SiPM=float(cs[0])
   ToFcorrection=abs((xpred-xref)/c_SiPM)
   if fixed_SiPM<8:
      if xpred >= xref: corrected_t = time - ToFcorrection
      else: corrected_t = time + ToFcorrection
   
   else: 
      if xpred >= xref: corrected_t = time + ToFcorrection
      else: corrected_t = time - ToFcorrection

   return (fixed_SiPM, corrected_t)

def correct_TW(iteration, fixed, rootfile, params, qdc, clock):
      
   highQDCthr=30.
   fixed_ch=MakeFixedCh(fixed)
   plane, bar, SiPM = fixed

   if qdc==-999. or clock==-999.:
      return -999.
   time=clock*6.25

   histname='dtvqdc_'+fixed_ch+'_iteration'+str(iteration)

   # Determine low QDC threshold as highest QDC bin with != 0 entries
   # rootfile e.g. /afs/cern.ch/work/a/aconsnd/Timing/rootfiles/timewalk_dists_24004_4.root
   res=DetermineLowQDCThreshold(rootfile, histname)
   if not isinstance(res, tuple):
      print('Error. result =', res)
      return res
   else:
      a, b =res
   if b==-999.:
      return -999.
   if b==-998.:
      return -998.

   lowQDCthr=b[0]

   if qdc>highQDCthr: corrected_time = time # No correction, append uncorrected times to the final list
   
   elif qdc>lowQDCthr and qdc<=highQDCthr: # Use log parameters computed for iteration i
      A, B, C = params
      TW_correction=A*ROOT.TMath.Log( (B+qdc)/(B+highQDCthr) )
      corrected_time = time + TW_correction
      # corrected_time=(SiPM, tcorr)
   
   else:
      A, B, C = params
      QDCmode, maxEntries=a
      QDC20, max20=b
      dQDC=QDCmode-QDC20
      y2=A*ROOT.TMath.Log(B+QDCmode)+C
      m, c = (y2-max20)/dQDC, y2-(y2-max20)*QDCmode/(dQDC)
      corrected_time = m * qdc + c
      # corrected_time=(SiPM, tcorr) # Investigate efficacy of linear extrapolation method between QDCmax and low threshold

   return (SiPM, corrected_time)

def correct_TW_poly(iteration, fixed, rootfile, polyparams, qdc, clock):
   highQDCthr=30.
   fixed_ch=MakeFixedCh(fixed)
   plane, bar, SiPM = fixed

   if qdc==-999. or clock==-999.:
      return -999.

   histname='dtvqdc_'+fixed_ch+'_iteration'+iteration
   res=DetermineLowQDCThreshold(rootfile, histname)
   if not isinstance(res, tuple):
      return res 

def DetermineLowQDCThreshold(rootfile, histname):
   iteration=rootfile[rootfile.find('iteration')+len('iteration')]

   infile=ROOT.TFile.Open(rootfile, 'read')
   hist=infile.Get(histname) # z.B dtvqdc_24004_4_iteration1

   detID, SiPM = histname.split('_')[1:3]
   lowQDCname='lowQDC_'+detID+'_'+SiPM

   xproj=hist.ProjectionX()
   modalqdc=xproj.GetMaximumBin()
   yproj_modalqdc=hist.ProjectionY('_py', modalqdc, modalqdc)
   meandt_modalqdc=yproj_modalqdc.GetMean()
   maxqdcentries=xproj.GetBinContent(modalqdc)

   a=(modalqdc, meandt_modalqdc)
   binwidth=xproj.GetBinWidth(1)

   b=-998.
   for qdc in range(modalqdc, 0, -1*int(binwidth)):
      tmp=hist.ProjectionY('_tmp'+str(qdc), qdc, qdc)
      tmp_mean=tmp.GetMean()
      if tmp.GetEntries()>0: b=(qdc, tmp_mean)
      else: break

   infile.Close()
   return a, b

def GoodFitFinder(plane, side):
   afswork='/afs/cern.ch/work/a/aconsnd/Timing/'
   SiPMsdict={'left':(3,4), 'right':(11, 12)}
   chi2s={}
   
   for bar in (3, 4, 5, 6):
      # if not bar in chi2s: chi2s[bar]=[]
      for SiPM in SiPMsdict[side]:
         fixed_ch = MakeFixedCh((plane, bar, SiPM))
         data=Getpolyparams_i(afswork, fixed_ch, 1)
         chi2pNDF=data[-1]
         # chi2s[bar].append((SiPM, chi2pNDF))
         chi2s[chi2pNDF]=(bar, SiPM)
            
   return chi2s
