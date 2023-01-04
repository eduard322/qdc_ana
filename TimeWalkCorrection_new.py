import ROOT, pickle, os, csv, subprocess, cProfile, io, pstats
from argparse import ArgumentParser
from pathlib import Path
import rootUtils as ut
import AnalysisFunctions as muAna
from time import perf_counter

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

parser = ArgumentParser()
parser.add_argument('-r', '--runNumber', dest='runNumber', help='runNumber', type=int, required=True)
parser.add_argument('-p', '--path', dest='path', help='path', type=str, required=True)
parser.add_argument('-t', '--track', dest='track', help='type of tracking', type=str, required=False)
parser.add_argument('-n', '--nEvents', dest='nEvents', help='number of events',type=int, required=False)
#parser.add_argument('-I', '--iterations', dest='iterations', help='number of iterations', required=False)
parser.add_argument('-s', '--startEvent', dest='startEvent', help='First event to read from eventChain', type=int, required=False)
parser.add_argument('-C', '--HTCondor', dest='HTCondor', help='int (0/1), on HTCondor?', default=0, type=int, required=False)
parser.add_argument('-A', '--AttenuationLength', dest='AttLength', help='int (0/1), running attenuation length determination?', default=0, type=int, required=False)
parser.add_argument('-m', '--mode', dest='mode', help='path', type=str, required=False)
#parser.add_argument('-F', '--fixed', dest='fixed', help='fixed channel', nargs=4, type=int, required=False)
parser.add_argument('-S', '--subsystem', dest='subsystem', help='fixed subsystem', type=int, required=False)
parser.add_argument('-P', '--plane', dest='plane', help='fixed plane', type=int, required=False)

options = parser.parse_args()

if options.HTCondor == 1: options.HTCondor=True 
else: options.HTCondor=False

"""
trackingdict={'H8':'DS', 'TI18':'Scifi'}
if options.track == None:
    options.track=trackingdict[options.path]
"""
options.track='DS' # Hard coded because Scifi acceptance is too low
afswork='/afs/cern.ch/work/a/aconsnd/Timing/'
afsuser='/afs/cern.ch/user/a/aconsnd/twfiles/'
if options.HTCondor: outpath=afswork
else: outpath=afsuser

datapaths={'H8':'/eos/experiment/sndlhc/convertedData/commissioning/TB_H8_october/', 'TI18':'/eos/experiment/sndlhc/convertedData/commissioning/TI18/'}

largeSiPMmap={0:0 ,1:1 ,3:2 ,4:3 ,6:4 ,7:5}
verticalBarDict={0:1, 1:3, 2:5, 3:6}
verticalPlanes=list(verticalBarDict.values())
v1, v2 = ROOT.TVector3(), ROOT.TVector3()

runNr = str(options.runNumber).zfill(6)
path = datapaths[options.path]
partitions = 0

if datapaths[options.path].find('eos')>0:
    path     = os.environ['EOSSHIP']+datapaths[options.path]
    # dirlist  = str( subprocess.check_output("xrdfs "+os.environ['EOSSHIP']+" ls "+datapaths[options.path],shell=True) )
# single file, before Feb'22
    if options.path == 'H8': 
        data = "sndsw_raw-"+runNr+".root"
    # if  dirlist.find(data)<0:
# check for partitions
    elif options.path == 'TI18':
        dirlist  = str( subprocess.check_output("xrdfs "+os.environ['EOSSHIP']+" ls "+datapaths[options.path]+"run_"+runNr,shell=True) )
        while 1>0:
            data = "sndsw_raw-"+str(partitions).zfill(4)
            if dirlist.find(data)>0:    partitions+=1
            else: break
else:
# check for partitions
    data = "sndsw_raw-"+runNr+".root"
    dirlist = os.listdir(datapaths[options.path])
    if not data in dirlist:
        dirlist  = os.listdir(datapaths[options.path]+"run_"+runNr)
    for x in dirlist:
        data = "data_"+ str(partitions).zfill(4)
        if x.find(data)>0:
            partitions+=1

print('\nFound '+str(partitions)+' partitions for run number '+str(runNr)+' in '+options.path+' configuration.\n')

import SndlhcGeo
# if (options.geoFile).find('../')<0: geo = SndlhcGeo.GeoInterface(path+options.geoFile)
if options.path=='H8':
    print(path)
    geo=SndlhcGeo.GeoInterface(path+'geofile_sndlhc_H6.root')
elif options.path=='TI18':
    print(path)
    geo=SndlhcGeo.GeoInterface(path+'geofile_sndlhc_TI18_V2_12July2022.root')
MuFilter = geo.modules['MuFilter']
zPos=muAna.BuildzPos(MuFilter)
Scifi = geo.modules['Scifi']
nav = ROOT.gGeoManager.GetCurrentNavigator()
print('\nGeometery constructed.\n')

A,B,locA,locB,globA,globB = ROOT.TVector3(),ROOT.TVector3(),ROOT.TVector3(),ROOT.TVector3(),ROOT.TVector3(),ROOT.TVector3()
latex = ROOT.TLatex()

if options.runNumber>0:
    eventChain = ROOT.TChain('rawConv')
    if partitions==0:
        eventChain.Add(path+'run_'+runNr+'/sndsw_raw-0000.root')
    else:
        for p in range(partitions):eventChain.Add(path+'run_'+runNr+'/sndsw_raw-'+str(p).zfill(4)+'.root')

else:
# for MC data
    f=ROOT.TFile.Open(options.fname)
    eventChain = f.cbmsim
eventChain.GetEvent(0)

run = ROOT.FairRunAna()
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
else:   eventTree = ioman.GetInTree()
# backward compatbility for early converted events
eventTree.GetEvent(0)
if eventTree.GetBranch('Digi_MuFilterHit'): eventTree.Digi_MuFilterHits = eventTree.Digi_MuFilterHit
OT = ioman.GetSink().GetOutTree()

highQDCthr=35.
systemAndPlanes = {1:2,2:5,3:4}
USbarlength=MuFilter.GetConfParF('MuFilter/UpstreamBarX')
MuFilter.GetPosition(20000, v1, v2)
xref = USbarlength/2.

cuts={'QDC/SiPM':False, 'nSiPMs':True, 'DSacceptance':False, 'slopes':False, 'yresidual':True}

def help():
    print("=="*15, "\niterative_cscint_determination(NumberOfIterations, Nev=-1). By default currently only plane 5, bar 5, SiPM index 4 is investigated.")
    print("=="*15, "\ncreateDistributionOf_cscintvalues(NumberOfIterations, fixed=(4, 4, 4)). By default currently only plane 5, bar 5, SiPM index 4 is investigated.")	
    print("=="*15,	"\nmakeDistributions(iteration, mode, Nev=-1, fixed=None, oneUShitperPlane=True, save=True)")
    print("=="*15,	"\ndetermine_cscint(iteration, fixed=None, save=True)")
    print("=="*15, "\ndetermineLogParams(iteration, fixed=None, save=True)")
    print('\n', "=="*15)

def makeDistributions(iteration, mode, Nev=None, save=True):

    t1_start=perf_counter()
    if not any( (mode=='ToF', mode=='TW', mode=='zeroth') ):
        print('Incorrect option for \'mode\' parameter passed to function and/or fixed==None.')
        from sys import exit
        exit()

    if os.path.exists(afswork+'SelectionCriteria/SelectionCriteria_run'+runNr+'.root'): cutNr=str(runNr)
    else: 
        if options.path=='TI18': cutNr='004208' # Might not work right now
        elif options.path=='H8': cutNr='000046' 
                
    IsCut=lambda x : x[1]==True
    cutdists2get=[i[0] for i in list(filter(IsCut, list(cuts.items())))]
    cutdists=muAna.GetCutDistributions(cutNr, cutdists2get)
    #return cutdists
    counters={'tracks':0, 'inresidual':0,'slopes':0, 'goodDST0':0}
	
    if mode == 'ToF':
        dtvqdc_dists={}
        
    elif mode == 'TW':
        dtvxpred_dists={}
        correctionfunction=lambda ps, qdc : 1/sum( [ ps[i]*qdc**i for i in range(len(ps)) ] )


    elif mode == 'zeroth':
        dtvxpred_dists={}
        
    if Nev==None:
        start=options.startEvent
        nEvents=options.nEvents
    else:
        start=0
        nEvents=int(Nev)

    deciles=[i/10 for i in range(11)]
    
    print('start event: {X}, number of events: {Y}'.format(X=start, Y=nEvents))

    for i in range(start, start+nEvents):
        eventTree.GetEvent(i)

        if ((i-start)/nEvents) in deciles: 
            progstr="Progress: {X}/{Y}, {Z}%".format( X=i, Y=(start+nEvents), Z=int(100*(i-start)/nEvents) )
            print("="*len(progstr)+'\n')
            print(progstr, '\n')

        if not muAna.OneHitPerUS(eventTree.Digi_MuFilterHits): continue

	######## Using SciFi tracks from TI18 data #######
        if options.track == 'DS':
            rc = trackTask.ExecuteTask('DS')
        else:
            rc = trackTask.ExecuteTask('Scifi')
        if not OT.Reco_MuonTracks.GetEntries()==1:
            continue
        theTrack=OT.Reco_MuonTracks[0]
        if not theTrack.getFitStatus().isFitConverged() and options.track!='DS':
            continue
        state=theTrack.getFittedState()
        pos=state.getPos()
        mom=state.getMom()
        slopeX, slopeY = mom.x()/mom.z(), mom.y()/mom.z()
        if cuts['slopes']:
            if abs(slopeX)>0.1 or abs(slopeY)>0.1: continue 
        counters['slopes']+=1
    	
        ###### Get T0 estimate from DS horizontal SiPMs used in single track event selection
        DST0=muAna.GetDSH_average(eventTree.Digi_MuFilterHits)
        if DST0 == -999.:
            continue # Multiple hits assigned to one bad 
        if DST0 == -998.: 			continue # Only 1 SiPM in a horizontal bar is firing	
        if DST0 == -6237.5:			continue # A bad value permeating the whole data set of run46
        DST0*=6.25
        if DST0 < 13.75 or DST0 > 15.45: continue # DST0 limit is run invariant
   
        ######## Determine average time measured by all fired bars in DS track fit. Store slope_x and y for histogram
        for j, aHit in enumerate(eventTree.Digi_MuFilterHits):
            if not aHit.isValid(): continue
            nSides=aHit.GetnSides()
            nSiPMs=aHit.GetnSiPMs()
        
            ######## cuts on nFired SiPMs and QDC/SiPM
            nFiredSiPMs_left, nFiredSiPMs_right = muAna.GetnFiredSiPMs(aHit)
            if cuts['QDC/SiPM']: # set to False
                if TotalQDC/TotalFiredSiPMs<15: continue
            if cuts['nSiPMs']: # set to True
                if nFiredSiPMs_left<5 or nFiredSiPMs_right<5: continue     # Minimise noise hits contribution
	    ########

            ########### dy-cut
            detID=aHit.GetDetectorID()
            if detID == '24008': print('here')
            subsystem, plane, bar = muAna.parseDetID(detID)
            tagged=False
            key=subsystem*10+plane
            
            z=zPos[key]
            lam = (z-pos.z())/mom.z()
            Ex = ROOT.TVector3(pos.x()+lam*mom.x(), pos.y()+lam*mom.y(), pos.z()+lam*mom.z())
            MuFilter.GetPosition(detID, v1, v2)
            dy=Ex.y() - v1.y()
            dy_min, dy_max = cutdists[str(key)+'_yresidual'].GetMean()-cutdists[str(key)+'_yresidual'].GetStdDev(), cutdists[str(key)+'_yresidual'].GetMean()+cutdists[str(key)+'_yresidual'].GetStdDev()
            if dy_min ==0 or dy_max==0: 
                dy_min=-3.
                dy_max=3.
            # distance between extrapolated-y and fired bar y-midpoint
            if subsystem==1 or subsystem==2: # Cut on xEx - bar midpoint from extrapolationCutInvestigation.py
                if dy>dy_min and dy<dy_max: tagged=True
            if not tagged: continue
            ########

            ######## x-pred is the extrapolated x-position in the bar's coordinate system
            xpred=v1.x() - Ex.x()
            xpred_R=-xpred
            
            # No point in cutting on xpred here. Just only perform the fit in the healthy range
            """ 
            if subsystem==1:
                if xpred<22 or xpred>40: continue
            if subsystem==2:
                if xpred<5 or xpred>75.: continue
            ########
            """
            ######## Arrays of all SiPM TDCs and QDCs
            channels_t = aHit.GetAllTimes()
            channels_qdc = aHit.GetAllSignals()
            for channel in channels_t:
                SiPM, clock=channel
                qdc=muAna.GetChannelVal(SiPM, channels_qdc)
                if qdc==-999.: continue # Only use hits with both QDC and a clock. Should never flag True
                fixed=(subsystem, plane, bar, SiPM)
                fixed_ch = muAna.MakeFixedCh(fixed)
                
                if mode == 'zeroth':
                    correctedtime=clock*6.25				
                    t_rel = DST0 - correctedtime
                    name=f'dtvxpred_{fixed_ch}_iteration{iteration}'
                    if name not in dtvxpred_dists:								
                        title = f'Uncorrected T_{{0}}^{{DS}}-t_{SiPM} v predicted x-position;x_{{predicted}} [cm];T_{{0}}^{{DS}}-t_{SiPM} [ns]'
                        if subsystem == 1: xmax=60.
                        elif subsystem == 2: xmax=100.
                        dtvxpred_dists[name] = ROOT.TH2F(name, title, 55, -10., xmax, 160, -20., 20.)
                    if 'DST0' not in dtvxpred_dists:
                        dtvxpred_dists['DST0']=ROOT.TH1F('DST0', 'DSH average;DSH average [ns];Counts',100,0.,25.)
                    name2=f'tSiPM_{fixed_ch}_iteration{iteration}'
                    if name2 not in dtvxpred_dists:
                        dtvxpred_dists[name2]=ROOT.TH1F(name2, f'Measured SiPM time {fixed_ch}; SiPM time [ns]; Counts', 100, 0., 25.)
                    name3=f'attlen_{fixed_ch}_iteration{iteration}'
                    if name3 not in dtvxpred_dists:
                        dtvxpred_dists[name3]=ROOT.TH2F(name3, f'Predicted position against QDC_{SiPM} {fixed_ch}; x_{{predicted}} [cm]; QDC_{SiPM} [a.u]', 55, -10., xmax, 100, 0., 25.)
                    dtvxpred_dists[name].Fill(xpred, t_rel)
                    dtvxpred_dists['DST0'].Fill(DST0)
                    dtvxpred_dists[name2].Fill(correctedtime)
                    dtvxpred_dists[name3].Fill(xpred, qdc) 

                elif mode=='ToF':
                    cdata=muAna.Getcscint_i(runNr, fixed_ch, iteration-1)
                    if cdata==-999.:continue
                    correctedtime = muAna.correct_ToF(SiPM, clock, xpred, cdata, xref)[1] # correct_ToF(detID, times, x_pred, cs, fixed)
                     
                    name=f'dtvqdc_{fixed_ch}_iteration{iteration}'
                    if not name in dtvqdc_dists:
                        title = f'T_{{0}}^{{DS}}-t_{SiPM} v QDC_{SiPM}. Iteration {iteration};QDC_{SiPM} [a.u];T_{{0}}^{{DS}}-t_{SiPM} [ns]'
                        dtvqdc_dists[name] = ROOT.TH2F(name, title,  200, 0., 200., 160, -20., 20.)
                    #if 'DST0' not in dtvqdc_dists:
                        #dtvqdc_dists['DST0']=ROOT.TH1F('DST0', 'DSH average;DSH average [ns];Counts',100,0.,25.)

                    t_rel = DST0 - correctedtime
                    dtvqdc_dists[name].Fill(qdc, t_rel)
                    #dtvqdc_dists['DST0'].Fill(DST0)

                elif mode=='TW':
                    time=clock*6.25
                    polyparams=muAna.GetPolyParams(runNr, fixed_ch, iteration)
                    
                    chi2=muAna.Getchi2pNDF(runNr, fixed_ch, iteration-1)
                    cdata=muAna.Getcscint_i(runNr, fixed_ch, iteration-1)
                    if polyparams == -999. or cdata == -999.: continue
                    if chi2>2: continue
                    if len(polyparams)==4:
                        to_apply=correctionfunction(polyparams[0:-1], qdc)	
                        correctedtime = time+to_apply+polyparams[-1]
                    elif len(polyparams)==3: 
                        to_apply=correctionfunction(polyparams, qdc)
                        correctedtime = time+to_apply
                    tmp_correctedtime2=muAna.correct_ToF(SiPM, clock, xpred, cdata, xref)[1]
                    correctedtime2 = tmp_correctedtime2+to_apply
                    t_rel = DST0 - correctedtime
                    t2_rel = DST0 - correctedtime2

                    name='dtvxpred_'+fixed_ch+'_iteration'+str(iteration)
                    if not name in dtvxpred_dists:
                        title = f'T_{{0}}^{{DS}}-t_{SiPM} v predicted x-position. Iteration {iteration};x_{{predicted}} [cm];T_{{0}}^{{DS}}-t_{{SiPM}} [ns]'
                        if subsystem == 1: xmax=60.
                        elif subsystem == 2: xmax=100.
                        dtvxpred_dists[name] = ROOT.TH2F(name, title, 55, 0., xmax, 160, -20., 20.)
                    
                    name2='dtcorrectedvqdc_'+fixed_ch+'_iteration'+str(iteration)
                    if not name2 in dtvxpred_dists:
                        title2 = f'T_{{0}}^{{DS}}-t^{{corrected}}_{SiPM} v QDC_{SiPM}. Iteration {iteration};QDC_{SiPM} [a.u];T_{{0}}^{{DS}}-t^{{corrected}}_{SiPM} [ns]'
                        dtvxpred_dists[name2] = ROOT.TH2F(name2, title2,  150, 0., 150., 160, -40., 40.)
                    
                    #if 'DST0' not in dtvxpred_dists:
                        #dtvxpred_dists['DST0']=ROOT.TH1F('DST0', 'DSH average;DSH average [ns];Counts',100,0.,25.)
                         
                    dtvxpred_dists[name].Fill(xpred, t_rel)
                    dtvxpred_dists[name2].Fill(qdc, t2_rel)
                    #dtvxpred_dists['DST0'].Fill(DST0)
                        
        theTrack.Delete()
	# return myhist
    if mode=='zeroth' or mode=='TW':
        hists=dtvxpred_dists
    elif mode=='ToF':
        hists=dtvqdc_dists
    #return hists
    print(len(hists))
    if not save: return hists
    else:
        for histname in hists:
            if len(histname.split('_'))==4:
        
                histkey, detID, SiPM, iterationNumber=histname.split('_')
                fixed_ch='_'.join((detID, SiPM))
                hist=hists[histname]
                outpath=f'{afswork}/splitfiles/run{runNr}/{fixed_ch}/'
                path_obj=Path(outpath)
                path_obj.mkdir(parents=True, exist_ok=True)
    
                outfile=f'timewalk_{detID}_{SiPM}_{start}.root'
                f=ROOT.TFile.Open(outpath+outfile, 'update')
                f.WriteObject(hist, hist.GetName(), 'kOverwrite')
                f.WriteObject(hists['DST0'], 'DST0', 'kOverwrite')
                f.Close()
    print(f'Files saved to {afswork}')
    
    t1_stop=perf_counter()
    print(f'Duration: {t1_stop-t1_start} s')

small=[2, 5, 10, 13]

def determinePolyParams(iteration, fixed, mode='cs2', save=True):

    ROOT.gROOT.SetBatch(True)
    t1_start=perf_counter()

    polymodes={'cs2':'1./( [p0] + [p1]*x +[p2]*x*x )', 'cs3':'1./( [p0] + [p1]*x + [p2]*x*x + [p3]*x*x*x)','6d':'pol6'}
    if mode not in polymodes:
        print('Mode must be either \'cs\', \'6d\'')
        from sys import exit
        exit()

    if not os.path.exists(afswork+'Polyparams/run'+str(runNr)): os.makedirs(afswork+'Polyparams/run'+str(runNr))	
    if not os.path.exists(afswork+'chi2s/run'+str(runNr)): os.makedirs(afswork+'chi2s/run'+str(runNr))	
        
    if path=='TI18': 
        nplanes=range(5)
        nbars=range(10)
        SiPMs=(0,1,3,4,6,7,8,9,11,12,14,15)
    elif path == 'H8':
        nplanes=range(5)
        nbars=range(10)
        SiPMs=(0,1,3,4,6,7,8,9,11,12,14,15)

    fixed_subsystem,fixed_plane,fixed_bar,fixed_SiPM = fixed
    fixed_ch=muAna.MakeFixedCh(fixed)  
    detID=int(fixed_ch.split('_')[0])
    
    xvqdc={}
    title='SiPM '+str(fixed_SiPM)
    canvname='twfitpoly_'+fixed_ch+'_'+str(iteration)
    xvqdc[canvname]=ROOT.TCanvas(canvname, title, 1200, 1200)
    hist_filename=afswork+'rootfiles/run'+str(runNr)+'/timewalk_'+fixed_ch+'.root'				
    hist_file = ROOT.TFile.Open(hist_filename, 'UPDATE')
    print('Number of keys in file:', len(hist_file.GetListOfKeys()))
    histname='dtvqdc_'+fixed_ch+'_iteration'+str(iteration)
    if not hasattr(hist_file, histname):
        print(hist_file, 'has no object', histname)
        return 0        
    hist=hist_file.Get(histname).Clone()
    print('Opened file', hist_file) 
    print('Retrieved histogram ', histname)
	
    #a, b = muAna.DetermineLowQDCThreshold(afswork+'rootfiles/run'+str(runNr)+'/timewalk_dists_'+fixed_ch+'_poly.root', 'dtvqdc_'+fixed_ch+'_iteration'+str(iteration))
    lowX=2
    f1=ROOT.TF1('fitfunction', polymodes[mode], lowX, 60)	
    f2=ROOT.TF1('fitfunction+offset', '1./([p0] + [p1]*x + [p2]*x*x) + [p3]', lowX, 60)
        
    if mode == 'cs2':
        sameplane=True # Use a constant offset parameter IF we use starting values from another plane
        if int(fixed_SiPM)<8:side='left'
        else: side='right'
        sameplane_chi2s=muAna.GoodFitFinder(runNr, int(fixed_subsystem), int(fixed_plane), side)
        print(f'len sameplane_chi2s:{len(sameplane_chi2s)}') 
        polyparamrun=runNr
        if len(sameplane_chi2s)<=40 and fixed_subsystem==2: 
            polyparamrun='000046'
            print('Using run 46')
            sameplane_chi2s=muAna.GoodFitFinder(polyparamrun, int(fixed_subsystem), int(fixed_plane), side)
        if len(sameplane_chi2s)<=38 and fixed_subsystem==1: 
            #polyparamrun='000046'
            #print('Using run 46')
            #sameplane_chi2s=muAna.GoodFitFinder(polyparamrun, int(fixed_subsystem), int(fixed_plane), side)
            first_sameplane=True
           
        else: 
            polyparamrun=runNr
            first_sameplane=False
        
        #first_sameplane=False

        startingparams=-999.  
        if not first_sameplane:
            #### Check if a SiPM on the same bar end has a chi2/NDF less than 1.5
            val=min(list(sameplane_chi2s.keys()))
            if val < 1.5:
                starting_bar, starting_SiPM = sameplane_chi2s.pop(val) # Pop best value!
                print("Taking the value from:", fixed_plane, starting_bar, starting_SiPM)
                startingchannel=muAna.MakeFixedCh((fixed_subsystem,fixed_plane, starting_bar, starting_SiPM))
                print('Starting parameters from', startingchannel)
                tmp=muAna.GetPolyParams(polyparamrun, startingchannel, iteration)
                #if len(tmp)==6: startingparams=tmp
                startingparams=tmp
                # else: startingparams==-999.
                #### Make a tuple of the two closest planes to the channel being looked at

            if startingparams==-999.: 
                    best_adjacent=muAna.adjacentplaneschi2s(polyparamrun, int(fixed_subsystem), int(fixed_plane), side)
                    best_effort=min(best_adjacent.keys())
                    #### If an adjacent plane contains a value less than 1.5, use this

                    if best_effort < 1.5:
                        plane2use=best_adjacent[best_effort]
                        adjacentPlane_chi2s=muAna.GoodFitFinder(polyparamrun,int(fixed_subsystem), plane2use, side) # chi2/NDF values for adjacent plane with lowest chi2/NDF
                        tmp=min(list(adjacentPlane_chi2s.keys()))
                        starting_bar, starting_SiPM = adjacentPlane_chi2s[tmp]
                        sameplane=False
                        print("sameplane false 1")
                        startingchannel=muAna.MakeFixedCh((fixed_subsystem,plane2use, starting_bar, starting_SiPM))
                        print('Starting parameters from', startingchannel)

                        tmp=muAna.GetPolyParams(polyparamrun, startingchannel, iteration)
                        # if len(tmp)==8: startingparams=tmp
                        startingparams=tmp
                        # else: startingparams=-999.
                                    
                        #### Check what is lower, the second best chi2/NDF from the same plane, or the best chi2/NDF of the 2 adjacent planes
                    else:
                        val2=min(sameplane_chi2s.keys()) # Best value is popped so the new min value is the second best result.
                        if val2 <= 1.5:
                            starting_bar, starting_SiPM = sameplane_chi2s.pop(val2)
                            # Taking the 2nd best value from same plane
                            startingchannel=muAna.MakeFixedCh((fixed_subsystem,fixed_plane, starting_bar, starting_SiPM))
                            print('Starting parameters from', startingchannel)
                            tmp=muAna.GetPolyParams(polyparamrun, startingchannel, iteration)
                            #if len(tmp)==6: startingparams=tmp
                            startingparams=tmp
                            # else: startingparams=-999.

                        elif best_effort <= 1.5:
                            plane2use=best_adjacent[best_effort]
                            adjacentPlane_chi2s=muAna.GoodFitFinder(polyparamrun, int(fixed_subsystem),plane2use, side)
                            tmp=min(adjacentPlane_chi2s.keys())
                            starting_bar, starting_SiPM = adjacentPlane_chi2s[tmp]
                            sameplane=False
                            print("sameplane false 2")
                            # Taking the best value from adjacent plane
                            startingchannel=muAna.MakeFixedCh((fixed_subsystem,plane2use, starting_bar, starting_SiPM))
                            print('Starting parameters from', startingchannel)
                            startingparams=muAna.GetPolyParams(polyparamrun, startingchannel, iteration)

                        else: 
                            # All values are shite. Taking from 24004_4.				
                            startingchannel=muAna.MakeFixedCh((2,4,4,4))
                            print('Starting parameters from', startingchannel)
                            tmp=muAna.GetPolyParams(polyparamrun, startingchannel, iteration)
                            if tmp != -999.:
                                if (fixed_plane != 4 and len(tmp)==8) or (fixed_plane==4 and len(tmp)==6): startingparams=tmp
                                if fixed_plane != 4: 
                                    sameplane=False
                                    print('sameplane false 3') 
                                # else: startingparams=-999.                            
            print('Starting params: {X}'.format(X=startingparams))
       
        #### Get histogram to make graph
        title='SiPM '+str(fixed_SiPM)
        canvname='twfitpoly_'+fixed_ch+'_'+str(iteration)
        xvqdc[canvname]=ROOT.TCanvas(canvname, title, 1200, 1200)

        hist_filename=afswork+'rootfiles/run'+str(runNr)+'/timewalk_'+fixed_ch+'.root'	    	
        hist_file = ROOT.TFile.Open(hist_filename, 'UPDATE')

        histname='dtvqdc_'+fixed_ch+'_iteration'+str(iteration)
        if not hasattr(hist_file, histname):
            print(hist_file, 'has no object', histname)
            return 0 
        hist=hist_file.Get(histname).Clone()
        if hist.GetEntries()<1000.:
            print(f'histogram has {hist.GetEntries()} entries. Not attempting fit')
            return 
        print('Opened file', hist_file)

        detID, SiPM = histname.split('_')[1:3]
        subsystem, plane, bar = muAna.parseDetID(int(detID))

        xvqdc[canvname].cd()
        hist.Draw('colz')

        # Making graph
        xvqdc['g_'+hist.GetName()]=ROOT.TGraphErrors()
        xvqdc['g_'+hist.GetName()].SetName('GraphFromtvqdc')
        g=xvqdc['g_'+hist.GetName()]
        xproj=hist.ProjectionX('tmpx')
        modalqdc=xproj.GetMaximumBin()
        binwidth=xproj.GetBinWidth(1)
        pointcounter=0

        for nx in range(lowX, int(highQDCthr)):
            tmp = hist.ProjectionY('tmp', nx, nx)
            if xproj.GetBinContent(nx) == 0: continue
            X = xproj.GetBinCenter(nx)
            dt=tmp.GetMean()
            std_error=tmp.GetStdDev()/ROOT.TMath.Sqrt(tmp.GetEntries())
            g.SetPoint(pointcounter, X, dt)
            g.SetPointError(pointcounter, 0.5*binwidth, std_error) 
            pointcounter+=1
        highX=min(g.GetPointX(g.GetN()-1), 35)
        g.SetLineColor(ROOT.kRed)
        g.SetLineWidth(2)
        g.Draw('same *')

        #### Make dictionaries to store parameters and chi2 + Ndf values for each fit
        chi2pNDF_dictionary={}
        chi2_dictionary={}

        #### First try with no starting values
        if sameplane: rc = g.Fit('fitfunction', 'S N Q', '', lowX, highX)
        else: rc = g.Fit('fitfunction+offset', 'S N Q', '', lowX, highX)
        result=rc.Get()
        if result.Ndf()!=0:
            zeroth_chi2pNDF=result.Chi2()/result.Ndf()
            chi2pNDF_dictionary[zeroth_chi2pNDF]=[[result.Parameter(i), result.ParError(i)] for i in range(len(result.Parameters()))]
            chi2_dictionary[zeroth_chi2pNDF]=[result.Chi2(), result.Ndf()]
            print("0th attempt chi2 per NDF:", zeroth_chi2pNDF)
        if not first_sameplane: print('Starting parameters:{X}'.format(X=startingparams))
        #### Set starting values if we have any 
        if startingparams != -999.: 
            for i in range(3): 
                # if sameplane: f1.SetParameter(i, startingparams[i*2])
                # else: f2.SetParameter(i, startingparams[i*2])		
                if sameplane: f1.SetParameter(i, startingparams[i])
                else: f2.SetParameter(i, startingparams[i])	
            if sameplane: rc = g.Fit('fitfunction', 'S N Q', '', lowX, highX)
            else: rc = g.Fit('fitfunction+offset', 'S N Q', '', lowX, highX)	

            result=rc.Get()
            if result.Ndf()!=0: 
                chi2pNDF=result.Chi2()/result.Ndf()
                print("1st attempt chi2 per NDF:", chi2pNDF)
            
                chi2pNDF_dictionary[chi2pNDF]=[[result.Parameter(i), result.ParError(i)] for i in range(len(result.Parameters()))]
                chi2_dictionary[chi2pNDF]=[result.Chi2(), result.Ndf()]

        if len(chi2pNDF_dictionary)!=0: best_chi2pNDF = min(list(chi2pNDF_dictionary.keys()))
        else: best_chi2pNDF = 5.

        #### If chi2pNDF greater than 1.5: try with all fits from the same plane
        sameplanetesting=False
        if best_chi2pNDF > 1.5: 
            sameplanetesting=True
            while len(sameplane_chi2s)>0:
                tmp_chi2 = min(list(sameplane_chi2s.keys()))
                tmp_bar, tmp_SiPM = sameplane_chi2s.pop(tmp_chi2)
                test_ch=muAna.MakeFixedCh((fixed_subsystem,fixed_plane, tmp_bar, tmp_SiPM))
                test_params=muAna.GetPolyParams(polyparamrun, test_ch, iteration)
                if test_params!=-999.:
                    # if len(test_params) == 6:
                    # for i in range(3):f1.SetParameter(i, test_params[i*2])
                    for i in range(3):f1.SetParameter(i, test_params[i])
                    rc=g.Fit('fitfunction', 'S N Q', '', lowX, highX)
                    res=rc.Get()
                    if res.Ndf()>0:
                        testing_chi2=res.Chi2()/res.Ndf()
                            # if testing_chi2 < min((chi2pNDF, zeroth_chi2pNDF)):
                        if testing_chi2 < zeroth_chi2pNDF:
                            chi2pNDF_dictionary[testing_chi2]=[[res.Parameter(i), res.ParError(i)] for i in range(len(res.Parameters()))]
                            chi2_dictionary[testing_chi2]=[res.Chi2(), res.Ndf()]

        if len(chi2pNDF_dictionary)!=0: best_chi2pNDF = min(list(chi2pNDF_dictionary.keys()))
        else: best_chi2pNDF=5.

        # If best chi2/ndf is still above 1.5 try with all fits in the US system from that run. 
        tempcounter=0
        if best_chi2pNDF > 1.5:
            remainingplanes=[i for i in range(5)]
            remainingplanes.pop(fixed_plane)
            for test_plane in remainingplanes:
                for test_bar in range(10):
                    for test_SiPM in (0,1,3,4,6,7,8,9,11,12,14,15):
                        test_ch=muAna.MakeFixedCh((fixed_subsystem,test_plane, test_bar, test_SiPM))
                        if test_ch == '21005_15': continue
                        test_params=muAna.GetPolyParams(polyparamrun, test_ch, iteration)
                        if test_params==-999.: continue
                        # if len(test_params) != 8: continue
                        # for i in range(3):  f2.SetParameter(i, test_params[i*2])
                        for i in range(3):  f2.SetParameter(i, test_params[i])
                        rc=g.Fit('fitfunction+offset', 'S N Q', '', lowX, highX)
                        res=rc.Get()
                        if res.Ndf()>0:
                            testing_chi2=res.Chi2()/res.Ndf()
                            if testing_chi2 < best_chi2pNDF: 
                                chi2pNDF_dictionary[testing_chi2]=[[res.Parameter(i), res.ParError(i)] for i in range(len(res.Parameters()))]
                                chi2_dictionary[testing_chi2]=[res.Chi2(), res.Ndf()]

                                # print("Best parameters from {X}".format(X=test_ch))
                                #print(len(res.Parameters()))
                                tempcounter+=1
                                if tempcounter==1: print('sameplane false last')	
                                if testing_chi2<1.5: break
        
        # return chi2pNDF_dictionary
        # if len(chi2pNDF_dictionary)!=0: best_chi2pNDF = min(list(chi2pNDF_dictionary.keys()))
        # else: best_chi2pNDF=5.
        if len(chi2pNDF_dictionary)==0:    
            print(f'All attempts at fit have failed. Channel {fixed_ch} is shit.')
            return 0

        best_params=chi2pNDF_dictionary[best_chi2pNDF]
        if len(best_params)==4: sameplane=False
        else: sameplane=True

        #### Testing with sequential removal of lowest QDC bin
        if best_chi2pNDF > 1.5:
            print('\n'+'***'*6+' Sequential removal of lowest QDC bin and re-evaluate chi2/NDF '+'***'*6+'\n')
            print('Current chi2/ndf:', best_chi2pNDF, '\n')
            print("best_params:", best_params, '\n')
            pointRemovalChi2s={}
            for i in range(3):
                # Remove lowest QDC point and retry with same starting values.
                if sameplane:
                    for j in range(3):	f1.SetParameter(j, best_params[j][0])
                    rc = g.Fit('fitfunction', 'S N Q', '', g.GetPointX(i+1), highX)
                elif not sameplane:
                    for j in range(3):	f2.SetParameter(j, best_params[j][0])
                    rc = g.Fit('fitfunction+offset', 'S N Q', '', g.GetPointX(i+1), highX)
   
                res=rc.Get()
                if res.Ndf()==0: continue
                testing_chi2=res.Chi2()/res.Ndf()
                print('chi2/NDF after removing point '+str(i)+': '+str(testing_chi2))
                if testing_chi2 < best_chi2pNDF:
                    pointRemovalChi2s[testing_chi2]=i
                    chi2pNDF_dictionary[testing_chi2]=[[res.Parameter(k), res.ParError(k)] for k in range(len(res.Parameters()))]
                    chi2_dictionary[testing_chi2]=[res.Chi2(), res.Ndf()]
                if testing_chi2 < 1.5: break

        best_chi2pNDF = min(list(chi2pNDF_dictionary.keys()))	
        best_params = chi2pNDF_dictionary[best_chi2pNDF]
        print('+++++'*10, '\n')
        print(best_chi2pNDF)
        print('+++++'*10, '\n')

        if best_chi2pNDF > 1.5:
            print('\n'+'***'*6+' Sequential removal of highest QDC bin and re-evaluate chi2/NDF '+'***'*6+'\n')
            print('Current chi2/ndf:', best_chi2pNDF, '\n')
            print("best_params:", best_params, '\n')
            if len(pointRemovalChi2s)==0: lowtmp=lowX
            else: lowtmp=min(pointRemovalChi2s.keys())
            for i in range(15):
                hightmp=highX-i
                if sameplane:
                    for j in range(3): f1.SetParameter(j, best_params[j][0])
                    rc=g.Fit('fitfunction', 'S N Q', '', lowtmp, hightmp)
                elif not sameplane:
                    for j in range(3): f2.SetParameter(j, best_params[j][0])   
                    rc=g.Fit('fitfunction+offset', 'S N Q', '', lowtmp, hightmp)
                
                res=rc.Get()
                if res.Ndf()==0: continue

                testing_chi2=res.Chi2()/res.Ndf()
                if testing_chi2 < best_chi2pNDF:
                    chi2pNDF_dictionary[testing_chi2]=[ [res.Parameter(k), res.ParError(k)] for k in range(len(res.Parameters())) ]
                    chi2_dictionary[testing_chi2]=[res.Chi2(), res.Ndf()]
                    print('chi2/NDF with {X} as the upper fit limit: {Y}'.format(X=hightmp, Y=testing_chi2))
                    if testing_chi2 < 1.5: break

        best_chi2pNDF = min(list(chi2pNDF_dictionary.keys()))	
        best_params = chi2pNDF_dictionary[best_chi2pNDF]
        #chi2_info=min(list(chi2_dictionary[best_chi2pNDF]))
        chi2_info=chi2_dictionary[best_chi2pNDF]
        print(f'chi2 info: {chi2_info}')
        print('+++++'*10, '\n')
        print(f'best chi2/NDF: {best_chi2pNDF}')
        print('+++++'*10, '\n')

        if sameplane:
            for i in range(f1.GetNpar()): f1.SetParameter(i, best_params[i][0])
            f1.Draw('same')
        else:
            for i in range(f2.GetNpar()): f2.SetParameter(i, best_params[i][0])
            f2.Draw('same')

        if save:
            writeoptions={1:'w', 2:'a'}
            opt=writeoptions[iteration]
            hist_file.WriteObject(xvqdc[canvname], canvname, 'kOverwrite')
            hist_file.Close()
            print('Polyfit of ', histname, ' written to', hist_filename)

            polyline=[item for sublist in best_params for item in sublist]
            polyline.insert(0,iteration)
            chi2_info.insert(0, iteration)
            #chi2_info.extend(best_chi2pNDF)
            chi2line=chi2_info
            #line.append(best_chi2pNDF)

            outfile=afswork+'Polyparams/run'+str(runNr)+'/polyparams_'+fixed_ch+'.csv'
            with open(outfile, opt) as polyparams_out:
                writer=csv.writer(polyparams_out)
                #print(polyline)
                writer.writerow(polyline)
                print('Polyparams written to', outfile)
            outfile=afswork+'chi2s/run'+str(runNr)+'/chi2s_'+fixed_ch+'.csv'
            with open(outfile, opt) as chi2s_out:
                writer=csv.writer(chi2s_out)
                #print(chi2line)
                writer.writerow(chi2line)
                print('chi2 written to', outfile)

        else:
            if sameplane: return g, f1
            else: return g,f2

def determine_cscint(iteration, fixed, save=True):
    
    ROOT.gROOT.SetBatch(True)
    xvt={}

    fixed_ch=muAna.MakeFixedCh(fixed)
    hist_filename=afswork+'rootfiles/run'+str(runNr)+'/timewalk_'+fixed_ch+'.root'
    fixed_subsystem, fixed_plane, fixed_bar, fixed_SiPM = fixed

    title='Linearly fitted T_0^{DS} - t_{SiPM} v QDC_{SiPM}'
    canvname='linearfit_'+fixed_ch+'_'+str(iteration)
    xvt[canvname]=ROOT.TCanvas(canvname, title, 1200, 1200)
    
    if not os.path.exists(afswork+'cscintvalues/run'+str(runNr)): os.makedirs(afswork+'cscintvalues/run'+str(runNr))

    hist_file = ROOT.TFile.Open(hist_filename, 'update')
    print('Number of keys in file at the start of this step:', len(hist_file.GetListOfKeys()))
    histname='dtvxpred_'+fixed_ch+'_iteration'+str(iteration)
    if not hasattr(hist_file, histname):
        print(hist_file, 'has no object', histname)
        return 0
    hist=hist_file.Get(histname)

    detID, SiPM = histname.split('_')[1:3]
    s, p, b = muAna.parseDetID(int(detID))
    key=detID+'_'+str(SiPM)

    xvt[histname]=hist.Clone()
    hist=xvt[histname]
    xvt[canvname].cd()
    hist.Draw('colz')
    xvt['g'+histname] = ROOT.TGraphErrors()
    g=xvt['g'+histname]
    xpred_proj=hist.ProjectionX() # TH1D of xpred 
    binwidth=xpred_proj.GetBinWidth(0)
    for xpred_bin in range(5, 75):
        tmp = hist.ProjectionY('tmp', xpred_bin, xpred_bin) # TH1D of (2x-L) for the given dt bin
        if tmp.GetEntries() < 20: continue
        x = xpred_proj.GetBinCenter(xpred_bin)
        dt = tmp.GetMean()
        std_error=tmp.GetStdDev()/ROOT.TMath.Sqrt(tmp.GetEntries())
        g.SetPoint(xpred_bin-1, x, dt) # SetPoint(i, x, y) - ith point with value (x, y)
        g.SetPointError(xpred_bin-1, binwidth*0.5, std_error)

    print('Number of entries in graph:', g.GetN())
    if g.GetN()<10:
        print('Too few points on graph for accurate fit')
        return 0
    g.SetLineColor(ROOT.kRed)
    g.SetLineWidth(2)
    g.SetTitle(hist.GetTitle()) # Keep title of histogram from unfitted version
    g.Draw('same *')
    if options.track=='Scifi' and s==1: xmax=35.
    elif options.track=='Scifi' and s==2: xmax=40.
    elif options.track=='DS' and s==1: xmax=35.
    elif options.track=='DS' and s==2: xmax=60.
    fitfunction=ROOT.TF1('fitfunction', '[0]+[1]*x', 5., xmax)
    if s<3 and int(SiPM)<8: fitfunction.SetParLimits(1, -1/5, -1/20)
    elif s<3 and int(SiPM)>7: fitfunction.SetParLimits(1, 1/20, 1/5)

    rc=g.Fit('fitfunction','S', '', 5., xmax)
    result=rc.Get()
    if result.Ndf() == 0: return 0
    chi2pNDF=result.Chi2()/result.Ndf()
    try:
        m=result.Parameter(1)
        intercept=result.Parameter(0)
    except ReferenceError:
        print('Cannot extract fit result')
        from sys import exit 
        exit()
    m_err=result.ParError(1)
    intercept=result.Parameter(1)
    intercept_err=result.ParError(0)
 
    if int(SiPM)<8: c, c_err = -1/m, (m_err/m)
    else: c, c_err = 1/m, (m_err/m)
    c_err=abs(c*m_err/m)
    print(key, " c_scint +/- err [cm/ns]: ", c, c_err )
    
    fitfunction.SetParameter(0, m)
    fitfunction.SetParameter(1, c)
    fitfunction.Draw('same')
    
    if save:
        hist_file.WriteObject(xvt[canvname], canvname, 'kOverwrite')
        hist_file.Close()
        print('Fitted dt v xpred written to', hist_filename)

        outfile=afswork+'cscintvalues/run'+str(runNr)+'/cscint_'+fixed_ch+'.csv'
        with open(outfile, 'a') as cs_out:
            writer=csv.writer(cs_out)
            line = (iteration, c, c_err, intercept, intercept_err, chi2pNDF)
            writer.writerow(line)
            print('c_scint value written to', outfile)

    elif not save: return xvt[canvname]

"""
def xQDC_dependence(Nev=-1, iteration=0, xwidth=2, fixed=tuple(options.fixed), oneUShitperPlane=True, overlayhists=False, QDCmodegraphs=False, save=True):
	
	modes={0: 'attenuationlength', 1:'overlayhists', 2:'QDCmodegraphs'}


	dy_hists={}

	xQDC={}
	graphs={}
	canvases={}

	fixed_plane, fixed_bar, fixed_SiPM = fixed
	fixed_ch=muAna.MakeFixedCh(fixed)

	if iteration==1:
		readparams=muAna.GetPolyparams_i(afswork, fixed_ch, iteration)
		polyparams=[ readparams[i] for i in tuple(filter( lambda i : i%2==0, range(len(readparams)) )) ]
		correctionfunction=lambda ps, qdc : 1./sum( [ps[i]*qdc**i for i in range(len(ps))] )

	try:	
		ut.readHists(dy_hists, "/afs/cern.ch/user/a/aconsnd/dy_distributions/dy_distributions_run"+str(options.runNumber)+".root")
	except OSError:
		print("Run extrapolationCutAnalysis.py for both x and y axes for run "+str(options.runNumber)+'.')

	if Nev==-1: Nev=eventTree.GetEntries()
	for i, event in enumerate(eventTree):

		if i%100000 == 0: print("============== \nevent: ", i,'/',int(Nev),'\n')
		if i>Nev:break

		if oneUShitperPlane: 
			if not muAna.OneHitPerUS(eventTree.Digi_MuFilterHits): continue

		###### Get DS track if possible. else, continue
		theTrack, stations=DS_track()
		if not hasattr(theTrack, "getFittedState"):
			continue
		if not theTrack.getFitStatus().isFitConverged() and options.track!='DS':
			theTrack.Delete()
			continue			
		# if not hasattr(theTrack, "getFittedState"): continue
		state=theTrack.getFittedState()
		pos=state.getPos()
		mom=state.getMom()
		########

		###### Get T0 estimate from DS horizontal SiPMs used in single track event selection
		DS_T0=muAna.GetDSH_average(stations)*6.25
		if DS_T0 == -999.: return stations # Multiple hits assigned to one bar
		if DS_T0 == -998.: continue # Only 1 SiPM in a horizontal bar is firing	

		######## Determine average time measured by all fired bars in DS track fit. Store slope_x and y for histogram
		# return stations
		for j, aHit in enumerate(eventTree.Digi_MuFilterHits):
			if not aHit.isValid(): continue
			nSides=aHit.GetnSides()
			nSiPMs=aHit.GetnSiPMs()
			if nSiPMs == 1: continue # Only using US hits

			######## cuts on nFired SiPMs and QDC/SiPM
			nFiredSiPMs_left, nFiredSiPMs_right = muAna.GetnFiredSiPMs(aHit)
			if cuts['QDC/SiPM']: # set to False
				if TotalQDC/TotalFiredSiPMs<15: continue
			if cuts['nSiPMs']: # set to True
				if nFiredSiPMs_left<5 or nFiredSiPMs_right<5: continue     # Minimise noise hits contribution
			########

			########### dy-cut
			detID=aHit.GetDetectorID()
			subsystem, plane, bar = muAna.parseDetID(detID)

			if fixed != None:
				fixed_plane, fixed_bar, SiPM = fixed
				if fixed_plane != plane: continue
				if fixed_bar != bar: continue

			tagged=False
			key=subsystem*10+plane
			z=zPos[key]
			lam = (z-pos.z())/mom.z()
			Ex = ROOT.TVector3(pos.x()+lam*mom.x(), pos.y()+lam*mom.y(), pos.z()+lam*mom.z())
			MuFilter.GetPosition(detID, v1, v2)
			dy=Ex.y() - v1.y()
			dy_cut_min, dy_cut_max = dy_hists[str(key)+'yEx-bary'].GetMean()-dy_hists[str(key)+'yEx-bary'].GetStdDev(), dy_hists[str(key)+'yEx-bary'].GetMean()+dy_hists[str(key)+'yEx-bary'].GetStdDev()
			# distance between extrapolated-y and fired bar y-midpoint
			if subsystem==1 or subsystem==2: # Cut on xEx - bar midpoint from extrapolationCutInvestigation.py
				if dy>dy_cut_min and dy<dy_cut_max: tagged=True
			if not tagged: continue
			########

			######## x-pred is the extrapolated x-position in the bar's coordinate system
			xpred, xpred_R=v1.x() - Ex.x(), Ex.x() - v2.x()
			if xpred<0. or xpred>USbarlength: continue
			if cuts['DSacceptance']: 
				if xpred>65:continue
			########

			######## Arrays of all SiPM TDCs and QDCs
			channels_t = aHit.GetAllTimes()
			channels_qdc = aHit.GetAllSignals()

			clock=muAna.GetChannelVal(SiPM, channels_t)
			qdc=muAna.GetChannelVal(SiPM, channels_qdc)

			# rootfile=afswork+'rootfiles/timewalk_dists_'+fixed_ch+'.root'

			if iteration==1:
				## if I have the offset I need to add this separately 
				to_apply=correctionfunction(polyparams, qdc)
				time=clock*6.25
				correctedtime=time+to_apply

			else: correctedtime=clock*6.25
			if correctedtime==-998.: continue
			# print(correctedtime)

			if muAna.IsSmallSiPMchannel(SiPM): continue # Only use large SiPMs
			qdc=muAna.GetChannelVal(SiPM, channels_qdc)
			if qdc==-999.: continue # Only use channels with both a TDC and QDC	

			tuncorrhist='T0-tuncorrsipm_'+str(detID)+'_'+str(SiPM)
			if not tuncorrhist in xQDC:
				title='T_{0}^{DS}-t_{'+str(SiPM)+'} uncorrected;T_{0}^{DS}-t_{'+str(SiPM)+'};Counts'
				xQDC[tuncorrhist]=ROOT.TH1F(tuncorrhist, title, 200, -10, 10)
			tuncorr_rel=DS_T0 - clock*6.25
			xQDC[tuncorrhist].Fill(tuncorr_rel)

			tcorrhist='T0-tcorrsipm_'+str(detID)+'_'+str(SiPM)
			if not tcorrhist in xQDC:
				title='T_{0}^{DS}-t_{'+str(SiPM)+'} corrected;T_{0}^{DS}-t_{'+str(SiPM)+'};Counts'
				xQDC[tcorrhist]=ROOT.TH1F(tcorrhist, title, 200, -10, 10)
			tcorr_rel=DS_T0 - correctedtime
			xQDC[tcorrhist].Fill(tcorr_rel)	

			qdchist=str(detID)+'_'+str(SiPM)+'_'+str(xpred//xwidth) 
			if not qdchist in xQDC:
				barrange=str(int(xpred)//xwidth*xwidth)+' - '+str((int(xpred)//xwidth+1)*xwidth)
				title='QDC distribution for '+str(detID)+' SiPM '+str(SiPM)+' '+barrange
				xQDC[qdchist]=ROOT.TH1F(qdchist, title, 60, 0., 60.)
			xQDC[qdchist].Fill(qdc)                                                                                                                                                                                                       

	# return xQDC
	
	if save:
		if xwidth==USbarlength/3.:timeoutfile='xQDC_'+fixed_ch+'_i'+str(iteration)+'.root'
		else: timeoutfile=afsuser+'qdcslices/xQDC'+str(int(xwidth))+'_'+fixed_ch+'.root'
		ut.writeHists(xQDC, timeoutfile)                                                     
		print('Slice hists written to', timeoutfile)

	if mode=='overlayhists':
		from os.path import exists
		legends={}
		entrieslist=[(i,xQDC[i].GetBinContent(xQDC[i].GetMaximumBin())) for i in xQDC]
		entrieslist.sort(key=lambda x: x[1])
		outfile_name='QDCslices_run'+str(options.runNumber)+'.root'
		if exists(outfile_name): canvoutfile=ROOT.TFile.Open(outfile_name, 'RECREATE')
		else: canvoutfile=ROOT.TFile.Open(outfile_name, 'new')
		colours = {0:ROOT.kBlack,1:ROOT.kRed,2:ROOT.kAzure,3:ROOT.kBlue,4:ROOT.kMagenta,5:ROOT.kCyan,
		6:ROOT.kAzure,7:ROOT.kPink,8:ROOT.kSpring}
		styles={0:ROOT.kFullCircle, 1:ROOT.kFullSquare, 2:ROOT.kFullTriangleDown, 3:ROOT.kFullCrossX}
		# return xQDC
		for idx, i in enumerate(entrieslist):
			histname, entries=i
			detID, SiPM, slicenumber = histname.split('_')
			canvname='_'.join( (detID, SiPM) )
			if canvname not in canvases:
				title='Overlayed QDC histograms for '+detID+' SiPM '+SiPM
				canvases[canvname]=ROOT.TCanvas(canvname, title, 1800, 1500)
				legends[canvname]=ROOT.TLegend()
			canvases[canvname].cd()
			xQDC[histname].SetMarkerSize(3)
			xQDC[histname].SetMarkerColor(colours[idx])
			xQDC[histname].SetLineColor(colours[idx])
			xQDC[histname].SetMarkerStyle(styles[idx])
			xQDC[histname].SetTitle(histname+';QDC value [a.u];Normalised counts')
			xQDC[histname].Scale(1./xQDC[histname].GetEntries())
			xQDC[histname].Draw('same')

		for i in range(int(USbarlength/xwidth)):
			legends[canvname].AddEntry(xQDC['_'.join( (fixed_ch, str(float(i))) )], '/'.join( (str(i), str(int(USbarlength/xwidth))) ))
		legends[canvname].Draw()

		canvoutfile.WriteObject(canvases[canvname], canvname, 'kOverwrite')
		print('canvases written to', canvoutfile.GetName())

	if mode=='QDCmodegraphs':
		graphsoutfile_name='x-QDC_graphs_run'+str(options.runNumber)+'.root'
		graphsoutfile=ROOT.TFile.Open(graphsoutfile_name, "update")
		counters={}
		for histname in xQDC:
			detID, SiPM, slicenumber = histname.split('_')
			graphname='_'.join( (detID, SiPM) )
			if graphname not in graphs:
				counters[graphname]=0
				title='Model QDC variation measured over 4 cm slices in x of '+str(detID)+' SiPM '+str(SiPM)+';x slice lower bound [cm]; Modal QDC'
				graphs[graphname] = ROOT.TGraphErrors()
				graphs[graphname].SetNameTitle(graphname, title)

			modalQDC=xQDC[histname].GetMaximumBin()
			graphs[graphname].SetPoint(counters[graphname], float(slicenumber), modalQDC)
			counters[graphname]+=1

		for graphname in graphs:
			graphsoutfile.WriteObject(graphs[graphname], graphname)

		print('Graphs written to', graphsoutfile_name)

	if save:
		timeoutfile='xQDC_'+fixed_ch+'_i'+str(iteration)+'.root'
		ut.writeHists(xQDC, timeoutfile)
		print('Time hists written to', timeoutfile)

def cscintComparison(fixed=tuple(options.fixed), save=True):
	
	fixed_ch=muAna.MakeFixedCh(fixed)
	# canv=ROOT.TCanvas('cscintComparison', 'cscintComparison', 1200, 1500)
	f=ROOT.TFile.Open(outpath+'rootfiles/timewalk_dists_'+fixed_ch+'.root', 'READ')
	hist=f.Get('linearfit_'+fixed_ch+'_0')
	legend=ROOT.TLegend(0.1,0.7,0.48,0.9)
	l1=f.Get('linearfit_'+fixed_ch+'_1').GetPrimitive('funcdtvxpred_'+fixed_ch+'_iteration1')
	canv.cd()
	l1.SetLineColor(ROOT.kBlack)
	l1.Draw('same')

	legend.AddEntry('funcdtvxpred_'+fixed_ch+'_iteration0', 'Uncorrected linear fit')
	legend.AddEntry('funcdtvxpred_'+fixed_ch+'_iteration1', 'Corrected linear fit')

	return canv

def cscintInstability(fixed=tuple(options.fixed), iteration=1, N=5):

	# Determine cscint for i*1e5 in range(N)
	fixed_ch=muAna.MakeFixedCh(fixed)
	fixed_plane, fixed_bar, fixed_SiPM = fixed

	h={}
	data={}
	if iteration==0: mode='zeroth'
	else: mode='TW'
	canvas=ROOT.TCanvas('cscintInstability', 'cscintInstability', 1500, 1800)
	# canvas.Divide(2, 3)
	for i in range(N):
		# canvas.cd(i+1)
		print('='*10,'\nMaking distributions with', int((i+1)*2e5), 'events.\n'+'='*10)
		
		dtvxpred_dists=makeDistributions(iteration, mode, Nev=(i+1)*2e5, fixed=tuple(options.fixed), oneUShitperPlane=True, save=False)
		histname='dtvxpred_'+fixed_ch+'_iteration'+str(iteration)
		h[histname+'_'+str((i+1)*2e5)]=dtvxpred_dists[histname].Clone()
		hist=h[histname+'_'+str((i+1)*2e5)]
		hist.Draw('colz')

		h['g'+histname+'_'+str((i+1)*2e5)] = ROOT.TGraph()
		g=h['g'+histname+'_'+str((i+1)*2e5)]
		xpred_proj=hist.ProjectionX() # TH1D of xpred 
		for xpred_bin in range(5, 75):
			tmp = hist.ProjectionY('tmp', xpred_bin, xpred_bin) # TH1D of (2x-L) for the given dt bin
			if tmp.GetEntries() < 20: continue
			x = xpred_proj.GetBinCenter(xpred_bin)
			dt = tmp.GetMean() #
			g.SetPoint(xpred_bin-1, x, dt) # SetPoint(i, x, y) - ith point with value (x, y)
		print('Number of entries in graph:', g.GetN())
		if g.GetN()<10:
			print('Too few points on graph for accurate fit')
			from sys import exit 
			exit()
		g.SetLineColor(ROOT.kRed)
		g.SetLineWidth(2)
		g.SetTitle(hist.GetTitle()) # Keep title of histogram from unfitted version

		# g.Draw('same')
		rc=g.Fit('pol1','S', '', 5., 75.)
		result=rc.Get()
		try:
			m=result.Parameter(1)
			intercept=result.Parameter(0)
		except ReferenceError:
			print('Cannot extract fit result')
			from sys import exit 
			exit()
		m_err=result.ParError(1)
		intercept_err=result.ParError(0)
		h['func'+histname]=ROOT.TF1('func'+histname, '[0]*x+[1]', 0, 80)
		h['func'+histname].SetParameter(0, m)
		h['func'+histname].SetParameter(1, intercept)
		h['func'+histname].Draw('same')
		c, cerr = -1/m, (m_err/m**2)
		data[str((i+1)*2e5)]=(c, cerr)
		print(histname, " c_scint +/- err [cm/ns]: ", c, cerr )

	cscintInstability=ROOT.TGraphErrors()
	name, title = 'cscintInstability_iteration'+str(iteration), 'c_{scint} values determined from iteration'+str(iteration)+' for N events;Number of events;c{scint} [cm/ns]'
	cscintInstability.SetNameTitle(name,title)
	for idx, entry in enumerate(data):
		c, cerr = data[entry]
		cscintInstability.SetPoint(idx, int(float(entry)), c)
		cscintInstability.SetPointError(idx, 0, cerr)
	canvas.cd()
	cscintInstability.Draw()
	canvas.Print('cscintInstability'+fixed_ch+'_i'+str(iteration)+'.root')

def dtoutliers(Nev=-1, fixed=tuple(options.fixed), oneUShitperPlane=True, save=True):

	d={}
	d[0]=ROOT.TLine(0., 0., 60., 0.)
	finalcanv=ROOT.TCanvas('finalcanv', 'finalcanv', 1800,2100)
	finalcanv.Divide(1,3)

	fixed_ch=muAna.MakeFixedCh(fixed)
	filename='rootfiles/timewalk_dists_'+fixed_ch+'.root'

	readparams=muAna.GetPolyparams_i(afswork, fixed_ch, iteration)
	polyparams=[ readparams[i] for i in tuple(filter( lambda i : i%2==0, range(len(readparams)) )) ]
	correctionfunction=lambda ps, qdc : 1./sum( [ps[i]*qdc**i for i in range(len(ps))] )	

	f=ROOT.TFile.Open(afswork+filename, 'READ')
	# hist=f.Get('linearfit_'+fixed_ch+'_0')
	legend=ROOT.TLegend(0.1,0.7,0.48,0.9)
	l1=f.Get('linearfit_'+fixed_ch+'_1').GetPrimitive('funcdtvxpred_'+fixed_ch+'_iteration1')
	m, intercept=l1.GetParameter(0), l1.GetParameter(1)
	
	dy_hists={}
	try:	
		ut.readHists(dy_hists, "/afs/cern.ch/user/a/aconsnd/dy_distributions/dy_distributions_run"+str(options.runNumber)+".root")
	except OSError:
		print("Run extrapolationCutAnalysis.py for both x and y axes for run "+str(options.runNumber)+'.')

	if Nev==-1: Nev=eventTree.GetEntries()
	for i, event in enumerate(eventTree):

		if i%100000 == 0: print("============== \nevent: ", i,'/',int(Nev),'\n')
		if i>Nev:break

		if oneUShitperPlane: 
			if not muAna.OneHitPerUS(eventTree.Digi_MuFilterHits): continue

		###### Get DS track if possible. else, continue
		theTrack, stations=DS_track()
		if not hasattr(theTrack, "getFittedState"):
			continue
		if not theTrack.getFitStatus().isFitConverged() and options.track!='DS':
			theTrack.Delete()
			continue			
		# if not hasattr(theTrack, "getFittedState"): continue
		state=theTrack.getFittedState()
		pos=state.getPos()
		mom=state.getMom()
		########

		###### Get T0 estimate from DS horizontal SiPMs used in single track event selection
		DS_T0=muAna.GetDSH_average(stations)*6.25
		if DS_T0 == -999.: return stations # Multiple hits assigned to one bar
		if DS_T0 == -998.: continue # Only 1 SiPM in a horizontal bar is firing	
		if DS_T0 == -6237.5: continue # A bad value permeating the whole data set of run46

		######## Determine average time measured by all fired bars in DS track fit. Store slope_x and y for histogram
		# return stations
		for j, aHit in enumerate(eventTree.Digi_MuFilterHits):
			if not aHit.isValid(): continue
			nSides=aHit.GetnSides()
			nSiPMs=aHit.GetnSiPMs()
			if nSiPMs == 1: continue # Only using US hits

			######## cuts on nFired SiPMs and QDC/SiPM
			nFiredSiPMs_left, nFiredSiPMs_right = muAna.GetnFiredSiPMs(aHit)
			if cuts['QDC/SiPM']: # set to False
				if TotalQDC/TotalFiredSiPMs<15: continue
			if cuts['nSiPMs']: # set to True
				if nFiredSiPMs_left<5 or nFiredSiPMs_right<5: continue     # Minimise noise hits contribution
			########

			########### dy-cut
			detID=aHit.GetDetectorID()
			subsystem, plane, bar = muAna.parseDetID(detID)

			if fixed != None:
				fixed_plane, fixed_bar, SiPM = fixed
				if fixed_plane != plane: continue
				if fixed_bar != bar: continue

			tagged=False
			key=subsystem*10+plane
			z=zPos[key]
			lam = (z-pos.z())/mom.z()
			Ex = ROOT.TVector3(pos.x()+lam*mom.x(), pos.y()+lam*mom.y(), pos.z()+lam*mom.z())
			MuFilter.GetPosition(detID, v1, v2)
			dy=Ex.y() - v1.y()
			dy_cut_min, dy_cut_max = dy_hists[str(key)+'yEx-bary'].GetMean()-dy_hists[str(key)+'yEx-bary'].GetStdDev(), dy_hists[str(key)+'yEx-bary'].GetMean()+dy_hists[str(key)+'yEx-bary'].GetStdDev()
			# distance between extrapolated-y and fired bar y-midpoint
			if subsystem==1 or subsystem==2: # Cut on xEx - bar midpoint from extrapolationCutInvestigation.py
				if dy>dy_cut_min and dy<dy_cut_max: tagged=True
			if not tagged: continue
			########

			######## x-pred is the extrapolated x-position in the bar's coordinate system
			xpred, xpred_R=v1.x() - Ex.x(), Ex.x() - v2.x()
			if xpred<0. or xpred>USbarlength: continue
			########

			######## Arrays of all SiPM TDCs and QDCs
			channels_t = aHit.GetAllTimes()
			channels_qdc = aHit.GetAllSignals()

			clock=muAna.GetChannelVal(SiPM, channels_t)
			qdc=muAna.GetChannelVal(SiPM, channels_qdc)

			if iteration==1:
				to_apply=correctionfunction(polyparams, qdc)
				time=clock*6.25
				correctedtime=time+to_apply
			else: correctedtime=clock*6.25

			if correctedtime==-998.: continue
			# print(correctedtime)

			if muAna.IsSmallSiPMchannel(SiPM): continue # Only use large SiPMs
			qdc=muAna.GetChannelVal(SiPM, channels_qdc)
			if qdc==-999.: continue # Only use channels with both a TDC and QDC	

			residualsname='linearfit-data_'+fixed_ch
			if not residualsname in d:
				title='linearfit - (T_{0}^{DS} - t_{'+str(SiPM)+'}) '+fixed_ch+';linearfit - (T_{0}^{DS} - t_{'+str(SiPM)+'}) [ns];Counts'
				d[residualsname]=ROOT.TH1D(residualsname, title, 200, -20, 20)

			resvqdcname='linearfit-datavqdc_'+fixed_ch
			if not resvqdcname in d:
				title='linearfit - (T_{0}^{DS} - t_{'+str(SiPM)+'}) '+fixed_ch+';QDC_{'+str(SiPM)+'};linearfit - (T_{0}^{DS} - t_{'+str(SiPM)+'}) [ns]'
				d[resvqdcname]=ROOT.TH2F(resvqdcname, title, 60, 0, 60, 200, -20, 20)

			DST0vtSiPMname='DST0vtSiPM_'+fixed_ch
			if not DST0vtSiPMname in d:
				title='(T_{0}^{DS} - t_{'+str(SiPM)+'}) '+fixed_ch+';t_{'+str(SiPM)+'} [ns];(T_{0}^{DS} [ns]'
				d[DST0vtSiPMname]=ROOT.TH2F(DST0vtSiPMname, title, 40, -20, 20, 40, -20, 20)

			t_rel = DS_T0 - correctedtime
			residual=(m*xpred+intercept) - t_rel 
			if abs(residual)>5: 
				# print(DS_T0, correctedtime)
				d[DST0vtSiPMname].Fill(correctedtime, DS_T0)

			d[residualsname].Fill(residual)
			d[resvqdcname].Fill(qdc,residual)

	if save:
		counter=0
		for hist in d:
			if hist==0: continue
			counter+=1
			finalcanv.cd(counter)
			if hist =='linearfit-datavqdc_'+fixed_ch: 
				d[hist].Draw('colz')	
				d[0].Draw('same')
				continue
			d[hist].Draw()
			
		outfile=ROOT.TFile.Open(afsuser+'T0tSiPM_residuals_'+fixed_ch+'.root', 'RECREATE')
		outfile.WriteObject(finalcanv, 'T0-tSiPM_OutliersInvestigation')
		print("Residuals saved to", outfile)
		outfile.Close()
	else: return d

def AttenuationLength(fixed=tuple(options.fixed), save=True):

	d={}
	
	fixed_ch=muAna.MakeFixedCh(fixed)
	SiPM=fixed[2]
	infile=afsuser+'qdcslices/xQDC2_'+fixed_ch+'.root'
	xwidth=2
	try:
		f=ROOT.TFile.Open(infile)
	except OSError:
		print('Nae file')
		return 69 

	ut.readHists(d, infile)
	langauparams={}

	params = {0:'Factor',1:'AttenuationLength',2:'Offset'}
	counter=0
	for histname in d:
		if histname.find(fixed_ch) != 0: continue

		hist=d[histname]
		Xrange=float(histname.split('_')[-1])*2
		if Xrange<5 or Xrange>60: continue
		if hist.GetEntries()<100:continue
		res=muAna.fit_langau(hist, '', 0, 60)
		#### check if fit was properly attempted
		# if res.Chi2()/res.Ndf() > 2: continue
		if res.ParName(1) != 'mostProbable': continue
		####

		langauparams[Xrange] = [[res.Parameter(i), res.ParError(i)] for i in range(4)]
	f.Close()
	xqdcgraph=ROOT.TGraphErrors()
	name, title='qdcvx_'+fixed_ch, 'QDC most probable value evaluated for '+str(xwidth)+' cm sections for '+fixed_ch+';x-position [cm];Log(QDC most probable value) [a.u]'
	xqdcgraph.SetNameTitle(name,title)

	fitfunction=ROOT.TF1('fitfunction', 'pol1', 5,60)
	for i,xbin in enumerate(langauparams):
		mpv,mpverror=langauparams[xbin][1][0], langauparams[xbin][1][1]
		xqdcgraph.SetPoint(i,xbin,ROOT.TMath.Log(mpv))
		xqdcgraph.SetPointError(i,1,mpverror/mpv)

	attlencanv=ROOT.TCanvas('AttenuationLengthCanvas', 'Determination of attenuation length of '+fixed_ch.split('_')[0]+' SiPM '+fixed_ch.split('_')[1]+'x-position [cm]; Most probable QDC value from fit [a.u]', 900, 600)
	attlencanv.cd()
	xqdcgraph.Draw('A *')
	if SiPM < 8:
		fitfunction.SetParLimits(1, 100, 500)
	else:
		fitfunction.SetParLimits(1, -500, -100)

	rc=xqdcgraph.Fit('fitfunction', 'S')
	res=rc.Get()
	failflag=False
	try:
		res.Print()
	except ReferenceError:
		failflag=True
	if not failflag:
		for idx, p in enumerate(res.Parameters()): fitfunction.SetParameter(idx, p)
		attlencanv.cd()
		fitfunction.Draw('same')

		m=res.Parameter(1)
		m_err=res.ParError(1)
		if SiPM < 8:
			attlength=-1./m
			attlength_err=abs(m_err/m)
		else:
			attlength=1./m
			attlength_err=abs(m_err/m)

		chi2pNDF = res.Chi2()/res.Ndf()

		print(attlength)
	if save: 
		attlencanv.Print(afswork+'attenuationlengths/rootfiles/attenuationlength_'+fixed_ch+'.root')
		filename = afswork+'attenuationlengths/csvfiles/attenuationlength_'+fixed_ch+'.csv'
		if not failflag:
			with open(filename, 'w') as h:
				writer=csv.writer(h)
				line = (attlength, attlength_err, chi2pNDF)
				print(line)
				writer.writerow(line)

def compareAttenuationLengths(wanted=((0,1,2,3,4), (3,4), (3,4,11,12)), save=True ):
	path=afswork+'attenuationlengths/'
	
	if wanted != None:
		wanted_planes, wanted_bars, wanted_SiPMs = wanted

	attlengraph=ROOT.TGraphErrors()

	for idx,filename in enumerate(os.listdir(path)):
		print(filename)
		key=filename[filename.find('2'):filename.find('.')]
		# detID, SiPM =key.split('_')[1:]
		detID, SiPM = key
		print(key)

		subsystem, plane, bar = muAna.parseDetID(int(detID))
		
		if plane not in wanted_planes or bar not in wanted_bars or int(SiPM) not in wanted_SiPMs: continue

		with open(path+filename, 'r') as handle:
			reader=csv.reader(handle)
			alldata=[r for r in reader]
			data=alldata[0]

		channelnumber=muAna.GetSiPMNumberInSystem_LandR(int(detID), int(SiPM))
		attlen, attlenerr, chi2pNDF = data 
		print(attlen, key)
		attlengraph.SetPoint(idx, channelnumber, float(attlen))
		attlengraph.SetPointError(idx, channelnumber, float(attlenerr))

	return attlengraph
"""
###
if options.HTCondor:
        modes={'c0':[0, 'zeroth'], 'getc0':[0], 'tw1':[0], 'gettw1':[0], 'c1':[1], 'getc1':[1], 'tw2':[1], 'gettw2':[1]}
        if options.mode not in modes.keys():
            from sys import exit
            exit()
        mode=options.mode
        if mode=='c0': makeDistributions(0, 'zeroth')
        elif mode=='getc0':

            fixed_subsystem=options.subsystem
            fixed_plane=options.plane
            if fixed_subsystem==1: 
                bars=range(7)
                SiPMs=range(12)
            elif fixed_subsystem==2:
                bars=range(10)
                SiPMs=(0,1,3,4,6,7,8,9,11,12,14,15)
            for bar in bars:
                for SiPM in SiPMs:
                    fixed=(fixed_subsystem, fixed_plane, bar, SiPM)
                    determine_cscint(0, fixed)

        elif mode=='tw1': makeDistributions(1, 'ToF')
        elif mode=='gettw1': 
 
            fixed_subsystem=options.subsystem
            fixed_plane=options.plane
            if fixed_subsystem==1: 
                bars=range(7)
                SiPMs=range(12)
            elif fixed_subsystem==2:
                bars=range(10)
                SiPMs=(0,1,3,4,6,7,8,9,11,12,14,15)
            for bar in bars:
                for SiPM in SiPMs:
                    fixed=(fixed_subsystem, fixed_plane, bar, SiPM)
                    determinePolyParams(1, fixed)

        elif mode=='c1': makeDistributions(1, 'TW')
        elif mode=='getc1':
            fixed_subsystem=options.subsystem
            fixed_plane=options.plane
            if fixed_subsystem==1: 
                bars=range(7)
                SiPMs=range(12)
            elif fixed_subsystem==2:
                bars=range(10)
                SiPMs=(0,1,3,4,6,7,8,9,11,12,14,15)
            for bar in bars:
                for SiPM in SiPMs:
                    fixed=(fixed_subsystem, fixed_plane, bar, SiPM)
                    determine_cscint(1, fixed)

#if options.Profiler:
#	makeDistributions(0, 'zeroth', Nev=int(options.nEvents), fixed=tuple(options.fixed), save=False)

#if options.AttLength:
#	xQDC_dependence(Nev=int(options.nEvents), iteration=0, xwidth=2, fixed=tuple(options.fixed), oneUShitperPlane=True, overlayhists=False, QDCmodegraphs=False, save=True)
#	AttenuationLength(fixed=tuple(options.fixed))
