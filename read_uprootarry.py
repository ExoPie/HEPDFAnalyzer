import uproot 
import numpy 
import pandas as pd 
import math 
import time

start = time.clock()

tree_ = uproot.open("Merged_WJetsInclusiveSkim.root")["outTree"]
#tree_ = uproot.open("Mergged_data_met.root")["outTree"]


def getpt_eta_phi(mupx, mupy,mupz):
    mupt = numpy.sqrt(mupx**2 + mupy**2)
    mup = numpy.sqrt(mupx**2 + mupy**2 + mupz**2)
    mueta = numpy.log((mup + mupz)/(mup - mupz))/2
    muphi = numpy.arctan2(mupy, mupx)
    return (mupt, mueta, muphi)

def Phi_mpi_pi(x):
    kPI=numpy.array(3.14159265)
    kPI = kPI.repeat(len(x))
    kTWOPI = 2 * kPI

    while ((x.any() >= kPI).any()): x = x - kTWOPI;
    while ((x.any() < -kPI).any()): x = x + kTWOPI;
    return x;

def DeltaPhi(phi1,phi2):
    phi = Phi_mpi_pi(phi1-phi2)
    return abs(phi)

def getrecoil(nEle,elept,elephi,elepx_,elepy_,met_,metphi_):
    dummy=-9999.0
    WenuRecoilPt=dummy; WenurecoilPhi=-10.0;  We_mass=dummy;
    if (nEle == 1).any():
        
        dphi = DeltaPhi(elephi,metphi_)

        MT = numpy.sqrt( 2 * elept * met * (1.0 - numpy.cos(dphi)) )
        #MT = numpy.sqrt( 2 * elept * met )

        #We_mass = MT(elept[0],met_, DeltaPhi(elephi[0],metphi_)) #transverse mass defined as sqrt{2pT*MET*(1-cos(dphi)}
        WenuRecoilPx = -( met_*numpy.cos(metphi_) + elepx_)
        WenuRecoilPy = -( met_*numpy.sin(metphi_) + elepy_)
        WenuRecoilPt = numpy.sqrt(WenuRecoilPx**2  +  WenuRecoilPy**2)
        WenurecoilPhi = numpy.arctan2(WenuRecoilPx,WenuRecoilPy)
    return WenuRecoilPt, WenurecoilPhi, MT




def getFirstElement(x):
    if len(x)>0: return x[0]

def getSecondElement(x):
    if len(x)>1: return x[1]
    
def getNthElement(x,n):
    if len(x)>n-1: return x[n-1]

def getMinimum(x):
    if len(x)>0: return min(x)
        
def countTrue(x):
    if len(x)>0: return numpy.sum(x)

nevent=100000

total_events = len(tree_ )
print ("total events: ", total_events)

entrystart=0
entrystop=0


steps = (total_events/nevent)+1
df_total_wmunu1b = pd.DataFrame()
df_total_wmunu2b = pd.DataFrame()
for i in range (steps): 
    entrystart=i*nevent
    entrystop=(i+1)*nevent

    ## event 
    run, lumi, event, met, metphi, mettrig, nmuloose, mupx, mupy, mupz, mue, muid, elepx, elepy, elepz, elee,elelooseid, \
        ntau, npho, njet, jetpx, jetpy, jetpz, jete, jetdeepcsv   = tree_.arrays(["st_runId", "st_lumiSection", "st_eventId", 
                                                                                  "st_pfMetCorrPt","st_pfMetCorrPhi","st_mettrigdecision",
                                                                                  "st_nMu","st_muPx", "st_muPy", "st_muPz", "st_muEnergy","st_isTightMuon",
                                                                                  "st_elePx", "st_elePy", "st_elePz", "st_eleEnergy","st_eleIsPassLoose",
                                                                                  "st_nTau_discBased_TightEleTightMuVeto","st_nPho",
                                                                                  "st_THINnJet", "st_THINjetPx", "st_THINjetPy", "st_THINjetPz", "st_THINjetEnergy", "st_THINjetDeepCSV"], 
                                                                                 entrystart = entrystart, entrystop=entrystop, outputtype=tuple)
    
    
    ## muons 
    mupt, mueta, muphi = getpt_eta_phi(mupx, mupy,mupz)
    mu_sel = (mupt>30) & (muid==True) 
    mu_sel_count = mu_sel.sum() ## this gives total number of true, which infact is the total number of selected muons 
    
    
    ## electrons 
    elept, eleeta, elephi = getpt_eta_phi(elepx, elepy, elepz)
    ele_sel = (elept>10.) & (elelooseid==True)
    ele_sel_count = ele_sel.sum()
    
    
    ## jets 
    jetpt, jeteta, jetphi = getpt_eta_phi(jetpx, jetpy, jetpz)
    jet_sel = (jetpt>50) & (numpy.abs(jeteta)<2.5) 
    jet_sel_count = jet_sel.sum()
    
    
    ## b-jets 
    bjet_m = (jetdeepcsv>0.6321)
    nbjet_m = bjet_m.sum()
    
    ## recoil 
    recoil_Wmu1b, recoilphi_Wmu1b, MT_Wmu1b =  getrecoil(mu_sel_count,mupt,muphi,mupx,mupy,met,metphi)
    
    
    ## dphi 
    dphi_jet_met = DeltaPhi(jetphi, metphi)
    
    
    
    Wmunu1b_df = pd.DataFrame({
        ## event
        'run':run,
        'lumi':lumi,
        'event':event,
        ## MET
        "mettrig":mettrig,
        "met":met,
        "metphi":metphi,
        ## muons 
        "NmuLoose":nmuloose,
        "Nmu":mu_sel_count,
        "mu_sel":mu_sel,
        "mupt":mupt, 
        "mueta":mueta,
        "muphi":muphi,
        "muid":muid,
        ## ele
        "Nele":ele_sel_count,
        "elept":elept,
        "eleeta":eleeta,
        "elephi":elephi,
        "elelooseid":elelooseid,
        ## recoil 
        "recoil_Wmu1b":recoil_Wmu1b,
        "recoilphi_Wmu1b":recoilphi_Wmu1b,
        "MT_Wmu1b":MT_Wmu1b,
        # jets 
        #    "nJet":njet,
        "nJet":jet_sel_count,
        "jetpt":jetpt,
        "jeteta":jeteta,
        "jetphi":jetphi,
        # b-jets 
         'nbjet_m':nbjet_m,
        #taus 
        "nTau":ntau,
        # photon 
        "nPho":npho,
        # dphi 
        "dphi_jet_met":dphi_jet_met,
        
    
        
        
    })
    
    print (len(mupt), len(mueta), len(muphi), len(mettrig), len(muid))
    
    
    
    ## this syntax will choose the required data from the tuple of each object in each event. 
    ## two function at the moment are getFirstElement and getMinimum
    ## in a similar way any function can be written and it should ideally work. 
    ## this method is usually slow and should be avoided if possible, 
    ## but this is still much faster than usual pythonic lopp so must be preferred over loop. 
    
    for ivar in ['recoil_Wmu1b','MT_Wmu1b','mueta','muid','muphi','mupt','recoilphi_Wmu1b', 'elept', 'eleeta', 'elephi', 'elelooseid', 'jetpt', 'jeteta', 'jetphi' ]:
        Wmunu1b_df[ivar] = Wmunu1b_df.apply(lambda x: getFirstElement(x[ivar]), axis=1)
        
    
    for ivar in ['dphi_jet_met']: 
        Wmunu1b_df[ivar] = Wmunu1b_df.apply(lambda x: getMinimum(x[ivar]), axis=1)
    
    
    ## W(mu)+1b CR     
    Wmunu1b_sel   =  Wmunu1b_df[(Wmunu1b_df.mettrig) & 
                                (Wmunu1b_df.Nele==0) & 
                                (Wmunu1b_df.nPho==0) & 
                                (Wmunu1b_df.nTau==0) & 
                                (Wmunu1b_df.Nmu==1) & 
                                (Wmunu1b_df.NmuLoose==1) &
                                (Wmunu1b_df.recoil_Wmu1b>200.)  &
                                (Wmunu1b_df.dphi_jet_met > 0.5) &
                                (Wmunu1b_df.MT_Wmu1b>0) & (Wmunu1b_df.MT_Wmu1b<160) & 
                                (Wmunu1b_df.nbjet_m==1) & 
                                (Wmunu1b_df.nJet==1) &
                                (met>100.)]
    
    
    Wmunu2b_sel   =  Wmunu1b_df[(Wmunu1b_df.mettrig) &
                                (Wmunu1b_df.Nele==0) &
                                (Wmunu1b_df.nPho==0) &
                                (Wmunu1b_df.nTau==0) &
                                (Wmunu1b_df.Nmu==1) &
                                (Wmunu1b_df.NmuLoose==1) &
                                (Wmunu1b_df.recoil_Wmu1b>200.)  &
                                (Wmunu1b_df.dphi_jet_met > 0.5) &
                                (Wmunu1b_df.MT_Wmu1b>0) & (Wmunu1b_df.MT_Wmu1b<160) &
                                (Wmunu1b_df.nbjet_m==2) &
                                (Wmunu1b_df.nJet==2) & 
                                (met>100.)]
    
    import root_pandas 
    
    df_total_wmunu1b = pd.concat([Wmunu1b_sel, df_total_wmunu1b])
    df_total_wmunu2b = pd.concat([Wmunu2b_sel, df_total_wmunu2b])
    
    print (Wmunu1b_sel)
end = time.clock()

print (df_total_wmunu1b)
df_total_wmunu1b.to_root("wmunu.root",key="crwmunu1b")
df_total_wmunu2b.to_root("wmunu.root",key="crwmunu2b",mode='a')




print("%.4gs" % (end-start))
