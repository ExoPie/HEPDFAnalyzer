import sys 
import uproot 
import numpy 
import pandas as pd 
import math 
import time

start = time.clock()

tree_ = uproot.open("Merged_WJetsInclusiveSkim.root")["outTree"]
#tree_ = uproot.open("Mergged_data_met.root")["outTree"]

from utils import *

ieve=0
def getCleanedNumber(eta1, eta2, phi1, phi2):
    print "inside getCleanedNumber", ieve
    ieve = ieve+1
    dr_=[]
    npho = len(eta1)
    njet = len(eta2)
    for ipho in range(npho):
        for ijet in range(njet):
            dr = Delta_R(eta1[ipho], eta2[ijet], phi1[ipho], phi2[ijet])
            dr_.append(dr)
    return len(dr_)

def deltaR(phoeta, jeteta, phophi, jetphi, cut_=0.4):
    
    phoeta_unzip, jeteta_unzip = phoeta.cross(jeteta, nested=True).unzip()
    phophi_unzip, jetphi_unzip = phophi.cross(jetphi, nested=True).unzip()
    
    deta_unzip = phoeta_unzip - jeteta_unzip
    dphi_unzip = phophi_unzip - jetphi_unzip
    
    dr_unzip = numpy.sqrt(deta_unzip**2 + dphi_unzip**2)
    dr_pho_jet_status = (dr_unzip<cut_).any() ## use axis in new version of awkward 
    return dr_pho_jet_status
    
    

df_total_wmunu1b = pd.DataFrame()
df_total_wmunu2b = pd.DataFrame()
df_total_wenu1b = pd.DataFrame()
df_total_wenu2b = pd.DataFrame()

nevent=1000000
total_events = len(tree_ )
entrystart=0
entrystop=0
steps = (total_events/nevent)
if total_events%nevent > 0: steps=steps+1
print ("total events: steps:", total_events, steps)

for i in range (steps): 
    entrystart=i*nevent
    entrystop=(i+1)*nevent

    ## event 
    run, lumi, event,\
        eletrig,\
        met, metphi, mettrig, \
        nmuloose, mupx, mupy, mupz, mue, muid, \
        elepx, elepy, elepz, elee,elelooseid, eletightid, elecharge, \
        ntau, npho, phopx, phopy, phopz, phoe, \
        njet, jetpx, jetpy, jetpz, jete, jetdeepcsv   = tree_.arrays(["st_runId", "st_lumiSection", "st_eventId",
                                                                      "st_eletrigdecision",
                                                                      "st_pfMetCorrPt","st_pfMetCorrPhi","st_mettrigdecision",
                                                                      "st_nMu","st_muPx", "st_muPy", "st_muPz", "st_muEnergy","st_isTightMuon",
                                                                      "st_elePx", "st_elePy", "st_elePz", "st_eleEnergy","st_eleIsPassLoose",'st_eleIsPassTight', 'st_eleCharge',
                                                                      "st_nTau_discBased_TightEleTightMuVeto","st_nPho", "st_phoPx","st_phoPy","st_phoPz","st_phoEnergy",
                                                                      "st_THINnJet", "st_THINjetPx", "st_THINjetPy", "st_THINjetPz", "st_THINjetEnergy", "st_THINjetDeepCSV"], 
                                                                     entrystart = entrystart, entrystop=entrystop, outputtype=tuple)
    
    
    #for ilumi in lumi: print ilumi
    ## muons 
    mupt, mueta, muphi = getpt_eta_phi(mupx, mupy,mupz)
    mu_sel = (mupt>30) & (muid==True) 
    mu_sel_tight_count = mu_sel.sum() ## this gives total number of true, which infact is the total number of selected muons 
    
    
    ## electrons 
    elept, eleeta, elephi = getpt_eta_phi(elepx, elepy, elepz)
    ele_sel = (elept>10.) & (elelooseid==True)
    ele_sel_count = ele_sel.sum()
    ele_sel_tight_count = ((elept>30.) & (eletightid==True)).sum()
    
    #print ("eletightid: ", sum(ele_sel_tight_count.tolist()))
    
    #print ("ele_sel_tight_count: ",ele_sel_tight_count.tolist())
    
    ## jets 
    jetpt, jeteta, jetphi = getpt_eta_phi(jetpx, jetpy, jetpz)
    jet_sel = (jetpt>50) & (numpy.abs(jeteta)<2.5) 
    jet_sel_count = jet_sel.sum()
    
    jet_sel_pt30 = (jetpt>30) & (numpy.abs(jeteta)<2.5)
    jet_sel_pt30_count = jet_sel_pt30.sum()

    ## b-jets 
    bjet_m = (jetdeepcsv>0.6321) & (numpy.abs(jeteta)<2.4)
    nbjet_m = bjet_m.sum()
    
    ## photons 
    phopt, phoeta, phophi = getpt_eta_phi(phopx, phopy, phopz)
    jeteta_sel_pt30 = jeteta[jet_sel_pt30]
    jetphi_sel_pt30 = jetphi[jet_sel_pt30]
    
    
    
    ## recoil 
    recoil_Wmu1b, recoilphi_Wmu1b, MT_Wmu1b =  getrecoil(mu_sel_tight_count,mupt,muphi,mupx,mupy,met,metphi)
    recoil_We1b, recoilphi_We1b, MT_We1b    =  getrecoil(ele_sel_tight_count,elept,elephi,elepx,elepy,met,metphi)

    
    ## dphi 
    dphi_jet_met = DeltaPhi(jetphi, metphi)
    
    
    
    common_df = pd.DataFrame({
        ## event
        'run':run,
        'lumi':lumi,
        'event':event,
        ## MET
        "mettrig":mettrig,
        "eletrig":eletrig,
        "met":met,
        "metphi":metphi,
        ## muons 
        "NmuLoose":nmuloose,
        "Nmu":mu_sel_tight_count,
        "mu_sel":mu_sel,
        "mupt":mupt, 
        "mueta":mueta,
        "muphi":muphi,
        "muid":muid,
        ## ele
        "NeleLoose":ele_sel_count,
        "Nele":ele_sel_tight_count,
        "elept":elept,
        "eleeta":eleeta,
        "elephi":elephi,
        "elelooseid":elelooseid,
        ## recoil 
        "recoil_Wmu1b":recoil_Wmu1b,
        "recoilphi_Wmu1b":recoilphi_Wmu1b,
        "MT_Wmu1b":MT_Wmu1b,
        
        "recoil_We1b":recoil_We1b,
        "recoilphi_We1b":recoilphi_We1b,
        "MT_We1b":MT_We1b,
        
        # jets 
        "nJet30":jet_sel_pt30_count,
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
        "phopt":phopt,
        "phoeta":phoeta,
        "phophi":phophi,
        # cleaningn jets 
        "jeteta_sel_pt30":jeteta_sel_pt30,
        "jetphi_sel_pt30":jetphi_sel_pt30,
        # dphi 
        "dphi_jet_met":dphi_jet_met,
        
    
        
        
    })
    
    
    #common_df["ncleanPho"] = 1
    '''
    phoeta_unzip, jeteta_unzip = phoeta.cross(jeteta, nested=True).unzip()
    phophi_unzip, jetphi_unzip = phophi.cross(jetphi, nested=True).unzip()
    
    deta_unzip = phoeta_unzip - jeteta_unzip
    dphi_unzip = phophi_unzip - jetphi_unzip
    
    dr_unzip = numpy.sqrt(deta_unzip**2 + dphi_unzip**2)
    dr_pho_jet_status = (dr_unzip<0.4).any() ## use axis in new version of awkward 
    '''
    
    mask_matched_pho_ele = (deltaR(phoeta, eleeta, phophi, elephi))    
    mask_matched_pho_mu  = deltaR(phoeta, mueta, phophi, muphi)

    mask_matched_pho = (mask_matched_pho_ele | mask_matched_pho_mu)
    npho_cleaned     = phoeta[~mask_matched_pho].sum()
    common_df["nphocleaned"] = npho_cleaned

    
    ## this syntax will choose the required data from the tuple of each object in each event. 
    ## two function at the moment are getFirstElement and getMinimum
    ## in a similar way any function can be written and it should ideally work. 
    ## this method is usually slow and should be avoided if possible, 
    ## but this is still much faster than usual pythonic lopp so must be preferred over loop. 
    
        
    for ivar in ['recoil_Wmu1b','MT_Wmu1b','mueta','muid','muphi','mupt','recoilphi_Wmu1b', 'elept', 'eleeta', 'elephi', 'elelooseid', 'jetpt', 'jeteta', 'jetphi', 'recoil_We1b', 'recoilphi_We1b', 'MT_We1b'  ]:
        common_df[ivar] = common_df.apply(lambda x: getFirstElement(x[ivar]), axis=1)
        
        
    #dphi_jet_met1 = dphi_jet_met.choose(1)
    #common_df['dphi_jet_met'] = dphi_jet_met1
    #3.031477, 3.239381, 3.369493, 4.027412 
    
    for ivar in ['dphi_jet_met']: 
        common_df[ivar] = common_df.apply(lambda x: getMinimum(x[ivar]), axis=1)
    
    #print ("common_df.nJet30: ", common_df.nJet30==2)
    
    
    ## W(mu)+1b CR     
    Wmunu1b_sel   =  common_df[(common_df.NeleLoose==0) & 
                               (common_df.nphocleaned == 0) &
                               (common_df.nTau==0) & 
                               (common_df.Nmu==1) & 
                               (common_df.NmuLoose==1) &
                               (common_df.recoil_Wmu1b>200.)  &
                               (common_df.dphi_jet_met > 0.5) &
                               (common_df.MT_Wmu1b>0) & (common_df.MT_Wmu1b<160) & 
                               (common_df.nbjet_m==1) & 
                               (common_df.nJet==1) &
                               (common_df.nJet30==1)
    ]
    
    


    Wmunu2b_sel   =  common_df[(common_df.NeleLoose==0) &
                               (common_df.nphocleaned == 0) &
                               (common_df.nTau==0) &
                               (common_df.Nmu==1) &
                               (common_df.NmuLoose==1) &
                               (common_df.recoil_Wmu1b>200.)  &
                               (common_df.dphi_jet_met > 0.5) &
                               (common_df.MT_Wmu1b>0) & (common_df.MT_Wmu1b<160) &
                               (common_df.nbjet_m==2) &
                               (common_df.nJet>=1) & 
                               (common_df.nJet30==2)
    ]
    
    
    Wenu1b_sel   =  common_df[(common_df.Nele==1) & 
                              (common_df.NeleLoose==1) &
                              (common_df.nphocleaned == 0) &
                              (common_df.nTau==0) & 
                              (common_df.NmuLoose==0) &
                              (common_df.recoil_We1b>200.)  &
                              (common_df.dphi_jet_met > 0.5) &
                              (common_df.MT_We1b>0) & (common_df.MT_We1b<160) & 
                              (common_df.nbjet_m==1) & 
                              (common_df.nJet==1) &
                              (common_df.nJet30==1)
                          ]


    Wenu2b_sel   =  common_df[(common_df.Nele==1) & 
                              (common_df.NeleLoose==1) &
                              (common_df.nphocleaned == 0) &
                              (common_df.nTau==0) & 
                              (common_df.NmuLoose==0) &
                              (common_df.recoil_We1b>200.)  &
                              (common_df.dphi_jet_met > 0.5) &
                              (common_df.MT_We1b>0) & (common_df.MT_We1b<160) & 
                              (common_df.nbjet_m==2) & 
                              (common_df.nJet>=1) &
                              (common_df.nJet30==2)
                          ]
    
    
    
    df_total_wmunu1b = pd.concat([Wmunu1b_sel, df_total_wmunu1b])
    df_total_wmunu2b = pd.concat([Wmunu2b_sel, df_total_wmunu2b])
    
    df_total_wenu1b  = pd.concat([Wenu1b_sel, df_total_wenu1b])
    df_total_wenu2b  = pd.concat([Wenu2b_sel, df_total_wenu2b])

end = time.clock()

df_total_wmunu1b['weight'] = df_total_wmunu1b.apply( lambda x: weight_( x.recoil_Wmu1b, x.Nele, x.Nmu, x.mupt, x.mueta ),
                                                     axis=1)

import root_pandas 

## We
df_total_wenu1b.to_root("wmunu.root",key="crwenu1b")
df_total_wenu2b.to_root("wmunu.root",key="crwenu2b",mode='a')

## Wmu 
df_total_wmunu1b.to_root("wmunu.root",key="crwmunu1b",mode='a')
df_total_wmunu2b.to_root("wmunu.root",key="crwmunu2b",mode='a')


## Tope 

## Topmu 

## Zee

## Zmumu

## SR1 and SR2 
print("%.4gs" % (end-start))
