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

nevent=1000000

total_events = len(tree_ )
print ("total events: ", total_events)

entrystart=0
entrystop=0

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

def deltaR(phoeta, jeteta, phophi, jetphi):
    
    phoeta_unzip, jeteta_unzip = phoeta.cross(jeteta, nested=True).unzip()
    phophi_unzip, jetphi_unzip = phophi.cross(jetphi, nested=True).unzip()
    
    deta_unzip = phoeta_unzip - jeteta_unzip
    dphi_unzip = phophi_unzip - jetphi_unzip
    
    dr_unzip = numpy.sqrt(deta_unzip**2 + dphi_unzip**2)
    dr_pho_jet_status = (dr_unzip<0.4).any() ## use axis in new version of awkward 
    return dr_pho_jet_status
    
    

steps = (total_events/nevent)+1
df_total_wmunu1b = pd.DataFrame()
df_total_wmunu2b = pd.DataFrame()
for i in range (steps): 
    entrystart=i*nevent
    entrystop=(i+1)*nevent

    ## event 
    run, lumi, event, met, metphi, mettrig, nmuloose, mupx, mupy, mupz, mue, muid, elepx, elepy, elepz, elee,elelooseid, \
        ntau, npho, phopx, phopy, phopz, phoe, njet, jetpx, jetpy, jetpz, jete, jetdeepcsv   = tree_.arrays(["st_runId", "st_lumiSection", "st_eventId", 
                                                                                    "st_pfMetCorrPt","st_pfMetCorrPhi","st_mettrigdecision",
                                                                                    "st_nMu","st_muPx", "st_muPy", "st_muPz", "st_muEnergy","st_isTightMuon",
                                                                                    "st_elePx", "st_elePy", "st_elePz", "st_eleEnergy","st_eleIsPassLoose",
                                                                                    "st_nTau_discBased_TightEleTightMuVeto","st_nPho", "st_phoPx","st_phoPy","st_phoPz","st_phoEnergy",
                                                                                    "st_THINnJet", "st_THINjetPx", "st_THINjetPy", "st_THINjetPz", "st_THINjetEnergy", "st_THINjetDeepCSV"], 
                                                                                    entrystart = entrystart, entrystop=entrystop, outputtype=tuple)
    
    
    #for ilumi in lumi: print ilumi
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
    recoil_Wmu1b, recoilphi_Wmu1b, MT_Wmu1b =  getrecoil(mu_sel_count,mupt,muphi,mupx,mupy,met,metphi)
    
    
    ## dphi 
    dphi_jet_met = DeltaPhi(jetphi, metphi)
    
    
    
    common_df = pd.DataFrame({
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
    
    print (len(mupt), len(mueta), len(muphi), len(mettrig), len(muid))
    
    common_df["ncleanPho"] = 1
    '''
    phoeta_unzip, jeteta_unzip = phoeta.cross(jeteta, nested=True).unzip()
    phophi_unzip, jetphi_unzip = phophi.cross(jetphi, nested=True).unzip()
    
    deta_unzip = phoeta_unzip - jeteta_unzip
    dphi_unzip = phophi_unzip - jetphi_unzip
    
    dr_unzip = numpy.sqrt(deta_unzip**2 + dphi_unzip**2)
    dr_pho_jet_status = (dr_unzip<0.4).any() ## use axis in new version of awkward 
    '''
    
    
    common_df["dr_pho_mu_status"] = (deltaR(phoeta, mueta, phophi, muphi)).sum()
    common_df["dr_pho_ele_status"] = (deltaR(phoeta, eleeta, phophi, elephi)).sum()
    
    
    #common_df["ncleanPho"] = common_df.apply(lambda x: getCleanedNumber(x["phoeta"], x["jeteta_sel_pt30"], x["phophi"], x["jetphi_sel_pt30"] ) )
    
    ## this syntax will choose the required data from the tuple of each object in each event. 
    ## two function at the moment are getFirstElement and getMinimum
    ## in a similar way any function can be written and it should ideally work. 
    ## this method is usually slow and should be avoided if possible, 
    ## but this is still much faster than usual pythonic lopp so must be preferred over loop. 
    
    for ivar in ['recoil_Wmu1b','MT_Wmu1b','mueta','muid','muphi','mupt','recoilphi_Wmu1b', 'elept', 'eleeta', 'elephi', 'elelooseid', 'jetpt', 'jeteta', 'jetphi' ]:
        common_df[ivar] = common_df.apply(lambda x: getFirstElement(x[ivar]), axis=1)
        
        
    #dphi_jet_met1 = dphi_jet_met.choose(1)
    #common_df['dphi_jet_met'] = dphi_jet_met1
    #3.031477, 3.239381, 3.369493, 4.027412 
    for ivar in ['dphi_jet_met']: 
        common_df[ivar] = common_df.apply(lambda x: getMinimum(x[ivar]), axis=1)
    
    
    
    ## W(mu)+1b CR     
    Wmunu1b_sel   =  common_df[(common_df.mettrig) & 
                                (common_df.Nele==0) & 
                                #(common_df.nPho==0) & 
                                (common_df.dr_pho_mu_status==0) &
                                (common_df.dr_pho_ele_status==0) &
                                (common_df.nTau==0) & 
                                (common_df.Nmu==1) & 
                                (common_df.NmuLoose==1) &
                                (common_df.recoil_Wmu1b>200.)  &
                                (common_df.dphi_jet_met > 0.5) &
                                (common_df.MT_Wmu1b>0) & (common_df.MT_Wmu1b<160) & 
                                (common_df.nbjet_m==1) & 
                                (common_df.nJet==1) &
                                (common_df.nJet30==1)
                                #(common_df.met>100.)
    ]
    
    
    Wmunu2b_sel   =  common_df[(common_df.mettrig) &
                                (common_df.Nele==0) &
                                (common_df.nPho==0) &
                                (common_df.nTau==0) &
                                (common_df.Nmu==1) &
                                (common_df.NmuLoose==1) &
                                (common_df.recoil_Wmu1b>200.)  &
                                (common_df.dphi_jet_met > 0.5) &
                                (common_df.MT_Wmu1b>0) & (common_df.MT_Wmu1b<160) &
                                (common_df.nbjet_m==2) &
                                (common_df.nJet==2) & 
                                (common_df.nJet30==1)
                                #(common_df.met>100.)
    ]
    
    
    
    Wenu1b_sel   =  common_df[(common_df.mettrig) & 
                              (common_df.Nele==1) & 
                              (common_df.dr_pho_mu_status==0) &
                              (common_df.dr_pho_ele_status==0) &
                              (common_df.nTau==0) & 
                              (common_df.NmuLoose==0) &
                              (common_df.recoil_Wmu1b>200.)  &
                              (common_df.dphi_jet_met > 0.5) &
                              (common_df.MT_Wmu1b>0) & (common_df.MT_Wmu1b<160) & 
                              (common_df.nbjet_m==1) & 
                              (common_df.nJet==1) &
                              (common_df.nJet30==1)
                              #(common_df.met>100.)
                           ]
    
    import root_pandas 
    
    df_total_wmunu1b = pd.concat([Wmunu1b_sel, df_total_wmunu1b])
    df_total_wmunu2b = pd.concat([Wmunu2b_sel, df_total_wmunu2b])
    
    print (Wmunu1b_sel)
end = time.clock()

df_total_wmunu1b['weight'] = df_total_wmunu1b.apply( lambda x: weight_( x.recoil_Wmu1b, x.Nele, x.Nmu, x.mupt, x.mueta ),#, x.recoil_Wmu1b, x.recoil_Wmu1b, x.Nele, x.elept, x.eleeta, x.elelooseid, x.Nmu, x.mupt, x.mueta, x.muid), 
                                                     axis=1)


#df_total_wmunu1b['weight'] = df_total_wmunu1b.apply( lambda x: weight_(common_weight, ep_pfMetCorrPt, ep_ZmumuRecoil, ep_WmunuRecoil, nEle, ep_elePt, ep_eleEta, ep_eleIsPTight, nMu, ep_muPt, ep_muEta, ep_isTightMuon), 
#                                                    axis=1)

print (df_total_wmunu1b)
df_total_wmunu1b.to_root("wmunu.root",key="crwmunu1b")
df_total_wmunu2b.to_root("wmunu.root",key="crwmunu2b",mode='a')


print("%.4gs" % (end-start))
