## Possible fixes needed 
## Delta Phi is not proper, it will be propagated to all vars require it, like Recoil , mt etc. 
## Add Photon veto using the new skimming files which has DR cleaning cut already applied 

import sys 
import uproot4 
import awkward1 as ak 
import numpy 
import math 
import time
start = time.clock()


inputfile=sys.argv[1]
print ("inputfile:", inputfile)
print ("outputfile: ", "output/"+inputfile.split("/")[-1])
outputfile = "output/"+inputfile.split("/")[-1]


from utils import getpt, geteta, getphi, getrecoil, DeltaPhi, getMT, getRecoilPhi, getrecoil1, getN


print ("opening rootfile")
mycache = uproot4.LRUArrayCache("10 MB")
tree_ = uproot4.open(inputfile)["outTree"].arrays(array_cache=mycache)

#tree_ = uproot4.open("Merged_WJetsInclusiveSkim.root")["outTree"].arrays(array_cache=mycache)
#tree_ = uproot4.open("Merged_WjetsInclusive_photveto.root")["outTree"].arrays(array_cache=mycache)
#tree_ = uproot4.open("/eos/cms/store/group/phys_exotica/bbMET/2016_SkimmedFiles/skim_setup_2016_v16_07-00/crab_DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_200918_215129_0000_0.root")["outTree"].arrays(array_cache=mycache)

#tree_ = uproot4.open("/eos/cms/store/group/phys_exotica/bbMET/2016_SkimmedFiles/skim_setup_2016_v16_07-00/crab_ttHTobb_M125_13TeV_powheg_pythia8_200918_215950_0000_0.root")["outTree"].arrays(array_cache=mycache)



#print ((tree_))

cms_events = ak.zip({   "run":tree_["st_runId"],"lumi":tree_["st_lumiSection"], "event": tree_["st_eventId"],
                        "jetpx":tree_["st_THINjetPx"], "jetpy":tree_["st_THINjetPy"], "jetpz":tree_["st_THINjetPz"], "jete":tree_["st_THINjetEnergy"],
                        "jetpt": getpt(tree_["st_THINjetPx"], tree_["st_THINjetPy"]), "jeteta":geteta(tree_["st_THINjetPx"], tree_["st_THINjetPy"], tree_["st_THINjetPz"]), 
                        "jetphi":getphi(tree_["st_THINjetPx"], tree_["st_THINjetPy"]), "jetcsv": tree_["st_THINjetDeepCSV"],
                        "metpt":tree_["st_pfMetCorrPt"], "metphi": tree_["st_pfMetCorrPhi"], "mettrig": tree_["st_mettrigdecision"],
                        "elepx":tree_["st_elePx"], "elepy":tree_["st_elePy"], "elepz":tree_["st_elePz"], "elee":tree_["st_eleEnergy"],
                        "eleidL":tree_["st_eleIsPassLoose"], "eleidT":tree_["st_eleIsPassTight"], "eleq":tree_["st_eleCharge"],
                        "elept": getpt(tree_["st_elePx"], tree_["st_elePy"]), "eleeta":geteta(tree_["st_elePx"], tree_["st_elePy"], tree_["st_elePz"]),
                        "elephi":getphi(tree_["st_elePx"], tree_["st_elePy"]), 
                        "mupx":tree_["st_muPx"], "mupy":tree_["st_muPy"], "mupz":tree_["st_muPz"], "mue":tree_["st_muEnergy"],
                        "muidT":tree_["st_isTightMuon"], "muq":tree_["st_muCharge"],
                        "mupt": getpt(tree_["st_muPx"], tree_["st_muPy"]), "mueta":geteta(tree_["st_muPx"], tree_["st_muPy"], tree_["st_muPz"]),
                        "muphi":getphi(tree_["st_muPx"], tree_["st_muPy"]),
                        "ntau":tree_["st_nTau_discBased_TightEleTightMuVeto"],
                        "npho":tree_["st_nPho"] ,
                        "phopx":tree_["st_phoPx"],"phopy":tree_["st_phoPy"],"phopz":tree_["st_phoPz"], "phoe":tree_["st_phoEnergy"], 
                        "phopt": getpt(tree_["st_phoPx"], tree_["st_phoPy"]), "phoeta":geteta(tree_["st_phoPx"], tree_["st_phoPy"], tree_["st_phoPz"])
                    },
                    depth_limit=1)


print ("event loading done")
print ("# of events: ",len(cms_events))

## add more columns/properties to the event 


'''
 another syntax to pick any element 
>>> a = ak.mask(cms_events.recoil_Wmunu, ak.num(cms_events.recoil_Wmunu, axis=1)>0, highlevel=False)[:,0]
>>> ak.to_list(a[:10])
[227.8722381591797, None, None, None, 196.6229705810547, 181.11541748046875, 250.32859802246094, None, 250.22610473632812, None]
'''
cms_events["mu_sel_tight"]  = (cms_events.mupt>30) & (cms_events.muidT==True) & (numpy.abs(cms_events.mueta)<2.4)
cms_events["mu_sel_tight0"] = ak.Array(getN(cms_events.mu_sel_tight,0))
cms_events["nMuTight"]      = ak.sum(cms_events.mu_sel_tight, axis=-1) 
cms_events["nMuLoose"]      = ak.sum( (cms_events.mupt>10), axis=-1) 
cms_events["mu_q0"]         = ak.Array(getN(cms_events.muq,0) )
cms_events["mu_q1"]         = ak.Array(getN(cms_events.muq,1)) 

cms_events["ele_sel_tight"] = (cms_events.eleidT==True) & (cms_events.elept>30) & (numpy.abs(cms_events.eleeta)<2.5) 
cms_events["ele_sel_tight0"] = ak.Array(getN(cms_events.ele_sel_tight,0))
cms_events["nEleTight"]     = ak.sum(cms_events.ele_sel_tight, axis=-1)
cms_events["nEleLoose"]     = ak.sum ( (cms_events.elept>10), axis = -1 )
cms_events["ele_q0"]        = ak.Array(getN(cms_events.eleq, 0))
cms_events["ele_q1"]        = ak.Array(getN(cms_events.eleq, 1))


cms_events["recoil_Wmunu"]  = getrecoil(cms_events.nMuTight,  cms_events.mupt, cms_events.muphi, cms_events.mupx, cms_events.mupy, cms_events.metpt, cms_events.metphi)
cms_events["recoil_Wmunu0"] = ak.firsts(cms_events.recoil_Wmunu)

cms_events["recoil_Wenu"]   = getrecoil(cms_events.nEleTight,  cms_events.elept, cms_events.elephi, cms_events.elepx, cms_events.elepy, cms_events.metpt, cms_events.metphi)
cms_events["recoil_Wenu0"]  = ak.firsts(cms_events.recoil_Wenu)

    

elepx0 = ak.Array(getN(cms_events.elepx,0))
elepx1 = ak.Array(getN(cms_events.elepx,1))

elepy0 = ak.Array(getN(cms_events.elepy,0))
elepy1 = ak.Array(getN(cms_events.elepy,1))

elepz0 = ak.Array(getN(cms_events.elepz,0))
elepz1 = ak.Array(getN(cms_events.elepz,1))

elee0 = ak.Array(getN(cms_events.elee,0))
elee1 = ak.Array(getN(cms_events.elee,1))

cms_events["Zee_mass"] = numpy.sqrt(  ( elee0+elee1)**2   - (elepx0+elepx1)**2 - (elepy0+elepy1)**2 - (elepz0+elepz1)**2 )
cms_events["Zee_pt"]      = numpy.sqrt ( (elepx0+elepx1)**2 + (elepy0+elepy1)**2 )
cms_events["Zee_recoil"]  = getrecoil1( (elepx0+elepx1), (elepy0+elepy1), cms_events.metpt, cms_events.metphi)

mupx0 = ak.Array(getN(cms_events.mupx,0))
mupx1 = ak.Array(getN(cms_events.mupx,1))

mupy0 = ak.Array(getN(cms_events.mupy,0))
mupy1 = ak.Array(getN(cms_events.mupy,1))

mupz0 = ak.Array(getN(cms_events.mupz,0))
mupz1 = ak.Array(getN(cms_events.mupz,1))

mue0 = ak.Array(getN(cms_events.mue,0))
mue1 = ak.Array(getN(cms_events.mue,1))

cms_events["Zmumu_mass"]    = numpy.sqrt(  ( mue0+mue1)**2   - (mupx0+mupx1)**2 - (mupy0+mupy1)**2 - (mupz0+mupz1)**2 )
cms_events["Zmumu_pt"]      = numpy.sqrt ( (mupx0+mupx1)**2 + (mupy0+mupy1)**2 )
cms_events["Zmumu_recoil"]  = getrecoil1( (mupx0+mupx1), (mupy0+mupy1), cms_events.metpt, cms_events.metphi)


#cms_events["recoil_Zmumu"] = getrecoil
cms_events["recoil_WmunuPhi"] = getRecoilPhi(cms_events.nMuTight,  cms_events.mupt, cms_events.muphi, cms_events.mupx, cms_events.mupy, cms_events.metpt, cms_events.metphi)
cms_events["recoil_WmunuPhi0"] = ak.firsts(cms_events.recoil_WmunuPhi)

cms_events["recoil_WenuPhi"] = getRecoilPhi(cms_events.nEleTight,  cms_events.elept, cms_events.elephi, cms_events.elepx, cms_events.elepy, cms_events.metpt, cms_events.metphi)
cms_events["recoil_WenuPhi0"] = ak.firsts(cms_events.recoil_WenuPhi)

cms_events["mt_Wmunu"]      = getMT(cms_events.nMuTight,  cms_events.mupt, cms_events.muphi, cms_events.mupx, cms_events.mupy, cms_events.metpt, cms_events.metphi)
cms_events["mt_Wmunu0"]     = ak.firsts(cms_events.mt_Wmunu)

cms_events["mt_Wenu"]      = getMT(cms_events.nEleTight,  cms_events.elept, cms_events.elephi, cms_events.elepx, cms_events.elepy, cms_events.metpt, cms_events.metphi)
cms_events["mt_Wenu0"]     = ak.firsts(cms_events.mt_Wenu)


cms_events["jet_sel_loose"] = (cms_events.jetpt > 30.0  ) & (numpy.abs(cms_events.jeteta)<2.5)
cms_events["jet_sel_tight"] = (cms_events.jetpt > 50.0  ) & (numpy.abs(cms_events.jeteta)<2.5) 
#cms_events["jet_sel_b"]     = (cms_events.jetcsv > 0.6321) & (numpy.abs(cms_events.jeteta)<2.4)
cms_events["jet_sel_b"]     = (cms_events.jetcsv[cms_events.jet_sel_loose==True] > 0.6321) & (numpy.abs(cms_events.jeteta[cms_events.jet_sel_loose==True])<2.4)

#cms_events["jet_sel_tight0"] = ak.Array(getN(cms_events.jet_sel_tight,0))
cms_events["jet_sel_tight0"] = ak.Array(getN(cms_events.jet_sel_tight[cms_events.jet_sel_loose==True],0))
cms_events["jet_sel_b_0"]   = ak.Array(getN(cms_events.jet_sel_b,0))
cms_events["jet_sel_b_1"]   = ak.Array(getN(cms_events.jet_sel_b,1))

cms_events["nJetLoose"] = ak.sum(cms_events.jet_sel_loose,axis=-1)
cms_events["nJetTight"] = ak.sum(cms_events.jet_sel_tight,axis=-1)
cms_events["nJetb"] = ak.sum(cms_events.jet_sel_b,axis=-1)

cms_events["dphi_jet_met"] = DeltaPhi(cms_events.jetphi[cms_events.jet_sel_loose==True], cms_events.metphi)
cms_events["min_dphi_jet_met"] = ak.min(cms_events.dphi_jet_met, axis=-1)

#--------------------------------------------------------------------------------------------------
## W --> lepton + nu 
#--------------------------------------------------------------------------------------------------

mask_wmunu1b = ( (cms_events.nEleLoose==0 ) &
                 (cms_events.npho==0) &
                 (cms_events.ntau==0) &
                 (cms_events.nMuLoose==1) &
                 (cms_events.nMuTight==1) &
                 (cms_events.recoil_Wmunu0>200.) & 
                 (cms_events.min_dphi_jet_met > 0.5) &
                 (cms_events.nJetLoose==1 ) & 
                 (cms_events.nJetTight==1 ) &
                 (cms_events.nJetb ==1 ) & 
                 (cms_events.mt_Wmunu0 >0 ) & (cms_events.mt_Wmunu0 < 160 )
             )

cms_events["mask_wmunu1b"] = mask_wmunu1b
wm = cms_events[mask_wmunu1b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]


mask_wmunu2b = ( (cms_events.nEleLoose==0 ) &
                 (cms_events.npho==0) &
                 (cms_events.ntau==0) &
                 (cms_events.nMuLoose==1) &
                 (cms_events.nMuTight==1) &
                 (cms_events.recoil_Wmunu0>200.) & 
                 (cms_events.min_dphi_jet_met > 0.5) &
                 (cms_events.nJetLoose==2 ) & 
                 (cms_events.nJetTight>=1 ) &
                 (cms_events.nJetb ==2 ) & 
                 (cms_events.mt_Wmunu0 >0 ) & (cms_events.mt_Wmunu0 < 160 )
             )


cms_events["mask_wmunu2b"] = mask_wmunu2b
wm = cms_events[mask_wmunu2b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]


mask_wenu1b = ( (cms_events.nEleTight==1 ) &
                (cms_events.nEleLoose==1 ) &
                (cms_events.npho==0) &
                (cms_events.ntau==0) &
                (cms_events.nMuLoose==0) &
                (cms_events.recoil_Wenu0>200.) & 
                (cms_events.min_dphi_jet_met > 0.5) &
                (cms_events.nJetLoose==1 ) & 
                (cms_events.nJetTight==1 ) &
                (cms_events.nJetb ==1 ) & 
                (cms_events.mt_Wenu0 >0 ) & (cms_events.mt_Wenu0 < 160 )
            )

cms_events["mask_wenu1b"] = mask_wenu1b
wm = cms_events[mask_wenu1b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]


mask_wenu2b = ( (cms_events.nEleTight==1 ) &
                (cms_events.nEleLoose==1 ) &
                (cms_events.npho==0) &
                (cms_events.ntau==0) &
                (cms_events.nMuLoose==0) &
                (cms_events.recoil_Wenu0>200.) & 
                (cms_events.min_dphi_jet_met > 0.5) &
                (cms_events.nJetLoose==2 ) & 
                (cms_events.nJetTight>=1 ) &
                (cms_events.nJetb ==2 ) & 
                (cms_events.mt_Wenu0 >0 ) & (cms_events.mt_Wenu0 < 160 )
            )

cms_events["mask_wenu2b"] = mask_wenu2b
wm = cms_events[mask_wenu2b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]
len(a)
#--------------------------------------------------------------------------------------------------
## Top --> semi-leptonic
#--------------------------------------------------------------------------------------------------
mask_topmunu1b = ( (cms_events.nEleLoose==0 ) &
                   (cms_events.npho==0) &
                   (cms_events.ntau==0) &
                   (cms_events.nMuLoose==1) &
                   (cms_events.nMuTight==1) &
                   (cms_events.recoil_Wmunu0>200.) &
                   (cms_events.min_dphi_jet_met > 0.5) &
                   ((cms_events.nJetLoose>=2 ) & (cms_events.nJetLoose<=600 )) &
                   (cms_events.nJetTight>=1 ) &
                   (cms_events.nJetb ==1 ) &
                   (cms_events.jet_sel_tight0 == True) &
                   (cms_events.jet_sel_b_0 == True) &
                   (cms_events.mt_Wmunu0 >0 ) & (cms_events.mt_Wmunu0 < 160 )
               )
cms_events["mask_topmunu1b"] = mask_topmunu1b
wm = cms_events[mask_topmunu1b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]
len(a)

mask_topmunu2b = ( (cms_events.nEleLoose==0 ) &
                   (cms_events.npho==0) &
                   (cms_events.ntau==0) &
                   (cms_events.nMuLoose==1) &
                   (cms_events.nMuTight==1) &
                   (cms_events.recoil_Wmunu0>200.) & 
                   (cms_events.min_dphi_jet_met > 0.5) &
                   ((cms_events.nJetLoose>=3 ) & (cms_events.nJetLoose<=6 )) &
                   (cms_events.nJetTight>=1 ) &
                   ((cms_events.nJetb ==2 ) & (cms_events.nJetLoose<=600 )) &
                   (cms_events.jet_sel_tight0==True) & 
                   (cms_events.jet_sel_b_0 == True) &
                   (cms_events.jet_sel_b_1 == True) &
                   (cms_events.mt_Wmunu0 >0 ) & (cms_events.mt_Wmunu0 < 160 )
               )

cms_events["mask_topmunu2b"] = mask_topmunu2b

wm = cms_events[mask_topmunu2b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]
len(a)

mask_topenu1b = ( (cms_events.nEleTight==1 ) &
                  (cms_events.nEleLoose==1 ) &
                  (cms_events.npho==0) &
                  (cms_events.ntau==0) &
                  (cms_events.nMuLoose==0) &
                  (cms_events.recoil_Wenu0>200.) & 
                  (cms_events.min_dphi_jet_met > 0.5) &
                  ((cms_events.nJetLoose>=2 ) & (cms_events.nJetLoose<=60 )) &
                  (cms_events.nJetTight>=1 ) &
                  (cms_events.nJetb ==1 ) &
                  (cms_events.jet_sel_tight0 == True) &
                  (cms_events.jet_sel_b_0 == True) &
                  (cms_events.mt_Wenu0 >0 ) & (cms_events.mt_Wenu0 < 160 )
              )

cms_events["mask_topenu1b"] = mask_topenu1b
wm = cms_events[mask_topenu1b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]
len(a)

mask_topenu2b = ( (cms_events.nEleTight==1 ) &
                  (cms_events.nEleLoose==1 ) &
                  (cms_events.npho==0) &
                  (cms_events.ntau==0) &
                  (cms_events.nMuLoose==0) &
                  (cms_events.recoil_Wenu0>200.) & 
                  (cms_events.min_dphi_jet_met > 0.5) &
                  ((cms_events.nJetLoose>=3 ) & (cms_events.nJetLoose<=60 )) &
                  (cms_events.nJetTight>=1 ) &
                  (cms_events.nJetb ==2 ) & 
                  (cms_events.jet_sel_tight0==True) &
                  (cms_events.jet_sel_b_0 == True) &
                  (cms_events.jet_sel_b_1 == True) &
                  (cms_events.mt_Wenu0 >0 ) & (cms_events.mt_Wenu0 < 160 )
              )


cms_events["mask_topenu2b"] = mask_topenu2b
wm = cms_events[mask_topenu2b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]
len(a)
#--------------------------------------------------------------------------------------------------
## Z --> ll
#--------------------------------------------------------------------------------------------------

mask_Zmumu1b = ( (cms_events.nEleLoose==0 ) &
                 (cms_events.npho==0) &
                 (cms_events.ntau==0) &
                 (cms_events.nMuLoose==2) &
                 (cms_events.nMuTight>=1) &
                 (cms_events.mu_sel_tight0==True) &
                 (cms_events.Zmumu_recoil>200.) & 
                 (cms_events.min_dphi_jet_met > 0.5) &
                 ((cms_events.nJetLoose>=1 ) &  (cms_events.nJetLoose<=2 )) & 
                 (cms_events.nJetTight>=1 ) &
                 (cms_events.nJetb ==1 ) & 
                 (cms_events.jet_sel_tight0 == True) &
                 (cms_events.jet_sel_b_0 == True) &
                 ( (cms_events.mu_q0 * cms_events.mu_q1 ) < 0) &
                 (cms_events.Zmumu_mass >=70 ) & (cms_events.Zmumu_mass <= 110 )
             )


cms_events["mask_Zmumu1b"] = mask_Zmumu1b

wm = cms_events[mask_Zmumu1b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]


mask_Zmumu2b = ( (cms_events.nEleLoose==0 ) &
                 (cms_events.npho==0) &
                 (cms_events.ntau==0) &
                 (cms_events.nMuLoose==2) &
                 (cms_events.nMuTight>=1) &
                 (cms_events.mu_sel_tight0==True) &
                 (cms_events.Zmumu_recoil>200.) & 
                 (cms_events.min_dphi_jet_met > 0.5) &
                 ((cms_events.nJetLoose>=2 ) &  (cms_events.nJetLoose<=3 )) & 
                 (cms_events.nJetTight>=1 ) &
                 (cms_events.nJetb ==2 ) & 
                 (cms_events.jet_sel_tight0==True) &
                 (cms_events.jet_sel_b_0 == True) &
                 (cms_events.jet_sel_b_1 == True) &
                 ( (cms_events.mu_q0 * cms_events.mu_q1 ) < 0) &
                 (cms_events.Zmumu_mass >=70 ) & (cms_events.Zmumu_mass <= 110 )
             )

cms_events["mask_Zmumu2b"] = mask_Zmumu2b
wm = cms_events[mask_Zmumu2b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]




mask_Zee1b = ( (cms_events.nEleLoose==2 ) &
               (cms_events.nEleTight>=1 ) &
               (cms_events.ele_sel_tight0==True)&
               (cms_events.npho==0) &
               (cms_events.ntau==0) &
               (cms_events.nMuLoose==0) &
               (cms_events.Zee_recoil>200.) & 
               (cms_events.min_dphi_jet_met > 0.5) &
               ((cms_events.nJetLoose>=1 ) &  (cms_events.nJetLoose<=2 )) & 
               (cms_events.nJetTight>=1 ) &
               (cms_events.nJetb ==1 ) & 
               (cms_events.jet_sel_tight0 == True) &
               (cms_events.jet_sel_b_0 == True) &
               ( (cms_events.ele_q0 * cms_events.ele_q1 ) < 0) &
               (cms_events.Zee_mass >=70 ) & (cms_events.Zee_mass <= 110 )
             )


cms_events["mask_Zee1b"] = mask_Zee1b

wm = cms_events[mask_Zee1b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]



mask_Zee2b = ( (cms_events.nEleLoose==2 ) &
               (cms_events.nEleTight>=1 ) &
               (cms_events.ele_sel_tight0==True) &
               (cms_events.npho==0) &
               (cms_events.ntau==0) &
               (cms_events.nMuLoose==0) &
               (cms_events.Zee_recoil>200.) &
               (cms_events.min_dphi_jet_met > 0.5) &
               ((cms_events.nJetLoose>=2 ) &  (cms_events.nJetLoose<=3 )) &
               (cms_events.nJetTight>=1 ) &
               (cms_events.nJetb ==2 ) &
               (cms_events.jet_sel_tight0==True) &
               (cms_events.jet_sel_b_0 == True) &
               (cms_events.jet_sel_b_1 == True) &
               ( (cms_events.ele_q0 * cms_events.ele_q1 ) < 0) &
               (cms_events.Zee_mass >=70 ) & (cms_events.Zee_mass <= 110 )
           )

cms_events["mask_Zee2b"] = mask_Zee2b

wm = cms_events[mask_Zee2b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]


mask_SR1b =  ( (cms_events.metpt > 200) &
               (cms_events.nEleLoose==0 ) &
               (cms_events.npho==0) &
               (cms_events.ntau==0) &
               (cms_events.nMuLoose==0) &
               (cms_events.min_dphi_jet_met > 0.5) &
               ((cms_events.nJetLoose>=1 ) &  (cms_events.nJetLoose<=2 )) & 
               (cms_events.nJetTight>=1 ) &
               (cms_events.nJetb ==1 ) & 
               (cms_events.jet_sel_tight0 == True) &
               (cms_events.jet_sel_b_0 == True) 
           )

cms_events["mask_SR1b"] = mask_SR1b
wm = cms_events[mask_SR1b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]




mask_SR2b =  ( (cms_events.metpt > 200) &
               (cms_events.nEleLoose==0 ) &
               (cms_events.npho==0) &
               (cms_events.ntau==0) &
               (cms_events.nMuLoose==0) &
               (cms_events.min_dphi_jet_met > 0.5) &
               ((cms_events.nJetLoose>=2 ) &  (cms_events.nJetLoose<=3 )) & 
               (cms_events.nJetTight>=1 ) &
               (cms_events.nJetb ==2 ) & 
               (cms_events.jet_sel_tight0 == True) &
               (cms_events.jet_sel_b_0 == True) &
               (cms_events.jet_sel_b_1 == True) 

           )

cms_events["mask_SR2b"] = mask_SR2b
wm = cms_events[mask_SR2b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]


end = time.clock()
print("%.4gs" % (end-start))
