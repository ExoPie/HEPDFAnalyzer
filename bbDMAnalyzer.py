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


from utils import getpt, geteta, getphi, getrecoil, DeltaPhi, getMT, getRecoilPhi, getrecoil1, getN, weight_W1mu_


print ("opening rootfile")
mycache = uproot4.LRUArrayCache("10 MB")

#tree_ = uproot4.open(inputfile)["outTree"].arrays(array_cache=mycache)

#tree_ = uproot4.open("Merged_WJetsInclusiveSkim.root")["outTree"].arrays(array_cache=mycache)
tree_ = uproot4.open("Merged_WjetsInclusive_photveto.root")["outTree"].arrays(array_cache=mycache)
#tree_ = uproot4.open("/eos/cms/store/group/phys_exotica/bbMET/2016_SkimmedFiles/skim_setup_2016_v16_07-00/crab_DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_200918_215129_0000_0.root")["outTree"].arrays(array_cache=mycache)

#tree_ = uproot4.open("/eos/cms/store/group/phys_exotica/bbMET/2016_SkimmedFiles/skim_setup_2016_v16_07-00/crab_ttHTobb_M125_13TeV_powheg_pythia8_200918_215950_0000_0.root")["outTree"].arrays(array_cache=mycache)



#print ((tree_))

cms_events = ak.zip({   "run":tree_["st_runId"],"lumi":tree_["st_lumiSection"], "event": tree_["st_eventId"],
                        "jetpx":tree_["st_THINjetPx"], "jetpy":tree_["st_THINjetPy"], "jetpz":tree_["st_THINjetPz"], "jete":tree_["st_THINjetEnergy"],
                        "jetpt": getpt(tree_["st_THINjetPx"], tree_["st_THINjetPy"]), "jeteta":geteta(tree_["st_THINjetPx"], tree_["st_THINjetPy"], tree_["st_THINjetPz"]), 
                        "jetphi":getphi(tree_["st_THINjetPx"], tree_["st_THINjetPy"]), "jetcsv": tree_["st_THINjetDeepCSV"],
                        "jetflav":tree_["st_THINjetHadronFlavor"],
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
                        "phopt": getpt(tree_["st_phoPx"], tree_["st_phoPy"]), "phoeta":geteta(tree_["st_phoPx"], tree_["st_phoPy"], tree_["st_phoPz"]),
                        "nTrueInt":tree_["st_pu_nTrueInt"], "nPUVert":tree_["st_pu_nPUVert"]
                    },
                    depth_limit=1)



out_events = ak.zip({"run":tree_["st_runId"],"lumi":tree_["st_lumiSection"], "event": tree_["st_eventId"] },
                    depth_limit=1)


print ("event loading done")
print ("# of events: ",len(cms_events))

## add more columns/properties to the event 
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

cms_events["jetptTight"]  = cms_events.jetpt[cms_events.jet_sel_tight==True] 
cms_events["jetetaTight"] = cms_events.jeteta[cms_events.jet_sel_tight==True] 
cms_events["jetphiTight"] = cms_events.jetphi[cms_events.jet_sel_tight==True] 

cms_events["jetptLoose"] = cms_events.jetpt[cms_events.jet_sel_loose==True] 
cms_events["jetetaLoose"] = cms_events.jeteta[cms_events.jet_sel_loose==True] 
cms_events["jetphiLoose"] = cms_events.jetphi[cms_events.jet_sel_loose==True] 


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
'''
wm = cms_events[mask_wmunu1b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]
'''

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
'''
wm = cms_events[mask_wmunu2b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]
'''

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
'''
wm = cms_events[mask_wenu1b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]
'''

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
'''
wm = cms_events[mask_wenu2b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]
len(a)
'''
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
'''
wm = cms_events[mask_topmunu1b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]
len(a)
'''
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
'''
wm = cms_events[mask_topmunu2b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]
len(a)
'''
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
'''
wm = cms_events[mask_topenu1b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]
len(a)
'''
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
'''
wm = cms_events[mask_topenu2b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]
len(a)
'''
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
'''
wm = cms_events[mask_Zmumu1b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]
'''

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
'''
wm = cms_events[mask_Zmumu2b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]
'''



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

'''
wm = cms_events[mask_Zee1b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]
'''


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
'''
wm = cms_events[mask_Zee2b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]
'''

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
'''
wm = cms_events[mask_SR1b]
event = ak.to_list(wm.event)
a=[i for i in event if i is not None]
'''



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
'''
wm = cms_events.event[mask_SR2b]
wm[~ak.is_none(wm)]
'''

###############
out_events["mu_sel_tight0"] = cms_events["mu_sel_tight0"]
out_events["nMuTight"]      = cms_events["nMuTight"]
out_events["nMuLoose"]      = cms_events["nMuLoose"]
out_events["mu_q0"]         =   cms_events["mu_q0"] 
out_events["mu_q1"]         =   cms_events["mu_q1"]
out_events["mupt0"] = ak.Array(getN(cms_events.mupt,0))
out_events["mupt1"] = ak.Array(getN(cms_events.mupt,1))
out_events["mueta0"] = ak.Array(getN(cms_events.mueta,0))
out_events["mueta1"] = ak.Array(getN(cms_events.mueta,1))
out_events["muphi0"] = ak.Array(getN(cms_events.muphi,0))
out_events["muphi1"] = ak.Array(getN(cms_events.muphi,1))

out_events["ele_sel_tight0"] = cms_events["ele_sel_tight0"]
out_events["nEleTight"] = cms_events["nEleTight"]
out_events["nEleLoose"] = cms_events["nEleLoose"]
out_events["ele_q0"] = cms_events["ele_q0"]
out_events["ele_q1"] = cms_events["ele_q1"]
out_events["elept0"] = ak.Array(getN(cms_events.elept,0))
out_events["elept1"] = ak.Array(getN(cms_events.elept,1))
out_events["eleeta0"] = ak.Array(getN(cms_events.eleeta,0))
out_events["eleeta1"] = ak.Array(getN(cms_events.eleeta,1))
out_events["elephi0"] = ak.Array(getN(cms_events.elephi,0))
out_events["elephi1"] = ak.Array(getN(cms_events.elephi,1))

out_events["recoil_Wmunu0"] = cms_events["recoil_Wmunu0"] 
out_events["recoil_Wenu0"]  = cms_events["recoil_Wenu0"]
out_events["recoil_WmunuPhi0"] =cms_events["recoil_WmunuPhi0"]
out_events["recoil_WenuPhi0"]  = cms_events["recoil_WenuPhi0"] 
out_events["mt_Wmunu0"] = cms_events["mt_Wmunu0"]
out_events["mt_Wenu0"] = cms_events["mt_Wenu0"]

out_events["Zee_mass"] = cms_events["Zee_mass"]
out_events["Zee_pt"] = cms_events["Zee_pt"]
out_events["Zee_recoil"]   = cms_events["Zee_recoil"]  
out_events["Zmumu_mass"] = cms_events["Zmumu_mass"]
out_events["Zmumu_pt"] = cms_events["Zmumu_pt"]
out_events["Zmumu_recoil"] = cms_events["Zmumu_recoil"]



out_events["nJetLoose"] = cms_events["nJetLoose"]
out_events["nJetTight"] = cms_events["nJetTight"]
out_events["nJetb"] = cms_events["nJetb"]
out_events["min_dphi_jet_met"] = cms_events["min_dphi_jet_met"]
cms_events["jet_sel_tight0"] = cms_events["jet_sel_tight0"] 
cms_events["jet_sel_b_0"] = cms_events["jet_sel_b_0"]  
cms_events["jet_sel_b_1"] = cms_events["jet_sel_b_1"]  






out_events["jetpt0"] = ak.Array(getN(cms_events.jetptTight,0))
out_events["jetpt1"] = ak.Array(getN(cms_events.jetptLoose,1))
out_events["jetpt2"] = ak.Array(getN(cms_events.jetptLoose,2))
out_events["jetpt3"] = ak.Array(getN(cms_events.jetptLoose,3))
out_events["jetpt4"] = ak.Array(getN(cms_events.jetptLoose,4))
out_events["jetpt5"] = ak.Array(getN(cms_events.jetptLoose,5))
out_events["jetpt6"] = ak.Array(getN(cms_events.jetptLoose,6))

out_events["jeteta0"] = ak.Array(getN(cms_events.jetetaTight,0))
out_events["jeteta1"] = ak.Array(getN(cms_events.jetetaLoose,1))
out_events["jeteta2"] = ak.Array(getN(cms_events.jetetaLoose,2))
out_events["jeteta3"] = ak.Array(getN(cms_events.jetetaLoose,3))
out_events["jeteta4"] = ak.Array(getN(cms_events.jetetaLoose,4))
out_events["jeteta5"] = ak.Array(getN(cms_events.jetetaLoose,5))
out_events["jeteta6"] = ak.Array(getN(cms_events.jetetaLoose,6))

out_events["jetphi0"] = ak.Array(getN(cms_events.jetphiTight,0))
out_events["jetphi1"] = ak.Array(getN(cms_events.jetphiLoose,1))
out_events["jetphi2"] = ak.Array(getN(cms_events.jetphiLoose,2))

out_events["jetflav0"] = ak.Array(getN(cms_events.jetflav[cms_events.jet_sel_tight==True],0))
out_events["jetflav1"] = ak.Array(getN(cms_events.jetflav[cms_events.jet_sel_loose==True],1))
out_events["jetflav2"] = ak.Array(getN(cms_events.jetflav[cms_events.jet_sel_loose==True],2))
out_events["jetflav3"] = ak.Array(getN(cms_events.jetflav[cms_events.jet_sel_loose==True],3))
out_events["jetflav4"] = ak.Array(getN(cms_events.jetflav[cms_events.jet_sel_loose==True],4))
out_events["jetflav5"] = ak.Array(getN(cms_events.jetflav[cms_events.jet_sel_loose==True],5))
out_events["jetflav6"] = ak.Array(getN(cms_events.jetflav[cms_events.jet_sel_loose==True],6))


out_events["csv0"]    = ak.Array(getN(cms_events.jetcsv[cms_events.jet_sel_tight==True],0))
out_events["csv1"]    = ak.Array(getN(cms_events.jetcsv[cms_events.jet_sel_loose==True],1))
out_events["csv2"]    = ak.Array(getN(cms_events.jetcsv[cms_events.jet_sel_loose==True],2))
out_events["csv3"]    = ak.Array(getN(cms_events.jetcsv[cms_events.jet_sel_loose==True],3))


out_events["mask_SR2b"] = cms_events["mask_SR2b"]
out_events["mask_SR1b"] = cms_events["mask_SR1b"]
out_events["mask_Zee2b"] =cms_events["mask_Zee2b"] 
out_events["mask_Zee1b"] =cms_events["mask_Zee1b"] 
out_events["mask_Zmumu2b"] =cms_events["mask_Zmumu2b"] 
out_events["mask_Zmumu1b"] =cms_events["mask_Zmumu1b"] 
out_events["mask_topenu2b"] = cms_events["mask_topenu2b"]
out_events["mask_topenu1b"] = cms_events["mask_topenu1b"]
out_events["mask_topmunu2b"] = cms_events["mask_topmunu2b"]
out_events["mask_topmunu1b"] = cms_events["mask_topmunu1b"]
out_events["mask_wenu1b"] = cms_events["mask_wenu1b"]
out_events["mask_wenu2b"] = cms_events["mask_wenu2b"]
out_events["mask_wmunu1b"] = cms_events["mask_wmunu1b"]
out_events["mask_wmunu2b"] = cms_events["mask_wmunu2b"]




####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  
####  ####   the following code to read the SFs can be moved in another python file ####  ####  ####  
####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  
## btagging SFs 
from coffea.btag_tools import BTagScaleFactor
btag_sf = BTagScaleFactor("data/DeepCSV_2016LegacySF_V1.csv.gz", "medium")
out_events["btagsf0"] = btag_sf.eval("central", out_events.jetflav0, abs(out_events.jeteta0), out_events.jetpt0)
out_events["btagsf1"] = btag_sf.eval("central", out_events.jetflav1, abs(out_events.jeteta1), out_events.jetpt1)
out_events["btagsf2"] = btag_sf.eval("central", out_events.jetflav2, abs(out_events.jeteta2), out_events.jetpt2)
out_events["btagsf3"] = btag_sf.eval("central", out_events.jetflav3, abs(out_events.jeteta3), out_events.jetpt3)
out_events["btagsf4"] = btag_sf.eval("central", out_events.jetflav4, abs(out_events.jeteta4), out_events.jetpt4)
out_events["btagsf5"] = btag_sf.eval("central", out_events.jetflav5, abs(out_events.jeteta5), out_events.jetpt5)
out_events["btagsf6"] = btag_sf.eval("central", out_events.jetflav6, abs(out_events.jeteta6), out_events.jetpt6)


## ele tight ID weights 
from coffea.lookup_tools import extractor
ext = extractor()

ext.add_weight_sets(["EGamma_SF2D_T EGamma_SF2D data/2016LegacyReReco_ElectronTight_Fall17V2.root",
                     "EGamma_SF2D_T_err EGamma_SF2D_error data/2016LegacyReReco_ElectronTight_Fall17V2.root"])

ext.add_weight_sets(["EGamma_SF2D_L EGamma_SF2D data/2016LegacyReReco_ElectronLoose_Fall17V2.root",
                     "EGamma_SF2D_L_err EGamma_SF2D_error data/2016LegacyReReco_ElectronLoose_Fall17V2.root"])

ext.add_weight_sets(["EGamma_SF2D_Trig hEffEtaPt data/electron_Trigger_eleTrig.root",
                     "EGamma_SF2D_Trig_err hEffEtaPt_error data/electron_Trigger_eleTrig.root" ])

ext.add_weight_sets(["EGamma_SF2D_Reco EGamma_SF2D data/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root",
                     "EGamma_SF2D_Reco_err EGamma_SF2D_error data/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root"])

ext.add_weight_sets(["EGamma_SF2D_Reco_lowpt EGamma_SF2D data/EGM2D_BtoH_low_RecoSF_Legacy2016.root",
                     "EGamma_SF2D_Reco_lowpt_err EGamma_SF2D_error data/EGM2D_BtoH_low_RecoSF_Legacy2016.root"])


ext.add_weight_sets(["met_trig hden_monojet_recoil_clone_passed data/metTriggerEfficiency_zmm_recoil_monojet_TH1F.root"])



ext.add_weight_sets(["muon_lowpt_BCDEF_LooseID NUM_LooseID_DEN_genTracks_pt_abseta data/Muon_low-pT_RunBCDEF_SF_ID.root", 
                     "muon_lowpt_BCDEF_TightID NUM_TightID_DEN_genTracks_pt_abseta data/Muon_low-pT_RunBCDEF_SF_ID.root",
                     "muon_lowpt_GH_LooseID NUM_LooseID_DEN_genTracks_pt_abseta data/Muon_low-pT_RunGH_SF_ID.root", 
                     "muon_lowpt_GH_TightID NUM_TightID_DEN_genTracks_pt_abseta data/Muon_low-pT_RunGH_SF_ID.root",
                     
                     "muon_highpt_BCDEF_LooseISO NUM_LooseRelIso_DEN_LooseID_eta_pt data/Muon_RunBCDEF_SF_ISO.root",
                     "muon_highpt_BCDEF_TightISO NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt data/Muon_RunBCDEF_SF_ISO.root",
                     "muon_highpt_GH_LooseISO NUM_LooseRelIso_DEN_LooseID_eta_pt data/Muon_RunGH_SF_ISO.root",
                     "muon_highpt_GH_TightISO NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt data/Muon_RunGH_SF_ISO.root",
                     
                     "muon_highpt_BCDEF_LooseID NUM_LooseID_DEN_genTracks_eta_pt data/Muon_RunBCDEF_SF_ID.root",
                     "muon_highpt_BCDEF_TightID NUM_TightID_DEN_genTracks_eta_pt data/Muon_RunBCDEF_SF_ID.root",
                     "muon_highpt_GH_LooseID NUM_LooseID_DEN_genTracks_eta_pt data/Muon_RunGH_SF_ID.root",
                     "muon_highpt_GH_TightID NUM_TightID_DEN_genTracks_eta_pt data/Muon_RunGH_SF_ID.root"])
                     


ext.add_weight_sets(["pu_weight puweight data/PU_Reweight_2016.root",
                     "pu_weight_Up puweight_Up data/PU_Reweight_2016.root",
                     "pu_weight_Down puweight_Down data/PU_Reweight_2016.root"])


ext.add_weight_sets(["btag_eff_lwp efficiency_btag_lwp data/bTagEffs_2016.root",
                     "btag_eff_mwp efficiency_btag_mwp data/bTagEffs_2016.root",
                     "ctag_eff_lwp efficiency_ctag_lwp data/bTagEffs_2016.root",
                     "ctag_eff_mwp efficiency_ctag_mwp data/bTagEffs_2016.root",
                     "ltag_eff_lwp efficiency_lighttag_lwp data/bTagEffs_2016.root",
                     "ltag_eff_mwp efficiency_lighttag_mwp data/bTagEffs_2016.root" ])
ext.finalize()
evaluator = ext.make_evaluator()

####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  #### 
####  ####   upto this should be in another python file ####  ####  ####  ####  ####  ####  ####  ####  ####  
####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  

out_events["btag_eff_lwp_0"] = evaluator["btag_eff_lwp"](out_events.jeteta0, out_events.jetpt0)
out_events["btag_eff_lwp_1"] = evaluator["btag_eff_lwp"](out_events.jeteta1, out_events.jetpt1) 

out_events["ctag_eff_lwp_0"] = evaluator["ctag_eff_lwp"](out_events.jeteta0, out_events.jetpt0)
out_events["ctag_eff_lwp_1"] = evaluator["ctag_eff_lwp"](out_events.jeteta1, out_events.jetpt1) 

out_events["ltag_eff_lwp_0"] = evaluator["ltag_eff_lwp"](out_events.jeteta0, out_events.jetpt0)
out_events["ltag_eff_lwp_1"] = evaluator["ltag_eff_lwp"](out_events.jeteta1, out_events.jetpt1) 

out_events["btag_eff_mwp_0"] = evaluator["btag_eff_mwp"](out_events.jeteta0, out_events.jetpt0)
out_events["btag_eff_mwp_1"] = evaluator["btag_eff_mwp"](out_events.jeteta1, out_events.jetpt1) 

out_events["ctag_eff_mwp_0"] = evaluator["ctag_eff_mwp"](out_events.jeteta0, out_events.jetpt0)
out_events["ctag_eff_mwp_1"] = evaluator["ctag_eff_mwp"](out_events.jeteta1, out_events.jetpt1) 

out_events["ltag_eff_mwp_0"] = evaluator["ltag_eff_mwp"](out_events.jeteta0, out_events.jetpt0)
out_events["ltag_eff_mwp_1"] = evaluator["ltag_eff_mwp"](out_events.jeteta1, out_events.jetpt1) 


out_events["eleTightSF0"]  = evaluator["EGamma_SF2D_T"](out_events.eleeta0, out_events.elept0)
out_events["eleLooseSF1"]  = evaluator["EGamma_SF2D_L"](out_events.eleeta1, out_events.elept1)
out_events["eleTrigSF0"]   = evaluator["EGamma_SF2D_Trig"](out_events.eleeta0, out_events.elept0)
out_events["eleRecoSF0"]   = evaluator["EGamma_SF2D_Reco"](out_events.eleeta0, out_events.elept0)

eleRecoSF1_hi   = evaluator["EGamma_SF2D_Reco"](out_events.eleeta1, out_events.elept1)
eleRecoSF1_lo   = evaluator["EGamma_SF2D_Reco_lowpt"](out_events.eleeta1, out_events.elept1)

eleRecoSF1_hi_ = ak.fill_none( ak.mask(  eleRecoSF1_hi , out_events.elept1>20.) ,0 )
eleRecoSF1_lo_ = ak.fill_none( ak.mask(  eleRecoSF1_lo , out_events.elept1>20.) ,0 )
out_events["eleRecoSF1"] = eleRecoSF1_hi_ + eleRecoSF1_lo_


bcdef_lumi  =19.554725529
gh_lumi = 16.224846377
total_lumi = bcdef_lumi + gh_lumi

##--------low pt Loose
muonLooseIDSF_lowpt1  = ( (bcdef_lumi*evaluator["muon_lowpt_BCDEF_LooseID"](out_events.mupt1, abs(out_events.mueta1) ))  +  
                                 (gh_lumi*evaluator["muon_lowpt_GH_LooseID"](out_events.mupt1, abs(out_events.mueta1) )) 
                             )/total_lumi


##----------- medium pt Loose
muonLooseIDSF1  = ( (bcdef_lumi*evaluator["muon_highpt_BCDEF_LooseID"](out_events.mueta1, out_events.mupt1))  +  
                   (gh_lumi*evaluator["muon_highpt_GH_LooseID"](out_events.mueta1, out_events.mupt1)) 
               )/total_lumi


muonLooseISOSF1  = ( (bcdef_lumi*evaluator["muon_highpt_BCDEF_LooseISO"](out_events.mueta1, out_events.mupt1))  +  
                    (gh_lumi*evaluator["muon_highpt_GH_LooseISO"](out_events.mueta1, out_events.mupt1)) 
                )/total_lumi

muon_loose_ID_low_SF_1  = ak.fill_none( ak.mask(  muonLooseIDSF_lowpt1 , out_events.mupt1<20.) ,0 )
muon_loose_ID_high_SF_1 = ak.fill_none( ak.mask(  muonLooseIDSF1       , out_events.mupt1>20.) ,0 )

muon_loose_ID_SF_1 = muon_loose_ID_low_SF_1 + muon_loose_ID_high_SF_1

out_events["muLooseSF1"]  = muon_loose_ID_SF_1 * muonLooseISOSF1 


##------------medium pt tight
muonTightIDSF0  = ( (bcdef_lumi*evaluator["muon_highpt_BCDEF_TightID"](out_events.mueta0, out_events.mupt0))  +  
                    (gh_lumi*evaluator["muon_highpt_GH_TightID"](out_events.mueta0, out_events.mupt0)) 
                )/total_lumi

muonTightISOSF0  = ( (bcdef_lumi*evaluator["muon_highpt_BCDEF_TightISO"](out_events.mueta0, out_events.mupt0))  +  
                                  (gh_lumi*evaluator["muon_highpt_GH_TightISO"](out_events.mueta0, out_events.mupt0)) 
                             )/total_lumi

out_events["muTightSF0"]  = muonTightIDSF0 * muonTightISOSF0

out_events["puweight"] = evaluator["pu_weight"](cms_events.nTrueInt)

out_events["mettrigWeight"] = evaluator["met_trig"](cms_events.metpt)
out_events["recoilWmunutrigWeight"] = evaluator["met_trig"](cms_events.recoil_Wmunu0)
out_events["recoilWenutrigWeight"] = evaluator["met_trig"](cms_events.recoil_Wenu0)
out_events["recoilZmumutrigWeight"] = evaluator["met_trig"](cms_events.Zmumu_recoil)
out_events["recoilZeetrigWeight"] = evaluator["met_trig"](cms_events.Zee_recoil)



## Fill weights for each CR so that we don't need to worry later 

out_events["weight_SR2b"]        = out_events.puweight * out_events.mettrigWeight 
out_events["weight_SR1b"]        = out_events.puweight * out_events.mettrigWeight 

out_events["weight_Zee2b"]       = out_events.puweight * out_events.eleTrigSF0 
out_events["weight_Zee1b"]       = out_events.puweight * out_events.eleTrigSF0 

out_events["weight_Zmumu2b"]     = out_events.puweight * out_events.recoilZmumutrigWeight 
out_events["weight_Zmumu1b"]     = out_events.puweight * out_events.recoilZmumutrigWeight 

out_events["weight_topenu2b"]    = out_events.puweight * out_events.eleTrigSF0 
out_events["weight_topenu1b"]    = out_events.puweight * out_events.eleTrigSF0 

out_events["weight_topmunu2b"]   = out_events.puweight * out_events.recoilWmunutrigWeight 
out_events["weight_topmunu1b"]   = out_events.puweight * out_events.recoilWmunutrigWeight 

out_events["weight_wenu1b"]      = out_events.puweight * out_events.eleTrigSF0 
out_events["weight_wenu2b"]      = out_events.puweight * out_events.eleTrigSF0 

out_events["weight_wmunu1b"]     = out_events.puweight * out_events.recoilWmunutrigWeight
out_events["weight_wmunu2b"]     = out_events.puweight * out_events.recoilWmunutrigWeight



## add EWK reweighting 
## add QCD reweighting 
## add missing weights 
## add up and down weights 
## get total events 
## get cross-section 
## 
## Fill Histograms 


print (ak.to_list(out_events[:10]))
ak.to_parquet(out_events, "analysis_wjets_allevents.parquet")


end = time.clock()
print("%.4gs" % (end-start))



'''
SF Links
-----------------
electrons: 
trigger: https://github.com/ExoPie/ExoPieUtils/blob/test_systematics/scalefactortools/data_2016/electron_Trigger_eleTrig.root
ele_reco: https://github.com/ExoPie/ExoPieUtils/blob/test_systematics/scalefactortools/data_2016/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root
ele_reco_lowpt: https://github.com/ExoPie/ExoPieUtils/blob/test_systematics/scalefactortools/data_2016/EGM2D_BtoH_low_RecoSF_Legacy2016.root
loose_id: https://github.com/ExoPie/ExoPieUtils/blob/test_systematics/scalefactortools/data_2016/2016LegacyReReco_ElectronLoose_Fall17V2.root
tight_id: https://github.com/ExoPie/ExoPieUtils/blob/test_systematics/scalefactortools/data_2016/2016LegacyReReco_ElectronTight_Fall17V2.root
Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Efficiencies_and_scale_factors

-----------------
trigger met: https://github.com/tiwariPC/ExoPieUtils/blob/test_systematics/scalefactortools/data_2016/metTriggerEfficiency_zmm_recoil_monojet_TH1F.root

------------------
muons: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs2016LegacyRereco


'''
