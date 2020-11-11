## add EWK reweighting 
## add QCD reweighting 
## add missing weights 



## add up and down weights 
## get total events 
## get cross-section 
## apply btag weights properly 
## no reweighting for data files 

import sys 
import uproot4 
import awkward1 as ak 
import numpy 
import math 
import time
from ROOT import TFile, TH1F
import copy 
import multiprocessing as mp
from utils import getpt, geteta, getphi, getrecoil, DeltaPhi, getMT, getRecoilPhi, getrecoil1, getN, VarToHist, SetHist, FileToList
from functools import partial

start = time.clock()
#trees_ = ['outTree']
inputfile=sys.argv[1]

def runOneFile(filename):
    #print ("filename: ", filename)
    #inputfile=filename
    outputfile = "output/"+inputfile.split("/")[-1]
    #outputfile = "tmp.root"
    mycache = uproot4.LRUArrayCache("1 MB")
    file_=uproot4.open(inputfile, num_workers=10)
    #print ("root file opened: ", filename)
    nevents=ak.to_list(file_["h_total_mcweight"].values())[2]
    #nevents = 1000000
    print ("histogram opened: ", nevents)
    
    #tree_ = uproot4.open(inputfile, num_workers=10)["outTree"].arrays(array_cache=mycache)
    tree_ = file_["outTree"].arrays(array_cache=mycache)
    print ("tree length", len(tree_))
    
    #tree_ = uproot4.open(inputfile)[trees[0]].arrays()
    #tree_ = uproot4.open(inputfile)["outTree"].arrays(array_cache=mycache)
    #tree_ = uproot4.open("Merged_WJetsInclusiveSkim.root")["outTree"].arrays(array_cache=mycache)
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
                            "nTrueInt":tree_["st_pu_nTrueInt"], "nPUVert":tree_["st_pu_nPUVert"],
                            "genpt":tree_["st_genParPt"]
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
    from regions import get_mask_wmunu1b, get_mask_wmunu2b, get_mask_wenu1b, get_mask_wenu2b, get_mask_topmunu1b, get_mask_topmunu2b, get_mask_topenu1b, get_mask_topenu2b, get_mask_Zmumu1b, get_mask_Zmumu2b, get_mask_Zee1b, get_mask_Zee2b, get_mask_SR1b, get_mask_SR2b 
    
    cms_events["mask_wmunu1b"]   = get_mask_wmunu1b  (cms_events)
    cms_events["mask_wmunu2b"]   = get_mask_wmunu2b  (cms_events)
    cms_events["mask_wenu1b"]    = get_mask_wenu1b   (cms_events)
    cms_events["mask_wenu2b"]    = get_mask_wenu2b   (cms_events)
    cms_events["mask_topmunu1b"] = get_mask_topmunu1b(cms_events)
    cms_events["mask_topmunu2b"] = get_mask_topmunu2b(cms_events)
    cms_events["mask_topenu1b"]  = get_mask_topenu1b (cms_events)
    cms_events["mask_topenu2b"]  = get_mask_topenu2b (cms_events)
    cms_events["mask_Zmumu1b"]   = get_mask_Zmumu1b  (cms_events)
    cms_events["mask_Zmumu2b"]   = get_mask_Zmumu2b  (cms_events)
    cms_events["mask_Zee1b"]     = get_mask_Zee1b    (cms_events)
    cms_events["mask_Zee2b"]     = get_mask_Zee2b    (cms_events)
    cms_events["mask_SR1b"]      = get_mask_SR1b     (cms_events)
    cms_events["mask_SR2b"]      = get_mask_SR2b     (cms_events)
    
    '''
    wm = cms_events.event[mask_SR2b]
    wm[~ak.is_none(wm)]
    '''
    
    ###############
    out_events["metpt"]             = cms_events["metpt"]
    out_events["metphi"]            = cms_events["metphi"]
    out_events["nTrueInt"]          = cms_events["nTrueInt"]
    out_events["nJetLoose"]         = cms_events["nJetLoose"]
                                    
    out_events["mu_sel_tight0"]     = cms_events["mu_sel_tight0"]
    out_events["nMuTight"]          = cms_events["nMuTight"]
    out_events["nMuLoose"]          = cms_events["nMuLoose"]
    out_events["mu_q0"]             = cms_events["mu_q0"] 
    out_events["mu_q1"]             = cms_events["mu_q1"]
    out_events["mupt0"]             = ak.Array(getN(cms_events.mupt,0))
    out_events["mupt1"]             = ak.Array(getN(cms_events.mupt,1))
    out_events["mueta0"]            = ak.Array(getN(cms_events.mueta,0))
    out_events["mueta1"]            = ak.Array(getN(cms_events.mueta,1))
    out_events["muphi0"]            = ak.Array(getN(cms_events.muphi,0))
    out_events["muphi1"]            = ak.Array(getN(cms_events.muphi,1))
                                    
    out_events["ele_sel_tight0"]    = cms_events["ele_sel_tight0"]
    out_events["nEleTight"]         = cms_events["nEleTight"]
    out_events["nEleLoose"]         = cms_events["nEleLoose"]
    out_events["ele_q0"]            = cms_events["ele_q0"]
    out_events["ele_q1"]            = cms_events["ele_q1"]
    out_events["elept0"]            = ak.Array(getN(cms_events.elept,0))
    out_events["elept1"]            = ak.Array(getN(cms_events.elept,1))
    out_events["eleeta0"]           = ak.Array(getN(cms_events.eleeta,0))
    out_events["eleeta1"]           = ak.Array(getN(cms_events.eleeta,1))
    out_events["elephi0"]           = ak.Array(getN(cms_events.elephi,0))
    out_events["elephi1"]           = ak.Array(getN(cms_events.elephi,1))
                                    
    out_events["recoil_Wmunu0"]     = cms_events["recoil_Wmunu0"] 
    out_events["recoil_Wenu0"]      = cms_events["recoil_Wenu0"]
    out_events["recoil_WmunuPhi0"]  = cms_events["recoil_WmunuPhi0"]
    out_events["recoil_WenuPhi0"]   = cms_events["recoil_WenuPhi0"] 
    out_events["mt_Wmunu0"]         = cms_events["mt_Wmunu0"]
    out_events["mt_Wenu0"]          = cms_events["mt_Wenu0"]
                                    
    out_events["Zee_mass"]          = cms_events["Zee_mass"]
    out_events["Zee_pt"]            = cms_events["Zee_pt"]
    out_events["Zee_recoil"]        = cms_events["Zee_recoil"]  
    out_events["Zmumu_mass"]        = cms_events["Zmumu_mass"]
    out_events["Zmumu_pt"]          = cms_events["Zmumu_pt"]
    out_events["Zmumu_recoil"]      = cms_events["Zmumu_recoil"]
                                    
    out_events["nJetLoose"]         = cms_events["nJetLoose"]
    out_events["nJetTight"]         = cms_events["nJetTight"]
    out_events["nJetb"]             = cms_events["nJetb"]
    out_events["min_dphi_jet_met"]  = cms_events["min_dphi_jet_met"]
    cms_events["jet_sel_tight0"]    = cms_events["jet_sel_tight0"] 
    cms_events["jet_sel_b_0"]       = cms_events["jet_sel_b_0"]  
    cms_events["jet_sel_b_1"]       = cms_events["jet_sel_b_1"]  
                                    
    out_events["jetpt0"]            = ak.Array(getN(cms_events.jetptTight,0))
    out_events["jetpt1"]            = ak.Array(getN(cms_events.jetptLoose,1))
    out_events["jetpt2"]            = ak.Array(getN(cms_events.jetptLoose,2))
    out_events["jetpt3"]            = ak.Array(getN(cms_events.jetptLoose,3))
    out_events["jetpt4"]            = ak.Array(getN(cms_events.jetptLoose,4))
    out_events["jetpt5"]            = ak.Array(getN(cms_events.jetptLoose,5))
    out_events["jetpt6"]            = ak.Array(getN(cms_events.jetptLoose,6))
                                   
    out_events["jeteta0"]           = ak.Array(getN(cms_events.jetetaTight,0))
    out_events["jeteta1"]           = ak.Array(getN(cms_events.jetetaLoose,1))
    out_events["jeteta2"]           = ak.Array(getN(cms_events.jetetaLoose,2))
    out_events["jeteta3"]           = ak.Array(getN(cms_events.jetetaLoose,3))
    out_events["jeteta4"]           = ak.Array(getN(cms_events.jetetaLoose,4))
    out_events["jeteta5"]           = ak.Array(getN(cms_events.jetetaLoose,5))
    out_events["jeteta6"]           = ak.Array(getN(cms_events.jetetaLoose,6))
                                   
    out_events["jetphi0"]           = ak.Array(getN(cms_events.jetphiTight,0))
    out_events["jetphi1"]           = ak.Array(getN(cms_events.jetphiLoose,1))
    out_events["jetphi2"]           = ak.Array(getN(cms_events.jetphiLoose,2))
    
    out_events["jetflav0"]          = ak.Array(getN(cms_events.jetflav[cms_events.jet_sel_tight==True],0))
    out_events["jetflav1"]          = ak.Array(getN(cms_events.jetflav[cms_events.jet_sel_loose==True],1))
    out_events["jetflav2"]          = ak.Array(getN(cms_events.jetflav[cms_events.jet_sel_loose==True],2))
    out_events["jetflav3"]          = ak.Array(getN(cms_events.jetflav[cms_events.jet_sel_loose==True],3))
    out_events["jetflav4"]          = ak.Array(getN(cms_events.jetflav[cms_events.jet_sel_loose==True],4))
    out_events["jetflav5"]          = ak.Array(getN(cms_events.jetflav[cms_events.jet_sel_loose==True],5))
    out_events["jetflav6"]          = ak.Array(getN(cms_events.jetflav[cms_events.jet_sel_loose==True],6))
    
    out_events["csv0"]              = ak.Array(getN(cms_events.jetcsv[cms_events.jet_sel_tight==True],0))
    out_events["csv1"]              = ak.Array(getN(cms_events.jetcsv[cms_events.jet_sel_loose==True],1))
    out_events["csv2"]              = ak.Array(getN(cms_events.jetcsv[cms_events.jet_sel_loose==True],2))
    out_events["csv3"]              = ak.Array(getN(cms_events.jetcsv[cms_events.jet_sel_loose==True],3))
    
    out_events["SR_2b"]             = cms_events["mask_SR2b"]
    out_events["SR_1b"]             = cms_events["mask_SR1b"]
    out_events["ZeeCR_2b"]          = cms_events["mask_Zee2b"] 
    out_events["ZeeCR_1b"]          = cms_events["mask_Zee1b"] 
    out_events["ZmumuCR_2b"]        = cms_events["mask_Zmumu2b"] 
    out_events["ZmumuCR_1b"]        = cms_events["mask_Zmumu1b"] 
    out_events["TopenuCR_2b"]       = cms_events["mask_topenu2b"]
    out_events["TopenuCR_1b"]       = cms_events["mask_topenu1b"]
    out_events["TopmunuCR_2b"]      = cms_events["mask_topmunu2b"]
    out_events["TopmunuCR_1b"]      = cms_events["mask_topmunu1b"]
    out_events["WenuCR_1b"]         = cms_events["mask_wenu1b"]
    out_events["WenuCR_2b"]         = cms_events["mask_wenu2b"]
    out_events["WmunuCR_1b"]        = cms_events["mask_wmunu1b"]
    out_events["WmunuCR_2b"]        = cms_events["mask_wmunu2b"]
    
    
    
    
    
    ## btagging SFs 
    from read_sfs import btag_sf
    from read_sfs import evaluator
    
    out_events["btagsf0"]           = btag_sf.eval("central", out_events.jetflav0, abs(out_events.jeteta0), out_events.jetpt0)
    out_events["btagsf1"]           = btag_sf.eval("central", out_events.jetflav1, abs(out_events.jeteta1), out_events.jetpt1)
    out_events["btagsf2"]           = btag_sf.eval("central", out_events.jetflav2, abs(out_events.jeteta2), out_events.jetpt2)
    out_events["btagsf3"]           = btag_sf.eval("central", out_events.jetflav3, abs(out_events.jeteta3), out_events.jetpt3)
    out_events["btagsf4"]           = btag_sf.eval("central", out_events.jetflav4, abs(out_events.jeteta4), out_events.jetpt4)
    out_events["btagsf5"]           = btag_sf.eval("central", out_events.jetflav5, abs(out_events.jeteta5), out_events.jetpt5)
    out_events["btagsf6"]           = btag_sf.eval("central", out_events.jetflav6, abs(out_events.jeteta6), out_events.jetpt6)
                                    
    ## btag efficiency              
    out_events["btag_eff_lwp_0"]    = evaluator["btag_eff_lwp"](out_events.jeteta0, out_events.jetpt0)
    out_events["btag_eff_lwp_1"]    = evaluator["btag_eff_lwp"](out_events.jeteta1, out_events.jetpt1) 
                                    
    out_events["ctag_eff_lwp_0"]    = evaluator["ctag_eff_lwp"](out_events.jeteta0, out_events.jetpt0)
    out_events["ctag_eff_lwp_1"]    = evaluator["ctag_eff_lwp"](out_events.jeteta1, out_events.jetpt1) 
                                    
    out_events["ltag_eff_lwp_0"]    = evaluator["ltag_eff_lwp"](out_events.jeteta0, out_events.jetpt0)
    out_events["ltag_eff_lwp_1"]    = evaluator["ltag_eff_lwp"](out_events.jeteta1, out_events.jetpt1) 
                                    
    out_events["btag_eff_mwp_0"]    = evaluator["btag_eff_mwp"](out_events.jeteta0, out_events.jetpt0)
    out_events["btag_eff_mwp_1"]    = evaluator["btag_eff_mwp"](out_events.jeteta1, out_events.jetpt1) 
                                    
    out_events["ctag_eff_mwp_0"]    = evaluator["ctag_eff_mwp"](out_events.jeteta0, out_events.jetpt0)
    out_events["ctag_eff_mwp_1"]    = evaluator["ctag_eff_mwp"](out_events.jeteta1, out_events.jetpt1) 
                                    
    out_events["ltag_eff_mwp_0"]    = evaluator["ltag_eff_mwp"](out_events.jeteta0, out_events.jetpt0)
    out_events["ltag_eff_mwp_1"]    = evaluator["ltag_eff_mwp"](out_events.jeteta1, out_events.jetpt1) 
                                    
                                    
    ## ele sfs                      
    out_events["eleTightSF0"]       = evaluator["EGamma_SF2D_T"](out_events.eleeta0, out_events.elept0)
    out_events["eleLooseSF1"]       = evaluator["EGamma_SF2D_L"](out_events.eleeta1, out_events.elept1)
    out_events["eleTrigSF0"]        = evaluator["EGamma_SF2D_Trig"](out_events.eleeta0, out_events.elept0)
    out_events["eleRecoSF0"]        = evaluator["EGamma_SF2D_Reco"](out_events.eleeta0, out_events.elept0)
                                    
    eleRecoSF1_hi                   = evaluator["EGamma_SF2D_Reco"](out_events.eleeta1, out_events.elept1)
    eleRecoSF1_lo                   = evaluator["EGamma_SF2D_Reco_lowpt"](out_events.eleeta1, out_events.elept1)
                                    
    eleRecoSF1_hi_                  = ak.fill_none( ak.mask(  eleRecoSF1_hi , out_events.elept1>20.) ,0 )
    eleRecoSF1_lo_                  = ak.fill_none( ak.mask(  eleRecoSF1_lo , out_events.elept1>20.) ,0 )
    out_events["eleRecoSF1"]        = eleRecoSF1_hi_ + eleRecoSF1_lo_
    
    
    ## muon sfs 
    bcdef_lumi  =19.554725529
    gh_lumi = 16.224846377
    total_lumi = bcdef_lumi + gh_lumi
    
    ##--------low pt Loose
    muonLooseIDSF_lowpt1                = ( (bcdef_lumi*evaluator["muon_lowpt_BCDEF_LooseID"](out_events.mupt1, abs(out_events.mueta1) ))  +  
                                            (gh_lumi*evaluator["muon_lowpt_GH_LooseID"](out_events.mupt1, abs(out_events.mueta1) ))        )/total_lumi
                                        
                                        
    ##----------- medium pt Loose       
    muonLooseIDSF1                      = ( (bcdef_lumi*evaluator["muon_highpt_BCDEF_LooseID"](out_events.mueta1, out_events.mupt1))  +  
                                            (gh_lumi*evaluator["muon_highpt_GH_LooseID"](out_events.mueta1, out_events.mupt1))              )/total_lumi
                                        
                                        
    muonLooseISOSF1                     = ( (bcdef_lumi*evaluator["muon_highpt_BCDEF_LooseISO"](out_events.mueta1, out_events.mupt1))  +  
                                            (gh_lumi*evaluator["muon_highpt_GH_LooseISO"](out_events.mueta1, out_events.mupt1))             )/total_lumi
                                        
    muon_loose_ID_low_SF_1              = ak.fill_none( ak.mask(  muonLooseIDSF_lowpt1 , out_events.mupt1<20.) ,0 )
    muon_loose_ID_high_SF_1             = ak.fill_none( ak.mask(  muonLooseIDSF1       , out_events.mupt1>20.) ,0 )
                                        
    muon_loose_ID_SF_1                  = muon_loose_ID_low_SF_1 + muon_loose_ID_high_SF_1
                                        
    out_events["muLooseSF1"]            = muon_loose_ID_SF_1 * muonLooseISOSF1 
                                        
                                        
    ##------------medium pt tight       
    muonTightIDSF0                      = ( (bcdef_lumi*evaluator["muon_highpt_BCDEF_TightID"](out_events.mueta0, out_events.mupt0))  +  
                                            (gh_lumi*evaluator["muon_highpt_GH_TightID"](out_events.mueta0, out_events.mupt0))               )/total_lumi
                                        
    muonTightISOSF0                     = ( (bcdef_lumi*evaluator["muon_highpt_BCDEF_TightISO"](out_events.mueta0, out_events.mupt0))  +  
                                            (gh_lumi*evaluator["muon_highpt_GH_TightISO"](out_events.mueta0, out_events.mupt0))               )/total_lumi
                                        
    out_events["muTightSF0"]            = muonTightIDSF0 * muonTightISOSF0
                                        
    out_events["puweight"]              = evaluator["pu_weight"](cms_events.nTrueInt)
    
    
    ## trigger sfs 
    out_events["mettrigWeight"]         = evaluator["met_trig"](cms_events.metpt)
    out_events["recoilWmunutrigWeight"] = evaluator["met_trig"](cms_events.recoil_Wmunu0)
    out_events["recoilWenutrigWeight"]  = evaluator["met_trig"](cms_events.recoil_Wenu0)
    out_events["recoilZmumutrigWeight"] = evaluator["met_trig"](cms_events.Zmumu_recoil)
    out_events["recoilZeetrigWeight"]   = evaluator["met_trig"](cms_events.Zee_recoil)
    
    
    
    ## Fill weights for each CR so that we don't need to worry later 
    
    out_events["weight_SR_2b"]          = out_events.puweight * out_events.mettrigWeight 
    out_events["weight_SR_1b"]          = out_events.puweight * out_events.mettrigWeight 
    
    out_events["weight_ZeeCR_2b"]       = out_events.puweight * out_events.eleTrigSF0 
    out_events["weight_ZeeCR_1b"]       = out_events.puweight * out_events.eleTrigSF0 
    
    out_events["weight_ZmumuCR_2b"]     = out_events.puweight * out_events.recoilZmumutrigWeight 
    out_events["weight_ZmumuCR_1b"]     = out_events.puweight * out_events.recoilZmumutrigWeight 
    
    out_events["weight_TopenuCR_2b"]    = out_events.puweight * out_events.eleTrigSF0 
    out_events["weight_TopenuCR_1b"]    = out_events.puweight * out_events.eleTrigSF0 
    
    out_events["weight_TopmunuCR_2b"]   = out_events.puweight * out_events.recoilWmunutrigWeight 
    out_events["weight_TopmunuCR_1b"]   = out_events.puweight * out_events.recoilWmunutrigWeight 
    
    out_events["weight_WenuCR_1b"]      = out_events.puweight * out_events.eleTrigSF0 
    out_events["weight_WenuCR_2b"]      = out_events.puweight * out_events.eleTrigSF0 
    
    out_events["weight_WmunuCR_1b"]     = out_events.puweight * out_events.recoilWmunutrigWeight
    out_events["weight_WmunuCR_2b"]     = out_events.puweight * out_events.recoilWmunutrigWeight
    
    
    ## Fill Histograms 
    from variables import vardict, regions, variables_common
    from binning import binning
    
    f = TFile(outputfile,"RECREATE")
    for ireg in regions:
        thisregion  = out_events[out_events[ireg]==True]
        thisregion_ = thisregion[~(ak.is_none(thisregion))]
        weight_ = "weight_"+ireg
        
        for ivar in variables_common[ireg]:
            hist_name_ = "h_reg_"+ireg+"_"+vardict[ivar]
            h = VarToHist(thisregion_[ivar], thisregion_[weight_], hist_name_, binning[ireg][ivar])
            f.cd()
            h.Write()
    
    h_total  = TH1F("h_total_mcweight", "h_total_mcweight", 2,0,2)
    h_total.SetBinContent(1,nevents)
    f.cd()
    h_total.Write()




    
    write_parquet=False
    if write_parquet:
        ak.to_parquet(out_events, "analysis_wjets_allevents.parquet")
    
    
    

def main():
    runOneFile("")

if __name__ == "__main__":
    main()
    end = time.clock()
    print("%.4gs" % (end-start))
    end = time.clock()
    print("%.4gs" % (end-start))
