import sys 
import uproot 
import numpy 
import pandas as pd 
import math 
import time

isCondor=False 

if isCondor:
    sys.path.append('ExoPieUtils/analysisutils/')
else:
    sys.path.append('../../ExoPieUtils/analysisutils/')


if isCondor:
    sys.path.append('ExoPieUtils/scalefactortools/')
else:
    sys.path.append('../../ExoPieUtils/scalefactortools/')


year='2016'
year_file = open("Year.py", "w")
if year == '2016':
    print('code is running for 2016')
    year_file.write('era="2016"')
elif year == '2017':
    print('code is running for 2017')
    year_file.write('era="2017"')
elif year == '2018':
    print('code is running for 2018')
    year_file.write('era="2018"')
else:
    print('please provide year')
    sys.exit()
year_file.close()

import ana_weight as wgt



def weight_( ep_WmunuRecoil, nEle, nMu, ep_muPt, ep_muEta):#, ep_ZmumuRecoil, ep_WmunuRecoil, nEle, ep_elePt, ep_eleEta, ep_eleIsPTight, nMu, ep_muPt, ep_muEta, ep_isTightMuon):
    total_weight=1.0
    if (nEle==0) & (nMu==1) :
        total_weight = weight_W1mu_(ep_WmunuRecoil, nEle, nMu, ep_muPt, ep_muEta)
        
    
    return total_weight

def weight_W1mu_(ep_WmunuRecoil, nEle, nMu, ep_muPt, ep_muEta):
    weightMET = 1.0
    weightMu = 1.0
    weightMET, weightMET_up, weightMET_down     = wgt.getMETtrig_First( ep_WmunuRecoil, 'R')
    weightMu, weightMu_up, weightMu_down        = wgt.mu_weight(ep_muPt, ep_muEta, 'T')
    return ( weightMET * weightMu)
    
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

        MT = numpy.sqrt( 2 * elept * met_ * (1.0 - numpy.cos(dphi)) )
        WenuRecoilPx = -( met_*numpy.cos(metphi_) + elepx_)
        WenuRecoilPy = -( met_*numpy.sin(metphi_) + elepy_)
        WenuRecoilPt = numpy.sqrt(WenuRecoilPx**2  +  WenuRecoilPy**2)
        WenurecoilPhi = numpy.arctan2(WenuRecoilPx,WenuRecoilPy)
    return WenuRecoilPt, WenurecoilPhi, MT
    
def Delta_R(eta1, eta2, phi1,phi2):
    deltaeta = eta1-eta2
    deltaphi = DeltaPhi(phi1,phi2)
    DR = numpy.sqrt ( deltaeta**2 + deltaphi**2 )
    return DR 

def jetcleaning(ak4eta, lepeta, ak4phi, lepphi, DRCut):
    ## usage: (obj_to_clean, obj_cleaned_against, so on
    dr_ = Delta_R(ak4eta, lepeta, ak4phi, lepphi)
    
    return (dr_ > DRCut)



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

