#!/usr/bin/env python

import os
import glob
import numpy as np
import pandas as pd
import ROOT
import sys


###################################################
##                                               ##
### COMPUTE POLARIZED FRACTIONS FOR WW ANALYSIS ### 
##                                               ##
###################################################

def compute_weights(GenPart_eta, GenPart_pt, GenPart_mass, GenPart_phi, pos_mother, status, particles):
  
  # Define the four-vectors
  
  Wp = ROOT.Math.PtEtaPhiEVector()
  Wm = ROOT.Math.PtEtaPhiEVector()
  genlp   = ROOT.Math.PtEtaPhiEVector()
  genlm = ROOT.Math.PtEtaPhiEVector()
  gennup = ROOT.Math.PtEtaPhiEVector()
  gennum = ROOT.Math.PtEtaPhiEVector()
  vector_lp = ROOT.TLorentzVector()
  vector_lm = ROOT.TLorentzVector()
  vector_nup = ROOT.TLorentzVector()
  vector_num = ROOT.TLorentzVector()
  
  genWp = ROOT.Math.PtEtaPhiEVector()
  genWm = ROOT.Math.PtEtaPhiEVector()
  vector_Wp = ROOT.TLorentzVector()
  vector_Wm = ROOT.TLorentzVector()
  
  number_elec = 0
  number_muon = 0
  number_tau = 0
  
  pos_wp = 999
  pos_wm = 999
  
  GenPart_status = np.array(status)
  GenPart_pdgId = np.array(particles)
  
  
  # Loop over generated particles
  # Save prompt electrons, muons, neutrinos and W bosons
  
  for p in range(len(particles)):
    
    if (particles[p]==11 and particles[pos_mother[p]]==-24 and particles[pos_mother[pos_mother[p]]]!=15):
        pos_wm = pos_mother[p]
        number_elec = number_elec + 1
        vector_lm.SetPtEtaPhiM(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], 0.0)
        genlm.SetCoordinates(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], vector_lm.E())
    
    elif (particles[p]==-11 and particles[pos_mother[p]]==24 and particles[pos_mother[pos_mother[p]]]!=-15):
        pos_wp = pos_mother[p]
        number_elec = number_elec + 1
        vector_lp.SetPtEtaPhiM(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], 0.0)
        genlp.SetCoordinates(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], vector_lp.E())
    
    elif (particles[p]==13 and particles[pos_mother[p]]==-24 and particles[pos_mother[pos_mother[p]]]!=15):
        pos_wm = pos_mother[p]
        number_muon = number_muon + 1
        vector_lm.SetPtEtaPhiM(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], 0.0)
        genlm.SetCoordinates(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], vector_lm.E())
    
    elif (particles[p]==-13 and particles[pos_mother[p]]==24 and particles[pos_mother[pos_mother[p]]]!=-15):
        pos_wp = pos_mother[p]
        number_muon = number_muon + 1
        vector_lp.SetPtEtaPhiM(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], 0.0)
        genlp.SetCoordinates(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], vector_lp.E())
        
    elif (particles[p]==15 and particles[pos_mother[p]]==-24):
        pos_wm = pos_mother[p]
        number_tau = number_tau + 1
        vector_lm.SetPtEtaPhiM(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], 0.0)
        genlm.SetCoordinates(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], vector_lm.E())
    
    elif (particles[p]==-15 and particles[pos_mother[p]]==24):
        pos_wp = pos_mother[p]
        number_tau = number_tau + 1
        vector_lp.SetPtEtaPhiM(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], 0.0)
        genlp.SetCoordinates(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], vector_lp.E())
         
    
    
    if (particles[p]==-12 and particles[pos_mother[p]]==-24 and particles[pos_mother[pos_mother[p]]]!=15):
        vector_num.SetPtEtaPhiM(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], 0.0)
        gennum.SetCoordinates(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], vector_num.E())
    
    elif (particles[p]==12 and particles[pos_mother[p]]==24 and particles[pos_mother[pos_mother[p]]]!=-15):
        vector_nup.SetPtEtaPhiM(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], 0.0)
        gennup.SetCoordinates(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], vector_nup.E())
    
    elif (particles[p]==-14 and particles[pos_mother[p]]==-24 and particles[pos_mother[pos_mother[p]]]!=15):
        vector_num.SetPtEtaPhiM(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], 0.0)
        gennum.SetCoordinates(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], vector_num.E())
    
    elif (particles[p]==14 and particles[pos_mother[p]]==24 and particles[pos_mother[pos_mother[p]]]!=-15):
        vector_nup.SetPtEtaPhiM(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], 0.0)
        gennup.SetCoordinates(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], vector_nup.E())
        
    elif (particles[p]==-16 and particles[pos_mother[p]]==-24):
        vector_num.SetPtEtaPhiM(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], 0.0)
        gennum.SetCoordinates(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], vector_num.E())
    
    elif (particles[p]==16 and particles[pos_mother[p]]==24):
        vector_nup.SetPtEtaPhiM(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], 0.0)
        gennup.SetCoordinates(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], vector_nup.E())

  vector_Wp.SetPtEtaPhiM(GenPart_pt.values[pos_wp], GenPart_eta.values[pos_wp], GenPart_phi.values[pos_wp], GenPart_mass.values[pos_wp]) # W plus
  genWp.SetCoordinates(GenPart_pt.values[pos_wp], GenPart_eta.values[pos_wp], GenPart_phi.values[pos_wp], vector_Wp.E())
      
  
  vector_Wm.SetPtEtaPhiM(GenPart_pt.values[pos_wm], GenPart_eta.values[pos_wm], GenPart_phi.values[pos_wm], GenPart_mass.values[pos_wm]) # W minus
  genWm.SetCoordinates(GenPart_pt.values[pos_wm], GenPart_eta.values[pos_wm], GenPart_phi.values[pos_wm], vector_Wm.E())
         
    
  if (((number_elec==0 or number_muon==0) and number_tau==0) or pos_wp==999 or pos_wm==999):
    return [-999, -999, -999, -999]
   
    
  # Boost over lepton from the Ws reference frame
  # Compute theta star
  
  wmRF = ROOT.Math.XYZVector()
  wpRF = ROOT.Math.XYZVector()
  
  wmRF = genWm.BoostToCM()
  wpRF = genWp.BoostToCM()
  
  leppWRF = ROOT.Math.XYZVector()
  lepmWRF = ROOT.Math.XYZVector()
  leppWRF = ROOT.Math.VectorUtil.boost(genlp, wpRF)
  lepmWRF = ROOT.Math.VectorUtil.boost(genlm, wmRF)
  
  theta_Wp_star = ROOT.Math.VectorUtil.Angle(leppWRF, genWp)
  theta_Wm_star = ROOT.Math.VectorUtil.Angle(lepmWRF, genWm)

  cos_Wp_theta_star = np.cos(theta_Wp_star)
  cos_Wm_theta_star = np.cos(theta_Wm_star)
  
  ###################################
  # Theoretical polarized fractions #
  ###################################
  
  # https://arxiv.org/pdf/2006.14867.pdf
  # https://arxiv.org/pdf/1204.6427.pdf
  
  f0_m = 0.26
  fL_m = 0.48
  fR_m = 0.25
  fT_m = fL_m + fR_m
  
  f0_p = 0.271
  fT_p = 0.729
  fL_p = fT_p/1.52
  fR_p = fT_p - fL_p
  
  # Compute single polarizations weights as a function of cos(theta_star)
  # Then, calculate doubly-polarized fractions
  
  # W minus
  
  weight_f0_Wm = (3/4) * f0_m * (1 - cos_Wm_theta_star*cos_Wm_theta_star)
  weight_fL_Wm = (3/8) * fL_m * (1 + cos_Wm_theta_star)*(1 + cos_Wm_theta_star)
  weight_fR_Wm = (3/8) * fR_m * (1 - cos_Wm_theta_star)*(1 - cos_Wm_theta_star)
  weight_fT_Wm = (3/8) * fL_m * (1 + cos_Wm_theta_star)*(1 + cos_Wm_theta_star) + (3/8) * fR_m * (1 - cos_Wm_theta_star)*(1 - cos_Wm_theta_star)
  weight_total_Wm = weight_f0_Wm + weight_fL_Wm + weight_fR_Wm

  # W plus
  
  weight_f0_Wp = (3/4) * f0_p * (1 - cos_Wp_theta_star*cos_Wp_theta_star)
  weight_fL_Wp = (3/8) * fL_p * (1 - cos_Wp_theta_star)*(1 - cos_Wp_theta_star)
  weight_fR_Wp = (3/8) * fR_p * (1 + cos_Wp_theta_star)*(1 + cos_Wp_theta_star)
  weight_fT_Wp = (3/8) * fL_p * (1 - cos_Wp_theta_star)*(1 - cos_Wp_theta_star) + (3/8) * fR_p * (1 + cos_Wp_theta_star)*(1 + cos_Wp_theta_star)
  weight_total_Wp = weight_f0_Wp + weight_fL_Wp + weight_fR_Wp

  # Doubly Polarized
  
  weight_LL = (weight_f0_Wp/weight_total_Wp)*(weight_f0_Wm/weight_total_Wm)
  weight_TL = (weight_fT_Wp/weight_total_Wp)*(weight_f0_Wm/weight_total_Wm)
  weight_LT = (weight_f0_Wp/weight_total_Wp)*(weight_fT_Wm/weight_total_Wm)
  weight_TT = (weight_fT_Wp/weight_total_Wp)*(weight_fT_Wm/weight_total_Wm)
  
  return [weight_LL, weight_TL, weight_LT, weight_TT]


if __name__ == "__main__":

  df = pd.read_pickle("dataset_ww.pkl")
  
  LL = []
  TL = []
  LT = []
  TT = []
  
  for event in range(0, len(df)):
    
    tmp = df.loc[event]

    # Get GenParticles as input 
  
    GenPart_eta = tmp.GenPart_eta
    GenPart_pt = tmp.GenPart_pt
    GenPart_mass = tmp.GenPart_mass
    GenPart_phi = tmp.GenPart_phi
    GenPartgenPartIdxMother = tmp.GenPartgenPartIdxMother
    GenPart_status = tmp.GenPart_status
    GenPart_pdgId = tmp.GenPart_pdgId
    
    weights_result = compute_weights(GenPart_eta, GenPart_pt, GenPart_mass, GenPart_phi, GenPartgenPartIdxMother, GenPart_status, GenPart_pdgId)
    
    LL.append(weights_result[0])
    TL.append(weights_result[1])
    LT.append(weights_result[2])
    TT.append(weights_result[3])
    
  df["weight_LL"] = LL
  df["weight_TL"] = TL
  df["weight_LT"] = LT
  df["weight_TT"] = TT

  df.to_pickle("dataset_ww_weighted.pkl")
  
