#include "LatinoAnalysis/MultiDraw/interface/TTreeFunction.h"
#include "LatinoAnalysis/MultiDraw/interface/FunctionLibrary.h"
#include "TSystem.h"
#include "iostream"
#include "vector"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TSystem.h"
#include <map>
#include "TString.h"
#include "momemta/ConfigurationReader.h"
#include "momemta/MoMEMta.h"
#include "momemta/Types.h"

using LorentzVectorM = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>;

using namespace momemta;

void normalizeInput(LorentzVector& p4) {
  if (p4.M() > 0) return;

  // Increase the energy until M is positive                                                                                                                                                              
  p4.SetE(p4.P());
  while (p4.M2() < 0) {
    double delta = p4.E() * 1e-5;
    p4.SetE(p4.E() + delta);
  };
}


float recoMoMEMta_WW(lep1pt, lep2pt, lep1phi, lep2phi, lep1eta, lep2eta, lep1Id, lep2Id, met, jet1pt, jet2pt, jet1phi, jet2phi, jet1eta, jet2eta, path){
  
  logging::set_level(logging::level::off);
  
  //Initializing 4-vectors
  TLorentzVector L1(0.,0.,0.,0.);
  TLorentzVector L2(0.,0.,0.,0.);
  TLorentzVector J1(0.,0.,0.,0.);
  TLorentzVector J2(0.,0.,0.,0.);
  
  
  L1.SetPtEtaPhiM(lep1pt, lep1eta, lep1phi, 0.0);
  L2.SetPtEtaPhiM(lep2pt, lep2eta, lep2phi, 0.0);
  J1.SetPtEtaPhiM(jet1pt, jet1eta, jet1phi, 0.0);
  J2.SetPtEtaPhiM(jet2pt, jet2eta, jet2phi, 0.0);
  
  
  ConfigurationReader configuration_top(path + "/TTbar_FullyLeptonic/TTbar_FullyLeptonic.lua");
  ConfigurationReader configuration_WW(path + "/TTbar_FullyLeptonic/WW_leptonic.lua");

  if (lep1Id < 0){
    configuration_top.getGlobalParameters().set("top_mass", 173.);
  }else{
    ConfigurationReader configuration_top(path + "/TTbar_FullyLeptonic/TTbar_FullyLeptonic_mue.lua");
    configuration_top.getGlobalParameters().set("top_mass", 173.);
    
    ConfigurationReader configuration_WW(path + "/TTbar_FullyLeptonic/WW_leptonic_mue.lua");
  }

  logging::set_level(logging::level::off);
  
  MoMEMta weight_top(configuration_top.freeze());
  MoMEMta weight_WW(configuration_WW.freeze());

  ParameterSet lua_parameters;
  lua_parameters.set("USE_TF", true);
  lua_parameters.set("USE_PERM", true);
  
  momemta::Particle lepton1 { "lepton1", LorentzVector(L1.Px(), L1.Py(), L1.Pz(), L1.E()), lep1 }; // muon                                                                                                                                
  momemta::Particle lepton2 { "lepton2", LorentzVector(L2.Px(), L2.Py(), L2.Pz(), L2.E()), lep2 }; // electron   
  momemta::Particle bjet1 { "bjet1", LorentzVector(J1.Px(), J1.Py(), J1.Pz(), J1.E()), 5 }; // Not necessary a bjet, but passed to MoMEMta as if it is                                                                                    
  momemta::Particle bjet2 { "bjet2", LorentzVector(J2.Px(), J2.Py(), J2.Pz(), J2.E()), -5 };
  
  LorentzVector met_p4 {NuNu.Px(), NuNu.Py(), NuNu.Pz(), NuNu.E()};
  
  // Normalize Input four-vectors for numerical estability
  normalizeInput(lepton1.p4);
  normalizeInput(lepton2.p4);
  normalizeInput(bjet1.p4);
  normalizeInput(bjet2.p4);
  
  std::vector<std::pair<double, double>> weights_top = weight_top.computeWeights({lepton1, bjet1, lepton2, bjet2}, met_p4);
  std::vector<std::pair<double, double>> weights_WW = weight_WW.computeWeights({lepton1, lepton2}, met_p4);
  
  double ME_top = (double)weights_top.back().first;
  double ME_WW =  (double)weights_WW.back().first;
  
  D_Top = abs(ME_WW) / (abs(ME_WW) + 1e4 * abs(ME_top));
  
  return D_Top
  
}
