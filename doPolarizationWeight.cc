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

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/Rotation3D.h"
#include "Math/EulerAngles.h"
#include "Math/AxisAngle.h"
#include "Math/Quaternion.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/RotationZ.h"
#include "Math/RotationZYX.h"
#include "Math/LorentzRotation.h"
#include "Math/Boost.h"
#include "Math/BoostX.h"
#include "Math/BoostY.h"
#include "Math/BoostZ.h"
#include "Math/Transform3D.h"
#include "Math/Plane3D.h"
#include "Math/VectorUtil.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMath.h"


/**
###################################################
##                                               ##
### COMPUTE POLARIZED FRACTIONS FOR WW ANALYSIS ### 
##                                               ##
###################################################
**/

class doPolarizationWeight : public multidraw::TTreeFunction {
public:
  //Class Constructor 
  doPolarizationWeight(char const* name);
  //Class Destructor 
  ~doPolarizationWeight() {
  }
  //Functions from Multidraw namespace (TTreeFunction class)
  char const* getName() const override {return "doPolarizationWeight"; }
  TTreeFunction* clone() const override {return new doPolarizationWeight(name_.c_str());}
  unsigned getNdata() override {return 1; }
  //This function will return the required value
  double evaluate(unsigned) override;

protected:
  void bindTree_(multidraw::FunctionLibrary&) override;
  std::string name_;
	
  FloatArrayReader* GenPart_pt{};
  FloatArrayReader* GenPart_eta{};
  FloatArrayReader* GenPart_phi{};
  FloatArrayReader* GenPart_mass{};
  FloatArrayReader* GenPart_pdgId{};
  FloatArrayReader* GenPart_status{};
  FloatArrayReader* GenPart_genPartIdxMother{};

private:

  Double_t LHCsqrts_= 13., mh_= 125.;
  
};

doPolarizationWeight::doPolarizationWeight(char const* name):
  TTreeFunction()
{
  name_ = name;
}

double
doPolarizationWeight::evaluate(unsigned)
{
   
  ROOT::Math::PtEtaPhiEVector genWp;
  ROOT::Math::PtEtaPhiEVector genWm;
  ROOT::Math::PtEtaPhiEVector genlp;
  ROOT::Math::PtEtaPhiEVector genlm;
  ROOT::Math::PtEtaPhiEVector gennup;
  ROOT::Math::PtEtaPhiEVector gennum;
  TLorentzVector vector_lp; 
  TLorentzVector vector_lm; 
  TLorentzVector vector_nup;
  TLorentzVector vector_num;
  TLorentzVector vector_Wp; 
  TLorentzVector vector_Wm;
    
  Int_t number_elec = 0;
  Int_t number_muon = 0;
  
  Int_t pos_wp = 999;
  Int_t pos_wm = 999;
  
  Int_t mother_pos = 0;
  
  Double_t nGen = GenPart_pt->GetSize();
  
  for (unsigned int p = 0; p < nGen; p++){
  
    mother_pos = GenPart_genPartIdxMother->At(p);
    if (GenPart_pdgId->At(p)==11 && GenPart_status->At(p)==1 && GenPart_pdgId->At(mother_pos)==-24){
      pos_wm = mother_pos;
      number_elec++;
      vector_lm.SetPtEtaPhiM(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), 0.0);
      genlm.SetCoordinates(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), vector_lm.E());
    }else if (GenPart_pdgId->At(p)==-11 && GenPart_status->At(p)==1 && GenPart_pdgId->At(mother_pos)==24){
      pos_wp = mother_pos;
      number_elec++;
      vector_lp.SetPtEtaPhiM(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), 0.0);
      genlp.SetCoordinates(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), vector_lp.E());    
    }else if (GenPart_pdgId->At(p)==13 && GenPart_status->At(p)==1 && GenPart_pdgId->At(mother_pos)==-24){
      pos_wm = mother_pos;
      number_muon++;
      vector_lm.SetPtEtaPhiM(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), 0.0);
      genlm.SetCoordinates(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), vector_lm.E());
    }if (GenPart_pdgId->At(p)==-13 && GenPart_status->At(p)==1 && GenPart_pdgId->At(mother_pos)==24){
      pos_wp = mother_pos;
      number_muon++;
      vector_lp.SetPtEtaPhiM(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), 0.0);
      genlp.SetCoordinates(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), vector_lp.E());    
    }
    
    if (GenPart_pdgId->At(p)==-12 && GenPart_status->At(p)==1 && GenPart_pdgId->At(mother_pos)==-24){
      vector_num.SetPtEtaPhiM(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), 0.0);
      gennum.SetCoordinates(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), vector_num.E());
    }else if (GenPart_pdgId->At(p)==12 && GenPart_status->At(p)==1 && GenPart_pdgId->At(mother_pos)==24){
      vector_nup.SetPtEtaPhiM(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), 0.0);
      gennup.SetCoordinates(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), vector_nup.E());    
    }else if (GenPart_pdgId->At(p)==-14 && GenPart_status->At(p)==1 && GenPart_pdgId->At(mother_pos)==-24){
      vector_num.SetPtEtaPhiM(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), 0.0);
      gennum.SetCoordinates(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), vector_num.E());
    }if (GenPart_pdgId->At(p)==14 && GenPart_status->At(p)==1 && GenPart_pdgId->At(mother_pos)==24){
      vector_nup.SetPtEtaPhiM(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), 0.0);
      gennup.SetCoordinates(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), vector_nup.E());    
    }
  
  } // End loop over particles
	
  // Assign the four-vector of the Ws
  vector_Wp.SetPtEtaPhiM(GenPart_pt->At(pos_wp), GenPart_eta->At(pos_wp), GenPart_phi->At(pos_wp), GenPart_mass->At(pos_wp)); // W plus
  genWp.SetCoordinates(GenPart_pt->At(pos_wp), GenPart_eta->At(pos_wp), GenPart_phi->At(pos_wp), vector_Wp.E());
      
  vector_Wm.SetPtEtaPhiM(GenPart_pt->At(pos_wm), GenPart_eta->At(pos_wm), GenPart_phi->At(pos_wm), GenPart_mass->At(pos_wm)); // W minus
  genWm.SetCoordinates(GenPart_pt->At(pos_wm), GenPart_eta->At(pos_wm), GenPart_phi->At(pos_wm), vector_Wm.E());
  
  std::vector<Double_t> weights;	
	
  if (number_elec!=1 || number_muon!=1 || pos_wp==999 || pos_wm==999){
	  weights = {-999, -999, -999, -999};
	  return  weights;
  }

  // Boost over lepton from the Ws reference frame
  // Compute theta star  

  ROOT::Math::XYZVector wpRF;
  ROOT::Math::XYZVector wmRF;
	  
  wmRF = genWm.BoostToCM();
  wpRF = genWp.BoostToCM();
  
  ROOT::Math::XYZVector leppWRF;
  ROOT::Math::XYZVector lepmWRF;	
	
  leppWRF = ROOT::Math::VectorUtil::boost(genlp, wpRF);
  lepmWRF = ROOT::Math::VectorUtil::boost(genlm, wmRF);
  
  Double_t theta_Wp_star = ROOT::Math::VectorUtil::Angle(leppWRF, genWp);
  Double_t theta_Wm_star = ROOT::Math::VectorUtil::Angle(lepmWRF, genWm);

  Double_t cos_Wp_theta_star = ROOT::Math::cos(theta_Wp_star);
  Double_t cos_Wm_theta_star = ROOT::Math::cos(theta_Wm_star);	

  	
  /////////////////////////////////////
  // Theoretical polarized fractions //
  /////////////////////////////////////
  
  // https://arxiv.org/pdf/2006.14867.pdf
  // https://arxiv.org/pdf/1204.6427.pdf
  
  Double_t f0_m = 0.26;
  Double_t fL_m = 0.48;
  Double_t fR_m = 0.25;
  Double_t fT_m = fL_m + fR_m;
  
  Double_t f0_p = 0.271;
  Double_t fT_p = 0.729;
  Double_t fL_p = fT_p/1.52;
  Double_t fR_p = fT_p - fL_p;	
	
  // W minus
  
  Double_t weight_f0_Wm = (3.0/4.0) * f0_m * (1 - cos_Wm_theta_star*cos_Wm_theta_star);
  Double_t weight_fL_Wm = (3.0/8.0) * fL_m * (1 + cos_Wm_theta_star)*(1 + cos_Wm_theta_star);
  Double_t weight_fR_Wm = (3.0/8.0) * fR_m * (1 - cos_Wm_theta_star)*(1 - cos_Wm_theta_star);
  Double_t weight_fT_Wm = (3.0/8.0) * fL_m * (1 + cos_Wm_theta_star)*(1 + cos_Wm_theta_star) + (3.0/8.0) * fR_m * (1 - cos_Wm_theta_star)*(1 - cos_Wm_theta_star);
  Double_t weight_total_Wm = weight_f0_Wm + weight_fL_Wm + weight_fR_Wm;

  // W plus
  
  Double_t weight_f0_Wp = (3.0/4.0) * f0_p * (1 - cos_Wp_theta_star*cos_Wp_theta_star);
  Double_t weight_fL_Wp = (3.0/8.0) * fL_p * (1 - cos_Wp_theta_star)*(1 - cos_Wp_theta_star);
  Double_t weight_fR_Wp = (3.0/8.0) * fR_p * (1 + cos_Wp_theta_star)*(1 + cos_Wp_theta_star);
  Double_t weight_fT_Wp = (3.0/8.0) * fL_p * (1 - cos_Wp_theta_star)*(1 - cos_Wp_theta_star) + (3.0/8.0) * fR_p * (1 + cos_Wp_theta_star)*(1 + cos_Wp_theta_star);
  Double_t weight_total_Wp = weight_f0_Wp + weight_fL_Wp + weight_fR_Wp;

  // Doubly Polarized
  
  Double_t weight_LL = (weight_f0_Wp/weight_total_Wp)*(weight_f0_Wm/weight_total_Wm);
  Double_t weight_TL = (weight_fT_Wp/weight_total_Wp)*(weight_f0_Wm/weight_total_Wm);
  Double_t weight_LT = (weight_f0_Wp/weight_total_Wp)*(weight_fT_Wm/weight_total_Wm);
  Double_t weight_TT = (weight_fT_Wp/weight_total_Wp)*(weight_fT_Wm/weight_total_Wm);
  
  weights = {weight_LL, weight_TL, weight_LT, weight_TT, cos_Wp_theta_star, cos_Wm_theta_star};
	
  if (name_=="LL"){
  	return weights[0];
  }else if (name_=="TL"){
  	return weights[1];
  }else if (name_=="LT"){
  	return weights[2];
  }else if (name_=="TT"){
  	return weights[3];
  }else if (name_=="cos_wp"){
  	return weights[4];
  }else if (name_=="cos_wm"){
  	return weights[5];
  }else{
  	return -999;
  }
	
  return weights;
	  
}	  
	  
void
doPolarizationWeight::bindTree_(multidraw::FunctionLibrary& _library)
{
  //GenPart
  _library.bindBranch(GenPart_eta, "GenPart_eta");
  _library.bindBranch(GenPart_pt, "GenPart_pt");
  _library.bindBranch(GenPart_mass, "GenPart_mass");
  _library.bindBranch(GenPart_phi, "GenPart_phi");
  _library.bindBranch(GenPart_status, "GenPart_status");
  _library.bindBranch(GenPart_pdgId, "GenPart_pdgId");
  _library.bindBranch(GenPart_genPartIdxMother, "GenPart_genPartIdxMother");
}
