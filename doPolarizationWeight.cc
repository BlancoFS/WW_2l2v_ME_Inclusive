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
  doPolarizationWeight();
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

doPolarizationWeight::doPolarizationWeight():
  TTreeFunction()
{

}

double
doPolarizationWeight::evaluate(unsigned)
{
   
  ROOT::Math::PtEtaPhiEVector Wp;
  ROOT::Math::PtEtaPhiEVector Wm;
  ROOT::Math::PtEtaPhiEVector genlp;
  ROOT::Math::PtEtaPhiEVector genlm;
  ROOT::Math::PtEtaPhiEVector gennup;
  ROOT::Math::PtEtaPhiEVector gennum;
  TLorentzVector vector_lp; 
  TLorentzVector vector_lm; 
  TLorentzVector vector_nup;
  TLorentzVector vector_num;
    
  Int_t number_elec = 0;
  Int_t number_muon = 0;
  
  Int_t pos_wp = 999;
  Int_t pos_wm = 999;
  
  Int_t mother_pos = 0;
  
  Double nGen{*GenPart_pt->size()};
  
  for (unsigned int p = 0; p < nGen; p++){
  
    mother_pos = GenPart_genPartIdxMother->At(p);
    if (GenPart_pdgId->At(p)==11 && GenPart_status->At(p)==1 and GenPart_pdgId->At(mother_pos)==-24){
      pos_wm = mother_pos;
      number_elec++;
      vector_lm.SetPtEtaPhiM(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), 0.0);
      genlm.SetCoordinates(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), vector_lm.E());
    }else if (GenPart_pdgId->At(p)==-11 && GenPart_status->At(p)==1 and GenPart_pdgId->At(mother_pos)==24){
      pos_wp = mother_pos;
      number_elec++;
      vector_lp.SetPtEtaPhiM(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), 0.0);
      genlp.SetCoordinates(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), vector_lp.E());    
    }else if (GenPart_pdgId->At(p)==13 && GenPart_status->At(p)==1 and GenPart_pdgId->At(mother_pos)==-24){
      pos_wm = mother_pos;
      number_muon++;
      vector_lm.SetPtEtaPhiM(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), 0.0);
      genlm.SetCoordinates(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), vector_lm.E());
    }if (GenPart_pdgId->At(p)==-13 && GenPart_status->At(p)==1 and GenPart_pdgId->At(mother_pos)==24){
      pos_wp = mother_pos;
      number_muon++;
      vector_lp.SetPtEtaPhiM(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), 0.0);
      genlp.SetCoordinates(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), vector_lp.E());    
    }
    
    if (GenPart_pdgId->At(p)==-12 && GenPart_status->At(p)==1 and GenPart_pdgId->At(mother_pos)==-24){
      vector_num.SetPtEtaPhiM(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), 0.0);
      gennum.SetCoordinates(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), vector_num.E());
    }else if (GenPart_pdgId->At(p)==12 && GenPart_status->At(p)==1 and GenPart_pdgId->At(mother_pos)==24){
      vector_nup.SetPtEtaPhiM(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), 0.0);
      gennup.SetCoordinates(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), vector_nup.E());    
    }else if (GenPart_pdgId->At(p)==-14 && GenPart_status->At(p)==1 and GenPart_pdgId->At(mother_pos)==-24){
      vector_num.SetPtEtaPhiM(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), 0.0);
      gennum.SetCoordinates(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), vector_num.E());
    }if (GenPart_pdgId->At(p)==14 && GenPart_status->At(p)==1 and GenPart_pdgId->At(mother_pos)==24){
      vector_nup.SetPtEtaPhiM(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), 0.0);
      gennup.SetCoordinates(GenPart_pt->At(p), GenPart_eta->At(p), GenPart_phi->At(p), vector_nup.E());    
    }
  
  } // End loop over particles
  
  
  //Initializing 4-vectors
  TLorentzVector L1(0.,0.,0.,0.);
  TLorentzVector L2(0.,0.,0.,0.);
  TLorentzVector LL(0.,0.,0.,0.);
  TLorentzVector NuNu(0.,0.,0.,0.);
  TLorentzVector Nu1(0.,0.,0.,0.);
  TLorentzVector Nu2(0.,0.,0.,0.);
  TLorentzVector W1(0.,0.,0.,0.);
  TLorentzVector W2(0.,0.,0.,0.);
  TLorentzVector Higgs(0.,0.,0.,0.);
  TLorentzVector J1(0.,0.,0.,0.);
  TLorentzVector J2(0.,0.,0.,0.);

  //Getting some values to select the events
  unsigned ncleanjet{*nCleanJet->Get()};
  unsigned nlep{*nLepton->Get()};
  float Pmet_pt{*PuppiMET_pt->Get()};
  float Pmet_phi{*PuppiMET_phi->Get()};

  //Conditions to select the event
  if(ncleanjet>=2 && nlep>1){
	 
    //STEP-1
    //4-vectors of the leptons
    //Select one electron and one muon
    int muons = 0;
    int electrons = 0;
    int lep1 = 0;
    int lep2 = 0;
	  
    // Loop over muons and electrons
    for (unsigned int ilep = 0; ilep<nlep; ilep++){
     if (abs(Lepton_pdgId->At(ilep)) == 13){
    	++muons;
    	if (muons == 1 && Lepton_pt->At(ilep) > 13){
    	  L1.SetPtEtaPhiM(Lepton_pt->At(ilep), Lepton_eta->At(ilep), Lepton_phi->At(ilep), 0.0); //Muon
        lep1 = Lepton_pdgId->At(ilep);
    	}
      }
      if (abs(Lepton_pdgId->At(ilep)) == 11){
    	++electrons;
    	if (electrons == 1 && Lepton_pt->At(ilep) > 13){
    	  L2.SetPtEtaPhiM(Lepton_pt->At(ilep), Lepton_eta->At(ilep), Lepton_phi->At(ilep), 0.0); //Electron
        lep2 = Lepton_pdgId->At(ilep);
    	}
     }
    }

    if (muons<1 || electrons<1){
 

    	std::vector<std::pair<double, double>> weights = weight.computeWeights({higgs, jet1, jet2});

    
    	return (double)weights.back().first;
    
    }
    
  }
  //End if(nCleanJet>=2 && nLepton>1)
  else return -9999; 
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
