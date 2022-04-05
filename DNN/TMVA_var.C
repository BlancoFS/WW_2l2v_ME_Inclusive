// To compile 
// root -l
// gSystem->Load("libLatinoAnalysisMultiDraw.so")
// .L hww_VBF_MYmvaBDTG.C+ 

#include <TMVA/Reader.h>
#include "TLorentzVector.h"
#include <TTree.h>
#include "TSystem.h"
#include "TROOT.h"
#include "TString.h"
#include "TFile.h"


using namespace std;

namespace multidraw {
  extern thread_local TTree* currentTree;
}

TMVA::Reader* reader = new TMVA::Reader();
bool initialized = false;
TString name_temp = "";
//TTree* latino = 0;

 float mll; 
 float lep1pt;
 float lep2pt;
 float lep1eta; 
 float lep2eta;
 float lep1phi; 
 float lep2phi; 
 float drll; 
 float ptll; 
 float mtw1; 
 float mtw2; 
 float dphill;
 float cos;
 float metpt;
 float metphi;
 float mpmet;
 float mth;  
 float rap1;  
 float rap2;  
 float dphilmet; 
 float dphilmet1; 
 float dphilmet2; 

 float loc_mll; 
 vector<float>* loc_leppt;
 vector<float>* loc_lepeta; 
 vector<float>* loc_lepphi; 
 float loc_drll; 
 float loc_ptll; 
 float loc_mtw1; 
 float loc_mtw2; 
 float loc_dphill;
 float loc_metpt;
 float loc_metphi;
 float loc_mpmet;
 float loc_mth;   
 float loc_dphilmet; 
 float loc_dphilmet1; 
 float loc_dphilmet2; 


vector<float>* loc0_ptj = 0;
vector<float>* loc0_etaj = 0;
vector<float>* loc0_phij = 0;
vector<float>* loc0_ptl = 0;
vector<float>* loc0_etal = 0;
vector<float>* loc0_phil = 0;
vector<float>* loc0_qgl = 0;


void init_TMVA_var(TTree* tree, String name_){

 tree->SetBranchAddress("mll", &loc_mll); 
 tree->SetBranchAddress("Lepton_pt", &loc_leppt);
 tree->SetBranchAddress("Lepton_eta", &loc_lepeta); 
 tree->SetBranchAddress("Lepton_phi", &loc_lepphi); 
 tree->SetBranchAddress("Lepton_pdgId", &loc_leppdgId); 
 tree->SetBranchAddress("drll", &loc_drll); 
 tree->SetBranchAddress("ptll", &loc_ptll); 
 tree->SetBranchAddress("mtw1", &loc_mtw1); 
 tree->SetBranchAddress("mtw2", &loc_mtw2); 
 tree->SetBranchAddress("dphill", &loc_dphill);
 tree->SetBranchAddress("MET_pt", &loc_metpt);
 tree->SetBranchAddress("MET_phi", &loc_metphi);
 tree->SetBranchAddress("mpmet", &loc_mpmet);
 tree->SetBranchAddress("mth", &loc_mth);   
 tree->SetBranchAddress("dphilmet", &loc_dphilmet); 
 tree->SetBranchAddress("dphilmet1", &loc_dphilmet1); 
 tree->SetBranchAddress("dphilmet2", &loc_dphilmet2); 
  
 reader->AddVariable("mll", &mll);
 reader->AddVariable("lep1pt", &lep1pt);
 reader->AddVariable("lep2pt", &lep2pt);
 reader->AddVariable("lep1eta", &lep1eta);
 reader->AddVariable("lep2eta", &lep2eta);
 reader->AddVariable("lep1phi", &lep1phi);
 reader->AddVariable("lep2phi", &lep2phi);
 reader->AddVariable("drll", &drll);
 reader->AddVariable("ptll", &ptll);
 reader->AddVariable("mtw1", &mtw1);
 reader->AddVariable("mtw2", &mtw2);
 reader->AddVariable("dphill", &dphill);
 reader->AddVariable("costheta", &cos);
 reader->AddVariable("metpt", &MET_pt);
 reader->AddVariable("metphi", &MET_phi);
 reader->AddVariable("mpmet", &mpmet);
 reader->AddVariable("mth", &mth);
 reader->AddVariable("rap1", &rap1);
 reader->AddVariable("rap2", &rap2);
 reader->AddVariable("dphilmet", &dphilmet);
 reader->AddVariable("dphilmet1", &dphilmet1);
 reader->AddVariable("dphilmet2", &dphilmet2);

  TString dir    = "/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/src/PlotsConfigurations/Configurations/WW/Full2016_v7/WW_helicity/DNN/dataset/weights/";
  TString prefix = "TMVAClassification";
  
  TString methodName = TString(name_) + TString(" method");
  TString weightfile = dir + prefix + TString("_") + TString(name_) + TString(".weights.xml");
  reader->BookMVA( methodName, weightfile );

}


float TMVA_var(int entry, String name_){
  
//	cout << multidraw::currentTree->GetCurrentFile()->GetName() << endl;
//        gSystem->Load("libLatinoAnalysisMultiDraw.so");

	
	if(name_temp != multidraw::currentTree->GetCurrentFile()->GetName()){
		cout << "name_temp = " << name_temp << endl;
		name_temp = multidraw::currentTree->GetCurrentFile()->GetName();
		cout << "name_temp = " << name_temp << endl;
		initialized = false;
	}

	
  if (!initialized){
		//latino = (TTree*)gDirectory->Get("latino");
		delete reader;
		reader = new TMVA::Reader();
		init_TMVA_var(multidraw::currentTree, String name_);
		cout << "check init" << endl;	
		initialized = true;		
  }

	multidraw::currentTree->GetEntry(entry);

  mll        = loc_mll;  
  lep1pt     = *(loc_leppt)[0]; 
  lep2pt     = *(loc_leppt)[1]; 
  lep1eta    = *(loc_lepeta)[0];  
  lep2eta    = *(loc_lepeta)[1]; 
  lep1phi    = *(loc_lepphi)[0];  
  lep2phi    = *(loc_lepphi)[1];  
  drll       = loc_drll;  
  ptll       = loc_ptll;  
  mtw1       = loc_mtw1;  
  mtw2       = loc_mtw2;  
  dphill     = loc_dphill; 
  metpt      = loc_metpt; 
  metphi     = loc_metphi; 
  mpmet      = loc_mpmet; 
  mth        = loc_mth;    
  dphilmet   = loc_dphilmet;  
  dphilmet1  = loc_dphilmet1;  
  dphilmet2  = loc_dphilmet2;  
  
  TLorentzVector L1(0., 0., 0., 0.);
  TLorentzVector L2(0., 0., 0., 0.);

  // Lepton 1 is positive                                                                                                                                                                                                                   
  if (*(loc_leppdgId)[0] < 0){
    L1.SetPtEtaPhiM(lep1pt, lep1eta, lep1phi, 0.0);
    L2.SetPtEtaPhiM(lep2pt, lep2eta, lep2phi, 0.0);
  }else{
    L2.SetPtEtaPhiM(lep1pt, lep1eta, lep1phi, 0.0);
    L1.SetPtEtaPhiM(lep2pt, lep2eta, lep2phi, 0.0);
  }

  //double theta_1 = L1.Theta();                                                                                                                                                                                                            
  //double theta_2 = L2.Theta();                                                                                                                                                                                                            

  ROOT::Math::PtEtaPhiEVector lep1;
  ROOT::Math::PtEtaPhiEVector lep2;

  lep1.SetCoordinates(L1.Pt(), L1.Phi(), L1.Eta(), L1.E());
  lep2.SetCoordinates(L2.Pt(), L2.Phi(), L2.Eta(), L2.E());


  Double_t theta_ll = ROOT::Math::VectorUtil::Angle(lep1, lep2);

  cos = (float)ROOT::Math::cos(theta_ll);
                                                                                                                                                                                                      
  rap1 = 0.5*log( (L1.E()+L1.Pz()) / (L1.E()-L1.Pz()) );
  rap2 = 0.5*log( (L2.E()+L2.Pz()) / (L2.E()-L2.Pz()) );

	TString methodName = TString(name_) + TString(" method");
	float classifier = reader->EvaluateMulticlass(methodName);
	//cout << entry << " " << classifier << endl;
	return classifier;

}
