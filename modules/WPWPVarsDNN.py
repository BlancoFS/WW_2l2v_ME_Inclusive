import ROOT

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class WPWPVarsDNN(Module):
    def __init__(self):

        self.bVetoCut = 0.2217
        self.metpt = 'event.MET_pt'
        self.metphi = 'event.MET_phi'
        self.lep1pt = 'event.Lepton_pt[0]'
        self.lep2pt = 'event.Lepton_pt[1]'
        self.lep1eta = 'event.Lepton_eta[0]'
        self.lep2eta = 'event.Lepton_eta[1]'
        self.lep1phi = 'event.Lepton_phi[0]'
        self.lep2phi = 'event.Lepton_phi[1]'
        

    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        
        self.inputFile = inputFile

        self.out = wrappedOutputTree        
        self.out.branch('metpt','F')
        self.out.branch('metphi','F')
        self.out.branch('lep1pt','F')
        self.out.branch('lep2pt','F')
        self.out.branch('lep1eta','F')
        self.out.branch('lep2eta','F')
        self.out.branch('lep1phi','F')
        self.out.branch('lep2phi','F')
        self.out.branch('mll','F')
        self.out.branch('mth','F')
        self.out.branch('dphill','F')
        self.out.branch('costheta','F')
        self.out.branch('mww','F')
        self.out.branch('ptll','F')
        self.out.branch('drll','F')
        self.out.branch('Weight_LL','F')
        self.out.branch('Weight_TL','F')
        self.out.branch('Weight_LT','F')
        self.out.branch('Weight_TT','F')

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        
        GenPart = Collection(event, 'GenPart')
        nParticles = GenPart._len
        
        PolWeight_LL = 0.0
        PolWeight_TL = 0.0
        PolWeight_LT = 0.0
        PolWeight_TT = 0.0

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
        number_tau  = 0
        
        for p in range(nParticles):
          
          if (GenPart[p].pdgId==11 and GenPart[GenPart[p].genPartIdxMother].pdgId==-24 and GenPart[GenPart[GenPart[p].genPartIdxMother].genPartIdxMother].pdgId==15):
              pos_wm = GenPart[p].genPartIdxMother
              number_elec = number_elec + 1
              vector_lm.SetPtEtaPhiM(GenPart[p].pt, GenPart[p].eta, GenPart[p].phi, 0.0)
              genlm.SetCoordinates(GenPart[p].pt, GenPart[p].eta, GenPart[p].phi, vector_lm.E())
          
          elif (GenPart[p].pdgId==-11 and GenPart[GenPart[p].genPartIdxMother].pdgId==24 and GenPart[GenPart[GenPart[p].genPartIdxMother].genPartIdxMother].pdgId==-15):
              pos_wp = GenPart[p].genPartIdxMother
              number_elec = number_elec + 1
              vector_lp.SetPtEtaPhiM(GenPart[p].pt, GenPart[p].eta, GenPart[p].phi, 0.0)
              genlp.SetCoordinates(GenPart[p].pt, GenPart[p].eta, GenPart[p].phi, vector_lp.E())
          
          elif (GenPart[p].pdgId==13 and GenPart[GenPart[p].genPartIdxMother].pdgId==-24 and GenPart[GenPart[GenPart[p].genPartIdxMother].genPartIdxMother].pdgId==15):
              pos_wm = GenPart[p].genPartIdxMother
              number_muon = number_muon + 1
              vector_lm.SetPtEtaPhiM(GenPart[p].pt, GenPart[p].eta, GenPart[p].phi, 0.0)
              genlm.SetCoordinates(GenPart[p].pt, GenPart[p].eta, GenPart[p].phi, vector_lm.E())
          
          elif (GenPart[p].pdgId==-13 and GenPart[GenPart[p].genPartIdxMother].pdgId==24 and GenPart[GenPart[GenPart[p].genPartIdxMother].genPartIdxMother].pdgId==-15):
              pos_wp = GenPart[p].genPartIdxMother
              number_muon = number_muon + 1
              vector_lp.SetPtEtaPhiM(GenPart[p].pt, GenPart[p].eta, GenPart[p].phi, 0.0)
              genlp.SetCoordinates(GenPart[p].pt, GenPart[p].eta, GenPart[p].phi, vector_lp.E())
                
          elif (GenPart[p].pdgId==15 and GenPart[GenPart[p].genPartIdxMother].pdgId==-24):
              pos_wm = GenPart[p].genPartIdxMother
              number_tau = number_tau + 1
              vector_lm.SetPtEtaPhiM(GenPart[p].pt, GenPart[p].eta, GenPart[p].phi, 0.0)
              genlm.SetCoordinates(GenPart[p].pt, GenPart[p].eta, GenPart[p].phi, vector_lm.E())
          
          elif (GenPart[p].pdgId==-15 and GenPart[GenPart[p].genPartIdxMother].pdgId==24):
              pos_wp = GenPart[p].genPartIdxMother
              number_tau = number_tau + 1
              vector_lp.SetPtEtaPhiM(GenPart[p].pt, GenPart[p].eta, GenPart[p].phi, 0.0)
              genlp.SetCoordinates(GenPart[p].pt, GenPart[p].eta, GenPart[p].phi, vector_lp.E())
          
          
          if (GenPart[p].pdgId==-12 nd GenPart[GenPart[p].genPartIdxMother].pdgId==-24 and GenPart[GenPart[GenPart[p].genPartIdxMother].genPartIdxMother].pdgId==15):
              vector_num.SetPtEtaPhiM(GenPart[p].pt, GenPart[p].eta, GenPart[p].phi, 0.0)
              gennum.SetCoordinates(GenPart[p].pt, GenPart[p].eta, GenPart[p].phi, vector_num.E())
          
          elif (GenPart[p].pdgId==12 and GenPart[GenPart[p].genPartIdxMother].pdgId==24 and GenPart[GenPart[GenPart[p].genPartIdxMother].genPartIdxMother].pdgId==-15):
              vector_nup.SetPtEtaPhiM(GenPart[p].pt, GenPart[p].eta, GenPart[p].phi, 0.0)
              gennup.SetCoordinates(GenPart[p].pt, GenPart[p].eta, GenPart[p].phi, vector_nup.E())
          
          elif (GenPart[p].pdgId==-14 and GenPart[GenPart[p].genPartIdxMother].pdgId==-24 and GenPart[GenPart[GenPart[p].genPartIdxMother].genPartIdxMother].pdgId==15):
              vector_num.SetPtEtaPhiM(GenPart[p].pt, GenPart[p].eta, GenPart[p].phi, 0.0)
              gennum.SetCoordinates(GenPart[p].pt, GenPart[p].eta, GenPart[p].phi, vector_num.E())
          
          elif (GenPart[p].pdgId==14 and GenPart[GenPart[p].genPartIdxMother].pdgId==24 and GenPart[GenPart[GenPart[p].genPartIdxMother].genPartIdxMother].pdgId==-15):
              vector_nup.SetPtEtaPhiM(GenPart[p].pt, GenPart[p].eta, GenPart[p].phi, 0.0)
              gennup.SetCoordinates(GenPart[p].pt, GenPart[p].eta, GenPart[p].phi, vector_nup.E())
                
          elif (GenPart[p].pdgId==-16 and GenPart[GenPart[p].genPartIdxMother].pdgId==-24):
              vector_num.SetPtEtaPhiM(GenPart[p].pt, GenPart[p].eta, GenPart[p].phi, 0.0)
              gennum.SetCoordinates(GenPart[p].pt, GenPart[p].eta, GenPart[p].phi, vector_num.E())
          
          elif (GenPart[p].pdgId==16 and GenPart[GenPart[p].genPartIdxMother].pdgId==24):
              vector_nup.SetPtEtaPhiM(GenPart[p].pt, GenPart[p].eta, GenPart[p].phi, 0.0)
              gennup.SetCoordinates(GenPart[p].pt, GenPart[p].eta, GenPart[p].phi, vector_nup.E())
          
          
        vector_Wp.SetPtEtaPhiM(GenPart[pos_wp].pt, GenPart[pos_wp].eta, GenPart[pos_wp].phi, GenPart[pos_wp].mass) # W plus
        genWp.SetCoordinates(GenPart[pos_wp].pt, GenPart[pos_wp].eta, GenPart[pos_wp].phi, vector_Wp.E())
            
        
        vector_Wm.SetPtEtaPhiM(GenPart[pos_wm].pt, GenPart[pos_wm].eta, GenPart[pos_wm].phi, GenPart[pos_wm].mass) # W minus
        genWm.SetCoordinates(GenPart[pos_wm].pt, GenPart[pos_wm].eta, GenPart[pos_wm].phi, vector_Wm.E())
          
        if (((number_elec==0 or number_muon==0) and number_tau==0) or pos_wp==999 or pos_wm==999):
            PolWeight_LL = 1.0
            PolWeight_TL = 1.0
            PolWeight_LT = 1.0
            PolWeight_TT = 1.0
        else:
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
            
            #f0_m = 0.26
            #fL_m = 0.48
            #fR_m = 0.25
            #fT_m = fL_m + fR_m
            
            f0_p = 0.271
            fT_p = 0.729
            fL_p = fT_p/1.52
            fR_p = fT_p - fL_p
            
            f0_m = 0.271
            fT_m = 0.729
            fL_m = fT_m/1.52
            fR_m = fT_m - fL_m
            
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
            
            PolWeight_LL = (weight_f0_Wp/weight_total_Wp)*(weight_f0_Wm/weight_total_Wm)
            PolWeight_TL = (weight_fT_Wp/weight_total_Wp)*(weight_f0_Wm/weight_total_Wm)
            PolWeight_LT = (weight_f0_Wp/weight_total_Wp)*(weight_fT_Wm/weight_total_Wm)
            PolWeight_TT = (weight_fT_Wp/weight_total_Wp)*(weight_fT_Wm/weight_total_Wm)
            
        theta_ll = ROOT.Math.VectorUtil.Angle(genlp, genlm)

        costhetall = ROOT.Math.cos(theta_ll)
        
        ### Fill Branch    
        
        self.out.fillBranch('metpt', self.metpt)
        self.out.fillBranch('metphi', self.metphi)
        self.out.fillBranch('lep1pt', self.lep1pt)
        self.out.fillBranch('lep2pt', self.lep2pt)
        self.out.fillBranch('lep1eta', self.lep1eta)
        self.out.fillBranch('lep2eta', self.lep2eta)
        self.out.fillBranch('lep1phi', self.lep1phi)
        self.out.fillBranch('lep2phi', self.lep2phi)
        self.out.fillBranch('mll', event.mll)
        self.out.fillBranch('mth', event.mth)
        self.out.fillBranch('dphill', event.dphill)
        self.out.fillBranch('costheta', costhetall)
        self.out.fillBranch('mww', event.Mww)
        self.out.fillBranch('ptll', event.ptll)
        self.out.fillBranch('drll', event.drll)
        self.out.fillBranch('Weight_LL', PolWeight_LL)
        self.out.fillBranch('Weight_TL', PolWeight_TL)
        self.out.fillBranch('Weight_LT', PolWeight_LT)
        self.out.fillBranch('Weight_TT', PolWeight_TT)
        

        return True
