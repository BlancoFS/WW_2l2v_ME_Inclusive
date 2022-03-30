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
        self.out.branch('metsphi','F')
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
        
        for iPart in range(nParticles):
          
          if (particles[p]==11 and particles[pos_mother[p]]==-24):
            pos_wm = pos_mother[p]
            number_elec = number_elec + 1
            vector_lm.SetPtEtaPhiM(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], 0.0)
            genlm.SetCoordinates(GenPart_pt.values[p], GenPart_eta.values[p], GenPart_phi.values[p], vector_lm.E())
          
          
          
          
    
    
    
            if (CleanJets[iJet].pt > 20. and ROOT.TMath.Abs(CleanJets[iJet].eta) < 2.5 and Jets[CleanJets[iJet].jetIdx].btagDeepB > self.bVetoCut):
                nbtagged += 1
            else:
                continue
                
        if 'TTTo2L2Nu' in str(self.inputFile):
            MCweight = (event.topGenPt * event.antitopGenPt > 0.) * (ROOT.TMath.Sqrt(ROOT.TMath.Exp(-0.158631 + 2.00214e-04*event.topGenPt - 3.09496e-07*event.topGenPt*event.topGenPt + 34.93/(event.topGenPt+135.633)) * ROOT.TMath.Exp(-0.158631 + 2.00214e-04*event.antitopGenPt - 3.09496e-07*event.antitopGenPt*event.antitopGenPt + 34.93/(event.antitopGenPt+135.633)))) + (event.topGenPt * event.antitopGenPt <= 0.)
        elif '_WWTo2L2Nu_' in str(self.inputFile):
            MCweight = event.nllW
        elif '_M-50' in str(self.inputFile): #DY high mass
            MCweight = (0.876979+event.gen_ptll*(4.11598e-03)-(2.35520e-05)*event.gen_ptll*event.gen_ptll)*(1.10211 * (0.958512 - 0.131835*ROOT.TMath.Erf((event.gen_ptll-14.1972)/10.1525)))*(event.gen_ptll<140)+0.891188*(event.gen_ptll>=140)
        elif '_M-10to50' in str(self.inputFile): #DY low mass
            MCweight = (8.61313e-01+event.gen_ptll*4.46807e-03-1.52324e-05*event.gen_ptll*event.gen_ptll)*(1.08683 * (0.95 - 0.0657370*ROOT.TMath.Erf((event.gen_ptll-11.)/5.51582)))*(event.gen_ptll<140)+1.141996*(event.gen_ptll>=140)
        elif '50_HT' in str(self.inputFile): #DY low mass
            MCweight = (8.61313e-01+event.gen_ptll*4.46807e-03-1.52324e-05*event.gen_ptll*event.gen_ptll)*(1.08683 * (0.95 - 0.0657370*ROOT.TMath.Erf((event.gen_ptll-11.)/5.51582)))*(event.gen_ptll<140)+1.141996*(event.gen_ptll>=140)
        elif (('GluGluWWTo2' in str(self.inputFile)) or ('GluGluToWW' in str(self.inputFile))):
            MCweight = 1.53/1.4
        elif (('ZZTo2L2Nu' in str(self.inputFile)) or ('ZZTo2L2Q' in str(self.inputFile)) or ('ZZTo4L' in str(self.inputFile)) or ('WZTo2L2Q' in str(self.inputFile))):
            MCweight = 1.11
        else:
            MCweight = 1.0

        self.out.fillBranch('nbtaggedJets',nbtagged)
        self.out.fillBranch('metpt',eval(self.metpt))
        self.out.fillBranch('metsig',eval(self.metsig))
        self.out.fillBranch('lep1pt',eval(self.lep1pt))
        self.out.fillBranch('lep2pt',eval(self.lep2pt))
        self.out.fillBranch('specialMCWeigths',MCweight)

        return True
