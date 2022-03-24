import ROOT
import math 
import numpy
import ctypes
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.modules.common.collectionMerger import 


import os.path



class WW_MEreweighter(Module):
    def __init__(self, sample):
        print '####################', sample
        self.sample = sample
        self.cmssw_base = os.getenv('CMSSW_BASE')
        self.cmssw_arch = os.getenv('SCRAM_ARCH')

        ROOT.gSystem.AddIncludePath("-I"+self.cmssw_base+"/interface/")
        ROOT.gSystem.AddIncludePath("-I"+self.cmssw_base+"/src/")
        ROOT.gSystem.Load(self.cmssw_base+"/lib/libmomemta.so.1.0.0")
        ROOT.gSystem.Load(self.cmssw_base+"/lib/libmomemta.so.1")
        ROOT.gSystem.Load(self.cmssw_base+"/lib/libmomemta.so")

        try:
            ROOT.gROOT.LoadMacro(self.cmssw_base+'/src/LatinoAnalysis/Gardener/python/variables/recoMoMEMta_WW.C+g')
        except RuntimeError:
            ROOT.gROOT.LoadMacro(self.cmssw_base+'/src/LatinoAnalysis/Gardener/python/variables/recoMoMEMta_WW.C++g')
            
      
    def beginJob(self):
        pass
    def endJob(self):
        pass
      
    
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        
        self.inputFile = inputFile

        self.out = wrappedOutputTree                
        self.out.branch('D_Top','F')
    
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
      
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        
        Jets = Collection(event, 'Jet')
        CleanJets = Collection(event, 'CleanJet')
        Leptons = Collection(event, 'Lepton')
        
        njet = 0
        nlep = Leptons._len
        
        for jet in range(CleanJets._len):
          if (CleanJets[jet].pt > 15):
            njet = njet + 1

          
        if (njet < 2 or nlep < 2):
          
          D_Top = 999
          
        else:
          
          D_Top = ROOT.recoMoMEMta_WW(Leptons[0].pt, Leptons[1].pt, Leptons[0].phi, Leptons[1].phi, Leptons[0].eta, Leptons[1].eta, event.MET_pt, CleanJets[0].pt, CleanJets[1].pt, CleanJets[0].phi, CleanJets[1].phi, CleanJets[0].eta, CleanJets[1].eta)

        self.out.fillBranch('D_Top', D_Top)
        
        return True
    
    
    
    
    
    
    
    




