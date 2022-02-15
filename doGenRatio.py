import uproot
import ROOT
import pandas as pd
import numpy as np



####### SAMPLE


path = '/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/src/PlotsConfigurations/Configurations/WW/Full2016_v7/WWj_Full2016_v7/rootFile/plots_WWj_2016_NNLOPS.root'


### Open file

file = uproot.open(path)


### Loop over cuts

for cutName in file.keys():
  
  ### Loop over variables
  
  for variableName in file[cutName].keys():
    
    cName = 'c_' + cutName + '_' + variableName
    
    rcanvas = ROOT.TCanvas("r"+cName, "r"+cName, 900, 800)
    rcanvas.SetRightMargin(0.24)
    
    rcanvas.Divide(1,2)
    rcanvas.cd(1)
    
    plotPad = rcanvas.GetPad(1)
    plotPad.SetPad(0.,0.2,1.,1.)
    
    hist.GetXaxis().SetLabelSize(0.)
    hist.GetXaxis().SetTitleSize(0.)
    
    hist_mc.GetXaxis().SetLabelSize(0.)
    hist_mc.GetXaxis().SetTitleSize(0.)
    
    








