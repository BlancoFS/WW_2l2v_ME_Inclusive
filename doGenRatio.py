import uproot
import ROOT
import pandas as pd
import numpy as np



####### SAMPLE

baseDir = '/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/src/PlotsConfigurations/Configurations/WW/Full2016_v7/WWj_Full2016_v7/plotWWj_2016'
path = '/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/src/PlotsConfigurations/Configurations/WW/Full2016_v7/WWj_Full2016_v7/rootFile/plots_WWj_2016_NNLOPS.root'


### Open file

file = uproot.open(path)


### Loop over cuts

for cutName in file.keys():
  
  ### Loop over variables
  
  for variableName in file[cutName].keys():
    
    cName = 'cratio_' + cutName + '_' + variableName
    
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
    
    hist_WW = file[cutName][variableName]['histo_WW']
    hist_WW_NNLOPS_1 = file[cutName][variableName]['histo_WW_NNLOPS;1']
    hist_WW_NNLOPS_2 = file[cutName][variableName]['histo_WW_NNLOPS;2']
    hist_WW_NNLOPS = hist_WW_NNLOPS_1 + hist_WW_NNLOPS_2
    hist_WW_NLO = file[cutName][variableName]['histo_WW_NLO']
    
    hist_WW.SetMarkerSize(0)
    hist_WW_NNLOPS.SetMarkerSize(0)
    hist_WW_NLO.SetMarkerSize(0)
    
    hist_WW.SetLineWidth(3)
    hist_WW_NNLOPS.SetLineWidth(3)
    hist_WW_NLO.SetLineWidth(3)
    
    hist_WW.SetLineStyle(0)
    hist_WW_NNLOPS.SetLineStyle(0)
    hist_WW_NLO.SetLineStyle(0)
    
    hist_WWSetLineColor(ROOT.kBlue)
    hist_WW_NNLOPSSetLineColor(ROOT.kBlack)
    hist_WW_NLOSetLineColor(ROOT.kBlue+3)

    hist_WW.Scale(1/hist_WW.Integral())
    hist_WW_NNLOPS.Scale(1/hist_WW_NNLOPS.Integral())
    hist_WW_NLO.Scale(1/hist_WW_NLO.Integral())

    hist_WW.Draw('HIST')
    hist_WW_NNLOPS.Draw('same')
    hist_WW_NLO.Draw('same')
    
    legend = ROOT.TLegend(0.95, 0.9, 0.75, 0.78)
    legend.AddEntry(hist_WW, "WW NLO to NNLO+NNLL")
    legend.AddEntry(hist_WW_NNLOPS, "WW NNLOPS")
    legend.AddEntry(hist_WW_NLO, "WW NLO")
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.Draw()
    
    rcanvas.cd(2)
    
    ratioPad = rcanvas.GetPad(2)
    ratioPad.SetPad(0.,0.,1.,0.31)
    
    ratioPad.SetFillStyle(4000)
    ratioPad.SetBottomMargin(0.2)
    
    hRatio = hist_WW.Clone()
    hRatio2 = hist_WW_NLO.Clone()
    hRatio.SetTitle(" ")
    
    hRatio.GetXaxis().SetLabelSize(0.1)
    hRatio.GetXaxis().SetTitleSize(0.1)
    hRatio.GetXaxis().SetTitleOffset(.85)
    hRatio.GetXaxis().SetTitle(variableName.split(';')[0])
    
    hRatio.GetYaxis().SetLabelSize(0.07)
    hRatio.GetYaxis().SetTitleSize(0.1)
    hRatio.GetYaxis().SetTitleOffset(0.5)
    hRatio.GetYaxis().SetTitle("MC/NNLOPS")
    hRatio.GetYaxis().SetRangeUser(0.5,1.5)

    
    hRatio.Divide(hist_WW_NNLOPS)
    hRatio2.Divide(hist_WW_NNLOPS)
    
    hRatio.Draw()
    hRatio2.Draw("same")
    
    Xmax = hRatio.GetXaxis().GetXmax()
    Xmin = hRatio.GetXaxis().GetXmin()
    
    l = ROOT.TLine(Xmin, 1, Xmax, 1)
    l.SetLineColor(1) 
    l.Draw("same") 
    
    rcanvas.Draw()
    rcanvas.SaveAs(baseDir + cName + '.png')
    
    rcanvas.cd(1)
                    
    plotPad.SetLogy()
    rcanvas.Update()
    rcanvas.Draw()
    
    rcanvas.SaveAs(baseDir + 'log_' + cName + '.png')
    
    




