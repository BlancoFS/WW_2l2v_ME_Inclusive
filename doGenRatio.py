#!/usr/bin/env python
import uproot
import ROOT
import pandas as pd
import numpy as np
import os
import collections
import tdrstyle
import CMS_lumi
from collections import OrderedDict




if __name__=='__main__':

  print("\n")
  print("###########################")
  print("###########################")
  print("### COMPARE NLO/NNLOPS  ###")
  print("###########################")
  print("###########################")
  print("\n")

  ####### SAMPLE
  
  baseDir = '/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/src/PlotsConfigurations/Configurations/WW/Full2016_v7/WWj_Full2016_v7/plotWWj_2016/'
  path = '/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/src/PlotsConfigurations/Configurations/WW/Full2016_v7/WWj_Full2016_v7/rootFile/plots_WWj_2016_NNLOPS.root'
  variablesPath = '/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/src/PlotsConfigurations/Configurations/WW/Full2016_v7/WWj_Full2016_v7/variables_CR.py'
  
  ROOT.gROOT.SetBatch()
  ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")
  tdrstyle.setTDRStyle()
  
  variables = collections.OrderedDict()
  if os.path.exists(variablesPath) :
    handle = open(variablesPath,'r')
    exec(handle)
    handle.close()
    
  ### Open file
  
  file = uproot.open(path)
  
  print("Processing file:\n")
  print(path)
  
  ### Loop over cuts
  
  for cutName in file.keys():
    
    print("New cut being processed: " + cutName.split(';')[0])
    print("\n")
    ### Loop over variables
    
    for variableName, variable in variables.iteritems():
      
      print("New variable added: " + variableName.split(';')[0])
      print("\n")
      
      cName = 'cratio_' + cutName.split(';')[0] + '_' + variableName
      
      rcanvas = ROOT.TCanvas("r"+cName, "r"+cName, 900, 800)
      rcanvas.SetRightMargin(0.24)
      
      ROOT.gStyle.SetPaintTextFormat("5.3f")

      rcanvas.Divide(1,2)
      rcanvas.cd(1)
      
      plotPad = rcanvas.GetPad(1)
      plotPad.SetPad(0.,0.2,1.,1.)
      
      WW = file[cutName][variableName]['histo_WW']
      WW_NNLOPS_1 = file[cutName][variableName]['histo_WW_NNLOPS;1']
      WW_NNLOPS_2 = file[cutName][variableName]['histo_WW_NNLOPS;2']
      WW_NLO = file[cutName][variableName]['histo_WW_NLO']
      

      bins = variable['range']

      hist_WW = ROOT.TH1F(cName, cName, bins[0], bins[1], bins[2])
      hist_WW_NNLOPS = ROOT.TH1F(cName + 'NNLOPS', cName + 'NNLOPS', bins[0], bins[1], bins[2])
      hist_WW_NLO = ROOT.TH1F(cName + 'NLO', cName + 'NLO', bins[0], bins[1], bins[2])
      
      for i in range(1, len(WW.values)+1):
        hist_WW.SetBinContent(i, WW.values[i-1])
        hist_WW_NNLOPS.SetBinContent(i, WW_NNLOPS_1.values[i-1] + WW_NNLOPS_2.values[i-1])
        hist_WW_NLO.SetBinContent(i, WW_NLO.values[i-1])
        
        
      hist_WW.GetXaxis().SetLabelSize(0.)
      hist_WW.GetXaxis().SetTitleSize(0.)
      
      hist_WW_NNLOPS.GetXaxis().SetLabelSize(0.)
      hist_WW_NNLOPS.GetXaxis().SetTitleSize(0.)
      
      hist_WW_NLO.GetXaxis().SetLabelSize(0.)
      hist_WW_NLO.GetXaxis().SetTitleSize(0.)
    
      hist_WW.GetYaxis().SetTitleSize(0.07)
      hist_WW.GetYaxis().SetTitleOffset(0.9)
      hist_WW.GetYaxis().SetTitle("a.u.")
      
      hist_WW.SetMarkerSize(0)
      hist_WW_NNLOPS.SetMarkerSize(0)
      hist_WW_NLO.SetMarkerSize(0)
      
      hist_WW.SetLineWidth(2)
      hist_WW_NNLOPS.SetLineWidth(3)
      hist_WW_NLO.SetLineWidth(2)
      
      hist_WW.SetLineStyle(1)
      hist_WW_NNLOPS.SetLineStyle(1)
      hist_WW_NLO.SetLineStyle(2)
      
      hist_WW.SetLineColor(ROOT.kRed+1)
      hist_WW_NNLOPS.SetLineColor(ROOT.kBlack)
      hist_WW_NLO.SetLineColor(ROOT.kRed-1)
      
      hist_WW.Scale(1/hist_WW.Integral())
      hist_WW_NNLOPS.Scale(1/hist_WW_NNLOPS.Integral())
      hist_WW_NLO.Scale(1/hist_WW_NLO.Integral())

      hist_WW.SetMaximum(hist_WW.GetMaximum() + 0.2*hist_WW.GetMaximum())
      
      hist_WW.Draw('HIST')
      hist_WW_NNLOPS.Draw('hist same')
      hist_WW_NLO.Draw('hist same')
      
      legend = ROOT.TLegend(0.95, 0.9, 0.75, 0.78)
      legend.AddEntry(hist_WW_NNLOPS, "WW NNLOPS")
      legend.AddEntry(hist_WW, "WW NLO to NNLO+NNLL")
      legend.AddEntry(hist_WW_NLO, "WW NLO")
      legend.SetTextFont(42)
      legend.SetBorderSize(0)
      legend.SetFillColor(0)
      legend.Draw()
      
      lumi = 35.9

      CMS_lumi.cmsText = 'CMS'
      CMS_lumi.writeExtraText = True
      CMS_lumi.extraText = 'Preliminary'
      if lumi!=-1:
        CMS_lumi.lumi_13TeV = "%0.1f fb^{-1}" % (lumi)
      else:
        CMS_lumi.lumi_13TeV = ""
      CMS_lumi.CMS_lumi(plotPad, 4, 11)

      rcanvas.Modified()
      rcanvas.Update()

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
      hRatio.GetYaxis().SetTitle("MC / NNLOPS")
      hRatio.GetYaxis().SetRangeUser(0.5,1.5)
      
      
      hRatio.Divide(hist_WW_NNLOPS)
      hRatio2.Divide(hist_WW_NNLOPS)
      
      hRatio.SetMarkerSize(0)
      hRatio2.SetMarkerSize(0)
      hRatio.SetLineStyle(1)
      hRatio2.SetLineStyle(2)
      hRatio.SetLineWidth(3)
      hRatio2.SetLineWidth(3)
      hRatio.SetLineColor(ROOT.kRed+1)
      hRatio2.SetLineColor(ROOT.kRed-1)

      hRatio.Draw("HIST")
      hRatio2.Draw("HIST same")
      
      Xmax = hRatio.GetXaxis().GetXmax()
      Xmin = hRatio.GetXaxis().GetXmin()
      
      l = ROOT.TLine(Xmin, 1, Xmax, 1)
      l.SetLineColor(1) 
      l.SetLineWidth(3)
      l.Draw("same") 
      
      rcanvas.Draw()
      rcanvas.SaveAs(baseDir + cName + '.png')
      
      rcanvas.cd(1)
      
      plotPad.SetLogy()
      rcanvas.Update()
      rcanvas.Draw()
      
      rcanvas.SaveAs(baseDir + 'log_' + cName + '.png')
      
