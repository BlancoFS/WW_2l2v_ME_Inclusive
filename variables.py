# variables

#variables['bdt'] = {      'name'  : 'VH2j_TMVAReader(Entry$)',# variable name
#                          'range' : (20, -1., 1.),            # variable range
#                          'xaxis' : 'BDT discriminant VH2j',  # x-axis name
#                          'fold'  : 3,                        # 0 = not fold (default), 1 = fold underflow bin, 2 = fold overflow bin, 3 = fold underflow and overflow
#                          'linesToAdd' : ['.L /afs/cern.ch/user/p/piedra/work/VH2jBDT/VH2j_TMVAReader.C+']}

#
# Centrality
#


variables['ml1j1']  = {  
    'name': 'mlj(CleanJet_pt[0], CleanJet_phi[0], CleanJet_eta[0], Lepton_pt[0], Lepton_phi[0], Lepton_eta[0])',
    'range': (40, 0.0, 500),
    'xaxis': 'm_{lj}',
    'linesToAdd': ['.L /afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/src/PlotsConfigurations/Configurations/WW/Full2016_v7/extended/mlj.C+'],
    'fold': 3
 }

variables['ml2j1']  = {
    'name': 'mlj(CleanJet_pt[0], CleanJet_phi[0], CleanJet_eta[0], Lepton_pt[1], Lepton_phi[1], Lepton_eta[1])',
    'range': (40, 0.0, 500),
    'xaxis': 'm_{lj}',
    'linesToAdd': ['.L /afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/src/PlotsConfigurations/Configurations/WW/Full2016_v7/extended/mlj.C+'],
    'fold': 3
 }

variables['ml1j2']  = {
    'name': 'mlj(CleanJet_pt[1], CleanJet_phi[1], CleanJet_eta[1], Lepton_pt[0], Lepton_phi[0], Lepton_eta[0])',
    'range': (40, 0.0, 500),
    'xaxis': 'm_{lj}',
    'linesToAdd': ['.L /afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/src/PlotsConfigurations/Configurations/WW/Full2016_v7/extended/mlj.C+'],
    'fold': 3
 }

variables['ml2j2']  = {
    'name': 'mlj(CleanJet_pt[1], CleanJet_phi[1], CleanJet_eta[1], Lepton_pt[1], Lepton_phi[1], Lepton_eta[1])',
    'range': (40, 0.0, 500),
    'xaxis': 'm_{lj}',
    'linesToAdd': ['.L /afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/src/PlotsConfigurations/Configurations/WW/Full2016_v7/extended/mlj.C+'],
    'fold': 3
 }


variables['DNN_0J'] = {      'name'  : 'dnn_0',
                           'range' : (30, 0., 1.),
                           'xaxis' : 'DNN OUTPUT with 0 jet',
                           'fold'  : 3
                    }

variables['DNN_1J'] = {      'name'  : 'dnn_1',                                                                                                                                                            
                            'range' : (30, 0., 1.),                                                                                                                                               
                           'xaxis' : 'DNN OUTPUT with 1 jet',                                                                                                                                              
                           'fold'  : 3                                                                                                                                                                     
                    } 

variables['DNN_2J'] = {      'name'  : 'dnn_2',                                                                                                                                                            
                           'range' : (30, 0., 1.),                                                                                                                                                         
                           'xaxis' : 'DNN OUTPUT with 2 jet',                                                                                                                                             
                           'fold'  : 3                                                                                                                                                                     
                    } 


variables['D_Top'] = {      'name'  : 'D_Top',
                           'range' : (30, 0., 1.),
                           'xaxis' : 'D_{Top}',
                           'fold'  : 3
                    }


#################


variables['ptWW'] = {      'name'  : 'pTWW',
                           'range' : (40, 0., 1000.),
                           'xaxis' : 'p_{T, WW} [GeV]',
                           'fold'  : 3}

variables['ptWW_400'] = {      'name'  : 'pTWW',
                               'range' : (40, 0., 400.),
                               'xaxis' : 'p_{T, WW} [GeV]',
                               'fold'  : 3}

variables['ptWW_200'] = {      'name'  : 'pTWW',
                               'range' : (30, 0., 200.),
                               'xaxis' : 'p_{T, WW} [GeV]',
                               'fold'  : 3}

variables['detall'] = { 'name'  : 'abs(Lepton_eta[0]-Lepton_eta[1])',
                          'range' : (20, 0., 5.),
                          'xaxis' : '|#Delta#eta_{l1l2}|',
                          'fold'  : 3}

variables['dphill'] = {   'name'  : 'abs(dphill)',     
                          'range' : (20, 0., 3.2),   
                          #'range' : (10, 0., 3.2),
                          'xaxis' : '#Delta#phi_{ll}',
                          'fold'  : 3}


variables['dphill_01'] = {   'name'  : 'abs(dphill)',
                          'range' : (10, 0., 1.),
                          'xaxis' : '#Delta#phi_{ll}',
                          'fold'  : 3}

variables['drll'] = {     'name'  : 'drll',
                          'range' : (30, 0., 5.),
                          #'range' : (10, 0., 2.5),
                          'xaxis' : '#DeltaR_{ll}',
                          'fold'  : 3}

variables['eta1'] = {     'name'  : 'Lepton_eta[0]',     
                          'range' : (40, -3.2, 3.2),   
                          #'range' : (15, -3.2, 3.2),
                          'xaxis' : '#eta 1st lepton',
                          'fold'  : 3}

variables['eta2'] = {     'name'  : 'Lepton_eta[1]',     
                          'range' : (40, -3.2, 3.2),
                          #'range' : (15, -3.2, 3.2),
                          'xaxis' : '#eta 2nd lepton',
                          'fold'  : 3}

variables['events'] = {   'name'  : '1',
                          'range' : (1, 0, 2),
                          'xaxis' : 'events',
                          'fold'  : 3}

variables['mll_top'] = {      'name'  : 'mll',
                          'range' : (20, 80., 200.),
                          'xaxis' : 'm_{ll} [GeV]',
                          'fold'  : 3}

variables['mll'] = { 'name'  : 'mll',
                     #'range' : (20, 0., 200.),
                     'range' : (40, 0., 400.),
                     'xaxis' : 'm_{ll} [GeV]',
                     'fold'  : 3}

variables['mll_DY'] = { 'name'  : 'mll',
                        #'range' : (20, 0., 200.),                                                                                                                                                
                        'range' : (40, 0., 80.),
                        'xaxis' : 'm_{ll} [GeV]',
                        'fold'  : 3}

variables['mpmet'] = {    'name'  : 'mpmet',      
                          'range' : (50, 0., 150.),  
                          'xaxis' : 'min. (proj. tk. E_{T}^{miss}, proj. E_{T}^{miss}) [GeV]', 
                          'fold'  : 3}

variables['mth'] = {      'name'  : 'mth',
                          #'range' : (40, 0., 200.),
                          'range' : (40, 0., 400.),
                          'xaxis' : 'm_{T}^{H} [GeV]',
                          'fold'  : 3}

variables['mth_top'] = {      'name'  : 'mth',
                              #'range' : (40, 0., 200.),                                                                                                                                                   
                              'range' : (40, 0., 200.),
                              'xaxis' : 'm_{T}^{H} [GeV]',
                              'fold'  : 3}

variables['mth_DY'] = {      'name'  : 'mth',
                             'range' : (20, 0., 40.),
                             'xaxis' : 'm_{T}^{H} [GeV]',
                             'fold'  : 3}

variables['mth_WW'] = {      'name'  : 'mth',
                             'range' : (12, 50., 200.),
                             'xaxis' : 'm_{T}^{H} [GeV]',
                             'fold'  : 3}

variables['mth_GGH'] = {      'name'  : 'mth',
                             'range' : (10, 60., 145.),
                             'xaxis' : 'm_{T}^{H} [GeV]',
                             'fold'  : 3}

variables['mtw1']  = {   'name': 'mtw1',            
                         'range' : (20,0,400),    
                         'xaxis' : 'm^{T}_{W1}  [GeV]',
                         'fold'  : 3
                     }

variables['mtw2']  = {   'name': 'mtw2',
                         'range' : (20,0,400),
                         'xaxis' : 'm^{T}_{W1}  [GeV]',
                         'fold'  : 3
                     }

variables['njet'] = {     'name'  : 'Sum$(CleanJet_pt>30)',     
                          'range' : (5, 0, 5),   
                          'xaxis' : 'number of jets',
                          'fold'  : 3}

variables['bjet'] = {     'name'  : 'Sum$(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.2217)',     
                          'range' : (5, 0, 5),   
                          'xaxis' : 'number of b jets',
                          'fold'  : 2}

variables['nvtx'] = {     'name'  : 'PV_npvsGood',      
                          'range' : (50, 0, 50),  
                          'xaxis' : 'number of vertices', 
                          'fold'  : 3}

variables['phi1'] = {     'name'  : 'Lepton_phi[0]',
                          'range' : (40, -3.2, 3.2),
                          'xaxis' : '#phi 1st lepton',
                          'fold'  : 3}

variables['phi2'] = {     'name'  : 'Lepton_phi[1]',
                          'range' : (40, -3.2, 3.2),
                          'xaxis' : '#phi 2nd lepton',
                          'fold'  : 3}

variables['pfmet'] = {    'name'  : 'MET_pt',     
                          'range' : (50, 0., 150.),   
                          'xaxis' : 'PF MET [GeV]',
                          'fold'  : 3}

variables['pt1'] = {      'name'  : 'Lepton_pt[0]',     
                          #'range' : (40, 0., 200.),   
                          'range' : (10, 20., 100.),
                          'xaxis' : 'p_{T} 1st lepton [GeV]',
                          'fold'  : 3}

variables['pt2'] = {      'name'  : 'Lepton_pt[1]',     
                          #'range' : (40, 0., 200.),   
                          'range' : (40, 13., 100.),
                          'xaxis' : 'p_{T} 2nd lepton [GeV]',
                          'fold'  : 3}

variables['ptll'] = {     'name'  : 'ptll',
                          #'range' : (40, 30., 200.),
                          'range' : (10, 30., 200.),
                          'xaxis' : 'p_{T}^{ll} [GeV]',
                          'fold'  : 3}

variables['puppimet'] = { 'name'  : 'PuppiMET_pt',
                          'range' : (50, 0., 150.),
                          'xaxis' : 'puppi MET [GeV]',
                          'fold'  : 3}

variables['rawmet'] = {   'name'  : 'RawMET_pt',     
                          'range' : (50, 0., 150.),   
                          'xaxis' : 'raw MET [GeV]',
                          'fold'  : 3}

variables['TkMET'] = {    'name'  : 'TkMET_pt',     
                          'range' : (50, 0., 150.),   
                          'xaxis' : 'tracker MET [GeV]',
                          'fold'  : 3}


## VARIABLES WITH 1 JET

variables['detal1j'] = {   'name'  : 'abs(Lepton_eta[0] - CleanJet_eta[0]) - 9999.0 * (Sum$(CleanJet_pt>30) == 0)',
                          'range' : (20, 0., 5.),
                          #'range' : (10, 3., 4.),                                                                                                                                                         
                          'xaxis' : '|#Delta#eta_{lj}|',
                          'fold'  : 2}

variables['dphil1j'] = {   'name'  : 'abs(Lepton_phi[0] - CleanJet_phi[0]) - 9999.0 * (Sum$(CleanJet_pt>30) == 0)',
                          'range' : (20, 0., 5.),
                          #'range' : (10, 3., 4.),                                                                                                                                                         
                          'xaxis' : '|#Delta#phi_{lj}|',
                          'fold'  : 2}

variables['detal2j'] = {   'name'  : 'abs(Lepton_eta[1] - CleanJet_eta[0]) - 9999.0 * (Sum$(CleanJet_pt>30) == 0)',
                          'range' : (20, 0., 5.),
                          #'range' : (10, 3., 4.),                                                                                                                                                         
                          'xaxis' : '|#Delta#eta_{lj}|',
                          'fold'  : 2}

variables['dphil2j'] = {   'name'  : 'abs(Lepton_phi[1] - CleanJet_phi[0]) - 9999.0 * (Sum$(CleanJet_pt>30) == 0)',
                          'range' : (20, 0., 5.),
                          #'range' : (10, 3., 4.),                                                                                                                                                         
                          'xaxis' : '|#Delta#phi_{lj}|',
                          'fold'  : 2}

## VARIABLES WITH 2 JETS

variables['detajj'] = {   'name'  : 'detajj - 9999.0 * (Sum$(CleanJet_pt>30) != 2)',
                          'range' : (20, 0., 5.),
                          #'range' : (10, 3., 4.),                                                                                                                                                         
                          'xaxis' : '|#Delta#eta_{jj}|',
                          'fold'  : 2}

variables['drjj'] = {     'name'  : 'sqrt((CleanJet_eta[0]-CleanJet_eta[1])**2 + (CleanJet_phi[0]-CleanJet_phi[1])**2)  - 9999.0 * (Sum$(CleanJet_pt>30) != 2)',
                          'range' : (20, 0., 8.),
                          'xaxis' : '#Delta R_{jj}',
                          'fold'  : 2}

variables['mjj'] = {   'name'  : 'mjj  - 9999.0 * (Sum$(CleanJet_pt>30) != 2)',
                       'range' : (40, 0., 400.),
                       #'range' : (10, 3., 4.),                                                                                                                                                         
                       'xaxis' : 'm_{jj}',
                       'fold'  : 2}

variables['dphijj'] = {   'name'  : 'CleanJet_phi[0] - CleanJet_phi[1] - 9999.0 * (Sum$(CleanJet_pt>30) != 2)',
                          'range' : (25, -3.2, 3.2),
                          #'range' : (10, 0., 3.2),                                                                                                                                                        
                          'xaxis' : '#Delta#phi_{jj}',
                          'fold'  : 2}

variables['Ctot'] = {     'name': 'log((abs(2*Lepton_eta[0]-CleanJet_eta[0]-CleanJet_eta[1])+abs(2*Lepton_eta[1]-CleanJet_eta[0]-CleanJet_eta[1]))/detajj) - 9999.0 * (Sum$(CleanJet_pt>30) != 2)',
                          'range' : (25, -4, 6),
                          #'range' : (10, -4, 2),                                                                                                                                                          
                          'xaxis' : 'Ctot',
                          'fold'  : 2}
