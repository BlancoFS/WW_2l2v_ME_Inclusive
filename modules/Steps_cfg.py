#!/usr/bin/env python

import os
from pprint import pprint



steps = {

'WPWPVarsDNN' : { 
                  'isChain'  : True ,
                  'do4MC'    : True  ,
                  'do4Data'  : True ,
                  'outputbranchsel': os.getenv('CMSSW_BASE') + '/src/LatinoAnalysis/NanoGardener/python/data/WPWPVarsDNN_branches.txt',
                  'selection': '"(nLepton >= 2 && \
                               Alt$(Lepton_pt[2],0) < 10. && \
                               Lepton_pt[0] > 25 && \
                               Lepton_pt[1] > 15)"',
               },

}
