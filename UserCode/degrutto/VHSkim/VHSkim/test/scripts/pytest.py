#! /usr/bin/env python


import sys
import os
import commands
import string

#os.system('mkdir -p %s' % destFolder)

def fib(n):
	os.system("sed -i 's/XXX/%s/g' ntupleZmm_cfg.py" % n)
	os.system('Ntupler ntupleZmm_cfg.py')
	os.system("sed -i 's/%s/XXX/g' ntupleZmm_cfg.py" % n)





L = [
"ZH115",
"Elev11v2",
"Muv2",
"TTbar",
"Z0Jets",
"Z1Jets0-100",
"Z1Jets100-300",
"Z1Jets300-800",
"Z1Jets800-1600",
"Z2Jets0-100",
"Z2Jets100-300",
"Z2Jets300-800",
"Z2Jets800-1600",
"Z3Jets0-100",
"Z3Jets100-300",
"Z3Jets300-800",
"Z3Jets800-1600",
"Z4Jets0-100",
"Z4Jets100-300",
"Z4Jets300-800",
"Z4Jets800-1600",
"Z5Jets0-100",
"Z5Jets100-300",
"Z5Jets300-800",
"Z5Jets800-1600",
"ZJetMad",
"ZZ",
"Zbb0Jets",
"Zbb1Jets",
"Zbb2Jets",
"Zbb3Jets",
"Zcc0Jets",
"Zcc1Jets",
"Zcc2Jets",
"Zcc3Jets",
"WW",
"WZ",
"Wjets",
"WlnHbb115",
"T-s",
"T-t",
"T-tW",
"ZllHbb"
]



for item in L:
        fib(item)

