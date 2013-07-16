#! /usr/bin/env python


import sys
import os
import commands
import string

#os.system('mkdir -p %s' % destFolder)

def fib(n):
	os.system("sed -i 's/XXX/%s/g' NewNtuplerZll.py" % n)
	os.system('./NewNtuplerZll.py')
	os.system("sed -i 's/%s/XXX/g' NewNtuplerZll.py" % n)
        print ("finished processing %s" % n)




L = [
"TTbar",
"ZJetsMad",
"ZZ",
"WJetsMad",
"WW",
"WZ",
"T-s",
"T-t",
"T-tWDR",
"Tbar-s",
"Tbar_t",
"Tbar_tWDR",
"ZllHbb105",
"ZllHbb110",
"ZllHbb115",
"ZllHbb120",
"ZllHbb125",
"ZllHbb130",
"ZllHbb135",
"ZJets100Mad",
"Zinv100",
"TTbarv4"
]

for item in L:
        fib(item)

