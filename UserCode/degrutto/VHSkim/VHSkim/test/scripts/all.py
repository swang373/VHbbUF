#! /usr/bin/env python


import sys
import os
import commands
import string

def fib(n):
	os.system("sed -i '/here/s/CCC/%s/g' CutProcessor.C micro_Zll.C plotBDT.C" % n)
	os.system('root -l -b -q CutProcessor.C')
	signal = open('/home/madfish/root/tmva/test/signal', 'r');
	background = open('/home/madfish/root/tmva/test/background', 'r');
	snum = int( signal.read() );
	bnum = int( background.read() );
	os.system("sed -i '/here/s/NSS/%i/g' Class_ZH.C" % snum)
	os.system("sed -i '/heres/NBB/%i/g' Class_ZH.C" % bnum)
	os.system('./Class')
	os.system('./Application')
	os.system("root -l -b -q micro_Zll.C")
	os.system("root -l -b -q plotBDT.C")
	os.system("sed -i '/here/s/%s/CCC/g' CutProcessor.C micro_Zll.C plotBDT.C" % n)
	os.system("sed -i '/here/s/%i/NSS/g' Class_ZH.C" % snum)
	os.system("sed -i '/here/s/%i/NBB/g' Class_ZH.C" % bnum)
	os.system("scp %sbb*.png BDT-%s*.png mfisher@alachua.ihepa.ufl.edu:~/public_html/twiki" % (n,n))
	




L = [
#"ttbar",
#"zjl",
#"zjh"
"relax4"
]

for item in L:
        fib(item)

