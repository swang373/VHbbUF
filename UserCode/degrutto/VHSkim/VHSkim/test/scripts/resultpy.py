#! /usr/bin/env python


import sys
import os
import commands
import string

higgsNo  = "110"

def fib(var):
    os.system('echo "combine vhbb_DC_ALL_CC_M%s_Aug16.txt -M HybridNew --frequentist --testStat LHC  --grid=H%s.root --expectedFromGrid %s >> %s-Results.txt"' % (higgsNo,higgsNo,var,higgsNo))
    os.system('echo " " >> %s-Results.txt' % higgsNo)
    os.system('echo " " >> %s-Results.txt' % higgsNo)
    os.system('echo "//////////////////////////////" >> %s-Results.txt' % higgsNo)
    os.system('echo "EXPECTED FROM GRID = %s" >> %s-Results.txt' % (var,higgsNo)) 
    os.system('echo "//////////////////////////////" >> %s-Results.txt' % higgsNo)
    os.system('echo " " >> %s-Results.txt' % higgsNo)
    os.system('echo " " >> %s-Results.txt' % higgsNo) 
    os.system('combine vhbb_DC_ALL_CC_M%s_Aug16.txt -M HybridNew --frequentist --testStat LHC  --grid=H%s.root --expectedFromGrid %s >> %s-Results.txt' % (higgsNo,higgsNo,var,higgsNo))



fib("0.5")
fib("0.16")
fib("0.84")
fib("0.0275")
fib("0.975")


os.system('echo " " >> %s-Results.txt' % higgsNo)
os.system('echo " " >> %s-Results.txt' % higgsNo)  
os.system('echo "///////////////////////////" >> %s-Results.txt' % higgsNo)
os.system('echo "EXPECTED FROM GRID = NONE" >> %s-Results.txt' % higgsNo) 
os.system('echo "///////////////////////////" >> %s-Results.txt' % higgsNo)
os.system('echo " " >> %s-Results.txt' % higgsNo)
os.system('echo " " >> %s-Results.txt' % higgsNo)  
os.system('combine vhbb_DC_ALL_CC_M%s_Aug16.txt -M HybridNew --frequentist --testStat LHC  --grid=H%s.root >> %s-Results.txt' % (higgsNo,higgsNo,higgsNo))



