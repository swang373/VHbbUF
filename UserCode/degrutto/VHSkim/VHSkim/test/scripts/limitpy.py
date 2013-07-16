#! /usr/bin/env python


import sys
import os
import commands
import string

higgsNo  = "135"

def fib(var):
    os.system('echo "combine vhbb_DC_ALL_CC_M%s_Aug16.txt -M HybridNew --frequentist --testStat LHC -s 22%s --singlePoint %s --saveToys --saveHybridResult > %s-%s.txt"' % (higgsNo,var,var,higgsNo,var))
    os.system('combine vhbb_DC_ALL_CC_M%s_Aug16.txt -M HybridNew --frequentist --testStat LHC -s 22%s --singlePoint %s --saveToys --saveHybridResult > %s-%s.txt' % (higgsNo,var,var,higgsNo,var))



fib(2)
fib(3)
fib(4)
fib(5)
fib(7)
fib(9)
fib(11)
fib(15)
fib(17)
fib(21)
fib(25)
fib(30)
fib(35)
fib(37)
fib(39)

