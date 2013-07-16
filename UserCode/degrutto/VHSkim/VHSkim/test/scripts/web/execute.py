#!/usr/bin/env python
# file: hypotenuse.py

import sys, math
import os
import commands
import string


#if len(sys.argv) != 3:  # the program name and the two arguments
  # stop the program and print an error message
#  sys.exit("Must provide two positive numbers")

# Convert the two arguments from strings into numbers
#x = float(sys.argv[1])


var = sys.argv[1]
xtitle = sys.argv[2]
cut = sys.argv[3]
xlow = sys.argv[4]
xhigh = sys.argv[5]
binsize = sys.argv[6]
dolog = sys.argv[7]

f = open("num", "r")
num = int(f.read())
num = num + 1
os.system('rm num')
os.system('echo %i > num' % num)
num = num - 1

os.system("echo 'makePlots(\"%s\",\"%s\",\"%s\",\"%i.png\",0.01,%s,%s,%s,Option,%s,0,0,9999,1);' > /home/madfish/Rizzi/%i.C" % (var,xtitle,cut,num,binsize,xlow,xhigh,dolog,num))
os.system('root -l -b -q /home/mfisher/Rizzi/%i.C > %i.log' % (num,num))

