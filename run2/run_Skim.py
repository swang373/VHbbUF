#!/usr/bin/env python
#filename = "./SkimTree.C"
# before run  python maketier2list.py -i 'arizzi/VHBBHeppyV11' -s pisa
import os 

from ROOT import gROOT
file = open('myfile.dat', 'r')
for f in file.readlines():
#    filename = 'root://xrootd.unl.edu//store/user/' + f.rstrip()
    filename = 'root://eoscms//eos/cms/store/user/'  + f.rstrip()
   #    print filename
    #    outname = f.rstrip() 
    #    print outname
    #    last = outname.split()[-1]
    last = f.rstrip().split("/")[-1]
    print 'last', last
    os.system("root -b -l -q Skim.C++\(\\\"%s\\\",\\\"%s\\\"\)" % (filename,f.rstrip().split('/')[-1] ))
#    os.system("root -b -l -q Skim.C++\(\\\"%s\\\"\)" % (filename) )
              #   root -b -l -q %s+\(\\\"%s\\\",\\\"%s\\\",%i,%i\)
#                print "root -b -l -q %s+\(\\\"%s\\\",\\\"%s\\\",%i,%i\)" % (filename,p,method,begin,end)
              
