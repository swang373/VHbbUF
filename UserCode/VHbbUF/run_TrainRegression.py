#!/usr/bin/env python

filename = "./TrainRegression.C"

from ROOT import gROOT
if not "5.34" in gROOT.GetVersion():
    print "Please use ROOT 5.34:"
    print "source startup_6_1_1.csh"
else:
    flag = gROOT.ProcessLine('.L %s++O' % filename)
    print '-'*40
    print
    print "root -b -l -q %s+\(\\\"BDTG\\\"\)" % (filename)
