#!/usr/bin/env python

filename = "./TrainClassification.C"

from ROOT import gROOT
flag = gROOT.ProcessLine('.L %s++O' % filename)

print '-'*40
print ''
print 'root -b -l -q %s+\(\\\"BDT\\\"\)' % (filename)

