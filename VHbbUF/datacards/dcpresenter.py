#!/usr/bin/env python

import os
import re

names = {
    1 : "high-pt",
    2 : "med-pt",
    3 : "low-pt",
    6 : "high+med",
    5 : "combo",
    }

for i in [1,2,3,6,5]:
    nA, na, nP, np, nM, nM1, nM2 = 0., 0., 0., 0., 0., 0., 0.
    
    with open("A%i.log" % i) as f:
        for line in f:
            m = re.match("Observed Limit: r < (.*)", line)
            if m:  nA = float(m.group(1))
            m = re.match("Expected 50.0%: r < (.*)", line)
            if m:  na = float(m.group(1))

    with open("P%i.log" % i) as f:
        for line in f:
            m = re.match("       \(Significance = (.*)\)", line)
            if m:  nP = float(m.group(1))
    
    with open("p%i.log" % i) as f:
        for line in f:
            m = re.match("       \(Significance = (.*)\)", line)
            if m:  np = float(m.group(1))
    
    with open("M%i.log" % i) as f:
        for line in f:
            m = re.match("Best fit r: (.*)  (.*)/(.*)  \(68% CL\)", line)
            if m:
                nM = float(m.group(1))
                nM1 = float(m.group(2))
                nM2 = float(m.group(3))

    if i == 1:
        print "          : limit          signif            mu"
        print "---------------------------------------------------------------"
    print "{0:10s}:  {1: .2f} ({2: .2f})   {3: .3f} ({4: .3f})   {5: .1f} [{6:+.1f} {7:+.1f}]".format(names[i], nA, na, nP, np, nM, nM1, nM2)
