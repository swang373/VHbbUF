#!/usr/bin/env python

filename = "./GrowTree.C"

from ROOT import gROOT
if not "5.34" in gROOT.GetVersion():
    print "Please use ROOT 5.34:"
    print "source startup_6_1_1.csh"
else:
    flag = gROOT.ProcessLine('.L %s++O' % filename)
    import os
    print '-'*40
    print
    maxentries = 1500000
    processes = [



         "ZnnH125", 
           "WlnH125", 
          "monoH",
#        "WW", "WZ", "ZZ",
        "T_t", "Tbar_t", "T_s", "Tbar_s", "T_tW", "Tbar_tW",
 
       "WJetsIncl", "WJetsHT100", "WJetsHT200", "WJetsHT400",  "WJetsHT600", 
         "ZJetsHT100", "ZJetsHT200", "ZJetsHT400",  "ZJetsHT600", 
        "TTPythia8", 
        "QCDHT100", "QCDHT250", "QCDHT500", "QCDHT1000", 
        ]
    method = "BDTG"
    outfiles = []
    for p in processes:
        if p[-1].isdigit() and ":" in p:
            outfiles.append([])
            p, njobs = p.split(":")
            njobs = int(njobs)
            for j in xrange(njobs):
                begin = j*maxentries
                end = (j+1)*maxentries if (j != njobs-1) else -1
                print "root -b -l -q %s+\(\\\"%s\\\",\\\"%s\\\",%i,%i\)" % (filename,p,method,begin,end)
                outfiles[-1].append("Step3_%s_%i_%i.root" % (p, begin, end))
        else:
            print "root -b -l -q %s+\(\\\"%s\\\",\\\"%s\\\"\)" % (filename,p,method)
    print 
    if outfiles:
        print "cd skim"
        for f in outfiles:
            tup = f[0].split("_")
            assert(len(tup)==4)
            print "hadd %s.root %s" % ("_".join([tup[0],tup[1]]), " ".join(f))
        print "cd -"
        print 
