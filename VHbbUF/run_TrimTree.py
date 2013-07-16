#!/usr/bin/env python

filename = "./TrimTree.C"

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
        "ZnnH110", "ZnnH115", "ZnnH120", "ZnnH125", "ZnnH130", "ZnnH135", "ZnnH140", "ZnnH145", "ZnnH150", 
        "WlnH110", "WlnH115", "WlnH120", "WlnH125", "WlnH130", "WlnH135", "WlnH140", "WlnH145", "WlnH150", 
        "WW", "WZ", "ZZ",
        "T_t", "Tbar_t", "T_s", "Tbar_s", "T_tW", "Tbar_tW",
        "WJetsPtW70", "WJetsPtW100", "WJetsPtW180", "WJetsHW", 
        "ZJetsHT50", "ZJetsHT100", "ZJetsHT200", "ZJetsHT400", "ZJetsHW", "ZJetsPtZ100",
        "TTPowheg", "TTMCatNLO", "TTFullLeptMG", "TTSemiLeptMG", "TTHadronicMG",
        "QCDPt50", "QCDPt80", "QCDPt120", "QCDPt170", "QCDPt300", "QCDPt470", "QCDPt600", "QCDPt800", "QCDPt1000", "QCDPt1400", "QCDPt1800",
        "Data_R", "Data_P",
        ]
    method = "BDT"
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
