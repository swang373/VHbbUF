import os
import re
from ROOT import gROOT, TTree, TFile
from array import array

### WARNING: make sure the loop over nEventsMin, MaxDepth, AdaBoostBeta 
###          are identical to the one in TrainBDT.C !!
###

# This specifies the set of NTrees parameters
r = [400,500,600,800]
#r = range(400,600,100)

w = "weights_V6FullB"
logs = "logs_V6FullB_20130329"

#w = "weights_HybDJ"
#logs = "logs_HybDJ_20130306"

#w = "weights_HybFJ"
#logs = "logs_HybFJ_20130306"

channel = "ZnunuHighPt"
loglimit = "a3"
logsignif = "b3"

#channel = "ZnunuMedPt"
#loglimit = "a4"
#logsignif = "b4"

#channel = "ZnunuLowPt"
#loglimit = "a5"
#logsignif = "b5"

#channel = "ZnunuLowCSV"
#loglimit = "a6"
#logsignif = "b6"

minlimit, minpvalue = 999, 999
maxsignif = -1
for i in r:
    store = {}
    ibdt = 0
    
    for j in xrange(5):
        nEventsMin_ = 150 + j*50
        for k in xrange(2):
            MaxDepth_ = 3 + k
            for l in xrange(1):
                AdaBoostBeta100_ = 10 + l*10
    
                limit = 100.
                with open("%s/%s_%i_%i" % (logs,loglimit,i,ibdt)) as f:
                    fr = f.read()
                    m = re.search("Expected 50.0\%: r < ([0-9.]+)", fr)
                    if m:
                        limit = float(m.group(1))
                        #print limit
                    else:
                        raise AttributeError("Limit not found!")
                pvalue = 100.
                signif = 0.
                with open("%s/%s_%i_%i" % (logs,logsignif,i,ibdt)) as f:
                    fr = f.read()
                    m = re.search("p-value of background: ([0-9.]+)", fr)
                    if m:
                        pvalue = float(m.group(1))
                        #print pvalue
                    else:
                        raise AttributeError("p-value not found!")
                    m = re.search("Significance = ([0-9.]+)", fr)
                    if m:
                        signif = float(m.group(1))
                        #print signif
                    else:
                        raise AttributeError("significance not found!")
                
                store[(i, nEventsMin_, MaxDepth_, AdaBoostBeta100_)] = (limit, pvalue, signif)
                ibdt += 1
    assert(ibdt<=10)
    
    f = TFile.Open("%sNTrees%i/125/TMVA_%s.root" % (w,i,channel))
    tree = f.Get("fomtree")
    nentries = tree.GetEntries()
    
    g = TFile.Open("%sNTrees%i/125/TMVA_%s_new.root" % (w,i,channel), "RECREATE")
    ctree = tree.CloneTree(0)

    limit  = array( 'f', [0.] )
    pvalue = array( 'f', [0.] )
    signif = array( 'f', [0.] )
    ctree.Branch( "COMB_limit", limit, "COMB_limit/F")
    ctree.Branch( "COMB_pvalue", pvalue, "COMB_pvalue/F")
    ctree.Branch( "COMB_signif", signif, "COMB_signif/F")

    
    for e in xrange(nentries):
        tree.GetEntry(e)
        tup = (tree.NTrees, tree.nEventsMin, tree.MaxDepth, tree.AdaBoostBeta100)
        limit[0] = store[tup][0]
        pvalue[0] = store[tup][1]
        signif[0] = store[tup][2]
        print tup, "%.3f  %.3f  %.3f " % store[tup], "%.3f  %.3f" % (tree.TMVA_kolS, tree.TMVA_kolB)
        ctree.Fill()
        
        if minlimit > limit[0]:
            minlimit = limit[0]
        if minpvalue > pvalue[0]:
            minpvalue = pvalue[0]
        if maxsignif < signif[0]:
            maxsignif = signif[0]
    ctree.Write()
    g.Close()
    f.Close()

print "%.3f, %.3f %.3f" % (minlimit, minpvalue, maxsignif)

print "hadd -f TMVA_%s_new.root %sNTrees*/125/TMVA_%s_new.root" % (channel,w,channel)

# fomtree->Scan("NTrees:nEventsMin:MaxDepth:COMB_limit:COMB_signif:TMVA_kolS:TMVA_kolB:USER_psig","COMB_limit<3.6 && COMB_signif>0.68 && NTrees<400")
