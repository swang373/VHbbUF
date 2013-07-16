#!/usr/bin/env python
import re
from sys import argv, stdout, stderr, exit
from optparse import OptionParser

# import ROOT with a fix to get batch mode (http://root.cern.ch/phpBB3/viewtopic.php?t=3198)
hasHelp = False
for X in ("-h", "-?", "--help"):
    if X in argv:
        hasHelp = True
        argv.remove(X)
argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
argv.remove( '-b-' )
if hasHelp: argv.append("-h")

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")

(options, args) = parser.parse_args()
if len(args) == 0:
    parser.print_usage()
    exit(1)

file = ROOT.TFile(args[0])
if file == None: raise RuntimeError, "Cannot open file %s" % args[0]
fit_s  = file.Get("fit_s")
fit_b  = file.Get("fit_b")
prefit = file.Get("nuisances_prefit")
if fit_s == None or fit_s.ClassName()   != "RooFitResult": raise RuntimeError, "File %s does not contain the output of the signal fit 'fit_s'"     % args[0]
if fit_b == None or fit_b.ClassName()   != "RooFitResult": raise RuntimeError, "File %s does not contain the output of the background fit 'fit_b'" % args[0]
if prefit == None or prefit.ClassName() != "RooArgSet":    raise RuntimeError, "File %s does not contain the prefit nuisances 'nuisances_prefit'"  % args[0]

isFlagged = {}
table = {}
fpf_b = fit_b.floatParsFinal()
fpf_s = fit_s.floatParsFinal()
pulls = []

for i in range(fpf_b.getSize()):
    nuis_b = fpf_b.at(i)
    name   = nuis_b.GetName();
    nuis_p = prefit.find(name)
    mean_p, sigma_p = (nuis_p.getVal(), nuis_p.getError())
    row = []
    
    m = re.match(r"CMS_vhbb_(.*)_SF.*", name)
    if m:
        row += [(nuis_b.getVal() - mean_p)/sigma_p, nuis_b.getError()/sigma_p]
        table[m.group(0)] = row
    m = re.match(r"CMS_vhbb_(WJSlope|ZJSlope)_Znunu_8TeV", name)
    if m:
        row += [(nuis_b.getVal() - mean_p)/sigma_p, nuis_b.getError()/sigma_p]
        table[m.group(0)] = row
    m = re.match(r"CMS_vhbb_(eff_b|fake_b_8TeV|res_j|scale_j)$", name)
    if m:
        row += [(nuis_b.getVal() - mean_p)/sigma_p, nuis_b.getError()/sigma_p]
        table[m.group(0)] = row

#print table

table2 = {}

for i in table.keys():
    for j in table.keys():
        row = []
        row += [fit_b.correlation(i,j)]
        table2[(i,j)] = row

#print table2

print "Scale factors"
orders = ["Wj0b", "Wj1b", "Wj2b", "WjLF", "WjHF", "Zj0b", "Zj1b", "Zj2b", "ZjLF", "ZjHF", "TT"]
#orders = ["ZH", "WH", "Wj0b", "Wj1b", "Wj2b", "WjLF", "WjHF", "Zj0b", "Zj1b", "Zj2b", "ZjLF", "ZjHF", "TT", "s_Top", "VV"]
results = []
for o in orders:
    for name in table.keys():
        m = re.match(r"CMS_vhbb_%s_SF.*" % o, name)
        if m:
            print "{0:<8}".format(o), "% .3f +/- %.3f" % tuple(table[name])
            results.append((o, name, table[name][0], table[name][1]))
print

print "Slopes"
results2 = []
for o in ["WJSlope", "ZJSlope"]:
    for name in table.keys():
        if (name == "CMS_vhbb_%s_Znunu_8TeV" % o):
            print "{0:<8}".format(o), "% .3f +/- %.3f" % tuple(table[name])
            results2.append((o, name, table[name][0], table[name][1]))
print

print "Shapes"
results3 = []
for o in ["eff_b", "fake_b_8TeV", "res_j", "scale_j"]:
    for name in table.keys():
        if (name == "CMS_vhbb_%s" % o):
            print "{0:<8}".format(o.replace("_8TeV","")), "% .3f +/- %.3f" % tuple(table[name])
            results3.append((o, name, table[name][0], table[name][1]))
print

print "Scale factor correlations"
for o, name, v, e in results:
    print "{0:<8}".format(o),
    for o2, name2, v2, e2 in results:
        print "% .3f  " % tuple(table2[name,name2]),
    print
print

print "Put scalefactors (1st row), scalefactors_lnN (2nd row)"
string = "".join(["% .3f, " % v for o, name, v, e in results])
print string
string = "".join(["% .3f, " % e for o, name, v, e in results])
print string
