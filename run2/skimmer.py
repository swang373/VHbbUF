#!/usr/bin/env python

from ROOT import TFile, TH1D, TTree, TTreeFormula, gROOT
from pyhelper import *
import os
import re

class Skimmer:
    """
    I skim the Step2 ntuples.
    """

    def __init__(self):
        self.reader = Reader()
        self.has_read = False
        #self.mandatory_settings = [
        #    "treename",
        #    "ftreename",
        #    "apply_bdt_regression",  # used by skim_for_classification
        #    ]

    def read(self, ini):
        self.reader.read(ini)
        if self.reader.has_read:
            self.has_read = True
        #for s in self.mandatory_settings:
        #    if not hasattr(self.reader, s):
        #        raise AttributeError("%s: Failed to find attribute: %s" % (self.__class__.__name__, s))
        return 0

    def copytree_with_friendtree(self, selection, ttree, friendttree):
        nentries = ttree.GetEntriesFast()
        assert(nentries == friendttree.GetEntriesFast())
        out_ttree = ttree.CloneTree(0)
        out_friendttree = friendttree.CloneTree(0)
        
        ttree.AddFriend(friendttree)
        ttreeformula = TTreeFormula("ttf", selection, ttree)
        if not ttreeformula.GetNdim():
            raise ValueError("%s: Failed to find any TTree variable from the selection: %s" % (self.__class__.__name__, selection))
        #warnings.filterwarnings( action='ignore', category=RuntimeWarning, message='creating converter for unknown type.*' )
        ttreeformula.SetQuickLoad(1)
        
        treenumber = -1
        for ientry in xrange(nentries):
            entrynumber = ttree.GetEntryNumber(ientry)
            friendentrynumber = friendttree.GetEntryNumber(ientry)
            if entrynumber < 0:  break
            localentry = ttree.LoadTree(entrynumber)
            if localentry < 0:  break
            if treenumber != ttree.GetTreeNumber():
                treenumber = ttree.GetTreeNumber()
                ttreeformula.UpdateFormulaLeaves()
            ndata = ttreeformula.GetNdata()
            keep = False
            for idata in xrange(ndata):
                keep |= (bool(ttreeformula.EvalInstance(idata)) != 0)
                if keep:  break
            if not keep:  continue

            ttree.GetEntry(entrynumber)
            out_ttree.Fill()
            friendttree.GetEntry(friendentrynumber)
            out_friendttree.Fill()
        return (out_ttree, out_friendttree)

    def splittree_with_friendtree(self, selection, ttree, friendttree, pass_tfile, fail_tfile):
        assert(pass_tfile != 0)
        assert(fail_tfile != 0)
        nentries = ttree.GetEntriesFast()
        assert(nentries == friendttree.GetEntriesFast())
        pass_tfile.cd()
        pass_ttree = ttree.CloneTree(0)
        pass_friendttree = friendttree.CloneTree(0)
        fail_tfile.cd()
        fail_ttree = ttree.CloneTree(0)
        fail_friendttree = friendttree.CloneTree(0)
        
        ttree.AddFriend(friendttree)
        ttreeformula = TTreeFormula("ttf", selection, ttree)
        if not ttreeformula.GetNdim():
            raise ValueError("%s: Failed to find any TTree variable from the selection: %s" % (self.__class__.__name__, selection))
        #warnings.filterwarnings( action='ignore', category=RuntimeWarning, message='creating converter for unknown type.*' )
        ttreeformula.SetQuickLoad(1)
        
        treenumber = -1
        for ientry in xrange(nentries):
            entrynumber = ttree.GetEntryNumber(ientry)
            friendentrynumber = friendttree.GetEntryNumber(ientry)
            if entrynumber < 0:  break
            localentry = ttree.LoadTree(entrynumber)
            if localentry < 0:  break
            if treenumber != ttree.GetTreeNumber():
                treenumber = ttree.GetTreeNumber()
                ttreeformula.UpdateFormulaLeaves()
            ndata = ttreeformula.GetNdata()
            keep = False
            for idata in xrange(ndata):
                keep |= (bool(ttreeformula.EvalInstance(idata)) != 0)
                if keep:  break

            if keep:
                ttree.GetEntry(entrynumber)
                pass_ttree.Fill()
                friendttree.GetEntry(friendentrynumber)
                pass_friendttree.Fill()
            else:
                ttree.GetEntry(entrynumber)
                fail_ttree.Fill()
                friendttree.GetEntry(friendentrynumber)
                fail_friendttree.Fill()                
        return (pass_ttree, pass_friendttree, fail_ttree, fail_friendttree)

    def get_friendtree(self, ttree, sample, isafriend=False):
        name, ext = os.path.splitext(sample)        
        if sample.startswith("dcache"):
            raise TypeError("%s: The sample cannot be in dCache: %s" % (self.__class__.__name__, sample))
        if not ext == ".root":
            raise ValueError("%s: Not a .root file: %s" % (self.__class__.__name__, sample))
        
        friendsample, treename = "", ""
        if isafriend:
            friendsample = name.replace("_friend","")+ext
            treename = self.reader.treename
        else:
            friendsample = name+"_friend"+ext
            treename = self.reader.ftreename
        if os.path.exists(friendsample):
            friendtfile = TFile.Open(friendsample, "READ")
            friendttree = friendtfile.Get(treename)
            if not isinstance(friendttree, TTree):
                raise TypeError("%s: Failed to get the friend tree! Got: %s" % (self.__class__.__name__, repr(friendttree)))
            return friendttree
        
        raise TypeError("%s: The sample has no friend: %s" % (self.__class__.__name__, sample))
        return TTree()

    def skim_for_classification(self, selection, indir, outdir, prefix=""):
        if not indir.endswith('/'):
            indir += '/'
        
        selections = self.reader.get_selections()
        if not selections:
            raise ValueError("%s: Found zero selection!" % (self.__class__.__name__))
        elif not selection in selections:
            raise ValueError("%s: Selection not found: \"%s\"" % (self.__class__.__name__, selection))

        # replace variables if apply_bdt_regression
        replacements = self.reader.get_regression_replacements()
        sel = selections[selection]
        for old, new in replacements.iteritems():
            sel = re.sub(r"(\A|\W)%s(\Z|\W)" % old, r"\1%s\2" % new, sel)
            #sel = sel.replace(old, new)

        for sample in os.listdir(indir):
            m = re.match(prefix+"(\S+).root", sample)
            if m:
                process = m.group(1)
                if "friend" in process:
                    continue
                if "Data" in process or "data" in process:
                    continue
                self.skim_sample(process, indir+sample, sel, outdir=outdir, prefix=prefix, addfriend=True)
        return 0

    def skim_sample(self, process, sample, selection, outdir, prefix="", addfriend=False, skimentries=1000000000):
        if not outdir.endswith('/'):
            outdir += '/'
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        tfile = TFile.Open(sample, "READ")
        treename = self.reader.treename
        ttree = tfile.Get(treename)
        if not isinstance(ttree, TTree):
            raise TypeError("%s: Failed to get the tree! Got: %s" % (self.__class__.__name__, repr(ttree)))
        
        if addfriend:
            friendttree = self.get_friendtree(ttree, sample)
            if skimentries != 1000000000:
                warnings.warn("%s: skimentries is not supported for addfriend=True" % (self.__class__.__name__))
        
        out = outdir + prefix + process + ".root"
        out_tfile = TFile.Open(out, "RECREATE")
        out_ttree_list = []
        if not addfriend:
            out_ttree = ttree.CopyTree(selection)
            if out_ttree.GetEntriesFast() > skimentries:
                out_ttree1 = out_ttree.CopyTree("", "", skimentries)
                out_ttree1.SetName("%s_1" % out_ttree1.GetName())
                out_ttree2 = out_ttree.CopyTree("", "", 1000000000, skimentries)
                out_ttree2.SetName("%s_2" % out_ttree2.GetName())
                out_ttree_list = [out_ttree1, out_ttree2]
            else:
                out_ttree_list = [out_ttree]
        else:
            (out_ttree, out_friendttree) = self.copytree_with_friendtree(selection, ttree, friendttree)
            out_ttree_list = [out_ttree, out_friendttree]
        
        entries = [ttree.GetEntriesFast(), out_ttree.GetEntriesFast()]
        assert(entries[0] > entries[1])
        
        if addfriend:
            friendentries = [friendttree.GetEntriesFast(), out_friendttree.GetEntriesFast()]    
            for i in xrange(2):
                assert(entries[i] == friendentries[i])
            assert(friendentries[0] > friendentries[1])

        if out_ttree_list:
            for o in out_ttree_list:
                o.Write()
        out_tfile.Close()
        
        print "%s: skimmed from %i to %i entries." % (process, entries[0], entries[1])
        tfile.Close()
        return 0

    def split_for_classification(self, indir="./", prefix=""):
        if not indir.endswith('/'):
            indir += '/'
        evendir = indir+"even/"
        odddir = indir+"odd/"
        if not os.path.exists(evendir):
            os.makedirs(evendir)
        if not os.path.exists(odddir):
            os.makedirs(odddir)
        
        #selection = "Entry$ %2"
        selection = "EVENT.event %2 == 0"
        for sample in os.listdir(indir):
            m = re.match(prefix+"(\S+).root", sample)
            if m:
                process = m.group(1)
                if "friend" in process:  continue
                
                tfile = TFile.Open(indir+sample, "READ")
                ttree = tfile.Get(self.reader.treename)
                if not isinstance(ttree, TTree):
                    raise TypeError("%s: Failed to get the tree! Got: %s" % (self.__class__.__name__, repr(ttree)))
                friendttree = tfile.Get(self.reader.ftreename)
                if not isinstance(friendttree, TTree):
                    raise TypeError("%s: Failed to get the friend tree! Got: %s" % (self.__class__.__name__, repr(friendttree)))

                even_tfile = TFile.Open(evendir+sample, "RECREATE")
                odd_tfile = TFile.Open(odddir+sample, "RECREATE")                
                (even_ttree, even_friendttree, odd_ttree, odd_friendttree) = self.splittree_with_friendtree(selection, ttree, friendttree, even_tfile, odd_tfile)
                entries = [ttree.GetEntriesFast(), even_ttree.GetEntriesFast(), odd_ttree.GetEntriesFast()]
                friendentries = [friendttree.GetEntriesFast(), even_friendttree.GetEntriesFast(), odd_friendttree.GetEntriesFast()]
                for i in xrange(3):
                    assert(entries[i] == friendentries[i])
                assert(entries[0] == entries[1] + entries[2])
                assert(friendentries[0] == friendentries[1] + friendentries[2])
                
                even_tfile.cd()
                even_ttree.Write()
                even_friendttree.Write()
                even_tfile.Close()
                
                odd_tfile.cd()
                odd_ttree.Write()
                odd_friendttree.Write()
                odd_tfile.Close()
                
                print "%s: %f = %f (even) + %f (odd)" % (process, entries[0], entries[1], entries[2])
                tfile.Close()
        return 0

    def process(self, indir, prefix="", processlist=[]):
        samples = {}
        for sample in os.listdir(indir):
            m = re.match(prefix+"(\S+).root", sample)
            if m:
                process = m.group(1)
                samples[process] = indir + sample

        if not processlist:
            processlist = self.reader.get_processes().keys()
            if not processlist:
                raise ValueError("%s: Found zero process!" % (self.__class__.__name__))
        
        for process in processlist:
            process_ = process
            m = re.match("(W|Z).+Jets(.*)", process)  # for (Z|W)(b|c|udsg)Jets, use (Z|W)Jets
            if m:
                process_ = m.group(1)+"Jets"+m.group(2)
            if not process_ in samples:
                warnings.warn("%s: Failed to find sample: \"%s\". Skipping..." % (self.__class__.__name__, process))
                continue
            
            tfile = TFile.Open(samples[process_], "READ")
            th1d = tfile.Get("Count")
            if not isinstance(th1d, TH1D):
                raise TypeError("%s: Failed to get the histogram! Got: %s" % (self.__class__.__name__, repr(th1d)))
            count = th1d.GetBinContent(1)
#            countwithpu = th1d.GetBinContent(2)
 #           if abs(count - countwithpu)/count > 0.2:
 #               warnings.warn("%s Count and CountWithPU disagree by more than 20%: %s != %s" % (self.__class__.__name__, count, countwithpu))
            
            value = self.reader.cfgparser.get("Process", process)
            self.reader.is_a_list(value, 4)
            values = [v.strip() for v in value.split(" , ")]
            
            print "{0:16s}= {1:18s} , {2:10.1f} , {3:13.4f} , {4:s}".format(process, values[0], count, count, values[3])
            tfile.Close()
        return 0

    def stitch(self, indir, prefix="", processlist=[]):
        samples = {}
        for sample in os.listdir(indir):
            m = re.match(prefix+"(\S+).root", sample)
            if m:
                process = m.group(1)
                samples[process] = indir + sample

        if not processlist:
            processlist = self.reader.get_processes().keys()
            if not processlist:
                raise ValueError("%s: Found zero process!" % (self.__class__.__name__))

        items = OrderedDict(self.reader.cfgparser.items("Stitch"))
        n = 0
        name = ""
        lhecut, lhenorm = [], []
        lhecounts = {}
        for (key, value) in items.iteritems():
            if n == 0:
                n = len(value.split(" , "))
                name = key[0:key.find("LHECUT")]
            elif not key.startswith(name):
                n = len(value.split(" , "))
                name = key[0:key.find("LHECUT")]
            self.reader.is_a_list(value, n)
            
            if "LHECUT" in key:
                lhecut = eval(value)
                lhecounts.clear()
            
            elif key in processlist:
                tfile = TFile.Open(samples[key], "READ")
                #ttree = tfile.Get("tree")
                #if not isinstance(ttree, TTree):
                #    raise TypeError("%s: Failed to get the tree! Got: %s" % (self.__class__.__name__, repr(ttree)))
                #lhecounts[key] = []
                #for i in xrange(n):
                #    count = ttree.GetEntries(lhecut[i])
                #    lhecounts[key].append(count)
                
                th1d = tfile.Get("LHECount")
                if not isinstance(th1d, TH1D):
                    raise TypeError("%s: Failed to get the histogram! Got: %s" % (self.__class__.__name__, repr(th1d)))
                lhecounts[key] = []
                for i in xrange(n):
                    count = th1d.GetBinContent(i+1)
                    lhecounts[key].append(count)
                
                items[key] = ""
                for i in xrange(n):
                    items[key] += "{0:10.1f}".format(lhecounts[key][i])
                    if i != n-1:  items[key] += " , "
        
            elif "LHENORM" in key:
                lhenorm = eval(value)
                for p in lhenorm:
                    assert(p in lhecounts)
            
            elif "LHECOUNT" in key:
                items[key] = ""
                for i in xrange(n):
                    totallhecounts = 0.
                    for p in lhecounts.keys():
                        totallhecounts += lhecounts[p][i]
                    items[key] += "{0:.1f} / {1:.1f}".format(lhecounts[lhenorm[i]][i], int(totallhecounts))
                    if i != n-1:  items[key] += " , "
        
        for (key, value) in items.iteritems():
            print "{0:16s}= {1}".format(key, value)
            
        return 0


## Main
## ----------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    
    indir = "/afs/cern.ch/work/d/degrutto/public/MiniAOD/ZnnHbb_Phys14_PU20bx25/skim/"
    outdir = "skim/"
    prefix = "skim_"


    argc = 1
    if len(sys.argv) != argc+1:
        print "ERROR: Expect exactly %i arguments: %s" % (argc, ' '.join(sys.argv[1:]))
        sys.exit(1)

    gROOT.LoadMacro("HelperFunctions.h")
    gROOT.SetBatch(True)
    
    # Call the class
    skimmer = Skimmer()
    skimmer.read(sys.argv[1])
    
    
    # Skim BDT classification selection
    #selection = "ZnunuHbb"
    #outdir = "skim_ZnunuHbb_classification/"
    #selection = "ZnunuHbbLowMET"
    #outdir = "skim_ZnunuHbb_classification_lowpt/"
    #selection = "ZnunuHbbLooseCSV"
    #outdir = "skim_ZnunuHbb_classification_loosecsv/"
    #skimmer.skim_for_classification(selection, indir=indir, outdir=outdir, prefix=prefix)
    
    # Split BDT classification even vs. odd
    #indir = "skim_ZnunuHbb_classification/"
    #indir = "skim_ZnunuHbb_classification_lowpt/"
    #indir = "skim_ZnunuHbb_classification_loosecsv/"
    #skimmer.split_for_classification(indir=indir, prefix=prefix)
    
    # Prepare counts, countswithpu in the .ini file
    #skimmer.process(indir=indir, prefix=prefix)
    #skimmer.process(indir=indir, prefix=prefix, processlist=["ZnnH125"])
    
    # Prepare LHE lhecut, lhenorm, lhecount in the .ini file
    skimmer.stitch(indir=indir, prefix=prefix)
