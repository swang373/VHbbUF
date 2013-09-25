#!/usr/bin/env python

import os
import re
import sys
import warnings
from math import sqrt

# import ROOT without screwing up argv (http://mike-is-a-geek.blogspot.com/2011/03/python-pyroot-optparse.html)
tmpargv = sys.argv[:]  # [:] for a copy, not reference
sys.argv = []
from ROOT import TFile, TH1, TH1F, TColor, TCanvas, TPad, TLegend, THStack, TLatex, TPaveText, TLine, RooFit, RooAbsData, gROOT, gPad
sys.argv = tmpargv

gROOT.ProcessLine(".L tdrstyle.C")
gROOT.ProcessLine("setTDRStyle()")
gROOT.ProcessLine(".L HelperBDTShape.h")
from ROOT import setHisto

# import Higgs combination pythons
from HiggsAnalysis.CombinedLimit.DatacardParser import *
from HiggsAnalysis.CombinedLimit.ShapeTools     import *


class NuisanceBox:
    """A container of channel, process, syst, up/down histograms, valshift, sigshift, correlation"""
    
    def __init__(self, channel, process, syst, hup, hdown, valshift, sigshift, correlation):
        self.channel = channel
        self.process = process
        self.syst = syst
        if not isinstance(hup, TH1F):
            raise ValueError("%s: hup must be a TH1F" % (self.__class__.__name__))
        self.hup = hup
        if not isinstance(hdown, TH1F):
            raise ValueError("%s: hdown must be a TH1F" % (self.__class__.__name__))
        self.hdown = hdown
        if not isinstance(valshift, float):
            raise ValueError("%s: valshift must be a float" % (self.__class__.__name__))
        self.valshift = valshift # Delta_x / sigma_in
        if not isinstance(sigshift, float):
            raise ValueError("%s: sigshift must be a float" % (self.__class__.__name__))
        self.sigshift = sigshift # sigma_out / sigma_in
        if not isinstance(correlation, float):
            raise ValueError("%s: correlation must be a float" % (self.__class__.__name__))
        self.correlation = correlation

    def __repr__(self):
        return self.get().__repr__()
    
    def get(self):
        return (self.channel, self.process, self.syst, self.hup, self.hdown, self.valshift, self.sigshift, self.correlation)


def read_datacard(nuisances, DC, MB, options):
    """Read nuisances from datacards"""
    
    TH1.AddDirectory(0)  # don't put TH1Fs in the TFile directory
    TH1.SetDefaultSumw2(1)
    
    for b in DC.bins:
        print "*** ", b , " ***"
        if (options.channel != None) and (options.channel != b): continue  # b is the channel

        # Taken from ShapeTools.py
        bentry = None
        if DC.shapeMap.has_key(b): bentry = DC.shapeMap[b]
        elif DC.shapeMap.has_key("*"): bentry = DC.shapeMap["*"]
        else: raise KeyError, "Shape map has no entry for channel '%s'" % (b)
        
        # Taken from datacardDump.py
        exps = {}  # not used
        for (p,e) in DC.exp[b].items(): # so that we get only self.DC.processes contributing to this bin
            exps[p] = [ e, [] ]
        for (lsyst,nofloat,pdf,pdfargs,errline) in DC.systs:
            if pdf in ('param', 'flatParam'): continue
            # begin skip systematics
            skipme = False
            for xs in options.excludeSyst:
                if re.search(xs, lsyst): 
                    skipme = True
            if skipme: continue
            # end skip systematics
            
            varname = ""
            
            for p in DC.exp[b].keys(): # so that we get only self.DC.processes contributing to this bin
                if errline[b][p] == 0: continue
                # begin skip processes
                skipme = False
                for xp in options.processVetos:
                    if re.search(xp, p):
                        skipme = True
                if skipme: continue
                # end skip processes
                
                # Taken from ShapeTools.py
                #names = []
                #if bentry.has_key(p): names = bentry[p]
                #elif bentry.has_key("*"): names = bentry["*"]
                #elif DC.shapeMap["*"].has_key(p): names = DC.shapeMap["*"][p]
                #elif DC.shapeMap["*"].has_key("*"): names = DC.shapeMap["*"]["*"]
                #else: raise KeyError, "Shape map has no entry for process '%s', channel '%s'" % (p,b)
                #if len(names) == 1 and names[0] == "FAKE": raise ValueError, "FAKE systematic is not allowed" % (p,b)
                #name = names[1].split(":")[0]
                #varname = options.varnames[name]
                varname = options.varnames[b]
                if options.doMJJ:
                    varname = varname.replace("BDT_", "MJJ_")
                
                s0 = MB.getShape(b,p)
                assert(s0.InheritsFrom("RooDataHist"))
                n0 = "%s_%s" % (b,p)
                h0 = s0.createHistogram(varname)
                h0.SetName(n0)
                
                nui = NuisanceBox(b, p, "", h0, h0, 0., 0., 0.)
                if n0 not in nuisances:
                    nuisances[n0] = nui
                
                if pdf == 'gmN':
                    exps[p][1].append(1/sqrt(pdfargs[0]+1))
                elif pdf == 'gmM':
                    exps[p][1].append(errline[b][p])
                elif type(errline[b][p]) == list: 
                    kmax = max(errline[b][p][0], errline[b][p][1], 1.0/errline[b][p][0], 1.0/errline[b][p][1])
                    exps[p][1].append(kmax-1.)
                elif pdf == 'lnN':
                    exps[p][1].append(max(errline[b][p], 1.0/errline[b][p])-1.)
                    
                    lnN = errline[b][p]
                    nup = n0+'_'+lsyst+'Up'
                    hup = h0.Clone(nup)
                    hup.Scale(lnN)
                    
                    ndown = n0+'_'+lsyst+'Down'
                    hdown = h0.Clone(ndown)
                    hdown.Scale(1.0 - (lnN-1.0))
                    
                    nui = NuisanceBox(b, p, lsyst, hup, hdown, 0., 0., 0.)
                    n1 = "%s_%s_%s" % (b,p,lsyst)
                    if n1 in nuisances:
                        raise ValueError("Duplicate: %s" % n1)
                    else:
                        nuisances[n1] = nui
                    
                elif ("shape" in pdf) and not options.norm:
                    #s0 = MB.getShape(b,p)
                    sUp = MB.getShape(b,p,lsyst+"Up")
                    sDown = MB.getShape(b,p,lsyst+"Down")
                    
                    if (s0.InheritsFrom("TH1")):
                        ratios = [sUp.Integral()/s0.Integral(), sDown.Integral()/s0.Integral()]
                        ratios += [1/ratios[0], 1/ratios[1]]
                        exps[p][1].append(max(ratios) - 1)
                    elif (s0.InheritsFrom("RooDataHist")):
                        
                        name = "%s_%s" % (b,p)
                        
                        nup = n0+'_'+lsyst+'Up'
                        hup = sUp.createHistogram(varname)
                        hup.SetName(nup)
                        
                        nup = n0+'_'+lsyst+'Down'
                        hdown = sDown.createHistogram(varname)
                        hdown.SetName(ndown)
                        
                        nui = NuisanceBox(b, p, lsyst, hup, hdown, 0., 0., 0.)
                        n1 = "%s_%s_%s" % (b,p,lsyst)
                        if n1 in nuisances:
                            raise ValueError("Duplicate: %s" % n1)
                        else:
                            nuisances[n1] = nui
                        
        # Real data
        s0 = MB.getShape(b,options.dataname)
        assert(s0.InheritsFrom("RooDataHist"))
        n0 = "%s_%s" % (b,options.dataname)
        h0 = s0.createHistogram(varname)
        h0.SetName(n0)
        # Reset the errors
        nbins = h0.GetNbinsX()
        for ibin in xrange(nbins+2):
            bincontent = h0.GetBinContent(ibin)
            h0.SetBinError(ibin, sqrt(bincontent))
        nui = NuisanceBox(b, options.dataname, "", h0, h0, 0., 0., 0.)
        if n0 in nuisances:
            raise ValueError("Duplicate: %s" % n0)
        else:
            nuisances[n0] = nui
    return 0


def read_mlfit(nuisances, options):
    """Read post-fit nuisances from mlfit.root"""
    
    # Taken from diffNuisances.py
    tfile = TFile.Open(options.infile, "READ")
    if tfile == None: raise RuntimeError, "Cannot open file %s" % rootfile
    fit_s  = tfile.Get("fit_s")
    fit_b  = tfile.Get("fit_b")
    prefit = tfile.Get("nuisances_prefit")
    if fit_s == None or fit_s.ClassName()   != "RooFitResult": raise RuntimeError, "File %s does not contain the output of the signal fit 'fit_s'"     % rootfile
    if fit_b == None or fit_b.ClassName()   != "RooFitResult": raise RuntimeError, "File %s does not contain the output of the background fit 'fit_b'" % rootfile
    if prefit == None or prefit.ClassName() != "RooArgSet":    raise RuntimeError, "File %s does not contain the prefit nuisances 'nuisances_prefit'"  % rootfile

    isFlagged = {}
    table = {}
    fpf_b = fit_b.floatParsFinal()
    fpf_s = fit_s.floatParsFinal()
    #pulls = []
    pulls = {}
    for i in range(fpf_s.getSize()):
        nuis_s = fpf_s.at(i)
        name = nuis_s.GetName()
        nuis_b = fpf_b.find(name)
        nuis_p = prefit.find(name)
        row = []
        flag = False
        mean_p, sigma_p = 0,0
        rho = fit_s.correlation(name, options.poi)
        if nuis_p != None:
            mean_p, sigma_p = (nuis_p.getVal(), nuis_p.getError())

        fit_name, nuis_x = 'b', nuis_b
        if options.fit_s:
            fit_name, nuis_x = 's',nuis_s            
        if nuis_x != None and nuis_p != None:
            valshift = (nuis_x.getVal() - mean_p)/sigma_p
            sigshift = nuis_x.getError()/sigma_p
            pulls[name] = (valshift, sigshift, rho)
    
    for n, nui in nuisances.iteritems():
        if nui.syst == "":  continue
        pull = pulls[nui.syst]
        nui.valshift = pull[0]
        nui.sigshift = pull[1]
        nui.rho = pull[2]
    return 0


# Order to add histogram to THStack
orders_b = ["ZmmHighPt_8TeV", "ZmmLowPt_8TeV", "ZeeHighPt_8TeV", "ZeeLowPt_8TeV", "WmnHighPt_8TeV", "WmnLowPt_8TeV", "WenHighPt_8TeV", "WenLowPt_8TeV", "ZnunuHighPt_8TeV", "ZnunuMedPt_8TeV", "ZnunuLowPt_8TeV"]
orders_p = ["ZH", "WH", "ZH_SM", "WH_SM", "ZbbHinv", "Wj0b", "Wj1b", "Wj2b", "WjLF", "WjHF", "Zj0b", "Zj1b", "Zj2b", "ZjLF", "ZjHF", "TT", "s_Top", "VVLF", "VVHF", "WZHF", "ZZHF", "WZ", "ZZ", "QCD"]

def get_postfit(nuisances, options):
    TH1.SetDefaultSumw2(1)
    g_upcol   = TColor.GetColor("#FF6600")
    g_downcol = TColor.GetColor("#0033FF")
    g_postcol = TColor.GetColor("#33FF33")
    
    c1 = TCanvas("c1", "c1", 700, 700)
    c1.SetLogy();
    c1.Print(options.pdfname+"[");
    
    histograms = []
    for b in orders_b:
        for p in orders_p:
            key = "%s_%s" % (b,p)
            if not key in nuisances:  continue
            h0 = nuisances[key].hup
            h1 = h0.Clone("%s_%s_postfit" % (b,p))
            h2 = h0.Clone("%s_%s_errpostfit" % (b,p))
            hh1 = h0.Clone("%s_%s_prefit" % (b,p))
            hh2 = h0.Clone("%s_%s_errprefit" % (b,p))
            h2.Reset()
            hh2.Reset()
            nbins = h0.GetNbinsX()
            
            for n, nui in nuisances.iteritems():
                if nui.channel == b and nui.process == p and nui.syst != "":
                    s0 = nui.hup if nui.valshift > 0 else nui.hdown
                    assert(nbins == s0.GetNbinsX())
                    for ibin in xrange(nbins+2):
                        binmean = h0.GetBinContent(ibin)
                        binsyst = s0.GetBinContent(ibin)
                        binsigma = 0.
                        if (binmean > 1e-9):
                            binsigma = (binsyst/binmean - 1.0)
                        
                        h1.SetBinContent(ibin, h1.GetBinContent(ibin) + abs(nui.valshift) * binsigma * binmean)
                        h2.SetBinContent(ibin, h2.GetBinContent(ibin) + pow(nui.sigshift * binsigma * binmean, 2))
                        hh1.SetBinContent(ibin, hh1.GetBinContent(ibin) + 0.0 * binsigma * binmean)
                        hh2.SetBinContent(ibin, hh2.GetBinContent(ibin) + pow(1.0 * binsigma * binmean, 2))
            
            for ibin in xrange(nbins+2):
                h2.SetBinContent(ibin, sqrt(h2.GetBinContent(ibin)))
                hh2.SetBinContent(ibin, sqrt(hh2.GetBinContent(ibin)))
                
                if h1.GetBinContent(ibin) < 0:
                    warnings.warn("negative bin content! h=%s, binc=%f" % (h1.GetName(), h1.GetBinContent(ibin)))
                    h1.SetBinContent(ibin, 0)
                #if h2.GetBinContent(ibin) > 1:
                #    raise ValueError("postfit error larger than prefit! h=%s, binc=%f" % (h2.GetName(), h2.GetBinContent(ibin)))
            
            histograms.append(h1)
            histograms.append(h2)
            histograms.append(hh1)
            histograms.append(hh2)
            
            for n, nui in nuisances.iteritems():
                if nui.channel == b and nui.process == p and nui.syst != "":
                    
                    h0.SetStats(0)
                    h0.SetTitle("; BDT")
                    h0.SetLineColor(1)
                    h0.SetLineWidth(2)
                    h0.SetFillColor(0)
                    #h0.SetMarkerStyle(20)
                    h0.SetMinimum(0.01)
                    h0.GetXaxis().CenterTitle()
                    
                    nui.hup.SetLineColor(g_upcol)
                    nui.hup.SetLineWidth(2)
                    nui.hup.SetFillColor(0)
                    
                    nui.hdown.SetLineColor(g_downcol)
                    nui.hdown.SetLineWidth(2)
                    nui.hdown.SetFillColor(0)
                    
                    h1.SetMarkerStyle(20)
                    h1.SetMarkerColor(g_postcol)
                    h1.SetLineColor(g_postcol)
                    h1.SetLineWidth(2)
                    h1.SetFillColor(0)
                    
                    h0.Draw("hist")
                    nui.hup.Draw("hist same")
                    nui.hdown.Draw("hist same")
                    h0.Draw("hist same")
                    h1.Draw("e1 same")
                    
                    leg = TLegend(0.35, 0.20, 0.92, 0.35)
                    leg.SetFillColor(0)
                    leg.SetFillStyle(0)
                    leg.SetLineColor(0)
                    leg.SetShadowColor(0)
                    leg.SetTextFont(62)
                    #leg.SetTextSize(0.015)
                    leg.AddEntry(h0, n, "l")
                    leg.AddEntry(nui.hup, nui.syst+"Up", "l")
                    leg.AddEntry(nui.hdown, nui.syst+"Down", "l")
                    leg.AddEntry(h1, "postfit", "p")
                    leg.Draw()

                    gPad.RedrawAxis()
                    gPad.Modified()
                    gPad.Update()
                    gPad.Print(options.pdfname)
                    
                    h1.SetLineWidth(1)
                    h1.SetLineColor(1)
        
        key = "%s_%s" % (b,options.dataname)
        if not key in nuisances:  continue
        histograms.append(nuisances[key].hup)
    
    c1.Print(options.pdfname+"]")
    return histograms


def put_histograms(histograms, options):
    tfile = TFile.Open(options.outname, "RECREATE")
    for h in histograms:
        h.Write()
    tfile.Close()
    return 0


def MakePlots(varname, whatfit, options):
    if whatfit not in ["prefit", "postfit"]:
        raise ValueError("whatfit must be either 'prefit' or 'postfit': %s" % whatfit)
    
    massH    = 125
    plotData = True
    plotLog  = True
    plotSig  = True
    plotdir  = "./"
    kRed     = 632
    kYellow  = 400
    kBlue    = 600
    kGray    = 920
    
    if options.doMJJ:
        plotLog = False
    
    tfile = TFile.Open(options.outname, "READ")
    hZH        = tfile.Get(varname+"_"+"ZH_SM"+"_%s" % whatfit)
    hWH        = tfile.Get(varname+"_"+"WH_SM"+"_%s" % whatfit)
    hZbbHinv   = tfile.Get(varname+"_"+"ZH"+"_%s" % whatfit)
    hWj0b      = tfile.Get(varname+"_"+"Wj0b"+"_%s" % whatfit)
    hWj1b      = tfile.Get(varname+"_"+"Wj1b"+"_%s" % whatfit)
    hWj2b      = tfile.Get(varname+"_"+"Wj2b"+"_%s" % whatfit)
    hZj0b      = tfile.Get(varname+"_"+"Zj0b"+"_%s" % whatfit)
    hZj1b      = tfile.Get(varname+"_"+"Zj1b"+"_%s" % whatfit)
    hZj2b      = tfile.Get(varname+"_"+"Zj2b"+"_%s" % whatfit)
    hTT        = tfile.Get(varname+"_"+"TT"+"_%s" % whatfit)
    hs_Top     = tfile.Get(varname+"_"+"s_Top"+"_%s" % whatfit)
    hVVLF      = tfile.Get(varname+"_"+"VVLF"+"_%s" % whatfit)
    hWZHF      = tfile.Get(varname+"_"+"WZ"+"_%s" % whatfit)
    hZZHF      = tfile.Get(varname+"_"+"ZZ"+"_%s" % whatfit)
    hQCD       = tfile.Get(varname+"_"+"QCD"+"_%s" % whatfit)
    hdata_obs  = tfile.Get(varname+"_"+"data_obs")
    
    eZH        = tfile.Get(varname+"_"+"ZH_SM"+"_err%s" % whatfit)
    eWH        = tfile.Get(varname+"_"+"WH_SM"+"_err%s" % whatfit)
    eZbbHinv   = tfile.Get(varname+"_"+"ZbbHinv"+"_err%s" % whatfit)
    eWj0b      = tfile.Get(varname+"_"+"Wj0b"+"_err%s" % whatfit)
    eWj1b      = tfile.Get(varname+"_"+"Wj1b"+"_err%s" % whatfit)
    eWj2b      = tfile.Get(varname+"_"+"Wj2b"+"_err%s" % whatfit)
    eZj0b      = tfile.Get(varname+"_"+"Zj0b"+"_err%s" % whatfit)
    eZj1b      = tfile.Get(varname+"_"+"Zj1b"+"_err%s" % whatfit)
    eZj2b      = tfile.Get(varname+"_"+"Zj2b"+"_err%s" % whatfit)
    eTT        = tfile.Get(varname+"_"+"TT"+"_err%s" % whatfit)
    es_Top     = tfile.Get(varname+"_"+"s_Top"+"_err%s" % whatfit)
    eVVLF      = tfile.Get(varname+"_"+"VVLF"+"_err%s" % whatfit)
    eWZHF      = tfile.Get(varname+"_"+"WZ"+"_err%s" % whatfit)
    eZZHF      = tfile.Get(varname+"_"+"ZZ"+"_err%s" % whatfit)
    eQCD       = tfile.Get(varname+"_"+"QCD"+"_err%s" % whatfit)
    
    varname = hdata_obs.GetXaxis().GetTitle()
    xtitle = "BDT"
    if options.doMJJ:  xtitle = "m_{T} [GeV]"
    nbins = hdata_obs.GetNbinsX()
    xlow = hdata_obs.GetXaxis().GetXmin()
    xup = hdata_obs.GetXaxis().GetXmax()
    hVH        = TH1F("VH"       , "", nbins, xlow, xup)
    hVVHF      = TH1F("VVHF"     , "", nbins, xlow, xup)
    hVV        = TH1F("VV"       , "", nbins, xlow, xup)
    hmc_exp    = TH1F("mc_exp"   , "", nbins, xlow, xup)
    
    eVH        = TH1F("VH"+"_err"   , "", nbins, xlow, xup)
    eVVHF      = TH1F("VVHF"+"_err" , "", nbins, xlow, xup)
    eVV        = TH1F("VV"+"_err"   , "", nbins, xlow, xup)
    
    # Make sum of histograms
    hVH.Add(hZH)
    hVH.Add(hWH)
    eVH.Add(eZH)
    eVH.Add(eWH)

    hVVHF.Add(hWZHF)
    hVVHF.Add(hZZHF)
    eVVHF.Add(eWZHF)
    eVVHF.Add(eZZHF)

    hVV.Add(hVVLF)
    hVV.Add(hVVHF)
    eVV.Add(eVVLF)
    eVV.Add(eVVHF)

    hmc_exp.Add(hWj0b)
    hmc_exp.Add(hWj1b)
    hmc_exp.Add(hWj2b)
    hmc_exp.Add(hZj0b)
    hmc_exp.Add(hZj1b)
    hmc_exp.Add(hZj2b)
    hmc_exp.Add(hTT)
    hmc_exp.Add(hs_Top)
    if not options.doVV:
        hmc_exp.Add(hVVLF)
        hmc_exp.Add(hVVHF)
        hmc_exp.Add(hVH)  # VH is counted as background
    else:
        hmc_exp.Add(hVVLF)
        hmc_exp.Add(hVH)  # VH is counted as background
    hmc_exp.Add(hQCD)
    
    if True:
        print "MakePlots(): Setting up histograms..."

        # Setup canvas and pads
        c1 = TCanvas("c1", "c1", 700, 700)
        #c1 = TCanvas("c1", "c1", 600, 600)

        pad1 = TPad("pad1", "top pad"   , 0.0, 0.3, 1.0, 1.0)
        pad1.SetBottomMargin(0.0)
        pad1.Draw()
        pad2 = TPad("pad2", "bottom pad", 0.0, 0.0, 1.0, 0.3)
        pad2.SetTopMargin(0.0)
        pad2.SetBottomMargin(0.35)
        pad2.Draw()
        pad1.cd()
        pad1.SetLogy(plotLog)

        # Setup histogram styles and stack the histograms.
        setHisto(hVH, "VH")
        setHisto(hZH, "VH")
        setHisto(hWH, "VH")
        setHisto(hZbbHinv, "VH_1")
        setHisto(hWj0b, "WjLF")
        setHisto(hWj1b, "WjHFc")
        setHisto(hWj2b, "WjHFb")
        setHisto(hZj0b, "ZjLF")
        setHisto(hZj1b, "ZjHFc")
        setHisto(hZj2b, "ZjHFb")
        setHisto(hTT, "TT")
        setHisto(hs_Top, "s_Top")
        setHisto(hVV, "VV")
        setHisto(hVVLF, "VV")
        setHisto(hVVHF, "VVHF")
        setHisto(hQCD, "QCD")
        
        hdata_test = hdata_obs.Clone("hdata_test")  # blinded plot
        hmc_test = hmc_exp.Clone("hmc_test")  # for chi2 and KS test
        hdata_test.Sumw2()
        hmc_test.Sumw2()
        nbins_plot = hdata_test.GetNbinsX()
        assert(nbins_plot == hmc_test.GetNbinsX())
        if not plotData:  # be blind to the most sensitive bins
            if not options.doMJJ:
                for i in xrange(max(int(nbins_plot*0.75), nbins_plot-5),nbins_plot+2):
                    hdata_test.SetBinContent(i, 0.)
                    hdata_test.SetBinError(i, 0.)
                    hmc_test.SetBinContent(i, 0.)
                    hmc_test.SetBinError(i, 0.)
            
            else:
                #for i in xrange(hdata_test.FindFixBin(105+1), hdata_test.FindFixBin(150-1)+1):
                #for i in xrange(hdata_test.FindFixBin(105+1)-2, hdata_test.FindFixBin(150-1)+1):  # hide VV as well
                for i in xrange(1, hdata_test.GetNbinsX()+1):  # all bins
                    hdata_test.SetBinContent(i, 0.)
                    hdata_test.SetBinError(i, 0.)
                    hmc_test.SetBinContent(i, 0.)
                    hmc_test.SetBinError(i, 0.)
        
        
        
        setHisto(hdata_test, "data_obs")

        print "MakePlots(): Setting up the stack..."
        hs = THStack("hs", "")
        if not options.doVV and not options.doMJJ:
            hs.Add(hVVHF)
            hs.Add(hVVLF)
        elif not options.doMJJ:
            hs.Add(hVVLF)
        
        hs.Add(hQCD)
        hs.Add(hs_Top)
        hs.Add(hTT)
        hs.Add(hWj0b)
        hs.Add(hWj1b)
        hs.Add(hWj2b)
        hs.Add(hZj0b)
        hs.Add(hZj1b)
        hs.Add(hZj2b)
        if not options.doVV and not options.doMJJ:
            if plotSig:  hs.Add(hVH)
            if plotSig:  hs.Add(hZbbHinv)
        elif not options.doMJJ:
            hs.Add(hVVHF)
            if plotSig:  hs.Add(hVH)
            if plotSig:  hs.Add(hZbbHinv)

        if options.doMJJ:
            hs.Add(hVVLF)
            hs.Add(hVVHF)
            if plotSig:  hs.Add(hVH)
            if plotSig:  hs.Add(hZbbHinv)


        ymax = max(hdata_test.GetMaximum(), hs.GetMaximum())
        hs.SetMaximum(ymax * 1.7 + (sqrt(ymax) if ymax>1 else 0))
        if plotLog:
            hs.SetMaximum(ymax * 200 + (sqrt(ymax) if ymax>1 else 0))
        hs.SetMinimum(0.01)
        
        # Setup auxiliary histograms
        print "MakePlots(): Setting up auxiliary histograms..."
        staterr = hmc_exp.Clone("staterr")
        staterr.Sumw2()
        #staterr.SetFillColor(kRed)
        staterr.SetFillColor(kGray+3)
        staterr.SetMarkerSize(0)
        staterr.SetFillStyle(3013)
        
        ratio = hdata_test.Clone("ratio")
        ratio.Sumw2()
        ratio.SetMarkerSize(0.8)
        #ratio.SetMarkerSize(0.5)
        ratio.Divide(hdata_test, hmc_exp, 1., 1., "")
        
        ratiostaterr = hmc_exp.Clone("ratiostaterr")
        ratiostaterr.Sumw2()
        ratiostaterr.SetStats(0)
        ratiostaterr.SetTitle("")
        #ratiostaterr.GetXaxis().SetTitle(xtitle)
        ratiostaterr.SetTitle(";"+xtitle)
        ratiostaterr.GetYaxis().SetTitle("Data/MC")
        ratiostaterr.SetMaximum(2.2)
        ratiostaterr.SetMinimum(0)
        ratiostaterr.SetMarkerSize(0)
        #ratiostaterr.SetFillColor(kRed)
        ratiostaterr.SetFillColor(kGray+3)
        ratiostaterr.SetFillStyle(3013)
        #ratiostaterr.SetFillStyle(3001)
        #ratiostaterr.GetXaxis().CenterTitle()
        ratiostaterr.GetXaxis().SetLabelSize(0.12)
        ratiostaterr.GetXaxis().SetTitleSize(0.14)
        ratiostaterr.GetXaxis().SetTitleOffset(1.10)
        ratiostaterr.GetYaxis().CenterTitle()
        ratiostaterr.GetYaxis().SetLabelSize(0.10)
        ratiostaterr.GetYaxis().SetTitleSize(0.12)
        #ratiostaterr.GetYaxis().SetTitleSize(0.10)
        ratiostaterr.GetYaxis().SetTitleOffset(0.6)
        ratiostaterr.GetYaxis().SetNdivisions(505)
        
        ratiounity = TLine(xlow,1,xup,1)
        ratiounity.SetLineStyle(2)
        
        for i in xrange(0, hmc_exp.GetNbinsX()+2):
            ratiostaterr.SetBinContent(i, 1.0)
            if (hmc_exp.GetBinContent(i) > 1e-6):  # use smaller tolerance?
                binerror = hmc_exp.GetBinError(i) / hmc_exp.GetBinContent(i)
                ratiostaterr.SetBinError(i, binerror)
            else:
                ratiostaterr.SetBinError(i, 999.)
           
            
            #if (!(hdata_test.GetBinContent(i) > 1e-6)):
            #    ratiostaterr.SetBinError(i, 0.)
        
        
        
        ratiosysterr = ratiostaterr.Clone("ratiosysterr")
        ratiosysterr.Sumw2()
        ratiosysterr.SetMarkerSize(0)
        #ratiosysterr.SetFillColor(kBlue)
        ratiosysterr.SetFillColor(kYellow-4)
        #ratiosysterr.SetFillStyle(3002)
        ratiosysterr.SetFillStyle(1001)
        
        for i in xrange(1, hmc_exp.GetNbinsX()+1):
            if hmc_exp.GetBinContent(i) > 1e-6:  # use smaller tolerance?

                binerror2        = (pow(hmc_exp.GetBinError(i), 2) +
                                    pow(eWj0b.GetBinContent(i), 2) +
                                    pow(eWj1b.GetBinContent(i), 2) +
                                    pow(eWj2b.GetBinContent(i), 2) +
                                    pow(eZj0b.GetBinContent(i), 2) +
                                    pow(eZj1b.GetBinContent(i), 2) +
                                    pow(eZj2b.GetBinContent(i), 2) +
                                    pow(eTT.GetBinContent(i), 2) +
                                    pow(es_Top.GetBinContent(i), 2))
                if not options.doVV:
                    binerror2 += pow(eVVLF.GetBinContent(i), 2)
                    binerror2 += pow(eVVHF.GetBinContent(i), 2)
                    binerror2 += pow(eVH.GetBinContent(i), 2)
                else:
                    binerror2 += pow(eVVLF.GetBinContent(i), 2)
                    binerror2 += pow(eVH.GetBinContent(i), 2)

                binerror = sqrt(binerror2)
                ratiosysterr.SetBinError(i, binerror / hmc_exp.GetBinContent(i))
                
                if hmc_test.GetBinContent(i) > 1e-6:  # use smaller tolerance?
                    hmc_test.SetBinError(i, binerror)
                if hdata_test.GetBinContent(i) > 1e-6:
                    hdata_test.SetBinError(i, sqrt(hdata_test.GetBinContent(i)))
            
        

        # Setup legends
        print "MakePlots(): Setting up legends..."
        #leg = TLegend(0.74, 0.56, 0.92, 0.92)
        #leg.SetFillColor(0)
        #leg.SetLineColor(0)
        #leg.SetShadowColor(0)
        #leg.SetTextFont(62)
        #leg.SetTextSize(0.03)
        #leg.AddEntry(hdata_test, "Data", "p")
        #if plotSig:  leg.AddEntry(hVH, "VH(%i)" % massH, "l")
        #if plotSig:  leg.AddEntry(hZbbHinv, "ZH(inv)", "l")
        #if options.doVV:
        #    leg.AddEntry(hVVHF, "VV(b#bar{b})", "l")
        #
        #leg.AddEntry(hZj2b, "Z + b#bar{b}", "f")
        #leg.AddEntry(hZj1b, "Z + b", "f")
        #leg.AddEntry(hZj0b, "Z + udscg", "f")
        #leg.AddEntry(hWj2b, "W + b#bar{b}", "f")
        #leg.AddEntry(hWj1b, "W + b", "f")
        #leg.AddEntry(hWj0b, "W + udscg", "f")
        #leg.AddEntry(hTT, "t#bar{t}", "f")
        #leg.AddEntry(hs_Top, "single top", "f")
        #leg.AddEntry(hQCD, "QCD", "f")
        #leg.AddEntry(hVVLF, "VV(udscg)", "f")
        #if not options.doVV:
        #    leg.AddEntry(hVVHF, "VV(b#bar{b})", "f")
        #
        #leg.AddEntry(staterr, "MC uncert. (stat)", "fl")

        #leg1 = TLegend(0.58, 0.68, 0.76, 0.92)
        leg1 = TLegend(0.50, 0.60, 0.72, 0.92)
        leg1.SetFillColor(0)
        leg1.SetLineColor(0)
        leg1.SetShadowColor(0)
        leg1.SetTextFont(62)
        leg1.SetTextSize(0.03)
        leg1.AddEntry(hdata_test, "Data", "p")
        if plotSig:  leg1.AddEntry(hVH, "VH(%i)" % massH, "l")
        if plotSig:  leg1.AddEntry(hZbbHinv, "ZH(inv)", "l")
        if options.doVV:
            leg1.AddEntry(hVVHF, "VV(b#bar{b})", "l")
        
        leg1.AddEntry(hTT, "t#bar{t}", "f")
        leg1.AddEntry(hs_Top, "single top", "f")
        leg1.AddEntry(hQCD, "QCD", "f")
        leg1.AddEntry(hVVLF, "VV(udscg)", "f")
        if not options.doVV:
            leg1.AddEntry(hVVHF, "VZ(b#bar{b})", "f")


        #leg2 = TLegend(0.72, 0.60, 0.94, 0.92)
        leg2 = TLegend(0.72, 0.60, 0.94, 0.88)
        leg2.SetFillColor(0)
        leg2.SetLineColor(0)
        leg2.SetShadowColor(0)
        leg2.SetTextFont(62)
        leg2.SetTextSize(0.03)
        leg2.AddEntry(hWj2b, "W + b#bar{b}", "f")
        leg2.AddEntry(hWj1b, "W + b", "f")
        leg2.AddEntry(hWj0b, "W + udscg", "f")
        leg2.AddEntry(hZj2b, "Z + b#bar{b}", "f")
        leg2.AddEntry(hZj1b, "Z + b", "f")
        leg2.AddEntry(hZj0b, "Z + udscg", "f")
        leg2.AddEntry(staterr, "MC uncert. (stat)", "f")

        ratioleg1 = TLegend(0.54, 0.88, 0.72, 0.96)
        #ratioleg1 = TLegend(0.50, 0.86, 0.69, 0.96)
        ratioleg1.AddEntry(ratiostaterr, "MC uncert. (stat)", "f")
        ratioleg1.SetFillColor(0)
        ratioleg1.SetLineColor(0)
        ratioleg1.SetShadowColor(0)
        ratioleg1.SetTextFont(62)
        ratioleg1.SetTextSize(0.06)
        ratioleg1.SetBorderSize(1)
        
        ratioleg2 = TLegend(0.72, 0.88, 0.95, 0.96)
        #ratioleg2 = TLegend(0.69, 0.86, 0.9, 0.96)
        ratioleg2.AddEntry(ratiosysterr, "MC uncert. (stat+syst)", "f")
        ratioleg2.SetFillColor(0)
        ratioleg2.SetLineColor(0)
        ratioleg2.SetShadowColor(0)
        ratioleg2.SetTextFont(62)
        ratioleg2.SetTextSize(0.06)
        ratioleg2.SetBorderSize(1)
        
        #ratioleg1 = TLegend(0.72, 0.88, 0.94, 0.96)
        ##ratioleg1 = TLegend(0.50, 0.86, 0.69, 0.96)
        #ratioleg1.AddEntry(ratiostaterr, "MC uncert. (stat)", "f")
        #ratioleg1.SetFillColor(0)
        #ratioleg1.SetLineColor(0)
        #ratioleg1.SetShadowColor(0)
        #ratioleg1.SetTextFont(62)
        ##ratioleg1.SetTextSize(0.06)
        #ratioleg1.SetTextSize(0.07)
        #ratioleg1.SetBorderSize(1)

        # Draw MC signal and background
        print "MakePlots(): Drawing..."
        pad1.cd()
        if plotLog: pad1.SetLogy()
        hs.Draw("hist")
        hs.GetXaxis().SetLabelSize(0)
        binwidth = (xup - xlow) / nbins_plot
        ytitle = "Events / %.3f" % binwidth
        hs.GetYaxis().SetTitle(ytitle)
        if " ; " in xtitle:
            hs.SetTitle(";"+xtitle)
        
        staterr.Draw("e2 same")
        if plotSig:
            hVH.SetLineWidth(3)
            hVH.SetFillColor(0)
            hVH.Draw("hist same")

        if plotSig:
            hZbbHinv.SetLineWidth(3)
            hZbbHinv.SetFillColor(0)
            hZbbHinv.Draw("hist same")

        if options.doVV:
            hVVHF.SetLineWidth(3)
            hVVHF.SetLineColor(kGray + 2)
            hVVHF.SetFillColor(0)
            hVVHF.Draw("hist same")

        
        # Draw data
        hdata_test.Draw("e1 same")
        
        # Draw legends
        #leg.Draw()
        leg1.Draw()
        leg2.Draw()
        latex = TLatex()
        latex.SetNDC()
        latex.SetTextAlign(12)
        latex.SetTextFont(62)
        latex.SetTextSize(0.052)
        latex.DrawLatex(0.19, 0.89, "CMS Preliminary")
        latex.SetTextSize(0.04)

        latex.DrawLatex(0.19, 0.84, "#sqrt{s} = 8 TeV, L = 18.9 fb^{-1}")




        #latex.DrawLatex(0.19, 0.79, "Z(#nu#bar{#nu})H(b#bar{b})")
        latex.DrawLatex(0.19, 0.79, "Z(b#bar{b})H(inv)")
        
        # Under/overflows a la TMVA
        #uoflow = "U/O-flow (Data,MC): (%.1f, %.1f) / (%.1f, %.1f)" % (hdata_test.GetBinContent(0), hmc_exp.GetBinContent(0), hdata_test.GetBinContent(nbins_plot+1), hmc_exp.GetBinContent(nbins_plot+1))
        #latex2 = TLatex(0.99, 0.1, uoflow)
        #latex2.SetNDC()
        #latex2.SetTextSize(0.02)
        #latex2.SetTextAngle(90)
        #latex2.AppendPad()
        
        # Draw ratio
        pad2.cd()
        pad2.SetGridy(0)
        ratiostaterr.Draw("e2")
        ratiosysterr.Draw("e2 same")
        ratiostaterr.Draw("e2 same")
        ratiounity.Draw()
        ratio.Draw("e1 same")
        
        # Draw ratio legends
        ratioleg1.Draw()
        ratioleg2.Draw()

        # Kolmogorov-Smirnov test and Chi2 test
        pave = TPaveText(0.18, 0.86, 0.35, 0.96, "brNDC")
        #pave = TPaveText(0.18, 0.86, 0.28, 0.96, "brNDC")
        pave.SetTextAlign(12)
        pave.SetLineColor(0)
        pave.SetFillColor(0)
        pave.SetShadowColor(0)
        pave.SetBorderSize(1)
        nchisq = hdata_test.Chi2Test(hmc_test, "UWCHI2/NDF")  # MC uncert. (stat)
        kolprob = hdata_test.KolmogorovTest(hmc_test)  # MC uncert. (stat)
        text = pave.AddText("#chi_{#nu}^{2} = %.3f, K_{s} = %.3f" % (nchisq, kolprob))
        #text = pave.AddText("#chi_{#nu}^{2} = %.3f" % (nchisq))
        text.SetTextFont(62)
        text.SetTextSize(0.06)
        #text.SetTextSize(0.07)
        pave.Draw()

        print "MakePlots(): Printing..."
        pad1.cd()
        gPad.RedrawAxis()
        gPad.Modified()
        gPad.Update()
        pad2.cd()
        gPad.RedrawAxis()
        gPad.Modified()
        gPad.Update()
        c1.cd()
        #plotname = TString("%s_%s_%s" % (channel, region, var))
        #FormatFileName(plotname)
        #gPad.Print(plotdir+plotname.Data()+".png")
        #gPad.Print(plotdir+plotname.Data()+".pdf")
        gPad.Print(varname+"_%s.png" % whatfit)
        gPad.Print(varname+"_%s.pdf" % whatfit)

        
    return 0


## Main
## ----------------------------------------------------------------------------

if __name__ == "__main__":
    # Note: Must have at least 1 shape systematic for each lnN systematic!
    
    from optparse import OptionParser

    # import ROOT with a fix to get batch mode (http://root.cern.ch/phpBB3/viewtopic.php?t=3198)
    sys.argv.append( '-b-' )
    gROOT.SetBatch(True)
    sys.argv.remove( '-b-' )

    # Setup OptionParser
    parser = OptionParser(usage="usage: %prog [options] datacard.txt -i mlfit.root \nrun with --help to get list of options")
    addDatacardParserOptions(parser)
    parser.add_option("-i", "--infile", type="string", dest="infile", help="input file with fit results [default: %default]")
    parser.add_option("--fit_s", action="store_true", dest="fit_s", default=True, help="use 'fit_s' [default: %default]")
    parser.add_option("--fit_b", action="store_false", dest="fit_s", help="use 'fit_b' instead of 'fit_s' [default: %default]")
    
    # Taken from datacardDump.py
    parser.add_option("-C", "--channel", type="string", dest="channel", default=None, help="Channel to dump")
    parser.add_option("-N", "--norm-only", dest="norm", default=False, action="store_true", help="Include only normalization uncertainties, not shape ones") 
    #parser.add_option("-f", "--format", type="string", dest="format", default="%8.3f +/- %6.3f", help="Format for output number")
    parser.add_option("--xs", "--exclude-syst", type="string", dest="excludeSyst", default=[], action="append", help="Systematic to exclude (regexp)")
    parser.add_option("--xp", "--exclude-process", type="string", dest="processVetos", default=[], action="append", help="Exclude processes that match this regexp; can specify multiple ones")
    #parser.add_option("-m", "--mass",     dest="mass",     default=0,  type="float",  help="Higgs mass to use. Will also be written in the Workspace as RooRealVar 'MH'.")
    #parser.add_option("-D", "--dataset",  dest="dataname", default="data_obs",  type="string",  help="Name of the observed dataset")
    parser.add_option("-p", "--poi", type="string", dest="poi", default="r", help="Name of signal strength parameter [default: %default]")
    parser.add_option("--mjj", action="store_true", dest="doMJJ", default=False, help="Do Mjj analysis [default: %default]")
    parser.add_option("--vv", action="store_true", dest="doVV", default=False, help="Do VV analysis [default: %default]")
    #parser.set_defaults(out="out.root")
    (options, args) = parser.parse_args()
    options.stat = False
    options.bin = True # fake that is a binary output, so that we parse shape lines

    argc = 1
    if len(sys.argv) < argc+1:
        parser.print_usage()
        sys.exit(1)

    # Taken from findRedundantSystematics.py
    options.fileName = args[0]
    if options.fileName.endswith(".gz"):
        import gzip
        datacardfile = gzip.open(options.fileName, "rb")
        options.fileName = options.fileName[:-3]
    elif options.fileName.endswith(".txt"):
        datacardfile = open(options.fileName, "r")
    else:
        raise RuntimeError("Cannot find the datacard, given: %s" % options.fileName)
    if "root" not in options.infile:
        raise RuntimeError("infile is not a ROOT file.")
    options.pdfname = options.fileName.replace("_8TeV.txt", "_postfit.pdf")
    options.pdfname = options.pdfname.replace(".txt", "_postfit.pdf")
    options.outname = options.pdfname.replace(".pdf", ".root")
    
    options.varnames = {  # Retarded RooFit never provides a way to access variable name without iostream
        #"WenLowPt_8TeV"         : "CMS_vhbb_BDT_Wln_8TeV",
        #"WenHighPt_8TeV"        : "CMS_vhbb_BDT_Wln_8TeV",
        #"WenHighPtLooseCSV_8TeV": "CMS_vhbb_BDT_Wln_8TeV",
        #"WmnLowPt_8TeV"         : "CMS_vhbb_BDT_Wln_8TeV",
        #"WmnHighPt_8TeV"        : "CMS_vhbb_BDT_Wln_8TeV",
        #"WmnHighPtLooseCSV_8TeV": "CMS_vhbb_BDT_Wln_8TeV",
        #"ZeeHighPt_8TeV"        : "CMS_vhbb_BDT_Zll_8TeV",
        #"ZeeLowPt_8TeV"         : "CMS_vhbb_BDT_Zll_8TeV",
        #"ZmmHighPt_8TeV"        : "CMS_vhbb_BDT_Zll_8TeV",
        #"ZmmLowPt_8TeV"         : "CMS_vhbb_BDT_Zll_8TeV",
        "ZnunuHighPt_8TeV"      : "CMS_vhbb_BDT_ZnunuHighPt_8TeV", 
        "ZnunuMedPt_8TeV"       : "CMS_vhbb_BDT_ZnunuMedPt_8TeV", 
        "ZnunuLowPt_8TeV"       : "CMS_vhbb_BDT_ZnunuLowPt_8TeV", 
        #"ZnunuLowCSV_8TeV"      : "CMS_vhbb_BDT_ZnunuLowCSV_8TeV", 
        }
    
    DC = parseCard(datacardfile, options)
    MB = ShapeBuilder(DC, options)
    
    
    nuisances = {}
    read_datacard(nuisances, DC, MB, options)
    read_mlfit(nuisances, options)
    
    histograms = get_postfit(nuisances, options)
    put_histograms(histograms, options)
    MakePlots("ZnunuHighPt_8TeV", "prefit", options)
    MakePlots("ZnunuHighPt_8TeV", "postfit", options)
    MakePlots("ZnunuMedPt_8TeV", "prefit", options)
    MakePlots("ZnunuMedPt_8TeV", "postfit", options)
    MakePlots("ZnunuLowPt_8TeV", "prefit", options)
    MakePlots("ZnunuLowPt_8TeV", "postfit", options)

#python bsmdcplotter.py vhbb_Znn_J14_bbb_combo_8TeV.txt -i mlfit_combo.root
#python bsmdcplotter.py vhbb_dcplotter_BDT.txt -i vhbb_dcplotter_mlfit_BDT.root
#python bsmdcplotter.py zhinv_dcplotter_BDT.txt -i zhinv_dcplotter_mlfit_BDT.root
