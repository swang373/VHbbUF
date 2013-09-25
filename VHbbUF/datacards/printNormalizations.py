#!/usr/bin/env python
import re
from sys import argv
from math import sqrt
from optparse import OptionParser
from ordereddict import OrderedDict

# import ROOT with a fix to get batch mode (http://root.cern.ch/phpBB3/viewtopic.php?t=3198)
hasHelp = False
for X in ("-h", "-?", "--help"):
    if X in argv:
        hasHelp = True
        argv.remove(X)
argv.append( '-b-' )
from ROOT import gROOT, gSystem, TFile
gROOT.SetBatch(True)
#gSystem.Load("libHiggsAnalysisCombinedLimit.so")
argv.remove( '-b-' )
if hasHelp: argv.append("-h")

parser = OptionParser(usage="usage: %prog [options] mlfit.root  \nrun with --help to get list of options")
parser.add_option("--xlow", dest="xlow", default=-999., type="float", help="Lower bound [default: %default]")
parser.add_option("--xup", dest="xup", default=999., type="float", help="Upper bound [default: %default]")
parser.add_option("-f", "--format", dest="format", default="text", type="string", help="Output format ('text', 'latex', 'twiki', 'html')")
parser.add_option("--plots", action="store_true", dest="plots", default=False, help="Draw stacked histograms")
parser.add_option("--workspace", dest="workspace", default="zhinv_Zbb_8TeV.root", type="string", help="RooWorkspace file [default: %default]")
parser.add_option("--mjj", action="store_true", dest="doMJJ", default=False, help="Do Mjj analysis [default: %default]")
parser.add_option("--vv", action="store_true", dest="doVV", default=False, help="Do VV analysis [default: %default]")

(options, args) = parser.parse_args()
if len(args) == 0:
    parser.print_usage()
    exit(1)

tfile = TFile.Open(args[0])
if tfile == None: raise RuntimeError, "Cannot open file %s" % args[0]

norms_all = []
for fit in ["prefit", "fit_b", "fit_s"]:
    norms = {}
    
    fitnorm = tfile.Get("norm_"+fit)
    fitdir = tfile.Get("shapes_"+fit)
    channels = [key.ReadObj() for key in fitdir.GetListOfKeys() if key.IsFolder()]
    for channel in channels:
        processes = [key.ReadObj() for key in channel.GetListOfKeys() if key.GetClassName() == "TH1F"]
        for process in processes:
            key = channel.GetName()+"/"+process.GetName()
            if process.GetName() != "data_obs":
                processnorm = fitnorm[key].getVal()
                processnormerr = fitnorm[key].getError()
            else:
                processnorm = process.Integral()
                processnormerr = sqrt(process.Integral())
            
            fraction = 1.0
            if options.xlow != -999. or options.xup != 999.:
                bxlow = process.FindFixBin(options.xlow)
                bxup = process.FindFixBin(options.xup)
                if process.Integral() > 1e-3:
                    fraction = process.Integral(bxlow, bxup) / process.Integral()
            norms[key] = (processnorm * fraction, processnormerr * fraction)
        # Create zero norms
        for process in ["QCD"]:
            key = channel.GetName()+"/"+process
            if key not in norms:
                norms[key] = (0., 0.)
    norms_all.append(norms)
    
dprocesses = OrderedDict([
    ("ZH", "\ZbbHinv"),
    ("ZH_SM", "\ZnnH (SM)"),
    ("WH_SM", "\WlnH (SM)"),
    ("ZZ", "ZZ(bb)"),
    ("WZ", "WZ(bb)"),
    ("VVLF", "VV(udscg)"),
    ("Zj2b", "\Zbb"),
    ("Zj1b", "\Zb"),
    ("Zj0b", "\Zudscg"),
    ("Wj2b", "\Wbb"),
    ("Wj1b", "\Wb"),
    ("Wj0b", "\Wudscg"),
    ("TT", "\\ttbar"),
    ("s_Top", "single top"),
    ("QCD", "QCD"),
    ("total_background", "Total background")
])

dchannels = OrderedDict([
    ("ZnunuHighPt_8TeV", "High \ptV"),
    ("ZnunuMedPt_8TeV", "Intermediate \ptV"),
    ("ZnunuLowPt_8TeV", "Low \ptV"),
])

for i in xrange(len(norms_all)):
    print "%-16s " % "Process",
    for channel in dchannels.keys():
        print " & %-19s " % ("  "+dchannels[channel]),
    print "\\\\"
    for process in dprocesses.keys():
        print "%-16s " % dprocesses[process],
        for channel in dchannels.keys():
            key = channel+"/"+process
            print " & %7.1f $\pm$ %5.1f " % norms_all[i][key],
        print "\\\\"


#_______________________________________________________________________________
def MakePlots(processes, plotname, whatfit, options):
    massH    = 125
    plotData = True
    plotLog  = True
    plotSig  = True
    plotdir  = "./"
    kRed     = 632
    kYellow  = 400
    kBlue    = 600
    kGray    = 920
    if whatfit == "prefit":
        whatfit_s = "pre-fit"
    else:
        whatfit_s = "post-fit"
    
    if options.doMJJ:
        plotLog = False
    
    hZH        = processes["ZH_SM"]
    hWH        = processes["WH_SM"]
    hZbbHinv   = processes["ZH"]
    hWj0b      = processes["Wj0b"]
    hWj1b      = processes["Wj1b"]
    hWj2b      = processes["Wj2b"]
    hZj0b      = processes["Zj0b"]
    hZj1b      = processes["Zj1b"]
    hZj2b      = processes["Zj2b"]
    hTT        = processes["TT"]
    hs_Top     = processes["s_Top"]
    hVVLF      = processes["VVLF"]
    hWZHF      = processes["WZ"]
    hZZHF      = processes["ZZ"]
    hQCD       = processes["QCD"]
    hdata_obs  = processes["data_obs"]
    hmc_exp    = processes["total_background"]
    
    eZH        = processes["ZH_SM"+"_err"]
    eWH        = processes["WH_SM"+"_err"]
    eZbbHinv   = processes["ZH"+"_err"]
    eWj0b      = processes["Wj0b"+"_err"]
    eWj1b      = processes["Wj1b"+"_err"]
    eWj2b      = processes["Wj2b"+"_err"]
    eZj0b      = processes["Zj0b"+"_err"]
    eZj1b      = processes["Zj1b"+"_err"]
    eZj2b      = processes["Zj2b"+"_err"]
    eTT        = processes["TT"+"_err"]
    es_Top     = processes["s_Top"+"_err"]
    eVVLF      = processes["VVLF"+"_err"]
    eWZHF      = processes["WZ"+"_err"]
    eZZHF      = processes["ZZ"+"_err"]
    eQCD       = processes["QCD"+"_err"]
    
    varname = hdata_obs.GetXaxis().GetTitle()
    xtitle = "BDT"
    if options.doMJJ:  xtitle = "m_{T} [GeV]"
    nbins = hdata_obs.GetNbinsX()
    xlow = hdata_obs.GetXaxis().GetXmin()
    xup = hdata_obs.GetXaxis().GetXmax()
    hVH        = TH1F("VH"       , "", nbins, xlow, xup)
    hVVHF      = TH1F("VVHF"     , "", nbins, xlow, xup)
    hVV        = TH1F("VV"       , "", nbins, xlow, xup)
    #hmc_exp    = TH1F("mc_exp"   , "", nbins, xlow, xup)
    
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

    #hmc_exp.Add(hWj0b)
    #hmc_exp.Add(hWj1b)
    #hmc_exp.Add(hWj2b)
    #hmc_exp.Add(hZj0b)
    #hmc_exp.Add(hZj1b)
    #hmc_exp.Add(hZj2b)
    #hmc_exp.Add(hTT)
    #hmc_exp.Add(hs_Top)
    #if not options.doVV:
    #    hmc_exp.Add(hVVLF)
    #    hmc_exp.Add(hVVHF)
    #    hmc_exp.Add(hVH)  # VH is counted as background
    #else:
    #    hmc_exp.Add(hVVLF)
    #    hmc_exp.Add(hVH)  # VH is counted as background
    #hmc_exp.Add(hQCD)
    
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
        leg2.AddEntry(staterr, "MC uncert. (%s)" % whatfit_s, "f")

        #ratioleg1 = TLegend(0.54, 0.88, 0.72, 0.96)
        ##ratioleg1 = TLegend(0.50, 0.86, 0.69, 0.96)
        #ratioleg1.AddEntry(ratiostaterr, "MC uncert. (stat)", "f")
        #ratioleg1.SetFillColor(0)
        #ratioleg1.SetLineColor(0)
        #ratioleg1.SetShadowColor(0)
        #ratioleg1.SetTextFont(62)
        #ratioleg1.SetTextSize(0.06)
        #ratioleg1.SetBorderSize(1)
        #
        #ratioleg2 = TLegend(0.72, 0.88, 0.95, 0.96)
        ##ratioleg2 = TLegend(0.69, 0.86, 0.9, 0.96)
        #ratioleg2.AddEntry(ratiosysterr, "MC uncert. (stat+syst)", "f")
        #ratioleg2.SetFillColor(0)
        #ratioleg2.SetLineColor(0)
        #ratioleg2.SetShadowColor(0)
        #ratioleg2.SetTextFont(62)
        #ratioleg2.SetTextSize(0.06)
        #ratioleg2.SetBorderSize(1)
        
        ratioleg1 = TLegend(0.72, 0.88, 0.94, 0.96)
        #ratioleg1 = TLegend(0.50, 0.86, 0.69, 0.96)
        ratioleg1.AddEntry(ratiostaterr, "MC uncert. (%s)" % whatfit_s, "f")
        ratioleg1.SetFillColor(0)
        ratioleg1.SetLineColor(0)
        ratioleg1.SetShadowColor(0)
        ratioleg1.SetTextFont(62)
        #ratioleg1.SetTextSize(0.06)
        ratioleg1.SetTextSize(0.07)
        ratioleg1.SetBorderSize(1)

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
        #ratiosysterr.Draw("e2 same")
        ratiostaterr.Draw("e2 same")
        ratiounity.Draw()
        ratio.Draw("e1 same")
        
        # Draw ratio legends
        ratioleg1.Draw()
        #ratioleg2.Draw()

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
        gPad.Print(plotdir+plotname+"_"+whatfit+".png")
        gPad.Print(plotdir+plotname+"_"+whatfit+".pdf")
        
    return 0

if options.plots:
    gROOT.ProcessLine(".L tdrstyle.C")
    gROOT.ProcessLine("setTDRStyle()")
    gROOT.ProcessLine(".L HelperBDTShape.h")
    from ROOT import TH1, TH1F, TColor, TCanvas, TPad, TLegend, THStack, TLatex, TPaveText, TLine, RooFit, RooAbsData, gPad, setHisto
    from array import array
    TH1.SetDefaultSumw2(1)
    
    wstfile = TFile.Open(options.workspace)
    if not options.doMJJ:
        varnameprefix = "zhinv_Zbb_BDT_"
    else:
        varnameprefix = "zhinv_Zbb_MJJ_"
    
    for whatfit in ["prefit", "fit_s"]:
        for channel in dchannels.keys():
            processes_plots = {}
            for process in dprocesses.keys():
                hist = tfile.Get("shapes_"+whatfit+"/"+channel+"/"+process)
                if hist:
                    hist.SetDirectory(0)
                    #key = channel+"/"+process
                    key = process
                    processes_plots[key] = hist
                    processes_plots[key+"_err"] = hist.Clone(process+"_err")
            
            # Create empty histograms
            for process in ["QCD"]:
                key = process
                if key not in processes_plots:
                    hist = processes_plots["total_background"].Clone(key)
                    hist.Reset()
                    hist.SetDirectory(0)
                    processes_plots[key] = hist
                    processes_plots[key+"_err"] = hist.Clone(process+"_err")
            
            # Need to get data_obs from somewhere else...
            process = "data_obs"
            ws = wstfile.Get(channel)
            rrv = ws.var(varnameprefix+channel)
            hist = ws.data(process).createHistogram(channel+"/"+process, rrv)
            if hist:
                # Set the proper errors and stats
                arr = array('d', [0., 0., 0., 0.])  # sumw, sumw2, sumwx, sumwx2
                for i in xrange(0,hist.GetNbinsX()+2):
                    hist.SetBinError(i, sqrt(hist.GetBinContent(i)))
                    x = hist.GetBinCenter(i)
                    w = hist.GetBinContent(i)
                    err = hist.GetBinError(i)
                    arr[0] += w
                    arr[1] += err*err
                    arr[2] += w*x
                    arr[3] += w*x*x
                hist.PutStats(arr)
                hist.SetDirectory(0)
                #key = channel+"/"+process
                key = process
                processes_plots[key] = hist
                processes_plots[key+"_err"] = hist.Clone(process+"_err")
            
            plotname = varnameprefix+channel
            MakePlots(processes_plots, plotname, whatfit, options)
    