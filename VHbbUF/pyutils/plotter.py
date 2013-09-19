from ROOT import TCanvas, TCut, TChain, TH1, TH1F, THStack, TLatex, TLegend, TPad, TPaveText, TString, gROOT, gPad, gStyle, gSystem
gROOT.ProcessLine(".L tdrstyle.C")
gROOT.ProcessLine("setTDRStyle()")
gROOT.ProcessLine(".L HelperBDTShape.h")
gROOT.ProcessLine(".L HelperFunctions.h")
from ROOT import setHisto, FormatFileName
from math import sqrt


channel  = "ZnunuHighPt"
region   = "TT"
massH    = 125
plotData = True
plotLog  = False
plotSig  = True
plotdir  = "plots_etc/"
kRed     = 632
kYellow  = 400
kBlue    = 600

TH1.SetDefaultSumw2(1)
gROOT.SetBatch(1)
if gSystem.AccessPathName(plotdir):
    gSystem.mkdir(plotdir)

class EventsJ11(object):
    pass


def Read(whatstep="Step4", cutallmc="", cutalldata=""):
    
    if whatstep == "Step4":
        indir    = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130404/stitch/"
        #indir    = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130314/stitch/"
        prefix   = "Step4_"
        suffix   = ".root"
        treename = "tree_%s_%s" % (channel, "ctrl")
        
        ZH = TChain(treename)
        ZH.Add(indir + prefix + "ZH%i" % massH + suffix)
        
        WH = TChain(treename)
        WH.Add(indir + prefix + "WH%i" % massH + suffix)
        
        Wj = TChain(treename)
        Wj.Add(indir + prefix + "Wj" + suffix)
        
        Zj = TChain(treename)
        Zj.Add(indir + prefix + "Zj" + suffix)
        
        TT = TChain(treename)
        TT.Add(indir + prefix + "TT" + suffix)
        
        s_Top = TChain(treename)
        s_Top.Add(indir + prefix + "s_Top" + suffix)
        
        VV = TChain(treename)
        VV.Add(indir + prefix + "VV" + suffix)
        
        QCD = TChain(treename)
        QCD.Add(indir + prefix + "QCD" + suffix)
        
        data_obs = TChain(treename)
        data_obs.Add(indir + prefix + "data_obs" + suffix)

        cutHF = TCut("eventFlav==5")
        cutLF = TCut("eventFlav!=5")
        
        cut2b = TCut("abs(hJet_flavour[0])==5 && abs(hJet_flavour[1])==5")
        cut1b = TCut("(abs(hJet_flavour[0])==5 && abs(hJet_flavour[1])!=5) || (abs(hJet_flavour[0])!=5 && abs(hJet_flavour[1])==5)")
        cut0b = TCut("abs(hJet_flavour[0])!=5 && abs(hJet_flavour[1])!=5")
    
    elif whatstep == "Step3":
        indir    = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/skim_ZnnH_Step3/"
        prefix   = "Step3_"
        suffix   = ".root"
        treename = "tree"
        
        ZH = TChain(treename)
        #ZH.Add(indir + prefix + "ZnnH%i" % massH + suffix)
        
        WH = TChain(treename)
        #WH.Add(indir + prefix + "WlnH%i" % massH + suffix)
        
        Wj = TChain(treename)
        Wj.Add(indir + prefix + "Wj" + suffix)
        #names = ["WJetsPtW70", "WJetsPtW100", "WJetsPtW180",]
        #for n in names:
        #    Wj.Add(indir + prefix + n + suffix)
        
        Zj = TChain(treename)
        #Zj.Add(indir + prefix + "Zj" + suffix)
        ##names = ["ZJetsHT50", "ZJetsHT100", "ZJetsHT200", "ZJetsHT400", "ZJetsPtZ100",]
        ##for n in names:
        ##    Zj.Add(indir + prefix + n + suffix)
        
        TT = TChain(treename)
        TT.Add(indir + prefix + "TT" + suffix)
        #names = ["TTFullLeptMG", "TTSemiLeptMG", "TTHadronicMG",]
        #for n in names:
        #    TT.Add(indir + prefix + n + suffix)
        
        s_Top = TChain(treename)
        s_Top.Add(indir + prefix + "s_Top" + suffix)
        #names = ["T_t", "Tbar_t", "T_s", "Tbar_s", "T_tW", "Tbar_tW",]
        #for n in names:
        #    s_Top.Add(indir + prefix + n + suffix)
        
        VV = TChain(treename)
        VV.Add(indir + prefix + "VV" + suffix)
        #names = ["WW", "WZ", "ZZ",]
        #for n in names:
        #    VV.Add(indir + prefix + n + suffix)
        
        QCD = TChain(treename)
        QCD.Add(indir + prefix + "QCD" + suffix)
        ##names = ["QCDPt50", "QCDPt80", "QCDPt120", "QCDPt170", "QCDPt300", "QCDPt470", "QCDPt600", "QCDPt800", "QCDPt1000", "QCDPt1400", "QCDPt1800",]
        ##for n in names:
        ##    QCD.Add(indir + prefix + n + suffix)
        
        data_obs = TChain(treename)
        data_obs.Add(indir + prefix + "data_obs" + suffix)
        #names = ["Data_R", "Data_P",]
        #for n in names:
        #    data_obs.Add(indir + prefix + n + suffix)

        cutHF = TCut("eventFlav==5")
        cutLF = TCut("eventFlav!=5")
        
        cut2b = TCut("abs(hJet_flavour[0])==5 && abs(hJet_flavour[1])==5")
        cut1b = TCut("(abs(hJet_flavour[0])==5 && abs(hJet_flavour[1])!=5) || (abs(hJet_flavour[0])!=5 && abs(hJet_flavour[1])==5)")
        cut0b = TCut("abs(hJet_flavour[0])!=5 && abs(hJet_flavour[1])!=5")
        
    ev = EventsJ11()
    #ev.ZH = ZH.CopyTree(cutallmc.GetTitle())
    #ev.WH = WH.CopyTree(cutallmc.GetTitle())
    ev.Wj0b = Wj.CopyTree((cutallmc + cut0b).GetTitle())
    ev.Wj1b = Wj.CopyTree((cutallmc + cut1b).GetTitle())
    ev.Wj2b = Wj.CopyTree((cutallmc + cut2b).GetTitle())
    #ev.Zj0b = Zj.CopyTree((cutallmc + cut0b).GetTitle())
    #ev.Zj1b = Zj.CopyTree((cutallmc + cut1b).GetTitle())
    #ev.Zj2b = Zj.CopyTree((cutallmc + cut2b).GetTitle())
    ev.TT = TT.CopyTree(cutallmc.GetTitle())
    ev.s_Top = s_Top.CopyTree(cutallmc.GetTitle())
    ev.VVLF = VV.CopyTree((cutallmc + cutLF).GetTitle())
    ev.VVHF = VV.CopyTree((cutallmc + cutHF).GetTitle())
    #ev.QCD = QCD.CopyTree(cutallmc.GetTitle())
    ev.data_obs = data_obs.CopyTree(cutalldata.GetTitle())
    
    return ev


def MakePlots(ev, var, cutmc, cutdata, xtitle, nbins, xlow, xup):
    
    hZH        = TH1F("ZH"       , "", nbins, xlow, xup)
    hWH        = TH1F("WH"       , "", nbins, xlow, xup)
    hWj0b      = TH1F("Wj0b"     , "", nbins, xlow, xup)
    hWj1b      = TH1F("Wj1b"     , "", nbins, xlow, xup)
    hWj2b      = TH1F("Wj2b"     , "", nbins, xlow, xup)
    hZj0b      = TH1F("Zj0b"     , "", nbins, xlow, xup)
    hZj1b      = TH1F("Zj1b"     , "", nbins, xlow, xup)
    hZj2b      = TH1F("Zj2b"     , "", nbins, xlow, xup)
    hTT        = TH1F("TT"       , "", nbins, xlow, xup)
    hs_Top     = TH1F("s_Top"    , "", nbins, xlow, xup)
    hVVLF      = TH1F("VVLF"     , "", nbins, xlow, xup)
    hVVHF      = TH1F("VVHF"     , "", nbins, xlow, xup)
    hQCD       = TH1F("QCD"      , "", nbins, xlow, xup)
    hVH        = TH1F("VH"       , "", nbins, xlow, xup)
    hVV        = TH1F("VV"       , "", nbins, xlow, xup)
    hmc_exp    = TH1F("mc_exp"   , "", nbins, xlow, xup)
    hdata_obs  = TH1F("data_obs" , "", nbins, xlow, xup)

    #ev.ZH.Project("ZH", var, cutmc.GetTitle())
    #ev.WH.Project("WH", var, cutmc.GetTitle())

    # Apply slope
    WJSlopeErr = 0.0020
    ZJSlopeErr = 0.0025
    WJSlope = -0.50
    ZJSlope = -1.00
    cutwjslope = TCut("apply_pt_slope(genW.pt, %f, 150)" % (WJSlope * WJSlopeErr))
    cutzjslope = TCut("apply_pt_slope(genZ.pt, %f, 130)" % (ZJSlope * ZJSlopeErr))

    # Apply scale factors
    scalefactors = (0.930, 2.124, 0.700, 1.175, 2.135, 1.125, 0.991)
    cutsf = TCut("%f" % scalefactors[0])
    ev.Wj0b.Project("Wj0b", var, (cutmc * cutsf * cutwjslope).GetTitle())
    cutsf = TCut("%f" % scalefactors[1])
    ev.Wj1b.Project("Wj1b", var, (cutmc * cutsf * cutwjslope).GetTitle())
    cutsf = TCut("%f" % scalefactors[2])
    ev.Wj2b.Project("Wj2b", var, (cutmc * cutsf * cutwjslope).GetTitle())
    cutsf = TCut("%f" % scalefactors[3])
    #ev.Zj0b.Project("Zj0b", var, (cutmc * cutsf * cutzjslope).GetTitle())
    cutsf = TCut("%f" % scalefactors[4])
    #ev.Zj1b.Project("Zj1b", var, (cutmc * cutsf * cutzjslope).GetTitle())
    cutsf = TCut("%f" % scalefactors[5])
    #ev.Zj2b.Project("Zj2b", var, (cutmc * cutsf * cutzjslope).GetTitle())
    cutsf = TCut("%f" % scalefactors[6])
    ev.TT.Project("TT", var, (cutmc * cutsf).GetTitle())
    
    ev.s_Top.Project("s_Top", var, cutmc.GetTitle())
    ev.VVLF.Project("VVLF", var, cutmc.GetTitle())
    ev.VVHF.Project("VVHF", var, cutmc.GetTitle())

    # Apply QCD scale factor
    scalefactorQCD = 0.005
    cutsf = TCut("%f" % scalefactorQCD)
    #ev.QCD.Project("QCD", var, (cutmc * cutsf).GetTitle())

    ev.data_obs.Project("data_obs", var, cutdata.GetTitle())
    
    hVH.Add(hZH)
    hVH.Add(hWH)

    hVV.Add(hVVLF)
    hVV.Add(hVVHF)

    hmc_exp.Add(hWj0b)
    hmc_exp.Add(hWj1b)
    hmc_exp.Add(hWj2b)
    hmc_exp.Add(hZj0b)
    hmc_exp.Add(hZj1b)
    hmc_exp.Add(hZj2b)
    hmc_exp.Add(hTT)
    hmc_exp.Add(hs_Top)
    hmc_exp.Add(hVV)
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
        
        hdata_plot = hdata_obs.Clone("hdata_plot")  # blinded plot
        hmc_test = hmc_exp.Clone("hmc_test")  # for chi2 and KS test
        hdata_plot.Sumw2()
        hmc_test.Sumw2()
        nbins_plot = hdata_plot.GetNbinsX()
        assert(nbins_plot == hmc_test.GetNbinsX())
        if not plotData:  # be blind to the most sensitive bins
            for i in xrange(max(int(nbins_plot*0.75), nbins_plot-5),nbins_plot+2):
                hdata_plot.SetBinContent(i, 0.)
                hdata_plot.SetBinError(i, 0.)
                hmc_test.SetBinContent(i, 0.)
                hmc_test.SetBinError(i, 0.)
                
        
        setHisto(hdata_plot, "data_obs")

        print "MakePlots(): Setting up the stack..."
        hs = THStack("hs", "")
#ifndef VVBDTANALYSIS
        hs.Add(hVV)
#else
        #hs.Add(hVVLF)
        #hs.Add(hVVHF)  // FIXME: should be plotted as signal?
#endif
        hs.Add(hQCD)
        hs.Add(hs_Top)
        hs.Add(hTT)
        hs.Add(hWj0b)
        hs.Add(hWj1b)
        hs.Add(hWj2b)
        hs.Add(hZj0b)
        hs.Add(hZj1b)
        hs.Add(hZj2b)
        if plotSig:
#ifdef VVBDTANALYSIS
        #hs.Add(hVVHF)
#endif
            hs.Add(hVH)
       
        
        ymax = max(hdata_plot.GetMaximum(), hs.GetMaximum())
        hs.SetMaximum(ymax + ymax / 2.0 + (sqrt(ymax) if ymax>1 else 0))
        if plotLog:
            hs.SetMaximum(ymax * 200 + (sqrt(ymax) if ymax>1 else 0))
        hs.SetMinimum(0.01)
        
        # Setup auxiliary histograms
        print "MakePlots(): Setting up auxiliary histograms..."
        staterr = hmc_exp.Clone("staterr")
        staterr.Sumw2()
        staterr.SetFillColor(kRed)
        staterr.SetMarkerSize(0)
        staterr.SetFillStyle(3013)
        
        ratio = hdata_plot.Clone("ratio")
        ratio.Sumw2()
        ratio.SetMarkerSize(0.8)
        #ratio.SetMarkerSize(0.5)
        ratio.Divide(hdata_plot, hmc_exp, 1., 1., "")
        
        ratiostaterr = hmc_exp.Clone("ratiostaterr")
        ratiostaterr.Sumw2()
        ratiostaterr.SetStats(0)
        ratiostaterr.SetTitle("")
        ratiostaterr.GetXaxis().SetTitle(xtitle)
        ratiostaterr.GetYaxis().SetTitle("Data/MC")
        ratiostaterr.SetMaximum(2.2)
        ratiostaterr.SetMinimum(0)
        ratiostaterr.SetMarkerSize(0)
        ratiostaterr.SetFillColor(kRed)
        ratiostaterr.SetFillStyle(3013)
        #ratiostaterr.SetFillStyle(3001)
        ratiostaterr.GetXaxis().CenterTitle()
        ratiostaterr.GetXaxis().SetLabelSize(0.12)
        ratiostaterr.GetXaxis().SetTitleSize(0.14)
        ratiostaterr.GetXaxis().SetTitleOffset(1.10)
        ratiostaterr.GetYaxis().CenterTitle()
        ratiostaterr.GetYaxis().SetLabelSize(0.10)
        ratiostaterr.GetYaxis().SetTitleSize(0.12)
        #ratiostaterr.GetYaxis().SetTitleSize(0.10)
        ratiostaterr.GetYaxis().SetTitleOffset(0.6)
        ratiostaterr.GetYaxis().SetNdivisions(505)
        
        for i in xrange(0, hmc_exp.GetNbinsX()+2):
            ratiostaterr.SetBinContent(i, 1.0)
            if (hmc_exp.GetBinContent(i) > 1e-6):  # use smaller tolerance?
                binerror = hmc_exp.GetBinError(i) / hmc_exp.GetBinContent(i)
                ratiostaterr.SetBinError(i, binerror)
            else:
                ratiostaterr.SetBinError(i, 999.)
           
            
            #if (!(hdata_plot.GetBinContent(i) > 1e-6)):
            #    ratiostaterr.SetBinError(i, 0.)
           
       
        
        ratiosysterr = ratiostaterr.Clone("ratiosysterr")
        ratiosysterr.Sumw2()
        ratiosysterr.SetMarkerSize(0)
        #ratiosysterr.SetFillColor(kBlue)
        ratiosysterr.SetFillColor(kYellow-4)
        #ratiosysterr.SetFillStyle(3002)
        ratiosysterr.SetFillStyle(1001)
        
        for i in xrange(0, hmc_exp.GetNbinsX()+2):
            if (hmc_exp.GetBinContent(i) > 1e-6):  # use smaller tolerance?
#ifndef VVBDTANALYSIS
                binerror2        = (pow(hmc_exp.GetBinError(i), 2) +
                                    pow(0.05 * hWj0b.GetBinContent(i), 2) +
                                    pow(0.25 * hWj1b.GetBinContent(i), 2) +
                                    pow(0.20 * hWj2b.GetBinContent(i), 2) +
                                    pow(0.05 * hZj0b.GetBinContent(i), 2) +
                                    pow(0.15 * hZj1b.GetBinContent(i), 2) +
                                    pow(0.10 * hZj2b.GetBinContent(i), 2) +
                                    pow(0.04 * hTT.GetBinContent(i), 2) +
                                    pow(0.25 * hs_Top.GetBinContent(i), 2) +
#ifndef VVBDTANALYSIS
                                    pow(0.30 * hVV.GetBinContent(i), 2))
#else
                                    #pow(0.30 * hVVLF.GetBinContent(i), 2))
#endif
                binerror = sqrt(binerror2) / hmc_exp.GetBinContent(i)
                ratiosysterr.SetBinError(i, binerror)
            
        

        # Setup legends
        print "MakePlots(): Setting up legends..."
        leg = TLegend(0.72, 0.62, 0.92, 0.92)
        leg.SetFillColor(0)
        leg.SetLineColor(0)
        leg.SetShadowColor(0)
        leg.SetTextFont(62)
        leg.SetTextSize(0.03)
        leg.AddEntry(hdata_plot, "Data", "p")
        if plotSig:
            leg.AddEntry(hVH, "VH(%i)" % massH, "l")
        leg.AddEntry(hZj2b, "Z + b#bar{b}", "f")
        leg.AddEntry(hZj1b, "Z + b", "f")
        leg.AddEntry(hZj0b, "Z + udscg", "f")
        leg.AddEntry(hWj2b, "W + b#bar{b}", "f")
        leg.AddEntry(hWj1b, "W + b", "f")
        leg.AddEntry(hWj0b, "W + udscg", "f")
        leg.AddEntry(hTT, "t#bar{t}", "f")
        leg.AddEntry(hs_Top, "single top", "f")
        leg.AddEntry(hQCD, "QCD", "f")
#ifndef VVBDTANALYSIS
        leg.AddEntry(hVV, "VV", "f")
#else
        #leg.AddEntry(hVVHF, "VV(b#bar{b})", "f")
        #leg.AddEntry(hVVLF, "VV(udscg)", "f")
#endif
        leg.AddEntry(staterr, "MC uncert. (stat)", "fl")

        ratioleg1 = TLegend(0.54, 0.86, 0.72, 0.96)
        #ratioleg1 = TLegend(0.50, 0.86, 0.69, 0.96)
        ratioleg1.AddEntry(ratiostaterr, "MC uncert. (stat)", "f")
        ratioleg1.SetFillColor(0)
        ratioleg1.SetLineColor(0)
        ratioleg1.SetShadowColor(0)
        ratioleg1.SetTextFont(62)
        ratioleg1.SetTextSize(0.06)
        ratioleg1.SetBorderSize(1)
        
        ratioleg2 = TLegend(0.72, 0.86, 0.95, 0.96)
        #ratioleg2 = TLegend(0.69, 0.86, 0.9, 0.96)
        ratioleg2.AddEntry(ratiosysterr, "MC uncert. (stat+syst)", "f")
        ratioleg2.SetFillColor(0)
        ratioleg2.SetLineColor(0)
        ratioleg2.SetShadowColor(0)
        ratioleg2.SetTextFont(62)
        ratioleg2.SetTextSize(0.06)
        ratioleg2.SetBorderSize(1)

        # Draw MC signal and background
        print "MakePlots(): Drawing..."
        pad1.cd()
        if plotLog: pad1.SetLogy()
        hs.Draw("hist")
        hs.GetXaxis().SetLabelSize(0)
        binwidth = (xup - xlow) / nbins_plot
        ytitle = "Events / %.3f" % binwidth
        hs.GetYaxis().SetTitle(ytitle)
        
        staterr.Draw("e2 same")
        if plotSig:
            hVH.SetLineWidth(3)
            hVH.SetFillColor(0)
            hVH.Draw("hist same")
       
        
        # Draw data
        hdata_plot.Draw("e1 same")
        
        # Draw legends
        leg.Draw()
        latex = TLatex()
        latex.SetNDC()
        latex.SetTextAlign(12)
        latex.SetTextFont(62)
        latex.SetTextSize(0.052)
        latex.DrawLatex(0.19, 0.89, "CMS Preliminary")
        latex.SetTextSize(0.04)
        latex.DrawLatex(0.19, 0.84, "#sqrt{s} = 8 TeV, L = 19.0 fb^{-1}")
        latex.DrawLatex(0.19, 0.79, "Z(#nu#bar{#nu})H(b#bar{b})")
        
        # Under/overflows a la TMVA
        uoflow = "U/O-flow (Data,MC): (%.1f, %.1f) / (%.1f, %.1f)" % (hdata_plot.GetBinContent(0), hmc_exp.GetBinContent(0), hdata_plot.GetBinContent(nbins_plot+1), hmc_exp.GetBinContent(nbins_plot+1))
        latex2 = TLatex(0.99, 0.1, uoflow)
        latex2.SetNDC()
        latex2.SetTextSize(0.02)
        latex2.SetTextAngle(90)
        latex2.AppendPad()
        
        # Draw ratio
        pad2.cd()
        pad2.SetGridy(0)
        ratiostaterr.Draw("e2")
        ratiosysterr.Draw("e2 same")
        ratiostaterr.Draw("e2 same")
        ratio.Draw("e1 same")
        
        # Draw ratio legends
        ratioleg1.Draw()
        ratioleg2.Draw()

        # Kolmogorov-Smirnov test and Chi2 test
        pave = TPaveText(0.18, 0.85, 0.35, 0.96, "brNDC")
        pave.SetTextAlign(12)
        pave.SetLineColor(0)
        pave.SetFillColor(0)
        pave.SetShadowColor(0)
        pave.SetBorderSize(1)
        nchisq = hdata_plot.Chi2Test(hmc_test, "UWCHI2/NDF")  # MC uncert. (stat)
        kolprob = hdata_plot.KolmogorovTest(hmc_test)  # MC uncert. (stat)
        text = pave.AddText("#chi_{#nu}^{2} = %.3f, K_{s} = %.3f" % (nchisq, kolprob))
        text.SetTextFont(62)
        text.SetTextSize(0.06)
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
        plotname = TString("%s_%s_%s" % (channel, region, var))
        FormatFileName(plotname)
        gPad.Print(plotdir+plotname.Data()+".png")
        gPad.Print(plotdir+plotname.Data()+".pdf")
        
    return 0
    

## Main
## ----------------------------------------------------------------------------

if __name__ == "__main__":
    
    #cutallmc    = TCut("")
    #cutalldata  = TCut("")
    #cutmc       = TCut("weightsHCP[0] * (selectFlags[5][0])")
    #cutdata     = TCut("selectFlags[5][0]")
    #ev = Read("Step4", cutallmc, cutalldata)
    #MakePlots(ev, "hJet_ptReg[0]", cutmc, cutdata, "H j1 p_{T} [GeV]", 20, 50, 350)

    cutallmc    = TCut("(Vtype==2||Vtype==3) && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.898 && METtype1corr.et>130 && HptReg>130 && max(hJet_ptReg[0],hJet_ptReg[1])>60 && min(hJet_ptReg[0],hJet_ptReg[1])>30 && HmassReg<250 && mindPhiMETJet_dPhi>0.5 && naJets_Znn>=1 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && !isBadHcalEvent && (FatH.FatHiggsFlag==1 && nfathFilterJets>0)")
    cutalldata  = cutallmc
    cutmc       = TCut("efflumi * PUweight * triggercorr2012ABCD( (triggerFlags[42]==1 || triggerFlags[39]==1), (triggerFlags[41]==1), METtype1corr.et, max(hJet_csv_nominal[0],hJet_csv_nominal[1]) ) * 19040/19624")
    cutdata     = TCut("EVENT.json && ( (190456<=EVENT.run && EVENT.run<=193752 && (triggerFlags[42]==1 || triggerFlags[49]==1 || triggerFlags[40]==1)) || (193752<=EVENT.run && EVENT.run<=208686 && (triggerFlags[42]==1 || triggerFlags[39]==1 || triggerFlags[41]==1)) )")
    ev = Read("Step3", cutallmc, cutalldata)
    #MakePlots(ev, "nPVs", cutmc, cutdata, "# of primary vertices (UNUSED)", 20, 0, 40)
    #MakePlots(ev, "rho25", cutmc, cutdata, "energy density in |#eta|<2.5 [GeV] (UNUSED)", 20, 0, 30)
    #MakePlots(ev, "hJet_ptRaw[0]", cutmc, cutdata, "H j1 raw p_{T} [GeV] (UNUSED)", 20, 50, 350)
    #MakePlots(ev, "smear_pt_res(hJet_ptRaw[0], hJet_genPt[0], hJet_eta[0])", cutmc, cutdata, "H j1 raw (JER corr.) p_{T} [GeV]", 20, 50, 350)
    #MakePlots(ev, "hJet_pt[0]", cutmc, cutdata, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakePlots(ev, "evalEt(hJet_pt[0], hJet_eta[0], hJet_phi[0], hJet_e[0])", cutmc, cutdata, "H j1 E_{T} [GeV]", 20, 50, 350)
    #MakePlots(ev, "evalMt(hJet_pt[0], hJet_eta[0], hJet_phi[0], hJet_e[0])", cutmc, cutdata, "H j1 m_{T} [GeV]", 20, 50, 350)
    #MakePlots(ev, "hJet_eta[0]", cutmc, cutdata, "H j1 #eta (UNUSED)", 24, -3.0, 3.0)
    #MakePlots(ev, "hJet_ptLeadTrack[0]", cutmc, cutdata, "H j1 leading track p_{T} [GeV]", 20, 0, 150)
    #MakePlots(ev, "max(0,hJet_vtx3dL[0])", cutmc, cutdata, "H j1 sec vtx 3-d flight length [cm]", 20, 0, 8)
    #MakePlots(ev, "max(0,hJet_vtx3deL[0])", cutmc, cutdata, "H j1 sec vtx 3-d flight length error [cm]", 20, 0, 0.8)
    #MakePlots(ev, "max(0,hJet_vtxMass[0])", cutmc, cutdata, "H j1 sec vtx mass [GeV]", 20, 0, 8)
    #MakePlots(ev, "max(0,hJet_vtxPt[0])", cutmc, cutdata, "H j1 sec vtx p_{T} [GeV]", 20, 0, 200)
    #MakePlots(ev, "hJet_chf[0]", cutmc, cutdata, "H j1 charged hadronic frac (UNUSED)", 21, 0, 1.05)
    #MakePlots(ev, "hJet_cef[0]", cutmc, cutdata, "H j1 charged electromagnetic frac", 21, 0, 1.05)
    #MakePlots(ev, "hJet_nhf[0]", cutmc, cutdata, "H j1 neutral hadronic frac (UNUSED)", 21, 0, 1.05)
    #MakePlots(ev, "hJet_nef[0]", cutmc, cutdata, "H j1 neutral electromagnetic frac(UNUSED)", 21, 0, 1.05)
    #MakePlots(ev, "hJet_nch[0]", cutmc, cutdata, "H j1 # charged particles (UNUSED)", 20, 0, 100)
    #MakePlots(ev, "hJet_nconstituents[0]", cutmc, cutdata, "H j1 # constituents [GeV]", 24, 0, 120)
    #MakePlots(ev, "hJet_JECUnc[0]", cutmc, cutdata, "H j1 JEC uncertainty [GeV]", 25, 0.005, 0.025)
    #MakePlots(ev, "max(0,hJet_SoftLeptptRel[0]*(hJet_SoftLeptIdlooseMu[0]==1 || hJet_SoftLeptId95[0]==1))", cutmc, cutdata, "H j1 soft lepton p_{T, rel} [GeV]", 20, 0, 80)
    #MakePlots(ev, "max(0,hJet_SoftLeptPt[0]*(hJet_SoftLeptIdlooseMu[0]==1 || hJet_SoftLeptId95[0]==1))", cutmc, cutdata, "H j1 soft lepton p_{T} [GeV]", 20, 0, 80)
    #MakePlots(ev, "max(0,hJet_SoftLeptdR[0]*(hJet_SoftLeptIdlooseMu[0]==1 || hJet_SoftLeptId95[0]==1))", cutmc, cutdata, "H j1 soft lepton #DeltaR", 12, 0, 0.6)
    #MakePlots(ev, "hJet_csv_nominal[0]", cutmc, cutdata, "H j1 CSV", 30, 0.0, 1.08)
    #MakePlots(ev, "hJet_ptReg[0]", cutmc, cutdata, "H j1 regressed p_{T} [GeV]", 20, 50, 350)
    
    
    MakePlots(ev, "fathFilterJets_ptRaw[0]", cutmc, cutdata, "H fj1 raw p_{T} [GeV] (UNUSED)", 20, 50, 350)
    MakePlots(ev, "smear_pt_res(fathFilterJets_ptRaw[0], fathFilterJets_genPt[0], fathFilterJets_eta[0])", cutmc, cutdata, "H fj1 raw (JER corr.) p_{T} [GeV]", 20, 50, 350)
    MakePlots(ev, "fathFilterJets_pt[0]", cutmc, cutdata, "H fj1 p_{T} [GeV]", 20, 50, 350)
    MakePlots(ev, "evalEt(fathFilterJets_pt[0], fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_e[0])", cutmc, cutdata, "H fj1 E_{T} [GeV]", 20, 50, 350)
    MakePlots(ev, "evalMt(fathFilterJets_pt[0], fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_e[0])", cutmc, cutdata, "H fj1 m_{T} [GeV]", 20, 50, 350)
    MakePlots(ev, "max(0,fathFilterJets_ptLeadTrack[0])", cutmc, cutdata, "H fj1 leading track p_{T} [GeV]", 20, 0, 150)
    MakePlots(ev, "max(0,fathFilterJets_vtx3dL[0])", cutmc, cutdata, "H fj1 sec vtx 3-d flight length [cm]", 20, 0, 8)
    MakePlots(ev, "max(0,fathFilterJets_vtx3deL[0])", cutmc, cutdata, "H fj1 sec vtx 3-d flight length error [cm]", 20, 0, 0.8)
    MakePlots(ev, "max(0,fathFilterJets_vtxMass[0])", cutmc, cutdata, "H fj1 sec vtx mass [GeV]", 20, 0, 8)
    MakePlots(ev, "max(0,fathFilterJets_vtxPt[0])", cutmc, cutdata, "H fj1 sec vtx p_{T} [GeV]", 20, 0, 200)
    MakePlots(ev, "fathFilterJets_cef[0]", cutmc, cutdata, "H fj1 charged electromagnetic frac", 21, 0, 1.05)
    MakePlots(ev, "fathFilterJets_nconstituents[0]", cutmc, cutdata, "H fj1 # constituents [GeV]", 24, 0, 120)
    MakePlots(ev, "fathFilterJets_JECUnc[0]", cutmc, cutdata, "H fj1 JEC uncertainty [GeV]", 25, 0.005, 0.025)
    MakePlots(ev, "fathFilterJets_jetArea[0]", cutmc, cutdata, "H fj1 area", 20, 0.0, 1.0)
    MakePlots(ev, "fathFilterJets_SoftLeptptRel[0]", cutmc, cutdata, "H fj1 soft lepton p_{T, rel} [GeV]", 20, 0, 80)
    MakePlots(ev, "fathFilterJets_SoftLeptPt[0]", cutmc, cutdata, "H fj1 soft lepton p_{T} [GeV]", 20, 0, 80)
    MakePlots(ev, "fathFilterJets_SoftLeptdR[0]", cutmc, cutdata, "H fj1 soft lepton #DeltaR", 12, 0, 0.6)
    MakePlots(ev, "fathFilterJets_csv_nominal[0]", cutmc, cutdata, "H fj1 CSV", 30, 0.0, 1.08)
    MakePlots(ev, "fathFilterJets_ptReg[0]", cutmc, cutdata, "H fj1 regressed p_{T} [GeV]", 20, 50, 350)
    MakePlots(ev, "FatH.filteredmass", cutmc, cutdata, "H_{fj} mass [GeV]", 17, 0, 255)
    MakePlots(ev, "FatHmassReg", cutmc, cutdata, "H_{fj} regressed mass [GeV]", 17, 0, 255)
    MakePlots(ev, "FatH.filteredpt", cutmc, cutdata, "H_{fj} p_{T} [GeV]", 20, 50, 350)
    MakePlots(ev, "FatHptReg", cutmc, cutdata, "H_{fj} regressed p_{T} [GeV]", 20, 50, 350)
