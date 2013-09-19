from ROOT import TCanvas, TCut, TChain, TH1, TH1F, THStack, TLatex, TLegend, TPad, TPaveText, TString, gROOT, gPad, gStyle, gSystem, TProfile
gROOT.ProcessLine(".L tdrstyle.C")
gROOT.ProcessLine("setTDRStyle()")
gROOT.ProcessLine(".L HelperBDTShape.h")
gROOT.ProcessLine(".L HelperFunctions.h")
from ROOT import setHisto, FormatFileName
from math import sqrt


channel  = "ZnunuHighPt"
region   = "Data"
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
        indir1   = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130314/stitch/"
        indir2   = indir1
        indir3   = indir1
        prefix   = "Step4_"
        suffix   = ".root"
        treename = "tree_%s_%s" % (channel, "ctrl")
        
        P1 = TChain(treename)
        P1.Add(indir1 + prefix + "TT" + suffix)
        
        P2 = TChain(treename)
        P2.Add(indir2 + prefix + "TT" + suffix)
        
        P3 = TChain(treename)
        P3.Add(indir3 + prefix + "TT" + suffix)

        cutHF = TCut("eventFlav==5")
        cutLF = TCut("eventFlav!=5")
        
        cut2b = TCut("abs(hJet_flavour[0])==5 && abs(hJet_flavour[1])==5")
        cut1b = TCut("(abs(hJet_flavour[0])==5 && abs(hJet_flavour[1])!=5) || (abs(hJet_flavour[0])!=5 && abs(hJet_flavour[1])==5)")
        cut0b = TCut("abs(hJet_flavour[0])!=5 && abs(hJet_flavour[1])!=5")
    
    elif whatstep == "Step3":
        indir1   = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/skim_ZnnH_Step3/"
        indir2   = indir1
        indir3   = indir1
        prefix   = "Step3_"
        suffix   = ".root"
        treename = "tree"
        
        P1 = TChain(treename)
        P1.Add(indir1 + prefix + "data_obs" + suffix)
        #names = ["Wj", "TT", "s_Top", "VV"]
        #for n in names:
        #    P1.Add(indir1 + prefix + n + suffix)

        P2 = TChain(treename)
        P2.Add(indir2 + prefix + "data_obs" + suffix)
        #names = ["Wj", "TT", "s_Top", "VV"]
        #for n in names:
        #    P2.Add(indir2 + prefix + n + suffix)
        
        P3 = TChain(treename)
        P3.Add(indir3 + prefix + "data_obs" + suffix)
        #names = ["Wj", "TT", "s_Top", "VV"]
        #for n in names:
        #    P3.Add(indir3 + prefix + n + suffix)

        cutHF = TCut("eventFlav==5")
        cutLF = TCut("eventFlav!=5")
        
        cut2b = TCut("abs(hJet_flavour[0])==5 && abs(hJet_flavour[1])==5")
        cut1b = TCut("(abs(hJet_flavour[0])==5 && abs(hJet_flavour[1])!=5) || (abs(hJet_flavour[0])!=5 && abs(hJet_flavour[1])==5)")
        cut0b = TCut("abs(hJet_flavour[0])!=5 && abs(hJet_flavour[1])!=5")
        
        cutABC  = TCut("EVENT.run<=203002")
        cutD1   = TCut("203002<EVENT.run && (!(207883<=EVENT.run && EVENT.run<=208307))")
        cutD2   = TCut("203002<EVENT.run && (207883<=EVENT.run && EVENT.run<=208307)")
        
    ev = EventsJ11()
    #ev.P1 = P1.CopyTree(cutallmc.GetTitle())
    #ev.P2 = P2.CopyTree(cutallmc.GetTitle())
    #ev.P3 = P3.CopyTree(cutallmc.GetTitle())
    ev.P1 = P1.CopyTree((cutalldata + cutABC).GetTitle())
    ev.P2 = P2.CopyTree((cutalldata + cutD1).GetTitle())
    ev.P3 = P3.CopyTree((cutalldata + cutD2).GetTitle())
    
    return ev


def MakePlots(ev, var, cutmc, cutdata, xtitle, nbins, xlow, xup):
    hP1        = TH1F("P1"       , "", nbins, xlow, xup)
    hP2        = TH1F("P2"       , "", nbins, xlow, xup)
    hP3        = TH1F("P3"       , "", nbins, xlow, xup)

    # Apply scale factors
    scalefactor1 = 1.0
    cutsf = TCut("%f" % scalefactor1)
    #ev.P1.Project("P1", var, (cutmc * cutsf).GetTitle())
    ev.P1.Project("P1", var, (cutdata).GetTitle())
    scalefactor2 = 1.0
    cutsf = TCut("%f" % scalefactor2)
    #ev.P2.Project("P2", var, (cutmc * cutsf).GetTitle())
    ev.P2.Project("P2", var, (cutdata).GetTitle())
    scalefactor3 = 1.0
    cutsf = TCut("%f" % scalefactor3)
    #ev.P3.Project("P3", var, (cutmc * cutsf).GetTitle())
    ev.P3.Project("P3", var, (cutdata).GetTitle())
    
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
        setHisto(hP1, "data_obs")
        setHisto(hP2, "VH")
        setHisto(hP3, "TT")

        hP1.Scale(1.0 / hP1.GetSumOfWeights())
        hP2.Scale(1.0 / hP2.GetSumOfWeights())
        hP3.Scale(1.0 / hP3.GetSumOfWeights())
        
        ymax = max(hP1.GetMaximum(), hP2.GetMaximum(), hP3.GetMaximum())
        hP1.SetMaximum(ymax + ymax / 2.0 + (sqrt(ymax) if ymax>1 else 0))
        if plotLog:
            hP1.SetMaximum(ymax * 200 + (sqrt(ymax) if ymax>1 else 0))
        #hP1.SetMinimum(0.01)
        
        # Setup auxiliary histograms
        print "MakePlots(): Setting up auxiliary histograms..."
        
        ratio = hP2.Clone("ratio")
        ratio.Sumw2()
        ratio.SetMarkerSize(0.8)
        #ratio.SetMarkerSize(0.5)
        ratio.Divide(hP2, hP1, 1., 1., "")
        
        ratioP3 = hP3.Clone("ratio3")
        ratioP3.Sumw2()
        ratioP3.SetMarkerSize(0.8)
        #ratioP3.SetMarkerSize(0.5)
        ratioP3.Divide(hP3, hP1, 1., 1., "")
        
        ratio.SetStats(0)
        ratio.SetTitle("")
        ratio.GetXaxis().SetTitle(xtitle)
        ratio.GetYaxis().SetTitle("ratio")
        ratio.SetMaximum(2.2)
        ratio.SetMinimum(0)
        #ratio.SetMarkerSize(0)
        #ratio.SetFillColor(kRed)
        #ratio.SetFillStyle(3013)
        ##ratio.SetFillStyle(3001)
        ratio.GetXaxis().CenterTitle()
        ratio.GetXaxis().SetLabelSize(0.12)
        ratio.GetXaxis().SetTitleSize(0.14)
        ratio.GetXaxis().SetTitleOffset(1.10)
        ratio.GetYaxis().CenterTitle()
        ratio.GetYaxis().SetLabelSize(0.10)
        ratio.GetYaxis().SetTitleSize(0.12)
        #ratio.GetYaxis().SetTitleSize(0.10)
        ratio.GetYaxis().SetTitleOffset(0.6)
        ratio.GetYaxis().SetNdivisions(505)

        # Setup legends
        print "MakePlots(): Setting up legends..."
        leg = TLegend(0.62, 0.62, 0.92, 0.92)
        leg.SetFillColor(0)
        leg.SetLineColor(0)
        leg.SetShadowColor(0)
        leg.SetTextFont(62)
        leg.SetTextSize(0.03)
        leg.AddEntry(hP1, "2012ABC", "lp")
        leg.AddEntry(hP2, "2012D ![207883:208307]", "lp")
        leg.AddEntry(hP3, "2012D [207883:208307]", "lp")
        
        # Draw MC signal and background
        print "MakePlots(): Drawing..."
        pad1.cd()
        if plotLog: pad1.SetLogy()
        hP1.Draw("e1")
        hP1.GetXaxis().SetLabelSize(0)
        binwidth = (xup - xlow) / nbins
        ytitle = "(normalized) / %.3f" % binwidth
        hP1.GetYaxis().SetTitle(ytitle)
        
        hP3.Draw("e1 same")
        hP2.Draw("e1 same")
        
        # Draw legends
        leg.Draw()
        latex = TLatex()
        latex.SetNDC()
        latex.SetTextAlign(12)
        latex.SetTextFont(62)
        latex.SetTextSize(0.052)
        latex.DrawLatex(0.19, 0.89, "CMS Preliminary")
        latex.SetTextSize(0.04)
        latex.DrawLatex(0.19, 0.84, "#sqrt{s} = 8 TeV, L = 19.6 fb^{-1}")
        latex.DrawLatex(0.19, 0.79, "Z(#nu#bar{#nu})H(b#bar{b})")
        
        # Under/overflows a la TMVA
        uoflow = "U/O-flow (Data,MC): (%.1f, %.1f) / (%.1f, %.1f)" % (hP2.GetBinContent(0), hP1.GetBinContent(0), hP2.GetBinContent(nbins+1), hP1.GetBinContent(nbins+1))
        latex2 = TLatex(0.99, 0.1, uoflow)
        latex2.SetNDC()
        latex2.SetTextSize(0.02)
        latex2.SetTextAngle(90)
        latex2.AppendPad()
        
        # Draw ratio
        pad2.cd()
        pad2.SetGridy(0)
        ratio.Draw("e1")
        
        ratioP3.Draw("e1 same")
        ratio.Draw("e1 same")
        
        # Kolmogorov-Smirnov test and Chi2 test
        pave = TPaveText(0.18, 0.85, 0.35, 0.96, "brNDC")
        pave.SetTextAlign(12)
        pave.SetLineColor(0)
        pave.SetFillColor(0)
        pave.SetShadowColor(0)
        pave.SetBorderSize(1)
        nchisq = hP2.Chi2Test(hP1, "UWCHI2/NDF")  # MC uncert. (stat)
        kolprob = hP2.KolmogorovTest(hP1)  # MC uncert. (stat)
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
        plotname = TString("%s_%s_diff_%s" % (channel, region, var))
        FormatFileName(plotname)
        gPad.Print(plotdir+plotname.Data()+".png")
        gPad.Print(plotdir+plotname.Data()+".pdf")

    return 0


def MakeProfiles(ev, var, cutmc, cutdata, xtitle, nbins, xlow, xup, ytitle, nbinsy, ylow, yup):
    hP1        = TProfile("P1"       , "", nbins, xlow, xup, ylow, yup, "")  # Data
    hP2        = TProfile("P2"       , "", nbins, xlow, xup, ylow, yup, "")  # MC
    hP3        = TProfile("P3"       , "", nbins, xlow, xup, ylow, yup, "")

    # Apply scale factors
    scalefactor1 = 1.0
    cutsf = TCut("%f" % scalefactor1)
    #ev.P1.Project("P1", var, (cutmc * cutsf).GetTitle(), "prof")
    ev.P1.Project("P1", var, (cutdata).GetTitle(), "prof")
    scalefactor2 = 1.0
    cutsf = TCut("%f" % scalefactor2)
    ev.P2.Project("P2", var, (cutmc * cutsf).GetTitle(), "prof")
    #ev.P2.Project("P2", var, (cutdata).GetTitle(), "prof")
    scalefactor3 = 1.0
    cutsf = TCut("%f" % scalefactor3)
    #ev.P3.Project("P3", var, (cutmc * cutsf).GetTitle(), "prof")
    #ev.P3.Project("P3", var, (cutdata).GetTitle(), "prof")
    
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

        hP1.SetMaximum(yup)
        hP1.SetMinimum(ylow)

        # Setup histogram styles and stack the histograms.
        hP2.SetFillColor(kBlue)
        hP2.SetMarkerSize(0)
        hP2.SetFillStyle(3013)
        
        # Setup auxiliary histograms
        print "MakePlots(): Setting up auxiliary histograms..."
        
        hP1p = hP1.ProjectionX()
        hP2p = hP2.ProjectionX()
        ratio = hP1p.Clone("ratio")
        ratio.Sumw2()
        ratio.SetMarkerSize(0.8)
        #ratio.SetMarkerSize(0.5)
        ratio.Divide(hP1p, hP2p, 1., 1., "")
        
        ratiostaterr = hP2p.Clone("ratiostaterr")
        ratiostaterr.Sumw2()
        ratiostaterr.SetStats(0)
        ratiostaterr.SetTitle("")
        ratiostaterr.GetXaxis().SetTitle(xtitle)
        ratiostaterr.GetYaxis().SetTitle("Data/MC")
        ratiostaterr.SetMaximum(2.2)
        ratiostaterr.SetMinimum(0)
        ratiostaterr.SetMarkerSize(0)
        ratiostaterr.SetFillColor(kBlue)
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
        
        for i in xrange(0, hP2p.GetNbinsX()+2):
            ratiostaterr.SetBinContent(i, 1.0)
            if (hP2p.GetBinContent(i) > 1e-6):  # use smaller tolerance?
                binerror = hP2p.GetBinError(i) / hP2p.GetBinContent(i)
                ratiostaterr.SetBinError(i, binerror)
            else:
                ratiostaterr.SetBinError(i, 999.)

        # Setup legends
        print "MakePlots(): Setting up legends..."
        leg = TLegend(0.72, 0.72, 0.92, 0.92)
        leg.SetFillColor(0)
        leg.SetLineColor(0)
        leg.SetShadowColor(0)
        leg.SetTextFont(62)
        leg.SetTextSize(0.03)
        leg.AddEntry(hP1, "Data", "p")
        leg.AddEntry(hP2, "MC", "f")
        
        ratioleg1 = TLegend(0.72, 0.86, 0.95, 0.96)
        #ratioleg1 = TLegend(0.69, 0.86, 0.9, 0.96)
        ratioleg1.AddEntry(ratiostaterr, "MC uncert. (stat)", "f")
        ratioleg1.SetFillColor(0)
        ratioleg1.SetLineColor(0)
        ratioleg1.SetShadowColor(0)
        ratioleg1.SetTextFont(62)
        ratioleg1.SetTextSize(0.06)
        ratioleg1.SetBorderSize(1)

        # Draw MC signal and background
        print "MakePlots(): Drawing..."
        pad1.cd()
        if plotLog: pad1.SetLogy()
        hP1.Draw("e1")
        hP1.GetXaxis().SetLabelSize(0)
        binwidth = (xup - xlow) / nbins
        hP1.GetYaxis().SetTitle(ytitle)
        
        hP2.Draw("e2 same")
        hP1.Draw("e1 same")
        
        # Draw legends
        leg.Draw()
        latex = TLatex()
        latex.SetNDC()
        latex.SetTextAlign(12)
        latex.SetTextFont(62)
        latex.SetTextSize(0.052)
        latex.DrawLatex(0.19, 0.89, "CMS Preliminary")
        latex.SetTextSize(0.04)
        latex.DrawLatex(0.19, 0.84, "#sqrt{s} = 8 TeV, L = 19.6 fb^{-1}")
        latex.DrawLatex(0.19, 0.79, "Z(#nu#bar{#nu})H(b#bar{b})")
        
        # Draw ratio
        pad2.cd()
        pad2.SetGridy(0)
        ratiostaterr.Draw("e2")
        ratio.Draw("e1 same")
        
        # Draw ratio legends
        ratioleg1.Draw()

        # Kolmogorov-Smirnov test and Chi2 test
        pave = TPaveText(0.18, 0.85, 0.35, 0.96, "brNDC")
        pave.SetTextAlign(12)
        pave.SetLineColor(0)
        pave.SetFillColor(0)
        pave.SetShadowColor(0)
        pave.SetBorderSize(1)
        nchisq = hP2p.Chi2Test(hP1p, "WWCHI2/NDF")  # MC uncert. (stat)
        kolprob = hP2p.KolmogorovTest(hP1p)  # MC uncert. (stat)
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
        plotname = TString("%s_%s_diff_%s" % (channel, region, var))
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

    cutallmc    = TCut("max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.898 && METtype1corr.et>130 && H.pt>130 && naJets_Znn>=1")
    cutalldata  = TCut("max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.898 && METtype1corr.et>130 && H.pt>130 && naJets_Znn>=1")
    #cutallmc    = TCut("max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.244 && naJets_Znn==0")
    #cutalldata  = TCut("max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.244 && naJets_Znn==0")
    cutmc       = TCut("efflumi * PUweight * (triggerFlags[42]==1 || triggerFlags[39]==1 || triggerFlags[41]==1)")
    cutdata     = TCut("EVENT.json && ( (190456<=EVENT.run && EVENT.run<=193752 && (triggerFlags[42]==1 || triggerFlags[49]==1 || triggerFlags[40]==1)) || (193752<=EVENT.run && EVENT.run<=208686 && (triggerFlags[42]==1 || triggerFlags[39]==1 || triggerFlags[41]==1)) )")
    ev = Read("Step3", cutallmc, cutalldata)
    MakePlots(ev, "nPVs", cutmc, cutdata, "# of primary vertices (UNUSED)", 20, 0, 40)
    MakePlots(ev, "rho25", cutmc, cutdata, "energy density in |#eta|<2.5 [GeV] (UNUSED)", 20, 0, 30)
    MakePlots(ev, "hJet_ptRaw[0]", cutmc, cutdata, "H j1 raw p_{T} [GeV] (UNUSED)", 20, 50, 350)
    MakePlots(ev, "smear_pt_res(hJet_ptRaw[0], hJet_genPt[0], hJet_eta[0])", cutmc, cutdata, "H j1 raw (JER corr.) p_{T} [GeV]", 20, 50, 350)
    MakePlots(ev, "hJet_pt[0]", cutmc, cutdata, "H j1 p_{T} [GeV]", 20, 50, 350)
    MakePlots(ev, "evalEt(hJet_pt[0], hJet_eta[0], hJet_phi[0], hJet_e[0])", cutmc, cutdata, "H j1 E_{T} [GeV]", 20, 50, 350)
    MakePlots(ev, "evalMt(hJet_pt[0], hJet_eta[0], hJet_phi[0], hJet_e[0])", cutmc, cutdata, "H j1 m_{T} [GeV]", 20, 50, 350)
    MakePlots(ev, "hJet_eta[0]", cutmc, cutdata, "H j1 #eta (UNUSED)", 24, -3.0, 3.0)
    MakePlots(ev, "hJet_ptLeadTrack[0]", cutmc, cutdata, "H j1 leading track p_{T} [GeV]", 20, 0, 150)
    MakePlots(ev, "max(0,hJet_vtx3dL[0])", cutmc, cutdata, "H j1 sec vtx 3-d flight length [cm]", 20, 0, 8)
    MakePlots(ev, "max(0,hJet_vtx3deL[0])", cutmc, cutdata, "H j1 sec vtx 3-d flight length error [cm]", 20, 0, 0.8)
    MakePlots(ev, "max(0,hJet_vtxMass[0])", cutmc, cutdata, "H j1 sec vtx mass [GeV]", 20, 0, 8)
    MakePlots(ev, "max(0,hJet_vtxPt[0])", cutmc, cutdata, "H j1 sec vtx p_{T} [GeV]", 20, 0, 200)
    MakePlots(ev, "hJet_chf[0]", cutmc, cutdata, "H j1 charged hadronic frac (UNUSED)", 21, 0, 1.05)
    MakePlots(ev, "hJet_cef[0]", cutmc, cutdata, "H j1 charged electromagnetic frac", 21, 0, 1.05)
    MakePlots(ev, "hJet_nhf[0]", cutmc, cutdata, "H j1 neutral hadronic frac (UNUSED)", 21, 0, 1.05)
    MakePlots(ev, "hJet_nef[0]", cutmc, cutdata, "H j1 neutral electromagnetic frac(UNUSED)", 21, 0, 1.05)
    MakePlots(ev, "hJet_nch[0]", cutmc, cutdata, "H j1 # charged particles (UNUSED)", 20, 0, 100)
    MakePlots(ev, "hJet_nconstituents[0]", cutmc, cutdata, "H j1 # constituents [GeV]", 24, 0, 120)
    MakePlots(ev, "hJet_JECUnc[0]", cutmc, cutdata, "H j1 JEC uncertainty [GeV]", 25, 0.005, 0.025)
    MakePlots(ev, "max(0,hJet_SoftLeptptRel[0]*(hJet_SoftLeptIdlooseMu[0]==1 || hJet_SoftLeptId95[0]==1))", cutmc, cutdata, "H j1 soft lepton p_{T, rel} [GeV]", 20, 0, 80)
    MakePlots(ev, "max(0,hJet_SoftLeptPt[0]*(hJet_SoftLeptIdlooseMu[0]==1 || hJet_SoftLeptId95[0]==1))", cutmc, cutdata, "H j1 soft lepton p_{T} [GeV]", 20, 0, 80)
    MakePlots(ev, "max(0,hJet_SoftLeptdR[0]*(hJet_SoftLeptIdlooseMu[0]==1 || hJet_SoftLeptId95[0]==1))", cutmc, cutdata, "H j1 soft lepton #DeltaR", 12, 0, 0.6)
    MakePlots(ev, "hJet_csv_nominal[0]", cutmc, cutdata, "H j1 CSV (UNUSED)", 30, 0.0, 1.08)
    MakePlots(ev, "hJet_ptReg[0]", cutmc, cutdata, "H j1 regressed p_{T} [GeV]", 20, 50, 350)
    
    
    #MakeProfiles(ev, "hJet_pt[0]:nPVs", cutmc, cutdata, "# of primary vertices (UNUSED)", 20, 0, 40, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakeProfiles(ev, "hJet_pt[0]:rho25", cutmc, cutdata, "energy density in |#eta|<2.5 [GeV] (UNUSED)", 20, 0, 30, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakeProfiles(ev, "hJet_pt[0]:hJet_ptRaw[0]", cutmc, cutdata, "H j1 raw p_{T} [GeV] (UNUSED)", 20, 50, 350, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakeProfiles(ev, "hJet_pt[0]:smear_pt_res(hJet_ptRaw[0], hJet_genPt[0], hJet_eta[0])", cutmc, cutdata, "H j1 raw (JER corr.) p_{T} [GeV]", 20, 50, 350, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakeProfiles(ev, "hJet_pt[0]:hJet_pt[0]", cutmc, cutdata, "H j1 p_{T} [GeV]", 20, 50, 350, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakeProfiles(ev, "hJet_pt[0]:evalEt(hJet_pt[0], hJet_eta[0], hJet_phi[0], hJet_e[0])", cutmc, cutdata, "H j1 E_{T} [GeV]", 20, 50, 350, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakeProfiles(ev, "hJet_pt[0]:evalMt(hJet_pt[0], hJet_eta[0], hJet_phi[0], hJet_e[0])", cutmc, cutdata, "H j1 m_{T} [GeV]", 20, 50, 350, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakeProfiles(ev, "hJet_pt[0]:hJet_eta[0]", cutmc, cutdata, "H j1 #eta (UNUSED)", 24, -3.0, 3.0, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakeProfiles(ev, "hJet_pt[0]:hJet_ptLeadTrack[0]", cutmc, cutdata, "H j1 leading track p_{T} [GeV]", 20, 0, 150, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakeProfiles(ev, "hJet_pt[0]:max(0,hJet_vtx3dL[0])", cutmc, cutdata, "H j1 sec vtx 3-d flight length [cm]", 20, 0, 8, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakeProfiles(ev, "hJet_pt[0]:max(0,hJet_vtx3deL[0])", cutmc, cutdata, "H j1 sec vtx 3-d flight length error [cm]", 20, 0, 0.8, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakeProfiles(ev, "hJet_pt[0]:max(0,hJet_vtxMass[0])", cutmc, cutdata, "H j1 sec vtx mass [GeV]", 20, 0, 8, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakeProfiles(ev, "hJet_pt[0]:max(0,hJet_vtxPt[0])", cutmc, cutdata, "H j1 sec vtx p_{T} [GeV]", 20, 0, 200, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakeProfiles(ev, "hJet_pt[0]:hJet_chf[0]", cutmc, cutdata, "H j1 charged hadronic frac (UNUSED)", 21, 0, 1.05, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakeProfiles(ev, "hJet_pt[0]:hJet_cef[0]", cutmc, cutdata, "H j1 charged electromagnetic frac", 21, 0, 1.05, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakeProfiles(ev, "hJet_pt[0]:hJet_nhf[0]", cutmc, cutdata, "H j1 neutral hadronic frac (UNUSED)", 21, 0, 1.05, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakeProfiles(ev, "hJet_pt[0]:hJet_nef[0]", cutmc, cutdata, "H j1 neutral electromagnetic frac(UNUSED)", 21, 0, 1.05, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakeProfiles(ev, "hJet_pt[0]:hJet_nch[0]", cutmc, cutdata, "H j1 # charged particles (UNUSED)", 20, 0, 100, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakeProfiles(ev, "hJet_pt[0]:hJet_nconstituents[0]", cutmc, cutdata, "H j1 # constituents [GeV]", 24, 0, 120, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakeProfiles(ev, "hJet_pt[0]:hJet_JECUnc[0]", cutmc, cutdata, "H j1 JEC uncertainty [GeV]", 25, 0.005, 0.025, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakeProfiles(ev, "hJet_pt[0]:max(0,hJet_SoftLeptptRel[0]*(hJet_SoftLeptIdlooseMu[0]==1 || hJet_SoftLeptId95[0]==1))", cutmc, cutdata, "H j1 soft lepton p_{T, rel} [GeV]", 20, 0, 80, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakeProfiles(ev, "hJet_pt[0]:max(0,hJet_SoftLeptPt[0]*(hJet_SoftLeptIdlooseMu[0]==1 || hJet_SoftLeptId95[0]==1))", cutmc, cutdata, "H j1 soft lepton p_{T} [GeV]", 20, 0, 80, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakeProfiles(ev, "hJet_pt[0]:max(0,hJet_SoftLeptdR[0]*(hJet_SoftLeptIdlooseMu[0]==1 || hJet_SoftLeptId95[0]==1))", cutmc, cutdata, "H j1 soft lepton #DeltaR", 12, 0, 0.6, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakeProfiles(ev, "hJet_pt[0]:hJet_csv_nominal[0]", cutmc, cutdata, "H j1 CSV (UNUSED)", 30, 0.0, 1.08, "H j1 p_{T} [GeV]", 20, 50, 350)
    #MakeProfiles(ev, "hJet_pt[0]:hJet_ptReg[0]", cutmc, cutdata, "H j1 regressed p_{T} [GeV]", 20, 50, 350, "H j1 p_{T} [GeV]", 20, 50, 350)
