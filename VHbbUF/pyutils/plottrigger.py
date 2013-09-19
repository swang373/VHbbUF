from ROOT import TCanvas, TCut, TChain, TH1, TH1F, THStack, TLatex, TLegend, TPad, TPaveText, TString, TLine, gROOT, gPad, gStyle, gSystem, TEfficiency, TGraphAsymmErrors, TColor
gROOT.ProcessLine(".L tdrstyle.C")
gROOT.ProcessLine("setTDRStyle()")
gROOT.ProcessLine(".L HelperBDTShape.h")
gROOT.ProcessLine(".L HelperFunctions.h")
from ROOT import setHisto, FormatFileName
from math import sqrt


channel  = "ZnunuHighPt"
#region   = "TT"
region   = "WjHF"
#region   = "WjLF"
whatdata = "SingleMu"
#whatdata = "SingleEl"
whatdata = "DijetMET"
massH    = 125
plotData = True
plotLog  = False
plotSig  = True
plotdir  = "plots_etc/"
kRed     = 632
kYellow  = 400
kBlue    = 600

colors = {
    "MET150"      : TColor.GetColor(193,39,45),
    "DijetMET100" : TColor.GetColor(253,185,19), 
    "MET80CSV"    : TColor.GetColor(0,106,68),
    "combo"       : TColor.GetColor(0,173,239),
    }

triggers = {
    # RunA, RunB+C+D
    "MET150"      : (42, 42), 
    "DijetMET100" : (49, 39),
    "MET80CSV"    : (40, 41),
    }

TH1.SetDefaultSumw2(1)
gROOT.SetBatch(1)
if gSystem.AccessPathName(plotdir):
    gSystem.mkdir(plotdir)

class EventsJ11(object):
    pass


def Read(whatdata="SingleMu", cutallmc="", cutalldata=""):
    
    if whatdata:
        indir    = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/skim_ZnnH_Step3/"
        prefix   = "Step3_"
        suffix   = ".root"
        treename = "tree"
        
        Data = TChain(treename)
        names = ["trig_%s" % whatdata]
        for n in names:
            Data.Add(indir + prefix + n + suffix)

        MC = TChain(treename)
        names = ["Wj", "TT", "s_Top", "VV"]  # FIXME: right now assume homogeneous background
        #names = ["TT"]
        for n in names:
            MC.Add(indir + prefix + n + suffix)

        cutHF = TCut("eventFlav==5")
        cutLF = TCut("eventFlav!=5")
        
        cut2b = TCut("abs(hJet_flavour[0])==5 && abs(hJet_flavour[1])==5")
        cut1b = TCut("(abs(hJet_flavour[0])==5 && abs(hJet_flavour[1])!=5) || (abs(hJet_flavour[0])!=5 && abs(hJet_flavour[1])==5)")
        cut0b = TCut("abs(hJet_flavour[0])!=5 && abs(hJet_flavour[1])!=5")
        
        cutABC  = TCut("EVENT.run<=203002")
        cutD1   = TCut("203002<EVENT.run && (!(207883<=EVENT.run && EVENT.run<=208307))")
        cutD2   = TCut("203002<EVENT.run && (207883<=EVENT.run && EVENT.run<=208307)")
    
    ev = EventsJ11()
    ev.Data = Data.CopyTree((cutallmc).GetTitle())
    ev.MC = MC.CopyTree((cutalldata).GetTitle())
    
    return ev


def MakePlots(ev, var, cutmc, cutdata, xtitle, nbins, xlow, xup, whatdata, whattrig):
    hData_P     = TH1F("Data_P" , "", nbins, xlow, xup)  # passed
    hData_T     = TH1F("Data_T" , "", nbins, xlow, xup)  # total
    hMC_P       = TH1F("MC_P"   , "", nbins, xlow, xup)  # passed
    hMC_T       = TH1F("MC_T"   , "", nbins, xlow, xup)  # total

    cutmc = "efflumi * PUweight * (triggerFlags[%i]==1)"
    cutdata = "EVENT.json && ( (190456<=EVENT.run && EVENT.run<=193752 && (triggerFlags[%i]==1)) || (193752<=EVENT.run && EVENT.run<=208686 && (triggerFlags[%i]==1)) )"

    if whatdata == "SingleMu" or whatdata == "SingleEl" or (whatdata == "DijetMET" and "METtype1corr.et" in var):
        reftrigA = 23
        reftrig = 23
        if whatdata == "SingleEl":
            reftrigA = 44
            reftrig = 44
        elif whatdata == "DijetMET":
            reftrigA = 49
            reftrig = 50
        if whattrig == "combo":
            ev.MC.Project("MC_P", var, (cutmc + " * (triggerFlags[%i]==1 || triggerFlags[%i]==1 || triggerFlags[%i]==1)") % (reftrig, triggers["MET150"][1], triggers["DijetMET100"][1], triggers["MET80CSV"][1]))
            ev.Data.Project("Data_P", var, (cutdata + " && ( (190456<=EVENT.run && EVENT.run<=193752 && (triggerFlags[%i]==1 || triggerFlags[%i]==1 || triggerFlags[%i]==1)) || (193752<=EVENT.run && EVENT.run<=208686 && (triggerFlags[%i]==1 || triggerFlags[%i]==1 || triggerFlags[%i]==1)) )") % (reftrigA, reftrig, triggers["MET150"][0], triggers["DijetMET100"][0], triggers["MET80CSV"][0], triggers["MET150"][1], triggers["DijetMET100"][1], triggers["MET80CSV"][1]))
        else:
            trigA = triggers[whattrig][0]
            trig = triggers[whattrig][1]
            ev.MC.Project("MC_P", var, (cutmc + " * (triggerFlags[%i]==1)") % (reftrig, trig))
            ev.Data.Project("Data_P", var, (cutdata + " && ( (190456<=EVENT.run && EVENT.run<=193752 && (triggerFlags[%i]==1)) || (193752<=EVENT.run && EVENT.run<=208686 && (triggerFlags[%i]==1)) )") % (reftrigA, reftrig, trigA, trig))
        ev.MC.Project("MC_T", var, cutmc % reftrig)
        ev.Data.Project("Data_T", var, cutdata % (reftrigA, reftrig))
    
    elif whatdata == "DijetMET":
        reftrigA = 49
        reftrig = 50
        if whattrig == "combo":
            ev.MC.Project("MC_P", var, (cutmc + " * triggercorrMET(METtype1corr.et) * (triggerFlags[%i]==1 || triggerFlags[%i]==1 || triggerFlags[%i]==1)") % (reftrig, triggers["MET150"][1], triggers["DijetMET100"][1], triggers["MET80CSV"][1]))
            ev.Data.Project("Data_P", var, (cutdata + " && (!(207883<=EVENT.run && EVENT.run<=208307)) && ( (190456<=EVENT.run && EVENT.run<=193752 && (triggerFlags[%i]==1 || triggerFlags[%i]==1 || triggerFlags[%i]==1)) || (193752<=EVENT.run && EVENT.run<=208686 && (triggerFlags[%i]==1 || triggerFlags[%i]==1 || triggerFlags[%i]==1)) )") % (reftrigA, reftrig, triggers["MET150"][0], triggers["DijetMET100"][0], triggers["MET80CSV"][0], triggers["MET150"][1], triggers["DijetMET100"][1], triggers["MET80CSV"][1]))
        else:
            trigA = triggers[whattrig][0]
            trig = triggers[whattrig][1]
            ev.MC.Project("MC_P", var, (cutmc + " * triggercorrMET(METtype1corr.et) * (triggerFlags[%i]==1)") % (reftrig, trig))
            ev.Data.Project("Data_P", var, (cutdata + " && (!(207883<=EVENT.run && EVENT.run<=208307)) && ( (190456<=EVENT.run && EVENT.run<=193752 && (triggerFlags[%i]==1)) || (193752<=EVENT.run && EVENT.run<=208686 && (triggerFlags[%i]==1)) )") % (reftrigA, reftrig, trigA, trig))
        ev.MC.Project("MC_T", var, (cutmc + " * triggercorrMET(METtype1corr.et)") % reftrig)
        ev.Data.Project("Data_T", var, (cutdata + " && (!(207883<=EVENT.run && EVENT.run<=208307))") % (reftrigA, reftrig))


    if whatdata == "SingleMu" or whatdata == "SingleEl" or (whatdata == "DijetMET" and "METtype1corr.et" in var):
        bini, binj = hData_P.FindFixBin(175), hData_P.FindFixBin(295)
        print hData_P.Integral(bini, binj), hData_T.Integral(bini, binj), hMC_P.Integral(bini, binj), hMC_T.Integral(bini, binj)
        print "high-pt eff (Data, MC): %.3f  %.3f" % (hData_P.Integral(bini, binj) / hData_T.Integral(bini, binj), hMC_P.Integral(bini, binj) / hMC_T.Integral(bini, binj))
        bini, binj = hData_P.FindFixBin(135), hData_P.FindFixBin(165)
        print "med-pt eff  (Data, MC): %.3f  %.3f" % (hData_P.Integral(bini, binj) / hData_T.Integral(bini, binj), hMC_P.Integral(bini, binj) / hMC_T.Integral(bini, binj))
        bini, binj = hData_P.FindFixBin(105), hData_P.FindFixBin(125)
        print "low-pt eff  (Data, MC): %.3f  %.3f" % (hData_P.Integral(bini, binj) / hData_T.Integral(bini, binj), hMC_P.Integral(bini, binj) / hMC_T.Integral(bini, binj))
    
    elif whatdata == "DijetMET":
        bini, binj = hData_P.FindFixBin(0.898), hData_P.FindFixBin(1.001)
        print "CSVT eff (Data, MC): %.3f  %.3f" % (hData_P.Integral(bini, binj) / hData_T.Integral(bini, binj), hMC_P.Integral(bini, binj) / hMC_T.Integral(bini, binj))
        bini, binj = hData_P.FindFixBin(0.679), hData_P.FindFixBin(0.898-0.036)
        print "CSVM eff (Data, MC): %.3f  %.3f" % (hData_P.Integral(bini, binj) / hData_T.Integral(bini, binj), hMC_P.Integral(bini, binj) / hMC_T.Integral(bini, binj))
        bini, binj = hData_P.FindFixBin(0.244), hData_P.FindFixBin(0.679-0.036)
        print "CSVL eff (Data, MC): %.3f  %.3f" % (hData_P.Integral(bini, binj) / hData_T.Integral(bini, binj), hMC_P.Integral(bini, binj) / hMC_T.Integral(bini, binj))
    
    #hData_E = TEfficiency(hData_P, hData_T)
    #hMC_E = TEfficiency(hMC_P, hMC_T)
    #hData_E.Draw("AP")
    #hMC_E.Draw("AP")
    #gData_E = hData_E.GetPaintedGraph()
    #gMC_E = hMC_E.GetPaintedGraph()
    
    gData_E = TGraphAsymmErrors(hData_P, hData_T, "cp")
    gMC_E = TGraphAsymmErrors(hMC_P, hMC_T, "n")
    gMC_E.SetLineColor(colors[whattrig])
    gMC_E.SetMarkerColor(colors[whattrig])
    
    
    if True:
        #print "MakePlots(): Setting up histograms..."

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
        
        h2 = hData_P.Clone("h2")
        h2.Divide(hData_P, hData_T)
        h1 = hMC_P.Clone("h1")
        h1.Divide(hMC_P, hMC_T)
        h1.SetMaximum(1.08)
        h1.SetMinimum(0.5)
        if region == "WjLF" and whattrig == "MET80CSV":
            h1.SetMaximum(0.58)
            h1.SetMinimum(0.0)
        
        h1.GetYaxis().SetNdivisions(505)
        
        ratio = h2.Clone("ratio")
        ratio.Sumw2()
        ratio.SetMarkerSize(0.8)
        #ratio.SetMarkerSize(0.5)
        ratio.Divide(h2, h1, 1., 1., "")
        
        ratio.SetStats(0)
        ratio.SetTitle("")
        ratio.GetXaxis().SetTitle(xtitle)
        ratio.GetYaxis().SetTitle("Data/MC")
        ratio.SetMaximum(1.2)
        ratio.SetMinimum(0.77)
        ratio.SetMarkerStyle(20)
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
        
        ratiounity = TLine(xlow,1,xup,1)
        ratiounity.SetLineStyle(2)
        
        
        # Print data/MC ratios
        if whatdata == "SingleMu" or whatdata == "SingleEl" or (whatdata == "DijetMET" and "METtype1corr.et" in var):
            n = 11
            for i in xrange(n):
                bini = ratio.FindFixBin(105+i*10)
                binc = ratio.GetBinContent(bini)
                if i==0:
                    print "if (met < 100.)  return 0.000;"
                print "if (met < %.0f.)  return %.3f;" % (100+(i+1)*10, binc)
                if i==n-1:
                    print "return 1.;"
        
        elif whatdata == "DijetMET":
            h1.SetMaximum(1.2)
            h1.SetMinimum(0)
            ratio.SetMaximum(1.6)
            ratio.SetMinimum(0.3)
            
            r = ratio.Fit("pol1", "FS" , "", 0.244, 1.001)
            print "if (csv>0.244) return (%.3f + (%.3f) * csv);" % (r.Parameter(0), r.Parameter(1))
            #r = ratio.Fit("pol1", "FS" , "", 0.679, 1.001)
            #print "if (csv>0.679) return (%.3f + (%.3f) * csv);" % (r.Parameter(0), r.Parameter(1))
            #r = ratio.Fit("pol1", "FS+", "", 0.244, 0.679)
            #print "if (csv>0.244) return (%.3f + (%.3f) * csv);" % (r.Parameter(0), r.Parameter(1))
            print "return (%.3f + (%.3f) * 0.244);" % (r.Parameter(0), r.Parameter(1))
        
        
        # Setup legends
        print "MakePlots(): Setting up legends..."
        leg = TLegend(0.67, 0.08, 0.92, 0.24)
        leg.SetFillColor(0)
        leg.SetLineColor(0)
        leg.SetShadowColor(0)
        leg.SetTextFont(62)
        leg.SetTextSize(0.05)
        #leg.AddEntry(gData_E, "Data", "lp")
        leg.AddEntry(gData_E, whatdata, "lp")
        leg.AddEntry(gMC_E, "MC", "lp")
        
        # Draw MC signal and background
        print "MakePlots(): Drawing..."
        pad1.cd()
        if plotLog:  pad1.SetLogy()
        h1.Reset()
        h1.Draw()
        h1.GetXaxis().SetLabelSize(0)
        binwidth = (xup - xlow) / nbins
        ytitle = "HLT efficiency"
        h1.GetYaxis().SetTitle(ytitle)
        
        gMC_E.Draw("P")
        gData_E.Draw("P")
        
        # Draw legends
        leg.Draw()
        latex = TLatex()
        latex.SetNDC()
        latex.SetTextAlign(12)
        latex.SetTextFont(62)
        latex.SetTextSize(0.052)
        latex.DrawLatex(0.19, 0.89, "CMS Preliminary")
        latex.SetTextSize(0.04)
        #latex.DrawLatex(0.19, 0.84, "#sqrt{s} = 8 TeV, L = 19.6 fb^{-1}")
        latex.DrawLatex(0.19, 0.84, "#sqrt{s} = 8 TeV, L = 19.0 fb^{-1}")
        #latex.DrawLatex(0.19, 0.79, "Z(#nu#bar{#nu})H(b#bar{b})")
        latex.DrawLatex(0.19, 0.74, whattrig + " trigger")
        
        # Draw ratio
        pad2.cd()
        pad2.SetGridy(0)
        ratio.Draw("e1")
        ratiounity.Draw()
        
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
        plotname = TString("%s_%s_%s_eff%s_%s" % (channel, region, whatdata, whattrig, var))
        FormatFileName(plotname)
        gPad.Print(plotdir+plotname.Data()+".png")
        gPad.Print(plotdir+plotname.Data()+".pdf")

    return 0


## Main
## ----------------------------------------------------------------------------

if __name__ == "__main__":
    
    vtype = " && Vtype==2 && abs(deltaPhi(vLepton_phi[0], METtype1corr.phi))<2.5"
    if whatdata == "SingleEl":
        vtype = " && Vtype==3"
    cutallmc    = TCut("max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.898 && H.pt>130 && 1<=naJets_Znn && naJets_Znn<=3" + vtype)
    cutalldata  = cutallmc
    if region == "WjLF":
        cutallmc    = TCut("max(hJet_csv_nominal[0],hJet_csv_nominal[1])<=0.898 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.244 && H.pt>130 && naJets_Znn<=0" + vtype)
        cutalldata  = cutallmc
    elif region == "WjHF":
        cutallmc    = TCut("max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.679 && H.pt>130 && naJets_Znn<=3" + vtype)
        cutalldata  = cutallmc
    #cutmc       = TCut("efflumi * PUweight")
    #cutdata     = TCut("EVENT.json && triggerFlags[23]")
    cutmc       = TCut("")  # useless
    cutdata     = TCut("")  # useless
    
#    var = "METtype1corr.et"
#    if whatdata == "SingleMu" or whatdata == "SingleEl" or (whatdata == "DijetMET" and "METtype1corr.et" in var):
#        ev = Read(whatdata, cutallmc, cutalldata)
#        whattrig = "MET150"
#        MakePlots(ev, var, cutmc, cutdata, "E_{T}^{miss} [GeV]", 35, 0, 350, whatdata=whatdata, whattrig=whattrig)
#        whattrig = "DijetMET100"
#        MakePlots(ev, var, cutmc, cutdata, "E_{T}^{miss} [GeV]", 35, 0, 350, whatdata=whatdata, whattrig=whattrig)
#        whattrig = "MET80CSV"
#        MakePlots(ev, var, cutmc, cutdata, "E_{T}^{miss} [GeV]", 35, 0, 350, whatdata=whatdata, whattrig=whattrig)
#        whattrig = "combo"
#        MakePlots(ev, var, cutmc, cutdata, "E_{T}^{miss} [GeV]", 35, 0, 350, whatdata=whatdata, whattrig=whattrig)
    
    var = "max(hJet_csv_nominal[0],hJet_csv_nominal[1])"
    if whatdata == "DijetMET" and not "METtype1corr.et" in var:
        cutallmc    = TCut("max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.244 && H.pt>90 && hJet_pt[0]>50 && hJet_pt[1]>30 && 2<=naJets_Znn && naJets_Znn<=3 && METtype1corr.et>115" + vtype)
        cutalldata  = cutallmc
        if region == "WjLF":
            cutallmc    = TCut("max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.244 && H.pt>90 && hJet_pt[0]>50 && hJet_pt[1]>30 && naJets_Znn<=0 && METtype1corr.et>115" + vtype)
            cutalldata  = cutallmc
        elif region == "WjHF":
            cutallmc    = TCut("max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.244 && H.pt>90 && hJet_pt[0]>50 && hJet_pt[1]>30 && naJets_Znn<=1 && METtype1corr.et>115" + vtype)
            cutalldata  = cutallmc
        ev = Read(whatdata, cutallmc, cutalldata)
        whattrig = "MET80CSV"
        MakePlots(ev, var, cutmc, cutdata, "CSV_{max}", 15, 0.0, 1.08, whatdata=whatdata, whattrig=whattrig)
        
        var = "max(Max$(hJet_csv_nominal),Max$(aJet_csv_nominal))"
        MakePlots(ev, var, cutmc, cutdata, "CSV_{max}(all jets)", 15, 0.0, 1.08, whatdata=whatdata, whattrig=whattrig)