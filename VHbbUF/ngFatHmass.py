import array
from ROOT import gROOT, AddressOf

gROOT.ProcessLine("typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVectorM;")
gROOT.ProcessLine("typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double> > LorentzVectorE;")

gROOT.LoadMacro("HelperFunctions.h")
gROOT.LoadMacro("HelperVHbbDataFormats.h")
gROOT.LoadMacro("$HOME/style-CMSTDR.C");
from ROOT import TFile, TTree, TTreeFormula, TH1F, TLatex, TLegend, TColor, gPad, gROOT, HiggsInfo, METInfo, genParticleInfo, LorentzVectorM, LorentzVectorE, setTDRStyle, deltaR, deltaPhi
setTDRStyle()
gROOT.SetBatch(1)

f = TFile.Open("skim_ZnnH_baseline/Step3_ZnnH125.root")
tree = f.Get("tree")

H = HiggsInfo()
tree.SetBranchAddress("H", AddressOf(H, "HiggsFlag"))
METtype1corr = METInfo()
tree.SetBranchAddress("METtype1corr", AddressOf(METtype1corr, "et"))
genZ = genParticleInfo()
tree.SetBranchAddress("genZ", AddressOf(genZ, "mass"))
genH = genParticleInfo()
tree.SetBranchAddress("genH", AddressOf(genH, "mass"))
genB = genParticleInfo()
tree.SetBranchAddress("genB", AddressOf(genB, "mass"))
genBbar = genParticleInfo()
tree.SetBranchAddress("genBbar", AddressOf(genBbar, "mass"))


#cut = "Vtype==4 && hJet_pt[0]>60 && hJet_pt[1]>30 && METtype1corr.et>170 && H.pt>130 && hJet_csv_nominal[0]>0.244 && hJet_csv_nominal[1]>0.244 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && H.mass>0"
cut = "Vtype==4 && hJet_pt[0]>60 && hJet_pt[1]>30 && hJet_csv_nominal[0]>0.244 && hJet_csv_nominal[1]>0.244 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && H.mass>0"
#fjcut = "Vtype==4 && fathFilterJets_pt[0]>60 && fathFilterJets_pt[1]>30 && METtype1corr.et>170 && FatH.pt>130 && fathFilterJets_csv[0]>0.244 && fathFilterJets_csv[1]>0.244 && abs(deltaPhi(FatH.phi,METtype1corr.phi))>2.0 && nfathFilterJets>0 && fathFilterJets_pt[0]>0 && fathFilterJets_pt[1]>0 && FatH.mass>0"
fjcut = "Vtype==4 && fathFilterJets_pt[0]>60 && fathFilterJets_pt[1]>30 && fathFilterJets_csv[0]>0.244 && fathFilterJets_csv[1]>0.244 && abs(deltaPhi(FatH.phi,METtype1corr.phi))>2.0 && nfathFilterJets>0 && fathFilterJets_pt[0]>0 && fathFilterJets_pt[1]>0 && FatH.mass>0"
fcut = TTreeFormula("ttfcut", cut, tree)
ffjcut = TTreeFormula("ttffjcut", fjcut, tree)

nbins, xmin, xmax = 70, 60., 200.
fathmassnorm2_5  = TH1F("FatHmassNorm2_5" , "", nbins, xmin, xmax)
fathmassnorm3_5  = TH1F("FatHmassNorm3_5" , "", nbins, xmin, xmax)
fathmassnorm2_10 = TH1F("FatHmassNorm2_10", "", nbins, xmin, xmax)
fathmassnorm3_10 = TH1F("FatHmassNorm3_10", "", nbins, xmin, xmax)
fathmassnorm2_15 = TH1F("FatHmassNorm2_15", "", nbins, xmin, xmax)
fathmassnorm3_15 = TH1F("FatHmassNorm3_15", "", nbins, xmin, xmax)
fathmassnorm2_20 = TH1F("FatHmassNorm2_20", "", nbins, xmin, xmax)
fathmassnorm3_20 = TH1F("FatHmassNorm3_20", "", nbins, xmin, xmax)
fathmassreg2_5  = TH1F("FatHmassReg2_5" , "", nbins, xmin, xmax)
fathmassreg3_5  = TH1F("FatHmassReg3_5" , "", nbins, xmin, xmax)
fathmassreg2_10 = TH1F("FatHmassReg2_10", "", nbins, xmin, xmax)
fathmassreg3_10 = TH1F("FatHmassReg3_10", "", nbins, xmin, xmax)
fathmassreg2_15 = TH1F("FatHmassReg2_15", "", nbins, xmin, xmax)
fathmassreg3_15 = TH1F("FatHmassReg3_15", "", nbins, xmin, xmax)
fathmassreg2_20 = TH1F("FatHmassReg2_20", "", nbins, xmin, xmax)
fathmassreg3_20 = TH1F("FatHmassReg3_20", "", nbins, xmin, xmax)

hmassnorm_pt1 = TH1F("HmassNorm_pt1" , "", nbins, xmin, xmax)
hmassnorm_pt2 = TH1F("HmassNorm_pt2" , "", nbins, xmin, xmax)
hmassnorm_pt3 = TH1F("HmassNorm_pt3" , "", nbins, xmin, xmax)
hmassnorm_pt4 = TH1F("HmassNorm_pt4" , "", nbins, xmin, xmax)
hmassnorm_pt5 = TH1F("HmassNorm_pt5" , "", nbins, xmin, xmax)
hmassnorm_pt6 = TH1F("HmassNorm_pt6" , "", nbins, xmin, xmax)
hmassreg_pt1 = TH1F("HmassReg_pt1" , "", nbins, xmin, xmax)
hmassreg_pt2 = TH1F("HmassReg_pt2" , "", nbins, xmin, xmax)
hmassreg_pt3 = TH1F("HmassReg_pt3" , "", nbins, xmin, xmax)
hmassreg_pt4 = TH1F("HmassReg_pt4" , "", nbins, xmin, xmax)
hmassreg_pt5 = TH1F("HmassReg_pt5" , "", nbins, xmin, xmax)
hmassreg_pt6 = TH1F("HmassReg_pt6" , "", nbins, xmin, xmax)
fathmassnorm_pt1 = TH1F("FatHmassNorm_pt1" , "", nbins, xmin, xmax)
fathmassnorm_pt2 = TH1F("FatHmassNorm_pt2" , "", nbins, xmin, xmax)
fathmassnorm_pt3 = TH1F("FatHmassNorm_pt3" , "", nbins, xmin, xmax)
fathmassnorm_pt4 = TH1F("FatHmassNorm_pt4" , "", nbins, xmin, xmax)
fathmassnorm_pt5 = TH1F("FatHmassNorm_pt5" , "", nbins, xmin, xmax)
fathmassnorm_pt6 = TH1F("FatHmassNorm_pt6" , "", nbins, xmin, xmax)
fathmassreg_pt1 = TH1F("FatHmassReg_pt1" , "", nbins, xmin, xmax)
fathmassreg_pt2 = TH1F("FatHmassReg_pt2" , "", nbins, xmin, xmax)
fathmassreg_pt3 = TH1F("FatHmassReg_pt3" , "", nbins, xmin, xmax)
fathmassreg_pt4 = TH1F("FatHmassReg_pt4" , "", nbins, xmin, xmax)
fathmassreg_pt5 = TH1F("FatHmassReg_pt5" , "", nbins, xmin, xmax)
fathmassreg_pt6 = TH1F("FatHmassReg_pt6" , "", nbins, xmin, xmax)

nevt = tree.GetEntries()
for ievt in xrange(0,nevt):
    tree.LoadTree(ievt)  # used by TTreeFormula
    tree.GetEntry(ievt)

    fcut.GetNdata()
    ffjcut.GetNdata()
    pass_cut = bool(fcut.EvalInstance())
    pass_fjcut = bool(ffjcut.EvalInstance())
    
    if(tree.hJet_pt[0]>0. and tree.hJet_pt[1]>0. and pass_cut):
        
        #if METtype1corr.et < 130:
        #    pass
        #elif METtype1corr.et < 150:
        #    hmassnorm_pt1.Fill(tree.HmassNorm)
        #    hmassreg_pt1.Fill(tree.HmassReg)
        #elif METtype1corr.et < 170:
        #    hmassnorm_pt2.Fill(tree.HmassNorm)
        #    hmassreg_pt2.Fill(tree.HmassReg)
        #elif METtype1corr.et < 200:
        #    hmassnorm_pt3.Fill(tree.HmassNorm)
        #    hmassreg_pt3.Fill(tree.HmassReg)
        #elif METtype1corr.et < 250:
        #    hmassnorm_pt4.Fill(tree.HmassNorm)
        #    hmassreg_pt4.Fill(tree.HmassReg)
        #else:
        #    hmassnorm_pt5.Fill(tree.HmassNorm)
        #    hmassreg_pt5.Fill(tree.HmassReg)
        
        if tree.HptNorm < 130:
            pass
        elif tree.HptNorm < 150:
            hmassnorm_pt1.Fill(tree.HmassNorm)
        elif tree.HptNorm < 170:
            hmassnorm_pt2.Fill(tree.HmassNorm)
        elif tree.HptNorm < 200:
            hmassnorm_pt3.Fill(tree.HmassNorm)
        elif tree.HptNorm < 250:
            hmassnorm_pt4.Fill(tree.HmassNorm)
        elif tree.HptNorm < 400:
            hmassnorm_pt5.Fill(tree.HmassNorm)
        else:
            hmassnorm_pt6.Fill(tree.HmassNorm)
        
        if tree.HptReg < 130:
            pass
        elif tree.HptReg < 150:
            hmassreg_pt1.Fill(tree.HmassReg)
        elif tree.HptReg < 170:
            hmassreg_pt2.Fill(tree.HmassReg)
        elif tree.HptReg < 200:
            hmassreg_pt3.Fill(tree.HmassReg)
        elif tree.HptReg < 250:
            hmassreg_pt4.Fill(tree.HmassReg)
        elif tree.HptReg < 400:
            hmassreg_pt5.Fill(tree.HmassReg)
        else:
            hmassreg_pt6.Fill(tree.HmassReg)

    if(tree.nfathFilterJets>0 and tree.fathFilterJets_pt[0]>0. and tree.fathFilterJets_pt[1]>0. and pass_fjcut):
        
        #if METtype1corr.et < 130:
        #    pass
        #elif METtype1corr.et < 150:
        #    fathmassnorm_pt1.Fill(tree.FatHmassNorm)
        #    fathmassreg_pt1.Fill(tree.FatHmassReg)
        #elif METtype1corr.et < 170:
        #    fathmassnorm_pt2.Fill(tree.FatHmassNorm)
        #    fathmassreg_pt2.Fill(tree.FatHmassReg)
        #elif METtype1corr.et < 200:
        #    fathmassnorm_pt3.Fill(tree.FatHmassNorm)
        #    fathmassreg_pt3.Fill(tree.FatHmassReg)
        #elif METtype1corr.et < 250:
        #    fathmassnorm_pt4.Fill(tree.FatHmassNorm)
        #    fathmassreg_pt4.Fill(tree.FatHmassReg)
        #else:
        #    fathmassnorm_pt5.Fill(tree.FatHmassNorm)
        #    fathmassreg_pt5.Fill(tree.FatHmassReg)
        
        if tree.FatHptNorm < 130:
            pass
        elif tree.FatHptNorm < 150:
            fathmassnorm_pt1.Fill(tree.FatHmassNorm)
        elif tree.FatHptNorm < 170:
            fathmassnorm_pt2.Fill(tree.FatHmassNorm)
        elif tree.FatHptNorm < 200:
            fathmassnorm_pt3.Fill(tree.FatHmassNorm)
        elif tree.FatHptNorm < 250:
            fathmassnorm_pt4.Fill(tree.FatHmassNorm)
        elif tree.FatHptNorm < 400:
            fathmassnorm_pt5.Fill(tree.FatHmassNorm)
        else:
            fathmassnorm_pt6.Fill(tree.FatHmassNorm)
        
        if tree.FatHptReg < 130:
            pass
        elif tree.FatHptReg < 150:
            fathmassreg_pt1.Fill(tree.FatHmassReg)
        elif tree.FatHptReg < 170:
            fathmassreg_pt2.Fill(tree.FatHmassReg)
        elif tree.FatHptReg < 200:
            fathmassreg_pt3.Fill(tree.FatHmassReg)
        elif tree.FatHptReg < 250:
            fathmassreg_pt4.Fill(tree.FatHmassReg)
        elif tree.FatHptReg < 400:
            fathmassreg_pt5.Fill(tree.FatHmassReg)
        else:
            fathmassreg_pt6.Fill(tree.FatHmassReg)
        
        if METtype1corr.et > 170:
            fj0 = LorentzVectorE(tree.fathFilterJets_pt[0], tree.fathFilterJets_eta[0], tree.fathFilterJets_phi[0], tree.fathFilterJets_e[0])
            fj1 = LorentzVectorE(tree.fathFilterJets_pt[1], tree.fathFilterJets_eta[1], tree.fathFilterJets_phi[1], tree.fathFilterJets_e[1])
            fj2 = LorentzVectorE()
            if (tree.nfathFilterJets>2 and tree.fathFilterJets_pt[2]>0.):
                fj2 = LorentzVectorE(tree.fathFilterJets_pt[2], tree.fathFilterJets_eta[2], tree.fathFilterJets_phi[2], tree.fathFilterJets_e[2])
            fh2 = fj0 + fj1
            fh3 = fj0 + fj1 + fj2
            fh3_5  = fj0 + fj1 + fj2 if (tree.nfathFilterJets>2 and tree.fathFilterJets_pt[2]> 5.) else fj0 + fj1
            fh3_10 = fj0 + fj1 + fj2 if (tree.nfathFilterJets>2 and tree.fathFilterJets_pt[2]>10.) else fj0 + fj1
            fh3_15 = fj0 + fj1 + fj2 if (tree.nfathFilterJets>2 and tree.fathFilterJets_pt[2]>15.) else fj0 + fj1
            fh3_20 = fj0 + fj1 + fj2 if (tree.nfathFilterJets>2 and tree.fathFilterJets_pt[2]>20.) else fj0 + fj1
            
            fj0r = LorentzVectorE(tree.fathFilterJets_ptReg[0], tree.fathFilterJets_eta[0], tree.fathFilterJets_phi[0], tree.fathFilterJets_e[0] * tree.fathFilterJets_ptReg[0] / tree.fathFilterJets_pt[0])
            fj1r = LorentzVectorE(tree.fathFilterJets_ptReg[1], tree.fathFilterJets_eta[1], tree.fathFilterJets_phi[1], tree.fathFilterJets_e[1] * tree.fathFilterJets_ptReg[1] / tree.fathFilterJets_pt[1])
            fj2r = LorentzVectorE()
            if (tree.nfathFilterJets>2 and tree.fathFilterJets_pt[2]>0.):
                fj2r = LorentzVectorE(tree.fathFilterJets_ptReg[2], tree.fathFilterJets_eta[2], tree.fathFilterJets_phi[2], tree.fathFilterJets_e[2] * tree.fathFilterJets_ptReg[2] / tree.fathFilterJets_pt[2])
            fh2r = fj0r + fj1r
            fh3r = fj0r + fj1r + fj2r
            fh3r_5  = fj0r + fj1r + fj2r if (tree.nfathFilterJets>2 and tree.fathFilterJets_ptReg[2]> 5.) else fj0r + fj1r
            fh3r_10 = fj0r + fj1r + fj2r if (tree.nfathFilterJets>2 and tree.fathFilterJets_ptReg[2]>10.) else fj0r + fj1r
            fh3r_15 = fj0r + fj1r + fj2r if (tree.nfathFilterJets>2 and tree.fathFilterJets_ptReg[2]>15.) else fj0r + fj1r
            fh3r_20 = fj0r + fj1r + fj2r if (tree.nfathFilterJets>2 and tree.fathFilterJets_ptReg[2]>20.) else fj0r + fj1r
            
            if (tree.nfathFilterJets>2 and tree.fathFilterJets_pt[2]> 5.):
                fathmassnorm3_5.Fill(fh3_5.mass())
            else:
                fathmassnorm2_5.Fill(fh3_5.mass())
            if (tree.nfathFilterJets>2 and tree.fathFilterJets_pt[2]>10.):
                fathmassnorm3_10.Fill(fh3_10.mass())
            else:
                fathmassnorm2_10.Fill(fh3_10.mass())
            if (tree.nfathFilterJets>2 and tree.fathFilterJets_pt[2]>15.):
                fathmassnorm3_15.Fill(fh3_15.mass())
            else:
                fathmassnorm2_15.Fill(fh3_15.mass())
            if (tree.nfathFilterJets>2 and tree.fathFilterJets_pt[2]>20.):
                fathmassnorm3_20.Fill(fh3_20.mass())
            else:
                fathmassnorm2_20.Fill(fh3_20.mass())
            
            
            if (tree.nfathFilterJets>2 and tree.fathFilterJets_ptReg[2]> 5.):
                fathmassreg3_5.Fill(fh3r_5.mass())
            else:
                fathmassreg2_5.Fill(fh3r_5.mass())
            if (tree.nfathFilterJets>2 and tree.fathFilterJets_ptReg[2]>10.):
                fathmassreg3_10.Fill(fh3r_10.mass())
            else:
                fathmassreg2_10.Fill(fh3r_10.mass())
            if (tree.nfathFilterJets>2 and tree.fathFilterJets_ptReg[2]>15.):
                fathmassreg3_15.Fill(fh3r_15.mass())
            else:
                fathmassreg2_15.Fill(fh3r_15.mass())
            if (tree.nfathFilterJets>2 and tree.fathFilterJets_ptReg[2]>20.):
                fathmassreg3_20.Fill(fh3r_20.mass())
            else:
                fathmassreg2_20.Fill(fh3r_20.mass())

            b0 = LorentzVectorM(genB.pt, genB.eta, genB.phi, genB.mass)
            b1 = LorentzVectorM(genBbar.pt, genBbar.eta, genBbar.phi, genBbar.mass)
            bb = b0 + b1

            verbose = False;
            if (verbose):
                print "ievt %6i" % ievt, "Nfj =", tree.nfathFilterJets, "fh2 = %6.2f, %6.2f" % (fh2.Pt(), fh2.M()), "fh3 = %6.2f, %6.2f" % (fh3.Pt(), fh3.M()) #, "fh3_5 = %6.2f, %6.2f" % (fh3_5.Pt(), fh3_5.M())
                print "  FJ        Norm: %6.2f Reg: %6.2f" % (tree.FatHmassNorm, tree.FatHmassReg), "fj1 = %6.2f" % (tree.fathFilterJets_pt[0]), "fj2 = %6.2f" % tree.fathFilterJets_pt[1], "fj3 = %6.2f" % tree.fathFilterJets_pt[2]
                print "  DJ        Norm: %6.2f Reg: %6.2f" % (tree.HmassNorm, tree.HmassReg), "dj1 = %6.2f" % (tree.hJet_pt[0]), "dj2 = %6.2f" % tree.hJet_pt[1]
                print "  MC        DR  : %6.2f pT : %6.2f" % (deltaR(genB.eta,genB.phi,genBbar.eta,genBbar.phi), bb.pt()), " b1  = %6.2f" % genB.pt, "b2  = %6.2f" % genBbar.pt


#-------------------------------------------------------------------------------
if True:
    text = TLatex(); text.SetNDC(); text.SetTextSize(0.025)
    kBlack = 1
    kRed = TColor.GetColor("#FF0000")
    kDRed = TColor.GetColor("#993333")
    kGreen = TColor.GetColor("#00FF00")
    kDGreen = TColor.GetColor("#339933")
    fxmin2r, fxmax2r = 105., 135.
    fxmin3r, fxmax3r = 105., 135.
    fxmin2, fxmax2 = 95., 130.
    fxmin3, fxmax3 = 95., 130.
    down = 0.12

    h2 = fathmassnorm2_5
    h3 = fathmassnorm3_5
    h2r = fathmassreg2_5
    h3r = fathmassreg3_5
    ymax = max(h2r.GetMaximum(), h3r.GetMaximum()) * 1.2
    h2r.Fit("gaus", "I", "", fxmin2r, fxmax2r)
    h3r.Fit("gaus", "I", "", fxmin3r, fxmax3r)
    h2.SetStats(0); h2.SetLineWidth(1); h2.SetMarkerSize(0); h2.SetLineColor(kDGreen)
    h3.SetStats(0); h3.SetLineWidth(1); h3.SetMarkerSize(0); h3.SetLineColor(kDRed)
    h2r.SetStats(0); h2r.SetLineWidth(2); h2r.SetMarkerSize(0); h2r.SetLineColor(kGreen)
    h3r.SetStats(0); h3r.SetLineWidth(2); h3r.SetMarkerSize(0); h3r.SetLineColor(kRed)
    h3r.SetMaximum(ymax)
    h3r.SetXTitle("H(125)^{fj} mass [GeV]"); h3r.Draw()
    h2.Draw("same")
    h3.Draw("same")
    h2r.Draw("same")
    h3r.Draw("same")

    leg = TLegend(0.58, 0.78, 0.92, 0.92)
    leg.AddEntry(h2r, "# filter jets = 2")
    leg.AddEntry(h3r, "# filter jets = 3")
    leg.Draw()
    text.SetTextColor(kBlack)
    text.DrawLatex(0.1, 0.97, "3rd filter jet p_{T} > 5 GeV")

    gaus = h2r.GetFunction("gaus")
    text.SetTextColor(kGreen)
    text.DrawLatex(0.2, 0.90, "gaus fit [%.0f,%.0f]" %(fxmin2r, fxmax2r))
    text.DrawLatex(0.2, 0.87, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    text.DrawLatex(0.2, 0.84, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gaus = h3r.GetFunction("gaus")
    text.SetTextColor(kRed)
    text.DrawLatex(0.2, 0.90-down, "gaus fit [%.0f,%.0f]" %(fxmin3r, fxmax3r))
    text.DrawLatex(0.2, 0.87-down, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    text.DrawLatex(0.2, 0.84-down, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gPad.Print("plots/breg_FatHmassReg2vs3_5.png")

    h2.Fit("gaus", "I", "", fxmin2, fxmax2)
    h3.Fit("gaus", "I", "", fxmin3, fxmax3)
    h3.SetMaximum(ymax)
    h3.SetXTitle("H(125)^{fj} mass [GeV]"); h3.Draw()
    h2.Draw("same")
    h3.Draw("same")

    leg = TLegend(0.58, 0.78, 0.92, 0.92)
    leg.AddEntry(h2, "# filter jets = 2")
    leg.AddEntry(h3, "# filter jets = 3")
    leg.Draw()
    text.SetTextColor(kBlack)
    text.DrawLatex(0.1, 0.97, "3rd filter jet p_{T} > 5 GeV")

    gaus = h2.GetFunction("gaus")
    text.SetTextColor(kDGreen)
    text.DrawLatex(0.2, 0.90, "gaus fit [%.0f,%.0f]" %(fxmin2, fxmax2))
    text.DrawLatex(0.2, 0.87, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    text.DrawLatex(0.2, 0.84, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gaus = h3.GetFunction("gaus")
    text.SetTextColor(kDRed)
    text.DrawLatex(0.2, 0.90-down, "gaus fit [%.0f,%.0f]" %(fxmin3, fxmax3))
    text.DrawLatex(0.2, 0.87-down, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    text.DrawLatex(0.2, 0.84-down, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gPad.Print("plots/breg_FatHmassNorm2vs3_5.png")

    #-------------------------------------------------------------------------------
    h2 = fathmassnorm2_10
    h3 = fathmassnorm3_10
    h2r = fathmassreg2_10
    h3r = fathmassreg3_10
    ymax = max(h2r.GetMaximum(), h3r.GetMaximum()) * 1.2
    h2r.Fit("gaus", "I", "", fxmin2r, fxmax2r)
    h3r.Fit("gaus", "I", "", fxmin3r, fxmax3r)
    h2.SetStats(0); h2.SetLineWidth(1); h2.SetMarkerSize(0); h2.SetLineColor(kDGreen)
    h3.SetStats(0); h3.SetLineWidth(1); h3.SetMarkerSize(0); h3.SetLineColor(kDRed)
    h2r.SetStats(0); h2r.SetLineWidth(2); h2r.SetMarkerSize(0); h2r.SetLineColor(kGreen)
    h3r.SetStats(0); h3r.SetLineWidth(2); h3r.SetMarkerSize(0); h3r.SetLineColor(kRed)
    h3r.SetMaximum(ymax)
    h3r.SetXTitle("H(125)^{fj} mass [GeV]"); h3r.Draw()
    h2.Draw("same")
    h3.Draw("same")
    h2r.Draw("same")
    h3r.Draw("same")

    leg = TLegend(0.58, 0.78, 0.92, 0.92)
    leg.AddEntry(h2r, "# filter jets = 2")
    leg.AddEntry(h3r, "# filter jets = 3")
    leg.Draw()
    text.SetTextColor(kBlack)
    text.DrawLatex(0.1, 0.97, "3rd filter jet p_{T} > 10 GeV")

    gaus = h2r.GetFunction("gaus")
    text.SetTextColor(kGreen)
    text.DrawLatex(0.2, 0.90, "gaus fit [%.0f,%.0f]" %(fxmin2r, fxmax2r))
    text.DrawLatex(0.2, 0.87, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    text.DrawLatex(0.2, 0.84, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gaus = h3r.GetFunction("gaus")
    text.SetTextColor(kRed)
    text.DrawLatex(0.2, 0.90-down, "gaus fit [%.0f,%.0f]" %(fxmin3r, fxmax3r))
    text.DrawLatex(0.2, 0.87-down, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    text.DrawLatex(0.2, 0.84-down, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gPad.Print("plots/breg_FatHmassReg2vs3_10.png")

    h2.Fit("gaus", "I", "", fxmin2, fxmax2)
    h3.Fit("gaus", "I", "", fxmin3, fxmax3)
    h3.SetMaximum(ymax)
    h3.SetXTitle("H(125)^{fj} mass [GeV]"); h3.Draw()
    h2.Draw("same")
    h3.Draw("same")

    leg = TLegend(0.58, 0.78, 0.92, 0.92)
    leg.AddEntry(h2, "# filter jets = 2")
    leg.AddEntry(h3, "# filter jets = 3")
    leg.Draw()
    text.SetTextColor(kBlack)
    text.DrawLatex(0.1, 0.97, "3rd filter jet p_{T} > 10 GeV")

    gaus = h2.GetFunction("gaus")
    text.SetTextColor(kDGreen)
    text.DrawLatex(0.2, 0.90, "gaus fit [%.0f,%.0f]" %(fxmin2, fxmax2))
    text.DrawLatex(0.2, 0.87, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    text.DrawLatex(0.2, 0.84, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gaus = h3.GetFunction("gaus")
    text.SetTextColor(kDRed)
    text.DrawLatex(0.2, 0.90-down, "gaus fit [%.0f,%.0f]" %(fxmin3, fxmax3))
    text.DrawLatex(0.2, 0.87-down, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    text.DrawLatex(0.2, 0.84-down, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gPad.Print("plots/breg_FatHmassNorm2vs3_10.png")

    #-------------------------------------------------------------------------------
    h2 = fathmassnorm2_15
    h3 = fathmassnorm3_15
    h2r = fathmassreg2_15
    h3r = fathmassreg3_15
    ymax = max(h2r.GetMaximum(), h3r.GetMaximum()) * 1.2
    h2r.Fit("gaus", "I", "", fxmin2r, fxmax2r)
    h3r.Fit("gaus", "I", "", fxmin3r, fxmax3r)
    h2.SetStats(0); h2.SetLineWidth(1); h2.SetMarkerSize(0); h2.SetLineColor(kDGreen)
    h3.SetStats(0); h3.SetLineWidth(1); h3.SetMarkerSize(0); h3.SetLineColor(kDRed)
    h2r.SetStats(0); h2r.SetLineWidth(2); h2r.SetMarkerSize(0); h2r.SetLineColor(kGreen)
    h3r.SetStats(0); h3r.SetLineWidth(2); h3r.SetMarkerSize(0); h3r.SetLineColor(kRed)
    h3r.SetMaximum(ymax)
    h3r.SetXTitle("H(125)^{fj} mass [GeV]"); h3r.Draw()
    h2.Draw("same")
    h3.Draw("same")
    h2r.Draw("same")
    h3r.Draw("same")

    leg = TLegend(0.58, 0.78, 0.92, 0.92)
    leg.AddEntry(h2r, "# filter jets = 2")
    leg.AddEntry(h3r, "# filter jets = 3")
    leg.Draw()
    text.SetTextColor(kBlack)
    text.DrawLatex(0.1, 0.97, "3rd filter jet p_{T} > 15 GeV")

    gaus = h2r.GetFunction("gaus")
    text.SetTextColor(kGreen)
    text.DrawLatex(0.2, 0.90, "gaus fit [%.0f,%.0f]" %(fxmin2r, fxmax2r))
    text.DrawLatex(0.2, 0.87, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    text.DrawLatex(0.2, 0.84, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gaus = h3r.GetFunction("gaus")
    text.SetTextColor(kRed)
    text.DrawLatex(0.2, 0.90-down, "gaus fit [%.0f,%.0f]" %(fxmin3r, fxmax3r))
    text.DrawLatex(0.2, 0.87-down, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    text.DrawLatex(0.2, 0.84-down, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gPad.Print("plots/breg_FatHmassReg2vs3_15.png")

    h2.Fit("gaus", "I", "", fxmin2, fxmax2)
    h3.Fit("gaus", "I", "", fxmin3, fxmax3)
    h3.SetMaximum(ymax)
    h3.SetXTitle("H(125)^{fj} mass [GeV]"); h3.Draw()
    h2.Draw("same")
    h3.Draw("same")

    leg = TLegend(0.58, 0.78, 0.92, 0.92)
    leg.AddEntry(h2, "# filter jets = 2")
    leg.AddEntry(h3, "# filter jets = 3")
    leg.Draw()
    text.SetTextColor(kBlack)
    text.DrawLatex(0.1, 0.97, "3rd filter jet p_{T} > 15 GeV")

    gaus = h2.GetFunction("gaus")
    text.SetTextColor(kDGreen)
    text.DrawLatex(0.2, 0.90, "gaus fit [%.0f,%.0f]" %(fxmin2, fxmax2))
    text.DrawLatex(0.2, 0.87, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    text.DrawLatex(0.2, 0.84, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gaus = h3.GetFunction("gaus")
    text.SetTextColor(kDRed)
    text.DrawLatex(0.2, 0.90-down, "gaus fit [%.0f,%.0f]" %(fxmin3, fxmax3))
    text.DrawLatex(0.2, 0.87-down, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    text.DrawLatex(0.2, 0.84-down, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gPad.Print("plots/breg_FatHmassNorm2vs3_15.png")

    #-------------------------------------------------------------------------------
    h2 = fathmassnorm2_20
    h3 = fathmassnorm3_20
    h2r = fathmassreg2_20
    h3r = fathmassreg3_20
    ymax = max(h2r.GetMaximum(), h3r.GetMaximum()) * 1.2
    h2r.Fit("gaus", "I", "", fxmin2r, fxmax2r)
    h3r.Fit("gaus", "I", "", fxmin3r, fxmax3r)
    h2.SetStats(0); h2.SetLineWidth(1); h2.SetMarkerSize(0); h2.SetLineColor(kDGreen)
    h3.SetStats(0); h3.SetLineWidth(1); h3.SetMarkerSize(0); h3.SetLineColor(kDRed)
    h2r.SetStats(0); h2r.SetLineWidth(2); h2r.SetMarkerSize(0); h2r.SetLineColor(kGreen)
    h3r.SetStats(0); h3r.SetLineWidth(2); h3r.SetMarkerSize(0); h3r.SetLineColor(kRed)
    h3r.SetMaximum(ymax)
    h3r.SetXTitle("H(125)^{fj} mass [GeV]"); h3r.Draw()
    h2.Draw("same")
    h3.Draw("same")
    h2r.Draw("same")
    h3r.Draw("same")

    leg = TLegend(0.58, 0.78, 0.92, 0.92)
    leg.AddEntry(h2r, "# filter jets = 2")
    leg.AddEntry(h3r, "# filter jets = 3")
    leg.Draw()
    text.SetTextColor(kBlack)
    text.DrawLatex(0.1, 0.97, "3rd filter jet p_{T} > 20 GeV")

    gaus = h2r.GetFunction("gaus")
    text.SetTextColor(kGreen)
    text.DrawLatex(0.2, 0.90, "gaus fit [%.0f,%.0f]" %(fxmin2r, fxmax2r))
    text.DrawLatex(0.2, 0.87, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    text.DrawLatex(0.2, 0.84, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gaus = h3r.GetFunction("gaus")
    text.SetTextColor(kRed)
    text.DrawLatex(0.2, 0.90-down, "gaus fit [%.0f,%.0f]" %(fxmin3r, fxmax3r))
    text.DrawLatex(0.2, 0.87-down, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    text.DrawLatex(0.2, 0.84-down, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gPad.Print("plots/breg_FatHmassReg2vs3_20.png")

    h2.Fit("gaus", "I", "", fxmin2, fxmax2)
    h3.Fit("gaus", "I", "", fxmin3, fxmax3)
    h3.SetMaximum(ymax)
    h3.SetXTitle("H(125)^{fj} mass [GeV]"); h3.Draw()
    h2.Draw("same")
    h3.Draw("same")

    leg = TLegend(0.58, 0.78, 0.92, 0.92)
    leg.AddEntry(h2, "# filter jets = 2")
    leg.AddEntry(h3, "# filter jets = 3")
    leg.Draw()
    text.SetTextColor(kBlack)
    text.DrawLatex(0.1, 0.97, "3rd filter jet p_{T} > 20 GeV")

    gaus = h2.GetFunction("gaus")
    text.SetTextColor(kDGreen)
    text.DrawLatex(0.2, 0.90, "gaus fit [%.0f,%.0f]" %(fxmin2, fxmax2))
    text.DrawLatex(0.2, 0.87, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    text.DrawLatex(0.2, 0.84, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gaus = h3.GetFunction("gaus")
    text.SetTextColor(kDRed)
    text.DrawLatex(0.2, 0.90-down, "gaus fit [%.0f,%.0f]" %(fxmin3, fxmax3))
    text.DrawLatex(0.2, 0.87-down, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    text.DrawLatex(0.2, 0.84-down, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gPad.Print("plots/breg_FatHmassNorm2vs3_20.png")


#-------------------------------------------------------------------------------
if True:
    text = TLatex(); text.SetNDC(); text.SetTextSize(0.025)
    kBlack = 1
    kRed = TColor.GetColor("#FF0000")
    kDRed = TColor.GetColor("#993333")
    kCyan = TColor.GetColor("#00FFFF")
    kDCyan = TColor.GetColor("#339999")
    widths = []
    down = 0.12

    h2 = hmassnorm_pt1
    h3 = fathmassnorm_pt1
    h2r = hmassreg_pt1
    h3r = fathmassreg_pt1
    fxmin2r, fxmax2r = 100., 140.
    fxmin3r, fxmax3r = 60., 100.
    fxmin2, fxmax2 = 100, 140.
    fxmin3, fxmax3 = 60., 100.
    ymax = max(h2r.GetMaximum(), h3r.GetMaximum()) * 1.2
    h2r.Fit("gaus", "I", "", fxmin2r, fxmax2r)
    h3r.Fit("gaus", "I", "", fxmin3r, fxmax3r)
    h2.SetStats(0); h2.SetLineWidth(1); h2.SetMarkerSize(0); h2.SetLineColor(kDCyan)
    h3.SetStats(0); h3.SetLineWidth(1); h3.SetMarkerSize(0); h3.SetLineColor(kDRed)
    h2r.SetStats(0); h2r.SetLineWidth(2); h2r.SetMarkerSize(0); h2r.SetLineColor(kCyan)
    h3r.SetStats(0); h3r.SetLineWidth(2); h3r.SetMarkerSize(0); h3r.SetLineColor(kRed)
    h3r.SetMaximum(ymax)
    h3r.SetXTitle("H(125) mass [GeV]"); h3r.Draw()
    h2.Draw("same")
    h3.Draw("same")
    h2r.Draw("same")
    h3r.Draw("same")

    leg = TLegend(0.58, 0.78, 0.92, 0.92)
    leg.AddEntry(h2r, "Dijet-system")
    leg.AddEntry(h3r, "Filterjet-system")
    leg.Draw()
    text.SetTextColor(kBlack)
    text.DrawLatex(0.1, 0.97, "130 < H p_{T} #leq 150 GeV")

    widths.append([])
    gaus = h2r.GetFunction("gaus")
    text.SetTextColor(kCyan)
    text.DrawLatex(0.2, 0.90, "gaus fit [%.0f,%.0f]" %(fxmin2r, fxmax2r))
    text.DrawLatex(0.2, 0.87, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    widths[-1].append(gaus.GetParameter(2)/gaus.GetParameter(1))
    text.DrawLatex(0.2, 0.84, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gaus = h3r.GetFunction("gaus")
    text.SetTextColor(kRed)
    text.DrawLatex(0.2, 0.90-down, "gaus fit [%.0f,%.0f]" %(fxmin3r, fxmax3r))
    text.DrawLatex(0.2, 0.87-down, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    widths[-1].append(gaus.GetParameter(2)/gaus.GetParameter(1))
    text.DrawLatex(0.2, 0.84-down, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gPad.Print("plots/breg_FatHmassRegDJvsFJ_1.png")

    h2.Fit("gaus", "I", "", fxmin2, fxmax2)
    h3.Fit("gaus", "I", "", fxmin3, fxmax3)
    h3.SetMaximum(ymax)
    h3.SetXTitle("H(125) mass [GeV]"); h3.Draw()
    h2.Draw("same")
    h3.Draw("same")

    leg = TLegend(0.58, 0.78, 0.92, 0.92)
    leg.AddEntry(h2, "Dijet-system")
    leg.AddEntry(h3, "Filterjet-system")
    leg.Draw()
    text.SetTextColor(kBlack)
    text.DrawLatex(0.1, 0.97, "130 < H p_{T} #leq 150 GeV")

    gaus = h2.GetFunction("gaus")
    text.SetTextColor(kDCyan)
    text.DrawLatex(0.2, 0.90, "gaus fit [%.0f,%.0f]" %(fxmin2, fxmax2))
    text.DrawLatex(0.2, 0.87, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    widths[-1].append(gaus.GetParameter(2)/gaus.GetParameter(1))
    text.DrawLatex(0.2, 0.84, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gaus = h3.GetFunction("gaus")
    text.SetTextColor(kDRed)
    text.DrawLatex(0.2, 0.90-down, "gaus fit [%.0f,%.0f]" %(fxmin3, fxmax3))
    text.DrawLatex(0.2, 0.87-down, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    widths[-1].append(gaus.GetParameter(2)/gaus.GetParameter(1))
    text.DrawLatex(0.2, 0.84-down, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gPad.Print("plots/breg_FatHmassNormDJvsFJ_1.png")

    #-------------------------------------------------------------------------------
    h2 = hmassnorm_pt2
    h3 = fathmassnorm_pt2
    h2r = hmassreg_pt2
    h3r = fathmassreg_pt2
    fxmin2r, fxmax2r = 105., 140.
    fxmin3r, fxmax3r = 80., 120.
    fxmin2, fxmax2 = 100, 140.
    fxmin3, fxmax3 = 75., 115.
    ymax = max(h2r.GetMaximum(), h3r.GetMaximum()) * 1.2
    h2r.Fit("gaus", "I", "", fxmin2r, fxmax2r)
    h3r.Fit("gaus", "I", "", fxmin3r, fxmax3r)
    h2.SetStats(0); h2.SetLineWidth(1); h2.SetMarkerSize(0); h2.SetLineColor(kDCyan)
    h3.SetStats(0); h3.SetLineWidth(1); h3.SetMarkerSize(0); h3.SetLineColor(kDRed)
    h2r.SetStats(0); h2r.SetLineWidth(2); h2r.SetMarkerSize(0); h2r.SetLineColor(kCyan)
    h3r.SetStats(0); h3r.SetLineWidth(2); h3r.SetMarkerSize(0); h3r.SetLineColor(kRed)
    h3r.SetMaximum(ymax)
    h3r.SetXTitle("H(125) mass [GeV]"); h3r.Draw()
    h2.Draw("same")
    h3.Draw("same")
    h2r.Draw("same")
    h3r.Draw("same")

    leg = TLegend(0.58, 0.78, 0.92, 0.92)
    leg.AddEntry(h2r, "Dijet-system")
    leg.AddEntry(h3r, "Filterjet-system")
    leg.Draw()
    text.SetTextColor(kBlack)
    text.DrawLatex(0.1, 0.97, "150 < H p_{T} #leq 170 GeV")

    widths.append([])
    gaus = h2r.GetFunction("gaus")
    text.SetTextColor(kCyan)
    text.DrawLatex(0.2, 0.90, "gaus fit [%.0f,%.0f]" %(fxmin2r, fxmax2r))
    text.DrawLatex(0.2, 0.87, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    widths[-1].append(gaus.GetParameter(2)/gaus.GetParameter(1))
    text.DrawLatex(0.2, 0.84, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gaus = h3r.GetFunction("gaus")
    text.SetTextColor(kRed)
    text.DrawLatex(0.2, 0.90-down, "gaus fit [%.0f,%.0f]" %(fxmin3r, fxmax3r))
    text.DrawLatex(0.2, 0.87-down, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    widths[-1].append(gaus.GetParameter(2)/gaus.GetParameter(1))
    text.DrawLatex(0.2, 0.84-down, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gPad.Print("plots/breg_FatHmassRegDJvsFJ_2.png")

    h2.Fit("gaus", "I", "", fxmin2, fxmax2)
    h3.Fit("gaus", "I", "", fxmin3, fxmax3)
    h3.SetMaximum(ymax)
    h3.SetXTitle("H(125) mass [GeV]"); h3.Draw()
    h2.Draw("same")
    h3.Draw("same")

    leg = TLegend(0.58, 0.78, 0.92, 0.92)
    leg.AddEntry(h2, "Dijet-system")
    leg.AddEntry(h3, "Filterjet-system")
    leg.Draw()
    text.SetTextColor(kBlack)
    text.DrawLatex(0.1, 0.97, "150 < H p_{T} #leq 170 GeV")

    gaus = h2.GetFunction("gaus")
    text.SetTextColor(kDCyan)
    text.DrawLatex(0.2, 0.90, "gaus fit [%.0f,%.0f]" %(fxmin2, fxmax2))
    text.DrawLatex(0.2, 0.87, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    widths[-1].append(gaus.GetParameter(2)/gaus.GetParameter(1))
    text.DrawLatex(0.2, 0.84, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gaus = h3.GetFunction("gaus")
    text.SetTextColor(kDRed)
    text.DrawLatex(0.2, 0.90-down, "gaus fit [%.0f,%.0f]" %(fxmin3, fxmax3))
    text.DrawLatex(0.2, 0.87-down, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    widths[-1].append(gaus.GetParameter(2)/gaus.GetParameter(1))
    text.DrawLatex(0.2, 0.84-down, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gPad.Print("plots/breg_FatHmassNormDJvsFJ_2.png")

    #-------------------------------------------------------------------------------
    h2 = hmassnorm_pt3
    h3 = fathmassnorm_pt3
    h2r = hmassreg_pt3
    h3r = fathmassreg_pt3
    fxmin2r, fxmax2r = 105., 140.
    fxmin3r, fxmax3r = 90., 130.
    fxmin2, fxmax2 = 100, 140.
    fxmin3, fxmax3 = 85., 125.
    ymax = max(h2r.GetMaximum(), h3r.GetMaximum()) * 1.2
    h2r.Fit("gaus", "I", "", fxmin2r, fxmax2r)
    h3r.Fit("gaus", "I", "", fxmin3r, fxmax3r)
    h2.SetStats(0); h2.SetLineWidth(1); h2.SetMarkerSize(0); h2.SetLineColor(kDCyan)
    h3.SetStats(0); h3.SetLineWidth(1); h3.SetMarkerSize(0); h3.SetLineColor(kDRed)
    h2r.SetStats(0); h2r.SetLineWidth(2); h2r.SetMarkerSize(0); h2r.SetLineColor(kCyan)
    h3r.SetStats(0); h3r.SetLineWidth(2); h3r.SetMarkerSize(0); h3r.SetLineColor(kRed)
    h3r.SetMaximum(ymax)
    h3r.SetXTitle("H(125) mass [GeV]"); h3r.Draw()
    h2.Draw("same")
    h3.Draw("same")
    h2r.Draw("same")
    h3r.Draw("same")

    leg = TLegend(0.58, 0.78, 0.92, 0.92)
    leg.AddEntry(h2r, "Dijet-system")
    leg.AddEntry(h3r, "Filterjet-system")
    leg.Draw()
    text.SetTextColor(kBlack)
    text.DrawLatex(0.1, 0.97, "170 < H p_{T} #leq 200 GeV")

    widths.append([])
    gaus = h2r.GetFunction("gaus")
    text.SetTextColor(kCyan)
    text.DrawLatex(0.2, 0.90, "gaus fit [%.0f,%.0f]" %(fxmin2r, fxmax2r))
    text.DrawLatex(0.2, 0.87, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    widths[-1].append(gaus.GetParameter(2)/gaus.GetParameter(1))
    text.DrawLatex(0.2, 0.84, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gaus = h3r.GetFunction("gaus")
    text.SetTextColor(kRed)
    text.DrawLatex(0.2, 0.90-down, "gaus fit [%.0f,%.0f]" %(fxmin3r, fxmax3r))
    text.DrawLatex(0.2, 0.87-down, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    widths[-1].append(gaus.GetParameter(2)/gaus.GetParameter(1))
    text.DrawLatex(0.2, 0.84-down, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gPad.Print("plots/breg_FatHmassRegDJvsFJ_3.png")

    h2.Fit("gaus", "I", "", fxmin2, fxmax2)
    h3.Fit("gaus", "I", "", fxmin3, fxmax3)
    h3.SetMaximum(ymax)
    h3.SetXTitle("H(125) mass [GeV]"); h3.Draw()
    h2.Draw("same")
    h3.Draw("same")

    leg = TLegend(0.58, 0.78, 0.92, 0.92)
    leg.AddEntry(h2, "Dijet-system")
    leg.AddEntry(h3, "Filterjet-system")
    leg.Draw()
    text.SetTextColor(kBlack)
    text.DrawLatex(0.1, 0.97, "170 < H p_{T} #leq 200 GeV")

    gaus = h2.GetFunction("gaus")
    text.SetTextColor(kDCyan)
    text.DrawLatex(0.2, 0.90, "gaus fit [%.0f,%.0f]" %(fxmin2, fxmax2))
    text.DrawLatex(0.2, 0.87, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    widths[-1].append(gaus.GetParameter(2)/gaus.GetParameter(1))
    text.DrawLatex(0.2, 0.84, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gaus = h3.GetFunction("gaus")
    text.SetTextColor(kDRed)
    text.DrawLatex(0.2, 0.90-down, "gaus fit [%.0f,%.0f]" %(fxmin3, fxmax3))
    text.DrawLatex(0.2, 0.87-down, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    widths[-1].append(gaus.GetParameter(2)/gaus.GetParameter(1))
    text.DrawLatex(0.2, 0.84-down, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gPad.Print("plots/breg_FatHmassNormDJvsFJ_3.png")

    #-------------------------------------------------------------------------------
    h2 = hmassnorm_pt4
    h3 = fathmassnorm_pt4
    h2r = hmassreg_pt4
    h3r = fathmassreg_pt4
    fxmin2r, fxmax2r = 105., 140.
    fxmin3r, fxmax3r = 100., 140.
    fxmin2, fxmax2 = 105, 140.
    fxmin3, fxmax3 = 95., 135.
    ymax = max(h2r.GetMaximum(), h3r.GetMaximum()) * 1.2
    h2r.Fit("gaus", "I", "", fxmin2r, fxmax2r)
    h3r.Fit("gaus", "I", "", fxmin3r, fxmax3r)
    h2.SetStats(0); h2.SetLineWidth(1); h2.SetMarkerSize(0); h2.SetLineColor(kDCyan)
    h3.SetStats(0); h3.SetLineWidth(1); h3.SetMarkerSize(0); h3.SetLineColor(kDRed)
    h2r.SetStats(0); h2r.SetLineWidth(2); h2r.SetMarkerSize(0); h2r.SetLineColor(kCyan)
    h3r.SetStats(0); h3r.SetLineWidth(2); h3r.SetMarkerSize(0); h3r.SetLineColor(kRed)
    h3r.SetMaximum(ymax)
    h3r.SetXTitle("H(125) mass [GeV]"); h3r.Draw()
    h2.Draw("same")
    h3.Draw("same")
    h2r.Draw("same")
    h3r.Draw("same")

    leg = TLegend(0.58, 0.78, 0.92, 0.92)
    leg.AddEntry(h2r, "Dijet-system")
    leg.AddEntry(h3r, "Filterjet-system")
    leg.Draw()
    text.SetTextColor(kBlack)
    text.DrawLatex(0.1, 0.97, "200 < H p_{T} #leq 250 GeV")

    widths.append([])
    gaus = h2r.GetFunction("gaus")
    text.SetTextColor(kCyan)
    text.DrawLatex(0.2, 0.90, "gaus fit [%.0f,%.0f]" %(fxmin2r, fxmax2r))
    text.DrawLatex(0.2, 0.87, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    widths[-1].append(gaus.GetParameter(2)/gaus.GetParameter(1))
    text.DrawLatex(0.2, 0.84, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gaus = h3r.GetFunction("gaus")
    text.SetTextColor(kRed)
    text.DrawLatex(0.2, 0.90-down, "gaus fit [%.0f,%.0f]" %(fxmin3r, fxmax3r))
    text.DrawLatex(0.2, 0.87-down, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    widths[-1].append(gaus.GetParameter(2)/gaus.GetParameter(1))
    text.DrawLatex(0.2, 0.84-down, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gPad.Print("plots/breg_FatHmassRegDJvsFJ_4.png")

    h2.Fit("gaus", "I", "", fxmin2, fxmax2)
    h3.Fit("gaus", "I", "", fxmin3, fxmax3)
    h3.SetMaximum(ymax)
    h3.SetXTitle("H(125) mass [GeV]"); h3.Draw()
    h2.Draw("same")
    h3.Draw("same")

    leg = TLegend(0.58, 0.78, 0.92, 0.92)
    leg.AddEntry(h2, "Dijet-system")
    leg.AddEntry(h3, "Filterjet-system")
    leg.Draw()
    text.SetTextColor(kBlack)
    text.DrawLatex(0.1, 0.97, "200 < H p_{T} #leq 250 GeV")

    gaus = h2.GetFunction("gaus")
    text.SetTextColor(kDCyan)
    text.DrawLatex(0.2, 0.90, "gaus fit [%.0f,%.0f]" %(fxmin2, fxmax2))
    text.DrawLatex(0.2, 0.87, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    widths[-1].append(gaus.GetParameter(2)/gaus.GetParameter(1))
    text.DrawLatex(0.2, 0.84, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gaus = h3.GetFunction("gaus")
    text.SetTextColor(kDRed)
    text.DrawLatex(0.2, 0.90-down, "gaus fit [%.0f,%.0f]" %(fxmin3, fxmax3))
    text.DrawLatex(0.2, 0.87-down, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    widths[-1].append(gaus.GetParameter(2)/gaus.GetParameter(1))
    text.DrawLatex(0.2, 0.84-down, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gPad.Print("plots/breg_FatHmassNormDJvsFJ_4.png")

    #-------------------------------------------------------------------------------
    h2 = hmassnorm_pt5
    h3 = fathmassnorm_pt5
    h2r = hmassreg_pt5
    h3r = fathmassreg_pt5
    fxmin2r, fxmax2r = 110., 145.
    fxmin3r, fxmax3r = 105., 145.
    fxmin2, fxmax2 = 110, 145.
    fxmin3, fxmax3 = 100., 135.
    ymax = max(h2r.GetMaximum(), h3r.GetMaximum()) * 1.2
    h2r.Fit("gaus", "I", "", fxmin2r, fxmax2r)
    h3r.Fit("gaus", "I", "", fxmin3r, fxmax3r)
    h2.SetStats(0); h2.SetLineWidth(1); h2.SetMarkerSize(0); h2.SetLineColor(kDCyan)
    h3.SetStats(0); h3.SetLineWidth(1); h3.SetMarkerSize(0); h3.SetLineColor(kDRed)
    h2r.SetStats(0); h2r.SetLineWidth(2); h2r.SetMarkerSize(0); h2r.SetLineColor(kCyan)
    h3r.SetStats(0); h3r.SetLineWidth(2); h3r.SetMarkerSize(0); h3r.SetLineColor(kRed)
    h3r.SetMaximum(ymax)
    h3r.SetXTitle("H(125) mass [GeV]"); h3r.Draw()
    h2.Draw("same")
    h3.Draw("same")
    h2r.Draw("same")
    h3r.Draw("same")

    leg = TLegend(0.58, 0.78, 0.92, 0.92)
    leg.AddEntry(h2r, "Dijet-system")
    leg.AddEntry(h3r, "Filterjet-system")
    leg.Draw()
    text.SetTextColor(kBlack)
    text.DrawLatex(0.1, 0.97, "250 < H p_{T} #leq 400")

    widths.append([])
    gaus = h2r.GetFunction("gaus")
    text.SetTextColor(kCyan)
    text.DrawLatex(0.2, 0.90, "gaus fit [%.0f,%.0f]" %(fxmin2r, fxmax2r))
    text.DrawLatex(0.2, 0.87, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    widths[-1].append(gaus.GetParameter(2)/gaus.GetParameter(1))
    text.DrawLatex(0.2, 0.84, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gaus = h3r.GetFunction("gaus")
    text.SetTextColor(kRed)
    text.DrawLatex(0.2, 0.90-down, "gaus fit [%.0f,%.0f]" %(fxmin3r, fxmax3r))
    text.DrawLatex(0.2, 0.87-down, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    widths[-1].append(gaus.GetParameter(2)/gaus.GetParameter(1))
    text.DrawLatex(0.2, 0.84-down, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gPad.Print("plots/breg_FatHmassRegDJvsFJ_5.png")

    h2.Fit("gaus", "I", "", fxmin2, fxmax2)
    h3.Fit("gaus", "I", "", fxmin3, fxmax3)
    h3.SetMaximum(ymax)
    h3.SetXTitle("H(125) mass [GeV]"); h3.Draw()
    h2.Draw("same")
    h3.Draw("same")

    leg = TLegend(0.58, 0.78, 0.92, 0.92)
    leg.AddEntry(h2, "Dijet-system")
    leg.AddEntry(h3, "Filterjet-system")
    leg.Draw()
    text.SetTextColor(kBlack)
    text.DrawLatex(0.1, 0.97, "250 < H p_{T} #leq 400")

    gaus = h2.GetFunction("gaus")
    text.SetTextColor(kDCyan)
    text.DrawLatex(0.2, 0.90, "gaus fit [%.0f,%.0f]" %(fxmin2, fxmax2))
    text.DrawLatex(0.2, 0.87, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    widths[-1].append(gaus.GetParameter(2)/gaus.GetParameter(1))
    text.DrawLatex(0.2, 0.84, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gaus = h3.GetFunction("gaus")
    text.SetTextColor(kDRed)
    text.DrawLatex(0.2, 0.90-down, "gaus fit [%.0f,%.0f]" %(fxmin3, fxmax3))
    text.DrawLatex(0.2, 0.87-down, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    widths[-1].append(gaus.GetParameter(2)/gaus.GetParameter(1))
    text.DrawLatex(0.2, 0.84-down, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gPad.Print("plots/breg_FatHmassNormDJvsFJ_5.png")

    #-------------------------------------------------------------------------------
    h2 = hmassnorm_pt6
    h3 = fathmassnorm_pt6
    h2r = hmassreg_pt6
    h3r = fathmassreg_pt6
    fxmin2r, fxmax2r = 120., 150.
    fxmin3r, fxmax3r = 110., 145.
    fxmin2, fxmax2 = 120, 150.
    fxmin3, fxmax3 = 105., 140.
    ymax = max(h2r.GetMaximum(), h3r.GetMaximum()) * 1.2
    h2r.Fit("gaus", "I", "", fxmin2r, fxmax2r)
    h3r.Fit("gaus", "I", "", fxmin3r, fxmax3r)
    h2.SetStats(0); h2.SetLineWidth(1); h2.SetMarkerSize(0); h2.SetLineColor(kDCyan)
    h3.SetStats(0); h3.SetLineWidth(1); h3.SetMarkerSize(0); h3.SetLineColor(kDRed)
    h2r.SetStats(0); h2r.SetLineWidth(2); h2r.SetMarkerSize(0); h2r.SetLineColor(kCyan)
    h3r.SetStats(0); h3r.SetLineWidth(2); h3r.SetMarkerSize(0); h3r.SetLineColor(kRed)
    h3r.SetMaximum(ymax)
    h3r.SetXTitle("H(125) mass [GeV]"); h3r.Draw()
    h2.Draw("same")
    h3.Draw("same")
    h2r.Draw("same")
    h3r.Draw("same")

    leg = TLegend(0.58, 0.78, 0.92, 0.92)
    leg.AddEntry(h2r, "Dijet-system")
    leg.AddEntry(h3r, "Filterjet-system")
    leg.Draw()
    text.SetTextColor(kBlack)
    text.DrawLatex(0.1, 0.97, "400 < H p_{T}")

    widths.append([])
    gaus = h2r.GetFunction("gaus")
    text.SetTextColor(kCyan)
    text.DrawLatex(0.2, 0.90, "gaus fit [%.0f,%.0f]" %(fxmin2r, fxmax2r))
    text.DrawLatex(0.2, 0.87, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    widths[-1].append(gaus.GetParameter(2)/gaus.GetParameter(1))
    text.DrawLatex(0.2, 0.84, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gaus = h3r.GetFunction("gaus")
    text.SetTextColor(kRed)
    text.DrawLatex(0.2, 0.90-down, "gaus fit [%.0f,%.0f]" %(fxmin3r, fxmax3r))
    text.DrawLatex(0.2, 0.87-down, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    widths[-1].append(gaus.GetParameter(2)/gaus.GetParameter(1))
    text.DrawLatex(0.2, 0.84-down, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gPad.Print("plots/breg_FatHmassRegDJvsFJ_6.png")

    h2.Fit("gaus", "I", "", fxmin2, fxmax2)
    h3.Fit("gaus", "I", "", fxmin3, fxmax3)
    h3.SetMaximum(ymax)
    h3.SetXTitle("H(125) mass [GeV]"); h3.Draw()
    h2.Draw("same")
    h3.Draw("same")

    leg = TLegend(0.58, 0.78, 0.92, 0.92)
    leg.AddEntry(h2, "Dijet-system")
    leg.AddEntry(h3, "Filterjet-system")
    leg.Draw()
    text.SetTextColor(kBlack)
    text.DrawLatex(0.1, 0.97, "400 < H p_{T}")

    gaus = h2.GetFunction("gaus")
    text.SetTextColor(kDCyan)
    text.DrawLatex(0.2, 0.90, "gaus fit [%.0f,%.0f]" %(fxmin2, fxmax2))
    text.DrawLatex(0.2, 0.87, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    widths[-1].append(gaus.GetParameter(2)/gaus.GetParameter(1))
    text.DrawLatex(0.2, 0.84, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gaus = h3.GetFunction("gaus")
    text.SetTextColor(kDRed)
    text.DrawLatex(0.2, 0.90-down, "gaus fit [%.0f,%.0f]" %(fxmin3, fxmax3))
    text.DrawLatex(0.2, 0.87-down, "#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f" %(gaus.GetParameter(1), gaus.GetParError(1), gaus.GetParameter(2), gaus.GetParError(2)))
    widths[-1].append(gaus.GetParameter(2)/gaus.GetParameter(1))
    text.DrawLatex(0.2, 0.84-down, "#chi^{2}/NDF = %.1f / %.0f" %(gaus.GetChisquare(), gaus.GetNDF()))
    gPad.Print("plots/breg_FatHmassNormDJvsFJ_6.png")

    print widths

    xbins = array.array('f', [130.,150.,170.,200.,250.,400.,700.])
    hh0 = TH1F("hh0", "", len(xbins)-1, xbins)
    hh1 = TH1F("hh1", "", len(xbins)-1, xbins)
    hh2 = TH1F("hh2", "", len(xbins)-1, xbins)
    hh3 = TH1F("hh3", "", len(xbins)-1, xbins)
    for i in xrange(len(widths)):
        hh0.SetBinContent(1+i, widths[i][0]*100.)
        hh1.SetBinContent(1+i, widths[i][1]*100.)
        hh2.SetBinContent(1+i, widths[i][2]*100.)
        hh3.SetBinContent(1+i, widths[i][3]*100.)

    hh0.SetStats(0); hh0.SetLineWidth(1); hh0.SetMarkerSize(1); hh0.SetLineColor(kCyan); hh0.SetMarkerColor(kCyan)
    hh1.SetStats(0); hh1.SetLineWidth(1); hh1.SetMarkerSize(1); hh1.SetLineColor(kRed); hh1.SetMarkerColor(kRed)
    hh2.SetStats(0); hh2.SetLineWidth(2); hh2.SetMarkerSize(1); hh2.SetLineColor(kDCyan); hh2.SetMarkerColor(kDCyan)
    hh3.SetStats(0); hh3.SetLineWidth(2); hh3.SetMarkerSize(1); hh3.SetLineColor(kDRed); hh3.SetMarkerColor(kDRed)
    hh0.GetYaxis().SetTitle("res = width/mean [%]"); hh0.GetXaxis().SetTitle("H(125) p_{T} [GeV]"); hh0.SetMinimum(0.); hh0.SetMaximum(25.); hh0.Draw("LP")
    hh3.Draw("LPsame")
    hh2.Draw("LPsame")
    hh1.Draw("LPsame")
    hh0.Draw("LPsame")
    
    leg = TLegend(0.48, 0.78, 0.92, 0.92)
    leg.AddEntry(hh2, "DJ before regression")
    leg.AddEntry(hh3, "FJ before regression")
    leg.AddEntry(hh0, "DJ after regression")
    leg.AddEntry(hh1, "FJ after regression")
    leg.Draw()
    text.SetTextColor(kBlack)
    gPad.Print("plots/breg_FatHmassResolution.png")
    
