from ROOT import TFile, TH1F, TH1, TLatex, gPad, gROOT, gStyle
import copy

file1 = "Step4_20130812/Step4_ZbbHinv125.root"
file2 = "skim/Step4_ZbbHinvZmmToZbb.root"
wait = True
#gROOT.SetBatch(1)

varexp = "BDTinvregular_125[0]"
selection = "selectFlags[0][0]"

#_______________________________________________________________________________
gROOT.LoadMacro("HelperFunctions.h")
gROOT.LoadMacro("tdrstyle.C")
gROOT.ProcessLine("setTDRStyle()")
TH1.SetDefaultSumw2()

gStyle.SetMarkerSize(0.7)

latex = TLatex()
latex.SetNDC()
latex.SetTextAlign(12)
latex.SetTextFont(42)
latex.SetTextSize(0.035)

tfile1 = TFile.Open(file1)
tfile2 = TFile.Open(file2)

hname = "BDT"
etc = ("; BDT; (normalized)", 100, -1, 0.4)

#_______________________________________________________________________________
for ih in xrange(3):
    h1 = TH1F("h1_%i" % ih, etc[0], etc[1], etc[2], etc[3])
    h2 = TH1F("h2_%i" % ih, etc[0], etc[1], etc[2], etc[3])
    h3 = TH1F("h3_%i" % ih, etc[0], etc[1], etc[2], etc[3])
    h1.SetLineWidth(2)
    h1.SetLineColor(1)
    h2.SetLineWidth(2)
    h2.SetLineColor(4)
    h3.SetLineWidth(2)
    h3.SetLineColor(6)

    if ih == 0:
        t1 = tfile1.tree_ZnunuHighPt_train
        t2 = tfile2.tree_ZnunuHighPt_train
    elif ih == 1:
        t1 = tfile1.tree_ZnunuMedPt_train
        t2 = tfile2.tree_ZnunuMedPt_train
    elif ih == 2:
        t1 = tfile1.tree_ZnunuLowPt_train
        t2 = tfile2.tree_ZnunuLowPt_train
    
    t1.Project("h1_%i" % ih, varexp, selection)
    t2.Project("h2_%i" % ih, varexp, selection)
    
    if ih == 0:
        t1 = tfile1.tree_ZnunuHighPt_test
        t2 = tfile2.tree_ZnunuHighPt_test
    elif ih == 1:
        t1 = tfile1.tree_ZnunuMedPt_test
        t2 = tfile2.tree_ZnunuMedPt_test
    elif ih == 2:
        t1 = tfile1.tree_ZnunuLowPt_test
        t2 = tfile2.tree_ZnunuLowPt_test
    
    t1.Project("+h1_%i" % ih, varexp, selection)
    t2.Project("+h2_%i" % ih, varexp, selection)

    h1.Scale(1.0/h1.Integral())
    h2.Scale(1.0/h2.Integral())
    h1.SetMaximum(h1.GetMaximum()*1.4)
    h1.Draw("hist")
    h2.Draw("samehist")
    #h3.Draw("samehist")
    
    latex.DrawLatex(0.19, 0.90, "Full BDT selection")
    latex.DrawLatex(0.19, 0.85, "#color[1]{FullSim}")
    latex.DrawLatex(0.19, 0.80, "#color[4]{Smear/Histo}")
    latex.DrawLatex(0.42, 0.85, "#color[1]{#mu=%.3f, #sigma=%.3f}" % (h1.GetMean(), h1.GetRMS()))
    latex.DrawLatex(0.42, 0.80, "#color[4]{#mu=%.3f, #sigma=%.3f}" % (h2.GetMean(), h2.GetRMS()))
    kolprob = h1.KolmogorovTest(h2)
    latex.DrawLatex(0.19, 0.75, "#color[2]{K_{s} = %.3f}" % kolprob)
    gPad.Print("plots_zmmtozbb/zmmtozbb_compare_ZbbHinv_%s_%i.png" % (hname, ih))
    gPad.Print("plots_zmmtozbb/zmmtozbb_compare_ZbbHinv_%s_%i.pdf" % (hname, ih))
    
    if wait:  raw_input("Press Enter to continue...")
