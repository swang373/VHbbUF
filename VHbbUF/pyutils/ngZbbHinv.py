from ROOT import gROOT, gPad
gROOT.LoadMacro("HelperFunctions.h")
gROOT.LoadMacro("tdrstyle.C")
gROOT.ProcessLine("setTDRStyle()")
from ROOT import TFile, TTree, TF1, TH1F, TLatex, TLegend, TColor, TCanvas, TPad, TLatex, TLine
from math import sqrt

#gROOT.SetBatch(1)


#_______________________________________________________________________________
def setRatioStyle(ratio):
    ratio.SetStats(0)
    ratio.GetXaxis().SetLabelSize(0.12)
    ratio.GetXaxis().SetTitleSize(0.14)
    ratio.GetXaxis().SetTitleOffset(1.10)
    ratio.GetYaxis().SetLabelSize(0.10)
    #ratio.GetYaxis().SetTitleSize(0.12)
    ratio.GetYaxis().SetTitleSize(0.10)
    ratio.GetYaxis().SetTitleOffset(0.6)
    ratio.GetYaxis().SetNdivisions(505)

#_______________________________________________________________________________
def drawCanvas():
    c1 = TCanvas("c1", "c1", 700, 700)
    c1.cd()
    return c1

#_______________________________________________________________________________
def drawCanvas2():
    c1 = TCanvas("c1", "c1", 700, 700)
    pad1 = TPad("pad1", "top pad"   , 0.0, 0.3, 1.0, 1.0)
    pad1.SetBottomMargin(0.0)
    pad1.SetNumber(1)
    pad1.Draw()
    pad1.cd()
    c1.cd()
    pad2 = TPad("pad2", "bottom pad", 0.0, 0.0, 1.0, 0.3)
    pad2.SetTopMargin(0.0)
    pad2.SetBottomMargin(0.35)
    pad2.SetNumber(2)
    pad2.Draw()
    pad2.cd()
    c1.cd()
    return c1

#_______________________________________________________________________________
def drawLatex():
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextAlign(12)
    latex.SetTextFont(42)
    latex.SetTextSize(0.032)
    return latex

def drawLegend(x1=0.74, y1=0.56, x2=0.92, y2=0.92):
    leg = TLegend(x1, y1, x2, y2)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    #leg.SetFillColor(0)
    #leg.SetLineColor(0)
    #leg.SetShadowColor(0)
    leg.SetTextFont(42)
    #leg.SetTextSize(0.03)
    return leg
    

#_______________________________________________________________________________
selection = "Vtype==4 && METtype1corr.et>100 && hJet_pt[0]>60 && hJet_pt[1]>30 && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5 && hJet_csv_nominal[0]>0.244 && hJet_csv_nominal[1]>0.244 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && H.mass>0"

class Foo:
    pass

def drawRegression(foo):
    c1 = drawCanvas()
    tree1 = TFile.Open(foo.infile1).tree
    tree2 = TFile.Open(foo.infile2).tree
    #h1 = TH1F("h1", "; m(jj) [GeV]", 100, 0, 250)
    h1 = TH1F("h1", "; m(jj) [GeV]; Events / 2.0", 60, 60, 180)
    h1.SetLineWidth(2)
    h1.SetStats(0)
    h2 = h1.Clone("h2")
    hr1 = h1.Clone("hr1")
    hr2 = h1.Clone("hr2")
    
    tree1.Project("h1", "H.mass", selection)
    tree2.Project("h2", "H.mass", selection)
    tree1.Project("hr1", "HmassReg", selection.replace("_pt", "_ptReg"))
    tree2.Project("hr2", "HmassReg", selection.replace("_pt", "_ptReg"))
    
    h1.SetMaximum(h1.GetMaximum() * 1.3)
    h1.Draw()
    #h2.Draw("same")
    hr1.SetLineColor(4)
    hr1.Draw("same")
    
    latex = drawLatex()
    latex.SetTextFont(62)
    latex.SetTextSize(0.052*0.7)
    latex.DrawLatex(0.64, 0.89, "CMS Simulation")
    latex.SetTextSize(0.04*0.7)
    latex.DrawLatex(0.64, 0.85, "#sqrt{s} = 8 TeV")
    
    leg = drawLegend(0.64, y1=0.66, x2=0.94, y2=0.82)
    leg.AddEntry(h1, foo.leg+" nominal" , "l")
    leg.AddEntry(hr1, foo.leg+" regressed", "l")
    leg.Draw()
    
    #gPad.Print("regression_"+foo.name+".png")
    #gPad.Print("regression_"+foo.name+".pdf")
    h1.Fit("gaus", "", "", foo.fitrange[0], foo.fitrange[1])
    hr1.Fit("gaus", "", "", foo.fitrange[0], foo.fitrange[1])
    
    h1.Draw()
    hr2.SetLineColor(6)
    hr2.Draw("same")
    
    latex = drawLatex()
    latex.SetTextFont(62)
    latex.SetTextSize(0.052*0.7)
    latex.DrawLatex(0.64, 0.89, "CMS Simulation")
    latex.SetTextSize(0.04*0.7)
    latex.DrawLatex(0.64, 0.85, "#sqrt{s} = 8 TeV")
    
    leg = drawLegend(0.64, y1=0.66, x2=0.94, y2=0.82)
    leg.AddEntry(h1, foo.leg+" nominal" , "l")
    leg.AddEntry(hr2, foo.leg+" regressed", "l")
    leg.Draw()
    
    #gPad.Print("regression_star_"+foo.name+".png")
    #gPad.Print("regression_star_"+foo.name+".pdf")
    hr2.Fit("gaus", "", "", foo.fitrange[0], foo.fitrange[1])

foo = Foo()
foo.selection = selection
foo.infile1 = "Step3_20130314/Step3_ZbbHinv125.root"
foo.infile2 = "skim/Step3_ZbbHinv125.root"
foo.name = "ZbbHinv"
foo.leg = "Z(bb)H(inv)"
foo.fitrange = (70, 115)
drawRegression(foo)

foo = Foo()
foo.selection = selection
foo.infile1 = "Step3_20130314/Step3_ZnnH125.root"
foo.infile2 = "skim/Step3_ZnnH125.root"
foo.name = "ZnnH"
foo.leg = "Z(#nu#nu)H(bb)"
foo.fitrange = (105, 145)
drawRegression(foo)

foo = Foo()
foo.selection = selection + " && abs(eventFlav)==5"
foo.infile1 = "Step3_20130314/Step3_ZZ.root"
foo.infile2 = "skim/Step3_ZZ.root"
foo.name = "ZZ"
foo.leg = "ZZ(bb)"
foo.fitrange = (70, 115)
#drawRegression(foo)

#_______________________________________________________________________________
foo = Foo()
foo.selection = selection
foo.infile1 = "Step3_20130314/Step3_ZbbHinv125.root"
foo.infile2 = "Step3_20130314/Step3_ZZ.root"
if False:
    c1 = drawCanvas()
    tree1 = TFile.Open(foo.infile1).tree
    tree2 = TFile.Open(foo.infile2).tree
    h1 = TH1F("h1", "; gen Z p_{T} [GeV]; (normalized)", 100, 80, 280)
    h1.SetLineWidth(2)
    h1.SetStats(0)
    h2 = h1.Clone("h2")
    tree1.Project("h1", "genH.pt", selection)
    tree2.Project("h2", "genZ.pt", selection + " && abs(eventFlav)==5")
    h1.Scale(1.0/h1.Integral(11,100))
    h2.Scale(1.0/h2.Integral(11,100))
    kPink = 900
    h1.SetLineColor(kPink -4)
    h2.SetLineColor(920)
    h2.SetFillColor(920)
    h2.SetMaximum(h2.GetMaximum() * 1.2)
    h2.Draw()
    h1.Draw("same")

    latex = drawLatex()
    latex.SetTextFont(62)
    latex.SetTextSize(0.052*0.7)
    latex.DrawLatex(0.64, 0.89, "CMS Simulation")
    latex.SetTextSize(0.04*0.7)
    latex.DrawLatex(0.64, 0.85, "#sqrt{s} = 8 TeV")
    
    leg = drawLegend(0.64, y1=0.74, x2=0.94, y2=0.82)
    leg.AddEntry(h1, "Z(bb)H(inv)", "l")
    leg.AddEntry(h2, "ZZ(bb)" , "lf")
    leg.Draw()
