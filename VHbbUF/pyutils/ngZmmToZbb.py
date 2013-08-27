from ROOT import TFile, TH1F, TH1, TLatex, gPad, gROOT, gStyle
from math import sqrt

infilenames = [
"/uscms_data/d3/degrutto/CMSSW_5_3_3_patch2/src/VHbbAnalysis/VHbbDataFormats/bin/ZH_ZToBB_HToInv_M-105_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1.root",
"/uscms_data/d3/degrutto/CMSSW_5_3_3_patch2/src/VHbbAnalysis/VHbbDataFormats/bin/ZH_ZToBB_HToInv_M-115_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1.root",
"/uscms_data/d3/degrutto/CMSSW_5_3_3_patch2/src/VHbbAnalysis/VHbbDataFormats/bin/ZH_ZToBB_HToInv_M-125_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1.root",
"/uscms_data/d3/degrutto/CMSSW_5_3_3_patch2/src/VHbbAnalysis/VHbbDataFormats/bin/ZH_ZToBB_HToInv_M-135_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1.root",
"/uscms_data/d3/degrutto/CMSSW_5_3_3_patch2/src/VHbbAnalysis/VHbbDataFormats/bin/ZH_ZToBB_HToInv_M-145_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1.root",
"/uscms_data/d3/degrutto/CMSSW_5_3_3_patch2/src/VHbbAnalysis/VHbbDataFormats/bin/ZH_ZToBB_HToInv_M-150_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1.root",
]

outfilename = "zmmtozbb_jetres.root"
useGenJet = False
useParton = not useGenJet


genzptbins = [
    (0,1000),
]
genptbins = [
    (20,30),
    (30,50),
    (50,70),
    (70,100),
    (100,160),
    (160,320),
    (320,640),
]
#genptbins = [
#    (30,60),
#    (60,100),
#    (100,200),
#    (200,2000),
#]
etabins = [
    (0.0,0.5),
    (0.5,1.0),
    (1.0,1.5),
    (1.5,2.5),
    (2.5,5.0),
]


gROOT.LoadMacro("HelperFunctions.h")
gROOT.LoadMacro("tdrstyle.C")
gROOT.ProcessLine("setTDRStyle()")
TH1.SetDefaultSumw2()

latex = TLatex()
latex.SetNDC()
latex.SetTextAlign(12)
latex.SetTextFont(42)
latex.SetTextSize(0.035)

gStyle.SetMarkerSize(0.7)

projappend_ = False
def projname(name):
    if projappend_:
        return "+"+name
    return name


histos1 = []
histos2 = []
hlabels = []

print "            norm   mean   RMS  "

nfiles = 0
for infilename in infilenames:
    infile = TFile.Open(infilename)
    intree = infile.tree
    nfiles += 1
    
    for etabin in etabins:
        for genptbin in genptbins:
            for genzptbin in genzptbins:
                
                hname1 = "h1_%i" % len(histos1)
                hname2 = "h2_%i" % len(histos2)
                hlabel = "hl_%i" % len(hlabels)
                #h1 = TH1F(hname1, "; p_{T}^{reco} / p_{T}^{gen}; Events / (0.02)", 150, 0, 3.0)
                #h2 = TH1F(hname2, "; p_{T}^{reco} / p_{T}^{gen}; Events / (0.02)", 150, 0, 3.0)
                h1 = TH1F(hname1, "; p_{T}^{reco} / p_{T}^{parton}; Events / (0.02)", 150, 0, 3.0)
                h2 = TH1F(hname2, "; p_{T}^{reco} / p_{T}^{parton}; Events / (0.02)", 150, 0, 3.0)
                hl = TH1F(hlabel, "; label", 10, 0, 10)
                h1.SetLineColor(1)
                h2.SetLineColor(4)
                h2.SetLineWidth(2)
                
                # N.B. genZ.pt is zero in Pythia sample
                
                if useGenJet:
                    projappend_ = False
                    intree.Project(projname(hname1), "hJet_pt[0] / hJet_genPt[0]", "hJet_id[0]==1 && hJet_puJetIdL[0]>0 && hJet_csv[0]>0 && hJet_csv_nominal[0]>0.244 && deltaR(hJet_eta[0],hJet_phi[0],hJet_genEta[0],hJet_genPhi[0])<0.3 && (%f<genH.pt && genH.pt<%f) && (%f<hJet_genPt[0] && hJet_genPt[0]<%f) && (%f<abs(hJet_genEta[0]) && abs(hJet_genEta[0])<%f)" % (genzptbin[0], genzptbin[1], genptbin[0], genptbin[1], etabin[0], etabin[1]))
                    print "%-10s  %6.0f  %.3f  %.3f" % (hname1, h1.Integral(0,1000), h1.GetMean(), h1.GetRMS())
                    projappend_ = True
                    intree.Project(projname(hname1), "hJet_pt[1] / hJet_genPt[1]", "hJet_id[1]==1 && hJet_puJetIdL[1]>0 && hJet_csv[1]>0 && hJet_csv_nominal[1]>0.244 && deltaR(hJet_eta[1],hJet_phi[1],hJet_genEta[1],hJet_genPhi[1])<0.3 && (%f<genH.pt && genH.pt<%f) && (%f<hJet_genPt[1] && hJet_genPt[1]<%f) && (%f<abs(hJet_genEta[1]) && abs(hJet_genEta[1])<%f)" % (genzptbin[0], genzptbin[1], genptbin[0], genptbin[1], etabin[0], etabin[1]))
                    print "%-10s  %6.0f  %.3f  %.3f" % (hname1, h1.Integral(0,1000), h1.GetMean(), h1.GetRMS())
                    
                    projappend_ = False
                    intree.Project(projname(hname2), "smear_gen_to_reco(hJet_pt[0], hJet_genPt[0], hJet_genEta[0]) / hJet_genPt[0]", "hJet_id[0]==1 && hJet_puJetIdL[0]>0 && hJet_csv[0]>0 && hJet_csv_nominal[0]>0.244 && deltaR(hJet_eta[0],hJet_phi[0],hJet_genEta[0],hJet_genPhi[0])<0.3 && (%f<genH.pt && genH.pt<%f) && (%f<hJet_genPt[0] && hJet_genPt[0]<%f) && (%f<abs(hJet_genEta[0]) && abs(hJet_genEta[0])<%f)" % (genzptbin[0], genzptbin[1], genptbin[0], genptbin[1], etabin[0], etabin[1]))
                    print "%-10s  %6.0f  %.3f  %.3f" % (hname2, h2.Integral(0,1000), h2.GetMean(), h2.GetRMS())
                    projappend_ = True
                    intree.Project(projname(hname2), "smear_gen_to_reco(hJet_pt[1], hJet_genPt[1], hJet_genEta[1]) / hJet_genPt[1]", "hJet_id[1]==1 && hJet_puJetIdL[1]>0 && hJet_csv[1]>0 && hJet_csv_nominal[1]>0.244 && deltaR(hJet_eta[1],hJet_phi[1],hJet_genEta[1],hJet_genPhi[1])<0.3 && (%f<genH.pt && genH.pt<%f) && (%f<hJet_genPt[1] && hJet_genPt[1]<%f) && (%f<abs(hJet_genEta[1]) && abs(hJet_genEta[1])<%f)" % (genzptbin[0], genzptbin[1], genptbin[0], genptbin[1], etabin[0], etabin[1]))
                    print "%-10s  %6.0f  %.3f  %.3f" % (hname2, h2.Integral(0,1000), h2.GetMean(), h2.GetRMS())
                
                elif useParton:
                    # reco jet 1
                    projappend_ = False
                    intree.Project(projname(hname1), "hJet_pt[0] / genB.pt", "hJet_id[0]==1 && hJet_puJetIdL[0]>0 && hJet_csv[0]>0 && hJet_csv_nominal[0]>0.244 && hJet_flavour[0]==5 && deltaR(hJet_genEta[0],hJet_genPhi[0],genB.eta,genB.phi)<0.1 && (%f<genH.pt && genH.pt<%f) && (%f<genB.pt && genB.pt<%f) && (%f<abs(genB.eta) && abs(genB.eta)<%f)" % (genzptbin[0], genzptbin[1], genptbin[0], genptbin[1], etabin[0], etabin[1]))
                    print "%-10s  %6.0f  %.3f  %.3f" % (hname1, h1.Integral(0,1000), h1.GetMean(), h1.GetRMS())
                    projappend_ = True
                    intree.Project(projname(hname1), "hJet_pt[0] / genBbar.pt", "hJet_id[0]==1 && hJet_puJetIdL[0]>0 && hJet_csv[0]>0 && hJet_csv_nominal[0]>0.244 && hJet_flavour[0]==-5 && deltaR(hJet_genEta[0],hJet_genPhi[0],genBbar.eta,genBbar.phi)<0.1 && (%f<genH.pt && genH.pt<%f) && (%f<genBbar.pt && genBbar.pt<%f) && (%f<abs(genBbar.eta) && abs(genBbar.eta)<%f)" % (genzptbin[0], genzptbin[1], genptbin[0], genptbin[1], etabin[0], etabin[1]))
                    print "%-10s  %6.0f  %.3f  %.3f" % (hname1, h1.Integral(0,1000), h1.GetMean(), h1.GetRMS())
                    # reco jet 2
                    projappend_ = True
                    intree.Project(projname(hname1), "hJet_pt[1] / genB.pt", "hJet_id[1]==1 && hJet_puJetIdL[1]>0 && hJet_csv[1]>0 && hJet_csv_nominal[1]>0.244 && hJet_flavour[1]==5 && deltaR(hJet_genEta[1],hJet_genPhi[1],genB.eta,genB.phi)<0.1 && (%f<genH.pt && genH.pt<%f) && (%f<genB.pt && genB.pt<%f) && (%f<abs(genB.eta) && abs(genB.eta)<%f)" % (genzptbin[0], genzptbin[1], genptbin[0], genptbin[1], etabin[0], etabin[1]))
                    print "%-10s  %6.0f  %.3f  %.3f" % (hname1, h1.Integral(0,1000), h1.GetMean(), h1.GetRMS())
                    projappend_ = True
                    intree.Project(projname(hname1), "hJet_pt[1] / genBbar.pt", "hJet_id[1]==1 && hJet_puJetIdL[1]>0 && hJet_csv[1]>0 && hJet_csv_nominal[1]>0.244 && hJet_flavour[1]==-5 && deltaR(hJet_genEta[1],hJet_genPhi[1],genBbar.eta,genBbar.phi)<0.1 && (%f<genH.pt && genH.pt<%f) && (%f<genBbar.pt && genBbar.pt<%f) && (%f<abs(genBbar.eta) && abs(genBbar.eta)<%f)" % (genzptbin[0], genzptbin[1], genptbin[0], genptbin[1], etabin[0], etabin[1]))
                    print "%-10s  %6.0f  %.3f  %.3f" % (hname1, h1.Integral(0,1000), h1.GetMean(), h1.GetRMS())
                    
                    # smear jet 1
                    projappend_ = False
                    intree.Project(projname(hname2), "smear_gen_to_reco(hJet_pt[0], genB.pt, genB.eta) / genB.pt", "hJet_id[0]==1 && hJet_puJetIdL[0]>0 && hJet_csv[0]>0 && hJet_csv_nominal[0]>0.244 && hJet_flavour[0]==5 && deltaR(hJet_genEta[0],hJet_genPhi[0],genB.eta,genB.phi)<0.1 && (%f<genH.pt && genH.pt<%f) && (%f<genB.pt && genB.pt<%f) && (%f<abs(genB.eta) && abs(genB.eta)<%f)" % (genzptbin[0], genzptbin[1], genptbin[0], genptbin[1], etabin[0], etabin[1]))
                    print "%-10s  %6.0f  %.3f  %.3f" % (hname2, h2.Integral(0,1000), h2.GetMean(), h2.GetRMS())
                    projappend_ = True
                    intree.Project(projname(hname2), "smear_gen_to_reco(hJet_pt[0], genBbar.pt, genBbar.eta) / genBbar.pt", "hJet_id[0]==1 && hJet_puJetIdL[0]>0 && hJet_csv[0]>0 && hJet_csv_nominal[0]>0.244 && hJet_flavour[0]==-5 && deltaR(hJet_genEta[0],hJet_genPhi[0],genBbar.eta,genBbar.phi)<0.1 && (%f<genH.pt && genH.pt<%f) && (%f<genBbar.pt && genBbar.pt<%f) && (%f<abs(genBbar.eta) && abs(genBbar.eta)<%f)" % (genzptbin[0], genzptbin[1], genptbin[0], genptbin[1], etabin[0], etabin[1]))
                    print "%-10s  %6.0f  %.3f  %.3f" % (hname2, h2.Integral(0,1000), h2.GetMean(), h2.GetRMS())
                    # smear jet 2
                    projappend_ = True
                    intree.Project(projname(hname2), "smear_gen_to_reco(hJet_pt[1], genB.pt, genB.eta) / genB.pt", "hJet_id[1]==1 && hJet_puJetIdL[1]>0 && hJet_csv[1]>0 && hJet_csv_nominal[1]>0.244 && hJet_flavour[1]==5 && deltaR(hJet_genEta[1],hJet_genPhi[1],genB.eta,genB.phi)<0.1 && (%f<genH.pt && genH.pt<%f) && (%f<genB.pt && genB.pt<%f) && (%f<abs(genB.eta) && abs(genB.eta)<%f)" % (genzptbin[0], genzptbin[1], genptbin[0], genptbin[1], etabin[0], etabin[1]))
                    print "%-10s  %6.0f  %.3f  %.3f" % (hname2, h2.Integral(0,1000), h2.GetMean(), h2.GetRMS())
                    projappend_ = True
                    intree.Project(projname(hname2), "smear_gen_to_reco(hJet_pt[1], genBbar.pt, genBbar.eta) / genBbar.pt", "hJet_id[1]==1 && hJet_puJetIdL[1]>0 && hJet_csv[1]>0 && hJet_csv_nominal[1]>0.244 && hJet_flavour[1]==-5 && deltaR(hJet_genEta[1],hJet_genPhi[1],genBbar.eta,genBbar.phi)<0.1 && (%f<genH.pt && genH.pt<%f) && (%f<genBbar.pt && genBbar.pt<%f) && (%f<abs(genBbar.eta) && abs(genBbar.eta)<%f)" % (genzptbin[0], genzptbin[1], genptbin[0], genptbin[1], etabin[0], etabin[1]))
                    print "%-10s  %6.0f  %.3f  %.3f" % (hname2, h2.Integral(0,1000), h2.GetMean(), h2.GetRMS())
                
                
                #h1.SetMaximum(h1.GetMaximum()*2e1)  # for log scale
                #h1.Draw("hist")
                #h2.Draw("samehist")
                #gPad.SetLogy()
                
                #latex.DrawLatex(0.19, 0.90, "%.0f<p_{T}^{Z}<%.0f, %.0f<p_{T}^{j}<%.0f, %.1f<|#eta^{j}|<%.1f" % (genzptbin[0], genzptbin[1], genptbin[0], genptbin[1], etabin[0], etabin[1]))
                #latex.DrawLatex(0.19, 0.85, "#color[4]{FullSim #mu=%.3f, #sigma=%.3f}" % (h1.GetMean(), h1.GetRMS()))
                #latex.DrawLatex(0.19, 0.80, "#color[2]{NGS param. #mu=%.3f, #sigma=%.3f}" % (h2.GetMean(), h2.GetRMS()))
                #gPad.Print("plots/zmmtozbb_%s.png" % h1.GetName())
                
                histos1.append(h1)
                histos2.append(h2)
                hl.SetBinContent(1, genzptbin[0])
                hl.SetBinContent(2, genzptbin[1])
                hl.SetBinContent(3, genptbin[0])
                hl.SetBinContent(4, genptbin[1])
                hl.SetBinContent(5, etabin[0])
                hl.SetBinContent(6, etabin[1])
                hlabels.append(hl)
                

hhistos1 = []
hhistos2 = []
nhistos = len(histos1)/nfiles
for i in xrange(nhistos):
    h1 = histos1[i]
    hh1 = h1.Clone(h1.GetName().replace("h1","hh1"))
    
    h2 = histos2[i]
    hh2 = h2.Clone(h2.GetName().replace("h2","hh2"))
    
    for j in xrange(nfiles-1):
        hh1.Add(histos1[(j+1)*nhistos+i])
        hh2.Add(histos2[(j+1)*nhistos+i])
    
    hh1.SetMaximum(hh1.GetMaximum()*3e1)  # for log scale
    hh1.Draw("")
    if hh2.Integral() > 1:
        hh2.Draw("samehist")
    gPad.SetLogy()
    
    hl = hlabels[i]
    latex.DrawLatex(0.19, 0.90, "%.0f<p_{T}^{Z}<%.0f, %.0f<p_{T}^{j}<%.0f, %.1f<|#eta^{j}|<%.1f" % tuple([hl.GetBinContent(l) for l in xrange(1,7)]))
    latex.DrawLatex(0.19, 0.85, "#color[1]{FullSim}")
    latex.DrawLatex(0.19, 0.80, "#color[4]{NGS gaussian}")
    latex.DrawLatex(0.42, 0.85, "#color[1]{#mu=%.3f, #sigma=%.3f}" % (hh1.GetMean(), hh1.GetRMS()))
    latex.DrawLatex(0.42, 0.80, "#color[4]{#mu=%.3f, #sigma=%.3f}" % (hh2.GetMean(), hh2.GetRMS()))
    gPad.Print("plots/zmmtozbb_%s.png" % hh1.GetName())
    gPad.Print("plots/zmmtozbb_%s.pdf" % hh1.GetName())
    
    hhistos1.append(hh1)
    hhistos2.append(hh2)
    

outfile = TFile.Open(outfilename, "RECREATE")
for h in histos1:
    h.Write()
for h in hhistos1:
    h.Write()
for h in hlabels:
    h.Write()
outfile.Close()
