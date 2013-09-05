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

infilenames_step3 = [
"Step3_20130314/Step3_ZbbHinv105.root",
"Step3_20130314/Step3_ZbbHinv115.root",
"Step3_20130314/Step3_ZbbHinv125.root",
"Step3_20130314/Step3_ZbbHinv135.root",
"Step3_20130314/Step3_ZbbHinv145.root",
"Step3_20130314/Step3_ZbbHinv150.root",
]

outfilename = "zmmtozbb_jetres.root"
outfilename_csv = "zmmtozbb_csv.root"

doCSV = False
useGenJet = False
useParton = not useGenJet
useStep3 = False  # should only be used for regressed pt
batch = True

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
ptbins = genptbins
etabins = [
    (0.0,0.5),
    (0.5,1.0),
    (1.0,1.5),
    (1.5,2.5),
#    (2.5,5.0),
]


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

if useStep3:
    infilenames = infilenames_step3
if doCSV:
    outfilename = outfilename_csv
if batch:
    gROOT.SetBatch(1)

projappend_ = False
def projname(name):
    if projappend_:
        return "+"+name
    return name

#_______________________________________________________________________________
histos1 = []  # pT resolution
histos2 = []  # pT resolution (smeared)
histos3 = []  # pT resolution (regressed)
hlabels = []
pulls = []

nfiles = 0
print "            norm   mean   RMS  "
for infilename in infilenames:
    infile = TFile.Open(infilename)
    intree = infile.tree
    nfiles += 1
    
    if not doCSV:
        for etabin in etabins:
            for genptbin in genptbins:
                for genzptbin in genzptbins:
                    
                    hname1 = "h1_%i" % len(histos1)
                    hname2 = "h2_%i" % len(histos2)
                    hname3 = "h3_%i" % len(histos3)
                    hlabel = "hl_%i" % len(hlabels)
                    #h1 = TH1F(hname1, "; p_{T}^{reco} / p_{T}^{gen}; Events / (0.02)", 150, 0, 3.0)
                    #h2 = TH1F(hname2, "; p_{T}^{reco} / p_{T}^{gen}; Events / (0.02)", 150, 0, 3.0)
                    h1 = TH1F(hname1, "; p_{T}^{reco} / p_{T}^{parton}; Events / (0.02)", 150, 0, 3.0)
                    h2 = TH1F(hname2, "; p_{T}^{reco} / p_{T}^{parton}; Events / (0.02)", 150, 0, 3.0)
                    h3 = TH1F(hname3, "; p_{T}^{regressed} / p_{T}^{parton}; Events / (0.02)", 150, 0, 3.0)
                    hl = TH1F(hlabel, "; label", 10, 0, 10)
                    h1.SetLineColor(1)
                    h2.SetLineColor(4)
                    h2.SetLineWidth(2)
                    h3.SetLineColor(418)  # kGreen+2
                    h3.SetMarkerStyle(24)
                    h3.SetMarkerColor(418)  # kGreen+2
                    
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
                        
                        # regressed jet 1
                        projappend_ = False
                        intree.Project(projname(hname3), "hJet_ptReg[0] / genB.pt", "hJet_id[0]==1 && hJet_puJetIdL[0]>0 && hJet_csv[0]>0 && hJet_csv_nominal[0]>0.244 && hJet_flavour[0]==5 && deltaR(hJet_genEta[0],hJet_genPhi[0],genB.eta,genB.phi)<0.1 && (%f<genH.pt && genH.pt<%f) && (%f<genB.pt && genB.pt<%f) && (%f<abs(genB.eta) && abs(genB.eta)<%f)" % (genzptbin[0], genzptbin[1], genptbin[0], genptbin[1], etabin[0], etabin[1]))
                        print "%-10s  %6.0f  %.3f  %.3f" % (hname3, h3.Integral(0,1000), h3.GetMean(), h3.GetRMS())
                        projappend_ = True
                        intree.Project(projname(hname3), "hJet_ptReg[0] / genBbar.pt", "hJet_id[0]==1 && hJet_puJetIdL[0]>0 && hJet_csv[0]>0 && hJet_csv_nominal[0]>0.244 && hJet_flavour[0]==-5 && deltaR(hJet_genEta[0],hJet_genPhi[0],genBbar.eta,genBbar.phi)<0.1 && (%f<genH.pt && genH.pt<%f) && (%f<genBbar.pt && genBbar.pt<%f) && (%f<abs(genBbar.eta) && abs(genBbar.eta)<%f)" % (genzptbin[0], genzptbin[1], genptbin[0], genptbin[1], etabin[0], etabin[1]))
                        print "%-10s  %6.0f  %.3f  %.3f" % (hname3, h3.Integral(0,1000), h3.GetMean(), h3.GetRMS())
                        # regressed jet 2
                        projappend_ = True
                        intree.Project(projname(hname3), "hJet_ptReg[1] / genB.pt", "hJet_id[1]==1 && hJet_puJetIdL[1]>0 && hJet_csv[1]>0 && hJet_csv_nominal[1]>0.244 && hJet_flavour[1]==5 && deltaR(hJet_genEta[1],hJet_genPhi[1],genB.eta,genB.phi)<0.1 && (%f<genH.pt && genH.pt<%f) && (%f<genB.pt && genB.pt<%f) && (%f<abs(genB.eta) && abs(genB.eta)<%f)" % (genzptbin[0], genzptbin[1], genptbin[0], genptbin[1], etabin[0], etabin[1]))
                        print "%-10s  %6.0f  %.3f  %.3f" % (hname3, h3.Integral(0,1000), h3.GetMean(), h3.GetRMS())
                        projappend_ = True
                        intree.Project(projname(hname3), "hJet_ptReg[1] / genBbar.pt", "hJet_id[1]==1 && hJet_puJetIdL[1]>0 && hJet_csv[1]>0 && hJet_csv_nominal[1]>0.244 && hJet_flavour[1]==-5 && deltaR(hJet_genEta[1],hJet_genPhi[1],genBbar.eta,genBbar.phi)<0.1 && (%f<genH.pt && genH.pt<%f) && (%f<genBbar.pt && genBbar.pt<%f) && (%f<abs(genBbar.eta) && abs(genBbar.eta)<%f)" % (genzptbin[0], genzptbin[1], genptbin[0], genptbin[1], etabin[0], etabin[1]))
                        print "%-10s  %6.0f  %.3f  %.3f" % (hname3, h3.Integral(0,1000), h3.GetMean(), h3.GetRMS())
                    
                    #h1.SetMaximum(h1.GetMaximum()*2e1)  # for log scale
                    #h1.Draw("hist")
                    #h2.Draw("samehist")
                    #gPad.SetLogy()
                    
                    #latex.DrawLatex(0.19, 0.90, "%.0f<p_{T}^{Z}<%.0f, %.0f<p_{T}^{j}<%.0f, %.1f<|#eta^{j}|<%.1f" % (genzptbin[0], genzptbin[1], genptbin[0], genptbin[1], etabin[0], etabin[1]))
                    #latex.DrawLatex(0.19, 0.85, "#color[4]{FullSim #mu=%.3f, #sigma=%.3f}" % (h1.GetMean(), h1.GetRMS()))
                    #latex.DrawLatex(0.19, 0.80, "#color[2]{NSC param. #mu=%.3f, #sigma=%.3f}" % (h2.GetMean(), h2.GetRMS()))
                    #gPad.Print("plots/zmmtozbb_%s.png" % h1.GetName())
                    
                    histos1.append(h1)
                    histos2.append(h2)
                    histos3.append(h3)
                    hl.SetBinContent(1, genzptbin[0])
                    hl.SetBinContent(2, genzptbin[1])
                    hl.SetBinContent(3, genptbin[0])
                    hl.SetBinContent(4, genptbin[1])
                    hl.SetBinContent(5, etabin[0])
                    hl.SetBinContent(6, etabin[1])
                    hlabels.append(hl)
    
    else:  # if doCSV
        for etabin in etabins:
            #for genptbin in genptbins:
            for ptbin in ptbins:
                for genzptbin in genzptbins:
                    
                    hname1 = "h1_%i" % len(histos1)
                    hlabel = "hl_%i" % len(hlabels)
                    h1 = TH1F(hname1, "; CSV; Events / (0.01)", 100, 0, 1.0)
                    hl = TH1F(hlabel, "; label", 10, 0, 10)
                    h1.SetLineColor(1)
                    
                    # N.B. genZ.pt is zero in Pythia sample
                    useRecoJet = True
                    if useRecoJet:
                        projappend_ = False
                        intree.Project(projname(hname1), "min(0.999,hJet_csv_nominal[0])", "hJet_id[0]==1 && hJet_puJetIdL[0]>0 && hJet_csv[0]>0 && abs(hJet_flavour[0])==5 && (%f<genH.pt && genH.pt<%f) && (%f<hJet_pt[0] && hJet_pt[0]<%f) && (%f<abs(hJet_eta[0]) && abs(hJet_eta[0])<%f)" % (genzptbin[0], genzptbin[1], ptbin[0], ptbin[1], etabin[0], etabin[1]))
                        print "%-10s  %6.0f  %.3f  %.3f" % (hname1, h1.Integral(0,1000), h1.GetMean(), h1.GetRMS())
                        projappend_ = True
                        intree.Project(projname(hname1), "min(0.999,hJet_csv_nominal[1])", "hJet_id[1]==1 && hJet_puJetIdL[1]>0 && hJet_csv[1]>0 && abs(hJet_flavour[1])==5 && (%f<genH.pt && genH.pt<%f) && (%f<hJet_pt[1] && hJet_pt[1]<%f) && (%f<abs(hJet_eta[1]) && abs(hJet_eta[1])<%f)" % (genzptbin[0], genzptbin[1], ptbin[0], ptbin[1], etabin[0], etabin[1]))
                        print "%-10s  %6.0f  %.3f  %.3f" % (hname1, h1.Integral(0,1000), h1.GetMean(), h1.GetRMS())
                        
                    histos1.append(h1)
                    hl.SetBinContent(1, genzptbin[0])
                    hl.SetBinContent(2, genzptbin[1])
                    hl.SetBinContent(3, ptbin[0])
                    hl.SetBinContent(4, ptbin[1])
                    hl.SetBinContent(5, etabin[0])
                    hl.SetBinContent(6, etabin[1])
                    hlabels.append(hl)
        
        h1 = TH1F("pull_%i" % len(pulls), "; deltaPullAngle/#pi; Events / (0.02)", 100, -1.0, 1.0)
        intree.Project("pull_%i" % len(pulls), "deltaPullAngle/pi", "hJet_id[0]==1 && hJet_puJetIdL[0]>0 && hJet_csv[0]>0 && abs(hJet_flavour[0])==5 && hJet_pt[0]>30 && abs(hJet_eta[0])<2.5 && hJet_id[1]==1 && hJet_puJetIdL[1]>0 && hJet_csv[1]>0 && abs(hJet_flavour[1])==5 && hJet_pt[1]>30 && abs(hJet_eta[1])<2.5")
        pulls.append(h1)

hhistos1 = []
hhistos2 = []
hhistos3 = []
nhistos = len(histos1)/nfiles
numbers = []

if not doCSV:
    def truncate(h, percent=4.550):
        begin, end = 0, h.GetNbinsX()+2
        low, up = begin, end
        integral = h.Integral(begin, end)
        for ib in xrange(begin, end):
            if (h.Integral(begin, low) < percent/2./100. * integral):
                low += 1
            if (h.Integral(up, end) < percent/2./100. * integral):
                up -= 1
        for ib in xrange(begin, low):
            h.SetBinContent(ib, 0)
            h.SetBinError(ib, 0)
        for ib in xrange(up, end):
            h.SetBinContent(ib, 0)
            h.SetBinError(ib, 0)
        print h.Integral(begin, end) / integral
    
    for i in xrange(nhistos):
        h1 = histos1[i]
        hh1 = h1.Clone(h1.GetName().replace("h1","hh1"))
        h2 = histos2[i]
        hh2 = h2.Clone(h2.GetName().replace("h2","hh2"))
        h3 = histos3[i]
        hh3 = h3.Clone(h3.GetName().replace("h3","hh3"))
        
        for j in xrange(nfiles-1):
            hh1.Add(histos1[(j+1)*nhistos+i])
            hh2.Add(histos2[(j+1)*nhistos+i])
            hh3.Add(histos3[(j+1)*nhistos+i])
        
        hh1.SetMaximum(hh1.GetMaximum()*3e1)  # for log scale
        hh1.Draw("")
        gPad.SetLogy()
        if hh2.Integral() > 1:
            hh2.Draw("samehist")
        
        hl = hlabels[i]
        latex.DrawLatex(0.19, 0.90, "%.0f<p_{T}^{Z}<%.0f, %.0f<p_{T}^{j}<%.0f, %.1f<|#eta^{j}|<%.1f" % tuple([hl.GetBinContent(l) for l in xrange(1,7)]))
        latex.DrawLatex(0.19, 0.85, "#color[1]{FullSim}")
        latex.DrawLatex(0.19, 0.80, "#color[4]{NSC gaussian}")
        latex.DrawLatex(0.42, 0.85, "#color[1]{#mu=%.3f, #sigma=%.3f}" % (hh1.GetMean(), hh1.GetRMS()))
        latex.DrawLatex(0.42, 0.80, "#color[4]{#mu=%.3f, #sigma=%.3f}" % (hh2.GetMean(), hh2.GetRMS()))
        gPad.Print("plots/zmmtozbb_%s.png" % hh1.GetName())
        gPad.Print("plots/zmmtozbb_%s.pdf" % hh1.GetName())
        
        hh1.Draw("")
        if hh3.Integral() > 1:
            hh3.Draw("same")
        
        latex.DrawLatex(0.19, 0.90, "%.0f<p_{T}^{Z}<%.0f, %.0f<p_{T}^{j}<%.0f, %.1f<|#eta^{j}|<%.1f" % tuple([hl.GetBinContent(l) for l in xrange(1,7)]))
        latex.DrawLatex(0.19, 0.85, "#color[1]{FullSim}")
        latex.DrawLatex(0.19, 0.80, "#color[418]{Regressed}")
        latex.DrawLatex(0.42, 0.85, "#color[1]{#mu=%.3f, #sigma=%.3f}" % (hh1.GetMean(), hh1.GetRMS()))
        latex.DrawLatex(0.42, 0.80, "#color[418]{#mu=%.3f, #sigma=%.3f (%.0f%%)}" % (hh3.GetMean(), hh3.GetRMS(), (hh3.GetRMS()/hh1.GetRMS()-1)*100))
        gPad.Print("plots/zmmtozbb_reg_%s.png" % hh1.GetName())
        gPad.Print("plots/zmmtozbb_reg_%s.pdf" % hh1.GetName())
        
        hh1_T = hh1.Clone(h1.GetName().replace("hh1","hh1_T"))
        hh3_T = hh3.Clone(h3.GetName().replace("hh3","hh3_T"))
        truncate(hh1_T)
        truncate(hh3_T)
        
        hh1.Draw("")
        if hh3.Integral() > 1:
            hh3.Draw("same")
        
        latex.DrawLatex(0.19, 0.90, "%.0f<p_{T}^{Z}<%.0f, %.0f<p_{T}^{j}<%.0f, %.1f<|#eta^{j}|<%.1f" % tuple([hl.GetBinContent(l) for l in xrange(1,7)]))
        latex.DrawLatex(0.19, 0.85, "#color[1]{FullSim}")
        latex.DrawLatex(0.19, 0.80, "#color[418]{Regressed}")
        latex.DrawLatex(0.42, 0.85, "#color[1]{#mu_{95%%}=%.3f, #sigma_{95%%}=%.3f}" % (hh1_T.GetMean(), hh1_T.GetRMS()))
        latex.DrawLatex(0.42, 0.80, "#color[418]{#mu_{95%%}=%.3f, #sigma_{95%%}=%.3f (%.0f%%)}" % (hh3_T.GetMean(), hh3_T.GetRMS(), (hh3_T.GetRMS()/hh1_T.GetRMS()-1)*100))
        gPad.Print("plots/zmmtozbb_reg95_%s.png" % hh1.GetName())
        gPad.Print("plots/zmmtozbb_reg95_%s.pdf" % hh1.GetName())
        
        hhistos1.append(hh1)
        hhistos2.append(hh2)
        hhistos3.append(hh3)
        numbers.append(hh3_T.GetRMS()/hh1_T.GetRMS()-1)

else:  # if doCSV
    for i in xrange(nhistos):
        h1 = histos1[i]
        hh1 = h1.Clone(h1.GetName().replace("h1","hh1"))
        
        for j in xrange(nfiles-1):
            hh1.Add(histos1[(j+1)*nhistos+i])
        
        hh1.SetMaximum(hh1.GetMaximum()*1.3)  # for linear scale
        hh1.Draw("")
        #gPad.SetLogy()
        
        hl = hlabels[i]
        latex.DrawLatex(0.19, 0.90, "%.0f<p_{T}^{Z}<%.0f, %.0f<p_{T}^{j}<%.0f, %.1f<|#eta^{j}|<%.1f" % tuple([hl.GetBinContent(l) for l in xrange(1,7)]))
        latex.DrawLatex(0.19, 0.85, "#color[1]{FullSim(reshaped)}")
        latex.DrawLatex(0.50, 0.85, "#color[1]{#mu=%.3f, #sigma=%.3f}" % (hh1.GetMean(), hh1.GetRMS()))
        gPad.Print("plots/zmmtozbb_csv_%s.png" % hh1.GetName())
        gPad.Print("plots/zmmtozbb_csv_%s.pdf" % hh1.GetName())
    
        hhistos1.append(hh1)
    
    # deltapullangle
    hh1 = pulls[0].Clone("deltaPullAngleOverPi")
    for h1 in pulls[1:]:
        hh1.Add(h1)
    hhistos1.append(hh1)

for n in numbers:
    print "%.3f," % n

outfile = TFile.Open(outfilename, "RECREATE")
for h in histos1:
    h.Write()
for h in hhistos1:
    h.Write()
for h in hlabels:
    h.Write()
outfile.Close()
