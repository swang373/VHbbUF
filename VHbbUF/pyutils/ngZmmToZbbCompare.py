from ROOT import TFile, TH1F, TH1, TLatex, gPad, gROOT, gStyle
import copy

file1 = "Step3_20130314/Step3_ZbbHinv125.root"
file2 = "skim/Step3_ZbbHinvZmmToZbb.root"
wait = False
gROOT.SetBatch(1)

#varexp = "deltaR(genB.eta, genB.phi, genBbar.eta, genBbar.phi)"
varexp = "hJet_pt[0]"
selection = "max(hJet_pt[0],hJet_pt[1])>60 && min(hJet_pt[0],hJet_pt[1])>30 && hJet_csv_nominal[0]>0.244 && hJet_csv_nominal[1]>0.244 && H.pt>80 && METtype1corr.et>80 && H.dR>0.5"
#varexp = "hJet_pt[0]/hJet_genPt[0]"
#selection = "hJet_id[0]==1 && hJet_puJetIdL[0]>0 && hJet_csv[0]>0 && hJet_csv_nominal[0]>0.244 && deltaR(hJet_eta[0],hJet_phi[0],hJet_genEta[0],hJet_genPhi[0])<0.3 && (%f<genH.pt && genH.pt<%f) && (%f<hJet_genPt[0] && hJet_genPt[0]<%f) && (%f<abs(hJet_genEta[0]) && abs(hJet_genEta[0])<%f)" % (0, 1000, 70, 100, 0, 0.5)
##varexp = "hJet_e[0]"
#varexp = "hJet_pt[0]/genB.pt"
#selection = "hJet_id[0]==1 && hJet_puJetIdL[0]>0 && hJet_csv[0]>0 && hJet_csv_nominal[0]>0. && hJet_flavour[0]==5 && deltaR(hJet_genEta[0],hJet_genPhi[0],genB.eta,genB.phi)<0.1 && (%f<genH.pt && genH.pt<%f) && (%f<genB.pt && genB.pt<%f) && (%f<abs(genB.eta) && abs(genB.eta)<%f)" % (0, 1000, 60, 80, 0.0, 2.5)
##varexp_ = "hJet_e[0]"
#varexp_ = "hJet_pt[0]/genBbar.pt"
#selection_ = "hJet_id[0]==1 && hJet_puJetIdL[0]>0 && hJet_csv[0]>0 && hJet_csv_nominal[0]>0. && hJet_flavour[0]==-5 && deltaR(hJet_genEta[0],hJet_genPhi[0],genBbar.eta,genBbar.phi)<0.1 && (%f<genH.pt && genH.pt<%f) && (%f<genBbar.pt && genBbar.pt<%f) && (%f<abs(genBbar.eta) && abs(genBbar.eta)<%f)" % (0, 60, 80, 10000, 0.0, 2.5)
#varexph3 = "smear_gen_by_double_gaussian(genB.pt, 0.995 ,  0.080 ,  0.876 ,  2.290 ,  0.389)/genB.pt"
#varexph3_ = "smear_gen_by_double_gaussian(genBbar.pt, 0.995 ,  0.080 ,  0.876 ,  2.290 ,  0.389)/genBbar.pt"

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
t1 = tfile1.tree
t2 = tfile2.tree

histos = []
if True:
    hname = "ptj1"
    varexp = "max(hJet_pt[0], hJet_pt[1])"
    etc = ("; p_{T}(j_{1}) [GeV]; (normalized)", 100, 0, 300)
    histos.append((hname, varexp, etc))

    hname = "ptj2"
    varexp = "min(hJet_pt[0], hJet_pt[1])"
    etc = ("; p_{T}(j_{2}) [GeV]; (normalized)", 100, 0, 300)
    histos.append((hname, varexp, etc))

    hname = "ptjj"
    varexp = "H.pt"
    etc = ("; p_{T}(jj) [GeV]; (normalized)", 100, 0, 300)
    histos.append((hname, varexp, etc))

    hname = "mjj"
    varexp = "H.mass"
    etc = ("; m(jj) [GeV]; (normalized)", 100, 0, 300)
    histos.append((hname, varexp, etc))

    hname = "drjj"
    varexp = "H.dR"
    etc = ("; #DeltaR(jj) [GeV]; (normalized)", 100, 0, 4)
    histos.append((hname, varexp, etc))

    hname = "detajj"
    varexp = "H.dEta"
    etc = ("; #Delta#eta(jj); (normalized)", 100, 0, 4)
    histos.append((hname, varexp, etc))

    hname = "dphijj"
    varexp = "H.dPhi"
    etc = ("; #Delta#phi(jj); (normalized)", 100, 0, 4)
    histos.append((hname, varexp, etc))

    hname = "met"
    varexp = "METtype1corr.et"
    etc = ("; MET [GeV]; (normalized)", 100, 0, 300)
    histos.append((hname, varexp, etc))

    hname = "jjmetdphi"
    varexp = "HMETdPhi"
    etc = ("; #Delta#phi(jj,MET); (normalized)", 100, 0, 4)
    histos.append((hname, varexp, etc))

    hname = "csvmax"
    varexp = "max(hJet_csv_nominal[0],hJet_csv_nominal[1])"
    etc = ("; CSV_{max} [GeV]; (normalized)", 110, 0, 1.1)
    histos.append((hname, varexp, etc))

    hname = "csvmin"
    varexp = "min(hJet_csv_nominal[0],hJet_csv_nominal[1])"
    etc = ("; CSV_{min} [GeV]; (normalized)", 110, 0, 1.1)
    histos.append((hname, varexp, etc))
    
    hname = "csvj1"
    varexp = "hJet_csv_nominal[0]"
    etc = ("; CSV(j_{1}) [GeV]; (normalized)", 110, 0, 1.1)
    histos.append((hname, varexp, etc))

    hname = "csvj2"
    varexp = "hJet_csv_nominal[0]"
    etc = ("; CSV(j_{2}) [GeV]; (normalized)", 110, 0, 1.1)
    histos.append((hname, varexp, etc))

    hname = "ptj1reg"
    varexp = "max(hJet_ptReg[0], hJet_ptReg[1])"
    etc = ("; regressed p_{T}(j_{1}) [GeV]; (normalized)", 100, 0, 300)
    histos.append((hname, varexp, etc))

    hname = "ptj2reg"
    varexp = "min(hJet_ptReg[0], hJet_ptReg[1])"
    etc = ("; regressed p_{T}(j_{2}) [GeV]; (normalized)", 100, 0, 300)
    histos.append((hname, varexp, etc))
    
    hname = "ptjjreg"
    varexp = "HptReg"
    etc = ("; regressed p_{T}(jj) [GeV]; (normalized)", 100, 0, 300)
    histos.append((hname, varexp, etc))

    hname = "mjjreg"
    varexp = "HmassReg"
    etc = ("; regressed m(jj) [GeV]; (normalized)", 100, 0, 300)
    histos.append((hname, varexp, etc))

    hname = "deltaPullAngle"
    varexp = "deltaPullAngle"
    etc = ("; #Delta#theta_{pull}; (normalized)", 100, -4, 4)
    histos.append((hname, varexp, etc))

#_______________________________________________________________________________
for ih, (hname, varexp, etc) in enumerate(histos):
    h1 = TH1F("h1_%i" % ih, etc[0], etc[1], etc[2], etc[3])
    h2 = TH1F("h2_%i" % ih, etc[0], etc[1], etc[2], etc[3])
    h3 = TH1F("h3_%i" % ih, etc[0], etc[1], etc[2], etc[3])
    h1.SetLineWidth(2)
    h1.SetLineColor(1)
    h2.SetLineWidth(2)
    h2.SetLineColor(4)
    h3.SetLineWidth(2)
    h3.SetLineColor(6)

    t1.Project("h1_%i" % ih, varexp, selection)
    t2.Project("h2_%i" % ih, varexp, selection)
    #t1.Project("h3", varexph3, selection)
    
    #t1.Project("+h1", varexp_, selection_)
    #t2.Project("+h2", varexp_, selection_)
    ##t1.Project("+h3", varexph3_, selection_)

    h1.Scale(1.0/h1.Integral())
    h2.Scale(1.0/h2.Integral())
    h1.SetMaximum(h1.GetMaximum()*1.3)
    h1.Draw("hist")
    h2.Draw("samehist")
    #h3.Draw("samehist")
    
    latex.DrawLatex(0.19, 0.90, "MET>80, p_{T}(jj)>80, p_{T}(j)>60,30, CSV(j)=CSVL")
    latex.DrawLatex(0.19, 0.85, "#color[1]{FullSim}")
    latex.DrawLatex(0.19, 0.80, "#color[4]{Smear/Histo}")
    gPad.Print("plots_zmmtozbb/zmmtozbb_compare_ZbbHinv_%s.png" % hname)
    gPad.Print("plots_zmmtozbb/zmmtozbb_compare_ZbbHinv_%s.pdf" % hname)
    
    if "mjj" in hname:
        print h1.GetMean(), h1.GetRMS(), h2.GetMean(), h2.GetRMS()
    
    if wait:  raw_input("Press Enter to continue...")
