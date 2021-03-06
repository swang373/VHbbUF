from ROOT import TFile, TH1F, TH1, TLatex, gPad, gROOT, gStyle
import copy

#file1 = "Step3_20130314/Step3_ZbbHinv125.root"
#file2 = "Step3_20130314/Step3_ZbbHinvZmmToZbb.root"
file1 = "skim_ZnnH_baseline/skim_DYJetsZmmToZbb.root"
file2 = "Step3_20130314/Step3_DYJetsZmmToZbb.root"
use_lep_for_t1 = True
wait = False
gROOT.SetBatch(1)

#varexp = "deltaR(genB.eta, genB.phi, genBbar.eta, genBbar.phi)"
varexp = "hJet_pt[0]"
selection = "max(hJet_pt[0],hJet_pt[1])>60 && min(hJet_pt[0],hJet_pt[1])>30 && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5 && hJet_csv_nominal[0]>0.244 && hJet_csv_nominal[1]>0.244 && H.pt>80 && METtype1corr.et>80 && H.dR>0.5 && !(vLepton_pt[0]<5 || abs(vLepton_eta[0])>5 || vLepton_pt[1]<5 || abs(vLepton_eta[1])>5)"
selectionlep = "max(vLepton_pt[0],vLepton_pt[1])>60 && min(vLepton_pt[0],vLepton_pt[1])>30 && abs(vLepton_eta[0])<2.5 && abs(vLepton_eta[1])<2.5 && V.pt>80 && METtype1corr.et>80 && V.dR>0.5 && !(vLepton_pt[0]<5 || abs(vLepton_eta[0])>5 || vLepton_pt[1]<5 || abs(vLepton_eta[1])>5)"
selection = "max(hJet_pt[0],hJet_pt[1])>60 && min(hJet_pt[0],hJet_pt[1])>30 && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5 && hJet_csv_nominal[0]>0.244 && hJet_csv_nominal[1]>0.244 && H.pt>80 && METtype1corr.et>80 && H.dR>0.5 && !(vLepton_pt[0]<5 || abs(vLepton_eta[0])>5 || vLepton_pt[1]<5 || abs(vLepton_eta[1])>5) && min(Min$(deltaPhiMETjets(METtype1corr.phi, hJet_phi, hJet_pt, hJet_eta, 25, 4.5)), Min$(deltaPhiMETjets(METtype1corr.phi, aJet_phi, aJet_pt, aJet_eta, 25, 4.5)))>0.5 && H.mass>0 && abs(deltaPhi(METtype1corr.phi, H.phi))>0 && H.mass>0"
selectionlep = "max(vLepton_pt[0],vLepton_pt[1])>60 && min(vLepton_pt[0],vLepton_pt[1])>30 && abs(vLepton_eta[0])<2.5 && abs(vLepton_eta[1])<2.5 && V.pt>80 && METtype1corr.et>80 && V.dR>0.5 && !(vLepton_pt[0]<5 || abs(vLepton_eta[0])>5 || vLepton_pt[1]<5 || abs(vLepton_eta[1])>5) && min(Min$(deltaPhiMETjets(METtype1corr.phi, vLepton_phi, vLepton_pt, vLepton_eta, 25, 4.5)), Min$(deltaPhiMETjets(METtype1corr.phi, aJet_phi, aJet_pt, aJet_eta, 25, 4.5)))>0.5 && V.mass>0 && abs(deltaPhi(METtype1corr.phi, V.phi))>0 && V.mass>0"

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
    etc = ("; p_{T}(j_{1}) [GeV]; (normalized)", 100, 0, 250)
    histos.append((hname, varexp, etc))

    hname = "ptj2"
    varexp = "min(hJet_pt[0], hJet_pt[1])"
    etc = ("; p_{T}(j_{2}) [GeV]; (normalized)", 100, 0, 250)
    histos.append((hname, varexp, etc))

    hname = "ptjj"
    varexp = "H.pt"
    etc = ("; p_{T}(jj) [GeV]; (normalized)", 100, 0, 250)
    histos.append((hname, varexp, etc))

    hname = "mjj"
    varexp = "H.mass"
    etc = ("; m(jj) [GeV]; (normalized)", 100, 0, 250)
    histos.append((hname, varexp, etc))

    hname = "drjj"
    varexp = "deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])"
    etc = ("; #DeltaR(jj) [GeV]; (normalized)", 100, 0, 4)
    histos.append((hname, varexp, etc))

    hname = "detajj"
    varexp = "abs(hJet_eta[0]-hJet_eta[1])"
    etc = ("; #Delta#eta(jj); (normalized)", 100, 0, 4)
    histos.append((hname, varexp, etc))

    hname = "dphijj"
    varexp = "deltaPhi(hJet_phi[0], hJet_phi[1])"
    etc = ("; #Delta#phi(jj); (normalized)", 100, 0, 4)
    histos.append((hname, varexp, etc))

    hname = "met"
    varexp = "METtype1corr.et"
    etc = ("; E_{T}^{miss} [GeV]; (normalized)", 100, 0, 250)
    histos.append((hname, varexp, etc))

    hname = "jjmetdphi"
    varexp = "deltaPhi(H.phi,METtype1corr.phi)"
    etc = ("; #Delta#phi(jj,E_{T}^{miss}); (normalized)", 100, 0, 4)
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
    varexp = "hJet_csv_nominal[1]"
    etc = ("; CSV(j_{2}) [GeV]; (normalized)", 110, 0, 1.1)
    histos.append((hname, varexp, etc))

#    hname = "deltaPullAngle"
#    varexp = "deltaPullAngle"
#    etc = ("; #Delta#theta_{pull}; (normalized)", 100, -4, 4)
#    histos.append((hname, varexp, etc))
    
    hname = "mindphi"
    #varexp = "mindPhiMETJet_dPhi"
    varexp = "min(Min$(deltaPhiMETjets(METtype1corr.phi, hJet_phi, hJet_pt, hJet_eta, 25, 4.5)), Min$(deltaPhiMETjets(METtype1corr.phi, aJet_phi, aJet_pt, aJet_eta, 25, 4.5)))"
    etc = ("; min #Delta #phi(E_{T}^{miss},jet); (normalized)", 100, 0, 4)
    histos.append((hname, varexp, etc))

#    hname = "ptj1reg"
#    varexp = "max(hJet_ptReg[0], hJet_ptReg[1])"
#    etc = ("; regressed p_{T}(j_{1}) [GeV]; (normalized)", 100, 0, 250)
#    histos.append((hname, varexp, etc))
#
#    hname = "ptj2reg"
#    varexp = "min(hJet_ptReg[0], hJet_ptReg[1])"
#    etc = ("; regressed p_{T}(j_{2}) [GeV]; (normalized)", 100, 0, 250)
#    histos.append((hname, varexp, etc))
#    
#    hname = "ptjjreg"
#    varexp = "HptReg"
#    etc = ("; regressed p_{T}(jj) [GeV]; (normalized)", 100, 0, 250)
#    histos.append((hname, varexp, etc))
#
#    hname = "mjjreg"
#    varexp = "HmassReg"
#    etc = ("; regressed m(jj) [GeV]; (normalized)", 100, 0, 250)
#    histos.append((hname, varexp, etc))


#_______________________________________________________________________________
for ih, (hname, varexp, etc) in enumerate(histos):
    h1 = TH1F("h1_%i" % ih, etc[0], etc[1]/4, etc[2], etc[3])
    h2 = TH1F("h2_%i" % ih, etc[0], etc[1]/4, etc[2], etc[3])
    h3 = TH1F("h3_%i" % ih, etc[0], etc[1]/4, etc[2], etc[3])
    h1.SetLineWidth(2)
    h1.SetLineColor(1)
    h1.SetFillColor(920)
    h1.SetFillStyle(3003)
    h2.SetLineWidth(2)
    h2.SetLineColor(4)
    h2.SetFillColor(591)  # kBlue-9
    h2.SetFillStyle(3005)
    h3.SetLineWidth(2)
    h3.SetLineColor(6)

    if use_lep_for_t1:
        t1.Project("h1_%i" % ih, varexp.replace("hJet_","vLepton_").replace("H.","V."), selectionlep)
    else:
        t1.Project("h1_%i" % ih, varexp, selection)
    t2.Project("h2_%i" % ih, varexp, selection)
    #t1.Project("h3", varexph3, selection)
    
    #t1.Project("+h1", varexp_, selection_)
    #t2.Project("+h2", varexp_, selection_)
    ##t1.Project("+h3", varexph3_, selection_)

    if h1.Integral()>0:  h1.Scale(1.0/h1.Integral())
    if h2.Integral()>0:  h2.Scale(1.0/h2.Integral())
    if h1.Integral()>0:
        h1.SetMaximum(max(h1.GetMaximum(),h2.GetMaximum())*1.3)
    else:
        h1.SetMaximum(h2.GetMaximum()*1.3)
    h1.Draw("hist")
    h2.Draw("samehist")
    #h3.Draw("samehist")
    
    latex.DrawLatex(0.19, 0.90, "MET>80, p_{T}(jj)>80, p_{T}(j)>60,30, CSV(j)=CSVL, min #Delta#phi>0.5")
    if "DYJets" in file1:
        latex.DrawLatex(0.19, 0.85, "#color[1]{FullSim DYJets}")
        latex.DrawLatex(0.19, 0.80, "#color[4]{Smear/Histo DYJets}")
        gPad.Print("plots_zmmtozbb/zmmtozbb_compare_DYJets_%s.png" % hname)
        gPad.Print("plots_zmmtozbb/zmmtozbb_compare_DYJets_%s.pdf" % hname)
    else:
        latex.DrawLatex(0.19, 0.85, "#color[1]{FullSim}")
        latex.DrawLatex(0.19, 0.80, "#color[4]{Smear/Histo}")
        gPad.Print("plots_zmmtozbb/zmmtozbb_compare_ZbbHinv_%s.png" % hname)
        gPad.Print("plots_zmmtozbb/zmmtozbb_compare_ZbbHinv_%s.pdf" % hname)
    
    if "mjj" in hname:
        print h1.GetMean(), h1.GetRMS(), h2.GetMean(), h2.GetRMS()
    if "mindphi" in hname:
        print h1.Integral(1,h1.FindFixBin(0.5))/max(1e-6,h1.Integral()), h2.Integral(1,h2.FindFixBin(0.5))/max(1e-6,h2.Integral())
    if "jjmetdphi" in hname:
        print h1.Integral(1,h1.FindFixBin(2.0))/max(1e-6,h1.Integral()), h2.Integral(1,h2.FindFixBin(2.0))/max(1e-6,h2.Integral())
    
    if wait:  raw_input("Press Enter to continue...")
