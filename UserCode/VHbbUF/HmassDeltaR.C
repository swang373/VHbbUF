{
    treename = TString("tree_ZnunuHighPt_train");
    dir = TString("Step4_20130214/stitch/Step4_");
    
    _fileTT = TFile::Open(dir + "TT.root");
    _fileZj = TFile::Open(dir + "Zj.root");
    _fileZH = TFile::Open(dir + "ZH125.root");
    
    _treeTT = (TTree*) _fileTT->Get(treename);
    _treeZj = (TTree*) _fileZj->Get(treename);
    gROOT->cd(0);
    _treeZb = (TTree*) _treeZj->CopyTree("eventFlav==5");
    _treeZH = (TTree*) _fileZH->Get(treename);
    
    // dR_vs_Mjj
    h2TT = new TH2F("h2TT", "; M(jj) [GeV]; #DeltaR(jj)", 25, 0., 250., 44, 0, 2.2);
    h2Zb = new TH2F("h2Zb", "; M(jj) [GeV]; #DeltaR(jj)", 25, 0., 250., 44, 0, 2.2);
    h2ZH = new TH2F("h2ZH", "; M(jj) [GeV]; #DeltaR(jj)", 25, 0., 250., 44, 0, 2.2);
    lat = new TLatex();
    lat->SetNDC();
    
    // DJ
    _treeTT->Project("h2TT", "deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1]):HmassReg", "2.0 * weightsHCP[0] * selectFlags[0][0]", "goff");
    h2TT->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[632]{DJ} t #bar{t}");
    gPad->Print("plots/dR_vs_Mjj_DJ_TT.png");
    
    _treeZb->Project("h2Zb", "deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1]):HmassReg", "2.0 * weightsHCP[0] * selectFlags[0][0]", "goff");
    h2Zb->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[632]{DJ} Z+b");
    gPad->Print("plots/dR_vs_Mjj_DJ_Zb.png");
    
    _treeZH->Project("h2ZH", "deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1]):HmassReg", "2.0 * weightsHCP[0] * selectFlags[0][0]", "goff");
    h2ZH->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[632]{DJ} ZH(125)");
    gPad->Print("plots/dR_vs_Mjj_DJ_ZH.png");
    
    // 2FJ
    _treeTT->Project("h2TT", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):FatH.filteredmass", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && (nfathFilterJets==2 || (nfathFilterJets==3 && Alt$(fathFilterJets_pt[2],0)<=15.)))", "goff");
    h2TT->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[616]{2FJ} t #bar{t}");
    gPad->Print("plots/dR_vs_Mjj_2FJ_TT.png");
    
    _treeZb->Project("h2Zb", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):FatH.filteredmass", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && (nfathFilterJets==2 || (nfathFilterJets==3 && Alt$(fathFilterJets_pt[2],0)<=15.)))", "goff");
    h2Zb->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[616]{2FJ} Z+b");
    gPad->Print("plots/dR_vs_Mjj_2FJ_Zb.png");
    
    _treeZH->Project("h2ZH", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):FatH.filteredmass", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && (nfathFilterJets==2 || (nfathFilterJets==3 && Alt$(fathFilterJets_pt[2],0)<=15.)))", "goff");
    h2ZH->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[616]{2FJ} ZH(125)");
    gPad->Print("plots/dR_vs_Mjj_2FJ_ZH.png");
    
    // 3FJ
    _treeTT->Project("h2TT", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):FatH.filteredmass", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && (nfathFilterJets==3 && fathFilterJets_pt[2]>15.))", "goff");
    h2TT->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[616]{3FJ} t #bar{t}");
    gPad->Print("plots/dR_vs_Mjj_3FJ_TT.png");
    
    _treeZb->Project("h2Zb", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):FatH.filteredmass", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && (nfathFilterJets==3 && fathFilterJets_pt[2]>15.))", "goff");
    h2Zb->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[616]{3FJ} Z+b");
    gPad->Print("plots/dR_vs_Mjj_3FJ_Zb.png");
    
    _treeZH->Project("h2ZH", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):FatH.filteredmass", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && (nfathFilterJets==3 && fathFilterJets_pt[2]>15.))", "goff");
    h2ZH->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[616]{3FJ} ZH(125)");
    gPad->Print("plots/dR_vs_Mjj_3FJ_ZH.png");
    
    // 2FJReg
    _treeTT->Project("h2TT", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):FatHmassReg", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && (nfathFilterJets==2 || (nfathFilterJets==3 && Alt$(fathFilterJets_pt[2],0)<=15.)))", "goff");
    h2TT->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[900]{2FJReg} t #bar{t}");
    gPad->Print("plots/dR_vs_Mjj_2FJReg_TT.png");
    
    _treeZb->Project("h2Zb", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):FatHmassReg", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && (nfathFilterJets==2 || (nfathFilterJets==3 && Alt$(fathFilterJets_pt[2],0)<=15.)))", "goff");
    h2Zb->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[900]{2FJReg} Z+b");
    gPad->Print("plots/dR_vs_Mjj_2FJReg_Zb.png");
    
    _treeZH->Project("h2ZH", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):FatHmassReg", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && (nfathFilterJets==2 || (nfathFilterJets==3 && Alt$(fathFilterJets_pt[2],0)<=15.)))", "goff");
    h2ZH->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[900]{2FJReg} ZH(125)");
    gPad->Print("plots/dR_vs_Mjj_2FJReg_ZH.png");
    
    // 3FJReg
    _treeTT->Project("h2TT", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):FatHmassReg", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && (nfathFilterJets==3 && fathFilterJets_pt[2]>15.))", "goff");
    h2TT->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[900]{3FJReg} t #bar{t}");
    gPad->Print("plots/dR_vs_Mjj_3FJReg_TT.png");
    
    _treeZb->Project("h2Zb", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):FatHmassReg", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && (nfathFilterJets==3 && fathFilterJets_pt[2]>15.))", "goff");
    h2Zb->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[900]{3FJReg} Z+b");
    gPad->Print("plots/dR_vs_Mjj_3FJReg_Zb.png");
    
    _treeZH->Project("h2ZH", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):FatHmassReg", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && (nfathFilterJets==3 && fathFilterJets_pt[2]>15.))", "goff");
    h2ZH->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[900]{3FJReg} ZH(125)");
    gPad->Print("plots/dR_vs_Mjj_3FJReg_ZH.png");
    
    
    // dR_vs_dR
    h2TT = new TH2F("h2TT", "; DJ #DeltaR(jj) [GeV]; FJ #DeltaR(jj)", 44, 0, 2.2., 44, 0, 2.2);
    h2Zb = new TH2F("h2Zb", "; DJ #DeltaR(jj) [GeV]; FJ #DeltaR(jj)", 44, 0, 2.2., 44, 0, 2.2);
    h2ZH = new TH2F("h2ZH", "; DJ #DeltaR(jj) [GeV]; FJ #DeltaR(jj)", 44, 0, 2.2., 44, 0, 2.2);
    
    // 2FJ
    _treeTT->Project("h2TT", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && (nfathFilterJets==2 || (nfathFilterJets==3 && Alt$(fathFilterJets_pt[2],0)<=15.)))", "goff");
    h2TT->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[616]{2FJ} t #bar{t}");
    gPad->Print("plots/dR_vs_dR_2FJ_TT.png");
    
    _treeZb->Project("h2Zb", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && (nfathFilterJets==2 || (nfathFilterJets==3 && Alt$(fathFilterJets_pt[2],0)<=15.)))", "goff");
    h2Zb->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[616]{2FJ} Z+b");
    gPad->Print("plots/dR_vs_dR_2FJ_Zb.png");
    
    _treeZH->Project("h2ZH", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && (nfathFilterJets==2 || (nfathFilterJets==3 && Alt$(fathFilterJets_pt[2],0)<=15.)))", "goff");
    h2ZH->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[616]{2FJ} ZH(125)");
    gPad->Print("plots/dR_vs_dR_2FJ_ZH.png");
    
    // 3FJ
    _treeTT->Project("h2TT", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && (nfathFilterJets==3 && fathFilterJets_pt[2]>15.))", "goff");
    h2TT->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[616]{3FJ} t #bar{t}");
    gPad->Print("plots/dR_vs_dR_3FJ_TT.png");
    
    _treeZb->Project("h2Zb", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && (nfathFilterJets==3 && fathFilterJets_pt[2]>15.))", "goff");
    h2Zb->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[616]{3FJ} Z+b");
    gPad->Print("plots/dR_vs_dR_3FJ_Zb.png");
    
    _treeZH->Project("h2ZH", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && (nfathFilterJets==3 && fathFilterJets_pt[2]>15.))", "goff");
    h2ZH->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[616]{3FJ} ZH(125)");
    gPad->Print("plots/dR_vs_dR_3FJ_ZH.png");
    
    // naJets_Znn==0
    _treeTT->Project("h2TT", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && naJets_Znn==0)", "goff");
    h2TT->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[432]{naJets==0} t #bar{t}");
    gPad->Print("plots/dR_vs_dR_naJets0_TT.png");
    
    _treeZb->Project("h2Zb", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && naJets_Znn==0)", "goff");
    h2Zb->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[432]{naJets==0} Z+b");
    gPad->Print("plots/dR_vs_dR_naJets0_Zb.png");
    
    _treeZH->Project("h2ZH", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && naJets_Znn==0)", "goff");
    h2ZH->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[432]{naJets==0} ZH(125)");
    gPad->Print("plots/dR_vs_dR_naJets0_ZH.png");
    
    // naJets_Znn>0
    _treeTT->Project("h2TT", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && naJets_Znn>0)", "goff");
    h2TT->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[432]{naJets>0} t #bar{t}");
    gPad->Print("plots/dR_vs_dR_naJets1_TT.png");
    
    _treeZb->Project("h2Zb", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && naJets_Znn>0)", "goff");
    h2Zb->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[432]{naJets>0} Z+b");
    gPad->Print("plots/dR_vs_dR_naJets1_Zb.png");
    
    _treeZH->Project("h2ZH", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && naJets_Znn>0)", "goff");
    h2ZH->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[432]{naJets>0} ZH(125)");
    gPad->Print("plots/dR_vs_dR_naJets1_ZH.png");
    
    // naJets_Znn==0, H pT < 250
    _treeZH->Project("h2ZH", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && HptReg<250 && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && naJets_Znn==0)", "goff");
    h2ZH->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[432]{naJets==0} #color[416]{H p_{T}<250} ZH(125)");
    gPad->Print("plots/dR_vs_dR_naJets0_HptReg_ZH.png");
    
    // naJets_Znn>0, H pT < 250
    _treeZH->Project("h2ZH", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1]):deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && HptReg<250 && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && fathFilterJets_pt[1]>15 && naJets_Znn>0)", "goff");
    h2ZH->Draw("COLZ");
    lat->DrawLatex(0.2, 0.88, "#color[432]{naJets>0} #color[416]{H p_{T}<250} ZH(125)");
    gPad->Print("plots/dR_vs_dR_naJets1_HptReg_ZH.png");
    
    
    // Mjj 1D
    h1TT_0 = new TH1F("h1TT_0", "; M(jj) [GeV]", 25, 0., 250.);
    h1Zb_0 = new TH1F("h1Zb_0", "; M(jj) [GeV]", 25, 0., 250.);
    h1ZH_0 = new TH1F("h1ZH_0", "; M(jj) [GeV]", 25, 0., 250.);
    h1TT_2 = new TH1F("h1TT_2", "; M(jj) [GeV]", 25, 0., 250.);
    h1Zb_2 = new TH1F("h1Zb_2", "; M(jj) [GeV]", 25, 0., 250.);
    h1ZH_2 = new TH1F("h1ZH_2", "; M(jj) [GeV]", 25, 0., 250.);
    h1TT_3 = new TH1F("h1TT_3", "; M(jj) [GeV]", 25, 0., 250.);
    h1Zb_3 = new TH1F("h1Zb_3", "; M(jj) [GeV]", 25, 0., 250.);
    h1ZH_3 = new TH1F("h1ZH_3", "; M(jj) [GeV]", 25, 0., 250.);
    h1TT_0->SetStats(0); h1TT_0->SetLineWidth(2); h1TT_0->SetLineColor(kBlack);
    h1TT_2->SetStats(0); h1TT_2->SetLineWidth(2); h1TT_2->SetLineColor(kMagenta);
    h1TT_3->SetStats(0); h1TT_3->SetLineWidth(2); h1TT_3->SetLineColor(kGreen);
    h1Zb_0->SetStats(0); h1Zb_0->SetLineWidth(2); h1Zb_0->SetLineColor(kBlack);
    h1Zb_2->SetStats(0); h1Zb_2->SetLineWidth(2); h1Zb_2->SetLineColor(kMagenta);
    h1Zb_3->SetStats(0); h1Zb_3->SetLineWidth(2); h1Zb_3->SetLineColor(kGreen);
    h1ZH_0->SetStats(0); h1ZH_0->SetLineWidth(2); h1ZH_0->SetLineColor(kBlack);
    h1ZH_2->SetStats(0); h1ZH_2->SetLineWidth(2); h1ZH_2->SetLineColor(kMagenta);
    h1ZH_3->SetStats(0); h1ZH_3->SetLineWidth(2); h1ZH_3->SetLineColor(kGreen);
    leg = new TLegend(0.7,0.8,0.9,0.9);
    leg->AddEntry(h1TT_0, "0FJ", "l");
    leg->AddEntry(h1TT_2, "2FJ", "l");
    leg->AddEntry(h1TT_3, "3FJ", "l");
    
    _treeTT->Project("h1TT_0", "HmassReg", "2.0 * weightsHCP[0] * (selectFlags[0][0] && !(nfathFilterJets>0 && FatH.FatHiggsFlag==1))", "goff");
    _treeTT->Project("h1TT_2", "HmassReg", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && (nfathFilterJets==2 || (nfathFilterJets==3 && Alt$(fathFilterJets_pt[2],0)<=15.)))", "goff");
    _treeTT->Project("h1TT_3", "HmassReg", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && (nfathFilterJets==3 && fathFilterJets_pt[2]>15.))", "goff");
    h1TT_0->SetMaximum(h1TT_0->GetMaximum()*1.5);
    h1TT_0->Draw();
    h1TT_2->Draw("same");
    h1TT_3->Draw("same");
    leg->Draw();
    lat->DrawLatex(0.2, 0.88, "#color[632]{DJ} t #bar{t}");
    gPad->Print("plots/Mjj_023FJ_DJ_TT.png");
    
    _treeZb->Project("h1Zb_0", "HmassReg", "2.0 * weightsHCP[0] * (selectFlags[0][0] && !(nfathFilterJets>0 && FatH.FatHiggsFlag==1))", "goff");
    _treeZb->Project("h1Zb_2", "HmassReg", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && (nfathFilterJets==2 || (nfathFilterJets==3 && Alt$(fathFilterJets_pt[2],0)<=15.)))", "goff");
    _treeZb->Project("h1Zb_3", "HmassReg", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && (nfathFilterJets==3 && fathFilterJets_pt[2]>15.))", "goff");
    h1Zb_0->SetMaximum(h1Zb_0->GetMaximum()*1.5);
    h1Zb_0->Draw();
    h1Zb_2->Draw("same");
    h1Zb_3->Draw("same");
    leg->Draw();
    lat->DrawLatex(0.2, 0.88, "#color[632]{DJ} Z+b");
    gPad->Print("plots/Mjj_023FJ_DJ_Zb.png");
    
    _treeZH->Project("h1ZH_0", "HmassReg", "2.0 * weightsHCP[0] * (selectFlags[0][0] && !(nfathFilterJets>0 && FatH.FatHiggsFlag==1))", "goff");
    _treeZH->Project("h1ZH_2", "HmassReg", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && (nfathFilterJets==2 || (nfathFilterJets==3 && Alt$(fathFilterJets_pt[2],0)<=15.)))", "goff");
    _treeZH->Project("h1ZH_3", "HmassReg", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && (nfathFilterJets==3 && fathFilterJets_pt[2]>15.))", "goff");
    h1ZH_0->SetMaximum(h1ZH_0->GetMaximum()*1.5);
    h1ZH_0->Draw();
    h1ZH_2->Draw("same");
    h1ZH_3->Draw("same");
    leg->Draw();
    lat->DrawLatex(0.2, 0.88, "#color[632]{DJ} ZH(125)");
    gPad->Print("plots/Mjj_023FJ_DJ_ZH.png");
    
    // filtered jets
    _treeTT->Project("h1TT_2", "FatH.filteredmass", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && (nfathFilterJets==2 || (nfathFilterJets==3 && Alt$(fathFilterJets_pt[2],0)<=15.)))", "goff");
    _treeTT->Project("h1TT_3", "FatH.filteredmass", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && (nfathFilterJets==3 && fathFilterJets_pt[2]>15.))", "goff");
    h1TT_2->SetMaximum(h1TT_0->GetMaximum());
    h1TT_2->Draw();
    h1TT_3->Draw("same");
    leg->Draw();
    lat->DrawLatex(0.2, 0.88, "#color[616]{FJ} t #bar{t}");
    gPad->Print("plots/Mjj_023FJ_FJ_TT.png");
    
    _treeZb->Project("h1Zb_2", "FatH.filteredmass", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && (nfathFilterJets==2 || (nfathFilterJets==3 && Alt$(fathFilterJets_pt[2],0)<=15.)))", "goff");
    _treeZb->Project("h1Zb_3", "FatH.filteredmass", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && (nfathFilterJets==3 && fathFilterJets_pt[2]>15.))", "goff");
    h1Zb_2->SetMaximum(h1Zb_0->GetMaximum());
    h1Zb_2->Draw();
    h1Zb_3->Draw("same");
    leg->Draw();
    lat->DrawLatex(0.2, 0.88, "#color[616]{FJ} Z+b");
    gPad->Print("plots/Mjj_023FJ_FJ_Zb.png");
    
    _treeZH->Project("h1ZH_2", "FatH.filteredmass", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && (nfathFilterJets==2 || (nfathFilterJets==3 && Alt$(fathFilterJets_pt[2],0)<=15.)))", "goff");
    _treeZH->Project("h1ZH_3", "FatH.filteredmass", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && (nfathFilterJets==3 && fathFilterJets_pt[2]>15.))", "goff");
    h1ZH_2->SetMaximum(h1ZH_0->GetMaximum());
    h1ZH_2->Draw();
    h1ZH_3->Draw("same");
    leg->Draw();
    lat->DrawLatex(0.2, 0.88, "#color[616]{FJ} ZH(125)");
    gPad->Print("plots/Mjj_023FJ_FJ_ZH.png");
    
    
    // dR 1D
    h1TT_0 = new TH1F("h1TT_0", "; #DeltaR(jj) [GeV]", 44, 0., 2.2);
    h1Zb_0 = new TH1F("h1Zb_0", "; #DeltaR(jj) [GeV]", 44, 0., 2.2);
    h1ZH_0 = new TH1F("h1ZH_0", "; #DeltaR(jj) [GeV]", 44, 0., 2.2);
    h1TT_2 = new TH1F("h1TT_2", "; #DeltaR(jj) [GeV]", 44, 0., 2.2);
    h1Zb_2 = new TH1F("h1Zb_2", "; #DeltaR(jj) [GeV]", 44, 0., 2.2);
    h1ZH_2 = new TH1F("h1ZH_2", "; #DeltaR(jj) [GeV]", 44, 0., 2.2);
    h1TT_3 = new TH1F("h1TT_3", "; #DeltaR(jj) [GeV]", 44, 0., 2.2);
    h1Zb_3 = new TH1F("h1Zb_3", "; #DeltaR(jj) [GeV]", 44, 0., 2.2);
    h1ZH_3 = new TH1F("h1ZH_3", "; #DeltaR(jj) [GeV]", 44, 0., 2.2);
    h1TT_0->SetStats(0); h1TT_0->SetLineWidth(2); h1TT_0->SetLineColor(kBlack);
    h1TT_2->SetStats(0); h1TT_2->SetLineWidth(2); h1TT_2->SetLineColor(kMagenta);
    h1TT_3->SetStats(0); h1TT_3->SetLineWidth(2); h1TT_3->SetLineColor(kGreen);
    h1Zb_0->SetStats(0); h1Zb_0->SetLineWidth(2); h1Zb_0->SetLineColor(kBlack);
    h1Zb_2->SetStats(0); h1Zb_2->SetLineWidth(2); h1Zb_2->SetLineColor(kMagenta);
    h1Zb_3->SetStats(0); h1Zb_3->SetLineWidth(2); h1Zb_3->SetLineColor(kGreen);
    h1ZH_0->SetStats(0); h1ZH_0->SetLineWidth(2); h1ZH_0->SetLineColor(kBlack);
    h1ZH_2->SetStats(0); h1ZH_2->SetLineWidth(2); h1ZH_2->SetLineColor(kMagenta);
    h1ZH_3->SetStats(0); h1ZH_3->SetLineWidth(2); h1ZH_3->SetLineColor(kGreen);
    leg = new TLegend(0.7,0.8,0.9,0.9);
    leg->AddEntry(h1TT_0, "0FJ", "l");
    leg->AddEntry(h1TT_2, "2FJ", "l");
    leg->AddEntry(h1TT_3, "3FJ", "l");
    
    _treeTT->Project("h1TT_0", "deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && !(nfathFilterJets>0 && FatH.FatHiggsFlag==1))", "goff");
    _treeTT->Project("h1TT_2", "deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && (nfathFilterJets==2 || (nfathFilterJets==3 && Alt$(fathFilterJets_pt[2],0)<=15.)))", "goff");
    _treeTT->Project("h1TT_3", "deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && (nfathFilterJets==3 && fathFilterJets_pt[2]>15.))", "goff");
    h1TT_0->SetMaximum(h1TT_0->GetMaximum()*1.5);
    h1TT_0->Draw();
    h1TT_2->Draw("same");
    h1TT_3->Draw("same");
    leg->Draw();
    lat->DrawLatex(0.2, 0.88, "#color[632]{DJ} t #bar{t}");
    gPad->Print("plots/dR_023FJ_DJ_TT.png");
    
    _treeZb->Project("h1Zb_0", "deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && !(nfathFilterJets>0 && FatH.FatHiggsFlag==1))", "goff");
    _treeZb->Project("h1Zb_2", "deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && (nfathFilterJets==2 || (nfathFilterJets==3 && Alt$(fathFilterJets_pt[2],0)<=15.)))", "goff");
    _treeZb->Project("h1Zb_3", "deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && (nfathFilterJets==3 && fathFilterJets_pt[2]>15.))", "goff");
    h1Zb_0->SetMaximum(h1Zb_0->GetMaximum()*1.5);
    h1Zb_0->Draw();
    h1Zb_2->Draw("same");
    h1Zb_3->Draw("same");
    leg->Draw();
    lat->DrawLatex(0.2, 0.88, "#color[632]{DJ} Z+b");
    gPad->Print("plots/dR_023FJ_DJ_Zb.png");
    
    _treeZH->Project("h1ZH_0", "deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && !(nfathFilterJets>0 && FatH.FatHiggsFlag==1))", "goff");
    _treeZH->Project("h1ZH_2", "deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && (nfathFilterJets==2 || (nfathFilterJets==3 && Alt$(fathFilterJets_pt[2],0)<=15.)))", "goff");
    _treeZH->Project("h1ZH_3", "deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && (nfathFilterJets==3 && fathFilterJets_pt[2]>15.))", "goff");
    h1ZH_0->SetMaximum(h1ZH_0->GetMaximum()*1.5);
    h1ZH_0->Draw();
    h1ZH_2->Draw("same");
    h1ZH_3->Draw("same");
    leg->Draw();
    lat->DrawLatex(0.2, 0.88, "#color[632]{DJ} ZH(125)");
    gPad->Print("plots/dR_023FJ_DJ_ZH.png");
    
    // filtered jets
    _treeTT->Project("h1TT_2", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && (nfathFilterJets==2 || (nfathFilterJets==3 && Alt$(fathFilterJets_pt[2],0)<=15.)))", "goff");
    _treeTT->Project("h1TT_3", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && (nfathFilterJets==3 && fathFilterJets_pt[2]>15.))", "goff");
    h1TT_2->SetMaximum(h1TT_0->GetMaximum());
    h1TT_2->Draw();
    h1TT_3->Draw("same");
    leg->Draw();
    lat->DrawLatex(0.2, 0.88, "#color[616]{FJ} t #bar{t}");
    gPad->Print("plots/dR_023FJ_FJ_TT.png");
    
    _treeZb->Project("h1Zb_2", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && (nfathFilterJets==2 || (nfathFilterJets==3 && Alt$(fathFilterJets_pt[2],0)<=15.)))", "goff");
    _treeZb->Project("h1Zb_3", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && (nfathFilterJets==3 && fathFilterJets_pt[2]>15.))", "goff");
    h1Zb_2->SetMaximum(h1Zb_0->GetMaximum());
    h1Zb_2->Draw();
    h1Zb_3->Draw("same");
    leg->Draw();
    lat->DrawLatex(0.2, 0.88, "#color[616]{FJ} Z+b");
    gPad->Print("plots/dR_023FJ_FJ_Zb.png");
    
    _treeZH->Project("h1ZH_2", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && (nfathFilterJets==2 || (nfathFilterJets==3 && Alt$(fathFilterJets_pt[2],0)<=15.)))", "goff");
    _treeZH->Project("h1ZH_3", "deltaR(fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_eta[1], fathFilterJets_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && nfathFilterJets>0 && FatH.FatHiggsFlag==1 && (nfathFilterJets==3 && fathFilterJets_pt[2]>15.))", "goff");
    h1ZH_2->SetMaximum(h1ZH_0->GetMaximum());
    h1ZH_2->Draw();
    h1ZH_3->Draw("same");
    leg->Draw();
    lat->DrawLatex(0.2, 0.88, "#color[616]{FJ} ZH(125)");
    gPad->Print("plots/dR_023FJ_FJ_ZH.png");
    
    // DJ (deltaR < 1.2)
    h1TT_0 = new TH1F("h1TT_0", "; M(jj) [GeV]", 25, 0., 250.);
    h1Zb_0 = new TH1F("h1Zb_0", "; M(jj) [GeV]", 25, 0., 250.);
    h1ZH_0 = new TH1F("h1ZH_0", "; M(jj) [GeV]", 25, 0., 250.);
    h1TT_2 = new TH1F("h1TT_2", "; M(jj) [GeV]", 25, 0., 250.);
    h1Zb_2 = new TH1F("h1Zb_2", "; M(jj) [GeV]", 25, 0., 250.);
    h1ZH_2 = new TH1F("h1ZH_2", "; M(jj) [GeV]", 25, 0., 250.);
    h1TT_0->SetStats(0); h1TT_0->SetLineWidth(2); h1TT_0->SetLineColor(kBlack);
    h1TT_2->SetStats(0); h1TT_2->SetLineWidth(2); h1TT_2->SetLineColor(kBlue); h1TT_2->SetFillColor(kBlue-7);
    h1Zb_0->SetStats(0); h1Zb_0->SetLineWidth(2); h1Zb_0->SetLineColor(kBlack);
    h1Zb_2->SetStats(0); h1Zb_2->SetLineWidth(2); h1Zb_2->SetLineColor(kBlue); h1Zb_2->SetFillColor(kBlue-7);
    h1ZH_0->SetStats(0); h1ZH_0->SetLineWidth(2); h1ZH_0->SetLineColor(kBlack);
    h1ZH_2->SetStats(0); h1ZH_2->SetLineWidth(2); h1ZH_2->SetLineColor(kBlue); h1ZH_2->SetFillColor(kBlue-7);
    leg = new TLegend(0.7,0.8,0.9,0.9);
    leg->AddEntry(h1TT_0, "all DJ", "lf");
    leg->AddEntry(h1TT_2, "DJ, #DeltaR < 1.2", "lf");
    
    _treeTT->Project("h1TT_0", "HmassReg", "2.0 * weightsHCP[0] * (selectFlags[0][0])", "goff");
    _treeTT->Project("h1TT_2", "HmassReg", "2.0 * weightsHCP[0] * (selectFlags[0][0] && deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])<1.2)", "goff");
    h1TT_0->Draw();
    h1TT_2->Draw("same");
    leg->Draw();
    lat->DrawLatex(0.2, 0.88, "#color[632]{DJ} t #bar{t}");
    gPad->Print("plots/dR_cut_DJ_TT.png");
    
    _treeZb->Project("h1Zb_0", "HmassReg", "2.0 * weightsHCP[0] * (selectFlags[0][0])", "goff");
    _treeZb->Project("h1Zb_2", "HmassReg", "2.0 * weightsHCP[0] * (selectFlags[0][0] && deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])<1.2)", "goff");
    h1Zb_0->Draw();
    h1Zb_2->Draw("same");
    leg->Draw();
    lat->DrawLatex(0.2, 0.88, "#color[632]{DJ} Z+b");
    gPad->Print("plots/dR_cut_DJ_Zb.png");
    
    _treeZH->Project("h1ZH_0", "HmassReg", "2.0 * weightsHCP[0] * (selectFlags[0][0])", "goff");
    _treeZH->Project("h1ZH_2", "HmassReg", "2.0 * weightsHCP[0] * (selectFlags[0][0] && deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])<1.2)", "goff");
    h1ZH_0->Draw();
    h1ZH_2->Draw("same");
    leg->Draw();
    lat->DrawLatex(0.2, 0.88, "#color[632]{DJ} ZH(125)");
    gPad->Print("plots/dR_cut_DJ_ZH.png");
    
    // DJ + (deltaR < 1.2)
    h1TT_0 = new TH1F("h1TT_0", "; #DeltaR(jj) [GeV]", 44, 0., 2.2);
    h1Zb_0 = new TH1F("h1Zb_0", "; #DeltaR(jj) [GeV]", 44, 0., 2.2);
    h1ZH_0 = new TH1F("h1ZH_0", "; #DeltaR(jj) [GeV]", 44, 0., 2.2);
    h1TT_2 = new TH1F("h1TT_2", "; #DeltaR(jj) [GeV]", 44, 0., 2.2);
    h1Zb_2 = new TH1F("h1Zb_2", "; #DeltaR(jj) [GeV]", 44, 0., 2.2);
    h1ZH_2 = new TH1F("h1ZH_2", "; #DeltaR(jj) [GeV]", 44, 0., 2.2);
    h1TT_0->SetStats(0); h1TT_0->SetLineWidth(2); h1TT_0->SetLineColor(kBlack);
    h1TT_2->SetStats(0); h1TT_2->SetLineWidth(2); h1TT_2->SetLineColor(kRed); h1TT_2->SetFillColor(kRed-7);
    h1Zb_0->SetStats(0); h1Zb_0->SetLineWidth(2); h1Zb_0->SetLineColor(kBlack);
    h1Zb_2->SetStats(0); h1Zb_2->SetLineWidth(2); h1Zb_2->SetLineColor(kRed); h1Zb_2->SetFillColor(kRed-7);
    h1ZH_0->SetStats(0); h1ZH_0->SetLineWidth(2); h1ZH_0->SetLineColor(kBlack);
    h1ZH_2->SetStats(0); h1ZH_2->SetLineWidth(2); h1ZH_2->SetLineColor(kRed); h1ZH_2->SetFillColor(kRed-7);
    leg = new TLegend(0.55,0.8,0.9,0.9);
    leg->AddEntry(h1TT_0, "all DJ", "lf");
    leg->AddEntry(h1TT_2, "DJ, # add. jets == 0", "lf");
    
    _treeTT->Project("h1TT_0", "deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0])", "goff");
    _treeTT->Project("h1TT_2", "deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && naJets_Znn==0)", "goff");
    h1TT_0->Draw();
    h1TT_2->Draw("same");
    leg->Draw();
    lat->DrawLatex(0.2, 0.88, "#color[632]{DJ} t #bar{t}");
    gPad->Print("plots/dR_naJets0_DJ_TT.png");
    
    _treeZb->Project("h1Zb_0", "deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0])", "goff");
    _treeZb->Project("h1Zb_2", "deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && naJets_Znn==0)", "goff");
    h1Zb_0->Draw();
    h1Zb_2->Draw("same");
    leg->Draw();
    lat->DrawLatex(0.2, 0.88, "#color[632]{DJ} Z+b");
    gPad->Print("plots/dR_naJets0_DJ_Zb.png");
    
    _treeZH->Project("h1ZH_0", "deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0])", "goff");
    _treeZH->Project("h1ZH_2", "deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", "2.0 * weightsHCP[0] * (selectFlags[0][0] && naJets_Znn==0)", "goff");
    h1ZH_0->Draw();
    h1ZH_2->Draw("same");
    leg->Draw();
    lat->DrawLatex(0.2, 0.88, "#color[632]{DJ} ZH(125)");
    gPad->Print("plots/dR_naJets0_DJ_ZH.png");
    
    
    
}
