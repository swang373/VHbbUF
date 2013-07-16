{
// from TMVA
const Int_t c_Canvas         = TColor::GetColor( "#f0f0f0" );
const Int_t c_FrameFill      = TColor::GetColor( "#fffffd" );
const Int_t c_TitleBox       = TColor::GetColor( "#5D6B7D" );
const Int_t c_TitleBorder    = TColor::GetColor( "#7D8B9D" );
const Int_t c_TitleText      = TColor::GetColor( "#FFFFFF" );
const Int_t c_SignalLine     = TColor::GetColor( "#0000ee" );
const Int_t c_SignalFill     = TColor::GetColor( "#7d99d1" );
const Int_t c_BackgroundLine = TColor::GetColor( "#ff0000" );
const Int_t c_BackgroundFill = TColor::GetColor( "#ff0000" );
const Int_t c_NovelBlue      = TColor::GetColor( "#2244a5" );
TLatex* text = new TLatex(); text->SetNDC(); text->SetTextSize(0.025);

//gROOT->SetBatch();

// Get trees
TFile* _file0 = TFile::Open("Step3_20130130/Step3_ZnnH125.root");
TTree* tree = (TTree*) _file0->Get("tree");

TFile* _fileData = TFile::Open("Step3_20130130/Step3_Data.root");
TTree* treeData = (TTree*) _fileData->Get("tree");

//TFile* _fileTT = TFile::Open("skim/Step3_TTPowheg.root");
TFile* _fileTT = TFile::Open("Step3_20130130/Step3_TTPowheg.root");
TTree* treeTT = (TTree*) _fileTT->Get("tree");

TFile* _fileZZ = TFile::Open("Step3_20130130/Step3_ZZ.root");
TTree* treeZZ = (TTree*) _fileZZ->Get("tree");

TChain* chainZJ = new TChain("tree");
chainZJ->Add("Step3_20130130/Step3_ZJetsHT100.root");
chainZJ->Add("Step3_20130130/Step3_ZJetsHT200.root");
chainZJ->Add("Step3_20130130/Step3_ZJetsHT400.root");

TChain* chainVV = new TChain("tree");
chainVV->Add("Step3_20130130/Step3_WW.root");
chainVV->Add("Step3_20130130/Step3_WZ.root");
chainVV->Add("Step3_20130130/Step3_ZZ.root");

TChain* chainST = new TChain("tree");
chainST->Add("Step3_20130130/Step3_T_s.root");
chainST->Add("Step3_20130130/Step3_T_t.root");
chainST->Add("Step3_20130130/Step3_T_tW.root");
chainST->Add("Step3_20130130/Step3_Tbar_s.root");
chainST->Add("Step3_20130130/Step3_Tbar_t.root");
chainST->Add("Step3_20130130/Step3_Tbar_tW.root");


//- H mass ---------------------------------------------------------------------
//TCut cut = "(Vtype==4||Vtype==3||Vtype==2) && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5 && hJet_id[0] && hJet_id[1] && hJet_pt[0]>20 && hJet_pt[1]>20 && hJet_csv[0]>0 && hJet_csv[1]>0 && H.pt>80 && METtype1corr.et>100";
TCut cut = "Vtype==4 && hJet_pt[0]>60 && hJet_pt[1]>30 && METtype1corr.et>170 && H.pt>130 && hJet_csv_nominal[0]>0.244 && hJet_csv_nominal[1]>0.244 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && H.mass>0";
//TCut cut = "Vtype==4 && hJet_ptReg[0]>80 && hJet_ptReg[1]>30 && HmassReg<250 && HptReg>160 && METtype1corr.et>160 && max(hJet_csvCor[0],hJet_csvCor[1])>0.679 && min(hJet_csvCor[0],hJet_csvCor[1])>0.244 && HVdPhi>2.0 && min(Min$(abs(deltaPhi(METtype1corr.phi,hJet_phi))),Min$(abs(deltaPhiMETjets(METtype1corr.phi,aJet_phi,aJet_pt,aJet_eta)))+999*(Sum$(aJet_pt>30 && abs(aJet_eta)<2.5)==0) )>0.5 && abs(deltaPhi(METtype1corr.phi,METnoPUCh.phi))<0.5";

/*
TH1::SetDefaultSumw2();
double xmin=110., xmax=140.;
double xxmin=80., xxmax=110.;
tree->Draw("H.mass >> h1(140,60,200)",cut);
tree->Draw("HmassReg >> h2(140,60,200)",cut);
treeZZ->Draw("H.mass >> h3(140,60,200)",cut);
treeZZ->Draw("HmassReg >> h4(140,60,200)",cut);
h1=h1;h2=h2;h3=h3;h4=h4;
h4->Scale(1.0/h4->Integral());
h3->Scale(1.0/h3->Integral());
h2->Scale(1.0/h2->Integral());
h1->Scale(1.0/h1->Integral());

double varh1[8], varh2[8], varh3[8], varh4[8];
h2->Fit("gaus", "I", "", xmin, xmax);
varh2[0] = gaus->GetParameter(0); varh2[1] = gaus->GetParameter(1); varh2[2] = gaus->GetParameter(2); varh2[3] = gaus->GetParError(0); varh2[4] = gaus->GetParError(1); varh2[5] = gaus->GetParError(2); varh2[6] = gaus->GetChisquare(); varh2[7] = gaus->GetNDF();
h4->Fit("gaus", "I", "", xxmin, xxmax);
varh4[0] = gaus->GetParameter(0); varh4[1] = gaus->GetParameter(1); varh4[2] = gaus->GetParameter(2); varh4[3] = gaus->GetParError(0); varh4[4] = gaus->GetParError(1); varh4[5] = gaus->GetParError(2); varh4[6] = gaus->GetChisquare(); varh4[7] = gaus->GetNDF();

h4->SetStats(0); h4->SetLineWidth(2); h4->SetMarkerSize(0); h4->SetLineColor(kGray+1);
h3->SetStats(0); h3->SetLineWidth(1); h3->SetMarkerSize(0); h3->SetLineColor(kBlack);
h2->SetStats(0); h2->SetLineWidth(2); h2->SetMarkerSize(0); h2->SetLineColor(kRed);
h1->SetStats(0); h1->SetLineWidth(1); h1->SetMarkerSize(0); h1->SetLineColor(kRed+2);
h3->SetTitle("; H(125) mass [GeV]; (normalized)"); h3->Draw("hist");
h4->Draw("histsame");
h2->Draw("histsame");
h1->Draw("histsame");
h4->GetFunction("gaus")->Draw("same");
h2->GetFunction("gaus")->Draw("same");

TLegend* leg = new TLegend(0.62, 0.78, 0.92, 0.92);
leg->AddEntry(h4, "ZZ regressed");
leg->AddEntry(h3, "ZZ");
leg->AddEntry(h2, "ZH(125) regressed");
leg->AddEntry(h1, "ZH(125)");
leg->Draw();

double down = 0.12;
TLatex* texth2 = new TLatex(); texth2->SetNDC(); texth2->SetTextSize(0.025); texth2->SetTextColor(kRed);
texth2->DrawLatex(0.61, 0.74-down, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth2->DrawLatex(0.61, 0.71-down, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh2[1], varh2[1+3], varh2[2], varh2[2+3]));
texth2->DrawLatex(0.61, 0.68-down, Form("#chi^{2}/NDF = %.1f / %.0f", varh2[6], varh2[7]));
TLatex* texth4 = new TLatex(); texth4->SetNDC(); texth4->SetTextSize(0.025); texth4->SetTextColor(kGray+1);
texth4->DrawLatex(0.61, 0.74, Form("gaus fit [%.0f,%.0f]", xxmin, xxmax));
texth4->DrawLatex(0.61, 0.71, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh4[1], varh4[1+3], varh4[2], varh4[2+3]));
texth4->DrawLatex(0.61, 0.68, Form("#chi^{2}/NDF = %.1f / %.0f", varh4[6], varh4[7]));
text->DrawLatex(0.015, 0.97, "H pT > 130 GeV, MET > 170 GeV, 2 CSVL");
gPad->Print("plots/breg_HmassRegSB.png");

h1->Fit("gaus", "I", "", xmin, xmax); 
varh1[0] = gaus->GetParameter(0); varh1[1] = gaus->GetParameter(1); varh1[2] = gaus->GetParameter(2); varh1[3] = gaus->GetParError(0); varh1[4] = gaus->GetParError(1); varh1[5] = gaus->GetParError(2); varh1[6] = gaus->GetChisquare(); varh1[7] = gaus->GetNDF();
h3->Fit("gaus", "I", "", xxmin, xxmax);
varh3[0] = gaus->GetParameter(0); varh3[1] = gaus->GetParameter(1); varh3[2] = gaus->GetParameter(2); varh3[3] = gaus->GetParError(0); varh3[4] = gaus->GetParError(1); varh3[5] = gaus->GetParError(2); varh3[6] = gaus->GetChisquare(); varh3[7] = gaus->GetNDF();
h3->SetStats(0); h3->SetLineWidth(2); h3->SetMarkerSize(0); h3->SetLineColor(kBlack);
h1->SetStats(0); h1->SetLineWidth(2); h1->SetMarkerSize(0); h1->SetLineColor(kRed+2);
h3->SetTitle("; H(125) mass [GeV]; (normalized)"); h3->Draw("hist");
h1->Draw("histsame");
h3->GetFunction("gaus")->Draw("same");
h1->GetFunction("gaus")->Draw("same");

TLegend* leg = new TLegend(0.62, 0.85, 0.92, 0.92);
leg->AddEntry(h3, "ZZ");
leg->AddEntry(h1, "ZH(125)");
leg->Draw();

TLatex* texth1 = new TLatex(); texth1->SetNDC(); texth1->SetTextSize(0.025); texth1->SetTextColor(kBlack);
texth1->DrawLatex(0.61, 0.74-down, Form("gaus fit [%.0f,%.0f]", xxmin, xxmax));
texth1->DrawLatex(0.61, 0.71-down, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh1[1], varh1[1+3], varh1[2], varh1[2+3]));
texth1->DrawLatex(0.61, 0.68-down, Form("#chi^{2}/NDF = %.1f / %.0f", varh1[6], varh1[7]));
TLatex* texth3 = new TLatex(); texth3->SetNDC(); texth3->SetTextSize(0.025); texth3->SetTextColor(kRed+2);
texth3->DrawLatex(0.61, 0.74, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth3->DrawLatex(0.61, 0.71, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh3[1], varh3[1+3], varh3[2], varh3[2+3]));
texth3->DrawLatex(0.61, 0.68, Form("#chi^{2}/NDF = %.1f / %.0f", varh3[6], varh3[7]));
text->DrawLatex(0.015, 0.97, "H pT > 130 GeV, MET > 170 GeV, 2 CSVL");
gPad->Print("plots/breg_HmassNormSB.png");
*/

//- FatH mass ------------------------------------------------------------------
//TCut fjcut = "Vtype==4 && hJet_pt[0]>60 && hJet_pt[1]>30 && METtype1corr.et>170 && H.pt>130 && hJet_csv_nominal[0]>0.244 && hJet_csv_nominal[1]>0.244 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && nfathFilterJets>0 && fathFilterJets_pt[0]>0 && fathFilterJets_pt[1]>0 && H.mass>0";
TCut fjcut = "Vtype==4 && fathFilterJets_pt[0]>60 && fathFilterJets_pt[1]>30 && abs(fathFilterJets-eta[0])<2.5 && abs(fathFilterJets-eta[1])<2.5 && METtype1corr.et>170 && FatH.pt>130 && fathFilterJets_csv[0]>0.244 && fathFilterJets_csv[1]>0.244 && abs(deltaPhi(FatH.phi,METtype1corr.phi))>2.0 && nfathFilterJets>0 && fathFilterJets_pt[0]>0 && fathFilterJets_pt[1]>0 && FatH.mass>0";
//TCut fjcut = "Vtype==4 && fathFilterJets_pt[0]>60 && fathFilterJets_pt[1]>30 && METtype1corr.et>170 && FatH.pt>130 && fathFilterJets_csv[0]>0.244 && fathFilterJets_csv[1]>0.244 && abs(deltaPhi(FatH.phi,METtype1corr.phi))>2.0 && nfathFilterJets>0 && fathFilterJets_pt[0]>0 && fathFilterJets_pt[1]>0 && FatH.mass>0 && nfathFilterJets>2 && fathFilterJets[2]_pt>15";

/*
TH1::SetDefaultSumw2();
double xmin=107., xmax=137.;
double xxmin=77., xxmax=107.;
tree->Draw("FatH.mass >> h1(140,60,200)",fjcut);
tree->Draw("FatHmassReg >> h2(140,60,200)",fjcut);
treeZZ->Draw("FatH.mass >> h3(140,60,200)",fjcut);
treeZZ->Draw("FatHmassReg >> h4(140,60,200)",fjcut);
h1=h1;h2=h2;h3=h3;h4=h4;
h4->Scale(1.0/h4->Integral());
h3->Scale(1.0/h3->Integral());
h2->Scale(1.0/h2->Integral());
h1->Scale(1.0/h1->Integral());

double varh1[8], varh2[8], varh3[8], varh4[8];
h2->Fit("gaus", "I", "", xmin, xmax);
varh2[0] = gaus->GetParameter(0); varh2[1] = gaus->GetParameter(1); varh2[2] = gaus->GetParameter(2); varh2[3] = gaus->GetParError(0); varh2[4] = gaus->GetParError(1); varh2[5] = gaus->GetParError(2); varh2[6] = gaus->GetChisquare(); varh2[7] = gaus->GetNDF();
h4->Fit("gaus", "I", "", xxmin, xxmax);
varh4[0] = gaus->GetParameter(0); varh4[1] = gaus->GetParameter(1); varh4[2] = gaus->GetParameter(2); varh4[3] = gaus->GetParError(0); varh4[4] = gaus->GetParError(1); varh4[5] = gaus->GetParError(2); varh4[6] = gaus->GetChisquare(); varh4[7] = gaus->GetNDF();

h4->SetStats(0); h4->SetLineWidth(2); h4->SetMarkerSize(0); h4->SetLineColor(kGray+1);
h3->SetStats(0); h3->SetLineWidth(1); h3->SetMarkerSize(0); h3->SetLineColor(kBlack);
h2->SetStats(0); h2->SetLineWidth(2); h2->SetMarkerSize(0); h2->SetLineColor(kRed);
h1->SetStats(0); h1->SetLineWidth(1); h1->SetMarkerSize(0); h1->SetLineColor(kRed+2);
h3->SetTitle("; H(125)^{fj} mass [GeV]; (normalized)"); h3->Draw("hist");
h4->Draw("histsame");
h2->Draw("histsame");
h1->Draw("histsame");
h4->GetFunction("gaus")->Draw("same");
h2->GetFunction("gaus")->Draw("same");

TLegend* leg = new TLegend(0.62, 0.78, 0.92, 0.92);
leg->AddEntry(h4, "ZZ regressed");
leg->AddEntry(h3, "ZZ");
leg->AddEntry(h2, "ZH(125) regressed");
leg->AddEntry(h1, "ZH(125)");
leg->Draw();

double down = 0.12;
TLatex* texth2 = new TLatex(); texth2->SetNDC(); texth2->SetTextSize(0.025); texth2->SetTextColor(kRed);
texth2->DrawLatex(0.61, 0.74-down, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth2->DrawLatex(0.61, 0.71-down, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh2[1], varh2[1+3], varh2[2], varh2[2+3]));
texth2->DrawLatex(0.61, 0.68-down, Form("#chi^{2}/NDF = %.1f / %.0f", varh2[6], varh2[7]));
TLatex* texth4 = new TLatex(); texth4->SetNDC(); texth4->SetTextSize(0.025); texth4->SetTextColor(kGray+1);
texth4->DrawLatex(0.61, 0.74, Form("gaus fit [%.0f,%.0f]", xxmin, xxmax));
texth4->DrawLatex(0.61, 0.71, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh4[1], varh4[1+3], varh4[2], varh4[2+3]));
texth4->DrawLatex(0.61, 0.68, Form("#chi^{2}/NDF = %.1f / %.0f", varh4[6], varh4[7]));
text->DrawLatex(0.015, 0.97, "H pT > 130 GeV, MET > 170 GeV, 2 CSVL");
gPad->Print("plots/breg_FatHmassRegSB.png");

double xmin=100., xmax=135.;
double xxmin=70., xxmax=105.;
h1->Fit("gaus", "I", "", xmin, xmax); 
varh1[0] = gaus->GetParameter(0); varh1[1] = gaus->GetParameter(1); varh1[2] = gaus->GetParameter(2); varh1[3] = gaus->GetParError(0); varh1[4] = gaus->GetParError(1); varh1[5] = gaus->GetParError(2); varh1[6] = gaus->GetChisquare(); varh1[7] = gaus->GetNDF();
h3->Fit("gaus", "I", "", xxmin, xxmax);
varh3[0] = gaus->GetParameter(0); varh3[1] = gaus->GetParameter(1); varh3[2] = gaus->GetParameter(2); varh3[3] = gaus->GetParError(0); varh3[4] = gaus->GetParError(1); varh3[5] = gaus->GetParError(2); varh3[6] = gaus->GetChisquare(); varh3[7] = gaus->GetNDF();
h3->SetStats(0); h3->SetLineWidth(2); h3->SetMarkerSize(0); h3->SetLineColor(kBlack);
h1->SetStats(0); h1->SetLineWidth(2); h1->SetMarkerSize(0); h1->SetLineColor(kRed+2);
h3->SetTitle("; H(125)^{fj} mass [GeV]; (normalized)"); h3->Draw("hist");
h1->Draw("histsame");
h3->GetFunction("gaus")->Draw("same");
h1->GetFunction("gaus")->Draw("same");

TLegend* leg = new TLegend(0.62, 0.85, 0.92, 0.92);
leg->AddEntry(h3, "ZZ");
leg->AddEntry(h1, "ZH(125)");
leg->Draw();

TLatex* texth1 = new TLatex(); texth1->SetNDC(); texth1->SetTextSize(0.025); texth1->SetTextColor(kBlack);
texth1->DrawLatex(0.61, 0.74-down, Form("gaus fit [%.0f,%.0f]", xxmin, xxmax));
texth1->DrawLatex(0.61, 0.71-down, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh1[1], varh1[1+3], varh1[2], varh1[2+3]));
texth1->DrawLatex(0.61, 0.68-down, Form("#chi^{2}/NDF = %.1f / %.0f", varh1[6], varh1[7]));
TLatex* texth3 = new TLatex(); texth3->SetNDC(); texth3->SetTextSize(0.025); texth3->SetTextColor(kRed+2);
texth3->DrawLatex(0.61, 0.74, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth3->DrawLatex(0.61, 0.71, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh3[1], varh3[1+3], varh3[2], varh3[2+3]));
texth3->DrawLatex(0.61, 0.68, Form("#chi^{2}/NDF = %.1f / %.0f", varh3[6], varh3[7]));
text->DrawLatex(0.015, 0.97, "H pT > 130 GeV, MET > 170 GeV, 2 CSVL");
gPad->Print("plots/breg_FatHmassNormSB.png");
*/

//- H mass difference ----------------------------------------------------------
/*
tree->Draw("HmassReg - FatHmassReg >> h1(50.,-50.,50.)", "efflumi" * (cut+"nfathFilterJets>0 && fathFilterJets_pt[0]>0 && fathFilterJets_pt[1]>0"), "goff");
treeTT->Draw("HmassReg - FatHmassReg >> h2(50.,-50.,50.)", "efflumi" * (cut+"nfathFilterJets>0 && fathFilterJets_pt[0]>0 && fathFilterJets_pt[1]>0"), "goff");
chainVV->Draw("HmassReg - FatHmassReg >> h3(50.,-50.,50.)", "efflumi" * (cut+"nfathFilterJets>0 && fathFilterJets_pt[0]>0 && fathFilterJets_pt[1]>0"), "goff");
//chainZJ->Draw("HmassReg - FatHmassReg >> h4(50.,-50.,50.)", "efflumi" * (cut+"nfathFilterJets>0 && fathFilterJets_pt[0]>0 && fathFilterJets_pt[1]>0 && eventFlav!=5"), "goff");
chainZJ->Draw("HmassReg - FatHmassReg >> h5(50.,-50.,50.)", "efflumi" * (cut+"nfathFilterJets>0 && fathFilterJets_pt[0]>0 && fathFilterJets_pt[1]>0 && eventFlav!=5"), "goff");

//h1=h1;h2=h2;h3=h3;h4=h4;h5=h5;
h1=h1;h2=h2;h3=h3;h5=h5;
h5->Scale(1.0/h5->Integral());
//h4->Scale(1.0/h4->Integral());
h3->Scale(1.0/h3->Integral());
h2->Scale(1.0/h2->Integral());
h1->Scale(1.0/h1->Integral());

h1->SetStats(0); h1->SetLineWidth(2); h1->SetMarkerSize(0); h1->SetLineColor(kRed);
h2->SetStats(0); h2->SetLineWidth(2); h2->SetMarkerSize(0); h2->SetLineColor(kBlue);
h3->SetStats(0); h3->SetLineWidth(2); h3->SetMarkerSize(0); h3->SetLineColor(kGray);
//h4->SetStats(0); h4->SetLineWidth(2); h4->SetMarkerSize(0); h4->SetLineColor(kYellow-6);
h5->SetStats(0); h5->SetLineWidth(2); h5->SetMarkerSize(0); h5->SetLineColor(kYellow);
h5->SetTitle("; H mass - H^{fj} mass [GeV]; (normalized)"); h5->SetMaximum(0.2); h5->Draw("hist");
//h4->Draw("histsame");
h3->Draw("histsame");
h2->Draw("histsame");
h1->Draw("histsame");

TLegend* leg = new TLegend(0.62, 0.75, 0.92, 0.92);
leg->AddEntry(h1, "ZH(125)");
leg->AddEntry(h2, "t #bar{t}");
leg->AddEntry(h3, "VV");
//leg->AddEntry(h4, "Z+udscg");
leg->AddEntry(h5, "Z+b#bar{b}");
leg->Draw();

text->DrawLatex(0.015, 0.97, "H pT > 130 GeV, MET > 170 GeV, 2 CSVL, regressed");
gPad->Print("plots/breg_HmassDiffSB.png");
*/

//- TopHad ---------------------------------------------------------------------

TCut tthighpt = "(Vtype==2||Vtype==3) && METtype1corr.et>170 && H.pt>130 && hJet_pt[0]>60 && hJet_pt[1]>30 &&  max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.898 && min(hJet_csv_nominal[0],hJet_csv_nominal[1])<0.244 && mindPhiMETJet_dPhi>0.5 && Max$(aJet_csv_nominal)>0.679 && naJets_Znn>=2 && naJets_Znn<=4 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && H.mass>0"; // FIXME: Vtype==4?
treeTT->Draw("TopHad_mass >> h1(40, 100., 300.)", tthighpt + "TopHad_pt>170.", "goff");
treeTT->Draw("TopHad_massReg >> h2(40, 100., 300.)", tthighpt + "TopHad_ptReg>170.", "goff");

TH1::SetDefaultSumw2();
double xmin=165., xmax=195.;
double varh1[8], varh2[8], varh3[8], varh4[8];
h1->Fit("gaus", "I", "", xmin, xmax); 
varh1[0] = gaus->GetParameter(0); varh1[1] = gaus->GetParameter(1); varh1[2] = gaus->GetParameter(2); varh1[3] = gaus->GetParError(0); varh1[4] = gaus->GetParError(1); varh1[5] = gaus->GetParError(2); varh1[6] = gaus->GetChisquare(); varh1[7] = gaus->GetNDF();
h2->Fit("gaus", "I", "", xmin, xmax);
varh2[0] = gaus->GetParameter(0); varh2[1] = gaus->GetParameter(1); varh2[2] = gaus->GetParameter(2); varh2[3] = gaus->GetParError(0); varh2[4] = gaus->GetParError(1); varh2[5] = gaus->GetParError(2); varh2[6] = gaus->GetChisquare(); varh2[7] = gaus->GetNDF();

h2->SetStats(0); h2->SetLineWidth(2); h2->SetMarkerSize(0); h2->SetLineColor(kRed);
h1->SetStats(0); h1->SetLineWidth(1); h1->SetMarkerSize(0); h1->SetLineColor(kBlack);
h2->SetTitle("; t_{had} mass [GeV];"); h2->Draw("hist");
h1->Draw("histsame");
h2->GetFunction("gaus")->Draw("same");
h1->GetFunction("gaus")->Draw("same");

TLegend* leg = new TLegend(0.62, 0.78, 0.92, 0.92);
leg->AddEntry(h2, "t#bar{t} regressed");
leg->AddEntry(h1, "t#bar{t}");
leg->Draw();

double down = 0.12;
TLatex* texth1 = new TLatex(); texth1->SetNDC(); texth1->SetTextSize(0.025); texth1->SetTextColor(kBlack);
texth1->DrawLatex(0.61, 0.74-down, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth1->DrawLatex(0.61, 0.71-down, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh1[1], varh1[1+3], varh1[2], varh1[2+3]));
texth1->DrawLatex(0.61, 0.68-down, Form("#chi^{2}/NDF = %.1f / %.0f", varh1[6], varh1[7]));
TLatex* texth2 = new TLatex(); texth2->SetNDC(); texth2->SetTextSize(0.025); texth2->SetTextColor(kRed);
texth2->DrawLatex(0.61, 0.74, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth2->DrawLatex(0.61, 0.71, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh2[1], varh2[1+3], varh2[2], varh2[2+3]));
texth2->DrawLatex(0.61, 0.68, Form("#chi^{2}/NDF = %.1f / %.0f", varh2[6], varh2[7]));
text->DrawLatex(0.015, 0.97, "t_{had} pT > 170 GeV, MET > 170 GeV, 1 CSVT, 1 !CSVL");
gPad->Print("plots/breg_TopHad_massReg_TT.png");


//- FatTopHad ---------------------------------------------------------------------
/*
//TCut fjtthighpt = "(Vtype==2||Vtype==3) && METtype1corr.et>170 && fathFilterJets_pt[0]>60 && fathFilterJets_pt[1]>30 && abs(fathFilterJets_eta[0])<2.5 && abs(fathFilterJets_eta[1])<2.5 && mindPhiMETJet_dPhi>0.5 && max(fathFilterJets_csv[0],fathFilterJets_csv[1])>0.898 && min(fathFilterJets_csv[0],fathFilterJets_csv[1])<0.244 && (Max$(aJetFat_csv)>0.679) && ((Sum$(fathFilterJets_pt>30 && abs(fathFilterJets_eta)<2.5) + Sum$(aJetFat_pt>30 && abs(aJetFat_eta)<4.5))>=4) && ((Sum$(fathFilterJets_pt>30 && abs(fathFilterJets_eta)<2.5) + Sum$(aJetFat_pt>30 && abs(aJetFat_eta)<4.5))<=6) && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && nfathFilterJets>0 && fathFilterJets_pt[0]>0 && fathFilterJets_pt[1]>0 && FatH.mass>0";
TCut fjtthighpt = "(Vtype==2||Vtype==3) && METtype1corr.et>170 && FatH.pt>130 && fathFilterJets_pt[0]>60 && fathFilterJets_pt[1]>30 && abs(fathFilterJets_eta[0])<2.5 && abs(fathFilterJets_eta[1])<2.5 && mindPhiMETJet_dPhi>0.5 && max(fathFilterJets_csv[0],fathFilterJets_csv[1])>0.898 && min(fathFilterJets_csv[0],fathFilterJets_csv[1])<0.244 && (Max$(aJetFat_csv)>0.679) && naJets_Znn>=2 && naJets_Znn<=4 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && nfathFilterJets>0 && fathFilterJets_pt[0]>0 && fathFilterJets_pt[1]>0 && FatH.mass>0";  // FIXME: Vtype==4

treeTT->Draw("FatTopHad_mass >> h1(40, 100., 300.)", fjtthighpt + "FatTopHad_pt>170. && FatTopHad_j3Pt>20", "goff");
treeTT->Draw("FatTopHad_massReg >> h2(40, 100., 300.)", fjtthighpt + "FatTopHad_ptReg>170. && FatTopHad_j3Pt>20", "goff");

TH1::SetDefaultSumw2();
double xmin=150., xmax=190.;
double xxmin=145., xxmax=185.;
double varh1[8], varh2[8], varh3[8], varh4[8];
h1->Fit("gaus", "I", "", xxmin, xxmax); 
varh1[0] = gaus->GetParameter(0); varh1[1] = gaus->GetParameter(1); varh1[2] = gaus->GetParameter(2); varh1[3] = gaus->GetParError(0); varh1[4] = gaus->GetParError(1); varh1[5] = gaus->GetParError(2); varh1[6] = gaus->GetChisquare(); varh1[7] = gaus->GetNDF();
h2->Fit("gaus", "I", "", xmin, xmax);
varh2[0] = gaus->GetParameter(0); varh2[1] = gaus->GetParameter(1); varh2[2] = gaus->GetParameter(2); varh2[3] = gaus->GetParError(0); varh2[4] = gaus->GetParError(1); varh2[5] = gaus->GetParError(2); varh2[6] = gaus->GetChisquare(); varh2[7] = gaus->GetNDF();

h2=h2;
h2->SetStats(0); h2->SetLineWidth(2); h2->SetMarkerSize(0); h2->SetLineColor(kRed);
h1->SetStats(0); h1->SetLineWidth(1); h1->SetMarkerSize(0); h1->SetLineColor(kBlack);
h2->SetTitle("; t_{had}^{fj} mass [GeV];"); h2->SetMaximum(h2->GetMaximum()*1.1); h2->Draw("hist");
h1->Draw("histsame");
h2->GetFunction("gaus")->Draw("same");
h1->GetFunction("gaus")->Draw("same");

TLegend* leg = new TLegend(0.62, 0.78, 0.92, 0.92);
leg->AddEntry(h2, "t#bar{t} regressed");
leg->AddEntry(h1, "t#bar{t}");
leg->Draw();

double down = 0.12;
TLatex* texth1 = new TLatex(); texth1->SetNDC(); texth1->SetTextSize(0.025); texth1->SetTextColor(kBlack);
texth1->DrawLatex(0.61, 0.74-down, Form("gaus fit [%.0f,%.0f]", xxmin, xxmax));
texth1->DrawLatex(0.61, 0.71-down, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh1[1], varh1[1+3], varh1[2], varh1[2+3]));
texth1->DrawLatex(0.61, 0.68-down, Form("#chi^{2}/NDF = %.1f / %.0f", varh1[6], varh1[7]));
TLatex* texth2 = new TLatex(); texth2->SetNDC(); texth2->SetTextSize(0.025); texth2->SetTextColor(kRed);
texth2->DrawLatex(0.61, 0.74, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth2->DrawLatex(0.61, 0.71, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh2[1], varh2[1+3], varh2[2], varh2[2+3]));
texth2->DrawLatex(0.61, 0.68, Form("#chi^{2}/NDF = %.1f / %.0f", varh2[6], varh2[7]));
text->DrawLatex(0.015, 0.97, "t_{had} pT > 170 GeV, MET > 170 GeV, 1 CSVT, 1 !CSVL");
gPad->Print("plots/breg_FatTopHad_massReg_TT.png");
*/

//- TopLep ---------------------------------------------------------------------
/*
TCut tthighpt = "(Vtype==2||Vtype==2) && METtype1corr.et>170 && H.pt>130 && hJet_pt[0]>60 && hJet_pt[1]>30 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.898 && min(hJet_csv_nominal[0],hJet_csv_nominal[1])<0.244 && mindPhiMETJet_dPhi>0.5 && V.mass<170 && naJets_Znn==1 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && H.mass>0";

//treeTT->Draw("TopLep_mass >> h1(40, 95., 295.)", tthighpt + "TopLep_pt>170.", "goff");
//treeTT->Draw("TopLep_massReg >> h2(40, 95., 295.)", tthighpt + "TopLep_ptReg>170.", "goff");
chainST->Draw("TopLep_mass >> h1(40, 95., 295.)", "efflumi"*(tthighpt + "TopLep_pt>150."), "goff");
chainST->Draw("TopLep_massReg >> h2(40, 95., 295.)", "efflumi"*(tthighpt + "TopLep_pt>150."), "goff");

TH1::SetDefaultSumw2();
double xmin=155., xmax=190.;
double varh1[8], varh2[8], varh3[8], varh4[8];
h1->Fit("gaus", "I", "", xmin, xmax); 
varh1[0] = gaus->GetParameter(0); varh1[1] = gaus->GetParameter(1); varh1[2] = gaus->GetParameter(2); varh1[3] = gaus->GetParError(0); varh1[4] = gaus->GetParError(1); varh1[5] = gaus->GetParError(2); varh1[6] = gaus->GetChisquare(); varh1[7] = gaus->GetNDF();
h2->Fit("gaus", "I", "", xmin, xmax);
varh2[0] = gaus->GetParameter(0); varh2[1] = gaus->GetParameter(1); varh2[2] = gaus->GetParameter(2); varh2[3] = gaus->GetParError(0); varh2[4] = gaus->GetParError(1); varh2[5] = gaus->GetParError(2); varh2[6] = gaus->GetChisquare(); varh2[7] = gaus->GetNDF();

h2->SetStats(0); h2->SetLineWidth(2); h2->SetMarkerSize(0); h2->SetLineColor(kRed);
h1->SetStats(0); h1->SetLineWidth(1); h1->SetMarkerSize(0); h1->SetLineColor(kBlack);
h2->SetTitle("; t_{lep} mass [GeV];"); h2->Draw("hist");
h1->Draw("histsame");
h2->GetFunction("gaus")->Draw("same");
h1->GetFunction("gaus")->Draw("same");

TLegend* leg = new TLegend(0.62, 0.78, 0.92, 0.92);
leg->AddEntry(h2, "t#bar{t} regressed");
leg->AddEntry(h1, "t#bar{t}");
leg->Draw();

double down = 0.12;
TLatex* texth1 = new TLatex(); texth1->SetNDC(); texth1->SetTextSize(0.025); texth1->SetTextColor(kBlack);
texth1->DrawLatex(0.61, 0.74-down, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth1->DrawLatex(0.61, 0.71-down, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh1[1], varh1[1+3], varh1[2], varh1[2+3]));
texth1->DrawLatex(0.61, 0.68-down, Form("#chi^{2}/NDF = %.1f / %.0f", varh1[6], varh1[7]));
TLatex* texth2 = new TLatex(); texth2->SetNDC(); texth2->SetTextSize(0.025); texth2->SetTextColor(kRed);
texth2->DrawLatex(0.61, 0.74, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth2->DrawLatex(0.61, 0.71, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh2[1], varh2[1+3], varh2[2], varh2[2+3]));
texth2->DrawLatex(0.61, 0.68, Form("#chi^{2}/NDF = %.1f / %.0f", varh2[6], varh2[7]));
text->DrawLatex(0.015, 0.97, "t_{lep} pT > 170 GeV, MET > 170 GeV, 1 CSVT, 1 !CSVL");
//gPad->Print("plots/breg_TopLep_massReg_TT.png");
gPad->Print("plots/breg_TopLep_massReg_ST.png");
*/

//------------------------------------------------------------------------------
/*
TCut cutData = "triggerFlags[42]||triggerFlags[39]||triggerFlags[49]";

treeData->Draw("TopLep_mass >> h1(20, 95., 295.)", tthighpt + cutData + "TopLep_pt>170.", "goff");
treeData->Draw("TopLep_massReg >> h2(20, 95., 295.)", tthighpt + cutData + "TopLep_ptReg>170.", "goff");

TH1::SetDefaultSumw2();
double xmin=150., xmax=200.;
double varh1[8], varh2[8], varh3[8], varh4[8];
h1->Fit("gaus", "I", "", xmin, xmax); 
varh1[0] = gaus->GetParameter(0); varh1[1] = gaus->GetParameter(1); varh1[2] = gaus->GetParameter(2); varh1[3] = gaus->GetParError(0); varh1[4] = gaus->GetParError(1); varh1[5] = gaus->GetParError(2); varh1[6] = gaus->GetChisquare(); varh1[7] = gaus->GetNDF();
h2->Fit("gaus", "I", "", xmin, xmax);
varh2[0] = gaus->GetParameter(0); varh2[1] = gaus->GetParameter(1); varh2[2] = gaus->GetParameter(2); varh2[3] = gaus->GetParError(0); varh2[4] = gaus->GetParError(1); varh2[5] = gaus->GetParError(2); varh2[6] = gaus->GetChisquare(); varh2[7] = gaus->GetNDF();

h2=h2;
h2->SetStats(0); h2->SetLineWidth(2); h2->SetMarkerSize(0); h2->SetLineColor(kRed);
h1->SetStats(0); h1->SetLineWidth(1); h1->SetMarkerSize(0); h1->SetLineColor(kBlack);
h2->SetTitle("; t_{lep} mass [GeV];"); h2->SetMaximum(h2->GetMaximum()*1.1); h2->Draw("hist");
h1->Draw("histsame");
h2->GetFunction("gaus")->Draw("same");
h1->GetFunction("gaus")->Draw("same");

TLegend* leg = new TLegend(0.62, 0.78, 0.92, 0.92);
leg->AddEntry(h2, "t#bar{t} regressed");
leg->AddEntry(h1, "t#bar{t}");
leg->Draw();

double down = 0.12;
TLatex* texth1 = new TLatex(); texth1->SetNDC(); texth1->SetTextSize(0.025); texth1->SetTextColor(kBlack);
texth1->DrawLatex(0.61, 0.74-down, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth1->DrawLatex(0.61, 0.71-down, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh1[1], varh1[1+3], varh1[2], varh1[2+3]));
texth1->DrawLatex(0.61, 0.68-down, Form("#chi^{2}/NDF = %.1f / %.0f", varh1[6], varh1[7]));
TLatex* texth2 = new TLatex(); texth2->SetNDC(); texth2->SetTextSize(0.025); texth2->SetTextColor(kRed);
texth2->DrawLatex(0.61, 0.74, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth2->DrawLatex(0.61, 0.71, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh2[1], varh2[1+3], varh2[2], varh2[2+3]));
texth2->DrawLatex(0.61, 0.68, Form("#chi^{2}/NDF = %.1f / %.0f", varh2[6], varh2[7]));
text->DrawLatex(0.015, 0.97, "t_{lep} pT > 170 GeV, MET > 170 GeV, 1 CSVT, 1 !CSVL");
gPad->Print("plots/breg_TopLep_massReg_Data.png");
*/

/*
gROOT->LoadMacro("HelperFunctions.h");
treeData->Draw("TopLep_pt >> h1(100, 0., 500.)", tthighpt + cutData + "TopLep_pt>100.", "goff");
treeTT->Draw("TopLep_pt >> h2(100, 0., 500.)", "efflumi * PUweight * triggerweight2012ABC(METtype1corr.et)"*(tthighpt + "TopLep_pt>100."), "goff");
chainST->Draw("TopLep_pt >> h3(100, 0., 500.)", "efflumi * PUweight * triggerweight2012ABC(METtype1corr.et)"*(tthighpt + "TopLep_pt>100."), "goff");
h3=h3;
//h2->Add(h3);
h2->SetLineColor(kCyan);
h1->Draw();
h2->Draw("same");
*/

}
