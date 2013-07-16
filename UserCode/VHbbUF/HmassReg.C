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

gROOT->SetBatch();

// Get tree
//TFile* _file0 = TFile::Open("skim_ZnnH_baseline/Step3_ZnnH125.root");
TFile* _file0 = TFile::Open("Step3_20130204/Step3_ZnnH125.root");
TTree* tree = (TTree*) _file0->Get("tree");


//- H mass ---------------------------------------------------------------------
//TCut cut = "(Vtype==4||Vtype==3||Vtype==2) && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5 && hJet_id[0] && hJet_id[1] && hJet_pt[0]>20 && hJet_pt[1]>20 && hJet_csv[0]>0 && hJet_csv[1]>0 && H.pt>80 && METtype1corr.et>100";
TCut cut = "Vtype==4 && hJet_pt[0]>60 && hJet_pt[1]>30 && METtype1corr.et>170 && H.pt>130 && hJet_csv_nominal[0]>0.244 && hJet_csv_nominal[1]>0.244 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && H.mass>0";
//TCut cut = "Vtype==4 && hJet_ptReg[0]>80 && hJet_ptReg[1]>30 && HmassReg<250 && HptReg>160 && METtype1corr.et>160 && max(hJet_csvCor[0],hJet_csvCor[1])>0.679 && min(hJet_csvCor[0],hJet_csvCor[1])>0.244 && HVdPhi>2.0 && min(Min$(abs(deltaPhi(METtype1corr.phi,hJet_phi))),Min$(abs(deltaPhiMETjets(METtype1corr.phi,aJet_phi,aJet_pt,aJet_eta)))+999*(Sum$(aJet_pt>30 && abs(aJet_eta)<2.5)==0) )>0.5 && abs(deltaPhi(METtype1corr.phi,METnoPUCh.phi))<0.5";

/*
double xmin=110., xmax=140.;
tree->Draw("H.mass >> h1(140,60,200)",cut);
tree->Draw("HmassReg >> h2(140,60,200)",cut);
tree->Draw("HmassGen >> h3(140,60,200)",cut);
double varh1[8], varh2[8], varh3[8];
h1->Fit("gaus", "I", "", xmin, xmax); 
varh1[0] = gaus->GetParameter(0); varh1[1] = gaus->GetParameter(1); varh1[2] = gaus->GetParameter(2); varh1[3] = gaus->GetParError(0); varh1[4] = gaus->GetParError(1); varh1[5] = gaus->GetParError(2); varh1[6] = gaus->GetChisquare(); varh1[7] = gaus->GetNDF();
h2->Fit("gaus", "I", "", xmin, xmax); double h2mean = gaus->GetParameter(1); double h2sigma = gaus->GetParameter(2);
varh2[0] = gaus->GetParameter(0); varh2[1] = gaus->GetParameter(1); varh2[2] = gaus->GetParameter(2); varh2[3] = gaus->GetParError(0); varh2[4] = gaus->GetParError(1); varh2[5] = gaus->GetParError(2); varh2[6] = gaus->GetChisquare(); varh2[7] = gaus->GetNDF();
h3->Fit("gaus", "I", "", xmin, xmax); double h3mean = gaus->GetParameter(1); double h3sigma = gaus->GetParameter(2);
varh3[0] = gaus->GetParameter(0); varh3[1] = gaus->GetParameter(1); varh3[2] = gaus->GetParameter(2); varh3[3] = gaus->GetParError(0); varh3[4] = gaus->GetParError(1); varh3[5] = gaus->GetParError(2); varh3[6] = gaus->GetChisquare(); varh3[7] = gaus->GetNDF();
h3->SetStats(0); h3->SetLineWidth(2); h3->SetLineColor(kBlue); h3->SetMarkerSize(0);
h3->SetTitle("; H(125) mass [GeV]"); h3->Draw("");
h2->SetStats(0); h2->SetLineWidth(2); h2->SetLineColor(kRed); h2->SetMarkerSize(0);
h2->Draw("same");
h1->SetStats(0); h1->SetLineWidth(2); h1->SetMarkerSize(0);
h1->Draw("same");

TLegend* leg = new TLegend(0.62, 0.78, 0.92, 0.92);
leg->AddEntry(h3, "gen");
leg->AddEntry(h2, "CMS regressed");
leg->AddEntry(h1, "CMS");
leg->Draw();

double down = 0.12;
TLatex* texth1 = new TLatex(); texth1->SetNDC(); texth1->SetTextSize(0.025);
texth1->DrawLatex(0.21, 0.90-down*2, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth1->DrawLatex(0.21, 0.87-down*2, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh1[1], varh1[1+3], varh1[2], varh1[2+3]));
texth1->DrawLatex(0.21, 0.84-down*2, Form("#chi^{2}/NDF = %.1f / %.0f", varh1[6], varh1[7]));
TLatex* texth2 = new TLatex(); texth2->SetNDC(); texth2->SetTextSize(0.025); texth2->SetTextColor(kRed);
texth2->DrawLatex(0.21, 0.90-down, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth2->DrawLatex(0.21, 0.87-down, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh2[1], varh2[1+3], varh2[2], varh2[2+3]));
texth2->DrawLatex(0.21, 0.84-down, Form("#chi^{2}/NDF = %.1f / %.0f", varh2[6], varh2[7]));
TLatex* texth3 = new TLatex(); texth3->SetNDC(); texth3->SetTextSize(0.025); texth3->SetTextColor(kBlue);
texth3->DrawLatex(0.21, 0.90, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth3->DrawLatex(0.21, 0.87, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh3[1], varh3[1+3], varh3[2], varh3[2+3]));
texth3->DrawLatex(0.21, 0.84, Form("#chi^{2}/NDF = %.1f / %.0f", varh3[6], varh3[7]));
text->DrawLatex(0.015, 0.97, "H pT > 130 GeV, MET > 170 GeV, 2 CSVL");
gPad->Print("plots/breg_HmassReg.png");

// repeat for low pT
TCut lowcut = "Vtype==4 && hJet_pt[0]>60 && hJet_pt[1]>30 && 130<METtype1corr.et && METtype1corr.et<170 && H.pt>130 && hJet_csv_nominal[0]>0.244 && hJet_csv_nominal[1]>0.244 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && H.mass>0";
//TCut lowcut = "Vtype==4 && hJet_pt[0]>60 && hJet_pt[1]>30 && 120<METtype1corr.et && METtype1corr.et<160 && H.pt>120 && hJet_csv[0]>0.244 && hJet_csv[1]>0.244"; // && Sum$(aJet_pt>30)==0";

double xmin=105., xmax=140.;
tree->Draw("H.mass >> h1(140,60,200)",lowcut);
tree->Draw("HmassReg >> h2(140,60,200)",lowcut);
tree->Draw("HmassGen >> h3(140,60,200)",lowcut);
double varh1[8], varh2[8], varh3[8];
h1->Fit("gaus", "I", "", xmin, xmax); 
varh1[0] = gaus->GetParameter(0); varh1[1] = gaus->GetParameter(1); varh1[2] = gaus->GetParameter(2); varh1[3] = gaus->GetParError(0); varh1[4] = gaus->GetParError(1); varh1[5] = gaus->GetParError(2); varh1[6] = gaus->GetChisquare(); varh1[7] = gaus->GetNDF();
h2->Fit("gaus", "I", "", xmin, xmax); double h2mean = gaus->GetParameter(1); double h2sigma = gaus->GetParameter(2);
varh2[0] = gaus->GetParameter(0); varh2[1] = gaus->GetParameter(1); varh2[2] = gaus->GetParameter(2); varh2[3] = gaus->GetParError(0); varh2[4] = gaus->GetParError(1); varh2[5] = gaus->GetParError(2); varh2[6] = gaus->GetChisquare(); varh2[7] = gaus->GetNDF();
h3->Fit("gaus", "I", "", xmin, xmax); double h3mean = gaus->GetParameter(1); double h3sigma = gaus->GetParameter(2);
varh3[0] = gaus->GetParameter(0); varh3[1] = gaus->GetParameter(1); varh3[2] = gaus->GetParameter(2); varh3[3] = gaus->GetParError(0); varh3[4] = gaus->GetParError(1); varh3[5] = gaus->GetParError(2); varh3[6] = gaus->GetChisquare(); varh3[7] = gaus->GetNDF();
h3->SetStats(0); h3->SetLineWidth(2); h3->SetLineColor(kBlue); h3->SetMarkerSize(0);
h3->SetTitle("; H(125) mass [GeV]"); h3->Draw("");
h2->SetStats(0); h2->SetLineWidth(2); h2->SetLineColor(kRed); h2->SetMarkerSize(0);
h2->Draw("same");
h1->SetStats(0); h1->SetLineWidth(2); h1->SetMarkerSize(0);
h1->Draw("same");

TLegend* leg = new TLegend(0.62, 0.78, 0.92, 0.92);
leg->AddEntry(h3, "gen");
leg->AddEntry(h2, "CMS regressed");
leg->AddEntry(h1, "CMS");
leg->Draw();

double down = 0.12;
TLatex* texth1 = new TLatex(); texth1->SetNDC(); texth1->SetTextSize(0.025);
texth1->DrawLatex(0.21, 0.90-down*2, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth1->DrawLatex(0.21, 0.87-down*2, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh1[1], varh1[1+3], varh1[2], varh1[2+3]));
texth1->DrawLatex(0.21, 0.84-down*2, Form("#chi^{2}/NDF = %.1f / %.0f", varh1[6], varh1[7]));
TLatex* texth2 = new TLatex(); texth2->SetNDC(); texth2->SetTextSize(0.025); texth2->SetTextColor(kRed);
texth2->DrawLatex(0.21, 0.90-down, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth2->DrawLatex(0.21, 0.87-down, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh2[1], varh2[1+3], varh2[2], varh2[2+3]));
texth2->DrawLatex(0.21, 0.84-down, Form("#chi^{2}/NDF = %.1f / %.0f", varh2[6], varh2[7]));
TLatex* texth3 = new TLatex(); texth3->SetNDC(); texth3->SetTextSize(0.025); texth3->SetTextColor(kBlue);
texth3->DrawLatex(0.21, 0.90, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth3->DrawLatex(0.21, 0.87, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh3[1], varh3[1+3], varh3[2], varh3[2+3]));
texth3->DrawLatex(0.21, 0.84, Form("#chi^{2}/NDF = %.1f / %.0f", varh3[6], varh3[7]));
text->DrawLatex(0.015, 0.97, "H pT > 130 GeV, 130 < MET < 170 GeV, 2 CSVL");
gPad->Print("plots/breg_HmassReg_lowmet.png");

// repeat for highest pT
double xmin=110., xmax=140.;
tree->Draw("H.mass >> h1(140,60,200)",cut+"METtype1corr.et>200 && HptReg>200");
tree->Draw("HmassReg >> h2(140,60,200)",cut+"METtype1corr.et>200 && HptReg>200");
tree->Draw("HmassGen >> h3(140,60,200)",cut+"METtype1corr.et>200 && HptReg>200");
double varh1[8], varh2[8], varh3[8];
h1->Fit("gaus", "I", "", xmin, xmax); 
varh1[0] = gaus->GetParameter(0); varh1[1] = gaus->GetParameter(1); varh1[2] = gaus->GetParameter(2); varh1[3] = gaus->GetParError(0); varh1[4] = gaus->GetParError(1); varh1[5] = gaus->GetParError(2); varh1[6] = gaus->GetChisquare(); varh1[7] = gaus->GetNDF();
h2->Fit("gaus", "I", "", xmin, xmax); double h2mean = gaus->GetParameter(1); double h2sigma = gaus->GetParameter(2);
varh2[0] = gaus->GetParameter(0); varh2[1] = gaus->GetParameter(1); varh2[2] = gaus->GetParameter(2); varh2[3] = gaus->GetParError(0); varh2[4] = gaus->GetParError(1); varh2[5] = gaus->GetParError(2); varh2[6] = gaus->GetChisquare(); varh2[7] = gaus->GetNDF();
h3->Fit("gaus", "I", "", xmin, xmax); double h3mean = gaus->GetParameter(1); double h3sigma = gaus->GetParameter(2);
varh3[0] = gaus->GetParameter(0); varh3[1] = gaus->GetParameter(1); varh3[2] = gaus->GetParameter(2); varh3[3] = gaus->GetParError(0); varh3[4] = gaus->GetParError(1); varh3[5] = gaus->GetParError(2); varh3[6] = gaus->GetChisquare(); varh3[7] = gaus->GetNDF();
h3->SetStats(0); h3->SetLineWidth(2); h3->SetLineColor(kBlue); h3->SetMarkerSize(0);
h3->SetTitle("; H(125) mass [GeV]"); h3->Draw("");
h2->SetStats(0); h2->SetLineWidth(2); h2->SetLineColor(kRed); h2->SetMarkerSize(0);
h2->Draw("same");
h1->SetStats(0); h1->SetLineWidth(2); h1->SetMarkerSize(0);
h1->Draw("same");

TLegend* leg = new TLegend(0.62, 0.78, 0.92, 0.92);
leg->AddEntry(h3, "gen");
leg->AddEntry(h2, "CMS regressed");
leg->AddEntry(h1, "CMS");
leg->Draw();

double down = 0.12;
TLatex* texth1 = new TLatex(); texth1->SetNDC(); texth1->SetTextSize(0.025);
texth1->DrawLatex(0.21, 0.90-down*2, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth1->DrawLatex(0.21, 0.87-down*2, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh1[1], varh1[1+3], varh1[2], varh1[2+3]));
texth1->DrawLatex(0.21, 0.84-down*2, Form("#chi^{2}/NDF = %.1f / %.0f", varh1[6], varh1[7]));
TLatex* texth2 = new TLatex(); texth2->SetNDC(); texth2->SetTextSize(0.025); texth2->SetTextColor(kRed);
texth2->DrawLatex(0.21, 0.90-down, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth2->DrawLatex(0.21, 0.87-down, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh2[1], varh2[1+3], varh2[2], varh2[2+3]));
texth2->DrawLatex(0.21, 0.84-down, Form("#chi^{2}/NDF = %.1f / %.0f", varh2[6], varh2[7]));
TLatex* texth3 = new TLatex(); texth3->SetNDC(); texth3->SetTextSize(0.025); texth3->SetTextColor(kBlue);
texth3->DrawLatex(0.21, 0.90, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth3->DrawLatex(0.21, 0.87, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh3[1], varh3[1+3], varh3[2], varh3[2+3]));
texth3->DrawLatex(0.21, 0.84, Form("#chi^{2}/NDF = %.1f / %.0f", varh3[6], varh3[7]));
text->DrawLatex(0.015, 0.97, "H pT > 200 GeV, MET > 200 GeV, 2 CSVL");
gPad->Print("plots/breg_HmassReg_met200.png");
*/

//- FatH mass ------------------------------------------------------------------
//TCut fjcut = "Vtype==4 && hJet_pt[0]>60 && hJet_pt[1]>30 && METtype1corr.et>170 && H.pt>130 && hJet_csv_nominal[0]>0.244 && hJet_csv_nominal[1]>0.244 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && nfathFilterJets>0 && fathFilterJets_pt[0]>0 && fathFilterJets_pt[1]>0 && H.mass>0";
TCut fjcut = "Vtype==4 && fathFilterJets_pt[0]>60 && fathFilterJets_pt[1]>30 && METtype1corr.et>170 && FatH.pt>130 && fathFilterJets_csv[0]>0.244 && fathFilterJets_csv[1]>0.244 && abs(deltaPhi(FatH.phi,METtype1corr.phi))>2.0 && nfathFilterJets>0 && fathFilterJets_pt[0]>0 && fathFilterJets_pt[1]>0 && FatH.mass>0";
//TCut fjcut = "Vtype==4 && fathFilterJets_pt[0]>60 && fathFilterJets_pt[1]>30 && METtype1corr.et>170 && FatH.pt>130 && fathFilterJets_csv[0]>0.244 && fathFilterJets_csv[1]>0.244 && abs(deltaPhi(FatH.phi,METtype1corr.phi))>2.0 && nfathFilterJets>0 && fathFilterJets_pt[0]>0 && fathFilterJets_pt[1]>0 && FatH.mass>0 && nfathFilterJets>2 && fathFilterJets[2]_pt>15";

/*
double xmin=107., xmax=137.;
tree->Draw("FatH.mass >> h1(140,60,200)",fjcut);
tree->Draw("FatHmassReg >> h2(140,60,200)",fjcut);
tree->Draw("FatHmassGen >> h3(140,60,200)",fjcut);
double varh1[8], varh2[8], varh3[8];
h1->Fit("gaus", "I", "", 100., 130.); 
varh1[0] = gaus->GetParameter(0); varh1[1] = gaus->GetParameter(1); varh1[2] = gaus->GetParameter(2); varh1[3] = gaus->GetParError(0); varh1[4] = gaus->GetParError(1); varh1[5] = gaus->GetParError(2); varh1[6] = gaus->GetChisquare(); varh1[7] = gaus->GetNDF();
h2->Fit("gaus", "I", "", xmin, xmax); double h2mean = gaus->GetParameter(1); double h2sigma = gaus->GetParameter(2);
varh2[0] = gaus->GetParameter(0); varh2[1] = gaus->GetParameter(1); varh2[2] = gaus->GetParameter(2); varh2[3] = gaus->GetParError(0); varh2[4] = gaus->GetParError(1); varh2[5] = gaus->GetParError(2); varh2[6] = gaus->GetChisquare(); varh2[7] = gaus->GetNDF();
h3->Fit("gaus", "I", "", xmin, xmax); double h3mean = gaus->GetParameter(1); double h3sigma = gaus->GetParameter(2);
varh3[0] = gaus->GetParameter(0); varh3[1] = gaus->GetParameter(1); varh3[2] = gaus->GetParameter(2); varh3[3] = gaus->GetParError(0); varh3[4] = gaus->GetParError(1); varh3[5] = gaus->GetParError(2); varh3[6] = gaus->GetChisquare(); varh3[7] = gaus->GetNDF();
h3->SetStats(0); h3->SetLineWidth(2); h3->SetLineColor(kBlue); h3->SetMarkerSize(0);
h3->SetTitle("; H(125)^{fj} mass [GeV]"); h3->Draw("");
h2->SetStats(0); h2->SetLineWidth(2); h2->SetLineColor(kRed); h2->SetMarkerSize(0);
h2->Draw("same");
h1->SetStats(0); h1->SetLineWidth(2); h1->SetMarkerSize(0);
h1->Draw("same");

TLegend* leg = new TLegend(0.62, 0.78, 0.92, 0.92);
leg->AddEntry(h3, "gen");
leg->AddEntry(h2, "CMS regressed");
leg->AddEntry(h1, "CMS");
leg->Draw();

double down = 0.12;
TLatex* texth1 = new TLatex(); texth1->SetNDC(); texth1->SetTextSize(0.025);
texth1->DrawLatex(0.21, 0.90-down*2, Form("gaus fit [%.0f,%.0f]", 100., 130.));
texth1->DrawLatex(0.21, 0.87-down*2, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh1[1], varh1[1+3], varh1[2], varh1[2+3]));
texth1->DrawLatex(0.21, 0.84-down*2, Form("#chi^{2}/NDF = %.1f / %.0f", varh1[6], varh1[7]));
TLatex* texth2 = new TLatex(); texth2->SetNDC(); texth2->SetTextSize(0.025); texth2->SetTextColor(kRed);
texth2->DrawLatex(0.21, 0.90-down, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth2->DrawLatex(0.21, 0.87-down, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh2[1], varh2[1+3], varh2[2], varh2[2+3]));
texth2->DrawLatex(0.21, 0.84-down, Form("#chi^{2}/NDF = %.1f / %.0f", varh2[6], varh2[7]));
TLatex* texth3 = new TLatex(); texth3->SetNDC(); texth3->SetTextSize(0.025); texth3->SetTextColor(kBlue);
texth3->DrawLatex(0.21, 0.90, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth3->DrawLatex(0.21, 0.87, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh3[1], varh3[1+3], varh3[2], varh3[2+3]));
texth3->DrawLatex(0.21, 0.84, Form("#chi^{2}/NDF = %.1f / %.0f", varh3[6], varh3[7]));
text->DrawLatex(0.015, 0.97, "H pT > 130 GeV, MET > 170 GeV, 2 CSVL");
gPad->Print("plots/breg_FatHmassReg.png");

// repeat for low pT
//TCut fjlowcut = "Vtype==4 && hJet_pt[0]>60 && hJet_pt[1]>30 && 130<METtype1corr.et && METtype1corr.et<170 && H.pt>130 && hJet_csv_nominal[0]>0.244 && hJet_csv_nominal[1]>0.244 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && nfathFilterJets>0 && fathFilterJets_pt[0]>0 && fathFilterJets_pt[1]>0 && H.mass>0";
TCut fjlowcut = "Vtype==4 && fathFilterJets_pt[0]>60 && fathFilterJets_pt[1]>30 && 130<METtype1corr.et && METtype1corr.et<170 && FatH.pt>130 && fathFilterJets_csv[0]>0.244 && fathFilterJets_csv[1]>0.244 && abs(deltaPhi(FatH.phi,METtype1corr.phi))>2.0 && nfathFilterJets>0 && fathFilterJets_pt[0]>0 && fathFilterJets_pt[1]>0 && FatH.mass>0";

double xmin=102., xmax=137.;
tree->Draw("H.mass >> h1(140,60,200)",fjlowcut);
tree->Draw("HmassReg >> h2(140,60,200)",fjlowcut);
tree->Draw("HmassGen >> h3(140,60,200)",fjlowcut);
double varh1[8], varh2[8], varh3[8];
h1->Fit("gaus", "I", "", 95., 130.); 
varh1[0] = gaus->GetParameter(0); varh1[1] = gaus->GetParameter(1); varh1[2] = gaus->GetParameter(2); varh1[3] = gaus->GetParError(0); varh1[4] = gaus->GetParError(1); varh1[5] = gaus->GetParError(2); varh1[6] = gaus->GetChisquare(); varh1[7] = gaus->GetNDF();
h2->Fit("gaus", "I", "", xmin, xmax); double h2mean = gaus->GetParameter(1); double h2sigma = gaus->GetParameter(2);
varh2[0] = gaus->GetParameter(0); varh2[1] = gaus->GetParameter(1); varh2[2] = gaus->GetParameter(2); varh2[3] = gaus->GetParError(0); varh2[4] = gaus->GetParError(1); varh2[5] = gaus->GetParError(2); varh2[6] = gaus->GetChisquare(); varh2[7] = gaus->GetNDF();
h3->Fit("gaus", "I", "", xmin, xmax); double h3mean = gaus->GetParameter(1); double h3sigma = gaus->GetParameter(2);
varh3[0] = gaus->GetParameter(0); varh3[1] = gaus->GetParameter(1); varh3[2] = gaus->GetParameter(2); varh3[3] = gaus->GetParError(0); varh3[4] = gaus->GetParError(1); varh3[5] = gaus->GetParError(2); varh3[6] = gaus->GetChisquare(); varh3[7] = gaus->GetNDF();
h3->SetStats(0); h3->SetLineWidth(2); h3->SetLineColor(kBlue); h3->SetMarkerSize(0);
h3->SetTitle("; H(125)^{fj} mass [GeV]"); h3->Draw("");
h2->SetStats(0); h2->SetLineWidth(2); h2->SetLineColor(kRed); h2->SetMarkerSize(0);
h2->Draw("same");
h1->SetStats(0); h1->SetLineWidth(2); h1->SetMarkerSize(0);
h1->Draw("same");

TLegend* leg = new TLegend(0.62, 0.78, 0.92, 0.92);
leg->AddEntry(h3, "gen");
leg->AddEntry(h2, "CMS regressed");
leg->AddEntry(h1, "CMS");
leg->Draw();

double down = 0.12;
TLatex* texth1 = new TLatex(); texth1->SetNDC(); texth1->SetTextSize(0.025);
texth1->DrawLatex(0.21, 0.90-down*2, Form("gaus fit [%.0f,%.0f]", 95., 130.));
texth1->DrawLatex(0.21, 0.87-down*2, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh1[1], varh1[1+3], varh1[2], varh1[2+3]));
texth1->DrawLatex(0.21, 0.84-down*2, Form("#chi^{2}/NDF = %.1f / %.0f", varh1[6], varh1[7]));
TLatex* texth2 = new TLatex(); texth2->SetNDC(); texth2->SetTextSize(0.025); texth2->SetTextColor(kRed);
texth2->DrawLatex(0.21, 0.90-down, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth2->DrawLatex(0.21, 0.87-down, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh2[1], varh2[1+3], varh2[2], varh2[2+3]));
texth2->DrawLatex(0.21, 0.84-down, Form("#chi^{2}/NDF = %.1f / %.0f", varh2[6], varh2[7]));
TLatex* texth3 = new TLatex(); texth3->SetNDC(); texth3->SetTextSize(0.025); texth3->SetTextColor(kBlue);
texth3->DrawLatex(0.21, 0.90, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth3->DrawLatex(0.21, 0.87, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh3[1], varh3[1+3], varh3[2], varh3[2+3]));
texth3->DrawLatex(0.21, 0.84, Form("#chi^{2}/NDF = %.1f / %.0f", varh3[6], varh3[7]));
text->DrawLatex(0.015, 0.97, "H pT > 130 GeV, 130 < MET < 170 GeV, 2 CSVL");
gPad->Print("plots/breg_FatHmassReg_lowmet.png");

// repeat for highest pt
double xmin=108., xmax=138.;
tree->Draw("FatH.mass >> h1(140,60,200)",fjcut+"METtype1corr.et>200 && FatHptReg>200");
tree->Draw("FatHmassReg >> h2(140,60,200)",fjcut+"METtype1corr.et>200 && FatHptReg>200");
tree->Draw("FatHmassGen >> h3(140,60,200)",fjcut+"METtype1corr.et>200 && FatHptReg>200");
double varh1[8], varh2[8], varh3[8];
h1->Fit("gaus", "I", "", 102., 132.); 
varh1[0] = gaus->GetParameter(0); varh1[1] = gaus->GetParameter(1); varh1[2] = gaus->GetParameter(2); varh1[3] = gaus->GetParError(0); varh1[4] = gaus->GetParError(1); varh1[5] = gaus->GetParError(2); varh1[6] = gaus->GetChisquare(); varh1[7] = gaus->GetNDF();
h2->Fit("gaus", "I", "", xmin, xmax); double h2mean = gaus->GetParameter(1); double h2sigma = gaus->GetParameter(2);
varh2[0] = gaus->GetParameter(0); varh2[1] = gaus->GetParameter(1); varh2[2] = gaus->GetParameter(2); varh2[3] = gaus->GetParError(0); varh2[4] = gaus->GetParError(1); varh2[5] = gaus->GetParError(2); varh2[6] = gaus->GetChisquare(); varh2[7] = gaus->GetNDF();
h3->Fit("gaus", "I", "", xmin, xmax); double h3mean = gaus->GetParameter(1); double h3sigma = gaus->GetParameter(2);
varh3[0] = gaus->GetParameter(0); varh3[1] = gaus->GetParameter(1); varh3[2] = gaus->GetParameter(2); varh3[3] = gaus->GetParError(0); varh3[4] = gaus->GetParError(1); varh3[5] = gaus->GetParError(2); varh3[6] = gaus->GetChisquare(); varh3[7] = gaus->GetNDF();
h3->SetStats(0); h3->SetLineWidth(2); h3->SetLineColor(kBlue); h3->SetMarkerSize(0);
h3->SetTitle("; H(125)^{fj} mass [GeV]"); h3->Draw("");
h2->SetStats(0); h2->SetLineWidth(2); h2->SetLineColor(kRed); h2->SetMarkerSize(0);
h2->Draw("same");
h1->SetStats(0); h1->SetLineWidth(2); h1->SetMarkerSize(0);
h1->Draw("same");

TLegend* leg = new TLegend(0.62, 0.78, 0.92, 0.92);
leg->AddEntry(h3, "gen");
leg->AddEntry(h2, "CMS regressed");
leg->AddEntry(h1, "CMS");
leg->Draw();

double down = 0.12;
TLatex* texth1 = new TLatex(); texth1->SetNDC(); texth1->SetTextSize(0.025);
texth1->DrawLatex(0.21, 0.90-down*2, Form("gaus fit [%.0f,%.0f]", 102., 132.));
texth1->DrawLatex(0.21, 0.87-down*2, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh1[1], varh1[1+3], varh1[2], varh1[2+3]));
texth1->DrawLatex(0.21, 0.84-down*2, Form("#chi^{2}/NDF = %.1f / %.0f", varh1[6], varh1[7]));
TLatex* texth2 = new TLatex(); texth2->SetNDC(); texth2->SetTextSize(0.025); texth2->SetTextColor(kRed);
texth2->DrawLatex(0.21, 0.90-down, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth2->DrawLatex(0.21, 0.87-down, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #psm %.1f", varh2[1], varh2[1+3], varh2[2], varh2[2+3]));
texth2->DrawLatex(0.21, 0.84-down, Form("#chi^{2}/NDF = %.1f / %.0f", varh2[6], varh2[7]));
TLatex* texth3 = new TLatex(); texth3->SetNDC(); texth3->SetTextSize(0.025); texth3->SetTextColor(kBlue);
texth3->DrawLatex(0.21, 0.90, Form("gaus fit [%.0f,%.0f]", xmin, xmax));
texth3->DrawLatex(0.21, 0.87, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh3[1], varh3[1+3], varh3[2], varh3[2+3]));
texth3->DrawLatex(0.21, 0.84, Form("#chi^{2}/NDF = %.1f / %.0f", varh3[6], varh3[7]));
text->DrawLatex(0.015, 0.97, "H pT > 200 GeV, MET > 200 GeV, 2 CSVL");
gPad->Print("plots/breg_FatHmassReg_met200.png");
*/

//- H jet pT resolution --------------------------------------------------------
/*
double xmin=0.8, xmax=1.2;
tree->Draw("hJet_pt[0] / hJet_genPt[0] >> h1(200,0,2)",cut+"hJet_genPt[0]>10");
tree->Draw("hJet_ptReg[0] / hJet_genPt[0] >> h2(200,0,2)",cut+"hJet_genPt[0]>10");
double varh1[8], varh2[8];
h1->Fit("gaus", "I", "", xmin, xmax); 
varh1[0] = gaus->GetParameter(0); varh1[1] = gaus->GetParameter(1); varh1[2] = gaus->GetParameter(2); varh1[3] = gaus->GetParError(0); varh1[4] = gaus->GetParError(1); varh1[5] = gaus->GetParError(2); varh1[6] = gaus->GetChisquare(); varh1[7] = gaus->GetNDF();
h2->Fit("gaus", "I", "", xmin, xmax); double h2mean = gaus->GetParameter(1); double h2sigma = gaus->GetParameter(2);
varh2[0] = gaus->GetParameter(0); varh2[1] = gaus->GetParameter(1); varh2[2] = gaus->GetParameter(2); varh2[3] = gaus->GetParError(0); varh2[4] = gaus->GetParError(1); varh2[5] = gaus->GetParError(2); varh2[6] = gaus->GetChisquare(); varh2[7] = gaus->GetNDF();

h2->SetStats(0); h2->SetLineWidth(2); h2->SetLineColor(kRed); h2->SetMarkerSize(0);
h2->SetTitle("; H jet 1 p_{T} / p_{T,gen} [GeV]"); h2->Draw("");
h1->SetStats(0); h1->SetLineWidth(2); h1->SetMarkerSize(0);
h1->Draw("same");

TLegend* leg = new TLegend(0.62, 0.78, 0.92, 0.92);
leg->AddEntry(h2, "CMS regressed");
leg->AddEntry(h1, "CMS");
leg->Draw();

double down = 0.12;
TLatex* texth1 = new TLatex(); texth1->SetNDC(); texth1->SetTextSize(0.025);
texth1->DrawLatex(0.21, 0.90-down, Form("gaus fit"));
texth1->DrawLatex(0.21, 0.87-down, Form("#mu = %.2f #pm %.2f, #sigma = %.2f #pm %.2f", varh1[1], varh1[1+3], varh1[2], varh1[2+3]));
texth1->DrawLatex(0.21, 0.84-down, Form("#chi^{2}/NDF = %.1f / %.0f", varh1[6], varh1[7]));
TLatex* texth2 = new TLatex(); texth2->SetNDC(); texth2->SetTextSize(0.025); texth2->SetTextColor(kRed);
texth2->DrawLatex(0.21, 0.90, Form("gaus fit"));
texth2->DrawLatex(0.21, 0.87, Form("#mu = %.2f #pm %.2f, #sigma = %.2f #pm %.2f", varh2[1], varh2[1+3], varh2[2], varh2[2+3]));
texth2->DrawLatex(0.21, 0.84, Form("#chi^{2}/NDF = %.1f / %.0f", varh2[6], varh2[7]));
text->DrawLatex(0.015, 0.97, "H pT > 130 GeV, MET > 170 GeV, 2 CSVL");
gPad->Print("plots/breg_hJet_ptReg0.png");

// repeat for jet 2
tree->Draw("hJet_pt[1] / hJet_genPt[1] >> h1(200,0,2)",cut+"hJet_genPt[1]>10");
tree->Draw("hJet_ptReg[1] / hJet_genPt[1] >> h2(200,0,2)",cut+"hJet_genPt[1]>10");
double varh1[8], varh2[8];
h1->Fit("gaus", "I", "", xmin, xmax); 
varh1[0] = gaus->GetParameter(0); varh1[1] = gaus->GetParameter(1); varh1[2] = gaus->GetParameter(2); varh1[3] = gaus->GetParError(0); varh1[4] = gaus->GetParError(1); varh1[5] = gaus->GetParError(2); varh1[6] = gaus->GetChisquare(); varh1[7] = gaus->GetNDF();
h2->Fit("gaus", "I", "", xmin, xmax); double h2mean = gaus->GetParameter(1); double h2sigma = gaus->GetParameter(2);
varh2[0] = gaus->GetParameter(0); varh2[1] = gaus->GetParameter(1); varh2[2] = gaus->GetParameter(2); varh2[3] = gaus->GetParError(0); varh2[4] = gaus->GetParError(1); varh2[5] = gaus->GetParError(2); varh2[6] = gaus->GetChisquare(); varh2[7] = gaus->GetNDF();

h2->SetStats(0); h2->SetLineWidth(2); h2->SetLineColor(kRed); h2->SetMarkerSize(0);
h2->SetTitle("; H jet 2 p_{T} / p_{T,gen} [GeV]"); h2->Draw("");
h1->SetStats(0); h1->SetLineWidth(2); h1->SetMarkerSize(0);
h1->Draw("same");

TLegend* leg = new TLegend(0.62, 0.78, 0.92, 0.92);
leg->AddEntry(h2, "CMS regressed");
leg->AddEntry(h1, "CMS");
leg->Draw();

double down = 0.12;
TLatex* texth1 = new TLatex(); texth1->SetNDC(); texth1->SetTextSize(0.025);
texth1->DrawLatex(0.21, 0.90-down, Form("gaus fit"));
texth1->DrawLatex(0.21, 0.87-down, Form("#mu = %.2f #pm %.2f, #sigma = %.2f #pm %.2f", varh1[1], varh1[1+3], varh1[2], varh1[2+3]));
texth1->DrawLatex(0.21, 0.84-down, Form("#chi^{2}/NDF = %.1f / %.0f", varh1[6], varh1[7]));
TLatex* texth2 = new TLatex(); texth2->SetNDC(); texth2->SetTextSize(0.025); texth2->SetTextColor(kRed);
texth2->DrawLatex(0.21, 0.90, Form("gaus fit"));
texth2->DrawLatex(0.21, 0.87, Form("#mu = %.2f #pm %.2f, #sigma = %.2f #pm %.2f", varh2[1], varh2[1+3], varh2[2], varh2[2+3]));
texth2->DrawLatex(0.21, 0.84, Form("#chi^{2}/NDF = %.1f / %.0f", varh2[6], varh2[7]));
text->DrawLatex(0.015, 0.97, "H pT > 130 GeV, MET > 170 GeV, 2 CSVL");
gPad->Print("plots/breg_hJet_ptReg1.png");
*/

//- FatH jet pT resolution -----------------------------------------------------
/*
double xmin=0.8, xmax=1.2;
tree->Draw("fathFilterJets_pt[0] / fathFilterJets_genPt[0] >> h1(200,0,2)",fjcut+"fathFilterJets_genPt[0]>10");
tree->Draw("fathFilterJets_ptReg[0] / fathFilterJets_genPt[0] >> h2(200,0,2)",fjcut+"fathFilterJets_genPt[0]>10");
double varh1[8], varh2[8];
h1->Fit("gaus", "I", "", xmin, xmax); 
varh1[0] = gaus->GetParameter(0); varh1[1] = gaus->GetParameter(1); varh1[2] = gaus->GetParameter(2); varh1[3] = gaus->GetParError(0); varh1[4] = gaus->GetParError(1); varh1[5] = gaus->GetParError(2); varh1[6] = gaus->GetChisquare(); varh1[7] = gaus->GetNDF();
h2->Fit("gaus", "I", "", xmin, xmax); double h2mean = gaus->GetParameter(1); double h2sigma = gaus->GetParameter(2);
varh2[0] = gaus->GetParameter(0); varh2[1] = gaus->GetParameter(1); varh2[2] = gaus->GetParameter(2); varh2[3] = gaus->GetParError(0); varh2[4] = gaus->GetParError(1); varh2[5] = gaus->GetParError(2); varh2[6] = gaus->GetChisquare(); varh2[7] = gaus->GetNDF();

h2->SetStats(0); h2->SetLineWidth(2); h2->SetLineColor(kRed); h2->SetMarkerSize(0);
h2->SetTitle("; H filter jet 1 p_{T} / p_{T,gen} [GeV]"); h2->Draw("");
h1->SetStats(0); h1->SetLineWidth(2); h1->SetMarkerSize(0);
h1->Draw("same");

TLegend* leg = new TLegend(0.62, 0.78, 0.92, 0.92);
leg->AddEntry(h2, "CMS regressed");
leg->AddEntry(h1, "CMS");
leg->Draw();

double down = 0.12;
TLatex* texth1 = new TLatex(); texth1->SetNDC(); texth1->SetTextSize(0.025);
texth1->DrawLatex(0.21, 0.90-down, Form("gaus fit"));
texth1->DrawLatex(0.21, 0.87-down, Form("#mu = %.2f #pm %.2f, #sigma = %.2f #pm %.2f", varh1[1], varh1[1+3], varh1[2], varh1[2+3]));
texth1->DrawLatex(0.21, 0.84-down, Form("#chi^{2}/NDF = %.1f / %.0f", varh1[6], varh1[7]));
TLatex* texth2 = new TLatex(); texth2->SetNDC(); texth2->SetTextSize(0.025); texth2->SetTextColor(kRed);
texth2->DrawLatex(0.21, 0.90, Form("gaus fit"));
texth2->DrawLatex(0.21, 0.87, Form("#mu = %.2f #pm %.2f, #sigma = %.2f #pm %.2f", varh2[1], varh2[1+3], varh2[2], varh2[2+3]));
texth2->DrawLatex(0.21, 0.84, Form("#chi^{2}/NDF = %.1f / %.0f", varh2[6], varh2[7]));
text->DrawLatex(0.015, 0.97, "H pT > 130 GeV, MET > 170 GeV, 2 CSVL");
gPad->Print("plots/breg_fathFilterJets_ptReg0.png");

// repeat for jet 2
tree->Draw("fathFilterJets_pt[1] / fathFilterJets_genPt[1] >> h1(200,0,2)",fjcut+"fathFilterJets_genPt[1]>10");
tree->Draw("fathFilterJets_ptReg[1] / fathFilterJets_genPt[1] >> h2(200,0,2)",fjcut+"fathFilterJets_genPt[1]>10");
double varh1[8], varh2[8];
h1->Fit("gaus", "I", "", xmin, xmax); 
varh1[0] = gaus->GetParameter(0); varh1[1] = gaus->GetParameter(1); varh1[2] = gaus->GetParameter(2); varh1[3] = gaus->GetParError(0); varh1[4] = gaus->GetParError(1); varh1[5] = gaus->GetParError(2); varh1[6] = gaus->GetChisquare(); varh1[7] = gaus->GetNDF();
h2->Fit("gaus", "I", "", xmin, xmax); double h2mean = gaus->GetParameter(1); double h2sigma = gaus->GetParameter(2);
varh2[0] = gaus->GetParameter(0); varh2[1] = gaus->GetParameter(1); varh2[2] = gaus->GetParameter(2); varh2[3] = gaus->GetParError(0); varh2[4] = gaus->GetParError(1); varh2[5] = gaus->GetParError(2); varh2[6] = gaus->GetChisquare(); varh2[7] = gaus->GetNDF();

h2->SetStats(0); h2->SetLineWidth(2); h2->SetLineColor(kRed); h2->SetMarkerSize(0);
h2->SetTitle("; H filter jet 2 p_{T} / p_{T,gen} [GeV]"); h2->Draw("");
h1->SetStats(0); h1->SetLineWidth(2); h1->SetMarkerSize(0);
h1->Draw("same");

TLegend* leg = new TLegend(0.62, 0.78, 0.92, 0.92);
leg->AddEntry(h2, "CMS regressed");
leg->AddEntry(h1, "CMS");
leg->Draw();

double down = 0.12;
TLatex* texth1 = new TLatex(); texth1->SetNDC(); texth1->SetTextSize(0.025);
texth1->DrawLatex(0.21, 0.90-down, Form("gaus fit"));
texth1->DrawLatex(0.21, 0.87-down, Form("#mu = %.2f #pm %.2f, #sigma = %.2f #pm %.2f", varh1[1], varh1[1+3], varh1[2], varh1[2+3]));
texth1->DrawLatex(0.21, 0.84-down, Form("#chi^{2}/NDF = %.1f / %.0f", varh1[6], varh1[7]));
TLatex* texth2 = new TLatex(); texth2->SetNDC(); texth2->SetTextSize(0.025); texth2->SetTextColor(kRed);
texth2->DrawLatex(0.21, 0.90, Form("gaus fit"));
texth2->DrawLatex(0.21, 0.87, Form("#mu = %.2f #pm %.2f, #sigma = %.2f #pm %.2f", varh2[1], varh2[1+3], varh2[2], varh2[2+3]));
texth2->DrawLatex(0.21, 0.84, Form("#chi^{2}/NDF = %.1f / %.0f", varh2[6], varh2[7]));
text->DrawLatex(0.015, 0.97, "H pT > 130 GeV, MET > 170 GeV, 2 CSVL");
gPad->Print("plots/breg_fathFilterJets_ptReg1.png");

// repeat for jet 3
tree->Draw("fathFilterJets_pt[2] / fathFilterJets_genPt[2] >> h1(100,0,2)",fjcut + "nfathFilterJets>2 && fathFilterJets_pt[2]>15 && fathFilterJets_genPt[0]>10");
tree->Draw("fathFilterJets_ptReg[2] / fathFilterJets_genPt[2] >> h2(100,0,2)",fjcut + "nfathFilterJets>2 && fathFilterJets_pt[2]>15 && fathFilterJets_genPt[0]>10");;
double varh1[8], varh2[8];
h1->Fit("gaus", "I", "", xmin, xmax); 
varh1[0] = gaus->GetParameter(0); varh1[1] = gaus->GetParameter(1); varh1[2] = gaus->GetParameter(2); varh1[3] = gaus->GetParError(0); varh1[4] = gaus->GetParError(1); varh1[5] = gaus->GetParError(2); varh1[6] = gaus->GetChisquare(); varh1[7] = gaus->GetNDF();
h2->Fit("gaus", "I", "", xmin, xmax); double h2mean = gaus->GetParameter(1); double h2sigma = gaus->GetParameter(2);
varh2[0] = gaus->GetParameter(0); varh2[1] = gaus->GetParameter(1); varh2[2] = gaus->GetParameter(2); varh2[3] = gaus->GetParError(0); varh2[4] = gaus->GetParError(1); varh2[5] = gaus->GetParError(2); varh2[6] = gaus->GetChisquare(); varh2[7] = gaus->GetNDF();

h2->SetStats(0); h2->SetLineWidth(2); h2->SetLineColor(kRed); h2->SetMarkerSize(0);
h2->SetTitle("; H filter jet 3 p_{T} / p_{T,gen} [GeV]"); h2->Draw("");
h1->SetStats(0); h1->SetLineWidth(2); h1->SetMarkerSize(0);
h1->Draw("same");

TLegend* leg = new TLegend(0.62, 0.78, 0.92, 0.92);
leg->AddEntry(h2, "CMS regressed");
leg->AddEntry(h1, "CMS");
leg->Draw();

double down = 0.12;
TLatex* texth1 = new TLatex(); texth1->SetNDC(); texth1->SetTextSize(0.025);
texth1->DrawLatex(0.21, 0.90-down, Form("gaus fit"));
texth1->DrawLatex(0.21, 0.87-down, Form("#mu = %.2f #pm %.2f, #sigma = %.2f #pm %.2f", varh1[1], varh1[1+3], varh1[2], varh1[2+3]));
texth1->DrawLatex(0.21, 0.84-down, Form("#chi^{2}/NDF = %.1f / %.0f", varh1[6], varh1[7]));
TLatex* texth2 = new TLatex(); texth2->SetNDC(); texth2->SetTextSize(0.025); texth2->SetTextColor(kRed);
texth2->DrawLatex(0.21, 0.90, Form("gaus fit"));
texth2->DrawLatex(0.21, 0.87, Form("#mu = %.2f #pm %.2f, #sigma = %.2f #pm %.2f", varh2[1], varh2[1+3], varh2[2], varh2[2+3]));
texth2->DrawLatex(0.21, 0.84, Form("#chi^{2}/NDF = %.1f / %.0f", varh2[6], varh2[7]));
text->DrawLatex(0.015, 0.97, "H pT > 130 GeV, MET > 170 GeV, 2 CSVL");
gPad->Print("plots/breg_fathFilterJets_ptReg2.png");
*/

//- Correlation scatter --------------------------------------------------------
/*
TCut cut = "Vtype==4 && hJet_pt[0]>20. && hJet_pt[1]>20. && hJet_genPt[0]>10. && hJet_genPt[1]>10. && METtype1corr.et>130 && hJet_csv_nominal[0]>0.244 && hJet_csv_nominal[1]>0.244";

pf1 = new TProfile("pf1","; H jet 1 or 2 p_{T} [GeV]; (p_{T} - p_{T,gen})/p_{T,gen}", 30, 0, 300);
pf2 = new TProfile("pf2","; H jet 1 or 2 p_{T,reg} [GeV]; (p_{T,reg} - p_{T,gen})/p_{T,gen}", 30, 0, 300);
sc1 = new TH2D("sc1","; H jet 1 or 2 p_{T} [GeV]; (p_{T} - p_{T,gen})/p_{T,gen}", 30, 0, 300, 60, -0.3, 0.3);
sc2 = new TH2D("sc2","; H jet 1 or 2 p_{T,reg} [GeV]; (p_{T,reg} - p_{T,gen})/p_{T,gen}", 30, 0, 300, 60, -0.3, 0.3);
tree->Draw("((hJet_pt - hJet_genPt)/hJet_genPt):hJet_pt >> pf1", cut, "goff");
tree->Draw("((hJet_pt - hJet_genPt)/hJet_genPt):hJet_pt >> sc1", cut, "goff");
tree->Draw("((hJet_ptReg - hJet_genPt)/hJet_genPt):hJet_ptReg >> pf2", cut, "goff");
tree->Draw("((hJet_ptReg - hJet_genPt)/hJet_genPt):hJet_ptReg >> sc2", cut, "goff");

fa1 = new TF1("fa1", "0", -9999., 9999.);
fa1->SetLineColor(kGray+2); fa1->SetLineStyle(7); fa1->SetLineWidth(2);
gStyle->SetOptStat("mr");
sc2->SetStats(0); sc2->Draw("COL");
fa1->Draw("lsame");
pf1->SetLineColor(kBlack); pf1->SetMarkerColor(kBlack); //pf1->SetMarkerSize(0.2);
pf1->Draw("sames e1");
pf2->SetLineColor(kRed); pf2->SetMarkerColor(kRed); //pf2->SetMarkerSize(0.2);
pf2->Draw("sames e1");
pf2->Fit("pol1", "NF", "", 30., 300.); pol1 = pol1;
text->SetTextSize(0.035); text->DrawLatex(0.19, 0.96, Form("linear fit: %.2g + %.2g * x (#chi^{2}=%.2f/%i)", pol1->GetParameter(0), pol1->GetParameter(1), pol1->GetChisquare(), pol1->GetNDF())); text->SetTextSize(0.025);
gPad->Print("plots/breg_correl_hJet_ptReg.png");

// repeat for H jet 1
pf1 = new TProfile("pf1_0","; H jet 1 p_{T} [GeV]; (p_{T} - p_{T,gen})/p_{T,gen}", 30, 0, 300);
pf2 = new TProfile("pf2_0","; H jet 1 p_{T,reg} [GeV]; (p_{T,reg} - p_{T,gen})/p_{T,gen}", 30, 0, 300);
sc1 = new TH2D("sc1_0","; H jet 1 p_{T} [GeV]; (p_{T} - p_{T,gen})/p_{T,gen}", 30, 0, 300, 60, -0.3, 0.3);
sc2 = new TH2D("sc2_0","; H jet 1 p_{T,reg} [GeV]; (p_{T,reg} - p_{T,gen})/p_{T,gen}", 30, 0, 300, 60, -0.3, 0.3);
tree->Draw("((hJet_pt[0] - hJet_genPt[0])/hJet_genPt[0]):hJet_pt[0] >> pf1_0", cut, "goff");
tree->Draw("((hJet_pt[0] - hJet_genPt[0])/hJet_genPt[0]):hJet_pt[0] >> sc1_0", cut, "goff");
tree->Draw("((hJet_ptReg[0] - hJet_genPt[0])/hJet_genPt[0]):hJet_ptReg[0] >> pf2_0", cut, "goff");
tree->Draw("((hJet_ptReg[0] - hJet_genPt[0])/hJet_genPt[0]):hJet_ptReg[0] >> sc2_0", cut, "goff");

sc2->SetStats(0); sc2->Draw("COL");
fa1->Draw("lsame");
pf1->SetLineColor(kBlack); pf1->SetMarkerColor(kBlack); //pf1->SetMarkerSize(0.2);
pf1->Draw("sames e1");
pf2->SetLineColor(kRed); pf2->SetMarkerColor(kRed); //pf2->SetMarkerSize(0.2);
pf2->Draw("sames e1");
pf2->Fit("pol1", "NF", "", 30., 300.); pol1 = pol1;
text->SetTextSize(0.035); text->DrawLatex(0.19, 0.96, Form("linear fit: %.2g + %.2g * x (#chi^{2}=%.2f/%i)", pol1->GetParameter(0), pol1->GetParameter(1), pol1->GetChisquare(), pol1->GetNDF())); text->SetTextSize(0.025);
gPad->Print("plots/breg_correl_hJet_ptReg0.png");

// repeat for H jet 2
pf1 = new TProfile("pf1_1","; H jet 2 p_{T} [GeV]; (p_{T} - p_{T,gen})/p_{T,gen}", 30, 0, 300);
pf2 = new TProfile("pf2_1","; H jet 2 p_{T,reg} [GeV]; (p_{T,reg} - p_{T,gen})/p_{T,gen}", 30, 0, 300);
sc1 = new TH2D("sc1_1","; H jet 2 p_{T} [GeV]; (p_{T} - p_{T,gen})/p_{T,gen}", 30, 0, 300, 60, -0.3, 0.3);
sc2 = new TH2D("sc2_1","; H jet 2 p_{T,reg} [GeV]; (p_{T,reg} - p_{T,gen})/p_{T,gen}", 30, 0, 300, 60, -0.3, 0.3);
tree->Draw("((hJet_pt[1] - hJet_genPt[1])/hJet_genPt[1]):hJet_pt[1] >> pf1_1", cut, "goff");
tree->Draw("((hJet_pt[1] - hJet_genPt[1])/hJet_genPt[1]):hJet_pt[1] >> sc1_1", cut, "goff");
tree->Draw("((hJet_ptReg[1] - hJet_genPt[1])/hJet_genPt[1]):hJet_ptReg[1] >> pf2_1", cut, "goff");
tree->Draw("((hJet_ptReg[1] - hJet_genPt[1])/hJet_genPt[1]):hJet_ptReg[1] >> sc2_1", cut, "goff");

sc2->SetStats(0); sc2->Draw("COL");
fa1->Draw("lsame");
pf1->SetLineColor(kBlack); pf1->SetMarkerColor(kBlack); //pf1->SetMarkerSize(0.2);
pf1->Draw("sames e1");
pf2->SetLineColor(kRed); pf2->SetMarkerColor(kRed); //pf2->SetMarkerSize(0.2);
pf2->Draw("sames e1");
pf2->Fit("pol1", "NF", "", 30., 300.); pol1 = pol1;
text->SetTextSize(0.035); text->DrawLatex(0.19, 0.96, Form("linear fit: %.2g + %.2g * x (#chi^{2}=%.2f/%i)", pol1->GetParameter(0), pol1->GetParameter(1), pol1->GetChisquare(), pol1->GetNDF())); text->SetTextSize(0.025);
gPad->Print("plots/breg_correl_hJet_ptReg1.png");
*/

//- Correlation scatter (filter jets) ------------------------------------------
/*
TCut fjcut = "Vtype==4 && nfathFilterJets>0 && fathFilterJets_pt[0]>20 && fathFilterJets_pt[1]>20 && fathFilterJets_genPt[0]>10 && fathFilterJets_genPt[1]>10 && METtype1corr.et>130 && fathFilterJets_csv[0]>0.244 && fathFilterJets_csv[1]>0.244";

pf1 = new TProfile("pf1","; H filter jet 1 or 2 p_{T} [GeV]; (p_{T} - p_{T,gen})/p_{T,gen}", 30, 0, 300);
pf2 = new TProfile("pf2","; H filter jet 1 or 2 p_{T,reg} [GeV]; (p_{T,reg} - p_{T,gen})/p_{T,gen}", 30, 0, 300);
sc1 = new TH2D("sc1","; H filter jet 1 or 2 p_{T} [GeV]; (p_{T} - p_{T,gen})/p_{T,gen}", 30, 0, 300, 60, -0.3, 0.3);
sc2 = new TH2D("sc2","; H filter jet 1 or 2 p_{T,reg} [GeV]; (p_{T,reg} - p_{T,gen})/p_{T,gen}", 30, 0, 300, 60, -0.3, 0.3);
tree->Draw("((fathFilterJets_pt[0] - fathFilterJets_genPt[0])/fathFilterJets_genPt[0]):fathFilterJets_pt[0] >> pf1", fjcut, "goff");
tree->Draw("((fathFilterJets_pt[1] - fathFilterJets_genPt[1])/fathFilterJets_genPt[1]):fathFilterJets_pt[1] >> +pf1", fjcut, "goff");
tree->Draw("((fathFilterJets_pt[0] - fathFilterJets_genPt[0])/fathFilterJets_genPt[0]):fathFilterJets_pt[0] >> sc1", fjcut, "goff");
tree->Draw("((fathFilterJets_pt[1] - fathFilterJets_genPt[1])/fathFilterJets_genPt[1]):fathFilterJets_pt[1] >> +sc1", fjcut, "goff");
tree->Draw("((fathFilterJets_ptReg[0] - fathFilterJets_genPt[0])/fathFilterJets_genPt[0]):fathFilterJets_ptReg[0] >> pf2", fjcut, "goff");
tree->Draw("((fathFilterJets_ptReg[1] - fathFilterJets_genPt[1])/fathFilterJets_genPt[1]):fathFilterJets_ptReg[1] >> +pf2", fjcut, "goff");
tree->Draw("((fathFilterJets_ptReg[0] - fathFilterJets_genPt[0])/fathFilterJets_genPt[0]):fathFilterJets_ptReg[0] >> sc2", fjcut, "goff");
tree->Draw("((fathFilterJets_ptReg[1] - fathFilterJets_genPt[1])/fathFilterJets_genPt[1]):fathFilterJets_ptReg[1] >> +sc2", fjcut, "goff");

fa1 = new TF1("fa1", "0", -9999., 9999.);
fa1->SetLineColor(kGray+2); fa1->SetLineStyle(7); fa1->SetLineWidth(2);
gStyle->SetOptStat("mr");
sc2->SetStats(0); sc2->Draw("COL");
fa1->Draw("lsame");
pf1->SetLineColor(kBlack); pf1->SetMarkerColor(kBlack); //pf1->SetMarkerSize(0.2);
pf1->Draw("sames e1");
pf2->SetLineColor(kRed); pf2->SetMarkerColor(kRed); //pf2->SetMarkerSize(0.2);
pf2->Draw("sames e1");
pf2->Fit("pol1", "NF", "", 20., 300.); pol1 = pol1;
text->SetTextSize(0.035); text->DrawLatex(0.21, 0.96, Form("linear fit: %.2g + %.2g * x (#chi^{2}=%.2f/%i)", pol1->GetParameter(0), pol1->GetParameter(1), pol1->GetChisquare(), pol1->GetNDF())); text->SetTextSize(0.025);
gPad->Print("plots/breg_correl_fathFilterJets_ptReg.png");

// repeat for H jet 1
pf1 = new TProfile("pf1_0","; H filter jet 1 p_{T} [GeV]; (p_{T} - p_{T,gen})/p_{T}", 30, 0, 300);
pf2 = new TProfile("pf2_0","; H filter jet 1 p_{T,reg} [GeV]; (p_{T,reg} - p_{T,gen})/p_{T,reg}", 30, 0, 300);
sc1 = new TH2D("sc1_0","; H filter jet 1 p_{T} [GeV]; (p_{T} - p_{T,gen})/p_{T}", 30, 0, 300, 60, -0.3, 0.3);
sc2 = new TH2D("sc2_0","; H filter jet 1 p_{T,reg} [GeV]; (p_{T,reg} - p_{T,gen})/p_{T,reg}", 30, 0, 300, 60, -0.3, 0.3);
tree->Draw("((fathFilterJets_pt[0] - fathFilterJets_genPt[0])/fathFilterJets_genPt[0]):fathFilterJets_pt[0] >> pf1_0", fjcut, "goff");
tree->Draw("((fathFilterJets_pt[0] - fathFilterJets_genPt[0])/fathFilterJets_genPt[0]):fathFilterJets_pt[0] >> sc1_0", fjcut, "goff");
tree->Draw("((fathFilterJets_ptReg[0] - fathFilterJets_genPt[0])/fathFilterJets_genPt[0]):fathFilterJets_ptReg[0] >> pf2_0", fjcut, "goff");
tree->Draw("((fathFilterJets_ptReg[0] - fathFilterJets_genPt[0])/fathFilterJets_genPt[0]):fathFilterJets_ptReg[0] >> sc2_0", fjcut, "goff");

sc2->SetStats(0); sc2->Draw("COL");
fa1->Draw("lsame");
pf1->SetLineColor(kBlack); pf1->SetMarkerColor(kBlack); //pf1->SetMarkerSize(0.2);
pf1->Draw("sames e1");
pf2->SetLineColor(kRed); pf2->SetMarkerColor(kRed); //pf2->SetMarkerSize(0.2);
pf2->Draw("sames e1");
pf2->Fit("pol1", "NF", "", 20., 300.); pol1 = pol1;
text->SetTextSize(0.035); text->DrawLatex(0.21, 0.96, Form("linear fit: %.2g + %.2g * x (#chi^{2}=%.2f/%i)", pol1->GetParameter(0), pol1->GetParameter(1), pol1->GetChisquare(), pol1->GetNDF())); text->SetTextSize(0.025);
gPad->Print("plots/breg_correl_fathFilterJets_ptReg0.png");

// repeat for H jet 2
pf1 = new TProfile("pf1_1","; H filter jet 2 p_{T} [GeV]; (p_{T} - p_{T,gen})/p_{T,gen}", 30, 0, 300);
pf2 = new TProfile("pf2_1","; H filter jet 2 p_{T,reg} [GeV]; (p_{T,reg} - p_{T,gen})/p_{T,gen}", 30, 0, 300);
sc1 = new TH2D("sc1_1","; H filter jet 2 p_{T} [GeV]; (p_{T} - p_{T,gen})/p_{T,gen}", 30, 0, 300, 60, -0.3, 0.3);
sc2 = new TH2D("sc2_1","; H filter jet 2 p_{T,reg} [GeV]; (p_{T,reg} - p_{T,gen})/p_{T,gen}", 30, 0, 300, 60, -0.3, 0.3);
tree->Draw("((fathFilterJets_pt[1] - fathFilterJets_genPt[1])/fathFilterJets_genPt[1]):fathFilterJets_pt[1] >> pf1_1", fjcut, "goff");
tree->Draw("((fathFilterJets_pt[1] - fathFilterJets_genPt[1])/fathFilterJets_genPt[1]):fathFilterJets_pt[1] >> sc1_1", fjcut, "goff");
tree->Draw("((fathFilterJets_ptReg[1] - fathFilterJets_genPt[1])/fathFilterJets_genPt[1]):fathFilterJets_ptReg[1] >> pf2_1", fjcut, "goff");
tree->Draw("((fathFilterJets_ptReg[1] - fathFilterJets_genPt[1])/fathFilterJets_genPt[1]):fathFilterJets_ptReg[1] >> sc2_1", fjcut, "goff");

sc2->SetStats(0); sc2->Draw("COL");
fa1->Draw("lsame");
pf1->SetLineColor(kBlack); pf1->SetMarkerColor(kBlack); //pf1->SetMarkerSize(0.2);
pf1->Draw("sames e1");
pf2->SetLineColor(kRed); pf2->SetMarkerColor(kRed); //pf2->SetMarkerSize(0.2);
pf2->Draw("sames e1");
pf2->Fit("pol1", "NF", "", 20., 300.); pol1 = pol1;
text->SetTextSize(0.035); text->DrawLatex(0.21, 0.96, Form("linear fit: %.2g + %.2g * x (#chi^{2}=%.2f/%i)", pol1->GetParameter(0), pol1->GetParameter(1), pol1->GetChisquare(), pol1->GetNDF())); text->SetTextSize(0.025);
gPad->Print("plots/breg_correl_fathFilterJets_ptReg1.png");

// repeat for H jet 3
fjcut += "nfathFilterJets>2 && fathFilterJets_pt[2]>20 && fathFilterJets_genPt[2]>10";
pf1 = new TProfile("pf1_2","; H filter jet 3 p_{T} [GeV]; (p_{T} - p_{T,gen})/p_{T,gen}", 30, 0, 300);
pf2 = new TProfile("pf2_2","; H filter jet 3 p_{T,reg} [GeV]; (p_{T,reg} - p_{T,gen})/p_{T,gen}", 30, 0, 300);
sc1 = new TH2D("sc1_2","; H filter jet 3 p_{T} [GeV]; (p_{T} - p_{T,gen})/p_{T,gen}", 30, 0, 300, 60, -0.3, 0.3);
sc2 = new TH2D("sc2_2","; H filter jet 3 p_{T,reg} [GeV]; (p_{T,reg} - p_{T,gen})/p_{T,gen}", 30, 0, 300, 60, -0.3, 0.3);
tree->Draw("((fathFilterJets_pt[2] - fathFilterJets_genPt[2])/fathFilterJets_genPt[2]):fathFilterJets_pt[2] >> pf1_2", fjcut, "goff");
tree->Draw("((fathFilterJets_pt[2] - fathFilterJets_genPt[2])/fathFilterJets_genPt[2]):fathFilterJets_pt[2] >> sc1_2", fjcut, "goff");
tree->Draw("((fathFilterJets_ptReg[2] - fathFilterJets_genPt[2])/fathFilterJets_genPt[2]):fathFilterJets_ptReg[2] >> pf2_2", fjcut, "goff");
tree->Draw("((fathFilterJets_ptReg[2] - fathFilterJets_genPt[2])/fathFilterJets_genPt[2]):fathFilterJets_ptReg[2] >> sc2_2", fjcut, "goff");

sc2->SetStats(0); sc2->Draw("COL");
fa1->Draw("lsame");
pf1->SetLineColor(kBlack); pf1->SetMarkerColor(kBlack); //pf1->SetMarkerSize(0.2);
pf1->Draw("sames e1");
pf2->SetLineColor(kRed); pf2->SetMarkerColor(kRed); //pf2->SetMarkerSize(0.2);
pf2->Draw("sames e1");
pf2->Fit("pol1", "NF", "", 20., 300.); pol1 = pol1;
text->SetTextSize(0.035); text->DrawLatex(0.21, 0.96, Form("linear fit: %.2g + %.2g * x (#chi^{2}=%.2f/%i)", pol1->GetParameter(0), pol1->GetParameter(1), pol1->GetChisquare(), pol1->GetNDF())); text->SetTextSize(0.025);
gPad->Print("plots/breg_correl_fathFilterJets_ptReg2.png");
*/

//- Test vs. Train -------------------------------------------------------------
/*
//gROOT->SetBatch(); const Int_t c_SignalLine = TColor::GetColor( "#0000ee" ); const Int_t c_SignalFill = TColor::GetColor( "#7d99d1" ); 
_file2 = TFile::Open("testTMVARegFJ3.root");
TrainTree->Draw("BDTG >> hbdt1(120,20,320)");
TestTree->Draw("BDTG >> hbdt2(120,20,320)");
hbdt1->SetMarkerSize(0.7); hbdt1->SetMarkerStyle(20); hbdt1->SetLineWidth(1); //hbdt1->SetLineColor(c_SignalLine); hbdt1->SetMarkerColor(c_SignalLine);
hbdt2->SetFillColor(c_SignalFill); hbdt2->SetFillStyle(1001); hbdt2->SetMarkerSize(0); hbdt2->SetLineColor(c_SignalLine); hbdt2->SetLineWidth(2);
hbdt1=hbdt1; hbdt2=hbdt2; hbdt1->Sumw2(); hbdt2->Sumw2(); hbdt1->Scale(1./hbdt1->Integral()); hbdt2->Scale(1./hbdt2->Integral()); 
c1->SetLogy(); hbdt2->SetTitle(";p_{T, regressed} [GeV]"); hbdt2->Draw("hist");
hbdt1->Draw("e1same");

int percent = 40;
TLegend* leg = new TLegend(0.21, 0.23, 0.64, 0.3);
leg->AddEntry(hbdt1, Form("Training (~%i%% MC)", percent));
leg->AddEntry(hbdt2, Form("Testing (~%i%% MC)", 100-percent));
leg->Draw();

hbdt1 = hbdt1;
double kol = hbdt2->KolmogorovTest( hbdt1 );
TLatex* text = new TLatex(); text->SetNDC(); text->SetTextSize(0.025);
text->DrawLatex(0.21, 0.2, Form("Kolmogorov-Smirnov test: signal probability = %5.3g", kol));
gPad->Print("plots/breg_fjovertrain_BDTG.png");
c1->SetLogy(0);
*/

//- genPt vs. quarks -----------------------------------------------------------
/*
tree->Draw("hJet_genPt[0] / genB.pt >> h1(200,0,2)", "70 < hJet_genPt[0] && hJet_genPt[0] < 100 && deltaR(hJet_genEta[0], hJet_genPhi[0], genB.eta, genB.phi)<0.3", "goff");
tree->Draw("hJet_genPt[0] / genBbar.pt >> h2(200,0,2)", "70 < hJet_genPt[0] && hJet_genPt[0] < 100 && deltaR(hJet_genEta[0], hJet_genPhi[0], genBbar.eta, genBbar.phi)<0.3", "goff");
tree->Draw("fathFilterJets_genPt[0] / genB.pt >> h3(200,0,2)", "70 < fathFilterJets_genPt[0] && fathFilterJets_genPt[0] < 100 && deltaR(fathFilterJets_genEta[0], fathFilterJets_genPhi[0], genB.eta, genB.phi)<0.3", "goff");
tree->Draw("fathFilterJets_genPt[0] / genBbar.pt >> h4(200,0,2)", "70 < fathFilterJets_genPt[0] && fathFilterJets_genPt[0] < 100 && deltaR(fathFilterJets_genEta[0], fathFilterJets_genPhi[0], genBbar.eta, genBbar.phi)<0.3", "goff");
h1->Add(h2);
h3->Add(h4);
h1->SetStats(0); h1->SetLineWidth(2); h1->SetLineColor(kCyan); h1->SetMarkerSize(0);
h3->SetStats(0); h3->SetLineWidth(2); h3->SetLineColor(kRed); h3->SetMarkerSize(0);
h1->SetTitle("70 < p_{T,gen-jet} < 100;p_{T,gen-jet 1}/p_{T,parton}"); h1->Draw();
h3->Draw("same");
TLegend* leg = new TLegend(0.58, 0.78, 0.94, 0.92);
leg->AddEntry(h1, "anti-k_{T} R=0.5 gen-jet");
leg->AddEntry(h3, "k_{T} R=0.3 gen-filtered jet");
leg->Draw();
gPad->Print("plots/breg_genPtResolution1.png");

tree->Draw("hJet_genPt[0] / genB.pt >> h1(200,0,2)", "100 < hJet_genPt[0] && hJet_genPt[0] < 200 && deltaR(hJet_genEta[0], hJet_genPhi[0], genB.eta, genB.phi)<0.3", "goff");
tree->Draw("hJet_genPt[0] / genBbar.pt >> h2(200,0,2)", "100 < hJet_genPt[0] && hJet_genPt[0] < 200 && deltaR(hJet_genEta[0], hJet_genPhi[0], genBbar.eta, genBbar.phi)<0.3", "goff");
tree->Draw("fathFilterJets_genPt[0] / genB.pt >> h3(200,0,2)", "100 < fathFilterJets_genPt[0] && fathFilterJets_genPt[0] < 200 && deltaR(fathFilterJets_genEta[0], fathFilterJets_genPhi[0], genB.eta, genB.phi)<0.3", "goff");
tree->Draw("fathFilterJets_genPt[0] / genBbar.pt >> h4(200,0,2)", "100 < fathFilterJets_genPt[0] && fathFilterJets_genPt[0] < 200 && deltaR(fathFilterJets_genEta[0], fathFilterJets_genPhi[0], genBbar.eta, genBbar.phi)<0.3", "goff");
h1->Add(h2);
h3->Add(h4);
h1->SetStats(0); h1->SetLineWidth(2); h1->SetLineColor(kCyan); h1->SetMarkerSize(0);
h3->SetStats(0); h3->SetLineWidth(2); h3->SetLineColor(kRed); h3->SetMarkerSize(0);
h1->SetTitle("100 < p_{T,gen-jet} < 200;p_{T,gen-jet 1}/p_{T,parton}"); h1->Draw();
h3->Draw("same");
TLegend* leg = new TLegend(0.58, 0.78, 0.94, 0.92);
leg->AddEntry(h1, "anti-k_{T} R=0.5 gen-jet");
leg->AddEntry(h3, "k_{T} R=0.3 gen-filtered jet");
leg->Draw();
gPad->Print("plots/breg_genPtResolution2.png");

tree->Draw("hJet_genPt[0] / genB.pt >> h1(200,0,2)", "200 < hJet_genPt[0] && deltaR(hJet_genEta[0], hJet_genPhi[0], genB.eta, genB.phi)<0.3", "goff");
tree->Draw("hJet_genPt[0] / genBbar.pt >> h2(200,0,2)", "200 < hJet_genPt[0] && deltaR(hJet_genEta[0], hJet_genPhi[0], genBbar.eta, genBbar.phi)<0.3", "goff");
tree->Draw("fathFilterJets_genPt[0] / genB.pt >> h3(200,0,2)", "200 < fathFilterJets_genPt[0] && deltaR(fathFilterJets_genEta[0], fathFilterJets_genPhi[0], genB.eta, genB.phi)<0.3", "goff");
tree->Draw("fathFilterJets_genPt[0] / genBbar.pt >> h4(200,0,2)", "200 < fathFilterJets_genPt[0] && deltaR(fathFilterJets_genEta[0], fathFilterJets_genPhi[0], genBbar.eta, genBbar.phi)<0.3", "goff");
h1->Add(h2);
h3->Add(h4);
h1->SetStats(0); h1->SetLineWidth(2); h1->SetLineColor(kCyan); h1->SetMarkerSize(0);
h3->SetStats(0); h3->SetLineWidth(2); h3->SetLineColor(kRed); h3->SetMarkerSize(0);
h1->SetTitle("200 < p_{T,gen-jet};p_{T,gen-jet 1}/p_{T,parton}"); h1->Draw();
h3->Draw("same");
TLegend* leg = new TLegend(0.58, 0.78, 0.94, 0.92);
leg->AddEntry(h1, "anti-k_{T} R=0.5 gen-jet");
leg->AddEntry(h3, "k_{T} R=0.3 gen-filtered jet");
leg->Draw();
gPad->Print("plots/breg_genPtResolution3.png");
*/

//- Test vs. Train Mjj ---------------------------------------------------------
/*
//gROOT->SetBatch(); const Int_t c_SignalLine = TColor::GetColor( "#0000ee" ); const Int_t c_SignalFill = TColor::GetColor( "#7d99d1" ); 
//_file0 = TFile::Open("skim_ZnunuHbb_regression/TMVARegApp/TMVARegApp_ZnunuHbb125_training.root"); _file1 = TFile::Open("skim_ZnunuHbb_regression/TMVARegApp/TMVARegApp_ZnunuHbb125_testing.root");
_file0 = TFile::Open("skim_ZnunuHbb_regression/TMVARegApp/TMVARegApp_ZnunuHbb125_traininglow.root"); _file1 = TFile::Open("skim_ZnunuHbb_regression/TMVARegApp/TMVARegApp_ZnunuHbb125_testinglow.root");
hbdt1 = (TH1F*) _file0->Get("breg_app_mjj_j1"); hbdt2 = (TH1F*) _file1->Get("breg_app_mjj_j1");
hbdt1->SetMarkerSize(0.7); hbdt1->SetMarkerStyle(20); hbdt1->SetLineWidth(1); //hbdt1->SetLineColor(c_SignalLine); hbdt1->SetMarkerColor(c_SignalLine); 
hbdt2->SetFillColor(c_SignalFill); hbdt2->SetFillStyle(1001); hbdt2->SetMarkerSize(0); hbdt2->SetLineColor(c_SignalLine); hbdt2->SetLineWidth(2);
hbdt1=hbdt1; hbdt2=hbdt2; hbdt1->Sumw2(); hbdt2->Sumw2(); hbdt1->Scale(1./hbdt1->Integral()); hbdt2->Scale(1./hbdt2->Integral()); 

double xmin=110., xmax=140.;
double varh1[8], varh2[8], varh3[8];
hbdt1->Fit("gaus", "I", "", xmin, xmax); 
varh1[0] = gaus->GetParameter(0); varh1[1] = gaus->GetParameter(1); varh1[2] = gaus->GetParameter(2); varh1[3] = gaus->GetParError(0); varh1[4] = gaus->GetParError(1); varh1[5] = gaus->GetParError(2); varh1[6] = gaus->GetChisquare(); varh1[7] = gaus->GetNDF();
hbdt2->Fit("gaus", "I", "", xmin, xmax); double h2mean = gaus->GetParameter(1); double h2sigma = gaus->GetParameter(2);
varh2[0] = gaus->GetParameter(0); varh2[1] = gaus->GetParameter(1); varh2[2] = gaus->GetParameter(2); varh2[3] = gaus->GetParError(0); varh2[4] = gaus->GetParError(1); varh2[5] = gaus->GetParError(2); varh2[6] = gaus->GetChisquare(); varh2[7] = gaus->GetNDF();

hbdt2->SetTitle("; H(125) mass [GeV]"); hbdt2->SetMaximum(0.2); hbdt2->Draw("hist");
hbdt1->Draw("e1same");

int percent = 20;
TLegend* leg = new TLegend(0.21, 0.23, 0.64, 0.3);
leg->AddEntry(hbdt1, Form("Training (~%i%% MC)", percent));
leg->AddEntry(hbdt2, Form("Testing (~%i%% MC)", 100-percent));
leg->Draw();

double down = 0.12;
TLatex* texth1 = new TLatex(); texth1->SetNDC(); texth1->SetTextSize(0.025);
texth1->DrawLatex(0.21, 0.90, Form("gaus fit (Training) [%.0f,%.0f]", xmin, xmax));
texth1->DrawLatex(0.21, 0.87, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh1[1], varh1[1+3], varh1[2], varh1[2+3]));
texth1->DrawLatex(0.21, 0.84, Form("#chi^{2}/NDF = %.1f / %.0f", varh1[6], varh1[7]));
TLatex* texth2 = new TLatex(); texth2->SetNDC(); texth2->SetTextSize(0.025); texth2->SetTextColor(c_SignalLine);
texth2->DrawLatex(0.21, 0.90-down, Form("gaus fit (Testing) [%.0f,%.0f]", xmin, xmax));
texth2->DrawLatex(0.21, 0.87-down, Form("#mu = %.1f #pm %.1f, #sigma = %.1f #pm %.1f", varh2[1], varh2[1+3], varh2[2], varh2[2+3]));
texth2->DrawLatex(0.21, 0.84-down, Form("#chi^{2}/NDF = %.1f / %.0f", varh2[6], varh2[7]));

hbdt1 = hbdt1;
double kol = hbdt2->KolmogorovTest( hbdt1 );
TLatex* text = new TLatex(); text->SetNDC(); text->SetTextSize(0.025);
text->DrawLatex(0.21, 0.2, Form("Kolmogorov-Smirnov test: signal probability = %5.3g", kol));
//text->DrawLatex(0.015, 0.97, "H pT > 160 GeV, MET > 160 GeV, 2 CSVL");
//gPad->Print("plots/breg_overtrain_mjj_BDTG.png");
text->DrawLatex(0.015, 0.97, "H pT > 120 GeV, 120 < MET < 160 GeV, 2 CSVL");
gPad->Print("plots/breg_overtrain_mjj_BDTG_lowmet.png");
*/

//- Correlation Matrix ---------------------------------------------------------
// from TMVA correlations.C
/*
_file2 = TFile::Open("testTMVARegFJ3.root");
TLatex* text = new TLatex(); text->SetNDC(); text->SetTextSize(0.025);
h2 = CorrelationMatrix;
c1 = new TCanvas("c1","c1",1300,700);
gPad->SetGrid();
gPad->SetTicks();
Float_t newMargin1 = 0.13;
Float_t newMargin2 = 0.13;
gPad->SetLeftMargin  ( newMargin2 );
gPad->SetBottomMargin( newMargin2 );
gPad->SetRightMargin ( newMargin1 );
gPad->SetTopMargin   ( newMargin1 );
gStyle->SetPalette( 1, 0 );
gStyle->SetPaintTextFormat( "3g" );

h2->SetStats(0);
h2->SetMarkerSize( 1.5 );
h2->SetMarkerColor( 0 );
Float_t labelSize = 0.040;
h2->GetXaxis()->SetLabelSize( labelSize );
h2->GetYaxis()->SetLabelSize( labelSize );
h2->LabelsOption( "d" );
h2->SetLabelOffset( 0.011 );// label offset on x axis
h2->Draw("colz"); // color pads   
gPad->Modify(); gPad->Update();

// modify properties of paletteAxis
TPaletteAxis* paletteAxis = (TPaletteAxis*)h2->GetListOfFunctions()->FindObject( "palette" );
paletteAxis->SetLabelSize( 0.03 );
paletteAxis->SetX1NDC( paletteAxis->GetX1NDC() + 0.02 );
h2->Draw("textsame");  // add text
text->DrawLatex(0.53, 0.905, "Linear correlation coefficients in %" );
gPad->Modified(); gPad->Update();
gPad->Print("plots/breg_fjcorrelationmatrix.png");
*/

}
