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

//gROOT->SetBatch(1);

//TFile* _file0 = TFile::Open("skim_ZnnH_baseline/Step3_ZnnH125.root");
TFile* _file0 = TFile::Open("dcache:/pnfs/cms/WAX/11/store/user/lpchbb/apana/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V4a_MC/DiJetPt_ZH_ZToNuNu_HToBB_M-125_8TeV-powheg-herwigpp.root");
TTree* tree = (TTree*) _file0->Get("tree");


//- # of filter jets -----------------------------------------------------------
//TCut cut = "Vtype==4 && hJet_pt[0]>60 && hJet_pt[1]>30 && METtype1corr.et>170 && H.pt>130 && hJet_csv[0]>0.244 && hJet_csv[1]>0.244 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && nfathFilterJets>0 && fathFilterJets_pt[0]>0 && fathFilterJets_pt[1]>0 && H.mass>0";
TCut cut = "Vtype==4 && hJet_pt[0]>60 && hJet_pt[1]>30 && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5 && METtype1corr.et>130 && hJet_csv_nominal[0]>0.244 && hJet_csv_nominal[1]>0.244 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && H.mass>0";
/*
tree->Draw("H.pt >> h1(40,0,400)",cut);
tree->Draw("H.pt >> h2(40,0,400)",cut + "FatH.FatHiggsFlag==1 && fathFilterJets_pt[0]>15 && fathFilterJets_pt[1]>15 &&  nfathFilterJets>=2");
tree->Draw("H.pt >> h3(40,0,400)",cut + "FatH.FatHiggsFlag==1 && fathFilterJets_pt[0]>15 && fathFilterJets_pt[1]>15 && (nfathFilterJets==2 || (nfathFilterJets==3 && Alt$(fathFilterJets_pt[2],0)<=15.))");
tree->Draw("H.pt >> h4(40,0,400)",cut + "FatH.FatHiggsFlag==1 && fathFilterJets_pt[0]>15 && fathFilterJets_pt[1]>15 && nfathFilterJets==3 && fathFilterJets_pt[2]>15");

h1=h1;
h1->SetStats(0); h1->SetLineWidth(2); h1->SetMarkerSize(0);
h2->SetStats(0); h2->SetLineWidth(2); h2->SetLineColor(kBlack); h2->SetMarkerSize(0);
h2->Divide(h1);
h2->SetTitle("; H(125) p_{T} [GeV]"); h2->SetMinimum(0.); h2->SetMaximum(1.15); h2->Draw();
h3->SetStats(0); h3->SetLineWidth(2); h3->SetLineColor(kRed); h3->SetMarkerSize(0);
h3->Divide(h1);
h3->Draw("same");
h4->SetStats(0); h4->SetLineWidth(2); h4->SetLineColor(kBlue); h4->SetMarkerSize(0);
h4->Divide(h1);
h4->Draw("same");

TLegend* leg = new TLegend(0.18, 0.75, 0.48, 0.92);
leg->SetFillStyle(0);
leg->AddEntry(h2, "# filter jets >= 2");
leg->AddEntry(h3, "# filter jets == 2");
leg->AddEntry(h4, "# filter jets == 3");
leg->Draw();

text->DrawLatex(0.015, 0.97, "MET > 130 , j1,j2^{ak5} > 60,30 , 2 CSVL");
gPad->Print("plots/sfj_nfathFilterJets.png");


h2=h2;
h3->Divide(h2);
h4->Divide(h2);
h4->SetTitle("; H(125) p_{T} [GeV]"); h4->SetMinimum(0.); h4->SetMaximum(1.15); h4->Draw();
h3->Draw("same");

TLegend* leg = new TLegend(0.18, 0.75, 0.48, 0.92);
leg->SetFillStyle(0);
leg->AddEntry(h3, "# filter jets == 2");
leg->AddEntry(h4, "# filter jets == 3");
leg->Draw();

text->DrawLatex(0.015, 0.97, "MET > 130 , j1,j2^{ak5} > 60,30 , 2 CSVL");
gPad->Print("plots/sfj_nfathFilterJets_1.png");
*/

//------------------------------------------------------------------------------

TCut gencut = "Vtype==4 && (genB.pt>60||genBbar.pt>60) && (genB.pt>30&&genBbar.pt>30) && abs(genB.eta)<2.5 && abs(genBbar.eta)<2.5 && METtype1corr.et>130 && abs(deltaPhi(genH.phi,METtype1corr.phi))>2.0";

tree->Draw("genH.pt >> h1(40,0,400)",gencut);
tree->Draw("genH.pt >> h2(40,0,400)",gencut + "FatH.FatHiggsFlag==1 && fathFilterJets_pt[0]>15 && fathFilterJets_pt[1]>15 && nfathFilterJets>=2");
tree->Draw("genH.pt >> h3(40,0,400)",gencut + "FatH.FatHiggsFlag==1 && fathFilterJets_pt[0]>15 && fathFilterJets_pt[1]>15 && (nfathFilterJets==2 || (nfathFilterJets==3 && Alt$(fathFilterJets_pt[2],0)<=15.))");
tree->Draw("genH.pt >> h4(40,0,400)",gencut + "FatH.FatHiggsFlag==1 && fathFilterJets_pt[0]>15 && fathFilterJets_pt[1]>15 && nfathFilterJets==3 && fathFilterJets_pt[2]>15");
tree->Draw("genH.pt >> h5(40,0,400)",gencut + "hJet_pt[0]>20 && hJet_pt[1]>20");

h1=h1;
h1->SetStats(0); h1->SetLineWidth(2); h1->SetMarkerSize(0);
h2->SetStats(0); h2->SetLineWidth(2); h2->SetLineColor(kBlack); h2->SetMarkerSize(0);
h2->Divide(h1);
h2->SetTitle("; H(125)^{gen} p_{T} [GeV]"); h2->SetMinimum(0.); h2->SetMaximum(1.15); h2->Draw();
h3->SetStats(0); h3->SetLineWidth(2); h3->SetLineColor(kRed); h3->SetMarkerSize(0);
h3->Divide(h1);
h3->Draw("same");
h4->SetStats(0); h4->SetLineWidth(2); h4->SetLineColor(kBlue); h4->SetMarkerSize(0);
h4->Divide(h1);
h4->Draw("same");
h5->SetStats(0); h5->SetLineWidth(2); h5->SetLineColor(kCyan); h5->SetMarkerSize(0);
h5->Divide(h1);
h5->Draw("same");

TLegend* leg = new TLegend(0.18, 0.75, 0.48, 0.92);
leg->SetFillStyle(0);
leg->AddEntry(h2, "# filter jets >= 2");
leg->AddEntry(h3, "# filter jets == 2");
leg->AddEntry(h4, "# filter jets == 3");
leg->AddEntry(h5, "# ak5 jets == 2");
leg->Draw();

text->DrawLatex(0.015, 0.97, "MET > 130 , b1,b2^{gen} > 60,30");
gPad->Print("plots/sfj_nfathFilterJets_gencut.png");

h2=h2;
h3->Divide(h2);
h4->Divide(h2);
h4->SetTitle("; H(125)^{gen} p_{T} [GeV]"); h4->SetMinimum(0.); h4->SetMaximum(1.15); h4->Draw();
h3->Draw("same");

TLegend* leg = new TLegend(0.18, 0.75, 0.48, 0.92);
leg->SetFillStyle(0);
leg->AddEntry(h3, "# filter jets == 2");
leg->AddEntry(h4, "# filter jets == 3");
leg->Draw();

text->DrawLatex(0.015, 0.97, "MET > 130 , b1,b2^{gen} > 60,30");
gPad->Print("plots/sfj_nfathFilterJets_gencut_1.png");


}
