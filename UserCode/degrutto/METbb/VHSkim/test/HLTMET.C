{
gROOT->SetStyle("Plain");
    TCanvas c1;
TChain noHLT("Events"); 
noHLT.Add("14AprnoHLTFilter_MET45_2CentralJet20_2BTagIP_SingleTrackTC_v1_matching.root");


//TChain MET45_2Btag("Events"); 
// MET45_2Btag.Add("2l2bMetEdmNtuples_MET45_v2.root");

TH1D * CSV1 = new TH1D("CSV1", "CSV1", 100, -1, 1);
TH1D * CSV1_2 = new TH1D("CSV1_2", "CSV1_2", 100, -1, 1);



TH1D * CSV2 = new TH1D("CSV2", "CSV2", 100, -1, 1);
TH1D * CSV2_2 = new TH1D("CSV2_2", "CSV2_2", 100, -1, 1);

TH1D *MET = new TH1D("MET", "MET", 100, 50, 250);

noHLT->Project("CSV1", "hjjJet1CSV",  "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 && pfMetEt>120 && hjjMass@.size()==1  &&  hjjJet1HLTBit==1");

noHLT->Project("CSV1_2", "hjjJet2CSV",  "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 && pfMetEt>120 && hjjMass@.size()==1 &&  hjjJet2HLTBit==1");


 CSV1->Add(CSV1_2);



CSV1->SetLineColor(kRed);
CSV1->SetLineWidth(3);
CSV2->SetLineWidth(3);

noHLT->Project("CSV2", "hjjJet1CSV",  "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4  && pfMetEt>120 && hjjMass@.size()==1 && ( hjjJet1HLTBit==0)");

noHLT->Project("CSV2_2", "hjjJet2CSV",  "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4  && pfMetEt>120 && hjjMass@.size()==1 && (hjjJet2HLTBit==0)");

 CSV2->Add(CSV2_2);
CSV2->SetLineColor(kBlue);

CSV1->SetTitle("CSV");
CSV2->SetTitle("CSV");



CSV1->DrawNormalized("");
CSV2->DrawNormalized("same");

 TLegend* leg=new TLegend(0.71,0.60,0.99,1.00);
 leg->SetLineColor(0);
 leg->SetFillColor(0);


 leg->AddEntry(CSV1,"HLT passed", "l");
 leg->AddEntry(CSV2,"HLT failed", "l");
 leg->Draw("same");

 //c1->SaveAs("CSV_HLT_MET45.eps");


 /*

 // now ranking the jets 
 


noHLT->Project("CSV1", "max(hjjJet1CSV, hjjJet2CSV)",  "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 && pfMetEt>120 && hjjMass@.size()==1  &&  ( ((hjjJet1CSV> hjjJet2CSV) && hjjJet1HLTBit==1)  || ((hjjJet2CSV> hjjJet1CSV) && hjjJet2HLTBit==1))");





CSV1->SetLineColor(kRed);
CSV1->SetLineWidth(3);
CSV2->SetLineWidth(3);

noHLT->Project("CSV2", "max(hjjJet1CSV, hjjJet2CSV)",  "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 && pfMetEt>120 && hjjMass@.size()==1  &&  ( ((hjjJet1CSV> hjjJet2CSV) && hjjJet1HLTBit==0)  || ((hjjJet2CSV> hjjJet1CSV) && hjjJet2HLTBit==0))");


CSV2->SetLineColor(kBlue);

CSV1->SetTitle("CSV");
CSV2->SetTitle("CSV");



CSV1->DrawNormalized("");
CSV2->DrawNormalized("same");

 TLegend* leg=new TLegend(0.71,0.60,0.99,1.00);
 leg->SetLineColor(0);
 leg->SetFillColor(0);


 leg->AddEntry(CSV1,"HLT passed", "l");
 leg->AddEntry(CSV2,"HLT failed", "l");
 leg->Draw("same");

 // c1->SaveAs("CSV_HLT_MET45_highRankingJet.eps");




 // low ranking jet

noHLT->Project("CSV1", "min(hjjJet1CSV, hjjJet2CSV)",  "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 && pfMetEt>120 && hjjMass@.size()==1  &&  ( ((hjjJet1CSV> hjjJet2CSV) && hjjJet1HLTBit==1)  || ((hjjJet2CSV> hjjJet1CSV) && hjjJet2HLTBit==1))");





CSV1->SetLineColor(kRed);
CSV1->SetLineWidth(3);
CSV2->SetLineWidth(3);

noHLT->Project("CSV2", "min(hjjJet1CSV, hjjJet2CSV)",  "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 && pfMetEt>120 && hjjMass@.size()==1  &&  ( ((hjjJet1CSV> hjjJet2CSV) && hjjJet1HLTBit==0)  || ((hjjJet2CSV> hjjJet1CSV) && hjjJet2HLTBit==0))");


CSV2->SetLineColor(kBlue);

CSV1->SetTitle("CSV");
CSV2->SetTitle("CSV");




CSV2->DrawNormalized("");
CSV1->DrawNormalized("same");

 TLegend* leg=new TLegend(0.71,0.60,0.99,1.00);
 leg->SetLineColor(0);
 leg->SetFillColor(0);


 leg->AddEntry(CSV1,"HLT passed", "l");
 leg->AddEntry(CSV2,"HLT failed", "l");
 leg->Draw("same");

 // c1->SaveAs("CSV_HLT_MET45_lowRankingJet.eps");



noHLT->Project("MET", "pfMetEt",  "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 &&  hjjMass@.size()==1 &&  hjjJet2HLTBit==1 && hjjJet1HLTBit==1");

 */

}

//TGraphAsymmErrors * eff = new TGraphAsymmErrors(CSV1, CSV2,   "cp");


    
    
    
//eff->GetXaxis()->SetTitle("pfMet");
//eff->GetYaxis()->SetTitle("HLT eff");
 
// eff->SetTitle("");
//       eff->Draw("ap");
