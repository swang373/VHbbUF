{
  gROOT->SetStyle("Plain");
    TCanvas c1;
TChain noHLT("Events"); 
noHLT.Add("2l2bMetEdmNtuples_noTrig.root");

TChain MET100("Events"); 
MET100.Add("2l2bMetEdmNtuples_MET100_v2.root");

TChain pfMHT150("Events"); 
pfMHT150.Add("2l2bMetEdmNtuples_PFMHT150_v2.root");
TChain MET80Jet55Jet35("Events"); 
 MET80Jet55Jet35.Add("2l2bMetEdmNtuples_MET80_v2.root");

TChain MET65_1Btag("Events"); 
 MET65_1Btag.Add("2l2bMetEdmNtuples_MET65_v2.root");

TChain MET45_2Btag("Events"); 
 MET45_2Btag.Add("2l2bMetEdmNtuples_MET45_v2.root");

TChain L1MET30("Events"); 
 L1MET30->Add("2l2bMetEdmNtuples_L1ETM30_L1Filter.root");

 TH1D * pfMet = new TH1D("pfMet", "pfMet", 26, -5, 255);
 TH1D * pfMet_met100 = new TH1D("pfMet_met100", "pfMet_met100", 26, -5, 255);
 TH1D * pfMet_mht150 = new TH1D("pfMet_mht150", "pfMet_mht150", 26, -5, 255);
 TH1D * pfMet_met80Jet55Jet35 = new TH1D("pfMet_met80Jet55Jet35", "pfMet_met80Jet55Jet35", 26, -5, 255);
 TH1D * pfMet_met651Btag = new TH1D("pfMet_met651Btag", "pfMet_met651Btag", 26, -5, 255);
 TH1D * pfMet_met452Btag = new TH1D("pfMet_met452Btag", "pfMet_met452Btag", 26, -5, 255);
 TH1D * pfMet_l1met30 = new TH1D("pfMet_l1met30", "pfMet_l1met30", 26, -5, 255);



 noHLT->Project("pfMet", "pfMetEt", "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 && ((hjjJet1CSV>0.9 && hjjJet2CSV>0.4) || (hjjJet2CSV>0.9 || hjjJet1CSV>0.4))");
 MET100->Project("pfMet_met100", "pfMetEt", "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 && ((hjjJet1CSV>0.9 && hjjJet2CSV>0.4) || (hjjJet2CSV>0.9 || hjjJet1CSV>0.4))");
 pfMet_met100->SetLineColor(kRed);
 pfMet_met100->SetLineWidth(2);

 pfMHT150->Project("pfMet_mht150", "pfMetEt", "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 && ((hjjJet1CSV>0.9 && hjjJet2CSV>0.4) || (hjjJet2CSV>0.9 || hjjJet1CSV>0.4))");


 pfMet_mht150->SetLineColor(kBlue);
 pfMet_mht150->SetLineWidth(2);

MET80Jet55Jet35->Project("pfMet_met80Jet55Jet35", "pfMetEt", "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 && ((hjjJet1CSV>0.9 && hjjJet2CSV>0.4) || (hjjJet2CSV>0.9 || hjjJet1CSV>0.4))");
 pfMet_met80Jet55Jet35->SetLineColor(kOrange-1);
 pfMet_met80Jet55Jet35->SetLineWidth(2);


MET65_1Btag->Project("pfMet_met651Btag", "pfMetEt", "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 && ((hjjJet1CSV>0.9 && hjjJet2CSV>0.4) || (hjjJet2CSV>0.9 || hjjJet1CSV>0.4))");
 pfMet_met651Btag->SetLineColor(kOrange-3);
 pfMet_met651Btag->SetLineWidth(2);

MET45_2Btag->Project("pfMet_met452Btag", "pfMetEt", "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 && ((hjjJet1CSV>0.9 && hjjJet2CSV>0.4) || (hjjJet2CSV>0.9 || hjjJet1CSV>0.4))");
 pfMet_met452Btag->SetLineColor(kYellow+2);
 pfMet_met452Btag->SetLineWidth(2);


L1MET30->Project("pfMet_l1met30", "pfMetEt", "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 && ((hjjJet1CSV>0.9 && hjjJet2CSV>0.4) || (hjjJet2CSV>0.9 || hjjJet1CSV>0.4))");
 pfMet_l1met30->SetLineColor(kYellow-4);
 pfMet_l1met30->SetLineWidth(2);

  
TGraphAsymmErrors * eff = new TGraphAsymmErrors(pfMet_met100, pfMet,   "cp");


    
    
    
 eff->GetXaxis()->SetTitle("pfMet");
 eff->GetYaxis()->SetTitle("HLT eff");
 
  eff->SetTitle("");
       eff->Draw("ap");

 TGraphAsymmErrors * eff2 = new TGraphAsymmErrors(pfMet_mht150, pfMet,   "cp");
 eff2->Draw("samep");

 TGraphAsymmErrors * eff3 = new TGraphAsymmErrors(pfMet_met80Jet55Jet35, pfMet,   "cp");
 eff3->Draw("samep");

 TGraphAsymmErrors * eff4 = new TGraphAsymmErrors(pfMet_met651Btag, pfMet,   "cp");
 eff4->Draw("samep");

 TGraphAsymmErrors * eff5 = new TGraphAsymmErrors(pfMet_met452Btag, pfMet,   "cp");
 eff5->Draw("samep");

 TGraphAsymmErrors * eff6 = new TGraphAsymmErrors(pfMet_l1met30, pfMet,   "cp");
 eff6->Draw("samep");

 TLegend* leg=new TLegend(0.105,0.50,0.34,0.93);
 // leg->SetLineColor(0);
 leg->SetFillColor(0);


 leg->AddEntry(eff,"CaloMET100", "l");
 leg->AddEntry(eff2,"pfMHT150", "l");
 leg->AddEntry(eff3,"CaloMET80_CeJ55CeJ35", "l");
 leg->AddEntry(eff4,"CaloMET65_2CeJ20_1BT", "l");
 leg->AddEntry(eff5,"CaloMET45_2CeJ20_2BT_1TC", "l");
 leg->AddEntry(eff6,"L1CaloMET30", "l");
 leg->Draw("same");
 
 c1.SaveAs("HLTEff_vs_pfMet.gif");
 c1.SaveAs("HLTEff_vs_pfMet.eps");
}
