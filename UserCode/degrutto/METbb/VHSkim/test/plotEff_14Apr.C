{

  gROOT->SetStyle("Plain");
    TCanvas c1;
TChain noHLT("Events"); 
noHLT.Add("14Apri_noHltFilter.root");

TChain MET100("Events"); 
MET100.Add("14Apr_MET100.root");

TChain pfMHT150("Events"); 
pfMHT150.Add("14Apr_HLT_pfMHT150_v1.root");

TChain MET80Jet55Jet35("Events"); 
 MET80Jet55Jet35.Add("14Apri_CentralJet55CentralJet35_MET80.root");

TChain MET802Jet20("Events"); 
 MET802Jet20.Add("14Apri_HLT_2CentralJet20_MET80_v1.root");

TChain MET65_1BtagIP("Events"); 
 MET65_1BtagIP.Add("14Apri_HLT_2CentralJet20_BtagIP_MET65_v1.root");


TChain MET65_1BtagTC("Events"); 
 MET65_1BtagTC.Add("14Apri_HLT_2CentralJet20_BtagIP_SingleTrackTC_MET65_v1.root");


TChain MET45_2BTagIP3D_1BTag6IP_SingleTrackTC("Events"); 
 MET45_2BTagIP3D_1BTag6IP_SingleTrackTC.Add("14Apri_HLT_MET45_2CentralJet20_2BTagIP3D_1BTag6IP_SingleTrackTC_v1.root");


TChain MET45_2CentralJet20_2BTagIP_1LooseTC_1Tight("Events"); 
 MET45_2CentralJet20_2BTagIP_1LooseTC_1Tight.Add("14Apri_HLT_MET45_2CentralJet20_2BTagIP_1LooseTC_1Tight_v1.root");

TChain MET45_2CentralJet20_2BTagIP_SingleTrackTC("Events"); 
 MET45_2CentralJet20_2BTagIP_SingleTrackTC.Add("14Apri_HLT_MET45_2CentralJet20_2BTagIP_SingleTrackTC_v1.root");

TChain MET45_2CentralJet20_2BTagIP3D_SingleTrackTC("Events"); 
 MET45_2CentralJet20_2BTagIP3D_SingleTrackTC.Add("14Apri_HLT_MET45_2CentralJet20_2BTagIP3D_SingleTrackTC_v1.root");



TChain MET45_2CentralJet20_1BTagIP3D_SingleTrackTC("Events"); 
 MET45_2CentralJet20_1BTagIP3D_SingleTrackTC.Add("14Apri_HLT_MET45_2CentralJet20_1BTagIP3D_SingleTrackTC_v1.root");

TChain MET45_2CentralJet20_1BTagIP_SingleTrackTC("Events"); 
 MET45_2CentralJet20_1BTagIP_SingleTrackTC.Add("14Apri_HLT_MET45_2CentralJet20_1BTagIP_SingleTrackTC_v1.root");




TChain L1MET30("Events"); 
 L1MET30->Add("14Apri_L1MET30.root");

 TH1D * pfMet = new TH1D("pfMet", "pfMet", 26, -5, 255);
 TH1D * pfMet_met100 = new TH1D("pfMet_met100", "pfMet_met100", 26, -5, 255);
 TH1D * pfMet_mht150 = new TH1D("pfMet_mht150", "pfMet_mht150", 26, -5, 255);

 TH1D * pfMet_met80Jet55Jet35 = new TH1D("pfMet_met80Jet55Jet35", "pfMet_met80Jet55Jet35", 26, -5, 255);
 TH1D * pfMet_met802Jet20 = new TH1D("pfMet_met802Jet20", "pfMet_met802Jet20", 26, -5, 255);

 TH1D * pfMet_met651BtagIP = new TH1D("pfMet_met651BtagIP", "pfMet_met651BtagIP", 26, -5, 255);
 TH1D * pfMet_met651BtagTC = new TH1D("pfMet_met651BtagTC", "pfMet_met651BtagTC", 26, -5, 255);



 TH1D * pfMet_met452Btag3D_1BTag6IP_SingleTrackTC = new TH1D("pfMet_met452Btag3D_1BTag6IP_SingleTrackTC", "pfMet_met452Btag3D_1BTag6IP_SingleTrackTC", 26, -5, 255);
 TH1D * pfMet_met452Btag3D_2BTag_SingleTrackTC = new TH1D("pfMet_met452Btag3D_2BTag_SingleTrackTC", "pfMet_met452Btag3D_2BTag_SingleTrackTC", 26, -5, 255);
 TH1D * pfMet_met452Btag_2BTag_1LooseTC_1Tight = new TH1D("pfMet_met452Btag_2BTag_1LooseTC_1Tight", "pfMet_met452Btag_2BTag_1LooseTC_1Tight", 26, -5, 255);
 TH1D * pfMet_met452Btag_2BTagIP_SingleTrackTC = new TH1D("pfMet_met452Btag_2BTagIP_SingleTrackTC", "pfMet_met452Btag_2BTagIP_SingleTrackTC", 26, -5, 255);
     TH1D * pfMet_met452Btag_2BTagIP3D_SingleTrackTC = new TH1D("pfMet_met452Btag_2BTagIP3D_SingleTrackTC", "pfMet_met452Btag_2BTagIP3D_SingleTrackTC", 26, -5, 255);
 TH1D * pfMet_met452Btag_1BTagIP3D_SingleTrackTC = new TH1D("pfMet_met452Btag_1BTagIP3D_SingleTrackTC", "pfMet_met452Btag_1BTagIP3D_SingleTrackTC", 26, -5, 255);
 TH1D * pfMet_met452Btag_1BTagIP_SingleTrackTC = new TH1D("pfMet_met452Btag_1BTagIP_SingleTrackTC", "pfMet_met452Btag_1BTagIP_SingleTrackTC", 26, -5, 255);




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

MET802Jet20->Project("pfMet_met802Jet20", "pfMetEt", "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 && ((hjjJet1CSV>0.9 && hjjJet2CSV>0.4) || (hjjJet2CSV>0.9 || hjjJet1CSV>0.4))");
 pfMet_met802Jet20->SetLineColor(kOrange-4);
 pfMet_met802Jet20->SetLineWidth(2);



MET65_1BtagIP->Project("pfMet_met651BtagIP", "pfMetEt", "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 && ((hjjJet1CSV>0.9 && hjjJet2CSV>0.4) || (hjjJet2CSV>0.9 || hjjJet1CSV>0.4))");
 pfMet_met651BtagIP->SetLineColor(kOrange-3);
 pfMet_met651BtagIP->SetLineWidth(2);

 //std::cout << "entries in 65 IP" << pfMet_met651BtagIP>GetIntegral()<<std::endl;

MET65_1BtagTC->Project("pfMet_met651BtagTC", "pfMetEt", "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 && ((hjjJet1CSV>0.9 && hjjJet2CSV>0.4) || (hjjJet2CSV>0.9 || hjjJet1CSV>0.4))");
 pfMet_met651BtagTC->SetLineColor(kOrange-5);
 pfMet_met651BtagTC->SetLineWidth(2);


 //std::cout << "entries in 65 TC" << pfMet_met651BtagTC->GetIntegral()<<std::endl;




MET45_2BTagIP3D_1BTag6IP_SingleTrackTC->Project("pfMet_met452Btag3D_1BTag6IP_SingleTrackTC", "pfMetEt", "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 && ((hjjJet1CSV>0.9 && hjjJet2CSV>0.4) || (hjjJet2CSV>0.9 || hjjJet1CSV>0.4))");
 pfMet_met452Btag3D_1BTag6IP_SingleTrackTC->SetLineColor(kYellow+2);
 pfMet_met452Btag3D_1BTag6IP_SingleTrackTC->SetLineWidth(2);



MET45_2CentralJet20_2BTagIP_SingleTrackTC->Project("pfMet_met452Btag_2BTagIP_SingleTrackTC", "pfMetEt", "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 && ((hjjJet1CSV>0.9 && hjjJet2CSV>0.4) || (hjjJet2CSV>0.9 || hjjJet1CSV>0.4))");
 pfMet_met452Btag_2BTagIP_SingleTrackTC->SetLineColor(kBlue);
 pfMet_met452Btag_2BTagIP_SingleTrackTC->SetLineWidth(2);

MET45_2CentralJet20_2BTagIP3D_SingleTrackTC->Project("pfMet_met452Btag_2BTagIP3D_SingleTrackTC", "pfMetEt", "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 && ((hjjJet1CSV>0.9 && hjjJet2CSV>0.4) || (hjjJet2CSV>0.9 || hjjJet1CSV>0.4))");
 pfMet_met452Btag_2BTagIP3D_SingleTrackTC->SetLineColor(kOrange);
 pfMet_met452Btag_2BTagIP3D_SingleTrackTC->SetLineWidth(2);


 MET45_2CentralJet20_2BTagIP_1LooseTC_1Tight->Project("pfMet_met452Btag_2BTag_1LooseTC_1Tight", "pfMetEt", "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 && ((hjjJet1CSV>0.9 && hjjJet2CSV>0.4) || (hjjJet2CSV>0.9 || hjjJet1CSV>0.4))");
 pfMet_met452Btag_2BTag_1LooseTC_1Tight->SetLineColor(kYellow+3);
pfMet_met452Btag_2BTag_1LooseTC_1Tight->SetLineWidth(2);



// MET45_2CentralJet20_2BTagIP3D_SingleTrackTC->Project("pfMet_met452Btag_2BTagIP3D_SingleTrackTC", "pfMetEt", "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 && ((hjjJet1CSV>0.9 && hjjJet2CSV>0.4) || (hjjJet2CSV>0.9 || hjjJet1CSV>0.4))");
// pfMet_met452Btag_2BTagIP3D_SingleTrackTC->SetLineColor(kViolet+3);
// pfMet_met452Btag_2BTagIP3D_SingleTrackTC->SetLineWidth(2);


MET45_2CentralJet20_1BTagIP3D_SingleTrackTC->Project("pfMet_met452Btag_1BTagIP3D_SingleTrackTC", "pfMetEt", "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 && ((hjjJet1CSV>0.9 && hjjJet2CSV>0.4) || (hjjJet2CSV>0.9 || hjjJet1CSV>0.4))");
 pfMet_met452Btag_1BTagIP3D_SingleTrackTC->SetLineColor(kYellow+6);
pfMet_met452Btag_1BTagIP3D_SingleTrackTC->SetLineWidth(2);



MET45_2CentralJet20_1BTagIP_SingleTrackTC->Project("pfMet_met452Btag_1BTagIP_SingleTrackTC", "pfMetEt", "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 && ((hjjJet1CSV>0.9 && hjjJet2CSV>0.4) || (hjjJet2CSV>0.9 || hjjJet1CSV>0.4))");
 pfMet_met452Btag_1BTagIP3D_SingleTrackTC->SetLineColor(kYellow+7);
pfMet_met452Btag_1BTagIP3D_SingleTrackTC->SetLineWidth(2);








L1MET30->Project("pfMet_l1met30", "pfMetEt", "hjjJet1Pt>30 && hjjJet2Pt>30 && abs(hjjJet1Eta)<2.4 && abs(hjjJet2Eta)<2.4 && ((hjjJet1CSV>0.9 && hjjJet2CSV>0.4) || (hjjJet2CSV>0.9 || hjjJet1CSV>0.4))");
 pfMet_l1met30->SetLineColor(kYellow-4);
 pfMet_l1met30->SetLineWidth(2);

  
TGraphAsymmErrors * eff100 = new TGraphAsymmErrors(pfMet_met100, pfMet,   "cp");


    
    
    
 eff100->GetXaxis()->SetTitle("pfMet");
 eff100->GetYaxis()->SetTitle("HLT eff");
 
  eff100->SetTitle("");
       eff100->Draw("ap");
       /*
 TGraphAsymmErrors * eff150 = new TGraphAsymmErrors(pfMet_mht150, pfMet,   "cp");
 eff150->Draw("samep");

 TGraphAsymmErrors * eff805535 = new TGraphAsymmErrors(pfMet_met80Jet55Jet35, pfMet,   "cp");
 eff805535->Draw("samep");

 TGraphAsymmErrors * eff80220 = new TGraphAsymmErrors(pfMet_met802Jet20, pfMet,   "cp");
 eff80220->Draw("samep");


 
 TGraphAsymmErrors * eff65_IP = new TGraphAsymmErrors(pfMet_met651BtagIP, pfMet,   "cp");
 eff65_IP->Draw("samep");
 

 TGraphAsymmErrors * eff65_TC = new TGraphAsymmErrors(pfMet_met651BtagTC, pfMet,   "cp");
 eff65_TC->Draw("samep");
 
 */

 TGraphAsymmErrors * eff4523D_16TC = new TGraphAsymmErrors(pfMet_met452Btag3D_1BTag6IP_SingleTrackTC, pfMet,   "cp");
 eff4523D_16TC->Draw("samep");

 
 // TGraphAsymmErrors * eff452_2TC = new TGraphAsymmErrors(pfMet_met452Btag_2BTag_SingleTrackTC, pfMet,   "cp");
 //eff452_2TC->Draw("samep");
 
 TGraphAsymmErrors * eff4523D_1TCL1TCT = new TGraphAsymmErrors(pfMet_met452Btag_2BTag_1LooseTC_1Tight, pfMet,"cp");
 eff4523D_1TCL1TCT->Draw("samep");

TGraphAsymmErrors * eff452_2TC = new TGraphAsymmErrors(pfMet_met452Btag_2BTagIP_SingleTrackTC, pfMet, "cp");
eff452_2TC->Draw("samep");

TGraphAsymmErrors * eff4523D_2TC = new TGraphAsymmErrors(pfMet_met452Btag_2BTagIP3D_SingleTrackTC, pfMet, "cp");
eff4523D_2TC->Draw("samep");


TGraphAsymmErrors * eff4513D_TC = new TGraphAsymmErrors(pfMet_met452Btag_1BTagIP3D_SingleTrackTC, pfMet,  "cp");
eff4513D_TC->Draw("samep");



 TGraphAsymmErrors * eff451_TC = new TGraphAsymmErrors(pfMet_met452Btag_1BTagIP_SingleTrackTC, pfMet,  "cp");
eff451_TC->Draw("samep");



 TGraphAsymmErrors * effl130 = new TGraphAsymmErrors(pfMet_l1met30, pfMet,   "cp");
 effl130->Draw("samep");



 TLegend* leg=new TLegend(0.105,0.50,0.34,0.93);
 // leg->SetLineColor(0);
 leg->SetFillColor(0);


 leg->AddEntry(eff100,"CaloMET100", "l");
   /*
 leg->AddEntry(eff150,"pfMHT150", "l");
 leg->AddEntry(eff805535,"CaloMET80_cJ55cJ35", "l");
 leg->AddEntry(eff80220,"CaloMET80_2cJ20", "l");
  leg->AddEntry(eff65_IP,"CaloMET65_2cJ20_1BTag2ndTC", "l");
   leg->AddEntry(eff65_TC,"CaloMET65_2cJ20_1BTag1stTC", "l");
   */

 leg->AddEntry(eff4523D_16TC,"CaloMET45_2cJ20_2BTag3D1stTC_one6IP", "l");


 leg->AddEntry(eff4523D_1TCL1TCT,"CaloMET45_2cJ20_2BTag3D_1stL2ndT", "l");

  leg->AddEntry(eff452_2TC,"CaloMET45_2cJ20_2BTag1stTC", "l");

 

  leg->AddEntry(eff4523D_2TC,"CaloMET45_2cJ20_2BTag3D1stTC", "l");

  leg->AddEntry(eff4513D_TC,"CaloMET45_2cJ20_1BTag3D1stTC", "l");

  leg->AddEntry(eff451_TC,"CaloMET45_2cJ20_1BTag1stTC", "l");
   

 leg->AddEntry(effl130,"L1CaloMET30", "l");
 leg->Draw("same");
 
 c1.SaveAs("HLTEff_vs_pfMet_14apr_metl65.gif");
 c1.SaveAs("HLTEff_vs_pfMet_14apr_metl65.eps");
}

