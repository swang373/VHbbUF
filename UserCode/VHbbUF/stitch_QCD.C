{
chain = new TChain("tree");
TString fname = "Step3_20130130/";
//chain->Add(fname+"Step3_QCDPt80.root");
chain->Add(fname+"Step3_QCDPt120.root");
chain->Add(fname+"Step3_QCDPt170.root");
chain->Add(fname+"Step3_QCDPt300.root");
chain->Add(fname+"Step3_QCDPt470.root");
chain->Add(fname+"Step3_QCDPt600.root");
chain->Add(fname+"Step3_QCDPt800.root");
chain->Add(fname+"Step3_QCDPt1000.root");
chain->Add(fname+"Step3_QCDPt1400.root");
chain->Add(fname+"Step3_QCDPt1800.root");
chain->Draw("hJet_pt[0] >> h0(250,0.,500.)", "efflumi","goff");
h0->Draw();
chain->Draw("hJet_genPt[0] >> h1(250,0.,500.)", "efflumi","goff");
h1->SetLineColor(kRed);
h1->Draw("same");
gPad->SetLogy();
gPad->Print("stitch_QCD.png");
}
