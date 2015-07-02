#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <vector>

#include "TSystem.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "THStack.h"
#include "TCut.h"
#include "TString.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TMath.h"

//#include "XSec_8TeV19invfb.h"
#include "XSec_13TeV20invfb.h"

#include "TLorentzVector.h"





 TCut cuttest = "Vtype==4 &&  Jet_btagCSV[hJCidx[0]]>0.814 && Jet_btagCSV[hJCidx[1]]>0.423   && HCSV_mass<250 &&  met_pt>170 && HCSV_pt>130 && ( (Jet_pt[hJCidx[0]]>80 && Jet_pt[hJCidx[1]]>30) || (Jet_pt[hJCidx[1]]>80 && Jet_pt[hJCidx[0]]>30 ))    &&  Sum$(aLeptons_pt>5  && aLeptons_relIso03<1. &&  ( deltaR(Jet_eta[hJCidx[0]], aLeptons_eta, Jet_phi[hJCidx[0]], aLeptons_phi)>0.4 &&  deltaR( Jet_eta[hJCidx[1]], aLeptons_eta, Jet_phi[hJCidx[1]], aLeptons_phi)>0.4))==0 && Sum$(Jet_pt>30 && abs(Jet_eta)<4.5 && Jet_puId==1)<=99 && MinIf$(abs(deltaPhi(met_phi,Jet_phi)) , Jet_pt>30 && abs(Jet_eta)<4.5 && Jet_puId==1)>0.5 && abs(deltaPhi(HCSV_phi, met_phi))>2.5";





// _____________________________________________________________________________
// This class reads the ntuples and applies basic cuts to reduce loading time.
class Events {
public:
    Events();
    ~Events();

    void read(TCut cutmc_all, TCut cutdata_all, TString processes="", TString treename = "test");

    void set_sf(const double sf[5]) {
        sf_WjLF = sf[0];
        sf_WjHF = sf[1];
        sf_ZjLF = sf[2];
        sf_ZjHF = sf[3];
        sf_TT   = sf[4];
        return;
    }

    // Reduced TTrees
    TTree * ZH;
    TTree * monoH;
    TTree * WH;
    TTree * WjLF;
    TTree * WjHF;
    TTree * ZjLF;
    TTree * ZjHF;

    TTree * ZjLF_HT100;
    TTree * ZjHF_HT100;
    TTree * ZjLF_HT200;
    TTree * ZjHF_HT200;
    TTree * ZjLF_HT600;
    TTree * ZjHF_HT600;


    TTree * TT;
    TTree * ST_T_s;
    TTree * ST_T_t;
    TTree * ST_T_tW;
    TTree * ST_Tbar_s;
    TTree * ST_Tbar_t;
    TTree * ST_Tbar_tW;
    TTree * s_Top;
    TTree * VV_WW;
    TTree * VV_WZ;
    TTree * VV_ZZ;
    TTree * QCD;
    TTree * data_obs;

    // Equivalent luminosities of processes
    float lumi_ZH;
    float lumi_monoH;
    float lumi_WH;
    float lumi_WjLF;
    float lumi_WjHF;
    float lumi_ZjLF;
    float lumi_ZjHF;

    float lumi_ZjLFHT100;
    float lumi_ZjHFHT100;
    float lumi_ZjLFHT200;
    float lumi_ZjHFHT200;
    float lumi_ZjLFHT600;
    float lumi_ZjHFHT600;

    float lumi_TT;
    float lumi_ST_T_s;
    float lumi_ST_T_t;
    float lumi_ST_T_tW;
    float lumi_ST_Tbar_s;
    float lumi_ST_Tbar_t;
    float lumi_ST_Tbar_tW;
    float lumi_s_Top;
    float lumi_VV_WW;
    float lumi_VV_WZ;
    float lumi_VV_ZZ;

    float lumi_QCD;

    // Scale factors
    float sf_WjLF;
    float sf_WjHF;
    float sf_ZjLF;
    float sf_ZjHF;
    float sf_TT;

    int massH;  // Higgs mass

};  // Events

// _____________________________________________________________________________
// These are functions to set the styles of TH1 and TLegend, respectively
void set_style(TH1 * h, const TString& p);
void set_style_norm(TH1 * h, const TString& p);
void set_style(TLegend * leg);

// _____________________________________________________________________________
// This function makes a single plot
void MakePlot(TTree * tree, TString var, TCut cut,
              TString title, int nbinsx, double xlow, double xup,
              TString plotname="plot", TString plotdir="",
              TString options="plotLog:plotNorm") {

    std::clog << "MakePlots(): Plot var: " << var << std::endl;
    std::clog << "MakePlots(): Using cut: " << cut << std::endl;

    bool plotLog      = options.Contains("plotLog")    && (!options.Contains("!plotLog"));
    bool plotNorm     = options.Contains("plotNorm")   && (!options.Contains("!plotNorm"));

    TH1F * h1         = new TH1F("h1"       , title, nbinsx, xlow, xup);
    tree->Project("h1", var, cut);

    TCanvas * c1 = new TCanvas("c1", "c1", 700, 700);
    c1->SetLogy(plotLog);
    h1->SetLineWidth  (2);
    h1->SetLineColor  (2);
    h1->SetMarkerSize (0);
    h1->GetYaxis()->SetTitle("Events");
    if (plotNorm) {
        h1->Scale(1.0 / h1->Integral());
        h1->GetYaxis()->SetTitle("arbitrary unit");
    }

    h1->Draw("hist");
    h1->Draw("same e1");

    // Setup legends
    //    TLegend * leg1 = new TLegend(0.50, 0.68, 0.72, 0.92);
    //set_style(leg1);
    //leg1->AddEntry(h1, tree->GetUserInfo()->At(0)->GetName() , "l");
    //    leg1->Draw();

    TLatex * latex = new TLatex();
    latex->SetNDC();
    latex->SetTextAlign(12);
    latex->SetTextFont(62);
    latex->SetTextSize(0.04);
    //latex->DrawLatex(0.19, 0.89, "CMS Preliminary");
    latex->DrawLatex(0.19, 0.89, "CMS Simulation 2015");
    latex->SetTextSize(0.04);
    latex->DrawLatex(0.19, 0.84, "#sqrt{s} = 13 TeV");


    gPad->Print(plotdir+plotname+".png");
    gPad->Print(plotdir+plotname+".pdf");

    // If not deleted, they will remain in the memory
    //delete c1;
    //delete h1;
    return;
}

// _____________________________________________________________________________
// This function makes two plots overlaid
void MakePlot2(TTree * tree1, TTree * tree2, TString var, TCut cut,
               TString title, int nbinsx, double xlow, double xup,
               TString plotname="plot", TString plotdir="",
               TString options="plotLog:plotNorm") {

    std::clog << "MakePlots(): Plot var: " << var << std::endl;
    std::clog << "MakePlots(): Using cut: " << cut << std::endl;

    bool plotLog      = options.Contains("plotLog")    && (!options.Contains("!plotLog"));
    bool plotNorm     = options.Contains("plotNorm")   && (!options.Contains("!plotNorm"));

    TH1F * h1         = new TH1F("h1"       , title, nbinsx, xlow, xup);
    TH1F * h2         = new TH1F("h2"       , title, nbinsx, xlow, xup);
    tree1->Project("h1", var, cut);
    tree2->Project("h2", var, cut);

    TCanvas * c1 = new TCanvas("c1", "c1", 700, 700);
    c1->SetLogy(plotLog);
    h1->SetLineWidth  (2);
    h1->SetLineColor  (2);
    h1->SetMarkerSize (0);
    h1->GetYaxis()->SetTitle("Events");
    h2->SetLineWidth  (2);
    h2->SetLineColor  (4);
    h2->SetMarkerSize (0);
    h2->GetYaxis()->SetTitle("Events");
    if (plotNorm) {
        h1->Scale(1.0 / h1->Integral());
        h1->GetYaxis()->SetTitle("arbitrary unit");
        h2->Scale(1.0 / h2->Integral());
        h2->GetYaxis()->SetTitle("arbitrary unit");
    }

    h1->Draw("hist");
    h1->Draw("same e1");
    h2->Draw("same hist");
    h2->Draw("same e1");


    // Setup legends
    //TLegend * leg1 = new TLegend(0.50, 0.68, 0.72, 0.92);
    //set_style(leg1);
    //leg1->AddEntry(h1, tree1->GetUserInfo()->At(0)->GetName() , "l");
    //leg1->AddEntry(h2, tree2->GetUserInfo()->At(0)->GetName() , "l");
    //    leg1->Draw();
    

    TLatex * latex = new TLatex();
    latex->SetNDC();
    latex->SetTextAlign(12);
    latex->SetTextFont(62);
    latex->SetTextSize(0.04);
    //latex->DrawLatex(0.19, 0.89, "CMS Preliminary");
    latex->DrawLatex(0.19, 0.89, "CMS Simulation 2015");
    latex->SetTextSize(0.04);
    latex->DrawLatex(0.19, 0.84, "#sqrt{s} = 13 TeV");


    gPad->Print(plotdir+plotname+".png");
    gPad->Print(plotdir+plotname+".pdf");



    // If not deleted, they will remain in the memory
    //delete c1;
    //delete h1;
    //delete h2;
    return;
}

// _____________________________________________________________________________
// This function makes stacked plots
void MakePlots(const Events * ev, TString var,
               TCut cutmc, TCut cutdata,
               TString title, int nbinsx, double xlow, double xup,
               TString plotname="plot", TString plotdir="",
               TString options="printStat:plotSig:plotData:plotLog:!plotNorm:!plotMonoH") {

    std::clog << "MakePlots(): Plot var: " << var << std::endl;
    std::clog << "MakePlots(): Using cutmc: " << cutmc << ", cutdata: " << cutdata << std::endl;

    // Parse options
    bool printStat    = options.Contains("printStat")  && (!options.Contains("!printStat"));
    //bool printCard    = options.Contains("printCard")  && (!options.Contains("!printCard"));
    //bool plotUOflow   = options.Contains("plotUOflow") && (!options.Contains("!plotUOflow"));  // FIXME: not yet implemented
    bool plotSig      = options.Contains("plotSig")    && (!options.Contains("!plotSig"));
    bool plotMonoH      = options.Contains("plotMonoH")    && (!options.Contains("!plotMonoH"));
    bool plotData     = options.Contains("plotData")   && (!options.Contains("!plotData"));
    bool plotLog      = options.Contains("plotLog")    && (!options.Contains("!plotLog"));

    bool plotNorm     = options.Contains("plotNorm")   && (!options.Contains("!plotNorm"));

    // Book histograms
    TH1F * hZH        = new TH1F("ZH"       , title, nbinsx, xlow, xup);
    TH1F * hmonoH        = new TH1F("monoH"       , title, nbinsx, xlow, xup);
    TH1F * hWH        = new TH1F("WH"       , title, nbinsx, xlow, xup);
    TH1F * hWjLF      = new TH1F("WjLF"     , title, nbinsx, xlow, xup);
    TH1F * hWjHF      = new TH1F("WjHF"     , title, nbinsx, xlow, xup);
    TH1F * hZjLF      = new TH1F("ZjLF"     , title, nbinsx, xlow, xup);
    TH1F * hZjHF      = new TH1F("ZjHF"     , title, nbinsx, xlow, xup);
    TH1F * hTT        = new TH1F("TT"       , title, nbinsx, xlow, xup);
    TH1F * hST        = new TH1F("ST"       , title, nbinsx, xlow, xup);
    TH1F * hVV        = new TH1F("VV"       , title, nbinsx, xlow, xup);
    TH1F * hQCD        = new TH1F("QCD"       , title, nbinsx, xlow, xup);
    TH1F * hdata_obs  = new TH1F("data_obs" , title, nbinsx, xlow, xup);
    TH1F * hVH        = new TH1F("VH"       , title, nbinsx, xlow, xup);
    TH1F * hmc_exp    = new TH1F("mc_exp"   , title, nbinsx, xlow, xup);


    // Draw histograms, apply cut and MC event weights (= equiv_lumi * PUweight)  // no PUweight for now...
    if (plotSig) {
      //        ev->ZH->Project("ZH", var, cutmc * Form("%5f * 1.", ev->lumi_ZH));
        ev->ZH->Project("ZH", var, cutmc );
        std::clog << "... DONE: project ZH." << std::endl;

	//        ev->WH->Project("WH", var, cutmc * Form("%5f * 1.", ev->lumi_WH));
        ev->WH->Project("WH", var, cutmc );
        std::clog << "... DONE: project WH." << std::endl;
    }

    if (plotMonoH) {
	//        ev->monoH->Project("monoH", var, cutmc * Form("%5f * 1.", ev->lumi_monoH));
        ev->monoH->Project("monoH", var, cutmc );
        std::clog << "... DONE: project monoH." << std::endl;
    }



    //    ev->WjLF->Project("WjLF", var, cutmc * Form("%5f * 1.", ev->lumi_WjLF) * Form("%f", ev->sf_WjLF));
    ev->WjLF->Project("WjLF", var, cutmc  * Form("%f", ev->sf_WjLF));
    std::clog << "... DONE: project WjLF, scaled by " << ev->sf_WjLF << "." << std::endl;

    //    ev->WjLF->Project("WjHF", var, cutmc * Form("%5f * 1", ev->lumi_WjHF) * Form("%f", ev->sf_WjHF));
    ev->WjHF->Project("WjHF", var, cutmc  * Form("%f", ev->sf_WjHF));
    std::clog << "... DONE: project WjHF, scaled by " << ev->sf_WjHF << "." << std::endl;

    //    ev->ZjLF->Project("ZjLF", var, cutmc * Form("%5f * 1", ev->lumi_ZjLF) * Form("%f", ev->sf_ZjLF));
       ev->ZjLF->Project("ZjLF", var, cutmc  * Form("%f", ev->sf_ZjLF));
     std::clog << "... DONE: project ZjLF, scaled by " << ev->sf_ZjLF << "." << std::endl;

     //    ev->ZjHF->Project("ZjHF", var, cutmc * Form("%5f * 1", ev->lumi_ZjHF) * Form("%f", ev->sf_ZjHF));
    ev->ZjHF->Project("ZjHF", var, cutmc  * Form("%f", ev->sf_ZjHF));
    std::clog << "... DONE: project ZjHF, scaled by " << ev->sf_ZjHF << "." << std::endl;    
    
    /*
    ev->ZjLF_HT100->Project("ZjLF", var, cutmc * Form("%5f * 1", ev->lumi_ZjLFHT100) * Form("%f", ev->sf_ZjLF)); 
    ev->ZjLF_HT200->Project("+ZjLF", var, cutmc * Form("%5f * 1", ev->lumi_ZjLFHT200) * Form("%f", ev->sf_ZjLF)); 
    ev->ZjLF_HT600->Project("+ZjLF", var, cutmc * Form("%5f * 1", ev->lumi_ZjLFHT600) * Form("%f", ev->sf_ZjLF)); 


    ev->ZjHF_HT100->Project("ZjHF", var, cutmc * Form("%5f * 1", ev->lumi_ZjHFHT100) * Form("%f", ev->sf_ZjHF)); 
    ev->ZjHF_HT200->Project("+ZjHF", var, cutmc * Form("%5f * 1", ev->lumi_ZjHFHT200) * Form("%f", ev->sf_ZjHF)); 
    ev->ZjHF_HT600->Project("+ZjHF", var, cutmc * Form("%5f * 1", ev->lumi_ZjHFHT600) * Form("%f", ev->sf_ZjHF)); 
    */


    //    ev->TT->Project("TT", var, cutmc * Form("%5f * 1.", ev->lumi_TT) * Form("%f", ev->sf_TT));
    ev->TT->Project("TT", var, cutmc  * Form("%f", ev->sf_TT));
    std::clog << "... DONE: project TT, scaled by " << ev->sf_TT << "." << std::endl;


    /*
    //    ev->ST_T_s    ->Project( "ST" , var, cutmc * Form("%5f * 1.", ev->lumi_ST_T_s));
    ev->ST_T_s    ->Project( "ST" , var, cutmc );

    //    ev->ST_T_t    ->Project("+ST" , var, cutmc * Form("%5f * 1.", ev->lumi_ST_T_t));
    ev->ST_T_t    ->Project("+ST" , var, cutmc );

    //    ev->ST_T_tW   ->Project("+ST" , var, cutmc * Form("%5f * 1.", ev->lumi_ST_T_tW));
    ev->ST_T_tW   ->Project("+ST" , var, cutmc );

    //    ev->ST_Tbar_s ->Project("+ST" , var, cutmc * Form("%5f * 1.", ev->lumi_ST_Tbar_s));
    ev->ST_Tbar_s ->Project("+ST" , var, cutmc );


    //    ev->ST_Tbar_t ->Project("+ST" , var, cutmc * Form("%5f * 1.", ev->lumi_ST_Tbar_t));
    ev->ST_Tbar_t ->Project("+ST" , var, cutmc );

    //    ev->ST_Tbar_tW->Project("+ST" , var, cutmc * Form("%5f * 1.", ev->lumi_ST_Tbar_tW));
    ev->ST_Tbar_tW->Project("+ST" , var, cutmc );
    */
    ev->s_Top->Project("ST" , var, cutmc );                                                                                                                                    



    std::clog << "... DONE: project ST." << std::endl;

    /*
    ev->VV_WW->Project( "VV", var, cutmc * Form("%5f * 1.", ev->lumi_VV_WW));
    ev->VV_WZ->Project("+VV", var, cutmc * Form("%5f * 1.", ev->lumi_VV_WZ));
    ev->VV_ZZ->Project("+VV", var, cutmc * Form("%5f * 1.", ev->lumi_VV_ZZ));
    */

    std::clog << "... DONE: project VV." << std::endl;


    ev->QCD ->Project("QCD" , var, cutmc );
    std::clog << "... DONE: project QCD." << std::endl;


    if (plotData) {
        ev->data_obs->Project("data_obs", var, cutdata);
    }
    std::clog << "... DONE: project data_obs." << std::endl;

    // Sum histograms
    hVH->Add(hZH);
    hVH->Add(hWH);
    std::clog << "... DONE: add ZH, WH to VH." << std::endl;
 


    if (plotNorm) {
       
        TCanvas * c1 = new TCanvas("c1", "c1", 700, 700);
	c1->SetLogy(plotLog);


        hVH->Scale(1.0 / hVH->Integral());
        hVH->GetYaxis()->SetTitle("arbitrary unit");
        hmonoH->Scale(1.0 / hmonoH->Integral());
        hmonoH->GetYaxis()->SetTitle("arbitrary unit");
        hST->Scale(1.0 / hST->Integral());
        hST->GetYaxis()->SetTitle("arbitrary unit");
        hTT->Scale(1.0 / hTT->Integral());
        hTT->GetYaxis()->SetTitle("arbitrary unit");
        hZjLF->Scale(1.0 / hZjLF->Integral());
        hZjLF->GetYaxis()->SetTitle("arbitrary unit");
        hZjHF->Scale(1.0 / hZjHF->Integral());
        hZjHF->GetYaxis()->SetTitle("arbitrary unit");
        hWjLF->Scale(1.0 / hWjLF->Integral());
        hWjLF->GetYaxis()->SetTitle("arbitrary unit");
        hWjHF->Scale(1.0 / hWjHF->Integral());
        hWjHF->GetYaxis()->SetTitle("arbitrary unit");
        hVV->Scale(1.0 / hVV->Integral());
        hVV->GetYaxis()->SetTitle("arbitrary unit");
        hQCD->Scale(1.0 / hQCD->Integral());
        hQCD->GetYaxis()->SetTitle("arbitrary unit");

    set_style_norm(hVH, "VH");
    set_style_norm(hZH, "VH");
    set_style_norm(hmonoH, "monoH");
    set_style_norm(hWH, "VH");
    set_style_norm(hWjLF, "WjLF");
    set_style_norm(hWjHF, "WjHFb");
    set_style_norm(hZjLF, "ZjLF");
    set_style_norm(hZjHF, "ZjHFb");
    set_style_norm(hTT, "TT");
    set_style_norm(hST, "ST");
    set_style_norm(hVV, "VVHF");
    set_style_norm(hQCD, "QCD");



    
	hVH->Draw("hist");
	hVH->Draw("same e1");
	hmonoH->Draw("same hist");

	hST->Draw("same hist");
	hST->Draw("same e1");
	hTT->Draw("same hist");
	hTT->Draw("same e1");
	hZjLF->Draw("same hist");
	hZjLF->Draw("same e1");
	hZjHF->Draw("same hist");
	hZjHF->Draw("same e1");
	hWjLF->Draw("same hist");
	hWjLF->Draw("same e1");
	hWjHF->Draw("same hist");
	hWjHF->Draw("same e1");
	hVV->Draw("same hist");
	hVV->Draw("same e1");
	hQCD->Draw("same hist");
	hQCD->Draw("same e1");

	// Setup legends
	TLegend * leg1 = new TLegend(0.50, 0.68, 0.72, 0.92);
	set_style(leg1);
	if (plotData) leg1->AddEntry(hdata_obs, "Data", "p");
	if (plotSig)  leg1->AddEntry(hVH, Form("VH(%i)", ev->massH), "l");
	if (plotMonoH)  leg1->AddEntry(hmonoH, Form("monoH(%i)", ev->massH), "l");
	leg1->AddEntry(hTT, "t#bar{t}", "f");
	leg1->AddEntry(hST, "single top", "f");
	leg1->AddEntry(hVV, "VV", "f");

	TLegend * leg2 = new TLegend(0.72, 0.68, 0.94, 0.92);
	set_style(leg2);
	leg2->AddEntry(hWjHF, "W + HF", "f");
	leg2->AddEntry(hWjLF, "W + LF", "f");
	leg2->AddEntry(hZjHF, "Z + HF", "f");
	leg2->AddEntry(hZjLF, "Z + LF", "f");
	leg2->AddEntry(hQCD, "QCD", "f");


        leg1->Draw("same");
        leg2->Draw("same");
 
	TLatex * latex = new TLatex();
	latex->SetNDC();
	latex->SetTextAlign(12);
	latex->SetTextFont(62);
	latex->SetTextSize(0.04);
	//latex->DrawLatex(0.19, 0.89, "CMS Preliminary");
	latex->DrawLatex(0.19, 0.89, "CMS Simulation 2015");
	latex->SetTextSize(0.04);
	latex->DrawLatex(0.19, 0.84, "#sqrt{s} = 13 TeV");
	
	
	gPad->Print(plotdir+plotname+".png");
	gPad->Print(plotdir+plotname+".pdf");



    // If not deleted, they will remain in the memory
    //delete c1;
    //delete h1;
    //delete h2;
	return;

   
    }



       


    hmc_exp->Add(hWjLF);
    hmc_exp->Add(hWjHF);
    hmc_exp->Add(hZjLF);
    hmc_exp->Add(hZjHF);
    hmc_exp->Add(hTT);
    hmc_exp->Add(hST);
    hmc_exp->Add(hVV);
    hmc_exp->Add(hQCD);
    std::clog << "... DONE: add all backgrounds to mc_exp." << std::endl;


    // Setting up histograms (boring stuff below)
    std::clog << "MakePlots(): Setting up histograms..." << std::endl;

    // Setup histogram stack
    THStack * hs = new THStack("hs", "");
    hs->Add(hST);
    hs->Add(hTT);
    hs->Add(hZjLF);
    hs->Add(hZjHF);
    hs->Add(hWjLF);
    hs->Add(hWjHF);
    hs->Add(hVV);
    hs->Add(hQCD);
    if (plotSig)  hs->Add(hVH);
    if (plotMonoH)  hs->Add(hmonoH);




    if (printStat) {
        double nS = hVH->Integral();
      double nS_err = -1; 
      if (plotMonoH) {
	 nS = hmonoH->IntegralAndError(0,hmonoH->GetNbinsX(), nS_err);
      } else {
	nS = hVH->IntegralAndError(0,hVH->GetNbinsX(), nS_err);
      }
      double nVH_err = -1;
      double nVH = hVH->IntegralAndError(0,hVH->GetNbinsX(), nVH_err);
        //double nB = ((TH1*) hs->GetStack()->At(hs->GetHists()->GetSize()-1))->Integral();  //< integral of all added histograms in the stack
      double nB_err =-1; 
      double nB = hmc_exp->IntegralAndError(0, hmc_exp->GetNbinsX(), nB_err);
        
      double nZjHF_err =-1; 
      double nZjHF = hZjHF->IntegralAndError(0, hZjHF->GetNbinsX(), nZjHF_err);

      double nZjLF_err =-1; 
      double nZjLF = hZjLF->IntegralAndError(0, hZjLF->GetNbinsX(), nZjLF_err);
      double nWjLF_err =-1; 
      double nWjLF = hWjLF->IntegralAndError(0, hWjLF->GetNbinsX(), nWjLF_err);

      double nWjHF_err =-1;  
      double nWjHF = hWjHF->IntegralAndError( 0, hWjHF->GetNbinsX(), nWjHF_err);
      double nTT_err =-1;
      double nTT = hTT->IntegralAndError(0, hTT->GetNbinsX(), nTT_err);

      double nST_err =-1;
      double nST = hST->IntegralAndError(0, hST->GetNbinsX(), nST_err);

        double signif_simple = nS / sqrt(nS + nB);  // only for reference
        double signif_punzi = nS / (3.0/2.0 + sqrt(nB) + 0.2*nB);
        std::cout << "MakePlots(): Printing FOM..." << std::endl;
        std::cout << Form("S                      = %.3f +- %.3f", nS, nS_err) << std::endl;
        std::cout << Form("B                      = %.3f +- %.3f", nB, nB_err) << std::endl;
        std::cout << Form("S/sqrt(S+B)            = %.3f ", signif_simple) << std::endl;
        std::cout << Form("S/(sqrt(B)+a/2+ 0.2*B) = %.3f", signif_punzi) << std::endl;
        std::cout << Form("VH                      = %.3f +- %.3f", nVH, nVH_err) << std::endl;
        std::cout << Form("ZjHF                      = %.3f +- %.3f", nZjHF, nZjHF_err) << std::endl;
        std::cout << Form("WjHF                      = %.3f +- %.3f", nWjHF, nWjHF_err) << std::endl;
        std::cout << Form("ZjLF                      = %.3f +- %.3f", nZjLF, nZjLF_err) << std::endl;
        std::cout << Form("WjLF                      = %.3f +- %.3f", nWjLF, nWjLF_err) << std::endl;
        std::cout << Form("TT                     = %.3f +- %.3f", nTT, nTT_err) << std::endl;
        std::cout << Form("ST                     = %.3f +- %.3f", nST, nST_err) << std::endl;

    }

    // Setup canvas and pads
    TCanvas * c1 = new TCanvas("c1", "c1", 700, 700);
    TPad * pad1 = new TPad("pad1", "top pad"   , 0.0, 0.3, 1.0, 1.0);
    pad1->SetBottomMargin(0.0);
    pad1->Draw();
    TPad * pad2 = new TPad("pad2", "bottom pad", 0.0, 0.0, 1.0, 0.3);
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.35);
    pad2->Draw();
    pad1->cd();
    pad1->SetLogy(plotLog);

    // Setup histogram styles
    set_style(hVH, "VH");
    set_style(hZH, "VH");
    set_style(hmonoH, "monoH");
    set_style(hWH, "VH");
    set_style(hWjLF, "WjLF");
    set_style(hWjHF, "WjHFb");
    set_style(hZjLF, "ZjLF");
    set_style(hZjHF, "ZjHFb");
    set_style(hTT, "TT");
    set_style(hST, "ST");
    set_style(hVV, "VVHF");
    set_style(hQCD, "QCD");
    set_style(hdata_obs, "data_obs");

    // Setup auxiliary histograms (ratios, errors, etc)
    TH1F * staterr = (TH1F *) hmc_exp->Clone("staterr");
    staterr->Sumw2();
    //staterr->SetFillColor(kRed);
    staterr->SetFillColor(kGray+3);
    staterr->SetMarkerSize(0);
    staterr->SetFillStyle(3013);

    TH1F * ratio = (TH1F *) hdata_obs->Clone("ratio");
    ratio->Sumw2();
    ratio->SetMarkerSize(0.8);
    //ratio->SetMarkerSize(0.5);
    ratio->Divide(hdata_obs, hmc_exp, 1., 1., "");

    TH1F * ratiostaterr = (TH1F *) hmc_exp->Clone("ratiostaterr");
    ratiostaterr->Sumw2();
    ratiostaterr->SetStats(0);
    ratiostaterr->SetTitle(title);
    ratiostaterr->GetYaxis()->SetTitle("Data/MC");
    ratiostaterr->SetMaximum(2.2);
    ratiostaterr->SetMinimum(0);
    ratiostaterr->SetMarkerSize(0);
    //ratiostaterr->SetFillColor(kRed);
    ratiostaterr->SetFillColor(kGray+3);
    ratiostaterr->SetFillStyle(3013);
    ratiostaterr->GetXaxis()->SetLabelSize(0.12);
    ratiostaterr->GetXaxis()->SetTitleSize(0.14);
    ratiostaterr->GetXaxis()->SetTitleOffset(1.10);
    ratiostaterr->GetYaxis()->SetLabelSize(0.10);
    ratiostaterr->GetYaxis()->SetTitleSize(0.12);
    ratiostaterr->GetYaxis()->SetTitleOffset(0.6);
    ratiostaterr->GetYaxis()->SetNdivisions(505);

    TLine* ratiounity = new TLine(xlow,1,xup,1);
    ratiounity->SetLineStyle(2);

    for (Int_t i = 0; i < hmc_exp->GetNbinsX()+2; i++) {
        ratiostaterr->SetBinContent(i, 1.0);
        if (hmc_exp->GetBinContent(i) > 1e-6) {  //< not empty
            double binerror = hmc_exp->GetBinError(i) / hmc_exp->GetBinContent(i);
            ratiostaterr->SetBinError(i, binerror);
        } else {
            ratiostaterr->SetBinError(i, 999.);
        }
    }

    TH1F * ratiosysterr = (TH1F *) ratiostaterr->Clone("ratiosysterr");
    ratiosysterr->Sumw2();
    ratiosysterr->SetMarkerSize(0);
    ratiosysterr->SetFillColor(kYellow-4);
    //ratiosysterr->SetFillStyle(3002);
    ratiosysterr->SetFillStyle(1001);

    for (Int_t i = 0; i < hmc_exp->GetNbinsX()+2; i++) {
        if (hmc_exp->GetBinContent(i) > 1e-6) {  //< not empty
            double binerror2 = (pow(hmc_exp->GetBinError(i), 2) +
                                pow(0.08 * hWjLF->GetBinContent(i), 2) +
                                pow(0.20 * hWjHF->GetBinContent(i), 2) +
                                pow(0.08 * hZjLF->GetBinContent(i), 2) +
                                pow(0.20 * hZjHF->GetBinContent(i), 2) +
                                pow(0.07 * hTT->GetBinContent(i), 2) +
                                pow(0.25 * hST->GetBinContent(i), 2) +
                                pow(0.25 * hVV->GetBinContent(i), 2)+
 	                        pow(0.50 * hQCD->GetBinContent(i), 2));

            double binerror = sqrt(binerror2);
            ratiosysterr->SetBinError(i, binerror / hmc_exp->GetBinContent(i));
        }
    }

    // Setup legends
    TLegend * leg1 = new TLegend(0.50, 0.68, 0.72, 0.92);
    set_style(leg1);
    if (plotData) leg1->AddEntry(hdata_obs, "Data", "p");
    if (plotSig)  leg1->AddEntry(hVH, Form("VH(%i)", ev->massH), "l");
    if (plotMonoH)  leg1->AddEntry(hmonoH, Form("monoH(%i)", ev->massH), "l");
    leg1->AddEntry(hTT, "t#bar{t}", "f");
    leg1->AddEntry(hST, "single top", "f");
    leg1->AddEntry(hVV, "VV", "f");

    TLegend * leg2 = new TLegend(0.72, 0.68, 0.94, 0.92);
    set_style(leg2);
    leg2->AddEntry(hWjHF, "W + HF", "f");
    leg2->AddEntry(hWjLF, "W + LF", "f");
    leg2->AddEntry(hZjHF, "Z + HF", "f");
    leg2->AddEntry(hZjLF, "Z + LF", "f");
    leg2->AddEntry(hQCD, "QCD", "f");
    leg2->AddEntry(staterr, "MC uncert. (stat)", "f");

    TLegend * ratioleg = new TLegend(0.72, 0.88, 0.94, 0.96);
    set_style(ratioleg);
    ratioleg->SetTextSize(0.07);
    ratioleg->AddEntry(ratiostaterr, "MC uncert. (stat)", "f");
    //   ratioleg->AddEntry(ratiosysterr, "MC uncert. (syst)", "f");

    TLegend * ratioleg2 = new TLegend(0.50, 0.88, 0.72, 0.96);
    set_style(ratioleg2);
    ratioleg2->SetTextSize(0.07);
    //    ratioleg->AddEntry(ratiostaterr, "MC uncert. (stat)", "f");
    ratioleg2->AddEntry(ratiosysterr, "MC uncert. (syst)", "f");


    // Draw MC signal and backgrounds
    std::clog << "MakePlots(): Drawing..." << std::endl;
    pad1->cd();
    if (plotLog) pad1->SetLogy();

    double ymax = TMath::Max(hdata_obs->GetMaximum(), hs->GetMaximum());
    hs->SetMaximum(ymax * 1.7 + (ymax>1 ? sqrt(ymax) : 0.));
    if (plotLog)  hs->SetMaximum(ymax * 200 + (ymax>1 ? sqrt(ymax) : 0.));
    hs->SetMinimum(0.01);
    hs->Draw("hist");
    hs->GetXaxis()->SetLabelSize(0);
    double binwidth = (xup - xlow) / nbinsx;
    TString ytitle = Form("Events / %.3f", binwidth);
    hs->GetYaxis()->SetTitle(ytitle);

    staterr->Draw("e2 same");
    if (plotSig) {
        hVH->SetLineColor(2);
        hVH->SetLineWidth(3);
        hVH->SetFillColor(0);
        hVH->Draw("hist same");


    }
    if (plotMonoH){
      hmonoH->SetLineColor(6);
      hmonoH->SetLineWidth(3);
      hmonoH->SetFillColor(0);
      hmonoH->Draw("hist same");
    } 
    if (plotData) {
        hdata_obs->Draw("e1 same");
    }

    // Draw legends
    leg1->Draw();
    leg2->Draw();
    TLatex * latex = new TLatex();
    latex->SetNDC();
    latex->SetTextAlign(12);
    latex->SetTextFont(62);
    latex->SetTextSize(0.04);
    //latex->DrawLatex(0.19, 0.89, "CMS Preliminary");
    latex->DrawLatex(0.19, 0.89, "CMS Simulation 2015");
    latex->SetTextSize(0.04);
    latex->DrawLatex(0.19, 0.84, "#sqrt{s} = 13 TeV, L = 20 fb^{-1}");
    // NOTE: change this to your channel
    //latex->DrawLatex(0.19, 0.79, "Z(#mu#bar{#mu})H(b#bar{b})");
    //latex->DrawLatex(0.19, 0.79, "Z(e#bar{e})H(b#bar{b})");
    //latex->DrawLatex(0.19, 0.79, "W(#mu#nu)H(b#bar{b})");
    //latex->DrawLatex(0.19, 0.79, "W(e#nu)H(b#bar{b})");
      latex->DrawLatex(0.19, 0.79, "Z(#nu#bar{#nu})H(b#bar{b})");
      if (plotMonoH)  latex->DrawLatex(0.19, 0.79, "monoH(b#bar{b})");

    // Draw ratio
    pad2->cd();
    pad2->SetGridy(0);
    ratiostaterr->Draw("e2");
    ratiosysterr->Draw("e2 same");
    ratiostaterr->Draw("e2 same");
    ratiounity->Draw();
    ratio->Draw("e1 same");
    ratioleg->Draw();
    ratioleg2->Draw();

    // Kolmogorov-Smirnov test and Chi2 test
    TPaveText * pave = new TPaveText(0.18, 0.86, 0.28, 0.96, "brNDC");
    if (plotData) {
        pave->SetLineColor(0);
        pave->SetFillColor(0);
        pave->SetShadowColor(0);
        pave->SetBorderSize(1);
        double nchisq = hdata_obs->Chi2Test(hmc_exp, "UWCHI2/NDF");  // MC uncert. (stat)
        //double kolprob = hdata_obs->KolmogorovTest(hmc_exp);  // MC uncert. (stat)
        TText * text = pave->AddText(Form("#chi_{#nu}^{2} = %.3f", nchisq));
        //TText * text = pave->AddText(Form("#chi_{#nu}^{2} = %.3f, K_{s} = %.3f", nchisq, kolprob));
        text->SetTextFont(62);
        text->SetTextSize(0.07);
        //text->SetTextSize(0.06);
        pave->Draw();
    }

    // Print
    std::clog << "MakePlots(): Printing..." << std::endl;
    pad1->cd();
    gPad->RedrawAxis();
    gPad->Modified();
    gPad->Update();
    pad2->cd();
    gPad->RedrawAxis();
    gPad->Modified();
    gPad->Update();
    c1->cd();

    gPad->Print(plotdir+plotname+".png");
    gPad->Print(plotdir+plotname+".pdf");

    TFile* outrootfile = TFile::Open(plotdir+plotname+".root", "RECREATE");
    hZH->Write();
    hWH->Write();
    hWjLF->Write();
    hWjHF->Write();
    hZjLF->Write();
    hZjHF->Write();
    hTT->Write();
    hST->Write();
    hVV->Write();
    hQCD->Write();
    hmc_exp->Write();
    hdata_obs->Write();
    outrootfile->Close();

    // Clean up
    delete staterr;
    delete ratio;
    delete ratiostaterr;
    delete ratiosysterr;
    delete leg1;
    delete leg2;
    delete ratioleg;
    delete ratioleg2;
    delete latex;
    delete pave;
    delete hs;
    delete pad1;
    delete pad2;
    delete c1;

    delete hZH;
    delete hWH;
    delete hWjLF;
    delete hWjHF;
    delete hZjLF;
    delete hZjHF;
    delete hTT;
    delete hST;
    delete hVV;
    delete hQCD;
    delete hdata_obs;
    delete hVH;
    delete hmonoH;
    delete hmc_exp;

    std::clog << "MakePlots(): DONE!" << std::endl;

    return;
}

//______________________________________________________________________________
void MakeDatacard(TString channel, TString dcname, TString wsname,
                  bool useshapes,
                  TString options="unblind:SplusB") {

    std::clog << "MakeDatacard(): Printing datacard..." << std::endl;

    // Parse options
    bool unblind    = options.Contains("unblind")    && (!options.Contains("!unblind"));
    bool SplusB     = options.Contains("SplusB")     && (!options.Contains("!SplusB"));
    channel += "_8TeV";

    if (useshapes && !unblind) {
        std::clog << "MakeDatacard(): I only support 'unblind' mode when using shapes." << std::endl;
        unblind = true;
    }

    std::ofstream dc;
    dc.setf(ios::fixed,ios::floatfield);
    dc.precision(3);
    dc.open(dcname.Data());

    TString wsname_new = dcname(0,dcname.Sizeof()-5)+".root";
    gSystem->Exec(Form("cp %s %s", wsname.Data(), wsname_new.Data()));
    TFile * wsfile = TFile::Open(wsname);

    int jmax = 8;
    dc << "imax 1 number of channels" << std::endl;
    dc << "jmax " << jmax << " number of processes minus 1 ('*' = automatic)" << std::endl;
    dc << "kmax * number of nuisance parameters (sources of systematical uncertainties)" << std::endl;
    dc << "-----------------------------------" << std::endl;
    if (useshapes) {
        //dc << "shapes * * " << wsname << " $CHANNEL:$PROCESS $CHANNEL:$PROCESS_$SYSTEMATIC" << std::endl;
        dc << "shapes * * " << wsname_new << " $PROCESS $PROCESS_$SYSTEMATIC" << std::endl;
        dc << "-----------------------------------" << std::endl;
    }
    dc << "bin         " << channel << std::endl;
    if (unblind) {
        dc << "observation " << ((TH1 *) wsfile->Get("data_obs"))->Integral() << std::endl;
    } else if (SplusB) {
        dc << "observation " << ((TH1 *) wsfile->Get("mc_exp"))->Integral() + ((TH1 *) wsfile->Get("ZH"))->Integral() + ((TH1 *) wsfile->Get("WH"))->Integral() << std::endl;
    } else {
        dc << "observation " << ((TH1 *) wsfile->Get("mc_exp"))->Integral() << std::endl;
    }
    dc << "#prediction " << ((TH1 *) wsfile->Get("mc_exp"))->Integral() << std::endl;
    dc << "-----------------------------------" << std::endl;

    dc << "bin         "; for (int j=0; j!=jmax+1; j++)  dc << channel << "   "; dc << std::endl;
    dc << "process     ZH         WH         WjLF       WjHF       ZjLF       ZjHF       TT         ST         VV         " << std::endl;
    dc << "process     -1         0          1          2          3          4          5          6          7          " << std::endl;
    dc << "rate        " << setw(10) << right << ((TH1 *) wsfile->Get("ZH"))->Integral() << " "
                         << setw(10) << right << ((TH1 *) wsfile->Get("WH"))->Integral() << " "
                         << setw(10) << right << ((TH1 *) wsfile->Get("WjLF"))->Integral() << " "
                         << setw(10) << right << ((TH1 *) wsfile->Get("WjHF"))->Integral() << " "
                         << setw(10) << right << ((TH1 *) wsfile->Get("ZjLF"))->Integral() << " "
                         << setw(10) << right << ((TH1 *) wsfile->Get("ZjHF"))->Integral() << " "
                         << setw(10) << right << ((TH1 *) wsfile->Get("TT"))->Integral() << " "
                         << setw(10) << right << ((TH1 *) wsfile->Get("ST"))->Integral() << " "
                         << setw(10) << right << ((TH1 *) wsfile->Get("VV"))->Integral() << std::endl;
    dc.precision(2);
    dc << "-----------------------------------" << std::endl;
    dc << "" << std::endl;
    dc << "# This datacard is for demonstration purposes only!" << std::endl;
    dc << "#################################### #####  ZH    WH    WjLF  WjHF  ZjLF  ZjHF  TT    ST    VV    " << std::endl;
    dc << "lumi_8TeV                            lnN    1.026 1.026 -     -     -     -     -     1.026 1.026 " << std::endl;
    dc << "pdf_qqbar                            lnN    1.01  1.01  -     -     -     -     -     -     1.01  " << std::endl;
    dc << "pdf_gg                               lnN    -     -     -     -     -     -     -     1.01  -     " << std::endl;
    dc << "QCDscale_VH                          lnN    1.04  1.04  -     -     -     -     -     -     -     " << std::endl;
    dc << "QCDscale_ttbar                       lnN    -     -     -     -     -     -     -     1.06  -     " << std::endl;
    dc << "QCDscale_VV                          lnN    -     -     -     -     -     -     -     -     1.04  " << std::endl;
    dc << "QCDscale_QCD                         lnN    -     -     -     -     -     -     -     -     -     " << std::endl;
    dc << "CMS_vhbb_boost_EWK_8TeV              lnN    1.05  1.10  -     -     -     -     -     -     -     " << std::endl;
    dc << "CMS_vhbb_boost_QCD_8TeV              lnN    1.10  1.10  -     -     -     -     -     -     -     " << std::endl;
    dc << "CMS_vhbb_ST                          lnN    -     -     -     -     -     -     -     1.25  -     " << std::endl;
    dc << "CMS_vhbb_VV                          lnN    -     -     -     -     -     -     -     -     1.25  " << std::endl;
    dc << "CMS_vhbb_eff_b                       lnN    1.07  1.07  1.07  1.07  1.07  1.07  1.07  1.07  1.07  " << std::endl;
    dc << "CMS_vhbb_fake_b_8TeV                 lnN    1.03  1.03  1.03  1.03  1.03  1.03  1.03  1.03  1.03  " << std::endl;
    dc << "CMS_vhbb_res_j                       lnN    1.05  1.05  1.05  1.05  1.05  1.05  1.05  1.05  1.05  " << std::endl;
    dc << "CMS_vhbb_scale_j                     lnN    1.05  1.05  1.05  1.05  1.05  1.05  1.05  1.05  1.05  " << std::endl;
    dc << "CMS_vhbb_WjLF_SF_" << channel << "            lnN    -     -     1.08  -     -     -     -     -     -     " << std::endl;
    dc << "CMS_vhbb_WjHF_SF_" << channel << "            lnN    -     -     -     1.20  -     -     -     -     -     " << std::endl;
    dc << "CMS_vhbb_ZjLF_SF_" << channel << "            lnN    -     -     -     -     1.08  -     -     -     -     " << std::endl;
    dc << "CMS_vhbb_ZjHF_SF_" << channel << "            lnN    -     -     -     -     -     1.20  -     -     -     " << std::endl;
    dc << "CMS_vhbb_TT_SF_" << channel << "              lnN    -     -     -     -     -     -     1.07  -     -     " << std::endl;
    if (channel == "Zmm_8TeV" || channel == "Wmn_8TeV")
        dc << "CMS_vhbb_eff_m                       lnN    1.01  1.01  -     -     -     -     -     1.01  1.01  " << std::endl;
    else if (channel == "Zee_8TeV" || channel == "Wen_8TeV")
        dc << "CMS_vhbb_eff_e                       lnN    1.01  1.01  -     -     -     -     -     1.01  1.01  " << std::endl;
    else if (channel == "Znn_8TeV")
        dc << "CMS_vhbb_trigger_MET                 lnN    1.03  1.03  -     -     -     -     -     1.03  1.03  " << std::endl;
    dc << "#################################### #####  ZH    WH    WjLF  WjHF  ZjLF  ZjHF  TT    ST    VV    " << std::endl;
    dc.close();

    wsfile->Close();

    std::clog << "MakeDatacard(): DONE!" << std::endl;

    return;
}




//______________________________________________________________________________
Events::Events()
  : ZH(0),
    WH(0),
    WjLF(0),
    WjHF(0),
    ZjLF(0),
    ZjHF(0),
    TT(0),
    ST_T_s(0),
    ST_T_t(0),
    ST_T_tW(0),
    ST_Tbar_s(0),
    ST_Tbar_t(0),
    ST_Tbar_tW(0),
    VV_WW(0),
    VV_WZ(0),
    VV_ZZ(0),
    data_obs(0),
    sf_WjLF(1.), sf_WjHF(1.), sf_ZjLF(1.), sf_ZjHF(1.), sf_TT(1.),
    massH(125) {}

Events::~Events() {
    delete ZH;
    delete WH;
    delete WjLF;
    delete WjHF;
    delete ZjLF;
    delete ZjHF;
    delete TT;
    delete ST_T_s;
    delete ST_T_t;
    delete ST_T_tW;
    delete ST_Tbar_s;
    delete ST_Tbar_t;
    delete ST_Tbar_tW;
    delete VV_WW;
    delete VV_WZ;
    delete VV_ZZ;
    delete data_obs;
}

void Events::read( TCut cutmc_all, TCut cutdata_all, TString processes , TString treename ) {
    //TString indir = "/uscms_data/d2/jiafu/CMSDAS2014/VHbbAnalysis/skim/";
    //TString indir = "/eos/uscms/store/user/cmsdas/2014/Hbb/Step2/";
  //    TString indir = "/afs/cern.ch/work/d/degrutto/public/MiniAOD/ZnnHbb_Phys14_PU20bx25/skimV11/step3/";
    //TString indir = "/afs/cern.ch/work/d/degrutto/public/MiniAOD/ZnnHbb_Phys14_PU20bx25/skimV11/step3/skim_ZnnH_classification/step4/"; //Step4_ZnnH125.root
    TString indir = "/afs/cern.ch/work/s/swang373/private/MiniAOD/ZnnHbb_Phys14_PU20bx25/skimV11/step3/skim_ZnnH_classification/step4/"; // swang373 01.07.2015

    TString prefix = "Step4_";
    TString suffix = ".root";
    //     TString treename = "tree_ZnunuHighPt_test";
     //TString treename = "tree";

    TCut cutHF = "(abs(Jet_mcFlavour[hJCidx[0]])==5 || abs(Jet_mcFlavour[hJCidx[1]])==5)";  //< for b quarks, pdgId = 5
    TCut cutLF = "(abs(Jet_mcFlavour[hJCidx[0]])!=5 && abs(Jet_mcFlavour[hJCidx[1]])!=5)";  //< for non-b quarks

    std::map<std::string, float> lumis = GetLumis();

    // parse processes
    bool loadAll  = processes == "" || processes.Contains("ALL") || processes.Contains("All") || processes.Contains("all");

    bool loadVH   = processes.Contains("VH");
    bool loadmonoH   = processes.Contains("monoH") ;
    bool loadZH   = processes.Contains("ZH") || loadVH || loadAll;
    bool loadWH   = processes.Contains("WH") || loadVH || loadAll;
    bool loadWJ   = processes.Contains("WJ") || processes.Contains("WjLF") || processes.Contains("WjHF") || loadAll;
    bool loadZJ   = processes.Contains("ZJ") ||   processes.Contains("ZjLF") ||  processes.Contains("ZjHF") ||   loadAll;
    bool loadTT   = processes.Contains("TT") || loadAll;
    bool loadST   = processes.Contains("ST") || loadAll;
    bool loadVV   = processes.Contains("VV") || loadAll;
    bool loadQCD   = processes.Contains("QCD") || loadAll;
    bool loadData = processes.Contains("DATA") || processes.Contains("Data") || processes.Contains("data") || loadAll;


    // Monte Carlo______________________________________________________________
    // NOTE: for Zll, use the default "ZllH"
    // NOTE: for Znn, change "ZllH" to "ZnnH"
    if (loadZH) {
        TChain ZH_(treename);
        //ZH_.Add(indir + prefix + Form("ZllH%i", massH) + suffix);
        //ZH = (TTree *) ZH_.CopyTree(cutmc_all);
        //lumi_ZH = lumis[Form("ZllH%i", massH)];
        ZH_.Add(indir + prefix + Form("ZnnH%i", massH) + suffix);
        ZH = (TTree *) ZH_.CopyTree(cutmc_all);
        lumi_ZH = lumis[Form("ZnnH%i", massH)];
        std::clog << "... DONE: ZH copy tree. N=" << ZH->GetEntries() << std::endl;
    }

    if (loadmonoH) {
        TChain monoH_(treename);
        //ZH_.Add(indir + prefix + Form("ZllH%i", massH) + suffix);
        //ZH = (TTree *) ZH_.CopyTree(cutmc_all);
        //lumi_ZH = lumis[Form("ZllH%i", massH)];
        monoH_.Add(indir + prefix + "monoH" + suffix);
        monoH = (TTree *) monoH_.CopyTree(cutmc_all);
        lumi_monoH = lumis["monoH"];
        std::clog << "... DONE: monoH copy tree. N=" << monoH->GetEntries() << std::endl;
    }



    if (loadWH) {
        TChain WH_(treename);
        WH_.Add(indir + prefix + Form("WlnH%i", massH) + suffix);
        WH = (TTree *) WH_.CopyTree(cutmc_all);
        lumi_WH = lumis[Form("WlnH%i", massH)];
        std::clog << "... DONE: WH copy tree. N=" << WH->GetEntries() << std::endl;
    }

    if (loadWJ) {
        TChain WjLF_(treename);
	//        WjLF_.Add(indir + prefix + "WJetsPtW100" + suffix);
        WjLF_.Add(indir + prefix + "WJets" + suffix);
        WjLF = (TTree *) WjLF_.CopyTree(cutmc_all + cutLF);
        lumi_WjLF = lumis["WJetsIncl"];
        std::clog << "... DONE: WjLF copy tree. N=" << WjLF->GetEntries() << std::endl;

        TChain WjHF_(treename);
        WjHF_.Add(indir + prefix + "WJets" + suffix);
        WjHF = (TTree *) WjHF_.CopyTree(cutmc_all + cutHF);
        lumi_WjHF = lumis["WJetsIncl"];
        std::clog << "... DONE: WjHF copy tree. N=" << WjHF->GetEntries() << std::endl;
    }

    if (loadZJ) {
 
        TChain ZjLF_(treename);
	//        ZjLF_.Add(indir + prefix + "WJetsPtW100" + suffix);
        ZjLF_.Add(indir + prefix + "ZJets" + suffix);
        ZjLF = (TTree *) ZjLF_.CopyTree(cutmc_all + cutLF);
        lumi_ZjLF = lumis["ZJetsHT600"];
        std::clog << "... DONE: ZjLF copy tree. N=" << ZjLF->GetEntries() << std::endl;

        TChain ZjHF_(treename);
        ZjHF_.Add(indir + prefix + "ZJets" + suffix);
        ZjHF = (TTree *) ZjHF_.CopyTree(cutmc_all + cutHF);
        lumi_ZjHF = lumis["ZJetsHT600"];
        std::clog << "... DONE: ZjHF copy tree. N=" << ZjHF->GetEntries() << std::endl;

      /*

        TChain ZjLF_HT100_(treename);
        ZjLF_HT100_.Add(indir + prefix + "ZJetsHT100" + suffix);
        ZjLF_HT100 = (TTree *) ZjLF_HT100_.CopyTree(cutmc_all + cutLF);
        lumi_ZjLFHT100 = lumis["ZJetsHT100"];
        std::clog << "... DONE: ZjLF HT100 copy tree. N=" << ZjLF_HT100->GetEntries() << std::endl;


        TChain ZjLF_HT200_(treename);
        ZjLF_HT200_.Add(indir + prefix + "ZJetsHT200" + suffix);
        ZjLF_HT200 = (TTree *) ZjLF_HT200_.CopyTree(cutmc_all + cutLF);
        lumi_ZjLFHT200 = lumis["ZJetsHT200"];
        std::clog << "... DONE: ZjLF HT200 copy tree. N=" << ZjLF_HT200->GetEntries() << std::endl;

	TChain ZjLF_HT600_(treename);
        ZjLF_HT600_.Add(indir + prefix + "ZJetsHT600" + suffix);
        ZjLF_HT600 = (TTree *) ZjLF_HT600_.CopyTree(cutmc_all + cutLF);
        lumi_ZjLFHT600 = lumis["ZJetsHT600"];
        std::clog << "... DONE: ZjLF HT600 copy tree. N=" << ZjLF_HT600->GetEntries() << std::endl;



        TChain ZjHF_HT100_(treename);
        ZjHF_HT100_.Add(indir + prefix + "ZJetsHT100" + suffix);
        ZjHF_HT100 = (TTree *) ZjHF_HT100_.CopyTree(cutmc_all + cutHF);
        lumi_ZjHFHT100 = lumis["ZJetsHT100"];
        std::clog << "... DONE: ZjHF HT100 copy tree. N=" << ZjHF_HT100->GetEntries() << std::endl;


        TChain ZjHF_HT200_(treename);
        ZjHF_HT200_.Add(indir + prefix + "ZJetsHT200" + suffix);
        ZjHF_HT200 = (TTree *) ZjHF_HT200_.CopyTree(cutmc_all + cutHF);
        lumi_ZjHFHT200 = lumis["ZJetsHT200"];
        std::clog << "... DONE: ZjHF HT200 copy tree. N=" << ZjHF_HT200->GetEntries() << std::endl;

	TChain ZjHF_HT600_(treename);
        ZjHF_HT600_.Add(indir + prefix + "ZJetsHT600" + suffix);
        ZjHF_HT600 = (TTree *) ZjHF_HT600_.CopyTree(cutmc_all + cutHF);
        lumi_ZjHFHT600 = lumis["ZJetsHT600"];
        std::clog << "... DONE: ZjHF HT600 copy tree. N=" << ZjHF_HT600->GetEntries() << std::endl;

      */



    }

    /* if (loadZjLF) { */
    /*     // NOTE: for Zll, use the default "DYJetsPtZ100" */
    /*     // NOTE: for Znn, change "DYJetsPtZ100" to "ZJetsPtZ100" */
    /*     TChain ZjLF_(treename); */
    /*     //ZjLF_.Add(indir + prefix + "DYJetsPtZ100" + suffix); */
    /*     //ZjLF = (TTree *) ZjLF_.CopyTree(cutmc_all + cutLF); */
    /*     //lumi_ZjLF = lumis["DYJetsPtZ100"]; */
    /*     ZjLF_.Add(indir + prefix + "ZJets" + suffix); */
    /*     ZjLF = (TTree *) ZjLF_.CopyTree(cutmc_all + cutLF); */
    /*     lumi_ZjLF = lumis["ZJetsHT600"]; */
    /*     std::clog << "... DONE: ZjLF copy tree. N=" << ZjLF->GetEntries() << std::endl; */

    /* } */


    /* if (loadZjHF) { */
    /*     // NOTE: for Zll, use the default "DYJetsPtZ100" */
    /*     // NOTE: for Znn, change "DYJetsPtZ100" to "ZJetsPtZ100" */
    /*     TChain ZjHF_(treename); */
    /*     //ZjLF_.Add(indir + prefix + "DYJetsPtZ100" + suffix); */
    /*     //ZjLF = (TTree *) ZjLF_.CopyTree(cutmc_all + cutLF); */
    /*     //lumi_ZjLF = lumis["DYJetsPtZ100"]; */
    /*     ZjHF_.Add(indir + prefix + "ZJets" + suffix); */
    /*     ZjHF = (TTree *) ZjHF_.CopyTree(cutmc_all + cutHF); */
    /*     lumi_ZjHF = lumis["ZJetsHT600"]; */
    /*     std::clog << "... DONE: ZjHF copy tree. N=" << ZjHF->GetEntries() << std::endl; */

    /* } */




    if (loadTT) {
        TChain TT_(treename);
        TT_.Add(indir + prefix + "TTPythia8" + suffix);
	//	    	TT_.Add(indir + prefix + "TTMadv2" + suffix);
        TT = (TTree *) TT_.CopyTree(cutmc_all);
       	lumi_TT = lumis["TTPythia8"];
	// lumi_TT = lumis["TTMadv2"];
        std::clog << "... DONE: TT copy tree. N=" << TT->GetEntries() << std::endl;

    }

    if (loadST) {
      /*
        TChain ST_T_s_(treename);
        ST_T_s_.Add(indir + prefix + "T_s" + suffix);
        ST_T_s = (TTree *) ST_T_s_.CopyTree(cutmc_all);
        lumi_ST_T_s = lumis["T_s"];
        std::clog << "... DONE: ST_T_s copy tree. N=" << ST_T_s->GetEntries() << std::endl;

        TChain ST_T_t_(treename);
        ST_T_t_.Add(indir + prefix + "T_t" + suffix);
        ST_T_t = (TTree *) ST_T_t_.CopyTree(cutmc_all);
        lumi_ST_T_t = lumis["T_t"];
        std::clog << "... DONE: ST_T_t copy tree. N="  << ST_T_t->GetEntries() << std::endl;

        TChain ST_T_tW_(treename);
        ST_T_tW_.Add(indir + prefix + "T_tW" + suffix);
        ST_T_tW = (TTree *) ST_T_tW_.CopyTree(cutmc_all);
        lumi_ST_T_tW = lumis["T_tW"];
        std::clog << "... DONE: ST_T_tW copy tree. N=" << ST_T_tW->GetEntries() << std::endl;

        TChain ST_Tbar_s_(treename);
        ST_Tbar_s_.Add(indir + prefix + "Tbar_s" + suffix);
        ST_Tbar_s = (TTree *) ST_Tbar_s_.CopyTree(cutmc_all);
        lumi_ST_Tbar_s = lumis["Tbar_s"];
        std::clog << "... DONE: ST_Tbar_s copy tree. N=" << ST_Tbar_s->GetEntries() << std::endl;

        TChain ST_Tbar_t_(treename);
        ST_Tbar_t_.Add(indir + prefix + "Tbar_t" + suffix);
        ST_Tbar_t = (TTree *) ST_Tbar_t_.CopyTree(cutmc_all);
        lumi_ST_Tbar_t = lumis["Tbar_t"];
        std::clog << "... DONE: ST_Tbar_t copy tree. N=" << ST_Tbar_t->GetEntries() << std::endl;

        TChain ST_Tbar_tW_(treename);
        ST_Tbar_tW_.Add(indir + prefix + "Tbar_tW" + suffix);
        ST_Tbar_tW = (TTree *) ST_Tbar_tW_.CopyTree(cutmc_all);
        lumi_ST_Tbar_tW = lumis["Tbar_tW"];
        std::clog << "... DONE: ST_Tbar_tW copy tree. N=" << ST_Tbar_tW->GetEntries() << std::endl; 
      */
                                                                                                                                                                                               
      TChain s_Top_(treename);                                                                                                                                                                
      s_Top_.Add(indir + prefix + "s_Top" + suffix);                                                                                                                                        
      s_Top = (TTree *) s_Top_.CopyTree(cutmc_all);                                                                                                                                      
      lumi_s_Top = lumis["Tbar_tW"];                                                                                                                                                          
      std::clog << "... DONE: s_Top copy tree. N=" << s_Top->GetEntries() << std::endl; 

    }

    if (loadVV) {
      /*
        TChain VV_WW_(treename);
        VV_WW_.Add(indir + prefix + "WW" + suffix);
        VV_WW = (TTree *) VV_WW_.CopyTree(cutmc_all);
        lumi_VV_WW = lumis["WW"];
        std::clog << "... DONE: VV_WW copy tree. N=" << VV_WW->GetEntries() << std::endl;

        TChain VV_WZ_(treename);
        VV_WZ_.Add(indir + prefix + "WZ" + suffix);
        VV_WZ = (TTree *) VV_WZ_.CopyTree(cutmc_all);
        lumi_VV_WZ = lumis["WZ"];
        std::clog << "... DONE: VV_WZ copy tree. N=" << VV_WZ->GetEntries() << std::endl;

        TChain VV_ZZ_(treename);
        VV_ZZ_.Add(indir + prefix + "ZZ" + suffix);
        VV_ZZ = (TTree *) VV_ZZ_.CopyTree(cutmc_all);
        lumi_VV_ZZ = lumis["ZZ"];
        std::clog << "... DONE: VV_ZZ copy tree. N=" << VV_ZZ->GetEntries() << std::endl;
      */

    }



    if (loadQCD) {
        TChain QCD_(treename);
	QCD_.Add(indir + prefix + "QCD" + suffix);
        QCD = (TTree *) QCD_.CopyTree(cutmc_all);
        lumi_QCD = lumis["QCD_HT1000"];
        std::clog << "... DONE: QCD copy tree. N=" << QCD->GetEntries() << std::endl;
    }





    // Data_____________________________________________________________________
    // NOTE: for Zmm and Wmn, use the default "SingleMu"
    // NOTE: for Wen, change both "SingleMu" to "SingleEl"
    // NOTE: for Zee, change both "SingleMu" to "DoubleEl"
    // NOTE: for Znn, change both "SingleMu" to "MET"
    if (loadData) {
        TChain data_obs_(treename);
        //data_obs_.Add(indir + prefix + "SingleMu_ReReco" + suffix);
        //data_obs_.Add(indir + prefix + "SingleMu_Prompt" + suffix);
        //data_obs_.Add(indir + prefix + "SingleEl_ReReco" + suffix);
        //data_obs_.Add(indir + prefix + "SingleEl_Prompt" + suffix);
        //data_obs_.Add(indir + prefix + "DoubleEl_ReReco" + suffix);
        //data_obs_.Add(indir + prefix + "DoubleEl_Prompt" + suffix);

	/*        data_obs_.Add(indir + prefix + "MET_ReReco" + suffix);
        data_obs_.Add(indir + prefix + "MET_Prompt" + suffix);
        data_obs = (TTree *) data_obs_.CopyTree(cutdata_all);
        std::clog << "... DONE: data_obs copy tree. N=" << data_obs->GetEntries() << std::endl;
	*/
    }

    return;
}


void set_style(TH1 * h, const TString& p) {
    if (p == "VH") {
        h->SetFillColor  (2);
        h->SetMarkerColor(2);
    } else if (p == "data_obs" || p == "Data") {
        h->SetMarkerSize(0.8);
        h->SetMarkerStyle(20);
    } else {
        h->SetLineColor(kBlack);
        if (p == "WjLF") {
            h->SetFillColor  (814);
            h->SetMarkerColor(814);
        } else if (p == "WjHFc") {
            h->SetFillColor  (816);
            h->SetMarkerColor(816);
        } else if (p == "WjHFb") {
            h->SetFillColor  (820);
            h->SetMarkerColor(820);
        } else if (p == "ZjLF") {
            h->SetFillColor  (401);
            h->SetMarkerColor(401);
        } else if (p == "ZjHFc") {
            h->SetFillColor  (41);
            h->SetMarkerColor(41);
        } else if (p == "ZjHFb") {
            h->SetFillColor  (5);
            h->SetMarkerColor(5);
        } else if (p == "TT") {
            h->SetFillColor  (596);
            h->SetMarkerColor(596);
        } else if (p == "ST") {
            h->SetFillColor  (840);
            h->SetMarkerColor(840);
        } else if (p == "VV") {
            h->SetFillColor  (922);
            h->SetMarkerColor(922);
        } else if (p == "VVHF") {
            h->SetFillColor  (920);
            h->SetMarkerColor(920);
        } else if (p == "QCD") {
            h->SetFillColor  (616);
            h->SetMarkerColor(616);
        } else if (p == "monoH") {
            h->SetFillColor  (6);
            h->SetMarkerColor(6);
        }
    }
    return;
}


void set_style_norm(TH1 * h, const TString& p) {
    if (p == "VH") {
        h->SetMaximum(10.0); 
        h->SetLineColor(2);
        h->SetLineWidth(3);
        h->SetMarkerColor(2);
    } else if (p == "data_obs" || p == "Data") {
        h->SetMaximum(10.0); 
        h->SetMarkerSize(0.8);
        h->SetLineWidth(3);
        h->SetMarkerStyle(20);
    } else {
        h->SetLineColor(kBlack);
        if (p == "WjLF") {
        h->SetMaximum(10.0); 
	  h->SetLineColor  (814);
        h->SetLineWidth(3);
            h->SetMarkerColor(814);
        } else if (p == "WjHFc") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
            h->SetLineColor  (816);
            h->SetMarkerColor(816);
        } else if (p == "WjHFb") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
            h->SetLineColor  (820);
            h->SetMarkerColor(820);
        } else if (p == "ZjLF") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
            h->SetLineColor  (401);
            h->SetMarkerColor(401);
        } else if (p == "ZjHFc") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
            h->SetLineColor  (41);
            h->SetMarkerColor(41);
        } else if (p == "ZjHFb") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
            h->SetLineColor  (5);
            h->SetMarkerColor(5);
        } else if (p == "TT") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
            h->SetLineColor  (596);
            h->SetMarkerColor(596);
        } else if (p == "ST") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
            h->SetLineColor  (840);
            h->SetMarkerColor(840);
        } else if (p == "VV") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
            h->SetLineColor  (922);
            h->SetMarkerColor(922);
        } else if (p == "VVHF") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
            h->SetLineColor  (920);
            h->SetMarkerColor(920);
        } else if (p == "QCD") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
            h->SetLineColor  (616);
            h->SetMarkerColor(616);
        } else if (p == "monoH") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
            h->SetLineColor  (6);
            h->SetMarkerColor(6);
        }

    }
    return;
}


void set_style(TLegend* leg) {
    leg->SetFillColor(0);
    leg->SetLineColor(0);
    leg->SetShadowColor(0);
    leg->SetTextFont(62);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(1);
    return;
}





