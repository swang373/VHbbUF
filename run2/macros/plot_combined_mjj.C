#include <iostream>

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "THStack.h"
#include "TCut.h"
#include "TString.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPaveText.h"


// TODO: add background-subtracted plot
// TODO: weighted by S/(S+B)?

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

void MakePlots(const std::vector<TString>& channels,
               TString plotname="plot", TString plotdir="") {

    std::clog << "MakePlots(): Using " << channels.size() << " channels." << std::endl;

    TString dir      = "datacards/";
    TString prefix   = "vhbb_shapes_";
    TString suffix   = "_8TeV.root";
    TString title    = "m(jj) [GeV]";

    bool plotSig      = true;
    bool plotData     = true;
    bool plotLog      = false;

    TH1F * hZH       = 0;
    TH1F * hWH       = 0;
    TH1F * hWjLF     = 0;
    TH1F * hWjHF     = 0;
    TH1F * hZjLF     = 0;
    TH1F * hZjHF     = 0;
    TH1F * hTT       = 0;
    TH1F * hST       = 0;
    TH1F * hVV       = 0;
    TH1F * hmc_exp   = 0;
    TH1F * hdata_obs = 0;


    for (unsigned int i=0; i < channels.size(); ++i) {
        TFile * infile = TFile::Open(dir+prefix+channels[i]+suffix);
        float weight = 1.0;

        if (i==0) {
            hZH       = (TH1F *) infile->Get("ZH");
            hWH       = (TH1F *) infile->Get("WH");
            hWjLF     = (TH1F *) infile->Get("WjLF");
            hWjHF     = (TH1F *) infile->Get("WjHF");
            hZjLF     = (TH1F *) infile->Get("ZjLF");
            hZjHF     = (TH1F *) infile->Get("ZjHF");
            hTT       = (TH1F *) infile->Get("TT");
            hST       = (TH1F *) infile->Get("ST");
            hVV       = (TH1F *) infile->Get("VV");
            hmc_exp   = (TH1F *) infile->Get("mc_exp");
            hdata_obs = (TH1F *) infile->Get("data_obs");
        } else {
            hZH       ->Add( (TH1F *) infile->Get("ZH"), weight);
            hWH       ->Add( (TH1F *) infile->Get("WH"), weight);
            hWjLF     ->Add( (TH1F *) infile->Get("WjLF"), weight);
            hWjHF     ->Add( (TH1F *) infile->Get("WjHF"), weight);
            hZjLF     ->Add( (TH1F *) infile->Get("ZjLF"), weight);
            hZjHF     ->Add( (TH1F *) infile->Get("ZjHF"), weight);
            hTT       ->Add( (TH1F *) infile->Get("TT"), weight);
            hST       ->Add( (TH1F *) infile->Get("ST"), weight);
            hVV       ->Add( (TH1F *) infile->Get("VV"), weight);
            hmc_exp   ->Add( (TH1F *) infile->Get("mc_exp"), weight);
            hdata_obs ->Add( (TH1F *) infile->Get("data_obs"), weight);
        }
    }

    TH1F * hVH = (TH1F *) hZH->Clone("VH");
    hVH->Add(hWH);

    int    nbinsx = hZH->GetNbinsX();
    double xlow   = hVH->GetBinLowEdge(1);
    double xup    = hVH->GetBinLowEdge(nbinsx+1);

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
    if (plotSig)  hs->Add(hVH);

    //if (printStat) {
    //    double nS = hVH->Integral();
    //    //double nB = ((TH1*) hs->GetStack()->At(hs->GetHists()->GetSize()-1))->Integral();  //< integral of all added histograms in the stack
    //    double nB = hmc_exp->Integral();
    //    double signif_simple = nS / sqrt(nS + nB);  // only for reference
    //    double signif_punzi = nS / (3.0/2.0 + sqrt(nB) + 0.2*nB);
    //    std::cout << "MakePlots(): Printing FOM..." << std::endl;
    //    std::cout << Form("S                      = %.3f", nS) << std::endl;
    //    std::cout << Form("B                      = %.3f", nB) << std::endl;
    //    std::cout << Form("S/sqrt(S+B)            = %.5f", signif_simple) << std::endl;
    //    std::cout << Form("S/(sqrt(B)+a/2+ 0.2*B) = %.5f", signif_punzi) << std::endl;
    //}

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
    set_style(hWH, "VH");
    set_style(hWjLF, "WjLF");
    set_style(hWjHF, "WjHFb");
    set_style(hZjLF, "ZjLF");
    set_style(hZjHF, "ZjHFb");
    set_style(hTT, "TT");
    set_style(hST, "ST");
    set_style(hVV, "VVHF");
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
                                pow(0.25 * hVV->GetBinContent(i), 2));

            double binerror = sqrt(binerror2);
            ratiosysterr->SetBinError(i, binerror / hmc_exp->GetBinContent(i));
        }
    }

    // Setup legends
    TLegend * leg1 = new TLegend(0.50, 0.68, 0.72, 0.92);
    set_style(leg1);
    if (plotData) leg1->AddEntry(hdata_obs, "Data", "p");
    //if (plotSig)  leg1->AddEntry(hVH, Form("VH(%i)", ev->massH), "l");
    if (plotSig)  leg1->AddEntry(hVH, "VH(125)", "l");
    leg1->AddEntry(hTT, "t#bar{t}", "f");
    leg1->AddEntry(hST, "single top", "f");
    leg1->AddEntry(hVV, "VV", "f");

    TLegend * leg2 = new TLegend(0.72, 0.68, 0.94, 0.92);
    set_style(leg2);
    leg2->AddEntry(hWjHF, "W + HF", "f");
    leg2->AddEntry(hWjLF, "W + LF", "f");
    leg2->AddEntry(hZjHF, "Z + HF", "f");
    leg2->AddEntry(hZjLF, "Z + LF", "f");
    leg2->AddEntry(staterr, "MC uncert. (stat)", "f");

    TLegend * ratioleg = new TLegend(0.72, 0.88, 0.94, 0.96);
    set_style(ratioleg);
    ratioleg->SetTextSize(0.07);
    ratioleg->AddEntry(ratiostaterr, "MC uncert. (stat)", "f");

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
    TString ytitle = Form("Events / %.0f GeV", binwidth);
    hs->GetYaxis()->SetTitle(ytitle);

    staterr->Draw("e2 same");
    if (plotSig) {
        hVH->SetLineColor(2);
        hVH->SetLineWidth(3);
        hVH->SetFillColor(0);
        hVH->Draw("hist same");
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
    latex->SetTextSize(0.052);
    latex->DrawLatex(0.19, 0.89, "CMS Preliminary");
    latex->SetTextSize(0.04);
    latex->DrawLatex(0.19, 0.84, "#sqrt{s} = 8 TeV, L = 18.9 fb^{-1}");
    // NOTE: change this to your channel
    latex->DrawLatex(0.19, 0.79, "Z(#mu#bar{#mu})H(b#bar{b})");
    //latex->DrawLatex(0.19, 0.79, "Z(e#bar{e})H(b#bar{b})");
    //latex->DrawLatex(0.19, 0.79, "W(#mu#nu)H(b#bar{b})");
    //latex->DrawLatex(0.19, 0.79, "W(e#nu)H(b#bar{b})");
    //latex->DrawLatex(0.19, 0.79, "Z(#nu#bar{#nu})H(b#bar{b})");

    // Draw ratio
    pad2->cd();
    pad2->SetGridy(0);
    ratiostaterr->Draw("e2");
    //ratiosysterr->Draw("e2 same");
    ratiostaterr->Draw("e2 same");
    ratiounity->Draw();
    ratio->Draw("e1 same");
    ratioleg->Draw();

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
    delete hdata_obs;
    delete hVH;
    delete hmc_exp;

    std::clog << "MakePlots(): DONE!" << std::endl;

    return;
}

void plot_combined_mjj() {

    const TString channels_[] = {
        "Zmm",
        "Zee",
        "Wmn",
        "Wen",
        "Znn"
    };

    const std::vector<TString> channels(channels_, channels_ + sizeof(channels_)/sizeof(channels_[0]));

    gROOT->LoadMacro("tdrstyle.C");
    gROOT->ProcessLine("setTDRStyle()");

    TH1::SetDefaultSumw2(1);
    TH1::AddDirectory(0);


    TString plotdir  = "plots/";
    TString plotname = "vhbb_combined_mjj";

    MakePlots(channels, plotname, plotdir);

}