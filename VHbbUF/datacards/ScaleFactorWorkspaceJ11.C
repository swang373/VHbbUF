#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"

#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooDataHist.h"


////////////////////////////////////////////////////////////////////////////////
/// Configuration                                                            ///
////////////////////////////////////////////////////////////////////////////////

const double g_xlow         = 0.00;
const double g_xup          = 2.16;
const Int_t  g_upcol        = TColor::GetColor("#FF6600");
const Int_t  g_downcol      = TColor::GetColor("#0033FF");

const Int_t jmax            = 11;
const Int_t nsyst           = 13 + 2 * 3;  // plus WjModel, ZjModel, TTModel
const TString g_realvarname = "CMS_vhbb_VAR_$CHANNEL_8TeV";  // $CHANNEL is replaced in MakeWorkSpace()
const TString g_dcname      = "vhbb_Znn_SF_J11_$CHANNEL_8TeV.txt";
const TString g_wsname      = "vhbb_Znn_SF_J11_$CHANNEL_8TeV.root";
const TString g_rootname    = "vhbb_Znn_SF_J11_$CHANNEL_TH1.root";
const TString g_pdfname     = "vhbb_Znn_SF_J11_$CHANNEL_TH1.pdf";

/// Search bins
const char* g_channels[] = {
    "ZnunuHighPt",
    "ZnunuMedPt",
    "ZnunuLowPt",
    //"ZnunuLowCSV",
};

const TString g_processes[jmax+1] = {
    "ZH", "WH", "Wj0b", "Wj1b", "Wj2b", "Zj0b", "Zj1b", "Zj2b", "TT", "s_Top", "VV", "QCD"
};

const TString g_systematics[nsyst] = {
    "NONE", 
    "CMS_vhbb_res_jUp", 
    "CMS_vhbb_res_jDown", 
    "CMS_vhbb_scale_jUp", 
    "CMS_vhbb_scale_jDown", 
    "CMS_vhbb_eff_bUp", 
    "CMS_vhbb_eff_bDown", 
    "CMS_vhbb_fake_b_8TeVUp", 
    "CMS_vhbb_fake_b_8TeVDown", 
    "UEPSUp", 
    "UEPSDown", 
    "CMS_vhbb_trigger_MET_Znunu_8TeVUp", 
    "CMS_vhbb_trigger_MET_Znunu_8TeVDown",
    //"CMS_vhbb_stat$PROCESS_$CHANNEL_8TeVUp",
    //"CMS_vhbb_stat$PROCESS_$CHANNEL_8TeVDown",
    "CMS_vhbb_WJModel_Znunu_8TeVUp",
    "CMS_vhbb_WJModel_Znunu_8TeVDown",
    "CMS_vhbb_ZJModel_Znunu_8TeVUp",
    "CMS_vhbb_ZJModel_Znunu_8TeVDown",
    "CMS_vhbb_TTModel_Znunu_8TeVUp",
    "CMS_vhbb_TTModel_Znunu_8TeVDown"
};

/// Regions
const TString g_regions[1+5] = {
    "VH", 
    "ZjLF", 
    "ZjHF", 
    "WjLF",
    "WjHF", 
    "TT"
};
////////////////////////////////////////////////////////////////////////////////
/// END Configuration                                                        ///
////////////////////////////////////////////////////////////////////////////////

using namespace RooFit;

void MakeSystPlot(const TString& channel, TFile * input, RooWorkspace * ws, const RooArgList * obs, 
                  const Int_t p, const Int_t up, const Int_t down)
{
    TString pdfname = g_pdfname;
    pdfname.ReplaceAll("$CHANNEL", channel);

    const TString& process = g_processes[p];
    TString systUp = g_systematics[up];
    systUp.ReplaceAll("$CHANNEL", channel);
    systUp.ReplaceAll("$PROCESS", process);
    
    TString systDown = g_systematics[down];
    systDown.ReplaceAll("$CHANNEL", channel);
    systDown.ReplaceAll("$PROCESS", process);

    if (process != "Wj0b" && process != "Wj1b" && process != "Wj2b") {
        if (systUp.Contains("WJModel") || systDown.Contains("WJModel"))
            return;
    }

    if (process != "Zj0b" && process != "Zj1b" && process != "Zj2b") {
        if (systUp.Contains("ZJModel") || systDown.Contains("ZJModel"))
            return;
    }
    
    if (process != "TT") {
        if (systUp.Contains("TTModel") || systDown.Contains("TTModel"))
            return;
    }
    
    TH1F * h = (TH1F *) input->Get(channel + "/" + process);
    TH1F * hUp = (TH1F *) input->Get(channel + "/" + process + "_" + systUp);
    TH1F * hDown = (TH1F *) input->Get(channel + "/" + process + "_" + systDown);
    if ((h->Integral() > 0. && hUp->Integral() <= 0.) || h->Integral() <= 0.) {
        TString clonename = hUp->GetName();
        delete hUp;
        hUp = (TH1F*) h->Clone(clonename);
    }
    if ((h->Integral() > 0. && hDown->Integral() <= 0.) || h->Integral() <= 0.) {
        TString clonename = hDown->GetName();
        delete hDown;
        hDown = (TH1F*) h->Clone(clonename);
    }
    if (process == "WH" || process == "ZH" || process == "VH") {
        if (systUp.Contains("eff_b"))
            systUp.ReplaceAll("eff_b", "eff_b_SIG");
        if (systDown.Contains("eff_b"))
            systDown.ReplaceAll("eff_b", "eff_b_SIG");
    }
    RooDataHist * dhUp = new RooDataHist(process + "_" + systUp, "", *obs, hUp);
    ws->import(*dhUp);
    RooDataHist * dhDown = new RooDataHist(process + "_" + systDown, "", *obs, hDown);
    ws->import(*dhDown);
    
    h->SetStats(0);
    h->SetTitle("; BDT");
    h->SetLineColor(1);
    h->SetLineWidth(2);
    h->SetFillColor(0);
    h->SetMarkerStyle(20);
    h->SetMinimum(0.01);
    h->GetXaxis()->CenterTitle();
    
    hUp->SetLineColor(g_upcol);
    hUp->SetLineWidth(2);
    hUp->SetFillColor(0);
    
    hDown->SetLineColor(g_downcol);
    hDown->SetLineWidth(2);
    hDown->SetFillColor(0);
    
    h->Draw("e1");
    hUp->Draw("hist same");
    hDown->Draw("hist same");
    h->Draw("e1 same");
    
    TLegend * leg = new TLegend(0.35, 0.20, 0.92, 0.35);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetLineColor(0);
    leg->SetShadowColor(0);
    leg->SetTextFont(62);
    //leg->SetTextSize(0.015);
    leg->AddEntry(h, process, "pl");
    leg->AddEntry(hUp, systUp, "l");
    leg->AddEntry(hDown, systDown, "l");
    leg->Draw();

    gPad->RedrawAxis();
    gPad->Modified();
    gPad->Update();
    gPad->Print(pdfname);
    
    delete dhUp;
    delete dhDown;
    delete leg;
    
    return;
}

void MakeWorkspace(const TString& channel, const TString& strategy)
{
    TString dcname       = g_dcname;
    TString wsname       = g_wsname;
    TString realvarname  = g_realvarname;
    TString rootname     = g_rootname;
    TString pdfname      = g_pdfname;
    dcname      .ReplaceAll("$CHANNEL", channel);
    wsname      .ReplaceAll("$CHANNEL", channel);
    realvarname .ReplaceAll("$CHANNEL", channel);
    rootname    .ReplaceAll("$CHANNEL", channel);
    pdfname     .ReplaceAll("$CHANNEL", channel);
    
    RooWorkspace * ws = new RooWorkspace(channel + "_8TeV", channel + "_8TeV");
    ws->factory(Form("%s[%.2f,%.2f]", realvarname.Data(), g_xlow, g_xup));
    RooRealVar * realvar = ws->var(realvarname.Data());
    RooArgList * obs = new RooArgList(*realvar);
    
    TFile * input = TFile::Open(rootname, "READ");
    TH1F * h(0);

    if (strategy == "blind") {
        h = (TH1F *) input->Get(channel + "/" + "mc_exp");
    } else if (strategy == "injectsignal") {
        h = (TH1F *) input->Get(channel + "/" + "mc_exp");
        TH1F * hZH = (TH1F *) input->Get(channel + "/" + "ZH_SM");
        h->Add(hZH);
        TH1F * hWH = (TH1F *) input->Get(channel + "/" + "WH_SM");
        h->Add(hWH);
        delete hZH;
        delete hWH;
    } else if (strategy == "real") {
        h = (TH1F *) input->Get(channel + "/" + "data_obs");
    } else {
        std::cerr << "Unknown strategy: " << strategy << std::endl;
    }

    if (strategy != "real")  gSystem->Exec(Form("sed -i 's/observation.*/observation %.3f/' %s ", h->Integral(), dcname.Data()));

    // data_obs
    RooDataHist * dh = new RooDataHist("data_obs", "", *obs, h);
    ws->import(*dh);
    delete dh;
    delete h;

    // no systematics
    for (Int_t p=0; p < jmax+1; p++) {
        const TString& process = g_processes[p];
        h = (TH1F *) input->Get(channel + "/" + process);
        dh = new RooDataHist(process, "", *obs, h);
        ws->import(*dh);
        delete dh;
        delete h;
    }

    // with systematics up/down
    TCanvas * c1 = new TCanvas("c1", "c1", 700, 700);
    c1->SetLogy();
    c1->Print(pdfname+"[");
    for (Int_t p=0; p < jmax+1; p++) {
        for (Int_t s=1; s < nsyst; s+=2) {
            MakeSystPlot(channel, input, ws, obs, p, s, s+1);
        }
    }
    c1->Print(pdfname+"]");
    delete c1;
    
    ws->writeToFile(wsname, kTRUE);
    delete input;
    delete ws;
    
    return;
}

////////////////////////////////////////////////////////////////////////////////
/// Main                                                                     ///
////////////////////////////////////////////////////////////////////////////////
void ScaleFactorWorkspaceJ11(TString strategy="real")  /// strategy: either one of "real", "blind", "injectsignal"
{
    gROOT->LoadMacro("tdrstyle.C");
    gROOT->ProcessLine("setTDRStyle()");
    
    TH1::SetDefaultSumw2(1);
    gROOT->SetBatch(1);

    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

    //if (!TString(gROOT->GetVersion()).Contains("5.34")) {
    //    std::cout << "INCORRECT ROOT VERSION! Please use 5.34:" << std::endl;
    //    std::cout << "source /uscmst1/prod/sw/cms/slc5_amd64_gcc472/lcg/root/5.34.03-cms4/bin/thisroot.csh" << std::endl;
    //    std::cout << "Return without doing anything." << std::endl;
    //    return;
    //}

    /// Set the channels
    const std::vector<TString> channels(g_channels, g_channels + sizeof(g_channels)/sizeof(g_channels[0]));

    for (UInt_t ichan = 0; ichan < channels.size(); ichan++) {
        const TString& channel = channels.at(ichan);
        std::clog << "\n*** " << channel << " ***\n" << std::endl;
        
        int begin = 1;
        int end   = 1+5;
        if (channel == "ZnunuLowCSV") {
            MakeWorkspace(channel + "_ZjHF", strategy);
            MakeWorkspace(channel + "_WjHF", strategy);
        } else {
            for (Int_t ireg=begin; ireg<end; ireg++) {
                MakeWorkspace(channel + "_" + g_regions[ireg], strategy);
            }
        }
    }
    
    return;
}
