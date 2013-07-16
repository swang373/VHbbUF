#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TMath.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"

#include "RooArgList.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooHistPdf.h"
#include "RooRandom.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"


//#define MJJANALYSIS

#define STATBINBYBIN

//#define GOODNUISHAPE
#ifdef GOODNUISHAPE
#include "HelperBDTShape.h"
#endif


////////////////////////////////////////////////////////////////////////////////
/// Configuration                                                            ///
////////////////////////////////////////////////////////////////////////////////


const int    g_upcol        = TColor::GetColor("#FF6600");
const int    g_downcol      = TColor::GetColor("#0033FF");
#ifndef MJJANALYSIS
const double g_xlow         = -1.;
const double g_xup          =  1.;
const TString g_realvarname = "CMS_vhbb_BDT_$CHANNEL_8TeV";  // $CHANNEL is replaced in MakeWorkSpace()
#else
const double g_xlow         =   0.;
const double g_xup          = 250.;
const TString g_realvarname = "CMS_vhbb_MJJ_$CHANNEL_8TeV";  // $CHANNEL is replaced in MakeWorkSpace()
#endif
const TString g_wsname      = "vhbb_Znn_J12_8TeV.root";
//const TString g_wsname      = "vhbb_Znn_J12_$CHANNEL_8TeV.root";
const TString g_dcname      = "vhbb_Znn_J12_$CHANNEL_8TeV.txt";
const TString g_dcbbbname   = "vhbb_Znn_J12_bbb_$CHANNEL_8TeV.txt";
const TString g_rootname    = "vhbb_Znn_J12_$CHANNEL_TH1.root";
const TString g_pdfname     = "vhbb_Znn_J12_$CHANNEL_TH1.pdf";

const int jmax              = 12;  // as in datacard
const int nsyst             = 1 + 2 * 7 + 2 * 6;  // plus stat, WJModel, ZJModel, TTModel, WJSlope, ZJSlope
//const int massH             = 125;

/// Search bins
const char* g_channels[] = {
    "ZnunuHighPt",
    "ZnunuMedPt",
    "ZnunuLowPt",
    //"ZnunuLowCSV",
};

const TString g_processes[jmax+1] = {
    "ZH", "WH", "Wj0b", "Wj1b", "Wj2b", "Zj0b", "Zj1b", "Zj2b", "TT", "s_Top", "VVLF", "VVHF", "QCD"
};

const TString g_systematics[nsyst] = {
    "NONE", 
    "CMS_vhbb_Znn_res_jUp", 
    "CMS_vhbb_Znn_res_jDown", 
    "CMS_vhbb_Znn_scale_jUp", 
    "CMS_vhbb_Znn_scale_jDown", 
    "CMS_vhbb_eff_bUp", 
    "CMS_vhbb_eff_bDown", 
    "CMS_vhbb_fake_b_8TeVUp", 
    "CMS_vhbb_fake_b_8TeVDown", 
    "UEPSUp", 
    "UEPSDown", 
    "CMS_vhbb_trigger_MET_Znn_8TeVUp", 
    "CMS_vhbb_trigger_MET_Znn_8TeVDown",
    "CMS_vhbb_trigger_CSV_Znn_8TeVUp", 
    "CMS_vhbb_trigger_CSV_Znn_8TeVDown",
    "CMS_vhbb_stat$PROCESS_$CHANNEL_8TeVUp",
    "CMS_vhbb_stat$PROCESS_$CHANNEL_8TeVDown",
    "CMS_vhbb_WJModel_Znn_8TeVUp",
    "CMS_vhbb_WJModel_Znn_8TeVDown",
    "CMS_vhbb_ZJModel_Znn_8TeVUp",
    "CMS_vhbb_ZJModel_Znn_8TeVDown",
    "CMS_vhbb_TTModel_Znn_8TeVUp",
    "CMS_vhbb_TTModel_Znn_8TeVDown",
    "CMS_vhbb_WJSlope_Znn_8TeVUp",
    "CMS_vhbb_WJSlope_Znn_8TeVDown",
    "CMS_vhbb_ZJSlope_Znn_8TeVUp",
    "CMS_vhbb_ZJSlope_Znn_8TeVDown",
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
    
    //std::cout << "VERBOSE: " << systUp << " " << systDown << std::endl;

    if (process != "Wj0b" && process != "Wj1b" && process != "Wj2b") {
        if (systUp.Contains("WJModel") || systDown.Contains("WJModel") || systUp.Contains("WJSlope") || systDown.Contains("WJSlope"))
            return;
    }
    if (process != "Zj0b" && process != "Zj1b" && process != "Zj2b") {
        if (systUp.Contains("ZJModel") || systDown.Contains("ZJModel") || systUp.Contains("ZJSlope") || systDown.Contains("ZJSlope"))
            return;
    }
    if (process != "TT") {
        if (systUp.Contains("TTModel") || systDown.Contains("TTModel"))
            return;
    }


    /// Get histograms
    TH1F * h = (TH1F *) input->Get(channel + "/" + process);
    TH1F * hUp = (TH1F *) input->Get(channel + "/" + process + "_" + systUp);
    TH1F * hDown = (TH1F *) input->Get(channel + "/" + process + "_" + systDown);

    /// hUp & hDown must be non-zero when h is non-zero
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


#if defined(GOODNUISHAPE) && !defined(MJJANALYSIS)
    if (process != "QCD" && process != "ZH" && process != "WH" && process != "VH") {  // FIXME
        if (systUp == "CMS_vhbb_eff_bUp" || 
            systUp == "CMS_vhbb_fake_b_8TeVUp" || 
            systUp == "CMS_vhbb_Znn_res_jUp" || 
            systUp == "CMS_vhbb_Znn_scale_jUp" || 
            systUp == "UEPSUp" || 
            systUp == "CMS_vhbb_trigger_MET_Znn_8TeVUp" || 
            systUp == "CMS_vhbb_trigger_CSV_Znn_8TeVUp" ||
            systUp == "CMS_vhbb_WJModel_Znn_8TeVUp" || 
            systUp == "CMS_vhbb_ZJModel_Znn_8TeVUp" ||
            systUp == "CMS_vhbb_TTModel_Znn_8TeVUp" ||
            systUp == "CMS_vhbb_WJSlope_Znn_8TeVUp" || 
            systUp == "CMS_vhbb_ZJSlope_Znn_8TeVUp" ) {
            
            UInt_t nbins = h->GetNbinsX();
            double errorf = 0.35; // should be the same as used in STATBINBYBIN
            const UInt_t begin = 1;
            const UInt_t end   = nbins+1;
            
            /// Group bins from the left
            UInt_t firstN = 1;  // starts with 1 bin
            findFirstN(firstN, begin, h, errorf);
            
            /// Group bins from the right
            UInt_t lastN = 1;  // starts with 1 bin
            findLastN(lastN, end, h, errorf);
            
            if (systUp == "CMS_vhbb_WJModel_Znn_8TeVUp" || 
                systUp == "CMS_vhbb_ZJModel_Znn_8TeVUp") {  // FIXME
                firstN = 1;  // starts with 1 bin
                findFirstN(firstN, begin, hUp, errorf);
                
                lastN = 1;  // starts with 1 bin
                findLastN(lastN, end, hUp, errorf);
            }
            
            //std::cout << process << " " << systUp << " " << firstN << " " << lastN << std::endl;
            assert(firstN > 0 && firstN < nbins && lastN > 0 && lastN < nbins);
            assert(begin+firstN-1 < end-1-lastN+1);
            
            /// Impose a maximum of non-empty bins
            UInt_t firstEmptyN = 0;
            for (UInt_t i=begin; i<end; i++) {
                if (h->GetBinContent(i) > 1E-7)  break;
                firstEmptyN++;
            }
            UInt_t lastEmptyN = 0;
            for (UInt_t i=end-1; i>begin-1; i--) {
                if (h->GetBinContent(i) > 1E-7)  break;
                lastEmptyN++;
            }
            
            if (systUp == "CMS_vhbb_WJModel_Znn_8TeVUp" || 
                systUp == "CMS_vhbb_ZJModel_Znn_8TeVUp") {  // FIXME
                firstN = firstN;
                lastN = lastN;
            } else {
                firstN = TMath::Min(firstN, (UInt_t) firstEmptyN + 3);
                lastN = TMath::Min(lastN, (UInt_t) lastEmptyN + 3);
            }
            
            double integralBegin     = h    ->Integral(begin, begin+firstN-1);
            double integralEnd       = h    ->Integral(end-1-lastN+1, end-1);
            double integralBeginUp   = hUp  ->Integral(begin, begin+firstN-1);
            double integralEndUp     = hUp  ->Integral(end-1-lastN+1, end-1);
            double integralBeginDown = hDown->Integral(begin, begin+firstN-1);
            double integralEndDown   = hDown->Integral(end-1-lastN+1, end-1);
            
            //std::cout << integralBegin << " " << integralBeginUp << " " << integralBeginDown << " " 
            //          << integralEnd << " " << integralEndUp << " " << integralEndDown << std::endl;
            //assert(integralBegin > 0 && integralBeginUp > 0 && integralBeginDown > 0 &&
            //       integralEnd > 0 && integralEndUp > 0 && integralEndDown > 0);  // FIXME
            
            for (UInt_t i=begin; i<begin+firstN-1+1; i++){
                const Double_t bincontent        = h->GetBinContent(i);
                const Double_t newbincontentUp   = bincontent * integralBeginUp / integralBegin;
                const Double_t newbincontentDown = bincontent * integralBeginDown / integralBegin;
                hUp->SetBinContent(i, newbincontentUp);
                hDown->SetBinContent(i, newbincontentDown);
            }
            
            for (UInt_t i=end-1-lastN+1; i<end-1+1; i++){
                const Double_t bincontent        = h->GetBinContent(i);
                const Double_t newbincontentUp   = bincontent * integralEndUp / integralEnd;
                const Double_t newbincontentDown = bincontent * integralEndDown / integralEnd;
                hUp->SetBinContent(i, newbincontentUp);
                hDown->SetBinContent(i, newbincontentDown);
            }
        }
    }
#endif  // defined(GOODNUISHAPE) && !defined(MJJANALYSIS)


    /// Fix normalizations for systematics that have control regions
    if (process == "Wj0b" || process == "Wj1b" || process == "Wj2b" || process == "Zj0b" || process == "Zj1b" || process == "Zj2b" || process == "TT" ) {
        if (systUp == "CMS_vhbb_eff_bUp" || 
            systUp == "CMS_vhbb_fake_b_8TeVUp" || 
            systUp == "CMS_vhbb_Znn_res_jUp" || 
            systUp == "CMS_vhbb_Znn_scale_jUp" || 
            systUp == "UEPSUp" || 
            systUp == "CMS_vhbb_trigger_MET_Znn_8TeVUp" || 
            systUp == "CMS_vhbb_trigger_CSV_Znn_8TeVUp" ) {
            hUp->Sumw2();
            hUp->Scale(h->Integral() / hUp->Integral());
            
            hDown->Sumw2();
            hDown->Scale(h->Integral() / hDown->Integral());
        }
    }
    
    /// Empty systUp and systDown bins when a norminal bin has zero content
    bool treatEmptyBins = true;
    if (treatEmptyBins) {
        UInt_t nbins = h->GetNbinsX();
        for (UInt_t ibin=1; ibin<nbins+1; ibin++) {
            if (h->GetBinContent(ibin) < 1E-7 && hUp->GetBinContent(ibin) > 1E-7) {
                std::cout << "WARNING: " << systUp << " " << process << " bin " << ibin << " has norminal: " << h->GetBinContent(ibin) << ", systUp: " << hUp->GetBinContent(ibin) << std::endl;
                std::cout << "WARNING: systUp is set to norminal!!" << std::endl;
                hUp->SetBinContent(ibin, h->GetBinContent(ibin));
            }
            if (h->GetBinContent(ibin) < 1E-7 && hDown->GetBinContent(ibin) > 1E-7) {
                std::cout << "WARNING: " << systDown << " " << process << " bin " << ibin << " has norminal: " << h->GetBinContent(ibin) << ", systDown: " << hDown->GetBinContent(ibin) << std::endl;
                std::cout << "WARNING: systDown is set to norminal!!" << std::endl;
                hDown->SetBinContent(ibin, h->GetBinContent(ibin));
            }
        }
    }


    /// Name changes
    // Decorrelate eff_b and eff_b_SIG
    //if (process == "WH" || process == "ZH" || process == "VH") {
    //    if (systUp.Contains("eff_b"))
    //        systUp.ReplaceAll("eff_b", "eff_b_SIG");
    //    if (systDown.Contains("eff_b"))
    //        systDown.ReplaceAll("eff_b", "eff_b_SIG");
    //}
    
    // Decorrelate trigger_CSV and trigger_CSV_fake
    if (process == "Wj0b" || process == "Zj0b" || process == "VVLF" || process == "QCD") {
        if (systUp.Contains("trigger_CSV"))
            systUp.ReplaceAll("trigger_CSV", "trigger_CSV_fake");
        if (systDown.Contains("trigger_CSV"))
            systDown.ReplaceAll("trigger_CSV", "trigger_CSV_fake");
    }
    
    // Decorrelate WJModel and ZJModel for 0b/1b/2b
    if (process == "Wj0b" || process == "Zj0b") {
        if (systUp.Contains("WJModel"))
            systUp.ReplaceAll("WJModel", "Wj0bModel");
        if (systDown.Contains("WJModel"))
            systDown.ReplaceAll("WJModel", "Wj0bModel");
        if (systUp.Contains("ZJModel"))
            systUp.ReplaceAll("ZJModel", "Zj0bModel");
        if (systDown.Contains("ZJModel"))
            systDown.ReplaceAll("ZJModel", "Zj0bModel");
    } else if (process == "Wj1b" || process == "Zj1b") {
        if (systUp.Contains("WJModel"))
            systUp.ReplaceAll("WJModel", "Wj1bModel");
        if (systDown.Contains("WJModel"))
            systDown.ReplaceAll("WJModel", "Wj1bModel");
        if (systUp.Contains("ZJModel"))
            systUp.ReplaceAll("ZJModel", "Zj1bModel");
        if (systDown.Contains("ZJModel"))
            systDown.ReplaceAll("ZJModel", "Zj1bModel");
    } else if (process == "Wj2b" || process == "Zj2b") {
        if (systUp.Contains("WJModel"))
            systUp.ReplaceAll("WJModel", "Wj2bModel");
        if (systDown.Contains("WJModel"))
            systDown.ReplaceAll("WJModel", "Wj2bModel");
        if (systUp.Contains("ZJModel"))
            systUp.ReplaceAll("ZJModel", "Zj2bModel");
        if (systDown.Contains("ZJModel"))
            systDown.ReplaceAll("ZJModel", "Zj2bModel");
    }


    /// Import
    RooDataHist * dhUp = new RooDataHist(process + "_" + systUp, "", *obs, hUp);
    ws->import(*dhUp);
    RooDataHist * dhDown = new RooDataHist(process + "_" + systDown, "", *obs, hDown);
    ws->import(*dhDown);


#ifdef STATBINBYBIN
    if (g_systematics[up] == "CMS_vhbb_stat$PROCESS_$CHANNEL_8TeVUp") {
        TString sedstring = Form("s/\\(%s.*\\)/", systUp.Data());
        sedstring.ReplaceAll("_8TeVUp", "\\)\\(_8TeV");
        
        TH1F * h_statUp(0);
        TH1F * h_statDown(0);
        UInt_t nbins = h->GetNbinsX();
        double errorf = 0.35;
        UInt_t nedgebins = TMath::Min(0.25 * nbins, 7.);  // decide who are the edge bins
        for (UInt_t ibin=1; ibin<nbins+1; ibin++) {

#ifndef MJJANALYSIS
            /// Skip if not an edge bin
            if (!(ibin < nedgebins || ibin > nbins-nedgebins))
                continue;
#endif
            
            h_statUp   = (TH1F *) h->Clone(Form("%s_CMS_vhbb_stat%s_%s_bin%i_8TeVUp"  , h->GetName(), h->GetName(), channel.Data(), ibin));
            h_statDown = (TH1F *) h->Clone(Form("%s_CMS_vhbb_stat%s_%s_bin%i_8TeVDown", h->GetName(), h->GetName(), channel.Data(), ibin));
            const Double_t bincontent        = h->GetBinContent(ibin);
            const Double_t binerror          = h->GetBinError(ibin);
            const Double_t newbincontentUp   = TMath::Max(0., bincontent + (binerror *  1.0));
            const Double_t newbincontentDown = TMath::Max(0., bincontent + (binerror * -1.0));
            h_statUp->SetBinContent(ibin, newbincontentUp);
            h_statDown->SetBinContent(ibin, newbincontentDown);
            
            if (process == "QCD") {  // inflate QCD errors as systematics are removed
                const Double_t newbincontentUp2  = TMath::Max(0., bincontent + (binerror *  1.25));
                const Double_t newbincontentDown2= TMath::Max(0., bincontent + (binerror * -1.25));
                h_statUp->SetBinContent(ibin, newbincontentUp2);
                h_statDown->SetBinContent(ibin, newbincontentDown2);
            }
            
            if ( !(TMath::AreEqualRel(bincontent, newbincontentUp, 1E-7) && TMath::AreEqualRel(bincontent, newbincontentDown, 1E-7)) &&
                  bincontent > 0.065 ) {  // only if different and bincontent > 0.065
                // Make process dependent threshold
                bool pass = false;
                if (process == "QCD") {
                    if (binerror/bincontent > errorf)  // error is already inflated
                        pass = true;
                } else if (process == "Wj0b" || process == "Zj0b" || process == "VVLF" || 
                           process == "Wj1b" || process == "Zj1b" ||
                           process == "s_Top") {  // include less from processes that don't contribute to the last bin
                    if (binerror/bincontent > errorf * 2)
                        pass = true;
                } else {
                    if (binerror/bincontent > errorf && ibin > nbins+1-1-nedgebins)  // only the right edge bins
                        pass = true;
                }
#ifdef MJJANALYSIS
                /// Always true for 105-150 GeV
                if (ibin >= h->FindFixBin(105+1) && ibin <= h->FindFixBin(150-1))
                    pass = true;
#endif
                
                if (pass) {
                    std::cout << "stat" << process << " " << ibin << " " << bincontent << " " << binerror/bincontent << std::endl;
                    sedstring += Form("\\1_bin%i\\2\\n", ibin);
                    
                    RooDataHist * dh_statUp = new RooDataHist(h_statUp->GetName(), "", *obs, h_statUp);
                    ws->import(*dh_statUp);
                    RooDataHist * dh_statDown = new RooDataHist(h_statDown->GetName(), "", *obs, h_statDown);
                    ws->import(*dh_statDown);
                    
                    delete dh_statUp;
                    delete dh_statDown;
                }
            }
            delete h_statUp;
            delete h_statDown;
        }
        sedstring += "/";
        
        TString dcname    = g_dcname;
        TString dcbbbname = g_dcbbbname;
        dcname   .ReplaceAll("$CHANNEL", channel);
        dcbbbname.ReplaceAll("$CHANNEL", channel);
        if (process == g_processes[0]) {
            gSystem->Exec(Form("cp %s %s", dcname.Data(), dcbbbname.Data()));
        }
        gSystem->Exec(Form("sed -i '%s' %s", sedstring.Data(), dcbbbname.Data()));
    }
#endif  // STATBINBYBIN


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
    TString dcbbbname    = g_dcbbbname;
    TString wsname       = g_wsname;
    TString realvarname  = g_realvarname;
    TString rootname     = g_rootname;
    TString pdfname      = g_pdfname;
    dcname      .ReplaceAll("$CHANNEL", channel);
    dcbbbname   .ReplaceAll("$CHANNEL", channel);
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
        
    } else if (strategy == "toy" || strategy == "toys") {
        h = (TH1F *) input->Get(channel + "/" + "mc_exp");
        int nbins = h->GetNbinsX();
        RooDataHist * dummy = new RooDataHist("dummy", "", *obs, h);
        RooHistPdf * toy = new RooHistPdf("toy", "", *realvar, *dummy, 0);
        RooDataSet * toyds = toy->generate(*realvar, int(h->Integral()));
        RooDataHist * dh = new RooDataHist("data_obs", "", *realvar, *toyds);
        delete h;
        h = (TH1F *) dh->createHistogram(realvarname.Data(), *realvar, Binning(nbins));
        delete dummy;
        delete toy;
        delete toyds;
        delete dh;
        
    } else if (strategy == "obs") {
        h = (TH1F *) input->Get(channel + "/" + "data_obs");
        
    } else {
        std::cerr << "Unknown strategy: " << strategy << std::endl;
    }

    // Import data_obs
    RooDataHist * dh = new RooDataHist("data_obs", "", *obs, h);
    ws->import(*dh);

    // Put the data_obs integral into the datacard
    gSystem->Exec(Form("sed -i 's/observation.*/observation %.3f/' %s ", h->Integral(), dcname.Data()));

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
    
    if (channel == TString(g_channels[0]))
        ws->writeToFile(wsname, kTRUE);
    else
        ws->writeToFile(wsname, kFALSE);
    delete input;
    delete ws;
    
    return;
}


////////////////////////////////////////////////////////////////////////////////
/// Main                                                                     ///
////////////////////////////////////////////////////////////////////////////////
void BDTShapeWorkspaceJ12(TString strategy="obs", int randomseed=0)  /// strategy: either one of "obs", "blind", "injectsignal", "toy"
{
    gROOT->LoadMacro("tdrstyle.C");
    gROOT->ProcessLine("setTDRStyle()");
    
    TH1::SetDefaultSumw2(1);
    gROOT->SetBatch(1);

    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    
    RooRandom::randomGenerator()->SetSeed(13+randomseed);

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
        MakeWorkspace(channel, strategy);
    }
        
    // Make copies when doing the signal injection
    if (strategy == "injectsignal") {
        TString wsname       = g_wsname;
        TString wsnameSI     = wsname;
        wsnameSI    .ReplaceAll("vhbb_Znn_", "vhbb_Znn_SI_");
        gSystem->Exec("cp " + wsname + " " + wsnameSI);
        
        for (UInt_t ichan = 0; ichan < channels.size(); ichan++) {
            const TString& channel = channels.at(ichan);
            TString dcname       = g_dcname;
            dcname      .ReplaceAll("$CHANNEL", channel);
            TString dcnameSI     = dcname;
            dcnameSI    .ReplaceAll("vhbb_Znn_", "vhbb_Znn_SI_");
            gSystem->Exec("cp " + dcname + " " + dcnameSI);
            gSystem->Exec(Form("sed -i 's/%s/%s/' %s ", wsname.Data(), wsnameSI.Data(), dcnameSI.Data()));
            
#ifdef STATBINBYBIN
            TString dcbbbname    = g_dcbbbname;
            dcbbbname   .ReplaceAll("$CHANNEL", channel);
            TString dcbbbnameSI  = dcbbbname;
            dcbbbnameSI .ReplaceAll("vhbb_Znn_", "vhbb_Znn_SI_");
            gSystem->Exec("cp " + dcbbbname + " " + dcbbbnameSI);
            gSystem->Exec(Form("sed -i 's/%s/%s/' %s ", wsname.Data(), wsnameSI.Data(), dcbbbnameSI.Data()));
#endif  // STATBINBYBIN
        }
    }
    
    return;
}
