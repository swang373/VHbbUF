////////////////////////////////////////////////////////////////////////////////
/// This macro runs on Step 4 with jmax = 10 in the following order
///   ZH         WH         WjLF       WjHF       ZjLF       ZjHF       TT         s_Top      VVLF       VVHF       QCD
///   9          10         1          2          3          4          5          6          7          0          8
///
/// It produces a datacard (set by dcname) and a root file with TH1F (set by rootname).
///
/// The following 1 + 2 x 9 systematic variations are considered:
///   NONE
///   CMS_vhbb_res_jUp
///   CMS_vhbb_res_jDown
///   CMS_vhbb_scale_jUp
///   CMS_vhbb_scale_jDown
///   CMS_vhbb_eff_bUp
///   CMS_vhbb_eff_bDown
///   CMS_vhbb_fake_b_8TeVUp
///   CMS_vhbb_fake_b_8TeVDown
///   UEPSUp
///   UEPSDown
///   CMS_vhbb_trigger_MET_Znunu_8TeVUp
///   CMS_vhbb_trigger_MET_Znunu_8TeVDown
///   CMS_vhbb_stat$PROCESS_$CHANNELUp
///   CMS_vhbb_stat$PROCESS_$CHANNELDown
///   CMS_vhbb_ZJModel_Znunu_8TeVUp
///   CMS_vhbb_ZJModel_Znunu_8TeVDown
///   CMS_vhbb_TTModel_Znunu_8TeVUp
///   CMS_vhbb_TTModel_Znunu_8TeVDown
////////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
//#include <map>

#include "TCanvas.h"
#include "TCut.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"


////////////////////////////////////////////////////////////////////////////////
/// Configuration                                                            ///
////////////////////////////////////////////////////////////////////////////////

const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130214/reload_20130216/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130214/stitch/";
//const TString indir     = "root://eoscms//eos/cms/store/caf/user/degrutto/2013/Step4_20130207/";
const TString prefix    = "Step4_";
const TString suffix    = ".root";
const int jmax          = 10;  // as in datacard
const int nsyst         = 13;
const int massH         = 125;

const TString g_var         = Form("BDTzzhfassignal_%i[ISYST]", massH);  // ISYST is replaced in the main function.
const TString g_varslice    = "";
//const TString g_varslice    = Form("slice(BDTttbar_%i[ISYST], BDTvjl_%i[ISYST], BDTzz_%i[ISYST], BDTregular_%i[ISYST], -0.5, -0.4, -0.4)", massH, massH, massH, massH);
const TString g_cutmc       = Form("2.0 * weightsHCP[ISYST] * (selectFlags[0][ISYST])");
const TString g_cutdata     = Form("selectFlags[0][ISYST]");
//const TString g_cutmc       = Form("2.0 * weightsHCP[ISYST] * (selectFlags[0][ISYST] && BDTttbar_125[ISYST]>-0.5 && BDTvjl_125[ISYST]>-0.4 && BDTzz_125[ISYST]>-0.4)");
//const TString g_cutdata     = Form("selectFlags[0][ISYST] && BDTttbar_125[ISYST]>-0.5 && BDTvjl_125[ISYST]>-0.4 && BDTzz_125[ISYST]>-0.4");
const TString g_dcname      = "vhbb_Znn_J10_$CHANNEL_8TeV.txt";
const TString g_wsname      = "vhbb_Znn_J10_$CHANNEL_8TeV.root";
const TString g_rootname    = "vhbb_Znn_J10_$CHANNEL_TH1.root";
const TCut    g_cutallmc    = "";  // used in CopyTree()
const TCut    g_cutalldata  = "";  // used in CopyTree()
const double  g_xlow        = -1.;
const double  g_xup         = 1.;

/// The channel
TString channel = "ZnunuHighPt";

/// Systematics
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
    "CMS_vhbb_trigger_MET_Znunu_8TeVDown"
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

/// Scale factors of WjLF, WjHF, ZjLF, ZjHF, TT (as in datacard)
/// for 1 + 2 x 6 systematics
const double g_scalefactors_ZnunuHighPt[5*nsyst] = {
    1.3, 1.0, 1.3, 1.0, 1.0, // NONE
    1.3, 1.0, 1.3, 1.0, 1.0, // CMS_vhbb_res_jUp
    1.3, 1.0, 1.3, 1.0, 1.0, // CMS_vhbb_res_jDown
    1.3, 1.0, 1.3, 1.0, 1.0, // CMS_vhbb_scale_jUp
    1.3, 1.0, 1.3, 1.0, 1.0, // CMS_vhbb_scale_jDown
    1.3, 1.0, 1.3, 1.0, 1.0, // CMS_vhbb_eff_bUp
    1.3, 1.0, 1.3, 1.0, 1.0, // CMS_vhbb_eff_bDown
    1.3, 1.0, 1.3, 1.0, 1.0, // CMS_vhbb_fake_b_8TeVUp
    1.3, 1.0, 1.3, 1.0, 1.0, // CMS_vhbb_fake_b_8TeVDown
    1.3, 1.0, 1.3, 1.0, 1.0, // UEPSUp
    1.3, 1.0, 1.3, 1.0, 1.0, // UEPSDown
    1.3, 1.0, 1.3, 1.0, 1.0, // CMS_vhbb_trigger_MET_Znunu_8TeVUp
    1.3, 1.0, 1.3, 1.0, 1.0, // CMS_vhbb_trigger_MET_Znunu_8TeVDown
};
const double g_scalefactors_lnN_ZnunuHighPt[5] = {
    1.04, 1.10, 1.06, 1.08, 1.06
};

const double g_scalefactors_ZnunuLowPt[5*nsyst] = {
    0.93, 1.00, 1.76, 1.00, 1.17, // NONE
    0.93, 1.00, 1.76, 1.00, 1.17, // CMS_vhbb_res_jUp
    0.93, 1.00, 1.76, 1.00, 1.17, // CMS_vhbb_res_jDown
    0.93, 1.00, 1.76, 1.00, 1.17, // CMS_vhbb_scale_jUp
    0.93, 1.00, 1.76, 1.00, 1.17, // CMS_vhbb_scale_jDown
    0.93, 1.00, 1.76, 1.00, 1.17, // CMS_vhbb_eff_bUp
    0.93, 1.00, 1.76, 1.00, 1.17, // CMS_vhbb_eff_bDown
    0.93, 1.00, 1.76, 1.00, 1.17, // CMS_vhbb_fake_b_8TeVUp
    0.93, 1.00, 1.76, 1.00, 1.17, // CMS_vhbb_fake_b_8TeVDown
    0.93, 1.00, 1.76, 1.00, 1.17, // UEPSUp
    0.93, 1.00, 1.76, 1.00, 1.17, // UEPSDown
    0.93, 1.00, 1.76, 1.00, 1.17, // CMS_vhbb_trigger_MET_Znunu_8TeVUp
    0.93, 1.00, 1.76, 1.00, 1.17, // CMS_vhbb_trigger_MET_Znunu_8TeVDown
};
const double g_scalefactors_lnN_ZnunuLowPt[5] = {
    1.04, 1.10, 1.06, 1.08, 1.06
};
 
const double g_scalefactors_ZnunuLowCSV[5*nsyst] = {
    0.93, 1.00, 1.33, 1.00, 1.17,  // NONE
    0.93, 1.00, 1.33, 1.00, 1.17,  // CMS_vhbb_res_jUp
    0.93, 1.00, 1.33, 1.00, 1.17,  // CMS_vhbb_res_jDown
    0.93, 1.00, 1.33, 1.00, 1.17,  // CMS_vhbb_scale_jUp
    0.93, 1.00, 1.33, 1.00, 1.17,  // CMS_vhbb_scale_jDown
    0.93, 1.00, 1.33, 1.00, 1.17,  // CMS_vhbb_eff_bUp
    0.93, 1.00, 1.33, 1.00, 1.17,  // CMS_vhbb_eff_bDown
    0.93, 1.00, 1.33, 1.00, 1.17,  // CMS_vhbb_fake_b_8TeVUp
    0.93, 1.00, 1.33, 1.00, 1.17,  // CMS_vhbb_fake_b_8TeVDown
    0.93, 1.00, 1.33, 1.00, 1.17,  // UEPSUp
    0.93, 1.00, 1.33, 1.00, 1.17,  // UEPSDown
    0.93, 1.00, 1.33, 1.00, 1.17,  // CMS_vhbb_trigger_MET_Znunu_8TeVUp
    0.93, 1.00, 1.33, 1.00, 1.17,  // CMS_vhbb_trigger_MET_Znunu_8TeVDown
};
const double g_scalefactors_lnN_ZnunuLowCSV[5] = {
    1.07, 1.10, 1.08, 1.10, 1.07
};
////////////////////////////////////////////////////////////////////////////////
/// END Configuration                                                        ///
////////////////////////////////////////////////////////////////////////////////

///_____________________________________________________________________________
/// Classes

class EventsJ10 {  // jmax = 10 in datacard, with 5 scale factors.
public:
    EventsJ10()
      : ZH(0),
        WH(0),
        WjLF(0),
        WjHF(0),
        ZjLF(0),
        ZjHF(0),
        TT(0),
        s_Top(0),
        VVLF(0),
        VVHF(0),
        QCD(0),
        ZjLFHW(0),
        ZjHFHW(0),
        TTMG(0),
        ZH_SM(0),
        WH_SM(0),
        //mc_exp(0),
        data_obs(0),
        sf_WjLF(1.), sf_WjHF(1.), sf_ZjLF(1.), sf_ZjHF(1.), sf_TT(1.) {}

    ~EventsJ10()
    {
        delete ZH;
        delete WH;
        delete WjLF;
        delete WjHF;
        delete ZjLF;
        delete ZjHF;
        delete TT;
        delete s_Top;
        delete VVLF;
        delete VVHF;
        delete QCD;
        delete ZjLFHW;
        delete ZjHFHW;
        delete TTMG;
        delete ZH_SM;
        delete WH_SM;
        //delete mc_exp;
        delete data_obs;
    }

    void check() const {
        assert(ZH       != 0 && ZH      ->GetEntriesFast() > 0);
        assert(WH       != 0 && WH      ->GetEntriesFast() > 0);
        assert(WjLF     != 0 && WjLF    ->GetEntriesFast() > 0);
        assert(WjHF     != 0 && WjHF    ->GetEntriesFast() > 0);
        assert(ZjLF     != 0 && ZjLF    ->GetEntriesFast() > 0);
        assert(ZjHF     != 0 && ZjHF    ->GetEntriesFast() > 0);
        assert(TT       != 0 && TT      ->GetEntriesFast() > 0);
        assert(s_Top    != 0 && s_Top   ->GetEntriesFast() > 0);
        assert(VVLF     != 0 && VVLF    ->GetEntriesFast() > 0);
        assert(VVHF     != 0 && VVHF    ->GetEntriesFast() > 0);
        //assert(QCD      != 0 && QCD     ->GetEntriesFast() > 0);
        assert(ZjLFHW   != 0 && ZjLFHW  ->GetEntriesFast() > 0);
        assert(ZjHFHW   != 0 && ZjHFHW  ->GetEntriesFast() > 0);
        assert(TTMG     != 0 && TTMG    ->GetEntriesFast() > 0);
        //assert(mc_exp   != 0 && mc_exp  ->GetEntriesFast() > 0);
        assert(ZH_SM    != 0 && ZH_SM   ->GetEntriesFast() > 0);
        assert(WH_SM    != 0 && WH_SM   ->GetEntriesFast() > 0);
        assert(data_obs != 0 && data_obs->GetEntriesFast() > 0);
        return;
    }
    
    void set_scale_factors(const Double_t sf[5]) {
        sf_WjLF = sf[0];
        sf_WjHF = sf[1];
        sf_ZjLF = sf[2];
        sf_ZjHF = sf[3];
        sf_TT   = sf[4];
        return;
    }

public:
    TTree * ZH;
    TTree * WH;
    TTree * WjLF;
    TTree * WjHF;
    TTree * ZjLF;
    TTree * ZjHF;
    TTree * TT;
    TTree * s_Top;
    TTree * VVLF;
    TTree * VVHF;
    TTree * QCD;  
    TTree * ZjLFHW; // for ZJModel
    TTree * ZjHFHW; // for ZJModel
    TTree * TTMG;   // for TTModel
    TTree * ZH_SM;  // for signal injection
    TTree * WH_SM;  // for signal injection
    //TTree * mc_exp;
    TTree * data_obs;
    Double_t sf_WjLF;
    Double_t sf_WjHF;
    Double_t sf_ZjLF;
    Double_t sf_ZjHF;
    Double_t sf_TT;
};  // EventsJ10


///_____________________________________________________________________________
/// Functions

double slice(const double bdt1, const double bdt2, const double bdt3, const double bdt4, 
             const double cut1, const double cut2, const double cut3,
             const double xlow=-1., const double xup=1.)
{
    // Assume original range is [-1,1], check if this is always true!
    if(!(-1. <= bdt1 && bdt1 <= 1.))  std::cerr << "ERROR: bdt1 = " << bdt1 << std::endl;
    if(!(-1. <= bdt2 && bdt2 <= 1.))  std::cerr << "ERROR: bdt2 = " << bdt2 << std::endl;
    if(!(-1. <= bdt3 && bdt3 <= 1.))  std::cerr << "ERROR: bdt3 = " << bdt3 << std::endl;
    if(!(-1. <= bdt4 && bdt4 <= 1.))  std::cerr << "ERROR: bdt4 = " << bdt4 << std::endl;
    double scale = (xup - xlow) / 2.0 / 4.0;
    if (bdt1 < cut1) {
        return ((bdt4+1.0)*scale + xlow);
    } else if (bdt2 < cut2) {
        return ((bdt4+3.0)*scale + xlow);
    } else if (bdt3 < cut3) {
        return ((bdt4+5.0)*scale + xlow);
    } else { // pass all 3 cuts
        return ((bdt4+7.0)*scale + xlow);
    }
}

void setHisto(TH1 * h, const TString& p)
{
    //h->Sumw2();
    if (p == "VH") {
        //h->SetLineColor(kBlack);
        h->SetLineColor(kRed);
        h->SetFillColor(kRed);
    } else if (p == "data_obs") {
        h->SetMarkerSize(0.8);
        h->SetMarkerStyle(20);
    } else {
        h->SetLineColor(kBlack);
        if (p == "WjLF") {
            h->SetFillColor(kSpring - 6);
        } else if (p == "WjHFc") {
            h->SetFillColor(kSpring - 3);
        } else if (p == "WjHFb") {
            h->SetFillColor(kSpring);
        } else if (p == "ZjLF") {
            h->SetFillColor(kYellow - 6);
        } else if (p == "ZjHFc") {
            h->SetFillColor(kYellow - 3);
        } else if (p == "ZjHFb") {
            h->SetFillColor(kYellow);
        } else if (p == "TT") {
            h->SetFillColor(kBlue);
        } else if (p == "s_Top") {
            h->SetFillColor(kTeal);
        } else if (p == "VV") {
            h->SetFillColor(kGray);
        } else if (p == "VVHF") {
            h->SetFillColor(kGray + 2);
        } else if (p == "QCD") {
            h->SetFillColor(kMagenta);
        }
    }
}

void FormatFileName(TString& str)
{
    str.ReplaceAll(" ", "");
    const char *specials = "[]{}()$+*?&";
    if (str.First(specials) != -1) {
        str.ReplaceAll("[","_");
        str.ReplaceAll("]","_");
        str.ReplaceAll("{","_");
        str.ReplaceAll("}","_");
        str.ReplaceAll("(","_");
        str.ReplaceAll(")","_");
        str.ReplaceAll(" ","_");
        str.ReplaceAll("$","_");
        str.ReplaceAll("+","_");
        str.ReplaceAll("*","_");
        str.ReplaceAll("?","_");
        str.ReplaceAll("&","_");
    }
    str.Remove(TString::kTrailing, '_');
}

/// newnbins        : the number of bins in the returned TH1F.
/// firstN, lastN   : they are passed by references. If they are set to 
///                   negative values, the function will find the proper 
///                   firstN and lastN and return them, so they can be 
///                   used to rebin other histograms.
///                   In short, they should be set once, then used by 
///                   others. They are supposed to be set by a histogram 
///                   of sum of MC backgrounds (no signal).
/// threshold       : threshold on the error fraction in the first and last bins.
/// newxlow, newxup : the boundaries of the returned TH1F. If they are
///                   the same, the boundaries of the input TH1F are used.
TH1F * RebinByStatErrors(const TH1F * h, Int_t newnbins, Int_t& firstN, Int_t& lastN, Double_t threshold=0.35, Double_t newxlow=0., Double_t newxup=0., const char * newname="")
{
    TH1F * hdummy = (TH1F *) h->Clone("dummy");
    hdummy->Sumw2();
    Int_t nbins = h->GetNbinsX();
    if (h->GetXaxis()->IsVariableBinSize()){
        hdummy->Error("Rebin", "Input histogram has variable bin sizes! h=%s", h->GetName());
        return hdummy;
    }
    if (threshold < 0.0 || threshold > 1.0){
        hdummy->Warning("Rebin", "threshold is not within [0,1]. Call TH1F::Rebin(). threshold=%f", threshold);
        hdummy->Rebin(float(hdummy->GetNbinsX())/float(newnbins));
        hdummy->SetName(newname);
        return hdummy;
    }
    if (nbins == newnbins){
        hdummy->Warning("Rebin", "newnbins is equal to current nbins. Return a clone. nbins=%i, newnbins=%i", nbins, newnbins);
        hdummy->SetName(newname);
        return hdummy;
    }
    if (newnbins < 3 || newnbins > nbins){
        hdummy->Error("Rebin", "newnbins must be within [3,%i]. newnbins=%i", nbins, newnbins);
        return hdummy;
    }
    if (firstN > nbins || lastN > nbins){
        hdummy->Error("Rebin", "firstN and lastN must be < %i. firstN=%i, lastN=%i", nbins, firstN, lastN);
        return hdummy;
    }
    if (h->GetBinContent(0) > 0. || h->GetBinContent(nbins+1) > 0.){
        hdummy->Error("Rebin", "Underflow and/or overflow is not empty! h=%s", h->GetName());
        return hdummy;
    }
    TH1F * hnew = (TH1F *) h->Clone(Form("%s_test1", newname));  // temp histogram with the same binning
    if (h->GetSumw2N() == 0){
        hdummy->Warning("Rebin", "Sumw2 is not created prior to this. h=%s", h->GetName());
        hnew->Sumw2();
    }
    UInt_t begin=1, end=nbins+1;
    // If firstN and lastN >= 1, don't change them.
    Bool_t isfixed = (firstN >= 1) && (lastN >= 1);
    // Group bins from the left
    if (firstN<1){
        Double_t firstbincontent  = hnew->GetBinContent(begin);
        Double_t firstbinerror    = hnew->GetBinError(begin);
        Double_t firstbinerror2   = (firstbinerror * firstbinerror);
        Double_t firstbinerrorf   = (firstbincontent > 1e-6) ? (TMath::Sqrt(firstbinerror2) / firstbincontent) : 1.0;
        firstN = 1;  // starts with 1 bin
        // Find firstN such that its error fraction becomes less than the threshold
        while ((firstbinerrorf > threshold) && (firstN <= nbins)){
            //hdummy->Info("Rebin", "firstN=%2i int. binc=%f int. bine=%f frac=%f", firstN, firstbincontent, sqrt(firstbinerror2), firstbinerrorf);
            
            Int_t i = begin + firstN;
            const Double_t bincontent     = hnew->GetBinContent(i);
            const Double_t binerror       = hnew->GetBinError(i);
            firstbincontent              += bincontent;
            firstbinerror2               += (binerror * binerror);
            firstbinerrorf                = (firstbincontent > 1e-6) ? (TMath::Sqrt(firstbinerror2) / firstbincontent) : 1.0;
            firstN ++;
        }
    }
    // Group bins from the right
    if (lastN<1){
        Double_t lastbincontent   = hnew->GetBinContent(end-1);
        Double_t lastbinerror     = hnew->GetBinError(end-1);
        Double_t lastbinerror2    = (lastbinerror * lastbinerror);
        Double_t lastbinerrorf    = (lastbincontent > 1e-6) ? (TMath::Sqrt(lastbinerror2) / lastbincontent) : 1.0;
        lastN = 1; // starts with 1 bin
        // Find lastN such that its error fraction becomes less than the threshold
        while ((lastbinerrorf > threshold) && (lastN <= nbins)){
            //hdummy->Info("Rebin", "lastN=%2i int. binc=%f int. bine=%f frac=%f", lastN, lastbincontent, sqrt(lastbinerror2), lastbinerrorf);
                    
            Int_t i = end-1 - lastN;
            const Double_t bincontent     = hnew->GetBinContent(i);
            const Double_t binerror       = hnew->GetBinError(i);
            lastbincontent               += bincontent;
            lastbinerror2                += (binerror * binerror);
            lastbinerrorf                 = (lastbincontent > 1e-6) ? (TMath::Sqrt(lastbinerror2) / lastbincontent) : 1.0;
            lastN ++;
        }
    }
    // Group bins in the middle
    Int_t nmiddlebins = nbins - firstN - lastN;
    Int_t ngroup = newnbins-2;
    Int_t nbg = nmiddlebins/ngroup;
    //hdummy->Info("Rebin", "firstN=%2i lastN=%2i nmiddlebins=%i before adjusting nmiddlebins", firstN, lastN, nmiddlebins);
    if (nmiddlebins < (newnbins-2)){
        hdummy->Error("Rebin", "Not enough number of bins! Perhaps the threshold is too low or newnbins is too many. threshold=%f, newnbins=%i", threshold, newnbins);
        return hdummy;
    } else if (!isfixed){
        // Group more bins from the left until nmiddlebins can be divided by ngroup
        while ((nbg*ngroup != nmiddlebins)){
            firstN++;  
            nmiddlebins = nbins - firstN - lastN;
            nbg = nmiddlebins/ngroup;
        }
    }
    //hdummy->Info("Rebin", "firstN=%2i lastN=%2i nmiddlebins=%i  after adjusting nmiddlebins", firstN, lastN, nmiddlebins);
    const TAxis* a = h->GetXaxis();
    TArrayD newlowedges(newnbins+1);                        // to make a histogram of newnbins, need newnbins+1 low edges
    newlowedges.SetAt(a->GetBinLowEdge(begin), 0);          // bin 0 low edge
    newlowedges.SetAt(a->GetBinLowEdge(end)  , newnbins);   // bin newnbins low edge
    for (Int_t i=1; i<newnbins; i++){
        UInt_t ibin = begin + firstN + (i-1)*nbg;
        newlowedges.SetAt(a->GetBinLowEdge(ibin), i);       // bin 1 ... newnbins-1 low edge
        if ((i == newnbins-1) && (ibin != (end - lastN))){
            hdummy->Error("Rebin", "Bin edges don't match! Perhaps the number of bins of input histogram has changed? firstN=%i, lastN=%i", firstN, lastN);
            return hdummy;
        }
    }
    //for (Int_t i=0; i<newnbins+1; i++)
    //    hdummy->Info("Rebin", "edge %2i=%f", i, newlowedges.GetAt(i));
    TH1F * hnewrebin = (TH1F *) hnew->Rebin(newnbins, Form("%s_test2", newname), newlowedges.GetArray());  // temp histogram with the new binning of variable size
    hnewrebin->Sumw2();
    // Sanity checks
    const TAxis* anew = hnewrebin->GetXaxis();
    if (! TMath::AreEqualRel(a->GetXmin(), anew->GetXmin(),1.E-10) ||
        ! TMath::AreEqualRel(a->GetXmax(), anew->GetXmax(),1.E-10) ) {
        hdummy->Error("Rebin", "Axis limits have changed! Please debug! h=%s", h->GetName());
        return hdummy;
    }
    if (TMath::AreEqualRel(newxlow, newxup, 1E-10)) {  // if newxlow == newxup, keep the boundaries
        newxlow = anew->GetXmin();
        newxup = anew->GetXmax();
    }
    TString newtitle = Form("%s;%s;%s", h->GetTitle(), h->GetXaxis()->GetTitle(), h->GetYaxis()->GetTitle());
    TH1F * hnewrebinfixed = new TH1F(newname, newtitle, newnbins, newxlow, newxup);  // final histogram with the new binning of fixed size
    hnewrebinfixed->Sumw2();
    for (Int_t i=0; i<newnbins+2; i++){  // include UOflow
        hnewrebinfixed->SetBinContent(i, hnewrebin->GetBinContent(i));
        hnewrebinfixed->SetBinError  (i, hnewrebin->GetBinError(i));
    }
    // Sanity checks
    if (! TMath::AreEqualRel(h->GetSumOfWeights(), hnewrebinfixed->GetSumOfWeights(), 1E-7)){  // they can differ slightly
        hdummy->Error("Rebin", "Sum of weights has changed! Please debug! h=%s, %f --> %f", h->GetName(), h->GetSumOfWeights(), hnewrebinfixed->GetSumOfWeights());
        return hdummy;
    }
    Double_t oldsumoferrors2 = 0., newsumoferrors2 = 0.;
    for (Int_t i=0; i<nbins+2; i++){  // include UOflow
        oldsumoferrors2 += (h->GetBinError(i) * h->GetBinError(i));
        if (i<newnbins+2)
            newsumoferrors2 += (hnewrebinfixed->GetBinError(i) * hnewrebinfixed->GetBinError(i));
    }
    if (! TMath::AreEqualRel(oldsumoferrors2, newsumoferrors2, 1E-10)){
        hdummy->Error("Rebin", "Sum of errors squared has changed! Please debug! h=%s, %f --> %f", h->GetName(), oldsumoferrors2, newsumoferrors2);
        return hdummy;
    }
    hdummy->Info("Rebin", "Transform %i bins into %i bins {|%i|,|%i|x%i,|%i|}", nbins, newnbins, firstN, nbg, ngroup, lastN);
    double firstbinerrorf = (hnewrebinfixed->GetBinError(1) / hnewrebinfixed->GetBinContent(1));
    double lastbinerrorf  = (hnewrebinfixed->GetBinError(newnbins) / hnewrebinfixed->GetBinContent(newnbins));
    //hdummy->Info("Rebin", "firstbinerrorf=%f, lastbinerrorf=%f", firstbinerrorf, lastbinerrorf);
    if (firstbinerrorf < 0 || lastbinerrorf < 0)
        hdummy->Error("Rebin", "Insensible error fractions! Please debug! h=%s, firstbinerrorf=%f, lastbinerrorf=%f", h->GetName(), firstbinerrorf, lastbinerrorf);
    delete hdummy;
    delete hnew;
    delete hnewrebin;
    return hnewrebinfixed;
}

// Overloaded for using background-specific BDT
TH1F * RebinByStatErrors(const TH1F * h, Int_t newnbins, Int_t& firstN1, Int_t& firstN2, Int_t& firstN3, Int_t& firstN4, Int_t& lastN1, Int_t& lastN2, Int_t& lastN3, Int_t& lastN4, Double_t threshold=0.35, Double_t newxlow=0., Double_t newxup=0., const char * newname="") {
    // newnbins is the # of bins in each individual partition.
    TH1F * hdummy = (TH1F *) h->Clone("dummy_split");
    hdummy->Sumw2();
    Int_t nbins = h->GetNbinsX();
    if (h->GetXaxis()->IsVariableBinSize()){
        hdummy->Error("Rebin", "Input histogram has variable bin sizes! h=%s", h->GetName());
        return hdummy;
    }
    if (nbins%4 != 0){
        hdummy->Error("Rebin", "nbins must be divisible by 4! nbins=%i", nbins);
        return hdummy;
    }
    if (nbins/4 == newnbins){
        hdummy->Warning("Rebin", "newnbins is equal to current nbins/4. Return a clone. nbins=%i/4, newnbins=%i", nbins, newnbins);
        hdummy->SetName(newname);
        return hdummy;
    }
    if (newnbins < 3 || newnbins > nbins/4){
        hdummy->Error("Rebin", "newnbins must be within [3,%i/4]. newnbins=%i", nbins, newnbins);
        return hdummy;
    }
    if (firstN1 > nbins/4 || firstN2 > nbins/4 || firstN3 > nbins/4 || firstN4 > nbins/4 || 
        lastN1 > nbins/4 || lastN2 > nbins/4 || lastN3 > nbins/4 || lastN4 > nbins/4){
        hdummy->Error("Rebin", "firstN and lastN must be < %i/4. firstN1=%i, firstN2=%i, firstN3=%i, firstN4=%i, lastN1=%i, lastN2=%i, lastN3=%i, lastN4=%i", nbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4);
        return hdummy;
    }
    TString newtitle = Form("%s;%s;%s", h->GetTitle(), h->GetXaxis()->GetTitle(), h->GetYaxis()->GetTitle());
    TH1F * hsplit11 = new TH1F(Form("%s_split11", newname), newtitle, nbins/4, newxlow, newxup);
    TH1F * hsplit12 = new TH1F(Form("%s_split12", newname), newtitle, nbins/4, newxlow, newxup);
    TH1F * hsplit13 = new TH1F(Form("%s_split13", newname), newtitle, nbins/4, newxlow, newxup);
    TH1F * hsplit14 = new TH1F(Form("%s_split14", newname), newtitle, nbins/4, newxlow, newxup);
    if (h->GetSumw2N() == 0){
        hdummy->Warning("Rebin", "Sumw2 is not created prior to this. h=%s", h->GetName());
        hsplit11->Sumw2();
        hsplit12->Sumw2();
        hsplit13->Sumw2();
        hsplit14->Sumw2();
    }
    
    for (Int_t i = 1; i < (1 + nbins); i++){  // ignore UOflow
        TH1F * hh;
        Int_t ii = i;
        if (i < (1 + nbins*1/4)) {
            hh = hsplit11;
            ii = i - nbins*0/4;
        } else if (i < (1 + nbins*2/4)) {
            hh = hsplit12;
            ii = i - nbins*1/4;
        } else if (i < (1 + nbins*3/4)) {
            hh = hsplit13;
            ii = i - nbins*2/4;
        } else if (i < (1 + nbins*4/4)) {
            hh = hsplit14;
            ii = i - nbins*3/4;
        }
        hh->SetBinContent(ii, h->GetBinContent(i));
        hh->SetBinError(ii, h->GetBinError(i));
    }

    TH1F * hsplit1 = RebinByStatErrors(hsplit11, newnbins, firstN1, lastN1, threshold, newxlow, newxup, Form("%s_split1", newname));
    TH1F * hsplit2 = RebinByStatErrors(hsplit12, newnbins, firstN2, lastN2, threshold, newxlow, newxup, Form("%s_split2", newname));
    TH1F * hsplit3 = RebinByStatErrors(hsplit13, newnbins, firstN3, lastN3, threshold, newxlow, newxup, Form("%s_split3", newname));
    TH1F * hsplit4 = RebinByStatErrors(hsplit14, newnbins, firstN4, lastN4, threshold, newxlow, newxup, Form("%s_split4", newname));
    
    TH1F * hnewrebinfixed = new TH1F(newname, newtitle, newnbins*4, newxlow, newxup);  // final histogram with the new binning of fixed size
    hnewrebinfixed->Sumw2();
    for (Int_t i = 1; i < (1 + newnbins); i++){  // ignore UOflow
        hnewrebinfixed->SetBinContent(i + newnbins*0, hsplit1->GetBinContent(i));
        hnewrebinfixed->SetBinContent(i + newnbins*1, hsplit2->GetBinContent(i));
        hnewrebinfixed->SetBinContent(i + newnbins*2, hsplit3->GetBinContent(i));
        hnewrebinfixed->SetBinContent(i + newnbins*3, hsplit4->GetBinContent(i));
        hnewrebinfixed->SetBinError(i + newnbins*0, hsplit1->GetBinError(i));
        hnewrebinfixed->SetBinError(i + newnbins*1, hsplit2->GetBinError(i));
        hnewrebinfixed->SetBinError(i + newnbins*2, hsplit3->GetBinError(i));
        hnewrebinfixed->SetBinError(i + newnbins*3, hsplit4->GetBinError(i));
    }
    
    // Sanity checks
    if (! TMath::AreEqualRel(h->GetSumOfWeights(), hnewrebinfixed->GetSumOfWeights(), 1E-7)){  // they can differ slightly
        hdummy->Error("Rebin", "Sum of weights has changed! Please debug! h=%s, %f --> %f", h->GetName(), h->GetSumOfWeights(), hnewrebinfixed->GetSumOfWeights());
        return hdummy;
    }
    Double_t oldsumoferrors2 = 0., newsumoferrors2 = 0.;
    for (Int_t i=0; i<nbins+2; i++){  // include UOflow
        oldsumoferrors2 += (h->GetBinError(i) * h->GetBinError(i));
    }
    for (Int_t i=0; i<(newnbins*4)+2; i++){  // include UOflow
        newsumoferrors2 += (hnewrebinfixed->GetBinError(i) * hnewrebinfixed->GetBinError(i));
    }
    if (! TMath::AreEqualRel(oldsumoferrors2, newsumoferrors2, 1E-10)){
        hdummy->Error("Rebin", "Sum of errors squared has changed! Please debug! h=%s, %f --> %f", h->GetName(), oldsumoferrors2, newsumoferrors2);
        return hdummy;
    }
    
    delete hdummy;
    delete hsplit11;
    delete hsplit12;
    delete hsplit13;
    delete hsplit14;
    delete hsplit1;
    delete hsplit2;
    delete hsplit3;
    delete hsplit4;
    return hnewrebinfixed;
}

// Check consistency between two histograms
bool CheckConsistency(const TH1 * h1, const TH1 * h2)
{
    if (h1->GetNbinsX() != h2->GetNbinsX()){
        h1->Error("CheckConsistency", "Two input histograms have different number of bins!");
        return false;
    }
    const TAxis * a1 = h1->GetXaxis();
    const TAxis * a2 = h2->GetXaxis();
    if (! TMath::AreEqualRel(a1->GetXmin(), a2->GetXmin(),1.E-12) ||
        ! TMath::AreEqualRel(a1->GetXmax(), a2->GetXmax(),1.E-12) ) {
        h1->Error("CheckConsistency", "Two input histograms have different axis limits!");
        return false;
    }
    const TArrayD * h1Array = a1->GetXbins(); 
    const TArrayD * h2Array = a2->GetXbins(); 
    Int_t fN = h1Array->fN;
    if ( fN != 0 ) {
        if ( h2Array->fN != fN ) {
            h1->Error("CheckConsistency", "Two input histograms have different bin limits!");
            return false;
        }
        else {
            for ( int i = 0; i < fN; ++i ) {
                if ( ! TMath::AreEqualRel( h1Array->GetAt(i), h2Array->GetAt(i), 1E-10 ) ) {
                    h1->Error("CheckConsistency", "Two input histograms have different bin limits!");
                    return false;
                }
            }
        }
    }
    return true;
}

TH1F * VaryStatErrors(const TH1F * h, const Double_t c, const char * newname="", bool doUOflow=false)
{
    TH1F * hdummy = (TH1F *) h->Clone("dummy");
    hdummy->Sumw2();
    Int_t nbins = h->GetNbinsX();
    if (TMath::Abs(c*c - 1.0) > 1e-10){
        hdummy->Error("Vary", "c must be +1 or -1!");
        return hdummy;
    }
    TH1F * hnew = (TH1F *) h->Clone(newname);
    if (h->GetSumw2N() == 0){
        hdummy->Warning("Vary", "Sumw2 is not created prior to this.");
        hnew->Sumw2();
    }
    if (h->GetSumOfWeights() < 1e-10){
        hdummy->Warning("Vary", "Empty histogram. Return a clone. h=%s", h->GetName());
        return hnew;
    }
    UInt_t begin=1, end=nbins+1;
    if (doUOflow){  // do underflow and overflow
        begin=0; end=nbins+2;
    } else if ((h->GetBinContent(0)>1e-6) || (h->GetBinContent(nbins+1)>1e-6)){
        hdummy->Warning("Vary", "Underflow and/or overflow are not empty and not modified.");
    }
    Double_t sumoferrors2 = 0.;
    for (UInt_t i=begin; i<end; i++){
        const Double_t bincontent     = h->GetBinContent(i);
        const Double_t binerror       = h->GetBinError(i);
        const Double_t newbincontent  = TMath::Max(0., bincontent + (binerror * c));
        // delta error = sqrt(old error)
        // new error = old error**2 +/- delta error**2
        //           = old error**2 +/- old error
        const Double_t newbinerror2   = TMath::Max(0., (binerror * binerror) + (binerror * c));
        const Double_t newbinerror    = TMath::Sqrt(newbinerror2);
        sumoferrors2                 += (binerror * binerror);

        hnew->SetBinContent(i, newbincontent);
        hnew->SetBinError  (i, newbinerror);
    }
    // new normalization = old normalization +/- sqrt of sum of old stat errors
    Double_t newsumofweights = h->GetSumOfWeights() + (TMath::Sqrt(sumoferrors2) * c);
    if (newsumofweights < 0.){
        hdummy->Error("Vary", "Negative new sum of weights!");
        return hdummy;
    }
    hnew->Scale(newsumofweights / (hnew->GetSumOfWeights()));
    delete hdummy;
    return hnew;
}

// FIXME: remove fabs
TH1F * VaryModelErrors(const TH1F * h, const TH1F * hmodel, Double_t c, const char* newname="", bool doUOflow=false)
{
    TH1F * hdummy = (TH1F *) h->Clone("dummy");
    hdummy->Sumw2();
    Int_t nbins = h->GetNbinsX();
    if (TMath::Abs(c*c - 1.0) > 1e-10){
        hdummy->Error("Vary", "c must be +1 or -1!");  // allow c to be any constant?
        return hdummy;
    }
    if (h == hmodel){
        hdummy->Error("Vary", "The input histograms are the same!");
        return hdummy;
    } else {
        if(!CheckConsistency(h, hmodel)) {
            hdummy->Error("Vary", "The input histograms are not consistent! h=%s", h->GetName());
            return hdummy;
        }
    }
    TH1F * hnew = (TH1F *) h->Clone(newname);
    if (h->GetSumw2N() == 0){
        hdummy->Warning("Vary", "Sumw2 is not created prior to this.");
        hnew->Sumw2();
    }
    if (h->GetSumOfWeights() < 1e-10){
        hdummy->Warning("Vary", "Empty histogram. Return a clone. h=%s", h->GetName());
        return hnew;
    }
    UInt_t begin=1, end=nbins+1;
    if (doUOflow){  // do underflow and overflow
        begin=0; end=nbins+2;
    } else if ((h->GetBinContent(0)>1e-6) || (h->GetBinContent(nbins+1)>1e-6)){
        hdummy->Warning("Vary", "Underflow and/or overflow are not empty and not modified.");
    }
    // Normalize hmodel to equal area
    TH1F * hmodelnorm = (TH1F *) hmodel->Clone(Form("%s_norm", hmodel->GetName()));
    if (hmodel->GetSumw2N() == 0){
        hdummy->Warning("Vary", "Sumw2 is not created prior to this.");
        hmodelnorm->Sumw2();
    }
    hmodelnorm->Scale(h->GetSumOfWeights() / (hmodelnorm->GetSumOfWeights()));
    for (UInt_t i=begin; i<end; i++){
        const Double_t bincontent     = h->GetBinContent(i);
        const Double_t binerror       = h->GetBinError(i);
        const Double_t bincontentm    = hmodelnorm->GetBinContent(i);
        const Double_t binerrorm      = hmodelnorm->GetBinError(i);
        const Double_t diff           = fabs(bincontentm - bincontent);
        const Double_t newbincontent  = TMath::Max(0., bincontent + (0.5 * diff * c));
        // the operation is h1 +/- 0.5 * (h2 - h1)
        // the bin error is simply sqrt(((1.0 +/- 0.5)**2 * e1**2) + (0.5**2 * e2**2))
        const Double_t newbinerror2   = TMath::Max(0., ((1.0+(0.5*c))*(1.0+(0.5*c)) * (binerror * binerror)) + (0.5*0.5 * (binerrorm * binerrorm)) );
        const Double_t newbinerror    = TMath::Sqrt(newbinerror2);

        hnew->SetBinContent(i, newbincontent);
        hnew->SetBinError  (i, newbinerror);
    }
    // Normalize hnew to equal area
    hnew->Scale(h->GetSumOfWeights() / (hnew->GetSumOfWeights()));
    delete hdummy;
    return hnew;
}

///_____________________________________________________________________________
/// The core function that makes plots, prints stats, and writes the datacard.
/// Note: QCD histograms are empty.

using namespace std;

/// Initialize rebinning "bookmarks" to negative values
/// These global variables that need to be reset when we switch channel.
Int_t firstN = -1, lastN = -1;  // for single BDT
Int_t firstN1 = -1, firstN2 = -1, firstN3 = -1, firstN4 = -1; // for PQRS BDT
Int_t lastN1  = -1, lastN2  = -1, lastN3  = -1, lastN4  = -1; // for PQRS BDT

void MakePlots(const EventsJ10 * ev, const TString var, 
               const TCut cutmc, const TCut cutdata, const TString syst, const TString region, 
               const char* axtitle, Int_t nbins, Double_t xlow, Double_t xup, 
               Int_t newnbins, Double_t rebinerrorf, 
               const Double_t scalefactors_lnN[5], 
               TString options="printStat:printCard:plotData:plotLog:plotSig")
{
    std::clog << "MakePlots(): Using cutmc = " << cutmc << ", cutdata = " << cutdata << std::endl;

    // Parse options
    options.ReplaceAll(" ", "");
    bool printStat = options.Contains("printStat") && (!options.Contains("!printStat"));
    bool printCard = options.Contains("printCard") && (!options.Contains("!printCard"));
    bool writeRoot = options.Contains("writeRoot") && (!options.Contains("!writeRoot"));
    bool plotData  = options.Contains("plotData")  && (!options.Contains("!plotData"));
    bool plotSig   = options.Contains("plotSig")   && (!options.Contains("!plotSig"));
    bool plotLog   = options.Contains("plotLog")   && (!options.Contains("!plotLog"));

    // Book histograms before rebinning
    TH1F * hZH_0        = new TH1F("ZH_0"      , "", nbins, xlow, xup);
    TH1F * hWH_0        = new TH1F("WH_0"      , "", nbins, xlow, xup);
    TH1F * hWjLF_0      = new TH1F("WjLF_0"    , "", nbins, xlow, xup);
    TH1F * hWjHF_0      = new TH1F("WjHF_0"    , "", nbins, xlow, xup);
    TH1F * hZjLF_0      = new TH1F("ZjLF_0"    , "", nbins, xlow, xup);
    TH1F * hZjHF_0      = new TH1F("ZjHF_0"    , "", nbins, xlow, xup);
    TH1F * hTT_0        = new TH1F("TT_0"      , "", nbins, xlow, xup);
    TH1F * hs_Top_0     = new TH1F("s_Top_0"   , "", nbins, xlow, xup);
    TH1F * hVVLF_0      = new TH1F("VVLF_0"    , "", nbins, xlow, xup);
    TH1F * hVVHF_0      = new TH1F("VVHF_0"    , "", nbins, xlow, xup);
    TH1F * hQCD_0       = new TH1F("QCD_0"     , "", nbins, xlow, xup);
    TH1F * hZjLFHW_0    = new TH1F("ZjLFHW_0"  , "", nbins, xlow, xup);
    TH1F * hZjHFHW_0    = new TH1F("ZjHFHW_0"  , "", nbins, xlow, xup);
    TH1F * hTTMG_0      = new TH1F("TTMG_0"    , "", nbins, xlow, xup);
    TH1F * hZH_SM_0     = new TH1F("ZH_SM_0"   , "", nbins, xlow, xup);
    TH1F * hWH_SM_0     = new TH1F("WH_SM_0"   , "", nbins, xlow, xup);
    TH1F * hVH_0        = new TH1F("VH_0"      , "", nbins, xlow, xup);
    TH1F * hmc_exp_0    = new TH1F("mc_exp_0"  , "", nbins, xlow, xup);
    TH1F * hdata_obs_0  = new TH1F("data_obs_0", "", nbins, xlow, xup);
    
    std::vector<TH1 *> histos_0;
    histos_0.push_back(hZH_0);
    histos_0.push_back(hWH_0);
    histos_0.push_back(hWjLF_0);
    histos_0.push_back(hWjHF_0);
    histos_0.push_back(hZjLF_0);
    histos_0.push_back(hZjHF_0);
    histos_0.push_back(hTT_0);
    histos_0.push_back(hs_Top_0);
    histos_0.push_back(hVVLF_0);
    histos_0.push_back(hVVHF_0);
    histos_0.push_back(hQCD_0);
    histos_0.push_back(hZjLFHW_0);
    histos_0.push_back(hZjHFHW_0);
    histos_0.push_back(hTTMG_0);
    histos_0.push_back(hZH_SM_0);
    histos_0.push_back(hWH_SM_0);
    histos_0.push_back(hVH_0);
    histos_0.push_back(hmc_exp_0);
    histos_0.push_back(hdata_obs_0);
    
    for (UInt_t ih = 0; ih < histos_0.size(); ih++)
        histos_0.at(ih)->Sumw2();
    
    // Project
    ev->ZH->Project("ZH_0", var, cutmc);
    std::clog << "... DONE: project ZH_0." << std::endl;
    ev->WH->Project("WH_0", var, cutmc);
    std::clog << "... DONE: project WH_0." << std::endl;
    
    // Apply scale factors
    TCut cutsf = "1.0";
    cutsf = Form("%f", ev->sf_WjLF);
    ev->WjLF->Project("WjLF_0", var, cutmc * cutsf);
    std::clog << "... DONE: project WjLF_0, scaled by " << cutsf << "." << std::endl;
    cutsf = Form("%f", ev->sf_WjHF);
    ev->WjHF->Project("WjHF_0", var, cutmc * cutsf);
    std::clog << "... DONE: project WjHF_0, scaled by " << cutsf << "." << std::endl;
    cutsf = Form("%f", ev->sf_ZjLF);
    ev->ZjLF->Project("ZjLF_0", var, cutmc * cutsf);
    std::clog << "... DONE: project ZjLF_0, scaled by " << cutsf << "." << std::endl;
    cutsf = Form("%f", ev->sf_ZjHF);
    ev->ZjHF->Project("ZjHF_0", var, cutmc * cutsf);
    std::clog << "... DONE: project ZjHF_0, scaled by " << cutsf << "." << std::endl;
    cutsf = Form("%f", ev->sf_TT);
    ev->TT->Project("TT_0", var, cutmc * cutsf);
    std::clog << "... DONE: project TT_0, scaled by " << cutsf << "." << std::endl;
    
    ev->s_Top->Project("s_Top_0", var, cutmc);
    std::clog << "... DONE: project s_Top_0." << std::endl;
    ev->VVLF->Project("VVLF_0", var, cutmc);
    std::clog << "... DONE: project VVLF_0." << std::endl;
    ev->VVHF->Project("VVHF_0", var, cutmc);
    std::clog << "... DONE: project VVHF_0." << std::endl;
    //ev->QCD->Project("QCD_0", var, cutmc);
    //std::clog << "... DONE: project QCD_0." << std::endl;
    
    ev->ZjLFHW->Project("ZjLFHW_0", var, cutmc);
    std::clog << "... DONE: project ZjLFHW_0." << std::endl;
    ev->ZjHFHW->Project("ZjHFHW_0", var, cutmc);
    std::clog << "... DONE: project ZjHFHW_0." << std::endl;
    ev->TTMG->Project("TTMG_0", var, cutmc);
    std::clog << "... DONE: project TTMG_0." << std::endl;
    
    ev->ZH_SM->Project("ZH_SM_0", var, cutmc);
    std::clog << "... DONE: project ZH_SM_0." << std::endl;
    ev->WH_SM->Project("WH_SM_0", var, cutmc);
    std::clog << "... DONE: project WH_SM_0." << std::endl;
    
    ev->data_obs->Project("data_obs_0", var, cutdata);
    std::clog << "... DONE: project data_obs_0." << std::endl;

    // Make sum of histograms
    hVH_0->Add(hZH_0);
    hVH_0->Add(hWH_0);
    std::clog << "... DONE: add MC signals to VH_0." << std::endl;

    hmc_exp_0->Add(hWjLF_0);
    hmc_exp_0->Add(hWjHF_0);
    hmc_exp_0->Add(hZjLF_0);
    hmc_exp_0->Add(hZjHF_0);
    hmc_exp_0->Add(hTT_0);
    hmc_exp_0->Add(hs_Top_0);
    hmc_exp_0->Add(hVVLF_0);
    //hmc_exp_0->Add(hVVHF_0);  // FIXME: add VH as background?
    //hmc_exp_0->Add(hQCD_0);
    std::clog << "... DONE: add MC backgrounds to mc_exp_0." << std::endl;

    // Do the rebinning
    TH1F * hZH, * hWH, * hWjLF, * hWjHF, * hZjLF, * hZjHF, * hTT, * hs_Top, * hVVLF, * hVVHF, * hQCD,
         * hZjLFHW, * hZjHFHW, * hTTMG, * hZH_SM, * hWH_SM, * hVH, * hmc_exp, * hdata_obs;

    if (!(var.Contains("slice"))) {
        hmc_exp   = RebinByStatErrors(hmc_exp_0  , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "mc_exp");
        hZH       = RebinByStatErrors(hZH_0      , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "ZH");
        hWH       = RebinByStatErrors(hWH_0      , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "WH");
        hWjLF     = RebinByStatErrors(hWjLF_0    , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "WjLF");
        hWjHF     = RebinByStatErrors(hWjHF_0    , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "WjHF");
        hZjLF     = RebinByStatErrors(hZjLF_0    , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "ZjLF");
        hZjHF     = RebinByStatErrors(hZjHF_0    , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "ZjHF");
        hTT       = RebinByStatErrors(hTT_0      , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "TT");
        hs_Top    = RebinByStatErrors(hs_Top_0   , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "s_Top");
        hVVLF     = RebinByStatErrors(hVVLF_0    , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "VVLF");
        hVVHF     = RebinByStatErrors(hVVHF_0    , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "VVHF");
        hQCD      = RebinByStatErrors(hQCD_0     , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "QCD");
        hZjLFHW   = RebinByStatErrors(hZjLFHW_0  , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "ZjLFHW");
        hZjHFHW   = RebinByStatErrors(hZjHFHW_0  , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "ZjHFHW");
        hTTMG     = RebinByStatErrors(hTTMG_0    , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "TTMG");
        hZH_SM    = RebinByStatErrors(hZH_SM_0   , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "ZH_SM");
        hWH_SM    = RebinByStatErrors(hWH_SM_0   , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "WH_SM");
        hVH       = RebinByStatErrors(hVH_0      , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "VH");
        hdata_obs = RebinByStatErrors(hdata_obs_0, newnbins, firstN, lastN, rebinerrorf, xlow, xup, "data_obs");
    } else {
        // newnbins decides the binning in each slice (as opposed to the total).
        hmc_exp   = RebinByStatErrors(hmc_exp_0  , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "mc_exp");
        hZH       = RebinByStatErrors(hZH_0      , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "ZH");
        hWH       = RebinByStatErrors(hWH_0      , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "WH");
        hWjLF     = RebinByStatErrors(hWjLF_0    , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "WjLF");
        hWjHF     = RebinByStatErrors(hWjHF_0    , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "WjHF");
        hZjLF     = RebinByStatErrors(hZjLF_0    , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "ZjLF");
        hZjHF     = RebinByStatErrors(hZjHF_0    , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "ZjHF");
        hTT       = RebinByStatErrors(hTT_0      , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "TT");
        hs_Top    = RebinByStatErrors(hs_Top_0   , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "s_Top");
        hVVLF     = RebinByStatErrors(hVVLF_0    , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "VVLF");
        hVVHF     = RebinByStatErrors(hVVHF_0    , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "VVHF");
        hQCD      = RebinByStatErrors(hQCD_0     , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "QCD");
        hZjLFHW   = RebinByStatErrors(hZjLFHW_0  , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "ZjLFHW");
        hZjHFHW   = RebinByStatErrors(hZjHFHW_0  , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "ZjHFHW");
        hTTMG     = RebinByStatErrors(hTTMG_0    , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "TTMG");
        hZH_SM    = RebinByStatErrors(hZH_SM_0   , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "ZH_SM");
        hWH_SM    = RebinByStatErrors(hWH_SM_0   , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "WH_SM");
        hVH       = RebinByStatErrors(hVH_0      , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "VH");
        hdata_obs = RebinByStatErrors(hdata_obs_0, newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "data_obs");
    }
    std::vector<TH1 *> histos;
    histos.push_back(hZH);
    histos.push_back(hWH);
    histos.push_back(hWjLF);
    histos.push_back(hWjHF);
    histos.push_back(hZjLF);
    histos.push_back(hZjHF);
    histos.push_back(hTT);
    histos.push_back(hs_Top);
    histos.push_back(hVVLF);
    histos.push_back(hVVHF);
    histos.push_back(hQCD);
    histos.push_back(hZjLFHW);
    histos.push_back(hZjHFHW);
    histos.push_back(hTTMG);
    histos.push_back(hZH_SM);
    histos.push_back(hWH_SM);
    histos.push_back(hVH);
    histos.push_back(hmc_exp);
    histos.push_back(hdata_obs);
    
    assert(histos.size() == histos_0.size());
    for (UInt_t ih = 0; ih < histos.size(); ih++)
        histos.at(ih)->Sumw2();
    
    if (printStat) {
        std::clog << "MakePlots(): Printing statistics..." << std::endl;
        
        const int p_bin1=0, p_bin2=9999;
        double p_integral=0., p_error=0.;
        double D=0., S=0., B=0.;
        std::cout.setf(ios::fixed,ios::floatfield);
        std::cout.precision(3);
        
        p_integral = hZH->IntegralAndError(p_bin1, p_bin2, p_error);
        S += p_integral;
        std::cout << setw(10) << left << "ZH " << "(" << p_bin1 << ", " << p_bin2 << ") = " << setw(8) << right << p_integral << " +/- " << setw(6) << right << p_error << std::endl;
        p_integral = hWH->IntegralAndError(p_bin1, p_bin2, p_error);
        S += p_integral;
        std::cout << setw(10) << left << "WH " << "(" << p_bin1 << ", " << p_bin2 << ") = " << setw(8) << right << p_integral << " +/- " << setw(6) << right << p_error << std::endl;
        
        p_integral = hWjLF->IntegralAndError(p_bin1, p_bin2, p_error);
        B += p_integral;
        std::cout << setw(10) << left << "WjLF " << "(" << p_bin1 << ", " << p_bin2 << ") = " << setw(8) << right << p_integral << " +/- " << setw(6) << right << p_error << std::endl;
        p_integral = hWjHF->IntegralAndError(p_bin1, p_bin2, p_error);
        B += p_integral;
        std::cout << setw(10) << left << "WjHF " << "(" << p_bin1 << ", " << p_bin2 << ") = " << setw(8) << right << p_integral << " +/- " << setw(6) << right << p_error << std::endl;
        p_integral = hZjLF->IntegralAndError(p_bin1, p_bin2, p_error);
        B += p_integral;
        std::cout << setw(10) << left << "ZjLF " << "(" << p_bin1 << ", " << p_bin2 << ") = " << setw(8) << right << p_integral << " +/- " << setw(6) << right << p_error << std::endl;
        p_integral = hZjHF->IntegralAndError(p_bin1, p_bin2, p_error);
        B += p_integral;
        std::cout << setw(10) << left << "ZjHF " << "(" << p_bin1 << ", " << p_bin2 << ") = " << setw(8) << right << p_integral << " +/- " << setw(6) << right << p_error << std::endl;
        p_integral = hTT->IntegralAndError(p_bin1, p_bin2, p_error);
        B += p_integral;
        std::cout << setw(10) << left << "TT " << "(" << p_bin1 << ", " << p_bin2 << ") = " << setw(8) << right << p_integral << " +/- " << setw(6) << right << p_error << std::endl;
        p_integral = hs_Top->IntegralAndError(p_bin1, p_bin2, p_error);
        B += p_integral;
        std::cout << setw(10) << left << "s_Top " << "(" << p_bin1 << ", " << p_bin2 << ") = " << setw(8) << right << p_integral << " +/- " << setw(6) << right << p_error << std::endl;
        p_integral = hVVLF->IntegralAndError(p_bin1, p_bin2, p_error);
        B += p_integral;
        std::cout << setw(10) << left << "VVLF " << "(" << p_bin1 << ", " << p_bin2 << ") = " << setw(8) << right << p_integral << " +/- " << setw(6) << right << p_error << std::endl;
        p_integral = hVVHF->IntegralAndError(p_bin1, p_bin2, p_error);
        B += p_integral;
        std::cout << setw(10) << left << "VVHF " << "(" << p_bin1 << ", " << p_bin2 << ") = " << setw(8) << right << p_integral << " +/- " << setw(6) << right << p_error << std::endl;
        p_integral = hQCD->IntegralAndError(p_bin1, p_bin2, p_error);
        B += p_integral;
        std::cout << setw(10) << left << "QCD " << "(" << p_bin1 << ", " << p_bin2 << ") = " << setw(8) << right << p_integral << " +/- " << setw(6) << right << p_error << std::endl;
        
        p_integral = hVH->IntegralAndError(p_bin1, p_bin2, p_error);
        std::cout << setw(10) << left << "VH " << "(" << p_bin1 << ", " << p_bin2 << ") = " << setw(8) << right << p_integral << " +/- " << setw(6) << right << p_error << std::endl;
        p_integral = hmc_exp->IntegralAndError(p_bin1, p_bin2, p_error);
        std::cout << setw(10) << left << "mc_exp " << "(" << p_bin1 << ", " << p_bin2 << ") = " << setw(8) << right << p_integral << " +/- " << setw(6) << right << p_error << std::endl;
        p_integral = hdata_obs->IntegralAndError(p_bin1, p_bin2, p_error);
        D += p_integral;
        std::cout << setw(10) << left << "data_obs " << "(" << p_bin1 << ", " << p_bin2 << ") = " << setw(8) << right << p_integral << " +/- " << setw(6) << right << p_error << std::endl;

        double Sover1p5B = S / (1.5 + sqrt(B) + (0.2 * B));
        double SoverSB = S / (sqrt(S + B));

        std::cout.precision(5);
        std::cout << "==> data, S, B, S/(1.5 + sqrt(B) + (0.2 * B)), S/sqrt(S+B): " << D << ", " << S << ", " << B << ", " << Sover1p5B << ", " << SoverSB << std::endl;
    }

    if (printCard && syst == "NONE") {
        std::clog << "MakePlots(): Printing the datacard..." << std::endl;

        TString dcname = g_dcname;
        TString wsname = g_wsname;
        dcname.ReplaceAll("$CHANNEL", channel);
        wsname.ReplaceAll("$CHANNEL", channel);
        const TString channel_8TeV = channel + "_8TeV";
        TString channel_8TeV_1     = channel_8TeV;
        if (channel == "ZnunuLowCSV")
            channel_8TeV_1 = "ZnunuHighPt_8TeV";
        
        std::ofstream dc;
        dc.setf(ios::fixed,ios::floatfield);
        dc.precision(3);
        dc.open(dcname);
        
        // FIXME: QCD histograms are empty and thus excluded.
        
        dc << "imax 1 number of channels" << std::endl;
        dc << "jmax " << jmax << " number of processes minus 1 ('*' = automatic)" << std::endl;
        dc << "kmax * number of nuisance parameters (sources of systematical uncertainties)" << std::endl;
        dc << "-----------------------------------" << std::endl;
        dc << "shapes * * " << wsname << " $CHANNEL:$PROCESS $CHANNEL:$PROCESS_$SYSTEMATIC" << std::endl;
        dc << "-----------------------------------" << std::endl;
        dc << "bin         " << channel_8TeV << std::endl;
        dc << "observation " << hdata_obs->Integral() << std::endl;
        dc << "#prediction " << hmc_exp->Integral() << std::endl;
        dc << "-----------------------------------" << std::endl;
        dc << "bin         "; for (Int_t j=0; j!=jmax+1; j++)  dc << channel_8TeV << " "; dc << std::endl;
        dc << "process     ZH         WH         WjLF       WjHF       ZjLF       ZjHF       TT         s_Top      VVLF       VVHF       QCD        " << std::endl;
        dc << "process     9          10         1          2          3          4          5          6          7          0          8          " << std::endl;
        dc << "rate        " << setw(10) << right << hZH->Integral() << " " << setw(10) << right << hWH->Integral() << " " << setw(10) << right << hWjLF->Integral() << " " << setw(10) << right << hWjHF->Integral() << " " << setw(10) << right << hZjLF->Integral() << " " << setw(10) << right << hZjHF->Integral() << " " << setw(10) << right << hTT->Integral() << " " << setw(10) << right << hs_Top->Integral() << " " << setw(10) << right << hVVLF->Integral() << " " << setw(10) << right << hVVHF->Integral() << " " << setw(10) << right << hQCD->Integral() << std::endl;
        dc.precision(2);
        dc << "-----------------------------------" << std::endl;
        dc << "" << std::endl;
        dc << "### Flat ########################## #####  ZH    WH    WjLF  WjHF  ZjLF  ZjHF  TT    s_Top VVLF  VVHF  QCD  " << std::endl;
        dc << "lumi_8TeV                           lnN    1.05  1.05  -     -     -     -     -     1.05  1.05  1.05  1.05 " << std::endl;
        dc << "pdf_qqbar                           lnN    1.01  1.01  -     -     -     -     -     -     1.01  1.01  -    " << std::endl;
        dc << "pdf_gg                              lnN    -     -     -     -     -     -     -     1.01  -     -     1.01 " << std::endl;
        dc << "QCDscale_VH                         lnN    1.04  1.04  -     -     -     -     -     -     -     -     -    " << std::endl;
        dc << "QCDscale_ttbar                      lnN    -     -     -     -     -     -     -     1.06  -     -     -    " << std::endl;
        dc << "QCDscale_VV                         lnN    -     -     -     -     -     -     -     -     1.04  1.04  -    " << std::endl;
        dc << "QCDscale_QCD                        lnN    -     -     -     -     -     -     -     -     -     -     1.30 " << std::endl;
        dc << "CMS_vhbb_boost_EWK                  lnN    1.05  1.10  -     -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_boost_QCD                  lnN    1.10  1.10  -     -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_ST                         lnN    -     -     -     -     -     -     -     1.25  -     -     -    " << std::endl;
        dc << "CMS_vhbb_VV                         lnN    -     -     -     -     -     -     -     -     1.30  1.30  -    " << std::endl;
        dc << "#CMS_vhbb_MET_nojets                lnN    1.03  1.03  -     -     -     -     -     1.03  1.03  1.03  -    " << std::endl;
        dc << "#CMS_vhbb_trigger_MET               lnN    1.03  1.03  -     -     -     -     -     1.03  1.03  1.03  -    " << std::endl;
        dc << "CMS_vhbb_WjLF_SF_" << channel_8TeV_1 << "   lnN    -     -     " << scalefactors_lnN[0] << "  -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_WjHF_SF_" << channel_8TeV_1 << "   lnN    -     -     -     " << scalefactors_lnN[1] << "  -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_ZjLF_SF_" << channel_8TeV_1 << "   lnN    -     -     -     -     " << scalefactors_lnN[2] << "  -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_ZjHF_SF_" << channel_8TeV_1 << "   lnN    -     -     -     -     -     " << scalefactors_lnN[3] << "  -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_TT_SF_" << channel_8TeV_1 << "     lnN    -     -     -     -     -     -     " << scalefactors_lnN[4] << "  -     -     -     -    " << std::endl;
        dc << "#CMS_vhbb_QCD_SF_" << channel_8TeV_1 << "   lnN    -     -     -     -     -     -     -     -     -     -    " << std::endl;
        dc << "### Shape ######################### #####  ZH    WH    WjLF  WjHF  ZjLF  ZjHF  TT    s_Top VVLF  VVHF  QCD  " << std::endl;
        dc << "UEPS                                shape  1.00  1.00  -     -     -     -     -     1.00  1.00  1.00  -    " << std::endl;
        dc << "CMS_vhbb_eff_b_SIG                  shape  1.00  1.00  -     -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_eff_b                      shape  -     -     1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  -    " << std::endl;
        dc << "CMS_vhbb_fake_b_8TeV                shape  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  -    " << std::endl;
        dc << "CMS_vhbb_res_j                      shape  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  -    " << std::endl;
        dc << "CMS_vhbb_scale_j                    shape  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  -    " << std::endl;
        dc << "CMS_vhbb_trigger_MET_Znunu_8TeV     shape  1.00  1.00  -     -     -     -     -     1.00  1.00  1.00  -    " << std::endl;
        dc << "CMS_vhbb_ZJModel_Znunu_8TeV         shape  -     -     -     -     1.00  1.00  -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_TTModel_Znunu_8TeV         shape  -     -     -     -     -     -     1.00  -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_statZH_" << channel_8TeV << "    shape  1.00  -     -     -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_statWH_" << channel_8TeV << "    shape  -     1.00  -     -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_statWjLF_" << channel_8TeV << "  shape  -     -     1.00  -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_statWjHF_" << channel_8TeV << "  shape  -     -     -     1.00  -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_statZjLF_" << channel_8TeV << "  shape  -     -     -     -     1.00  -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_statZjHF_" << channel_8TeV << "  shape  -     -     -     -     -     1.00  -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_statTT_" << channel_8TeV << "    shape  -     -     -     -     -     -     1.00  -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_stats_Top_" << channel_8TeV << " shape  -     -     -     -     -     -     -     1.00  -     -     -    " << std::endl;
        dc << "CMS_vhbb_statVVLF_" << channel_8TeV << "  shape  -     -     -     -     -     -     -     -     1.00  -     -    " << std::endl;
        dc << "CMS_vhbb_statVVHF_" << channel_8TeV << "  shape  -     -     -     -     -     -     -     -     -     1.00  -    " << std::endl;
        dc << "CMS_vhbb_statQCD_" << channel_8TeV << "   shape  -     -     -     -     -     -     -     -     -     -     -    " << std::endl;
        dc << "################################### #####  ZH    WH    WjLF  WjHF  ZjLF  ZjHF  TT    s_Top VVLF  VVHF  QCD  " << std::endl;

        dc.close();
        
        std::clog << "MakePlots(): The datacard is written." << std::endl;
    } // end if syst == NONE


    if (syst == "NONE") {
        std::clog << "MakePlots(): Setting up histograms..." << std::endl;
        
        // Setup canvas and pads
        TCanvas * c1 = new TCanvas("c1", "c1", 700, 700);
        //TCanvas * c1 = new TCanvas("c1", "c1", 600, 600);

        TPad * pad1 = new TPad("pad1", "top pad"   , 0.0, 0.3, 1.0, 1.0);
        pad1->SetBottomMargin(0.0);
        pad1->Draw();
        TPad * pad2 = new TPad("pad2", "bottom pad", 0.0, 0.0, 1.0, 0.3);
        pad2->SetTopMargin(0.0);
        pad2->SetBottomMargin(0.35);
        pad2->Draw();
        pad1->cd();
        pad1->SetLogy(plotLog);

        // Setup histogram styles and stack the histograms.
        setHisto(hVH, "VH");
        setHisto(hZH, "VH");
        setHisto(hWH, "VH");
        setHisto(hWjLF, "WjLF");
        setHisto(hWjHF, "WjHFb");
        setHisto(hZjLF, "ZjLF");
        setHisto(hZjHF, "ZjHFb");
        setHisto(hTT, "TT");
        setHisto(hs_Top, "s_Top");
        setHisto(hVVLF, "VV");
        setHisto(hVVHF, "VVHF");
        setHisto(hQCD, "QCD");
        setHisto(hZjLFHW, "ZjLF");
        setHisto(hZjHFHW, "ZjHFb");
        setHisto(hTTMG, "TT");
        
        TH1F * hdata_plot = (TH1F *) hdata_obs->Clone("hdata_plot");  // blinded plot
        TH1F * hmc_test = (TH1F *) hmc_exp->Clone("hmc_test");  // for chi2 and KS test
        hdata_plot->Sumw2();
        hmc_test->Sumw2();
        Int_t nbins_plot = hdata_plot->GetNbinsX();
        assert(nbins_plot == hmc_test->GetNbinsX());
        if (!plotData) {  // be blind to the most sensitive bins
            for (Int_t i = TMath::Max(nbins_plot*0.85, nbins_plot-5.); i < nbins_plot+2; i++) {
                hdata_plot->SetBinContent(i, 0.);
                hdata_plot->SetBinError(i, 0.);
                hmc_test->SetBinContent(i, 0.);
                hmc_test->SetBinError(i, 0.);
            }
        }
        setHisto(hdata_plot, "data_obs");

        std::clog << "MakePlots(): Setting up the stack..." << std::endl;
        THStack * hs = new THStack("hs", "");
        //hs->Add(hQCD);
        hs->Add(hVVLF);
        hs->Add(hVVHF);  // FIXME: should be plotted as signal?
        hs->Add(hs_Top);
        hs->Add(hTT);
        hs->Add(hWjLF);
        hs->Add(hWjHF);
        hs->Add(hZjLF);
        hs->Add(hZjHF);
        if (plotSig) {
            hs->Add(hZH);
            hs->Add(hWH);
        }
        
        double ymax = TMath::Max(hdata_plot->GetMaximum(), hs->GetMaximum());
        hs->SetMaximum(ymax + ymax / 2.0 + sqrt(ymax));
        if (plotLog)
            hs->SetMaximum(ymax * 200 + sqrt(ymax));
        hs->SetMinimum(0.01);
        
        // Setup auxiliary histograms
        std::clog << "MakePlots(): Setting up auxiliary histograms..." << std::endl;
        TH1F * staterr = (TH1F *) hmc_exp->Clone("staterr");
        staterr->Sumw2();
        staterr->SetFillColor(kRed);
        staterr->SetMarkerSize(0);
        staterr->SetFillStyle(3013);
        
        TH1F * ratio = (TH1F *) hdata_plot->Clone("ratio");
        ratio->Sumw2();
        ratio->SetMarkerSize(0.8);
        //ratio->SetMarkerSize(0.5);
        ratio->Divide(hdata_plot, hmc_exp, 1., 1., "");
        
        TH1F * ratiostaterr = (TH1F *) hmc_exp->Clone("ratiostaterr");
        ratiostaterr->Sumw2();
        ratiostaterr->SetStats(0);
        ratiostaterr->SetTitle("");
        ratiostaterr->GetXaxis()->SetTitle(axtitle);
        ratiostaterr->GetYaxis()->SetTitle("Data/MC");
        ratiostaterr->SetMinimum(0);
        ratiostaterr->SetMaximum(2.2);
        ratiostaterr->SetMarkerSize(0);
        ratiostaterr->SetFillColor(kRed);
        ratiostaterr->SetFillStyle(3013);
        //ratiostaterr->SetFillStyle(3001);
        ratiostaterr->GetXaxis()->CenterTitle();
        ratiostaterr->GetXaxis()->SetLabelSize(0.12);
        ratiostaterr->GetXaxis()->SetTitleSize(0.14);
        ratiostaterr->GetXaxis()->SetTitleOffset(1.10);
        ratiostaterr->GetYaxis()->CenterTitle();
        ratiostaterr->GetYaxis()->SetLabelSize(0.10);
        ratiostaterr->GetYaxis()->SetTitleSize(0.12);
        //ratiostaterr->GetYaxis()->SetTitleSize(0.10);
        ratiostaterr->GetYaxis()->SetTitleOffset(0.6);
        ratiostaterr->GetYaxis()->SetNdivisions(505);
        
        for (Int_t i = 0; i < hmc_exp->GetNbinsX()+2; i++) {
            ratiostaterr->SetBinContent(i, 1.0);
            if (hmc_exp->GetBinContent(i) > 1e-6) {  // use smaller tolerance?
                double binerror = hmc_exp->GetBinError(i) / hmc_exp->GetBinContent(i);
                ratiostaterr->SetBinError(i, binerror);
            } else {
                ratiostaterr->SetBinError(i, 999.);
            }
            
            //if (!(hdata_plot->GetBinContent(i) > 1e-6)) {
            //    ratiostaterr->SetBinError(i, 0.);
            //}
        }
        
        TH1F * ratiosysterr = (TH1F *) ratiostaterr->Clone("ratiosysterr");
        ratiosysterr->Sumw2();
        ratiosysterr->SetMarkerSize(0);
        //ratiosysterr->SetFillColor(kBlue);
        ratiosysterr->SetFillColor(kYellow-4);
        //ratiosysterr->SetFillStyle(3002);
        ratiosysterr->SetFillStyle(1001);
        
        for (Int_t i = 0; i < hmc_exp->GetNbinsX()+2; i++) {
            if (hmc_exp->GetBinContent(i) > 1e-6) {  // use smaller tolerance?
                double binerror2 = (pow(hmc_exp->GetBinError(i), 2) +
                                    pow(0.07 * hWjLF->GetBinContent(i), 2) +
                                    pow(0.15 * hWjHF->GetBinContent(i), 2) +
                                    pow(0.07 * hZjLF->GetBinContent(i), 2) +
                                    pow(0.15 * hZjHF->GetBinContent(i), 2) +
                                    pow(0.08 * hTT->GetBinContent(i), 2) +
                                    pow(0.25 * hs_Top->GetBinContent(i), 2) +
                                    pow(0.30 * hVVLF->GetBinContent(i), 2));  // FIXME: add VH as background?
                double binerror = sqrt(binerror2) / hmc_exp->GetBinContent(i);
                ratiosysterr->SetBinError(i, binerror);
            }
        }

        // Setup legends
        std::clog << "MakePlots(): Setting up legends..." << std::endl;
        TLegend * leg = new TLegend(0.72, 0.62, 0.92, 0.92);
        leg->SetFillColor(0);
        leg->SetLineColor(0);
        leg->SetShadowColor(0);
        leg->SetTextFont(62);
        leg->SetTextSize(0.03);
        leg->AddEntry(hdata_plot, "Data", "p");
        if (plotSig)
            leg->AddEntry(hVH, Form("VH(%i)", massH), "l");
        leg->AddEntry(hZjHF, "Z + b#bar{b}", "f");
        leg->AddEntry(hZjLF, "Z + udscg", "f");
        leg->AddEntry(hWjHF, "W + b#bar{b}", "f");
        leg->AddEntry(hWjLF, "W + udscg", "f");
        leg->AddEntry(hTT, "t#bar{t}", "f");
        leg->AddEntry(hs_Top, "single top", "f");
        leg->AddEntry(hVVHF, "VV(b#bar{b})", "f");
        leg->AddEntry(hVVLF, "VV(udscg)", "f");
        //leg->AddEntry(hQCD, "QCD", "f");
        leg->AddEntry(staterr, "MC uncert. (stat)", "fl");

        TLegend * ratioleg1 = new TLegend(0.54, 0.86, 0.72, 0.96);
        //TLegend * ratioleg1 = new TLegend(0.50, 0.86, 0.69, 0.96);
        ratioleg1->AddEntry(ratiostaterr, "MC uncert. (stat)", "f");
        ratioleg1->SetFillColor(0);
        ratioleg1->SetLineColor(0);
        ratioleg1->SetShadowColor(0);
        ratioleg1->SetTextFont(62);
        ratioleg1->SetTextSize(0.06);
        ratioleg1->SetBorderSize(1);
        
        TLegend * ratioleg2 = new TLegend(0.72, 0.86, 0.95, 0.96);
        //TLegend * ratioleg2 = new TLegend(0.69, 0.86, 0.9, 0.96);
        ratioleg2->AddEntry(ratiosysterr, "MC uncert. (stat+syst)", "f");
        ratioleg2->SetFillColor(0);
        ratioleg2->SetLineColor(0);
        ratioleg2->SetShadowColor(0);
        ratioleg2->SetTextFont(62);
        ratioleg2->SetTextSize(0.06);
        ratioleg2->SetBorderSize(1);

        // Draw MC signal and background
        std::clog << "MakePlots(): Drawing..." << std::endl;
        pad1->cd();
        if (plotLog) pad1->SetLogy();
        hs->Draw("hist");
        hs->GetXaxis()->SetLabelSize(0);
        double binwidth = (xup - xlow) / nbins_plot;
        TString ytitle = Form("Events / %.3f", binwidth);
        hs->GetYaxis()->SetTitle(ytitle);
        
        staterr->Draw("e2 same");
        if (plotSig) {
            hVH->SetLineWidth(3);
            hVH->SetFillColor(0);
            hVH->Draw("hist same");
        }
        
        // Draw data
        hdata_plot->Draw("e1 same");
        
        // Draw legends
        leg->Draw();
        TLatex * latex = new TLatex();
        latex->SetNDC();
        latex->SetTextAlign(12);
        latex->SetTextFont(62);
        latex->SetTextSize(0.052);
        latex->DrawLatex(0.19, 0.89, "CMS Preliminary");
        latex->SetTextSize(0.04);
        latex->DrawLatex(0.19, 0.84, "#sqrt{s} = 8 TeV, L = 12.1 fb^{-1}");
        latex->DrawLatex(0.19, 0.79, "Z(#nu#bar{#nu})H(b#bar{b})");
        
        // Under/overflows a la TMVA
        TString uoflow = Form("U/O-flow (Data,MC): (%.1f, %.1f) / (%.1f, %.1f)", hdata_plot->GetBinContent(0), hmc_exp->GetBinContent(0), hdata_plot->GetBinContent(nbins_plot+1), hmc_exp->GetBinContent(nbins_plot+1));
        TLatex * latex2 = new TLatex(0.99, 0.1, uoflow);
        latex2->SetNDC();
        latex2->SetTextSize(0.02);
        latex2->SetTextAngle(90);
        latex2->AppendPad();
        
        // Draw ratio
        pad2->cd();
        pad2->SetGridy(0);
        ratiostaterr->Draw("e2");
        ratiosysterr->Draw("e2 same");
        ratiostaterr->Draw("e2 same");
        ratio->Draw("e1 same");
        
        // Draw ratio legends
        ratioleg1->Draw();
        ratioleg2->Draw();

        // Kolmogorov-Smirnov test and Chi2 test
        TPaveText * pave = new TPaveText(0.18, 0.85, 0.35, 0.96, "brNDC");
        pave->SetTextAlign(12);
        pave->SetLineColor(0);
        pave->SetFillColor(0);
        pave->SetShadowColor(0);
        pave->SetBorderSize(1);
        double nchisq = hdata_plot->Chi2Test(hmc_test, "UWCHI2/NDF");  // MC uncert. (stat)
        double kolprob = hdata_plot->KolmogorovTest(hmc_test);  // MC uncert. (stat)
        TText * text = pave->AddText(Form("#chi_{#nu}^{2} = %.3f, K_{s} = %.3f", nchisq, kolprob));
        text->SetTextFont(62);
        text->SetTextSize(0.06);
        pave->Draw();

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
        TString plotname = Form("plotsJ10/%s_%s_%s", channel.Data(), region.Data(), var.Data());
        FormatFileName(plotname);
        gPad->Print(plotname+".png");
        gPad->Print(plotname+".pdf");
        //gPad->Print(plotname+".eps");
        
        delete hdata_plot;
        delete hmc_test;
        delete staterr;
        delete ratio;
        delete ratiostaterr;
        delete ratiosysterr;
        delete leg;
        delete ratioleg1;
        delete ratioleg2;
        delete latex;
        delete pave;
        delete hs;
        delete pad1;
        delete pad2;
        delete c1;
    }  // end if syst == NONE

    if (writeRoot) {
        std::clog << "MakePlots(): Writing out histograms..." << std::endl;
        TFile * out(0);
        TString rootname = g_rootname;
        rootname.ReplaceAll("$CHANNEL", channel);
        if (syst == "NONE") {
            out = TFile::Open(rootname, "RECREATE");
            out->mkdir(channel);
        } else {
            out = TFile::Open(rootname, "UPDATE");
        }
        out->cd(channel);
        
        // Write the histograms
        for (UInt_t ih = 0; ih < histos.size(); ih++) {
            if (syst == "NONE") {
                histos.at(ih)->Write(histos.at(ih)->GetName());
            } else {
                const TString p = histos.at(ih)->GetName();
                if (p != "data_obs") {
                    histos.at(ih)->Write(histos.at(ih)->GetName() + ("_" + syst));
                }
            }
        }

        // Add histograms for stat, ZjModel, TTModel
        if (syst == "NONE") {
            std::clog << "MakePlots(): Writing out stat, model syst histograms..." << std::endl;
            TString name;
            for (UInt_t ih = 0; ih < histos.size(); ih++) {
                const TString p = histos.at(ih)->GetName();
                if (p != "data_obs") {
                    name = Form("%s_CMS_vhbb_stat%s_%s_8TeVUp", p.Data(), p.Data(), channel.Data());
                    TH1F * h_statUp = (TH1F *) VaryStatErrors((TH1F *) histos.at(ih), 1.0, name);
                    name = Form("%s_CMS_vhbb_stat%s_%s_8TeVDown", p.Data(), p.Data(), channel.Data());
                    TH1F * h_statDown = (TH1F *) VaryStatErrors((TH1F *) histos.at(ih), -1.0, name);
                    h_statUp->Write(h_statUp->GetName());
                    h_statDown->Write(h_statDown->GetName());
                    
                    delete h_statUp;
                    delete h_statDown;
                }
            }
            name = "ZjLF_CMS_vhbb_ZJModel_Znunu_8TeVUp";
            TH1F * hZjLF_modelUp = VaryModelErrors((TH1F *) hZjLF, (TH1F *) hZjLFHW, 1.0, name);
            name = "ZjHF_CMS_vhbb_ZJModel_Znunu_8TeVUp";
            TH1F * hZjHF_modelUp = VaryModelErrors((TH1F *) hZjHF, (TH1F *) hZjHFHW, 1.0, name);
            name = "ZjLF_CMS_vhbb_ZJModel_Znunu_8TeVDown";
            TH1F * hZjLF_modelDown = VaryModelErrors((TH1F *) hZjLF, (TH1F *) hZjLFHW, -1.0, name);
            name = "ZjHF_CMS_vhbb_ZJModel_Znunu_8TeVDown";
            TH1F * hZjHF_modelDown = VaryModelErrors((TH1F *) hZjHF, (TH1F *) hZjHFHW, -1.0, name);
            hZjLF_modelUp->Write(hZjLF_modelUp->GetName());
            hZjHF_modelUp->Write(hZjHF_modelUp->GetName());
            hZjLF_modelDown->Write(hZjLF_modelDown->GetName());
            hZjHF_modelDown->Write(hZjHF_modelDown->GetName());
            
            delete hZjLF_modelUp;
            delete hZjHF_modelUp;
            delete hZjLF_modelDown;
            delete hZjHF_modelDown;
            
            name = "TT_CMS_vhbb_TTModel_Znunu_8TeVUp";
            TH1F * hTT_modelUp = VaryModelErrors((TH1F *) hTT, (TH1F *) hTTMG, 1.0, name);
            name = "TT_CMS_vhbb_TTModel_Znunu_8TeVDown";
            TH1F * hTT_modelDown = VaryModelErrors((TH1F *) hTT, (TH1F *) hTTMG, -1.0, name);
            
            hTT_modelUp->Write(hTT_modelUp->GetName());
            hTT_modelDown->Write(hTT_modelDown->GetName());
            
            delete hTT_modelUp;
            delete hTT_modelDown;
        }
        out->Close();
        delete out;
    }

    for (UInt_t ih = 0; ih < histos_0.size(); ih++)
        delete histos_0.at(ih);
    for (UInt_t ih = 0; ih < histos.size(); ih++)
        delete histos.at(ih);

    std::clog << "MakePlots(): DONE!" << std::endl;

    return;
}


EventsJ10 * Read(const TString region, const TCut cutallmc, const TCut cutalldata)
{
    std::clog << "Read(): Using cutallmc = " << cutallmc << ", cutalldata = " << cutalldata << std::endl;
    std::clog << "Read():       indir = " << indir << std::endl;
    
    TString treename = Form("tree_%s_ctrl", channel.Data());
    if (region == "VH")
        treename = Form("tree_%s_test", channel.Data());
    
    TChain ZH(treename);
    ZH.Add(indir + prefix + Form("ZH%i", massH) + suffix);
    
    TChain WH(treename);
    WH.Add(indir + prefix + Form("WH%i", massH) + suffix);
    
    TChain Wj(treename);
    Wj.Add(indir + prefix + "Wj" + suffix);
    
    TChain Zj(treename);
    Zj.Add(indir + prefix + "Zj" + suffix);
    
    TChain TT(treename);
    TT.Add(indir + prefix + "TT" + suffix);
    
    TChain s_Top(treename);
    s_Top.Add(indir + prefix + "s_Top" + suffix);
    
    TChain VV(treename);
    VV.Add(indir + prefix + "VV" + suffix);
    
    TChain QCD(treename);
    QCD.Add(indir + prefix + "QCD" + suffix);
    
    TChain ZjHW(treename);
    ZjHW.Add(indir + prefix + "ZjHW" + suffix);
    
    TChain TTMG(treename);
    TTMG.Add(indir + prefix + "TTMG" + suffix);
    
    TChain ZH_SM(treename);
    ZH_SM.Add(indir + prefix + "ZH125" + suffix);
    
    TChain WH_SM(treename);
    WH_SM.Add(indir + prefix + "WH125" + suffix);
    
    TChain data_obs(treename);
    data_obs.Add(indir + prefix + "data_obs" + suffix);

    // Start copying trees
    const TCut cutHF = "eventFlav==5";
    const TCut cutLF = "eventFlav!=5";
    
    EventsJ10 * ev = new EventsJ10();
    ev->ZH = (TTree *) ZH.CopyTree(cutallmc);
    std::clog << "... DONE: ZH copy tree." << std::endl;
    ev->WH = (TTree *) WH.CopyTree(cutallmc);
    std::clog << "... DONE: WH copy tree." << std::endl;
    ev->WjLF = (TTree *) Wj.CopyTree(cutallmc + cutLF);
    std::clog << "... DONE: WjLF copy tree." << std::endl;
    ev->WjHF = (TTree *) Wj.CopyTree(cutallmc + cutHF);
    std::clog << "... DONE: WjHF copy tree." << std::endl;
    ev->ZjLF = (TTree *) Zj.CopyTree(cutallmc + cutLF);
    std::clog << "... DONE: ZjLF copy tree." << std::endl;
    ev->ZjHF = (TTree *) Zj.CopyTree(cutallmc + cutHF);
    std::clog << "... DONE: ZjHF copy tree." << std::endl;
    ev->TT = (TTree *) TT.CopyTree(cutallmc);
    std::clog << "... DONE: TT copy tree." << std::endl;
    ev->s_Top = (TTree *) s_Top.CopyTree(cutallmc);
    std::clog << "... DONE: s_Top copy tree." << std::endl;
    ev->VVLF = (TTree *) VV.CopyTree(cutallmc + cutLF);
    std::clog << "... DONE: VVLF copy tree." << std::endl;
    ev->VVHF = (TTree *) VV.CopyTree(cutallmc + cutHF);
    std::clog << "... DONE: VVHF copy tree." << std::endl;
    ev->QCD = (TTree *) QCD.CopyTree(cutallmc);
    std::clog << "... DONE: QCD copy tree." << std::endl;
    ev->ZjLFHW = (TTree *) ZjHW.CopyTree(cutallmc + cutLF);
    std::clog << "... DONE: ZjLFHW copy tree." << std::endl;
    ev->ZjHFHW = (TTree *) ZjHW.CopyTree(cutallmc + cutHF);
    std::clog << "... DONE: ZjHFHW copy tree." << std::endl;
    ev->TTMG = (TTree *) TTMG.CopyTree(cutallmc);
    std::clog << "... DONE: TTMG copy tree." << std::endl;
    ev->ZH_SM = (TTree *) ZH_SM.CopyTree(cutallmc);
    std::clog << "... DONE: ZH_SM copy tree." << std::endl;
    ev->WH_SM = (TTree *) WH_SM.CopyTree(cutallmc);
    std::clog << "... DONE: WH_SM copy tree." << std::endl;
    
    ev->data_obs = (TTree *) data_obs.CopyTree(cutalldata);
    std::cout << "... DONE: data_obs copy tree." << std::endl;

    return ev;
}

////////////////////////////////////////////////////////////////////////////////
/// Main                                                                     ///
////////////////////////////////////////////////////////////////////////////////

void BDTShapeJ10(Int_t nbins=500, Int_t newnbins=16, Double_t rebinerrorf=0.35, bool pqrs=false, bool isControlRegion=false)
{
    gROOT->LoadMacro("tdrstyle.C");
    gROOT->ProcessLine("setTDRStyle()");
    
    TH1::SetDefaultSumw2(1);
    gROOT->SetBatch(1);
    
    assert(newnbins <= nbins);  // nbins is used in plots before rebinning
    
    if (!TString(gROOT->GetVersion()).Contains("5.32")) {
        std::cout << "INCORRECT ROOT VERSION! Please use 5.32:" << std::endl;
        std::cout << "source /uscmst1/prod/sw/cms/slc5_amd64_gcc462/lcg/root/5.32.00-cms14/bin/thisroot.csh" << std::endl;
        std::cout << "Return without doing anything." << std::endl;
        return;
    }
    
    if (gSystem->AccessPathName("plotsJ10"))
        gSystem->mkdir("plotsJ10");
    
    
    ///-- ZnunuHighPt ----------------------------------------------------------
    channel = "ZnunuHighPt";
    std::clog << "\n*** " << channel << " ***\n" << std::endl;

    int begin = (isControlRegion) ?   1 : 0;
    int end   = (isControlRegion) ? 1+5 : 1;
    for (Int_t ireg=begin; ireg<end; ireg++) {
        EventsJ10 * ev = Read(g_regions[ireg], g_cutallmc, g_cutalldata);
        ev->check();
    
        /// Select scale factors
        const Double_t * scalefactors = g_scalefactors_ZnunuHighPt;
        const Double_t * scalefactors_lnN = g_scalefactors_lnN_ZnunuHighPt;
        
        /// Initialize rebinning "bookmarks" to negative values
        /// These global variables that need to be reset when we switch channel.
        firstN = -1, lastN = -1;  // for single BDT
        firstN1 = -1, firstN2 = -1, firstN3 = -1, firstN4 = -1; // for PQRS BDT
        lastN1  = -1, lastN2  = -1, lastN3  = -1, lastN4  = -1; // for PQRS BDT

        /// Make BDT plots
        /// Loop over systematics
        for (Int_t isyst = 0; isyst < nsyst; isyst++) {
            ev->set_scale_factors(&scalefactors[isyst*5]);
            std::clog << "isyst " << isyst << ": scale factors: " << ev->sf_WjLF << " (WjLF), " << ev->sf_WjHF << " (WjHF), " << ev->sf_ZjLF << " (ZjLF), " << ev->sf_ZjHF << " (ZjHF), " << ev->sf_TT << " (TT)." << std::endl;
            
            /// Replace "ISYST" with the index of systematic
            TString var = (!pqrs) ? g_var : g_varslice;
            TString cutmc=g_cutmc, cutdata=g_cutdata, ssyst=Form("%i",isyst);
            var.ReplaceAll("ISYST",ssyst);
            cutmc.ReplaceAll("ISYST",ssyst);
            cutdata.ReplaceAll("ISYST",ssyst);
            
            MakePlots(ev, var, cutmc.Data(), cutdata.Data(), g_systematics[isyst], g_regions[ireg], "BDT", nbins, g_xlow, g_xup, newnbins, rebinerrorf, scalefactors_lnN, "printStat:printCard:writeRoot:!plotData:plotLog:plotSig");
            std::clog << "isyst " << isyst << ": DONE!\n" << std::endl;
        }
        
        /// Make plots of 13 input variables
        ev->set_scale_factors(&scalefactors[0*5]);
        TCut cutmc_ctrl = Form("weightsHCP[0] * (selectFlags[%i][0])", ireg);
        TCut cutdata_ctrl = Form("(selectFlags[%i][0])", ireg);
        TCut cutmass_ctrl = "";
        if (ireg==0) {  // signal region
            cutmc_ctrl = Form("2.0 * weightsHCP[0] * (selectFlags[%i][0])", ireg);
            cutmass_ctrl = "(HmassReg<100 || 140<HmassReg)";
        }
        
        if (ireg!=0) {  // control region
            MakePlots(ev, Form("BDTregular_%i[0]", massH), cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "BDT", nbins, g_xlow, g_xup, newnbins, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:plotLog:plotSig");  // apply mass veto
        }
        MakePlots(ev, "HmassReg", cutmc_ctrl*cutmass_ctrl, cutdata_ctrl*cutmass_ctrl, g_systematics[0], g_regions[ireg], "H mass [GeV]", 25, 0, 250, 25, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");  // apply mass veto
        MakePlots(ev, "HptReg", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H p_{T} [GeV]", 20, 50, 450, 20, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
        MakePlots(ev, "hJet_ptReg[0]", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 p_{T} [GeV]", 20, 50, 350, 20, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
        MakePlots(ev, "hJet_ptReg[1]", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j2 p_{T} [GeV]", 10, 20, 170, 10, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
        MakePlots(ev, "deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "#Delta R(j1,j2)", 16, 0, 3.2, 16, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
        MakePlots(ev, "abs(hJet_eta[0]-hJet_eta[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "|#Delta #eta|(j1,j2)", 16, 0, 3.2, 16, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
        MakePlots(ev, "min(abs(deltaPullAngle),pi)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "|#Delta #theta_{pull}|(j1,j2)", 16, 0, 3.2, 16, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
        MakePlots(ev, "METtype1corr.et", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "Type-1 corr. pfMET [GeV]", 20, 50, 450, 20, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
        MakePlots(ev, "HMETdPhi", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "|#Delta #varphi|(H, pfMET)", 16, 0, 3.2, 16, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
        MakePlots(ev, "mindPhiMETJet_dPhi", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "min |#Delta #varphi|(pfMET, cj30)", 16, 0, 3.2, 16, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
        MakePlots(ev, "max(hJet_csv_nominal[0],hJet_csv_nominal[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "CSV_{max}(j1,j2)", 20, 0, 1, 20, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
        MakePlots(ev, "min(hJet_csv_nominal[0],hJet_csv_nominal[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "CSV_{min}(j1,j2)", 20, 0, 1, 20, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
        MakePlots(ev, "naJets_Znn", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "# of add. jets", 8, 0, 8, 8, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
        MakePlots(ev, "MaxIf$(max(aJet_csv_nominal,0), aJet_pt>20 && abs(aJet_eta)<2.5)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "CSV_{max}(add. cj20)", 20, 0, 1, 20, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
        MakePlots(ev, "min(mindRHJet_dR,5.653)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "min #DeltaR(H, add. j30)", 20, 0, 6, 20, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
        
        // Angle (need the definition of the functions)
        //MakePlots(ev, "abs(evalCosThetaHbb(hJet_ptReg[0], hJet_pt[0], hJet_eta[0], hJet_phi[0], hJet_e[0], hJet_ptReg[1], hJet_pt[1], hJet_eta[1], hJet_phi[1], hJet_e[1]))", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "cos #theta* (H,j)", 20, 0, 1, 20, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
        //MakePlots(ev, "evalHMETMt(HptReg, H.phi, METtype1corr.et, METtype1corr.phi)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "mass_{T}(H,pfMET)", 25, 250, 1000, 20, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
        //MakePlots(ev, "evalHMETMassiveMt(HptReg, H.phi, METtype1corr.et, METtype1corr.phi)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "massive mass_{T}(H,pfMET)", 25, 250, 1000, 20, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
        
        // FatH.filteredmass
        cutmass_ctrl = "(FatH.filteredmass<100 || 140<FatH.filteredmass)";
        MakePlots(ev, "FatH.filteredmass * (nfathFilterJets>0 && FatH.FatHiggsFlag==1)", cutmc_ctrl*cutmass_ctrl, cutdata_ctrl*cutmass_ctrl, g_systematics[0], g_regions[ireg], "H_{fj} mass [GeV]", 25, 0, 250, 25, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
        MakePlots(ev, "FatH.filteredpt * (nfathFilterJets>0 && FatH.FatHiggsFlag==1)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H_{fj} p_{T} [GeV]", 20, 50, 450, 20, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
        
        // FatHmassReg
        cutmass_ctrl = "(FatHmassReg<100 || 140<FatHmassReg)";
        MakePlots(ev, "FatHmassReg * (nfathFilterJets>0 && FatH.FatHiggsFlag==1)", cutmc_ctrl*cutmass_ctrl, cutdata_ctrl*cutmass_ctrl, g_systematics[0], g_regions[ireg], "H_{fj} mass [GeV]", 25, 0, 250, 25, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
        MakePlots(ev, "FatHptReg * (nfathFilterJets>0 && FatH.FatHiggsFlag==1)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H_{fj} p_{T} [GeV]", 20, 50, 450, 20, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
        MakePlots(ev, "fathFilterJets_ptReg[0] * (nfathFilterJets>0 && FatH.FatHiggsFlag==1)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H fj1 p_{T} [GeV]", 20, 50, 350, 20, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
        MakePlots(ev, "fathFilterJets_ptReg[1] * (nfathFilterJets>0 && FatH.FatHiggsFlag==1)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H fj2 p_{T} [GeV]", 10, 20, 170, 10, rebinerrorf, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");

        delete ev;
    } // end loop over regions
    
    ///-- ZnunuLowPt -----------------------------------------------------------
    channel = "ZnunuLowPt";
    std::clog << "\n*** " << channel << " ***\n" << std::endl;
    
    EventsJ10 * ev = Read(g_regions[0], g_cutallmc, g_cutalldata);
    ev->check();
    
    /// Select scale factors
    const Double_t * scalefactors = g_scalefactors_ZnunuLowPt;
    const Double_t * scalefactors_lnN = g_scalefactors_lnN_ZnunuLowPt;
    
    /// Initialize rebinning "bookmarks" to negative values
    /// These global variables that need to be reset when we switch channel.
    firstN = -1, lastN = -1;  // for single BDT
    firstN1 = -1, firstN2 = -1, firstN3 = -1, firstN4 = -1; // for PQRS BDT
    lastN1  = -1, lastN2  = -1, lastN3  = -1, lastN4  = -1; // for PQRS BDT
    
    /// Loop over systematics
    for (Int_t isyst = 0; isyst < nsyst; isyst++) {
        newnbins = 9;  // FIXME: currently it's hardcoded
        rebinerrorf = 0.25;
        
        ev->set_scale_factors(&scalefactors[isyst*5]);
        std::clog << "isyst " << isyst << ": scale factors: " << ev->sf_WjLF << " (WjLF), " << ev->sf_WjHF << " (WjHF), " << ev->sf_ZjLF << " (ZjLF), " << ev->sf_ZjHF << " (ZjHF), " << ev->sf_TT << " (TT)." << std::endl;
        
        /// Replace "ISYST" with the index of systematic
        TString var = (!pqrs) ? g_var : g_varslice;
        TString cutmc=g_cutmc, cutdata=g_cutdata, ssyst=Form("%i",isyst);
        var.ReplaceAll("ISYST",ssyst);
        cutmc.ReplaceAll("ISYST",ssyst);
        cutdata.ReplaceAll("ISYST",ssyst);
        
        MakePlots(ev, var, cutmc.Data(), cutdata.Data(), g_systematics[isyst], g_regions[0], "BDT", nbins, g_xlow, g_xup, newnbins, rebinerrorf, scalefactors_lnN, "printStat:printCard:writeRoot:!plotData:plotLog:plotSig");
        std::clog << "isyst " << isyst << ": DONE!\n" << std::endl;
    }
    delete ev;
    
    ///-- ZnunuLowCSV ----------------------------------------------------------
    channel = "ZnunuLowCSV";
    std::clog << "\n*** " << channel << " ***\n" << std::endl;
    
    ev = Read(g_regions[0], g_cutallmc, g_cutalldata);
    ev->check();
    
    /// Select scale factors
    scalefactors = g_scalefactors_ZnunuLowCSV;
    scalefactors_lnN = g_scalefactors_lnN_ZnunuLowCSV;
    
    /// Initialize rebinning "bookmarks" to negative values
    /// These global variables that need to be reset when we switch channel.
    firstN = -1, lastN = -1;  // for single BDT
    firstN1 = -1, firstN2 = -1, firstN3 = -1, firstN4 = -1; // for PQRS BDT
    lastN1  = -1, lastN2  = -1, lastN3  = -1, lastN4  = -1; // for PQRS BDT
    
    /// Loop over systematics
    for (Int_t isyst = 0; isyst < nsyst; isyst++) {
        newnbins = 8;
        rebinerrorf = 0.35;
        
        ev->set_scale_factors(&scalefactors[isyst*5]);
        std::clog << "isyst " << isyst << ": scale factors: " << ev->sf_WjLF << " (WjLF), " << ev->sf_WjHF << " (WjHF), " << ev->sf_ZjLF << " (ZjLF), " << ev->sf_ZjHF << " (ZjHF), " << ev->sf_TT << " (TT)." << std::endl;
        
        /// Replace "ISYST" with the index of systematic
        TString var = (!pqrs) ? g_var : g_varslice;
        TString cutmc=g_cutmc, cutdata=g_cutdata, ssyst=Form("%i",isyst);
        var.ReplaceAll("ISYST",ssyst);
        cutmc.ReplaceAll("ISYST",ssyst);
        cutdata.ReplaceAll("ISYST",ssyst);
        
        MakePlots(ev, var, cutmc.Data(), cutdata.Data(), g_systematics[isyst], g_regions[0],"BDT", nbins, g_xlow, g_xup, newnbins, rebinerrorf, scalefactors_lnN, "printStat:printCard:writeRoot:!plotData:plotLog:plotSig");
        std::clog << "isyst " << isyst << ": DONE!\n" << std::endl;
    }
    delete ev;
}
