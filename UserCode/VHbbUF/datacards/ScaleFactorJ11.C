////////////////////////////////////////////////////////////////////////////////
/// This macro runs on Step 4 with jmax = 9 in the following order
///   ZH         WH         Wj0b       Wj1b       Wj2b       Zj0b       Zj1b       Zj2b       TT         s_Top      VV         QCD
///   -1         0          1          2          3          4          5          6          7          8          9          10
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
///   CMS_vhbb_WJModel_Znunu_8TeVUp
///   CMS_vhbb_WJModel_Znunu_8TeVDown
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

#include "HelperBDTShape.h"
#include "HelperFunctions.h"  // FIXME

////////////////////////////////////////////////////////////////////////////////
/// Configuration                                                            ///
////////////////////////////////////////////////////////////////////////////////

const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130314/stitch/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130302/stitch/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130221/stitch/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130214/stitch/";
//const TString indir     = "root://eoscms//eos/cms/store/caf/user/degrutto/2013/Step4_20130207/";
const TString prefix    = "Step4_";
const TString suffix    = ".root";

const TCut    g_cutallmc    = "";          // used in CopyTree()
const TCut    g_cutalldata  = g_cutallmc;  // used in CopyTree()
const double  g_xlow        = -1.;
const double  g_xup         = 1.;

const int jmax              = 11;  // as in datacard
const int nsyst             = 13;
const int massH             = 125;
const TString g_dcname      = "vhbb_Znn_SF_J11_$CHANNEL_8TeV.txt";
const TString g_wsname      = "vhbb_Znn_SF_J11_$CHANNEL_8TeV.root";
const TString g_rootname    = "vhbb_Znn_SF_J11_$CHANNEL_TH1.root";
const TString g_plotdir     = "plotsJ11/";

/// The channel
TString channel = "ZnunuHighPt";
TString channel_CR = channel;  // this is changed in MakePlots()

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

/// Scale factors of Wj0b, Wj1b, Wj2b, Zj0b, Zj1b, Zj2b, TT (as in datacard)
const double g_scalefactors_prefit[7] = {
    1.00, 2.00, 1.00, 1.00, 2.00, 1.00, 1.00
};
const double g_scalefactors_lnN_prefit[7] = {
    1.10, 1.50, 1.40, 1.10, 1.25, 1.20, 1.08
};
////////////////////////////////////////////////////////////////////////////////
/// END Configuration                                                        ///
////////////////////////////////////////////////////////////////////////////////

///_____________________________________________________________________________
/// Classes

class EventsJ11 {  // jmax = 11 in datacard, with 7 scale factors.
public:
    EventsJ11()
      : ZH(0),
        WH(0),
        Wj0b(0),
        Wj1b(0),
        Wj2b(0),
        Zj0b(0),
        Zj1b(0),
        Zj2b(0),
        TT(0),
        s_Top(0),
        VVLF(0),
        VVHF(0),
        QCD(0),
        Wj0b_syst(0),
        Wj1b_syst(0),
        Wj2b_syst(0),
        Zj0b_syst(0),
        Zj1b_syst(0),
        Zj2b_syst(0),
        TT_syst(0),
        ZH_SM(0),
        WH_SM(0),
        data_obs(0),
        sf_Wj0b(1.), sf_Wj1b(1.), sf_Wj2b(1.), sf_Zj0b(1.), sf_Zj1b(1.), sf_Zj2b(1.), sf_TT(1.) {}

    ~EventsJ11()
    {
        delete ZH;
        delete WH;
        delete Wj0b;
        delete Wj1b;
        delete Wj2b;
        delete Zj0b;
        delete Zj1b;
        delete Zj2b;
        delete TT;
        delete s_Top;
        delete VVLF;
        delete VVHF;
        delete QCD;
        delete Wj0b_syst;
        delete Wj1b_syst;
        delete Wj2b_syst;
        delete Zj0b_syst;
        delete Zj1b_syst;
        delete Zj2b_syst;
        delete TT_syst;
        delete ZH_SM;
        delete WH_SM;
        delete data_obs;
    }

    void check() const {
        assert(ZH        != 0 && ZH       ->GetEntriesFast() > 0);
        assert(WH        != 0 && WH       ->GetEntriesFast() > 0);
        assert(Wj0b      != 0 && Wj0b     ->GetEntriesFast() > 0);
        assert(Wj1b      != 0 && Wj1b     ->GetEntriesFast() > 0);
        assert(Wj2b      != 0 && Wj2b     ->GetEntriesFast() > 0);
        assert(Zj0b      != 0 && Zj0b     ->GetEntriesFast() > 0);
        assert(Zj1b      != 0 && Zj1b     ->GetEntriesFast() > 0);
        assert(Zj2b      != 0 && Zj2b     ->GetEntriesFast() > 0);
        assert(TT        != 0 && TT       ->GetEntriesFast() > 0);
        assert(s_Top     != 0 && s_Top    ->GetEntriesFast() > 0);
        assert(VVLF      != 0 && VVLF     ->GetEntriesFast() > 0);
        assert(VVHF      != 0 && VVHF     ->GetEntriesFast() > 0);
        //assert(QCD       != 0 && QCD      ->GetEntriesFast() > 0);
        assert(Wj0b_syst != 0 && Wj0b_syst->GetEntriesFast() > 0);
        assert(Wj1b_syst != 0 && Wj1b_syst->GetEntriesFast() > 0);
        assert(Wj2b_syst != 0 && Wj2b_syst->GetEntriesFast() > 0);
        assert(Zj0b_syst != 0 && Zj0b_syst->GetEntriesFast() > 0);
        assert(Zj1b_syst != 0 && Zj1b_syst->GetEntriesFast() > 0);
        assert(Zj2b_syst != 0 && Zj2b_syst->GetEntriesFast() > 0);
        assert(TT_syst   != 0 && TT_syst  ->GetEntriesFast() > 0);
        //assert(ZH_SM     != 0 && ZH_SM    ->GetEntriesFast() > 0);
        //assert(WH_SM     != 0 && WH_SM    ->GetEntriesFast() > 0);
        assert(data_obs  != 0 && data_obs ->GetEntriesFast() > 0);
        return;
    }
    
    void set_scale_factors(const double sf[7]) {
        sf_Wj0b = sf[0];
        sf_Wj1b = sf[1];
        sf_Wj2b = sf[2];
        sf_Zj0b = sf[3];
        sf_Zj1b = sf[4];
        sf_Zj2b = sf[5];
        sf_TT   = sf[6];
        return;
    }

public:
    TTree * ZH;
    TTree * WH;
    TTree * Wj0b;
    TTree * Wj1b;
    TTree * Wj2b;
    TTree * Zj0b;
    TTree * Zj1b;
    TTree * Zj2b;
    TTree * TT;
    TTree * s_Top;
    TTree * VVLF;
    TTree * VVHF;
    TTree * QCD;
    TTree * Wj0b_syst;  // for WJModel
    TTree * Wj1b_syst;  // for WJModel
    TTree * Wj2b_syst;  // for WJModel
    TTree * Zj0b_syst;  // for ZJModel
    TTree * Zj1b_syst;  // for ZJModel
    TTree * Zj2b_syst;  // for ZJModel
    TTree * TT_syst;    // for TTModel
    TTree * ZH_SM;      // for signal injection
    TTree * WH_SM;      // for signal injection
    TTree * data_obs;
    double sf_Wj0b;
    double sf_Wj1b;
    double sf_Wj2b;
    double sf_Zj0b;
    double sf_Zj1b;
    double sf_Zj2b;
    double sf_TT;
};  // EventsJ11

///_____________________________________________________________________________
/// Functions

///_____________________________________________________________________________
/// The core function that makes plots, prints stats, and writes the datacard.
/// Note: QCD histograms are empty.

using namespace std;

/// Declare rebinner
Rebinner* rebinner = 0;

void MakePlots(const EventsJ11 * ev, const TString var, const TString var2, 
               const TCut cutmc, const TCut cutdata, const TString syst, const TString region, 
               const char* xtitle, int nbins, double xlow, double xup, 
               long long newnbins, double rebinerrorf, 
               const double scalefactors_lnN[7], 
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
    TH1F * hZH_0        = new TH1F("ZH_0"       , "", nbins, xlow, xup);
    TH1F * hWH_0        = new TH1F("WH_0"       , "", nbins, xlow, xup);
    TH1F * hWj0b_0      = new TH1F("Wj0b_0"     , "", nbins, xlow, xup);
    TH1F * hWj1b_0      = new TH1F("Wj1b_0"     , "", nbins, xlow, xup);
    TH1F * hWj2b_0      = new TH1F("Wj2b_0"     , "", nbins, xlow, xup);
    TH1F * hZj0b_0      = new TH1F("Zj0b_0"     , "", nbins, xlow, xup);
    TH1F * hZj1b_0      = new TH1F("Zj1b_0"     , "", nbins, xlow, xup);
    TH1F * hZj2b_0      = new TH1F("Zj2b_0"     , "", nbins, xlow, xup);
    TH1F * hTT_0        = new TH1F("TT_0"       , "", nbins, xlow, xup);
    TH1F * hs_Top_0     = new TH1F("s_Top_0"    , "", nbins, xlow, xup);
    TH1F * hVVLF_0      = new TH1F("VVLF_0"     , "", nbins, xlow, xup);
    TH1F * hVVHF_0      = new TH1F("VVHF_0"     , "", nbins, xlow, xup);
    TH1F * hQCD_0       = new TH1F("QCD_0"      , "", nbins, xlow, xup);
    TH1F * hWj0b_syst_0 = new TH1F("Wj0b_syst_0", "", nbins, xlow, xup);
    TH1F * hWj1b_syst_0 = new TH1F("Wj1b_syst_0", "", nbins, xlow, xup);
    TH1F * hWj2b_syst_0 = new TH1F("Wj2b_syst_0", "", nbins, xlow, xup);
    TH1F * hZj0b_syst_0 = new TH1F("Zj0b_syst_0", "", nbins, xlow, xup);
    TH1F * hZj1b_syst_0 = new TH1F("Zj1b_syst_0", "", nbins, xlow, xup);
    TH1F * hZj2b_syst_0 = new TH1F("Zj2b_syst_0", "", nbins, xlow, xup);
    TH1F * hTT_syst_0   = new TH1F("TT_syst_0"  , "", nbins, xlow, xup);
    //TH1F * hZH_SM_0     = new TH1F("ZH_SM_0"    , "", nbins, xlow, xup);
    //TH1F * hWH_SM_0     = new TH1F("WH_SM_0"    , "", nbins, xlow, xup);
    TH1F * hVH_0        = new TH1F("VH_0"       , "", nbins, xlow, xup);
    TH1F * hVV_0        = new TH1F("VV_0"       , "", nbins, xlow, xup);
    TH1F * hmc_exp_0    = new TH1F("mc_exp_0"   , "", nbins, xlow, xup);
    TH1F * hdata_obs_0  = new TH1F("data_obs_0" , "", nbins, xlow, xup);
    
    std::vector<TH1 *> histos_0;
    histos_0.push_back(hZH_0);
    histos_0.push_back(hWH_0);
    histos_0.push_back(hWj0b_0);
    histos_0.push_back(hWj1b_0);
    histos_0.push_back(hWj2b_0);
    histos_0.push_back(hZj0b_0);
    histos_0.push_back(hZj1b_0);
    histos_0.push_back(hZj2b_0);
    histos_0.push_back(hTT_0);
    histos_0.push_back(hs_Top_0);
    histos_0.push_back(hVVLF_0);
    histos_0.push_back(hVVHF_0);
    histos_0.push_back(hQCD_0);
    histos_0.push_back(hWj0b_syst_0);
    histos_0.push_back(hWj1b_syst_0);
    histos_0.push_back(hWj2b_syst_0);
    histos_0.push_back(hZj0b_syst_0);
    histos_0.push_back(hZj1b_syst_0);
    histos_0.push_back(hZj2b_syst_0);
    histos_0.push_back(hTT_syst_0);
    //histos_0.push_back(hZH_SM_0);
    //histos_0.push_back(hWH_SM_0);
    histos_0.push_back(hVH_0);
    histos_0.push_back(hVV_0);
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
    cutsf = Form("%f", ev->sf_Wj0b);
    ev->Wj0b->Project("Wj0b_0", var, cutmc * cutsf);
    std::clog << "... DONE: project Wj0b_0, scaled by " << cutsf << "." << std::endl;
    cutsf = Form("%f", ev->sf_Wj1b);
    ev->Wj1b->Project("Wj1b_0", var, cutmc * cutsf);
    std::clog << "... DONE: project Wj1b_0, scaled by " << cutsf << "." << std::endl;
    cutsf = Form("%f", ev->sf_Wj2b);
    ev->Wj2b->Project("Wj2b_0", var, cutmc * cutsf);
    std::clog << "... DONE: project Wj2b_0, scaled by " << cutsf << "." << std::endl;
    cutsf = Form("%f", ev->sf_Zj0b);
    ev->Zj0b->Project("Zj0b_0", var, cutmc * cutsf);
    std::clog << "... DONE: project Zj0b_0, scaled by " << cutsf << "." << std::endl;
    cutsf = Form("%f", ev->sf_Zj1b);
    ev->Zj1b->Project("Zj1b_0", var, cutmc * cutsf);
    std::clog << "... DONE: project Zj1b_0, scaled by " << cutsf << "." << std::endl;
    cutsf = Form("%f", ev->sf_Zj2b);
    ev->Zj2b->Project("Zj2b_0", var, cutmc * cutsf);
    std::clog << "... DONE: project Zj2b_0, scaled by " << cutsf << "." << std::endl;
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

    ev->Wj0b_syst->Project("Wj0b_syst_0", var, cutmc);
    std::clog << "... DONE: project Wj0b_syst_0." << std::endl;
    ev->Wj1b_syst->Project("Wj1b_syst_0", var, cutmc);
    std::clog << "... DONE: project Wj1b_syst_0." << std::endl;
    ev->Wj2b_syst->Project("Wj2b_syst_0", var, cutmc);
    std::clog << "... DONE: project Wj2b_syst_0." << std::endl;
    ev->Zj0b_syst->Project("Zj0b_syst_0", var, cutmc);
    std::clog << "... DONE: project Zj0b_syst_0." << std::endl;
    ev->Zj1b_syst->Project("Zj1b_syst_0", var, cutmc);
    std::clog << "... DONE: project Zj1b_syst_0." << std::endl;
    ev->Zj2b_syst->Project("Zj2b_syst_0", var, cutmc);
    std::clog << "... DONE: project Zj2b_syst_0." << std::endl;
    ev->TT_syst->Project("TT_syst_0", var, cutmc);
    std::clog << "... DONE: project TT_syst_0." << std::endl;
    
    //ev->ZH_SM->Project("ZH_SM_0", var, cutmc);
    //std::clog << "... DONE: project ZH_SM_0." << std::endl;
    //ev->WH_SM->Project("WH_SM_0", var, cutmc);
    //std::clog << "... DONE: project WH_SM_0." << std::endl;
    
    ev->data_obs->Project("data_obs_0", var, cutdata);
    std::clog << "... DONE: project data_obs_0." << std::endl;
    
    // Project var2
    ev->ZH->Project("+ZH_0", var2, cutmc);
    std::clog << "... DONE: project ZH_0." << std::endl;
    ev->WH->Project("+WH_0", var2, cutmc);
    std::clog << "... DONE: project WH_0." << std::endl;
    
    // Apply scale factors
    cutsf = Form("%f", ev->sf_Wj0b);
    ev->Wj0b->Project("+Wj0b_0", var2, cutmc * cutsf);
    std::clog << "... DONE: project Wj0b_0, scaled by " << cutsf << "." << std::endl;
    cutsf = Form("%f", ev->sf_Wj1b);
    ev->Wj1b->Project("+Wj1b_0", var2, cutmc * cutsf);
    std::clog << "... DONE: project Wj1b_0, scaled by " << cutsf << "." << std::endl;
    cutsf = Form("%f", ev->sf_Wj2b);
    ev->Wj2b->Project("+Wj2b_0", var2, cutmc * cutsf);
    std::clog << "... DONE: project Wj2b_0, scaled by " << cutsf << "." << std::endl;
    cutsf = Form("%f", ev->sf_Zj0b);
    ev->Zj0b->Project("+Zj0b_0", var2, cutmc * cutsf);
    std::clog << "... DONE: project Zj0b_0, scaled by " << cutsf << "." << std::endl;
    cutsf = Form("%f", ev->sf_Zj1b);
    ev->Zj1b->Project("+Zj1b_0", var2, cutmc * cutsf);
    std::clog << "... DONE: project Zj1b_0, scaled by " << cutsf << "." << std::endl;
    cutsf = Form("%f", ev->sf_Zj2b);
    ev->Zj2b->Project("+Zj2b_0", var2, cutmc * cutsf);
    std::clog << "... DONE: project Zj2b_0, scaled by " << cutsf << "." << std::endl;
    cutsf = Form("%f", ev->sf_TT);
    ev->TT->Project("+TT_0", var2, cutmc * cutsf);
    std::clog << "... DONE: project TT_0, scaled by " << cutsf << "." << std::endl;
    
    ev->s_Top->Project("+s_Top_0", var2, cutmc);
    std::clog << "... DONE: project s_Top_0." << std::endl;
    ev->VVLF->Project("+VVLF_0", var2, cutmc);
    std::clog << "... DONE: project VVLF_0." << std::endl;
    ev->VVHF->Project("+VVHF_0", var2, cutmc);
    std::clog << "... DONE: project VVHF_0." << std::endl;
    //ev->QCD->Project("+QCD_0", var2, cutmc);
    //std::clog << "... DONE: project QCD_0." << std::endl;
    
    ev->Wj0b_syst->Project("+Wj0b_syst_0", var2, cutmc);
    std::clog << "... DONE: project Wj0b_syst_0." << std::endl;
    ev->Wj1b_syst->Project("+Wj1b_syst_0", var2, cutmc);
    std::clog << "... DONE: project Wj1b_syst_0." << std::endl;
    ev->Wj2b_syst->Project("+Wj2b_syst_0", var2, cutmc);
    std::clog << "... DONE: project Wj2b_syst_0." << std::endl;
    ev->Zj0b_syst->Project("+Zj0b_syst_0", var2, cutmc);
    std::clog << "... DONE: project Zj0b_syst_0." << std::endl;
    ev->Zj1b_syst->Project("+Zj1b_syst_0", var2, cutmc);
    std::clog << "... DONE: project Zj1b_syst_0." << std::endl;
    ev->Zj2b_syst->Project("+Zj2b_syst_0", var2, cutmc);
    std::clog << "... DONE: project Zj2b_syst_0." << std::endl;
    ev->TT_syst->Project("+TT_syst_0", var2, cutmc);
    std::clog << "... DONE: project TT_syst_0." << std::endl;
    
    //ev->ZH_SM->Project("+ZH_SM_0", var2, cutmc);
    //std::clog << "... DONE: project ZH_SM_0." << std::endl;
    //ev->WH_SM->Project("+WH_SM_0", var2, cutmc);
    //std::clog << "... DONE: project WH_SM_0." << std::endl;
    
    ev->data_obs->Project("+data_obs_0", var2, cutdata);
    std::clog << "... DONE: project data_obs_0." << std::endl;

    // Make sum of histograms
    hVH_0->Add(hZH_0);
    hVH_0->Add(hWH_0);
    std::clog << "... DONE: add ZH_0, WH_0 to VH_0." << std::endl;

    hVV_0->Add(hVVLF_0);
    hVV_0->Add(hVVHF_0);
    std::clog << "... DONE: add VVLF_0, VVHF_0 to VV_0." << std::endl;

    hmc_exp_0->Add(hWj0b_0);
    hmc_exp_0->Add(hWj1b_0);
    hmc_exp_0->Add(hWj2b_0);
    hmc_exp_0->Add(hZj0b_0);
    hmc_exp_0->Add(hZj1b_0);
    hmc_exp_0->Add(hZj2b_0);
    hmc_exp_0->Add(hTT_0);
    hmc_exp_0->Add(hs_Top_0);
    hmc_exp_0->Add(hVV_0);
    //hmc_exp_0->Add(hQCD_0);
    std::clog << "... DONE: add MC backgrounds to mc_exp_0." << std::endl;

    if (rebinner == 0) {
        rebinner = new Rebinner(newnbins, rebinerrorf, xlow, xup);
        rebinner->setSignalAndBackground(hVH_0, hmc_exp_0);
    }

    TH1F * hZH        = rebinner->rebin(hZH_0         , newnbins, "ZH"        );
    TH1F * hWH        = rebinner->rebin(hWH_0         , newnbins, "WH"        );
    TH1F * hWj0b      = rebinner->rebin(hWj0b_0       , newnbins, "Wj0b"      );
    TH1F * hWj1b      = rebinner->rebin(hWj1b_0       , newnbins, "Wj1b"      );
    TH1F * hWj2b      = rebinner->rebin(hWj2b_0       , newnbins, "Wj2b"      );
    TH1F * hZj0b      = rebinner->rebin(hZj0b_0       , newnbins, "Zj0b"      );
    TH1F * hZj1b      = rebinner->rebin(hZj1b_0       , newnbins, "Zj1b"      );
    TH1F * hZj2b      = rebinner->rebin(hZj2b_0       , newnbins, "Zj2b"      );
    TH1F * hTT        = rebinner->rebin(hTT_0         , newnbins, "TT"        );
    TH1F * hs_Top     = rebinner->rebin(hs_Top_0      , newnbins, "s_Top"     );
    TH1F * hVVLF      = rebinner->rebin(hVVLF_0       , newnbins, "VVLF"      );
    TH1F * hVVHF      = rebinner->rebin(hVVHF_0       , newnbins, "VVHF"      );
    TH1F * hQCD       = rebinner->rebin(hQCD_0        , newnbins, "QCD"       );
    TH1F * hWj0b_syst = rebinner->rebin(hWj0b_syst_0  , newnbins, "Wj0b_syst" );
    TH1F * hWj1b_syst = rebinner->rebin(hWj1b_syst_0  , newnbins, "Wj1b_syst" );
    TH1F * hWj2b_syst = rebinner->rebin(hWj2b_syst_0  , newnbins, "Wj2b_syst" );
    TH1F * hZj0b_syst = rebinner->rebin(hZj0b_syst_0  , newnbins, "Zj0b_syst" );
    TH1F * hZj1b_syst = rebinner->rebin(hZj1b_syst_0  , newnbins, "Zj1b_syst" );
    TH1F * hZj2b_syst = rebinner->rebin(hZj2b_syst_0  , newnbins, "Zj2b_syst" );
    TH1F * hTT_syst   = rebinner->rebin(hTT_syst_0    , newnbins, "TT_syst"   );
    //TH1F * hZH_SM     = rebinner->rebin(hZH_SM_0      , newnbins, "ZH_SM"     );
    //TH1F * hWH_SM     = rebinner->rebin(hWH_SM_0      , newnbins, "WH_SM"     );
    TH1F * hVH        = rebinner->rebin(hVH_0         , newnbins, "VH"        );
    TH1F * hVV        = rebinner->rebin(hVV_0         , newnbins, "VV"        );
    TH1F * hmc_exp    = rebinner->rebin(hmc_exp_0     , newnbins, "mc_exp"    );
    TH1F * hdata_obs  = rebinner->rebin(hdata_obs_0   , newnbins, "data_obs"  );

    std::vector<TH1 *> histos;
    histos.push_back(hZH);
    histos.push_back(hWH);
    histos.push_back(hWj0b);
    histos.push_back(hWj1b);
    histos.push_back(hWj2b);
    histos.push_back(hZj0b);
    histos.push_back(hZj1b);
    histos.push_back(hZj2b);
    histos.push_back(hTT);
    histos.push_back(hs_Top);
    histos.push_back(hVVLF);
    histos.push_back(hVVHF);
    histos.push_back(hQCD);
    histos.push_back(hWj0b_syst);
    histos.push_back(hWj1b_syst);
    histos.push_back(hWj2b_syst);
    histos.push_back(hZj0b_syst);
    histos.push_back(hZj1b_syst);
    histos.push_back(hZj2b_syst);
    histos.push_back(hTT_syst);
    //histos.push_back(hZH_SM);
    //histos.push_back(hWH_SM);
    histos.push_back(hVH);
    histos.push_back(hVV);
    histos.push_back(hmc_exp);
    histos.push_back(hdata_obs);
    
    assert(histos.size() == histos_0.size());
    for (UInt_t ih = 0; ih < histos.size(); ih++)
        histos.at(ih)->Sumw2();

    if (printCard && syst == "NONE") {
        std::clog << "MakePlots(): Printing the datacard..." << std::endl;

        TString dcname = g_dcname;
        TString wsname = g_wsname;
        dcname.ReplaceAll("$CHANNEL", channel_CR);
        wsname.ReplaceAll("$CHANNEL", channel_CR);
        const TString channel_8TeV = channel_CR + "_8TeV";
        TString channel_8TeV_1 = channel + "_8TeV";  // scale factors
        if (channel == "ZnunuLowCSV")
            channel_8TeV_1.ReplaceAll("ZnunuLowCSV", "ZnunuHighPt");
        
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
        dc << "bin         "; for (int j=0; j!=jmax+1; j++)  dc << channel_8TeV << " "; dc << std::endl;
        dc << "process     ZH         WH         Wj0b       Wj1b       Wj2b       Zj0b       Zj1b       Zj2b       TT         s_Top      VV         QCD        " << std::endl;
        dc << "process     -1         0          1          2          3          4          5          6          7          8          9          10         " << std::endl;
        dc << "rate        " << setw(10) << right << hZH->Integral() << " " << setw(10) << right << hWH->Integral() << " " << setw(10) << right << hWj0b->Integral() << " " << setw(10) << right << hWj1b->Integral() << " " << setw(10) << right << hWj2b->Integral() << " " << setw(10) << right << hZj0b->Integral() << " " << setw(10) << right << hZj1b->Integral() << " " << setw(10) << right << hZj2b->Integral() << " " << setw(10) << right << hTT->Integral() << " " << setw(10) << right << hs_Top->Integral() << " " << setw(10) << right << hVV->Integral() << " " << setw(10) << right << hQCD->Integral() << std::endl;
        dc.precision(2);
        dc << "-----------------------------------" << std::endl;
        dc << "" << std::endl;
        dc << "### Flat ########################## #####  ZH    WH    Wj0b  Wj1b  Wj2b  Zj0b  Zj1b  Zj2b  TT    s_Top VV    QCD  " << std::endl;
        dc << "CMS_vhbb_ZH_SF_" << channel_8TeV_1 << "     lnN    1.30  -     -     -     -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_WH_SF_" << channel_8TeV_1 << "     lnN    -     1.30  -     -     -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_Wj0b_SF_" << channel_8TeV_1 << "   lnN    -     -     " << scalefactors_lnN[0] << "  -     -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_Wj1b_SF_" << channel_8TeV_1 << "   lnN    -     -     -     " << scalefactors_lnN[1] << "  -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_Wj2b_SF_" << channel_8TeV_1 << "   lnN    -     -     -     -     " << scalefactors_lnN[2] << "  -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_Zj0b_SF_" << channel_8TeV_1 << "   lnN    -     -     -     -     -     " << scalefactors_lnN[3] << "  -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_Zj1b_SF_" << channel_8TeV_1 << "   lnN    -     -     -     -     -     -     " << scalefactors_lnN[4] << "  -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_Zj2b_SF_" << channel_8TeV_1 << "   lnN    -     -     -     -     -     -     -     " << scalefactors_lnN[5] << "  -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_TT_SF_" << channel_8TeV_1 << "     lnN    -     -     -     -     -     -     -     -     " << scalefactors_lnN[6] << "  -     -     -    " << std::endl;
        dc << "CMS_vhbb_s_Top_SF_" << channel_8TeV_1 << "  lnN    -     -     -     -     -     -     -     -     -     1.25  -     -    " << std::endl;
        dc << "CMS_vhbb_VV_SF_" << channel_8TeV_1 << "     lnN    -     -     -     -     -     -     -     -     -     -     1.30  -    " << std::endl;
        dc << "#CMS_vhbb_QCD_SF_" << channel_8TeV_1 << "   lnN    -     -     -     -     -     -     -     -     -     -     -     2.00 " << std::endl;
        dc << "### Shape ######################### #####  ZH    WH    Wj0b  Wj1b  Wj2b  Zj0b  Zj1b  Zj2b  TT    s_Top VV    QCD  " << std::endl;
        dc << "UEPS                                shape  1.00  1.00  -     -     -     -     -     -     -     1.00  1.00  -    " << std::endl;
        dc << "CMS_vhbb_eff_b_SIG                  shape  1.00  1.00  -     -     -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_eff_b                      shape  -     -     1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  -    " << std::endl;
        dc << "CMS_vhbb_fake_b_8TeV                shape  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  -    " << std::endl;
        dc << "CMS_vhbb_res_j                      shape  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  -    " << std::endl;
        dc << "CMS_vhbb_scale_j                    shape  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  -    " << std::endl;
        dc << "CMS_vhbb_trigger_MET_Znunu_8TeV     shape  1.00  1.00  -     -     -     -     -     -     -     1.00  1.00  -    " << std::endl;
        if (channel != "ZnunuLowPt") {
        dc << "CMS_vhbb_WJModel_Znunu_8TeV         shape  -     -     1.00  1.00  1.00  -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_ZJModel_Znunu_8TeV         shape  -     -     -     -     -     1.00  1.00  1.00  -     -     -     -    " << std::endl;
        }
        dc << "CMS_vhbb_TTModel_Znunu_8TeV         shape  -     -     -     -     -     -     -     -     1.00  -     -     -    " << std::endl;
        dc << "################################### #####  ZH    WH    Wj0b  Wj1b  Wj2b  Zj0b  Zj1b  Zj2b  TT    s_Top VV    QCD  " << std::endl;

        dc.close();
        
        std::clog << "MakePlots(): The datacard is written." << std::endl;
    } // end if syst == NONE

    if (writeRoot) {
        std::clog << "MakePlots(): Writing out histograms..." << std::endl;
        TFile * out(0);
        TString rootname = g_rootname;
        rootname.ReplaceAll("$CHANNEL", channel_CR);
        if (syst == "NONE") {
            out = TFile::Open(rootname, "RECREATE");
            out->mkdir(channel_CR);
        } else {
            out = TFile::Open(rootname, "UPDATE");
        }
        out->cd(channel_CR);

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

        // Add histograms for stat, WjModel, ZjModel, TTModel
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
            
            name = "Wj0b_CMS_vhbb_WJModel_Znunu_8TeVUp";
            TH1F * hWj0b_modelUp = VaryModelErrors((TH1F *) hWj0b, (TH1F *) hWj0b_syst, 1.0, name);
            name = "Wj1b_CMS_vhbb_WJModel_Znunu_8TeVUp";
            TH1F * hWj1b_modelUp = VaryModelErrors((TH1F *) hWj1b, (TH1F *) hWj1b_syst, 1.0, name);
            name = "Wj2b_CMS_vhbb_WJModel_Znunu_8TeVUp";
            TH1F * hWj2b_modelUp = VaryModelErrors((TH1F *) hWj2b, (TH1F *) hWj2b_syst, 1.0, name);
            name = "Wj0b_CMS_vhbb_WJModel_Znunu_8TeVDown";
            TH1F * hWj0b_modelDown = VaryModelErrors((TH1F *) hWj0b, (TH1F *) hWj0b_syst, -1.0, name);
            name = "Wj1b_CMS_vhbb_WJModel_Znunu_8TeVDown";
            TH1F * hWj1b_modelDown = VaryModelErrors((TH1F *) hWj1b, (TH1F *) hWj1b_syst, -1.0, name);
            name = "Wj2b_CMS_vhbb_WJModel_Znunu_8TeVDown";
            TH1F * hWj2b_modelDown = VaryModelErrors((TH1F *) hWj2b, (TH1F *) hWj2b_syst, -1.0, name);
            hWj0b_modelUp->Write(hWj0b_modelUp->GetName());
            hWj1b_modelUp->Write(hWj1b_modelUp->GetName());
            hWj2b_modelUp->Write(hWj2b_modelUp->GetName());
            hWj0b_modelDown->Write(hWj0b_modelDown->GetName());
            hWj1b_modelDown->Write(hWj1b_modelDown->GetName());
            hWj2b_modelDown->Write(hWj2b_modelDown->GetName());
            
            delete hWj0b_modelUp;
            delete hWj1b_modelUp;
            delete hWj2b_modelUp;
            delete hWj0b_modelDown;
            delete hWj1b_modelDown;
            delete hWj2b_modelDown;
            
            name = "Zj0b_CMS_vhbb_ZJModel_Znunu_8TeVUp";
            TH1F * hZj0b_modelUp = VaryModelErrors((TH1F *) hZj0b, (TH1F *) hZj0b_syst, 1.0, name);
            name = "Zj1b_CMS_vhbb_ZJModel_Znunu_8TeVUp";
            TH1F * hZj1b_modelUp = VaryModelErrors((TH1F *) hZj1b, (TH1F *) hZj1b_syst, 1.0, name);
            name = "Zj2b_CMS_vhbb_ZJModel_Znunu_8TeVUp";
            TH1F * hZj2b_modelUp = VaryModelErrors((TH1F *) hZj2b, (TH1F *) hZj2b_syst, 1.0, name);
            name = "Zj0b_CMS_vhbb_ZJModel_Znunu_8TeVDown";
            TH1F * hZj0b_modelDown = VaryModelErrors((TH1F *) hZj0b, (TH1F *) hZj0b_syst, -1.0, name);
            name = "Zj1b_CMS_vhbb_ZJModel_Znunu_8TeVDown";
            TH1F * hZj1b_modelDown = VaryModelErrors((TH1F *) hZj1b, (TH1F *) hZj1b_syst, -1.0, name);
            name = "Zj2b_CMS_vhbb_ZJModel_Znunu_8TeVDown";
            TH1F * hZj2b_modelDown = VaryModelErrors((TH1F *) hZj2b, (TH1F *) hZj2b_syst, -1.0, name);
            hZj0b_modelUp->Write(hZj0b_modelUp->GetName());
            hZj1b_modelUp->Write(hZj1b_modelUp->GetName());
            hZj2b_modelUp->Write(hZj2b_modelUp->GetName());
            hZj0b_modelDown->Write(hZj0b_modelDown->GetName());
            hZj1b_modelDown->Write(hZj1b_modelDown->GetName());
            hZj2b_modelDown->Write(hZj2b_modelDown->GetName());
            
            delete hZj0b_modelUp;
            delete hZj1b_modelUp;
            delete hZj2b_modelUp;
            delete hZj0b_modelDown;
            delete hZj1b_modelDown;
            delete hZj2b_modelDown;
            
            name = "TT_CMS_vhbb_TTModel_Znunu_8TeVUp";
            TH1F * hTT_modelUp = VaryModelErrors((TH1F *) hTT, (TH1F *) hTT_syst, 1.0, name);
            name = "TT_CMS_vhbb_TTModel_Znunu_8TeVDown";
            TH1F * hTT_modelDown = VaryModelErrors((TH1F *) hTT, (TH1F *) hTT_syst, -1.0, name);
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


EventsJ11 * Read(bool isControlRegion, const TCut cutallmc, const TCut cutalldata)
{
    std::clog << "Read(): Using cutallmc = " << cutallmc << ", cutalldata = " << cutalldata << std::endl;
    std::clog << "Read():       indir = " << indir << std::endl;
    
    const TString treename = Form("tree_%s_%s", channel.Data(), (isControlRegion ? "ctrl" : "test"));
    
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
    
    TChain Wj_syst(treename);
    Wj_syst.Add(indir + prefix + "WjHW" + suffix);
    
    TChain Zj_syst(treename);
    Zj_syst.Add(indir + prefix + "ZjHW" + suffix);
    
    TChain TT_syst(treename);
    TT_syst.Add(indir + prefix + "TTPowheg" + suffix);
    
    TChain ZH_SM(treename);
    ZH_SM.Add(indir + prefix + "ZH125" + suffix);
    
    TChain WH_SM(treename);
    WH_SM.Add(indir + prefix + "WH125" + suffix);
    
    TChain data_obs(treename);
    data_obs.Add(indir + prefix + "data_obs" + suffix);

    // Start copying trees
    const TCut cutHF = "eventFlav==5";
    const TCut cutLF = "eventFlav!=5";

    const TCut cut2b = "abs(hJet_flavour[0])==5 && abs(hJet_flavour[1])==5";
    const TCut cut1b = "(abs(hJet_flavour[0])==5 && abs(hJet_flavour[1])!=5) || (abs(hJet_flavour[0])!=5 && abs(hJet_flavour[1])==5)";
    const TCut cut0b = "abs(hJet_flavour[0])!=5 && abs(hJet_flavour[1])!=5";
    
    EventsJ11 * ev = new EventsJ11();
    ev->ZH = (TTree *) ZH.CopyTree(cutallmc);
    std::clog << "... DONE: ZH copy tree." << std::endl;
    ev->WH = (TTree *) WH.CopyTree(cutallmc);
    std::clog << "... DONE: WH copy tree." << std::endl;
    ev->Wj0b = (TTree *) Wj.CopyTree(cutallmc + cut0b);
    std::clog << "... DONE: Wj0b copy tree." << std::endl;
    ev->Wj1b = (TTree *) Wj.CopyTree(cutallmc + cut1b);
    std::clog << "... DONE: Wj1b copy tree." << std::endl;
    ev->Wj2b = (TTree *) Wj.CopyTree(cutallmc + cut2b);
    std::clog << "... DONE: Wj2b copy tree." << std::endl;
    ev->Zj0b = (TTree *) Zj.CopyTree(cutallmc + cut0b);
    std::clog << "... DONE: Zj0b copy tree." << std::endl;
    ev->Zj1b = (TTree *) Zj.CopyTree(cutallmc + cut1b);
    std::clog << "... DONE: Zj1b copy tree." << std::endl;
    ev->Zj2b = (TTree *) Zj.CopyTree(cutallmc + cut2b);
    std::clog << "... DONE: Zj2b copy tree." << std::endl;
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
    ev->Wj0b_syst = (TTree *) Wj_syst.CopyTree(cutallmc + cut0b);
    std::clog << "... DONE: Wj0b_syst copy tree." << std::endl;
    ev->Wj1b_syst = (TTree *) Wj_syst.CopyTree(cutallmc + cut1b);
    std::clog << "... DONE: Wj1b_syst copy tree." << std::endl;
    ev->Wj2b_syst = (TTree *) Wj_syst.CopyTree(cutallmc + cut2b);
    std::clog << "... DONE: Wj2b_syst copy tree." << std::endl;
    ev->Zj0b_syst = (TTree *) Zj_syst.CopyTree(cutallmc + cut0b);
    std::clog << "... DONE: Zj0b_syst copy tree." << std::endl;
    ev->Zj1b_syst = (TTree *) Zj_syst.CopyTree(cutallmc + cut1b);
    std::clog << "... DONE: Zj1b_syst copy tree." << std::endl;
    ev->Zj2b_syst = (TTree *) Zj_syst.CopyTree(cutallmc + cut2b);
    std::clog << "... DONE: Zj2b_syst copy tree." << std::endl;
    ev->TT_syst = (TTree *) TT_syst.CopyTree(cutallmc);
    std::clog << "... DONE: TT_syst copy tree." << std::endl;
    //ev->ZH_SM = (TTree *) ZH_SM.CopyTree(cutallmc);
    //std::clog << "... DONE: ZH_SM copy tree." << std::endl;
    //ev->WH_SM = (TTree *) WH_SM.CopyTree(cutallmc);
    //std::clog << "... DONE: WH_SM copy tree." << std::endl;
    
    ev->data_obs = (TTree *) data_obs.CopyTree(cutalldata);
    std::clog << "... DONE: data_obs copy tree." << std::endl;

    return ev;
}

void calc_scalefactors(int n, int niterations, double * scalefactors, double * scalefactors_lnN, double * fitresults, double * fitresults_lnN) {
    niterations += 1;
    for (int i=0; i<n; i++) {
        scalefactors[i]     = g_scalefactors_prefit[i];
        scalefactors_lnN[i] = g_scalefactors_lnN_prefit[i] - 1.0;
        for (int j=0; j<niterations; j++) {
            // scale factors are updated first before their uncertainties;
            scalefactors[i]     += (fitresults[n*j + i] * scalefactors_lnN[i]);
            scalefactors_lnN[i] *= fitresults_lnN[n*j + i];
        }
        scalefactors_lnN[i] += 1.0;
    }
}


////////////////////////////////////////////////////////////////////////////////
/// Main                                                                     ///
////////////////////////////////////////////////////////////////////////////////

void ScaleFactorJ11(int nbins=500, long long newnbins=15, double rebinerrorf=0.35, bool freerebin=false, bool isControlRegion=true)
{
    gROOT->LoadMacro("tdrstyle.C");
    gROOT->ProcessLine("setTDRStyle()");
    gROOT->LoadMacro("HelperFunctions.h");  // FIXME
    
    TH1::SetDefaultSumw2(1);
    gROOT->SetBatch(1);
    
    //if (!TString(gROOT->GetVersion()).Contains("5.34")) {
    //    std::cout << "INCORRECT ROOT VERSION! Please use 5.34:" << std::endl;
    //    std::cout << "source /uscmst1/prod/sw/cms/slc5_amd64_gcc472/lcg/root/5.34.03-cms4/bin/thisroot.csh" << std::endl;
    //    std::cout << "Return without doing anything." << std::endl;
    //    return;
    //}
    
    double fitresults_ZnunuHighPt[]     = {
         0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
        -0.917,  0.040, -0.897, -1.268, -0.310, -0.505, -1.698, 
    };
    double fitresults_lnN_ZnunuHighPt[] = {
         1.000,  1.000,  1.000,  1.000,  1.000,  1.000,  1.000,
         0.334,  0.637,  0.735,  0.380,  0.347,  0.526,  0.450, 
    };
    
    double fitresults_ZnunuMedPt[]      = {
         0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
        -0.796,  0.614,  0.591,  1.801,  1.311,  0.146, -1.141, 
    };
    double fitresults_lnN_ZnunuMedPt[]  = {
         1.000,  1.000,  1.000,  1.000,  1.000,  1.000,  1.000,
         0.261,  0.627,  0.617,  0.415,  0.346,  0.629,  0.349, 
    };
    
    double fitresults_ZnunuLowPt[]      = {
         0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
        -3.389,  2.070,  1.110,  1.440,  0.177,  0.514, -1.345, 
    };
    double fitresults_lnN_ZnunuLowPt[]  = {
         1.000,  1.000,  1.000,  1.000,  1.000,  1.000,  1.000,
         0.229,  0.300,  0.652,  0.468,  0.396,  0.562,  0.345, 
    };
    
    double fitresults_ZnunuLowCSV[]     = {
         0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
    };
    double fitresults_lnN_ZnunuLowCSV[] = { 
         1.000,  1.000,  1.000,  1.000,  1.000,  1.000,  1.000,
    };
    
    double scalefactors_ZnunuHighPt[7], scalefactors_lnN_ZnunuHighPt[7];
    double scalefactors_ZnunuMedPt [7], scalefactors_lnN_ZnunuMedPt [7];
    double scalefactors_ZnunuLowPt [7], scalefactors_lnN_ZnunuLowPt [7];
    double scalefactors_ZnunuLowCSV[7], scalefactors_lnN_ZnunuLowCSV[7];
    // 0 is the first iteration with prefit scale factors and uncertainties
    calc_scalefactors(7, 1, scalefactors_ZnunuHighPt, scalefactors_lnN_ZnunuHighPt, fitresults_ZnunuHighPt, fitresults_lnN_ZnunuHighPt);
    calc_scalefactors(7, 1, scalefactors_ZnunuMedPt , scalefactors_lnN_ZnunuMedPt , fitresults_ZnunuMedPt , fitresults_lnN_ZnunuMedPt );
    calc_scalefactors(7, 1, scalefactors_ZnunuLowPt , scalefactors_lnN_ZnunuLowPt , fitresults_ZnunuLowPt , fitresults_lnN_ZnunuLowPt );
    calc_scalefactors(7, 0, scalefactors_ZnunuLowCSV, scalefactors_lnN_ZnunuLowCSV, fitresults_ZnunuLowCSV, fitresults_lnN_ZnunuLowCSV);
    

    if (gSystem->AccessPathName(g_plotdir))
        gSystem->mkdir(g_plotdir);

    const bool makeslice = (newnbins > 100);
    if (!makeslice)  assert(nbins > newnbins);
    const int begin = (isControlRegion) ?   1 : 0;
    const int end   = (isControlRegion) ? 1+5 : 1;

    const TString maxCSV        = "max(hJet_csv_nominal[0],hJet_csv_nominal[1])";
    const TString maxCSV_upBC   = "max(hJet_csv_upBC[0],hJet_csv_upBC[1])";
    const TString maxCSV_downBC = "max(hJet_csv_downBC[0],hJet_csv_downBC[1])";
    const TString maxCSV_upL    = "max(hJet_csv_upL[0],hJet_csv_upL[1])";
    const TString maxCSV_downL  = "max(hJet_csv_downL[0],hJet_csv_downL[1])";
    const TString minCSV        = "1.08+min(hJet_csv_nominal[0],hJet_csv_nominal[1])";
    const TString minCSV_upBC   = "1.08+min(hJet_csv_upBC[0],hJet_csv_upBC[1])";
    const TString minCSV_downBC = "1.08+min(hJet_csv_downBC[0],hJet_csv_downBC[1])";
    const TString minCSV_upL    = "1.08+min(hJet_csv_upL[0],hJet_csv_upL[1])";
    const TString minCSV_downL  = "1.08+min(hJet_csv_downL[0],hJet_csv_downL[1])";


    ///-- ZnunuHighPt ----------------------------------------------------------
    channel = "ZnunuHighPt";
    std::clog << "\n*** " << channel << " ***\n" << std::endl;

    EventsJ11 * ev = Read(isControlRegion, g_cutallmc, g_cutalldata);
    ev->check();

    /// Select scale factors
    const double * scalefactors     = scalefactors_ZnunuHighPt;
    const double * scalefactors_lnN = scalefactors_lnN_ZnunuHighPt;

    for (int ireg=begin; ireg<end; ireg++) {
        channel_CR = channel + "_" + g_regions[ireg];

        /// Reset rebinner
        delete rebinner; rebinner = 0;

        /// Loop over systematics
        for (int isyst = 0; isyst < nsyst; isyst++) {
            ev->set_scale_factors(scalefactors);
            std::clog << "isyst " << isyst << ": scale factors: " << ev->sf_Wj0b << " (Wj0b), " << ev->sf_Wj1b << " (Wj1b), " << ev->sf_Wj2b << " (Wj2b), " << ev->sf_Zj0b << " (Zj0b), " << ev->sf_Zj1b << " (Zj1b), " << ev->sf_Zj2b << " (Zj2b), " << ev->sf_TT << " (TT)." << std::endl;
            
            TString var = maxCSV;
            if (isyst == 5)       var = maxCSV_upBC;
            else if (isyst == 6)  var = maxCSV_downBC;
            else if (isyst == 7)  var = maxCSV_upL;
            else if (isyst == 8)  var = maxCSV_downL;
            
            TString var2 = minCSV;
            if (isyst == 5)       var2 = minCSV_upBC;
            else if (isyst == 6)  var2 = minCSV_downBC;
            else if (isyst == 7)  var2 = minCSV_upL;
            else if (isyst == 8)  var2 = minCSV_downL;
            
            TCut cutmc_ctrl = Form("weightsHCP[%i] * (selectFlags[%i][%i]) * triggercorr2012ABCD( (triggerFlags[42]==1 || triggerFlags[39]==1), (triggerFlags[41]==1), METtype1corr.et, max(hJet_csv_nominal[0],hJet_csv_nominal[1]) )", isyst, ireg, isyst);  // FIXME
            //TCut cutmc_ctrl = Form("weightsHCP[%i] * (selectFlags[%i][%i])", isyst, ireg, isyst);
            TCut cutdata_ctrl = Form("(selectFlags[%i][%i])", ireg, isyst);
            
            MakePlots(ev, var, var2, cutmc_ctrl, cutdata_ctrl, g_systematics[isyst], g_regions[ireg], var.Data(), 60, 0.00, 2.16, 60, rebinerrorf, scalefactors_lnN, "!printStat:printCard:writeRoot:!plotData:plotLog:plotSig");
            std::clog << "isyst " << isyst << ": DONE!\n" << std::endl;
        }
    }  // end loop over regions
    delete ev;


    ///-- ZnunuMedPt -----------------------------------------------------------
    channel = "ZnunuMedPt";
    std::clog << "\n*** " << channel << " ***\n" << std::endl;

    ev = Read(isControlRegion, g_cutallmc, g_cutalldata);
    ev->check();

    /// Select scale factors
    scalefactors     = scalefactors_ZnunuMedPt;
    scalefactors_lnN = scalefactors_lnN_ZnunuMedPt;

    for (int ireg=begin; ireg<end; ireg++) {
        channel_CR = channel + "_" + g_regions[ireg];

        /// Reset rebinner
        delete rebinner; rebinner = 0;

        /// Loop over systematics
        for (int isyst = 0; isyst < nsyst; isyst++) {
            ev->set_scale_factors(scalefactors);
            std::clog << "isyst " << isyst << ": scale factors: " << ev->sf_Wj0b << " (Wj0b), " << ev->sf_Wj1b << " (Wj1b), " << ev->sf_Wj2b << " (Wj2b), " << ev->sf_Zj0b << " (Zj0b), " << ev->sf_Zj1b << " (Zj1b), " << ev->sf_Zj2b << " (Zj2b), " << ev->sf_TT << " (TT)." << std::endl;
            
            TString var = maxCSV;
            if (isyst == 5)       var = maxCSV_upBC;
            else if (isyst == 6)  var = maxCSV_downBC;
            else if (isyst == 7)  var = maxCSV_upL;
            else if (isyst == 8)  var = maxCSV_downL;
            
            TString var2 = minCSV;
            if (isyst == 5)       var2 = minCSV_upBC;
            else if (isyst == 6)  var2 = minCSV_downBC;
            else if (isyst == 7)  var2 = minCSV_upL;
            else if (isyst == 8)  var2 = minCSV_downL;
            
            TCut cutmc_ctrl = Form("weightsHCP[%i] * (selectFlags[%i][%i]) * triggercorr2012ABCD( (triggerFlags[42]==1 || triggerFlags[39]==1), (triggerFlags[41]==1), METtype1corr.et, max(hJet_csv_nominal[0],hJet_csv_nominal[1]) )", isyst, ireg, isyst);  // FIXME
            //TCut cutmc_ctrl = Form("weightsHCP[%i] * (selectFlags[%i][%i])", isyst, ireg, isyst);
            TCut cutdata_ctrl = Form("(selectFlags[%i][%i])", ireg, isyst);
            
            MakePlots(ev, var, var2, cutmc_ctrl, cutdata_ctrl, g_systematics[isyst], g_regions[ireg], var.Data(), 60, 0.00, 2.16, 60, rebinerrorf, scalefactors_lnN, "!printStat:printCard:writeRoot:!plotData:plotLog:plotSig");
            std::clog << "isyst " << isyst << ": DONE!\n" << std::endl;
        }
    }  // end loop over regions
    delete ev;


    ///-- ZnunuLowPt -----------------------------------------------------------
    channel = "ZnunuLowPt";
    std::clog << "\n*** " << channel << " ***\n" << std::endl;

    ev = Read(isControlRegion, g_cutallmc, g_cutalldata);
    ev->check();

    /// Select scale factors
    scalefactors     = scalefactors_ZnunuLowPt;
    scalefactors_lnN = scalefactors_lnN_ZnunuLowPt;

    for (int ireg=begin; ireg<end; ireg++) {
        channel_CR = channel + "_" + g_regions[ireg];

        /// Reset rebinner
        delete rebinner; rebinner = 0;

        /// Loop over systematics
        for (int isyst = 0; isyst < nsyst; isyst++) {
            ev->set_scale_factors(scalefactors);
            std::clog << "isyst " << isyst << ": scale factors: " << ev->sf_Wj0b << " (Wj0b), " << ev->sf_Wj1b << " (Wj1b), " << ev->sf_Wj2b << " (Wj2b), " << ev->sf_Zj0b << " (Zj0b), " << ev->sf_Zj1b << " (Zj1b), " << ev->sf_Zj2b << " (Zj2b), " << ev->sf_TT << " (TT)." << std::endl;
            
            TString var = maxCSV;
            if (isyst == 5)       var = maxCSV_upBC;
            else if (isyst == 6)  var = maxCSV_downBC;
            else if (isyst == 7)  var = maxCSV_upL;
            else if (isyst == 8)  var = maxCSV_downL;
            
            TString var2 = minCSV;
            if (isyst == 5)       var2 = minCSV_upBC;
            else if (isyst == 6)  var2 = minCSV_downBC;
            else if (isyst == 7)  var2 = minCSV_upL;
            else if (isyst == 8)  var2 = minCSV_downL;
            
            TCut cutmc_ctrl = Form("weightsHCP[%i] * (selectFlags[%i][%i]) * triggercorr2012ABCD( (triggerFlags[42]==1 || triggerFlags[39]==1), (triggerFlags[41]==1), METtype1corr.et, max(hJet_csv_nominal[0],hJet_csv_nominal[1]) )", isyst, ireg, isyst);  // FIXME
            //TCut cutmc_ctrl = Form("weightsHCP[%i] * (selectFlags[%i][%i])", isyst, ireg, isyst);
            TCut cutdata_ctrl = Form("(selectFlags[%i][%i])", ireg, isyst);
            
            MakePlots(ev, var, var2, cutmc_ctrl, cutdata_ctrl, g_systematics[isyst], g_regions[ireg], var.Data(), 60, 0.00, 2.16, 60, rebinerrorf, scalefactors_lnN, "!printStat:printCard:writeRoot:!plotData:plotLog:plotSig");
            std::clog << "isyst " << isyst << ": DONE!\n" << std::endl;
        }
    }  // end loop over regions
    delete ev;

/*
    ///-- ZnunuLowCSV ----------------------------------------------------------
    channel = "ZnunuLowCSV";
    std::clog << "\n*** " << channel << " ***\n" << std::endl;

    ev = Read(isControlRegion, g_cutallmc, g_cutalldata);
    ev->check();

    /// Select scale factors
    scalefactors     = scalefactors_ZnunuLowCSV;
    scalefactors_lnN = scalefactors_lnN_ZnunuLowCSV;

    for (int ireg=begin; ireg<end; ireg++) {
        channel_CR = channel + "_" + g_regions[ireg];
        if (ireg == 1)  // only two control regions for ZnunuLowCSV
            channel_CR = channel + "_" + g_regions[2];
        else if (ireg == 2)
            channel_CR = channel + "_" + g_regions[4];
        else
            continue;

        /// Reset rebinner
        delete rebinner; rebinner = 0;

        /// Loop over systematics
        for (int isyst = 0; isyst < nsyst; isyst++) {
            ev->set_scale_factors(scalefactors);
            std::clog << "isyst " << isyst << ": scale factors: " << ev->sf_Wj0b << " (Wj0b), " << ev->sf_Wj1b << " (Wj1b), " << ev->sf_Wj2b << " (Wj2b), " << ev->sf_Zj0b << " (Zj0b), " << ev->sf_Zj1b << " (Zj1b), " << ev->sf_Zj2b << " (Zj2b), " << ev->sf_TT << " (TT)." << std::endl;
            
            TString var = maxCSV;
            if (isyst == 5)       var = maxCSV_upBC;
            else if (isyst == 6)  var = maxCSV_downBC;
            else if (isyst == 7)  var = maxCSV_upL;
            else if (isyst == 8)  var = maxCSV_downL;
            
            TString var2 = minCSV;
            if (isyst == 5)       var2 = minCSV_upBC;
            else if (isyst == 6)  var2 = minCSV_downBC;
            else if (isyst == 7)  var2 = minCSV_upL;
            else if (isyst == 8)  var2 = minCSV_downL;
            
            TCut cutmc_ctrl = Form("weightsHCP[%i] * (selectFlags[%i][%i])", isyst, ireg, isyst);
            TCut cutdata_ctrl = Form("(selectFlags[%i][%i])", ireg, isyst);
            
            MakePlots(ev, var, var2, cutmc_ctrl, cutdata_ctrl, g_systematics[isyst], g_regions[ireg], var.Data(), 60, 0.00, 2.16, 60, rebinerrorf, scalefactors_lnN, "!printStat:printCard:writeRoot:!plotData:plotLog:plotSig");
            std::clog << "isyst " << isyst << ": DONE!\n" << std::endl;
        }
    }  // end loop over regions
    delete ev;
*/

    std::cout << "Final scale factors: " << std::endl;
    std::cout << fixed << setprecision(3) << "ZnunuHighPt: " << std::endl;
    std::cout << scalefactors_ZnunuHighPt[0] << ", " << scalefactors_ZnunuHighPt[1] << ", " << scalefactors_ZnunuHighPt[2] << ", " << scalefactors_ZnunuHighPt[3] << ", " << scalefactors_ZnunuHighPt[4] << ", " << scalefactors_ZnunuHighPt[5] << ", " << scalefactors_ZnunuHighPt[6] << std::endl;
    std::cout << fixed << setprecision(2) << "ZnunuHighPt lnN: " << std::endl;
    std::cout << scalefactors_lnN_ZnunuHighPt[0] << ", " << scalefactors_lnN_ZnunuHighPt[1] << ", " << scalefactors_lnN_ZnunuHighPt[2] << ", " << scalefactors_lnN_ZnunuHighPt[3] << ", " << scalefactors_lnN_ZnunuHighPt[4] << ", " << scalefactors_lnN_ZnunuHighPt[5] << ", " << scalefactors_lnN_ZnunuHighPt[6] << std::endl;
    std::cout << fixed << setprecision(3) << "ZnunuMedPt: " << std::endl;
    std::cout << scalefactors_ZnunuMedPt[0] << ", " << scalefactors_ZnunuMedPt[1] << ", " << scalefactors_ZnunuMedPt[2] << ", " << scalefactors_ZnunuMedPt[3] << ", " << scalefactors_ZnunuMedPt[4] << ", " << scalefactors_ZnunuMedPt[5] << ", " << scalefactors_ZnunuMedPt[6] << std::endl;
    std::cout << fixed << setprecision(2) << "ZnunuMedPt lnN: " << std::endl;
    std::cout << scalefactors_lnN_ZnunuMedPt[0] << ", " << scalefactors_lnN_ZnunuMedPt[1] << ", " << scalefactors_lnN_ZnunuMedPt[2] << ", " << scalefactors_lnN_ZnunuMedPt[3] << ", " << scalefactors_lnN_ZnunuMedPt[4] << ", " << scalefactors_lnN_ZnunuMedPt[5] << ", " << scalefactors_lnN_ZnunuMedPt[6] << std::endl;
    std::cout << fixed << setprecision(3) << "ZnunuLowPt: " << std::endl;
    std::cout << scalefactors_ZnunuLowPt[0] << ", " << scalefactors_ZnunuLowPt[1] << ", " << scalefactors_ZnunuLowPt[2] << ", " << scalefactors_ZnunuLowPt[3] << ", " << scalefactors_ZnunuLowPt[4] << ", " << scalefactors_ZnunuLowPt[5] << ", " << scalefactors_ZnunuLowPt[6] << std::endl;
    std::cout << fixed << setprecision(2) << "ZnunuLowPt lnN: " << std::endl;
    std::cout << scalefactors_lnN_ZnunuLowPt[0] << ", " << scalefactors_lnN_ZnunuLowPt[1] << ", " << scalefactors_lnN_ZnunuLowPt[2] << ", " << scalefactors_lnN_ZnunuLowPt[3] << ", " << scalefactors_lnN_ZnunuLowPt[4] << ", " << scalefactors_lnN_ZnunuLowPt[5] << ", " << scalefactors_lnN_ZnunuLowPt[6] << std::endl;
    std::cout << fixed << setprecision(3) << "ZnunuLowCSV: " << std::endl;
    std::cout << scalefactors_ZnunuLowCSV[0] << ", " << scalefactors_ZnunuLowCSV[1] << ", " << scalefactors_ZnunuLowCSV[2] << ", " << scalefactors_ZnunuLowCSV[3] << ", " << scalefactors_ZnunuLowCSV[4] << ", " << scalefactors_ZnunuLowCSV[5] << ", " << scalefactors_ZnunuLowCSV[6] << std::endl;
    std::cout << fixed << setprecision(2) << "ZnunuLowCSV lnN: " << std::endl;
    std::cout << scalefactors_lnN_ZnunuLowCSV[0] << ", " << scalefactors_lnN_ZnunuLowCSV[1] << ", " << scalefactors_lnN_ZnunuLowCSV[2] << ", " << scalefactors_lnN_ZnunuLowCSV[3] << ", " << scalefactors_lnN_ZnunuLowCSV[4] << ", " << scalefactors_lnN_ZnunuLowCSV[5] << ", " << scalefactors_lnN_ZnunuLowCSV[6] << std::endl;
}
