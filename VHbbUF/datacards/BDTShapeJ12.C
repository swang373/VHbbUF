////////////////////////////////////////////////////////////////////////////////
/// This macro runs on Step 4 with jmax = 12 in the following order
///   ZH         WH         Wj0b       Wj1b       Wj2b       Zj0b       Zj1b       Zj2b       TT         s_Top      VVLF       VVHF       QCD
///   -1         0          1          2          3          4          5          6          7          8          9          10         11
///
/// It produces a datacard (set by dcname) and a root file with TH1F (set by rootname).
///
/// The following 1 + 2 x 14 systematic variations are considered:
///   NONE
///   CMS_vhbb_Znn_res_jUp
///   CMS_vhbb_Znn_res_jDown
///   CMS_vhbb_Znn_scale_jUp
///   CMS_vhbb_Znn_scale_jDown
///   CMS_vhbb_eff_bUp
///   CMS_vhbb_eff_bDown
///   CMS_vhbb_fake_b_8TeVUp
///   CMS_vhbb_fake_b_8TeVDown
///   UEPSUp
///   UEPSDown
///   CMS_vhbb_trigger_MET_Znn_8TeVUp
///   CMS_vhbb_trigger_MET_Znn_8TeVDown
///   CMS_vhbb_trigger_CSV_Znn_8TeVUp
///   CMS_vhbb_trigger_CSV_Znn_8TeVDown
///   CMS_vhbb_stat$PROCESS_$CHANNELUp
///   CMS_vhbb_stat$PROCESS_$CHANNELDown
///   CMS_vhbb_WJModel_Znn_8TeVUp
///   CMS_vhbb_WJModel_Znn_8TeVDown
///   CMS_vhbb_ZJModel_Znn_8TeVUp
///   CMS_vhbb_ZJModel_Znn_8TeVDown
///   CMS_vhbb_TTModel_Znn_8TeVUp
///   CMS_vhbb_TTModel_Znn_8TeVDown
///   CMS_vhbb_WJSlope_Znn_8TeVUp
///   CMS_vhbb_WJSlope_Znn_8TeVDown
///   CMS_vhbb_ZJSlope_Znn_8TeVUp
///   CMS_vhbb_ZJSlope_Znn_8TeVDown
///   CMS_vhbb_boost_QCD_8TeVUp
///   CMS_vhbb_boost_QCD_8TeVDown
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
#include "TPRegexp.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

#include "HelperBDTShape.h"
#include "HelperFunctions.h"

//#define MJJANALYSIS

//#define VVANALYSIS

//#define HCPANALYSIS

#define QCDSHAPE

#define VHEWKCORRECTION

#define VHQCDCORRECTION


////////////////////////////////////////////////////////////////////////////////
/// Configuration                                                            ///
////////////////////////////////////////////////////////////////////////////////

const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130404/reload_20130401/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130326/reload_20130401/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130326/reload_20130401_VV/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130326/reload_20130327/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130314/reload_20130323/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130302/reload_20130306/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130302/reload_20130302/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130221/reload_20130228/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130214/reload_20130216/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130404/stitch/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130326/stitch/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130314/stitch/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130302/stitch/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130221/stitch/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130214/stitch/";
//const TString indir     = "root://eoscms//eos/cms/store/caf/user/degrutto/2013/Step4_20130207/";
const TString prefix    = "Step4_";
const TString suffix    = ".root";

const int jmax              = 12;  // as in datacard
const int nsyst             = 1 + 2 * 7;
const int massH             = 125;

//const TCut    g_cutallmc    = "";          // used in CopyTree()
//const TCut    g_cutalldata  = g_cutallmc;  // used in CopyTree()
const TCut    g_cutallmc    = "VtypeWithTau==4";  // FIXME
const TCut    g_cutalldata  = "VtypeWithTau==4 && (!(207883<=EVENT.run && EVENT.run<=208307))";  // FIXME
//const TCut    g_cutallmc    = "(HmassReg<110 || HmassReg>140) && VtypeWithTau==4";  // FIXME
//const TCut    g_cutalldata  = "(HmassReg<110 || HmassReg>140) && VtypeWithTau==4 && (!(207883<=EVENT.run && EVENT.run<=208307))";  // FIXME
//const TCut    g_cutallmc    = "(BDTtt_125[0]>-0.5 && BDTvjlf_125[0]>-0.5 && BDTzz_125[0]>-0.3) && VtypeWithTau==4";  // FIXME
//const TCut    g_cutalldata  = "(BDTtt_125[0]>-0.5 && BDTvjlf_125[0]>-0.5 && BDTzz_125[0]>-0.3) && VtypeWithTau==4 && (!(207883<=EVENT.run && EVENT.run<=208307))";  // FIXME
const double  g_xlow        = -1.;
const double  g_xup         = 1.;
const bool    g_manyplots   = false;
const TString g_wsname      = "vhbb_Znn_J12_8TeV.root";
//const TString g_wsname      = "vhbb_Znn_J12_$CHANNEL_8TeV.root";
const TString g_dcname      = "vhbb_Znn_J12_$CHANNEL_8TeV.txt";
const TString g_rootname    = "vhbb_Znn_J12_$CHANNEL_TH1.root";
const TString g_plotdir     = "plotsJ12/";

#ifndef VVANALYSIS
const TString g_var         = Form("BDTregular_%i[ISYST]", massH);  // ISYST is replaced in the main function.
//const TString g_var         = Form("BDTregular_fj_%i[ISYST]", massH);  // ISYST is replaced in the main function.
const TString g_varslice    = Form("slice(BDTregular_%i[ISYST], (BDTtt_%i[ISYST]>-0.5), (BDTvjlf_%i[ISYST]>-0.5), (BDTzz_%i[ISYST]>-0.3))", massH, massH, massH, massH);  // BDTzz > -0.35?
//const TString g_varslice    = Form("slice(BDTregular_fj_%i[ISYST], (BDTtt_%i[ISYST]>-0.6), (BDTvjlf_%i[ISYST]>-0.5), (BDTzz_%i[ISYST]>-0.3))", massH, massH, massH, massH);
//const TString g_varslice    = Form("hybridslice(BDTregularhybrid_%i[ISYST], BDTregularhybrid_fj_%i[ISYST], (nfathFilterJets>0 && FatH.FatHiggsFlag==1))", massH, massH);
//const TString g_varslice    = Form("slice(BDTttbar_%i[ISYST], BDTvjl_%i[ISYST], BDTzz_%i[ISYST], BDTregular_%i[ISYST], -0.5, -0.4, -0.4)", massH, massH, massH, massH);
#else
const TString g_var         = Form("BDTregular__sigzzhf_%i[ISYST]", massH);  // ISYST is replaced in the main function.
const TString g_varslice    = Form("slice(BDTregular__sigzzhf_%i[ISYST], (BDTtt__sigzzhf_%i[ISYST]>-0.4), (BDTvjlf__sigzzhf_%i[ISYST]>-0.35))", massH, massH, massH);
#endif


/// The channel
TString channel = "ZnunuHighPt";

/// Systematics
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
const double g_scalefactors_ZnunuHighPt[7] = {
  //1.154, 2.374, 0.736, 1.131, 2.198, 1.136, 0.954
  //0.825, 1.994, 0.717, 1.165, 2.158, 1.101, 0.933
    0.930, 2.124, 0.700, 1.175, 2.135, 1.125, 0.991
};
const double g_scalefactors_lnN_ZnunuHighPt[7] = {
  //1.03, 1.30, 1.40, 1.03, 1.07, 1.10, 1.03
  //1.03, 1.30, 1.38, 1.04, 1.09, 1.12, 1.04
    1.03, 1.15, 1.15, 1.03, 1.09, 1.11, 1.03
};

const double g_scalefactors_ZnunuMedPt[7] = {
  //1.002, 2.175, 1.132, 1.239, 2.502, 1.114, 0.959
  //0.868, 1.923, 1.022, 1.306, 2.523, 1.102, 0.944
    0.910, 2.080, 0.750, 1.246, 2.300, 1.113, 1.015
};
const double g_scalefactors_lnN_ZnunuMedPt[7] = {
  //1.03, 1.28, 1.29, 1.03, 1.07, 1.10, 1.03
  //1.03, 1.34, 1.35, 1.04, 1.09, 1.15, 1.04
    1.03, 1.15, 1.15, 1.03, 1.10, 1.12, 1.03
};

const double g_scalefactors_ZnunuLowPt[7] = {
  //0.865, 2.588, 0.783, 1.259, 2.190, 1.422, 0.939
  //0.726, 2.797, 0.946, 1.270, 2.159, 1.426, 1.000
    0.850, 2.224, 0.683, 1.203, 2.480, 1.420, 1.050
};
const double g_scalefactors_lnN_ZnunuLowPt[7] = {
  //1.05, 1.34, 1.39, 1.03, 1.07, 1.09, 1.04
  //1.03, 1.31, 1.40, 1.04, 1.09, 1.12, 1.05
    1.04, 1.15, 1.15, 1.04, 1.12, 1.12, 1.04
};

const double g_scalefactors_ZnunuLowCSV[7] = {
    1.000, 2.000, 1.000, 1.000, 2.000, 1.000, 1.000
};
const double g_scalefactors_lnN_ZnunuLowCSV[7] = {
    1.15, 1.60, 1.50, 1.12, 1.30, 1.25, 1.10
};

const double g_scalefactors_lnN_prefit[7] = {
    1.15, 1.60, 1.50, 1.12, 1.30, 1.25, 1.10
};

//const double g_scalefactors_QCD[3] = {0.02, 0.2, 0.5};
const double g_scalefactors_QCD[3] = {0.004, 0.04, 0.06};
////////////////////////////////////////////////////////////////////////////////
/// END Configuration                                                        ///
////////////////////////////////////////////////////////////////////////////////

///_____________________________________________________________________________
/// Classes

class EventsJ12 {  // jmax = 12 in datacard, with 7 scale factors.
public:
    EventsJ12()
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

    ~EventsJ12()
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
#ifdef QCDSHAPE
        assert(QCD       != 0 && QCD      ->GetEntriesFast() > 0);
#endif
        assert(Wj0b_syst != 0 && Wj0b_syst->GetEntriesFast() > 0);
        assert(Wj1b_syst != 0 && Wj1b_syst->GetEntriesFast() > 0);
        assert(Wj2b_syst != 0 && Wj2b_syst->GetEntriesFast() > 0);
        assert(Zj0b_syst != 0 && Zj0b_syst->GetEntriesFast() > 0);
        assert(Zj1b_syst != 0 && Zj1b_syst->GetEntriesFast() > 0);
        assert(Zj2b_syst != 0 && Zj2b_syst->GetEntriesFast() > 0);
        assert(TT_syst   != 0 && TT_syst  ->GetEntriesFast() > 0);
        assert(ZH_SM     != 0 && ZH_SM    ->GetEntriesFast() > 0);
        assert(WH_SM     != 0 && WH_SM    ->GetEntriesFast() > 0);
        assert(data_obs  != 0 && data_obs ->GetEntriesFast() > 0);
        return;
    }
    
    void set_scalefactors(const double sf[7]) {
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
};  // EventsJ12

///_____________________________________________________________________________
/// Functions

float slice(const float var1, int cut1)
{
    const float xlow=-1., xup=1.;
    // Assume original range is [xlow,xup], check if this is always true!
    if(!(xlow <= var1 && var1 <= xup))  std::cerr << "ERROR: var1 = " << var1 << std::endl;
    const float scale = (xup - xlow) / 2.0 / 2.0;
    if (!cut1) {
        return ((var1+1.0)*scale + xlow);
    } else { // pass all 1 cuts
        return ((var1+3.0)*scale + xlow);
    }
}

float slice(const float var1, int cut1, int cut2)
{
    const float xlow=-1., xup=1.;
    // Assume original range is [xlow,xup], check if this is always true!
    if(!(xlow <= var1 && var1 <= xup))  std::cerr << "ERROR: var1 = " << var1 << std::endl;
    const float scale = (xup - xlow) / 2.0 / 3.0;
    if (!cut1) {
        return ((var1+1.0)*scale + xlow);
    } else if (!cut2) {
        return ((var1+3.0)*scale + xlow);
    } else { // pass all 2 cuts
        return ((var1+5.0)*scale + xlow);
    }
}

float slice(const float var1, int cut1, int cut2, int cut3)
{
    const float xlow=-1., xup=1.;
    // Assume original range is [xlow,xup], check if this is always true!
    if(!(xlow <= var1 && var1 <= xup))  std::cerr << "ERROR: var1 = " << var1 << std::endl;
    const float scale = (xup - xlow) / 2.0 / 4.0;
    if (!cut1) {
        return ((var1+1.0)*scale + xlow);
    } else if (!cut2) {
        return ((var1+3.0)*scale + xlow);
    } else if (!cut3) {
        return ((var1+5.0)*scale + xlow);
    } else { // pass all 3 cuts
        return ((var1+7.0)*scale + xlow);
    }
}

float hybridslice(const float var1, const float var2, int cut1)
{
    const float xlow=-1., xup=1.;
    // Assume original range is [xlow,xup], check if this is always true!
    if(!(xlow <= var1 && var1 <= xup))  std::cerr << "ERROR: var1 = " << var1 << std::endl;
    const float scale = (xup - xlow) / 2.0 / 2.0;
    if (!cut1) {
        return ((var1+1.0)*scale + xlow);
    } else { // pass all 1 cuts
        return ((var2+3.0)*scale + xlow);
    }
}

///_____________________________________________________________________________
/// The core function that makes plots, prints stats, and writes the datacard.

using namespace std;

/// Declare rebinner
Rebinner* rebinner = 0;

void MakePlots(const EventsJ12 * ev, TString var_, 
               TCut cutmc_, TCut cutdata_, const TString syst, const TString region, 
               const char* xtitle, int nbins, double xlow, double xup, 
               long long newnbins, double errorffirst, double errorflast, 
               const double scalefactors_lnN[7], 
               TString options="printStat:printCard:plotData:plotLog:plotSig")
{
#ifdef HCPANALYSIS
    cutmc_   *= TCut("12350/19040 * PUweightAB/PUweight");
    cutdata_ *= TCut("EVENT.run<=203002");
#endif

#ifdef MJJANALYSIS
    TPRegexp regexp("selectFlags\\[0\\]\\[([0-9]+)\\]");
    TString str = "";
    str       = cutmc_.GetTitle();
    regexp.Substitute(str, "selectFlags[0][$1] * selectMjj[$1]", "g");  // "g" for global
    cutmc_    = TCut(str);
    cutmc_   *= "0.5";
    str       = cutdata_.GetTitle();
    regexp.Substitute(str, "selectFlags[0][$1] * selectMjj[$1]", "g");  // "g" for global
    cutdata_  = TCut(str);
    
    if (TString(xtitle) == "BDT") {
        var_      = "HmassReg+0";
        if (syst == "CMS_vhbb_Znn_res_jUp")
            var_  = "HmassReg_res_j_up";
        else if (syst == "CMS_vhbb_Znn_res_jDown")
            var_  = "HmassReg_res_j_down";
        else if (syst == "CMS_vhbb_Znn_scale_jUp")
            var_  = "HmassReg_scale_j_up";
        else if (syst == "CMS_vhbb_Znn_scale_jDown")
            var_  = "HmassReg_scale_j_down";
        
        // FIXME post mortem fix
        if (channel == "ZnunuHighPt") {
            if (syst == "CMS_vhbb_Znn_res_jUp") {
                cutmc_   *= "HptReg_res_j_up>190";
                cutdata_ *= "HptReg_res_j_up>190";
            } else if (syst == "CMS_vhbb_Znn_res_jDown") {
                cutmc_   *= "HptReg_res_j_down>190";
                cutdata_ *= "HptReg_res_j_down>190";
            } else if (syst == "CMS_vhbb_Znn_scale_jUp") {
                cutmc_   *= "HptReg_scale_j_up>190";
                cutdata_ *= "HptReg_scale_j_up>190";
            } else if (syst == "CMS_vhbb_Znn_scale_jDown") {
                cutmc_   *= "HptReg_scale_j_down>190";
                cutdata_ *= "HptReg_scale_j_down>190";
            } else {
                cutmc_   *= "HptReg>190";
                cutdata_ *= "HptReg>190";
            }
        } else if (channel == "ZnunuMedPt") {
            //if (syst == "CMS_vhbb_Znn_res_jUp") {
            //    cutmc_   *= "HptReg_res_j_up>140 && hJet_ptReg_res_j_up[0]>80";
            //    cutdata_ *= "HptReg_res_j_up>140 && hJet_ptReg_res_j_up[0]>80";
            //} else if (syst == "CMS_vhbb_Znn_res_jDown") {
            //    cutmc_   *= "HptReg_res_j_down>140 && hJet_ptReg_res_j_down[0]>80";
            //    cutdata_ *= "HptReg_res_j_down>140 && hJet_ptReg_res_j_down[0]>80";
            //} else if (syst == "CMS_vhbb_Znn_scale_jUp") {
            //    cutmc_   *= "HptReg_scale_j_up>140 && hJet_ptReg_scale_j_up[0]>80";
            //    cutdata_ *= "HptReg_scale_j_up>140 && hJet_ptReg_scale_j_up[0]>80";
            //} else if (syst == "CMS_vhbb_Znn_scale_jDown") {
            //    cutmc_   *= "HptReg_scale_j_down>140 && hJet_ptReg_scale_j_down[0]>80";
            //    cutdata_ *= "HptReg_scale_j_down>140 && hJet_ptReg_scale_j_down[0]>80";
            //} else {
            //    cutmc_   *= "HptReg>140 && hJet_ptReg[0]>80";
            //    cutdata_ *= "HptReg>140 && hJet_ptReg[0]>80";
            //}
            
            if (syst == "CMS_vhbb_Znn_res_jUp") {
                cutmc_   *= "HptReg_res_j_up > 140";
                cutdata_ *= "HptReg_res_j_up > 140";
            } else if (syst == "CMS_vhbb_Znn_res_jDown") {
                cutmc_   *= "HptReg_res_j_down > 140";
                cutdata_ *= "HptReg_res_j_down > 140";
            } else if (syst == "CMS_vhbb_Znn_scale_jUp") {
                cutmc_   *= "HptReg_scale_j_up > 140";
                cutdata_ *= "HptReg_scale_j_up > 140";
            } else if (syst == "CMS_vhbb_Znn_scale_jDown") {
                cutmc_   *= "HptReg_scale_j_down > 140";
                cutdata_ *= "HptReg_scale_j_down > 140";
            } else {
                cutmc_   *= "HptReg > 140";
                cutdata_ *= "HptReg > 140";
            }

        }
        
        xtitle    = "M_{b#bar{b}} [GeV]";
        nbins     = 17;
        xlow      = 0.;
        xup       = 255.;
        newnbins  = 17;
        
        if (!options.Contains("!plotLog"))
            options.ReplaceAll("plotLog", "!plotLog");
    }
#endif

    if (channel == "ZnunuLowPt") {  // FIXME
        cutmc_   *= TCut("nPVs<=25");
        cutdata_ *= TCut("nPVs<=25");
    }

    const TCut    cutmc   = cutmc_;
    const TCut    cutdata = cutdata_;
    const TString var     = var_;
    
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
    TH1F * hZH_SM_0     = new TH1F("ZH_SM_0"    , "", nbins, xlow, xup);
    TH1F * hWH_SM_0     = new TH1F("WH_SM_0"    , "", nbins, xlow, xup);
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
    histos_0.push_back(hZH_SM_0);
    histos_0.push_back(hWH_SM_0);
    histos_0.push_back(hVH_0);
    histos_0.push_back(hVV_0);
    histos_0.push_back(hmc_exp_0);
    histos_0.push_back(hdata_obs_0);
    
    for (UInt_t ih = 0; ih < histos_0.size(); ih++)
        histos_0.at(ih)->Sumw2();

    TCut cutvhewk = "weightSignalEWKNew";
    TCut cutvhqcd = "weightSignalQCD";
#if !defined(VHEWKCORRECTION)
    ev->ZH->Project("ZH_0", var, cutmc);
    std::clog << "... DONE: project ZH_0." << std::endl;
    ev->WH->Project("WH_0", var, cutmc);
    std::clog << "... DONE: project WH_0." << std::endl;
#elif !defined(VHQCDCORRECTION)
    ev->ZH->Project("ZH_0", var, cutmc * cutvhewk);
    std::clog << "... DONE: project ZH_0." << std::endl;
    ev->WH->Project("WH_0", var, cutmc * cutvhewk);
    std::clog << "... DONE: project WH_0." << std::endl;
#else
    ev->ZH->Project("ZH_0", var, cutmc * cutvhewk * cutvhqcd);
    std::clog << "... DONE: project ZH_0." << std::endl;
    ev->WH->Project("WH_0", var, cutmc * cutvhewk * cutvhqcd);
    std::clog << "... DONE: project WH_0." << std::endl;
#endif


    // Apply slope
    double WJSlopeErr = 0.0020;
    double ZJSlopeErr = 0.0025;
    double WJSlope = -0.50;
    double ZJSlope = -1.00;
    TCut cutwjslope = Form("apply_pt_slope(genW.pt, %f, 150)", WJSlope * WJSlopeErr);
    TCut cutzjslope = Form("apply_pt_slope(genZ.pt, %f, 130)", ZJSlope * ZJSlopeErr);
    
    // Apply scale factors
    TCut cutsf = "1.0";
    cutsf = Form("%f", ev->sf_Wj0b);
    ev->Wj0b->Project("Wj0b_0", var, cutmc * cutsf * cutwjslope);
    std::clog << "... DONE: project Wj0b_0, scaled by " << cutsf * cutwjslope << "." << std::endl;
    cutsf = Form("%f", ev->sf_Wj1b);
    ev->Wj1b->Project("Wj1b_0", var, cutmc * cutsf * cutwjslope);
    std::clog << "... DONE: project Wj1b_0, scaled by " << cutsf * cutwjslope << "." << std::endl;
    cutsf = Form("%f", ev->sf_Wj2b);
    ev->Wj2b->Project("Wj2b_0", var, cutmc * cutsf * cutwjslope);
    std::clog << "... DONE: project Wj2b_0, scaled by " << cutsf * cutwjslope << "." << std::endl;
    cutsf = Form("%f", ev->sf_Zj0b);
    ev->Zj0b->Project("Zj0b_0", var, cutmc * cutsf * cutzjslope);
    std::clog << "... DONE: project Zj0b_0, scaled by " << cutsf * cutzjslope << "." << std::endl;
    cutsf = Form("%f", ev->sf_Zj1b);
    ev->Zj1b->Project("Zj1b_0", var, cutmc * cutsf * cutzjslope);
    std::clog << "... DONE: project Zj1b_0, scaled by " << cutsf * cutzjslope << "." << std::endl;
    cutsf = Form("%f", ev->sf_Zj2b);
    ev->Zj2b->Project("Zj2b_0", var, cutmc * cutsf * cutzjslope);
    std::clog << "... DONE: project Zj2b_0, scaled by " << cutsf * cutzjslope << "." << std::endl;
    cutsf = Form("%f", ev->sf_TT);
    ev->TT->Project("TT_0", var, cutmc * cutsf);
    std::clog << "... DONE: project TT_0, scaled by " << cutsf << "." << std::endl;
    
    ev->s_Top->Project("s_Top_0", var, cutmc);
    std::clog << "... DONE: project s_Top_0." << std::endl;
    ev->VVLF->Project("VVLF_0", var, cutmc);
    std::clog << "... DONE: project VVLF_0." << std::endl;
    ev->VVHF->Project("VVHF_0", var, cutmc);
    std::clog << "... DONE: project VVHF_0." << std::endl;

#ifdef QCDSHAPE
    // Apply QCD scale factor
    if (channel == "ZnunuHighPt")
        cutsf = Form("%f", g_scalefactors_QCD[0]);
    else if (channel == "ZnunuMedPt")
        cutsf = Form("%f", g_scalefactors_QCD[1]);
    else if (channel == "ZnunuLowPt")
        cutsf = Form("%f", g_scalefactors_QCD[2]);
    ev->QCD->Project("QCD_0", var, cutmc * cutsf);
    std::clog << "... DONE: project QCD_0, scaled by " << cutsf << "." << std::endl;
    //ev->QCD->Project("QCD_0", var, cutmc);
    //std::clog << "... DONE: project QCD_0." << std::endl;
#endif

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

    ev->ZH_SM->Project("ZH_SM_0", var, cutmc);
    std::clog << "... DONE: project ZH_SM_0." << std::endl;
    ev->WH_SM->Project("WH_SM_0", var, cutmc);
    std::clog << "... DONE: project WH_SM_0." << std::endl;
    
    ev->data_obs->Project("data_obs_0", var, cutdata);
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
#ifndef VVANALYSIS
    hmc_exp_0->Add(hVVLF_0);
    hmc_exp_0->Add(hVVHF_0);
#else
    hmc_exp_0->Add(hVVLF_0);
    //hmc_exp_0->Add(hVH_0);  // VH is not counted as background
#endif
    //hmc_exp_0->Add(hQCD_0);  // QCD is added after rebinning
    std::clog << "... DONE: add MC backgrounds to mc_exp_0." << std::endl;

    if (rebinner == 0) {
        rebinner = new Rebinner(newnbins, errorffirst, errorflast, xlow, xup);
        rebinner->set_signal_backgr(hVH_0, hmc_exp_0);
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
    TH1F * hZH_SM     = rebinner->rebin(hZH_SM_0      , newnbins, "ZH_SM"     );
    TH1F * hWH_SM     = rebinner->rebin(hWH_SM_0      , newnbins, "WH_SM"     );
    TH1F * hVH        = rebinner->rebin(hVH_0         , newnbins, "VH"        );
    TH1F * hVV        = rebinner->rebin(hVV_0         , newnbins, "VV"        );
    TH1F * hmc_exp    = rebinner->rebin(hmc_exp_0     , newnbins, "mc_exp"    );
    TH1F * hdata_obs  = rebinner->rebin(hdata_obs_0   , newnbins, "data_obs"  );
    
    // Add QCD after rebinning
    hmc_exp->Add(hQCD);

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
    histos.push_back(hZH_SM);
    histos.push_back(hWH_SM);
    histos.push_back(hVH);
    histos.push_back(hVV);
    histos.push_back(hmc_exp);
    histos.push_back(hdata_obs);
    
    assert(histos.size() == histos_0.size());
    for (UInt_t ih = 0; ih < histos.size(); ih++)
        histos.at(ih)->Sumw2();
    
    if (printStat && syst == "NONE") {
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
        
        p_integral = hWj0b->IntegralAndError(p_bin1, p_bin2, p_error);
        B += p_integral;
        std::cout << setw(10) << left << "Wj0b " << "(" << p_bin1 << ", " << p_bin2 << ") = " << setw(8) << right << p_integral << " +/- " << setw(6) << right << p_error << std::endl;
        p_integral = hWj1b->IntegralAndError(p_bin1, p_bin2, p_error);
        B += p_integral;
        std::cout << setw(10) << left << "Wj1b " << "(" << p_bin1 << ", " << p_bin2 << ") = " << setw(8) << right << p_integral << " +/- " << setw(6) << right << p_error << std::endl;
        p_integral = hWj2b->IntegralAndError(p_bin1, p_bin2, p_error);
        B += p_integral;
        std::cout << setw(10) << left << "Wj2b " << "(" << p_bin1 << ", " << p_bin2 << ") = " << setw(8) << right << p_integral << " +/- " << setw(6) << right << p_error << std::endl;
        p_integral = hZj0b->IntegralAndError(p_bin1, p_bin2, p_error);
        B += p_integral;
        std::cout << setw(10) << left << "Zj0b " << "(" << p_bin1 << ", " << p_bin2 << ") = " << setw(8) << right << p_integral << " +/- " << setw(6) << right << p_error << std::endl;
        p_integral = hZj1b->IntegralAndError(p_bin1, p_bin2, p_error);
        B += p_integral;
        std::cout << setw(10) << left << "Zj1b " << "(" << p_bin1 << ", " << p_bin2 << ") = " << setw(8) << right << p_integral << " +/- " << setw(6) << right << p_error << std::endl;
        p_integral = hZj2b->IntegralAndError(p_bin1, p_bin2, p_error);
        B += p_integral;
        std::cout << setw(10) << left << "Zj2b " << "(" << p_bin1 << ", " << p_bin2 << ") = " << setw(8) << right << p_integral << " +/- " << setw(6) << right << p_error << std::endl;
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
        p_integral = hVV->IntegralAndError(p_bin1, p_bin2, p_error);
        std::cout << setw(10) << left << "VV " << "(" << p_bin1 << ", " << p_bin2 << ") = " << setw(8) << right << p_integral << " +/- " << setw(6) << right << p_error << std::endl;
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
        TString channel_8TeV_1 = channel_8TeV;  // scale factors
        if (channel == "ZnunuLowCSV")
            channel_8TeV_1.ReplaceAll("ZnunuLowCSV", "ZnunuHighPt");
        
        std::ofstream dc;
        dc.setf(ios::fixed,ios::floatfield);
        dc.precision(3);
        dc.open(dcname);

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
        dc << "process     ZH         WH         Wj0b       Wj1b       Wj2b       Zj0b       Zj1b       Zj2b       TT         s_Top      VVLF       VVHF       QCD        " << std::endl;
#ifndef VVANALYSIS
        dc << "process     -1         0          1          2          3          4          5          6          7          8          9          10         11         " << std::endl;
#else
        dc << "process     11         12         1          2          3          4          5          6          7          8          9          0          10         " << std::endl;
#endif
        dc << "rate        " << setw(10) << right << hZH->Integral() << " " << setw(10) << right << hWH->Integral() << " " << setw(10) << right << hWj0b->Integral() << " " << setw(10) << right << hWj1b->Integral() << " " << setw(10) << right << hWj2b->Integral() << " " << setw(10) << right << hZj0b->Integral() << " " << setw(10) << right << hZj1b->Integral() << " " << setw(10) << right << hZj2b->Integral() << " " << setw(10) << right << hTT->Integral() << " " << setw(10) << right << hs_Top->Integral() << " " << setw(10) << right << hVVLF->Integral() << " " << setw(10) << right << hVVHF->Integral() << " " << setw(10) << right << hQCD->Integral() << std::endl;
        dc.precision(2);
        dc << "-----------------------------------" << std::endl;
        dc << "" << std::endl;
        dc << "### Flat ########################### #####  ZH    WH    Wj0b  Wj1b  Wj2b  Zj0b  Zj1b  Zj2b  TT    s_Top VVLF  VVHF  QCD  " << std::endl;
        dc << "lumi_8TeV                            lnN    1.05  1.05  -     -     -     -     -     -     -     1.05  1.05  1.05  1.05 " << std::endl;
        dc << "pdf_qqbar                            lnN    1.01  1.01  -     -     -     -     -     -     -     -     1.01  1.01  -    " << std::endl;
        dc << "pdf_gg                               lnN    -     -     -     -     -     -     -     -     -     1.01  -     -     1.01 " << std::endl;
        dc << "QCDscale_VH                          lnN    1.04  1.04  -     -     -     -     -     -     -     -     -     -     -    " << std::endl;
        dc << "QCDscale_ttbar                       lnN    -     -     -     -     -     -     -     -     -     1.06  -     -     -    " << std::endl;
        dc << "QCDscale_VV                          lnN    -     -     -     -     -     -     -     -     -     -     1.04  1.04  -    " << std::endl;
        dc << "QCDscale_QCD                         lnN    -     -     -     -     -     -     -     -     -     -     -     -     1.30 " << std::endl;
#ifndef VHEWKCORRECTION
        dc << "CMS_vhbb_boost_EWK_8TeV              lnN    1.05  1.10  -     -     -     -     -     -     -     -     -     -     -    " << std::endl;
#else
        dc << "CMS_vhbb_boost_EWK_8TeV              lnN    1.02  1.02  -     -     -     -     -     -     -     -     -     -     -    " << std::endl;
#endif
#ifndef VHQCDCORRECTION
        dc << "CMS_vhbb_boost_QCD_8TeV              lnN    1.10  1.10  -     -     -     -     -     -     -     -     -     -     -    " << std::endl;
#else
        dc << "CMS_vhbb_boost_QCD_8TeV              lnN    1.05  1.05  -     -     -     -     -     -     -     -     -     -     -    " << std::endl;
#endif
        dc << "CMS_vhbb_ST                          lnN    -     -     -     -     -     -     -     -     -     1.15  -     -     -    " << std::endl;
#ifndef VVANALYSIS
        dc << "CMS_vhbb_VV                          lnN    -     -     -     -     -     -     -     -     -     -     1.25  1.25  -    " << std::endl;
#else
        dc << "CMS_vhbb_VV                          lnN    -     -     -     -     -     -     -     -     -     -     1.15  1.15  -    " << std::endl;
        dc << "CMS_vhbb_VH                          lnN    1.50  1.50  -     -     -     -     -     -     -     -     -     -     -    " << std::endl;
#endif
        dc << "#CMS_vhbb_MET_nojets                 lnN    1.03  1.03  -     -     -     -     -     -     -     1.03  1.03  1.03  1.03 " << std::endl;
        dc << "#CMS_vhbb_trigger_MET                lnN    1.03  1.03  -     -     -     -     -     -     -     1.03  1.03  1.03  1.03 " << std::endl;
        dc << "CMS_vhbb_Wj0b_SF_" << channel_8TeV_1 << "    lnN    -     -     " << scalefactors_lnN[0] << "  -     -     -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_Wj1b_SF_" << channel_8TeV_1 << "    lnN    -     -     -     " << scalefactors_lnN[1] << "  -     -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_Wj2b_SF_" << channel_8TeV_1 << "    lnN    -     -     -     -     " << scalefactors_lnN[2] << "  -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_Zj0b_SF_" << channel_8TeV_1 << "    lnN    -     -     -     -     -     " << scalefactors_lnN[3] << "  -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_Zj1b_SF_" << channel_8TeV_1 << "    lnN    -     -     -     -     -     -     " << scalefactors_lnN[4] << "  -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_Zj2b_SF_" << channel_8TeV_1 << "    lnN    -     -     -     -     -     -     -     " << scalefactors_lnN[5] << "  -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_TT_SF_" << channel_8TeV_1 << "      lnN    -     -     -     -     -     -     -     -     " << scalefactors_lnN[6] << "  -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_QCD_SF_" << channel_8TeV_1 << "     lnN    -     -     -     -     -     -     -     -     -     -     -     -     1.60 " << std::endl;
        dc << "### Shape ########################## #####  ZH    WH    Wj0b  Wj1b  Wj2b  Zj0b  Zj1b  Zj2b  TT    s_Top VVLF  VVHF  QCD  " << std::endl;
        //dc << "CMS_vhbb_boost_QCD_8TeV                   shape  1.00  1.00  -     -     -     -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_eff_b                       shape  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  -    " << std::endl;
        dc << "CMS_vhbb_fake_b_8TeV                 shape  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  -    " << std::endl;
        dc << "CMS_vhbb_Znn_res_j                   shape  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  -    " << std::endl;
        dc << "CMS_vhbb_Znn_scale_j                 shape  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  -    " << std::endl;
        dc << "#UEPS                                shape  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  -    " << std::endl;
        //dc << "CMS_vhbb_trigger_MET_Znn_8TeV        shape  1.00  1.00  -     -     -     -     -     -     -     1.00  1.00  1.00  -    " << std::endl;
        //dc << "CMS_vhbb_trigger_CSV_Znn_8TeV        shape  1.00  1.00  -     -     -     -     -     -     -     1.00  -     1.00  -    " << std::endl;
        //dc << "CMS_vhbb_trigger_CSV_fake_Znn_8TeV   shape  -     -     -     -     -     -     -     -     -     -     1.00  -     -    " << std::endl;
        dc << "CMS_vhbb_trigger_MET_Znn_8TeV        shape  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  -    " << std::endl;
        dc << "CMS_vhbb_trigger_CSV_Znn_8TeV        shape  1.00  1.00  -     1.00  1.00  -     1.00  1.00  1.00  1.00  -     1.00  -    " << std::endl;
        dc << "CMS_vhbb_trigger_CSV_fake_Znn_8TeV   shape  -     -     1.00  -     -     1.00  -     -     -     -     1.00  -     -    " << std::endl;
        if (channel != "ZnunuLowPt") {
        //dc << "CMS_vhbb_WJModel_Znn_8TeV            shape  -     -     1.00  1.00  1.00  -     -     -     -     -     -     -     -    " << std::endl;
        //dc << "CMS_vhbb_ZJModel_Znn_8TeV            shape  -     -     -     -     -     1.00  1.00  1.00  -     -     -     -     -    " << std::endl;
        dc << "#CMS_vhbb_Wj0bModel_Znn_8TeV          shape  -     -     1.00  -     -     -     -     -     -     -     -     -     -    " << std::endl;
        dc << "#CMS_vhbb_Wj1bModel_Znn_8TeV          shape  -     -     -     1.00  -     -     -     -     -     -     -     -     -    " << std::endl;
        dc << "#CMS_vhbb_Wj2bModel_Znn_8TeV          shape  -     -     -     -     1.00  -     -     -     -     -     -     -     -    " << std::endl;
        dc << "#CMS_vhbb_Zj0bModel_Znn_8TeV          shape  -     -     -     -     -     1.00  -     -     -     -     -     -     -    " << std::endl;
        dc << "#CMS_vhbb_Zj1bModel_Znn_8TeV          shape  -     -     -     -     -     -     1.00  -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_Zj2bModel_Znn_8TeV          shape  -     -     -     -     -     -     -     1.00  -     -     -     -     -    " << std::endl;
        }
        dc << "CMS_vhbb_TTModel_Znn_8TeV            shape  -     -     -     -     -     -     -     -     1.00  -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_WJSlope_Znn_8TeV            shape  -     -     1.00  1.00  1.00  -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_ZJSlope_Znn_8TeV            shape  -     -     -     -     -     1.00  1.00  1.00  -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_statZH_" << channel_8TeV << "     shape  1.00  -     -     -     -     -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_statWH_" << channel_8TeV << "     shape  -     1.00  -     -     -     -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_statWj0b_" << channel_8TeV << "   shape  -     -     1.00  -     -     -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_statWj1b_" << channel_8TeV << "   shape  -     -     -     1.00  -     -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_statWj2b_" << channel_8TeV << "   shape  -     -     -     -     1.00  -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_statZj0b_" << channel_8TeV << "   shape  -     -     -     -     -     1.00  -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_statZj1b_" << channel_8TeV << "   shape  -     -     -     -     -     -     1.00  -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_statZj2b_" << channel_8TeV << "   shape  -     -     -     -     -     -     -     1.00  -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_statTT_" << channel_8TeV << "     shape  -     -     -     -     -     -     -     -     1.00  -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_stats_Top_" << channel_8TeV << "  shape  -     -     -     -     -     -     -     -     -     1.00  -     -     -    " << std::endl;
        dc << "CMS_vhbb_statVVLF_" << channel_8TeV << "   shape  -     -     -     -     -     -     -     -     -     -     1.00  -     -    " << std::endl;
        dc << "CMS_vhbb_statVVHF_" << channel_8TeV << "   shape  -     -     -     -     -     -     -     -     -     -     -     1.00  -    " << std::endl;
        dc << "CMS_vhbb_statQCD_" << channel_8TeV << "    shape  -     -     -     -     -     -     -     -     -     -     -     -     1.00 " << std::endl;
        dc << "#################################### #####  ZH    WH    Wj0b  Wj1b  Wj2b  Zj0b  Zj1b  Zj2b  TT    s_Top VVLF  VVHF  QCD  " << std::endl;
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
        setHisto(hWj0b, "WjLF");
        setHisto(hWj1b, "WjHFc");
        setHisto(hWj2b, "WjHFb");
        setHisto(hZj0b, "ZjLF");
        setHisto(hZj1b, "ZjHFc");
        setHisto(hZj2b, "ZjHFb");
        setHisto(hTT, "TT");
        setHisto(hs_Top, "s_Top");
        setHisto(hVV, "VV");
        setHisto(hVVLF, "VV");
        setHisto(hVVHF, "VVHF");
        setHisto(hQCD, "QCD");
        
        TH1F * hdata_test = (TH1F *) hdata_obs->Clone("hdata_test");  // blinded plot
        TH1F * hmc_test = (TH1F *) hmc_exp->Clone("hmc_test");  // for chi2 and KS test
        hdata_test->Sumw2();
        hmc_test->Sumw2();
        Int_t nbins_plot = hdata_test->GetNbinsX();
        assert(nbins_plot == hmc_test->GetNbinsX());
        if (!plotData) {  // be blind to the most sensitive bins
#ifndef MJJANALYSIS
            for (Int_t i = TMath::Max(Int_t(nbins_plot*0.75), nbins_plot-5); i < nbins_plot+2; i++) {
                hdata_test->SetBinContent(i, 0.);
                hdata_test->SetBinError(i, 0.);
                hmc_test->SetBinContent(i, 0.);
                hmc_test->SetBinError(i, 0.);
            }
#else
            for (Int_t i = hdata_test->FindFixBin(105+1); i < hdata_test->FindFixBin(150-1)+1; i++) {
            //for (Int_t i = hdata_test->FindFixBin(105+1)-2; i < hdata_test->FindFixBin(150-1)+1; i++) {  // hide VV as well
                hdata_test->SetBinContent(i, 0.);
                hdata_test->SetBinError(i, 0.);
                hmc_test->SetBinContent(i, 0.);
                hmc_test->SetBinError(i, 0.);
            }
#endif
        }
        setHisto(hdata_test, "data_obs");

        std::clog << "MakePlots(): Setting up the stack..." << std::endl;
        THStack * hs = new THStack("hs", "");
#if !defined(VVANALYSIS) && !defined(MJJANALYSIS)
        hs->Add(hVVHF);
        hs->Add(hVVLF);
#elif !defined(MJJANALYSIS)
        hs->Add(hVVLF);
#endif
        hs->Add(hQCD);
        hs->Add(hs_Top);
        hs->Add(hTT);
        hs->Add(hWj0b);
        hs->Add(hWj1b);
        hs->Add(hWj2b);
        hs->Add(hZj0b);
        hs->Add(hZj1b);
        hs->Add(hZj2b);
#if !defined(VVANALYSIS) && !defined(MJJANALYSIS)
        if (plotSig)  hs->Add(hVH);
#elif !defined(MJJANALYSIS)
        hs->Add(hVVHF);
        if (plotSig)  hs->Add(hVH);
#endif
#if defined(MJJANALYSIS)
        hs->Add(hVVLF);
        hs->Add(hVVHF);
        if (plotSig)  hs->Add(hVH);
#endif
        
        double ymax = TMath::Max(hdata_test->GetMaximum(), hs->GetMaximum());
        hs->SetMaximum(ymax * 1.7 + (ymax>1 ? sqrt(ymax) : 0.));
        if (plotLog)
            hs->SetMaximum(ymax * 200 + (ymax>1 ? sqrt(ymax) : 0.));
        hs->SetMinimum(0.01);
        
        // Setup auxiliary histograms
        std::clog << "MakePlots(): Setting up auxiliary histograms..." << std::endl;
        TH1F * staterr = (TH1F *) hmc_exp->Clone("staterr");
        staterr->Sumw2();
        staterr->SetFillColor(kRed);
        staterr->SetMarkerSize(0);
        staterr->SetFillStyle(3013);
        
        TH1F * ratio = (TH1F *) hdata_test->Clone("ratio");
        ratio->Sumw2();
        ratio->SetMarkerSize(0.8);
        //ratio->SetMarkerSize(0.5);
        ratio->Divide(hdata_test, hmc_exp, 1., 1., "");
        
        TH1F * ratiostaterr = (TH1F *) hmc_exp->Clone("ratiostaterr");
        ratiostaterr->Sumw2();
        ratiostaterr->SetStats(0);
        ratiostaterr->SetTitle("");
        //ratiostaterr->GetXaxis()->SetTitle(xtitle);
        ratiostaterr->SetTitle(TString(";")+xtitle);
        ratiostaterr->GetYaxis()->SetTitle("Data/MC");
        ratiostaterr->SetMaximum(2.2);
        ratiostaterr->SetMinimum(0);
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
            
            //if (!(hdata_test->GetBinContent(i) > 1e-6)) {
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
                                    pow(0.05 * hWj0b->GetBinContent(i), 2) +
                                    pow(0.20 * hWj1b->GetBinContent(i), 2) +
                                    pow(0.20 * hWj2b->GetBinContent(i), 2) +
                                    pow(0.05 * hZj0b->GetBinContent(i), 2) +
                                    pow(0.15 * hZj1b->GetBinContent(i), 2) +
                                    pow(0.15 * hZj2b->GetBinContent(i), 2) +
                                    pow(0.04 * hTT->GetBinContent(i), 2) +
                                    pow(0.25 * hs_Top->GetBinContent(i), 2) +
#ifndef VVANALYSIS
                                    pow(0.25 * hVVLF->GetBinContent(i), 2) +
                                    pow(0.25 * hVVHF->GetBinContent(i), 2));
#else                                    
                                    pow(0.25 * hVVLF->GetBinContent(i), 2));
#endif
                double binerror = sqrt(binerror2);
                ratiosysterr->SetBinError(i, binerror / hmc_exp->GetBinContent(i));
                
                if (hmc_test->GetBinContent(i) > 1e-6)  // use smaller tolerance?
                    hmc_test->SetBinError(i, binerror);
                if (hdata_test->GetBinContent(i) > 1e-6)
                    hdata_test->SetBinError(i, sqrt(hdata_test->GetBinContent(i)));
            }
        }

        // Setup legends
        std::clog << "MakePlots(): Setting up legends..." << std::endl;
//        TLegend * leg = new TLegend(0.74, 0.56, 0.92, 0.92);
//        leg->SetFillColor(0);
//        leg->SetLineColor(0);
//        leg->SetShadowColor(0);
//        leg->SetTextFont(62);
//        leg->SetTextSize(0.03);
//        leg->AddEntry(hdata_test, "Data", "p");
//        if (plotSig)  leg->AddEntry(hVH, Form("VH(%i)", massH), "l");
//#ifdef VVANALYSIS
//        leg->AddEntry(hVVHF, "VV(b#bar{b})", "l");
//#endif
//        leg->AddEntry(hZj2b, "Z + b#bar{b}", "f");
//        leg->AddEntry(hZj1b, "Z + b", "f");
//        leg->AddEntry(hZj0b, "Z + udscg", "f");
//        leg->AddEntry(hWj2b, "W + b#bar{b}", "f");
//        leg->AddEntry(hWj1b, "W + b", "f");
//        leg->AddEntry(hWj0b, "W + udscg", "f");
//        leg->AddEntry(hTT, "t#bar{t}", "f");
//        leg->AddEntry(hs_Top, "single top", "f");
//        leg->AddEntry(hQCD, "QCD", "f");
//        leg->AddEntry(hVVLF, "VV(udscg)", "f");
//#ifndef VVANALYSIS
//        leg->AddEntry(hVVHF, "VV(b#bar{b})", "f");
//#endif
//        leg->AddEntry(staterr, "MC uncert. (stat)", "fl");

        TLegend * leg1 = new TLegend(0.58, 0.68, 0.76, 0.92);
        leg1->SetFillColor(0);
        leg1->SetLineColor(0);
        leg1->SetShadowColor(0);
        leg1->SetTextFont(62);
        leg1->SetTextSize(0.03);
        leg1->AddEntry(hdata_test, "Data", "p");
        if (plotSig)  leg1->AddEntry(hVH, Form("VH(%i)", massH), "l");
#ifdef VVANALYSIS
        leg1->AddEntry(hVVHF, "VV(b#bar{b})", "l");
#endif
        leg1->AddEntry(hTT, "t#bar{t}", "f");
        leg1->AddEntry(hs_Top, "single top", "f");
        leg1->AddEntry(hQCD, "QCD", "f");
        leg1->AddEntry(hVVLF, "VV(udscg)", "f");
#ifndef VVANALYSIS
        leg1->AddEntry(hVVHF, "VZ(b#bar{b})", "f");
#endif

        TLegend * leg2 = new TLegend(0.76, 0.68, 0.94, 0.92);
        leg2->SetFillColor(0);
        leg2->SetLineColor(0);
        leg2->SetShadowColor(0);
        leg2->SetTextFont(62);
        leg2->SetTextSize(0.03);
        leg2->AddEntry(hWj2b, "W + b#bar{b}", "f");
        leg2->AddEntry(hWj1b, "W + b", "f");
        leg2->AddEntry(hWj0b, "W + udscg", "f");
        leg2->AddEntry(hZj2b, "Z + b#bar{b}", "f");
        leg2->AddEntry(hZj1b, "Z + b", "f");
        leg2->AddEntry(hZj0b, "Z + udscg", "f");
        leg2->AddEntry(staterr, "MC uncert. (stat)", "f");
        

        //TLegend * ratioleg1 = new TLegend(0.54, 0.86, 0.72, 0.96);
        TLegend * ratioleg1 = new TLegend(0.54, 0.88, 0.72, 0.96);
        //TLegend * ratioleg1 = new TLegend(0.50, 0.86, 0.69, 0.96);
        ratioleg1->AddEntry(ratiostaterr, "MC uncert. (stat)", "f");
        ratioleg1->SetFillColor(0);
        ratioleg1->SetLineColor(0);
        ratioleg1->SetShadowColor(0);
        ratioleg1->SetTextFont(62);
        ratioleg1->SetTextSize(0.06);
        ratioleg1->SetBorderSize(1);
        
        //TLegend * ratioleg2 = new TLegend(0.72, 0.86, 0.95, 0.96);
        TLegend * ratioleg2 = new TLegend(0.72, 0.88, 0.95, 0.96);
        //TLegend * ratioleg2 = new TLegend(0.69, 0.86, 0.9, 0.96);
        ratioleg2->AddEntry(ratiosysterr, "MC uncert. (stat+syst)", "f");
        ratioleg2->SetFillColor(0);
        ratioleg2->SetLineColor(0);
        ratioleg2->SetShadowColor(0);
        ratioleg2->SetTextFont(62);
        ratioleg2->SetTextSize(0.06);
        ratioleg2->SetBorderSize(1);

//        //TLegend * ratioleg1 = new TLegend(0.72, 0.86, 0.94, 0.96);
//        TLegend * ratioleg1 = new TLegend(0.76, 0.88, 0.94, 0.96);
//        //TLegend * ratioleg1 = new TLegend(0.50, 0.86, 0.69, 0.96);
//        ratioleg1->AddEntry(ratiostaterr, "MC uncert. (stat)", "f");
//        ratioleg1->SetFillColor(0);
//        ratioleg1->SetLineColor(0);
//        ratioleg1->SetShadowColor(0);
//        ratioleg1->SetTextFont(62);
//        //ratioleg1->SetTextSize(0.06);
//        ratioleg1->SetTextSize(0.07);
//        ratioleg1->SetBorderSize(1);

        // Draw MC signal and background
        std::clog << "MakePlots(): Drawing..." << std::endl;
        pad1->cd();
        if (plotLog) pad1->SetLogy();
        hs->Draw("hist");
        hs->GetXaxis()->SetLabelSize(0);
        double binwidth = (xup - xlow) / nbins_plot;
        TString ytitle = Form("Events / %.3f", binwidth);
        hs->GetYaxis()->SetTitle(ytitle);
        if (TString(xtitle).Contains(" ; "))
            hs->SetTitle(TString(";")+xtitle);
        
        staterr->Draw("e2 same");
        if (plotSig) {
            hVH->SetLineWidth(3);
            hVH->SetFillColor(0);
            hVH->Draw("hist same");
        }
#ifdef VVANALYSIS
        hVVHF->SetLineWidth(3);
        hVVHF->SetLineColor(kGray + 2);
        hVVHF->SetFillColor(0);
        hVVHF->Draw("hist same");
#endif
        
        // Draw data
        hdata_test->Draw("e1 same");
        
        // Draw legends
        //leg->Draw();
        leg1->Draw();
        leg2->Draw();
        TLatex * latex = new TLatex();
        latex->SetNDC();
        latex->SetTextAlign(12);
        latex->SetTextFont(62);
        latex->SetTextSize(0.052);
        latex->DrawLatex(0.19, 0.89, "CMS Preliminary");
        latex->SetTextSize(0.04);
#ifndef HCPANALYSIS
        latex->DrawLatex(0.19, 0.84, "#sqrt{s} = 8 TeV, L = 19.0 fb^{-1}");
#else
        latex->DrawLatex(0.19, 0.84, "#sqrt{s} = 8 TeV, L = 12.3 fb^{-1}");
#endif
        
        if (region == "VH")
            latex->DrawLatex(0.19, 0.79, "Z(#nu#bar{#nu})H(b#bar{b})");
        else if (region == "ZjLF")
            latex->DrawLatex(0.19, 0.79, "Z(#nu#bar{#nu})H(b#bar{b}), Z + udscg enriched");
        else if (region == "ZjHF")
            latex->DrawLatex(0.19, 0.79, "Z(#nu#bar{#nu})H(b#bar{b}), Z + b#bar{b} enriched");
        else if (region == "WjLF")
            latex->DrawLatex(0.19, 0.79, "Z(#nu#bar{#nu})H(b#bar{b}), W + udscg enriched");
        else if (region == "WjHF")
            latex->DrawLatex(0.19, 0.79, "Z(#nu#bar{#nu})H(b#bar{b}), W + b#bar{b} enriched");
        else if (region == "TT")
            latex->DrawLatex(0.19, 0.79, "Z(#nu#bar{#nu})H(b#bar{b}), t#bar{t} enriched");
        
        // Under/overflows a la TMVA
        TString uoflow = Form("U/O-flow (Data,MC): (%.1f, %.1f) / (%.1f, %.1f)", hdata_test->GetBinContent(0), hmc_exp->GetBinContent(0), hdata_test->GetBinContent(nbins_plot+1), hmc_exp->GetBinContent(nbins_plot+1));
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
        //TPaveText * pave = new TPaveText(0.18, 0.86, 0.35, 0.96, "brNDC");
        TPaveText * pave = new TPaveText(0.18, 0.86, 0.28, 0.96, "brNDC");
        pave->SetTextAlign(12);
        pave->SetLineColor(0);
        pave->SetFillColor(0);
        pave->SetShadowColor(0);
        pave->SetBorderSize(1);
        double nchisq = hdata_test->Chi2Test(hmc_test, "UWCHI2/NDF");  // MC uncert. (stat)
        double kolprob = hdata_test->KolmogorovTest(hmc_test);  // MC uncert. (stat)
        //TText * text = pave->AddText(Form("#chi_{#nu}^{2} = %.3f, K_{s} = %.3f", nchisq, kolprob));
        TText * text = pave->AddText(Form("#chi_{#nu}^{2} = %.3f", nchisq));
        text->SetTextFont(62);
        //text->SetTextSize(0.06);
        text->SetTextSize(0.07);
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
        TString plotname = Form("%s_%s_%s", channel.Data(), region.Data(), var.Data());
        FormatFileName(plotname);
        gPad->Print(g_plotdir+plotname+".png");
        gPad->Print(g_plotdir+plotname+".pdf");
        
        delete hdata_test;
        delete hmc_test;
        delete staterr;
        delete ratio;
        delete ratiostaterr;
        delete ratiosysterr;
        //delete leg;
        delete leg1;
        delete leg2;
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

        // Add histograms for stat, WJModel, ZJModel, TTModel, WJSlope, ZJSlope, boost_QCD
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
            
            name = "Wj0b_CMS_vhbb_WJModel_Znn_8TeVUp";
            TH1F * hWj0b_modelUp = VaryModelErrors((TH1F *) hWj0b, (TH1F *) hWj0b_syst, 1.0, name);
            name = "Wj1b_CMS_vhbb_WJModel_Znn_8TeVUp";
            TH1F * hWj1b_modelUp = VaryModelErrors((TH1F *) hWj1b, (TH1F *) hWj1b_syst, 1.0, name);
            name = "Wj2b_CMS_vhbb_WJModel_Znn_8TeVUp";
            TH1F * hWj2b_modelUp = VaryModelErrors((TH1F *) hWj2b, (TH1F *) hWj2b_syst, 1.0, name);
            name = "Wj0b_CMS_vhbb_WJModel_Znn_8TeVDown";
            TH1F * hWj0b_modelDown = VaryModelErrors((TH1F *) hWj0b, (TH1F *) hWj0b_syst, -1.0, name);
            name = "Wj1b_CMS_vhbb_WJModel_Znn_8TeVDown";
            TH1F * hWj1b_modelDown = VaryModelErrors((TH1F *) hWj1b, (TH1F *) hWj1b_syst, -1.0, name);
            name = "Wj2b_CMS_vhbb_WJModel_Znn_8TeVDown";
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
            
            name = "Zj0b_CMS_vhbb_ZJModel_Znn_8TeVUp";
            TH1F * hZj0b_modelUp = VaryModelErrors((TH1F *) hZj0b, (TH1F *) hZj0b_syst, 1.0, name);
            name = "Zj1b_CMS_vhbb_ZJModel_Znn_8TeVUp";
            TH1F * hZj1b_modelUp = VaryModelErrors((TH1F *) hZj1b, (TH1F *) hZj1b_syst, 1.0, name);
            name = "Zj2b_CMS_vhbb_ZJModel_Znn_8TeVUp";
            TH1F * hZj2b_modelUp = VaryModelErrors((TH1F *) hZj2b, (TH1F *) hZj2b_syst, 1.0, name);
            name = "Zj0b_CMS_vhbb_ZJModel_Znn_8TeVDown";
            TH1F * hZj0b_modelDown = VaryModelErrors((TH1F *) hZj0b, (TH1F *) hZj0b_syst, -1.0, name);
            name = "Zj1b_CMS_vhbb_ZJModel_Znn_8TeVDown";
            TH1F * hZj1b_modelDown = VaryModelErrors((TH1F *) hZj1b, (TH1F *) hZj1b_syst, -1.0, name);
            name = "Zj2b_CMS_vhbb_ZJModel_Znn_8TeVDown";
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
            
            name = "TT_CMS_vhbb_TTModel_Znn_8TeVUp";
            TH1F * hTT_modelUp = VaryModelErrors((TH1F *) hTT, (TH1F *) hTT_syst, 1.0, name);
            name = "TT_CMS_vhbb_TTModel_Znn_8TeVDown";
            TH1F * hTT_modelDown = VaryModelErrors((TH1F *) hTT, (TH1F *) hTT_syst, -1.0, name);
            hTT_modelUp->Write(hTT_modelUp->GetName());
            hTT_modelDown->Write(hTT_modelDown->GetName());
            
            delete hTT_modelUp;
            delete hTT_modelDown;
            
            TH1F * hWj0b_slopeUp_0    = new TH1F("Wj0b_slopeUp_0"   , "", nbins, xlow, xup);
            TH1F * hWj1b_slopeUp_0    = new TH1F("Wj1b_slopeUp_0"   , "", nbins, xlow, xup);
            TH1F * hWj2b_slopeUp_0    = new TH1F("Wj2b_slopeUp_0"   , "", nbins, xlow, xup);
            TH1F * hZj0b_slopeUp_0    = new TH1F("Zj0b_slopeUp_0"   , "", nbins, xlow, xup);
            TH1F * hZj1b_slopeUp_0    = new TH1F("Zj1b_slopeUp_0"   , "", nbins, xlow, xup);
            TH1F * hZj2b_slopeUp_0    = new TH1F("Zj2b_slopeUp_0"   , "", nbins, xlow, xup);
            TH1F * hWj0b_slopeDown_0  = new TH1F("Wj0b_slopeDown_0" , "", nbins, xlow, xup);
            TH1F * hWj1b_slopeDown_0  = new TH1F("Wj1b_slopeDown_0" , "", nbins, xlow, xup);
            TH1F * hWj2b_slopeDown_0  = new TH1F("Wj2b_slopeDown_0" , "", nbins, xlow, xup);
            TH1F * hZj0b_slopeDown_0  = new TH1F("Zj0b_slopeDown_0" , "", nbins, xlow, xup);
            TH1F * hZj1b_slopeDown_0  = new TH1F("Zj1b_slopeDown_0" , "", nbins, xlow, xup);
            TH1F * hZj2b_slopeDown_0  = new TH1F("Zj2b_slopeDown_0" , "", nbins, xlow, xup);
            
            cutwjslope = Form("apply_pt_slope(genW.pt, %f, 150)", (WJSlope + 1.0) * WJSlopeErr);
            cutzjslope = Form("apply_pt_slope(genZ.pt, %f, 130)", (ZJSlope + 1.0) * ZJSlopeErr);
            
            TCut cutslope = "1.0";
            cutsf = Form("%f", ev->sf_Wj0b);
            ev->Wj0b->Project("Wj0b_slopeUp_0", var, cutmc * cutsf * cutwjslope);
            cutsf = Form("%f", ev->sf_Wj1b);
            ev->Wj1b->Project("Wj1b_slopeUp_0", var, cutmc * cutsf * cutwjslope);
            cutsf = Form("%f", ev->sf_Wj2b);
            ev->Wj2b->Project("Wj2b_slopeUp_0", var, cutmc * cutsf * cutwjslope);
            cutsf = Form("%f", ev->sf_Zj0b);
            ev->Zj0b->Project("Zj0b_slopeUp_0", var, cutmc * cutsf * cutzjslope);
            cutsf = Form("%f", ev->sf_Zj1b);
            ev->Zj1b->Project("Zj1b_slopeUp_0", var, cutmc * cutsf * cutzjslope);
            cutsf = Form("%f", ev->sf_Zj2b);
            ev->Zj2b->Project("Zj2b_slopeUp_0", var, cutmc * cutsf * cutzjslope);
            
            cutwjslope = Form("apply_pt_slope(genW.pt, %f, 150)", (WJSlope - 1.0) * WJSlopeErr);
            cutzjslope = Form("apply_pt_slope(genZ.pt, %f, 130)", (ZJSlope - 1.0) * ZJSlopeErr);
            
            cutsf = Form("%f", ev->sf_Wj0b);
            ev->Wj0b->Project("Wj0b_slopeDown_0", var, cutmc * cutsf * cutwjslope);
            cutsf = Form("%f", ev->sf_Wj1b);
            ev->Wj1b->Project("Wj1b_slopeDown_0", var, cutmc * cutsf * cutwjslope);
            cutsf = Form("%f", ev->sf_Wj2b);
            ev->Wj2b->Project("Wj2b_slopeDown_0", var, cutmc * cutsf * cutwjslope);
            cutsf = Form("%f", ev->sf_Zj0b);
            ev->Zj0b->Project("Zj0b_slopeDown_0", var, cutmc * cutsf * cutzjslope);
            cutsf = Form("%f", ev->sf_Zj1b);
            ev->Zj1b->Project("Zj1b_slopeDown_0", var, cutmc * cutsf * cutzjslope);
            cutsf = Form("%f", ev->sf_Zj2b);
            ev->Zj2b->Project("Zj2b_slopeDown_0", var, cutmc * cutsf * cutzjslope);
            
            TH1F * hWj0b_slopeUp    = rebinner->rebin(hWj0b_slopeUp_0   , newnbins, "Wj0b_CMS_vhbb_WJSlope_Znn_8TeVUp"    );
            TH1F * hWj1b_slopeUp    = rebinner->rebin(hWj1b_slopeUp_0   , newnbins, "Wj1b_CMS_vhbb_WJSlope_Znn_8TeVUp"    );
            TH1F * hWj2b_slopeUp    = rebinner->rebin(hWj2b_slopeUp_0   , newnbins, "Wj2b_CMS_vhbb_WJSlope_Znn_8TeVUp"    );
            TH1F * hZj0b_slopeUp    = rebinner->rebin(hZj0b_slopeUp_0   , newnbins, "Zj0b_CMS_vhbb_ZJSlope_Znn_8TeVUp"    );
            TH1F * hZj1b_slopeUp    = rebinner->rebin(hZj1b_slopeUp_0   , newnbins, "Zj1b_CMS_vhbb_ZJSlope_Znn_8TeVUp"    );
            TH1F * hZj2b_slopeUp    = rebinner->rebin(hZj2b_slopeUp_0   , newnbins, "Zj2b_CMS_vhbb_ZJSlope_Znn_8TeVUp"    );
            TH1F * hWj0b_slopeDown  = rebinner->rebin(hWj0b_slopeDown_0 , newnbins, "Wj0b_CMS_vhbb_WJSlope_Znn_8TeVDown"  );
            TH1F * hWj1b_slopeDown  = rebinner->rebin(hWj1b_slopeDown_0 , newnbins, "Wj1b_CMS_vhbb_WJSlope_Znn_8TeVDown"  );
            TH1F * hWj2b_slopeDown  = rebinner->rebin(hWj2b_slopeDown_0 , newnbins, "Wj2b_CMS_vhbb_WJSlope_Znn_8TeVDown"  );
            TH1F * hZj0b_slopeDown  = rebinner->rebin(hZj0b_slopeDown_0 , newnbins, "Zj0b_CMS_vhbb_ZJSlope_Znn_8TeVDown"  );
            TH1F * hZj1b_slopeDown  = rebinner->rebin(hZj1b_slopeDown_0 , newnbins, "Zj1b_CMS_vhbb_ZJSlope_Znn_8TeVDown"  );
            TH1F * hZj2b_slopeDown  = rebinner->rebin(hZj2b_slopeDown_0 , newnbins, "Zj2b_CMS_vhbb_ZJSlope_Znn_8TeVDown"  );
            
            hWj0b_slopeUp->Scale(hWj0b->GetSumOfWeights() / hWj0b_slopeUp->GetSumOfWeights());
            hWj1b_slopeUp->Scale(hWj1b->GetSumOfWeights() / hWj1b_slopeUp->GetSumOfWeights());
            hWj2b_slopeUp->Scale(hWj2b->GetSumOfWeights() / hWj2b_slopeUp->GetSumOfWeights());
            hZj0b_slopeUp->Scale(hZj0b->GetSumOfWeights() / hZj0b_slopeUp->GetSumOfWeights());
            hZj1b_slopeUp->Scale(hZj1b->GetSumOfWeights() / hZj1b_slopeUp->GetSumOfWeights());
            hZj2b_slopeUp->Scale(hZj2b->GetSumOfWeights() / hZj2b_slopeUp->GetSumOfWeights());
            hWj0b_slopeDown->Scale(hWj0b->GetSumOfWeights() / hWj0b_slopeDown->GetSumOfWeights());
            hWj1b_slopeDown->Scale(hWj1b->GetSumOfWeights() / hWj1b_slopeDown->GetSumOfWeights());
            hWj2b_slopeDown->Scale(hWj2b->GetSumOfWeights() / hWj2b_slopeDown->GetSumOfWeights());
            hZj0b_slopeDown->Scale(hZj0b->GetSumOfWeights() / hZj0b_slopeDown->GetSumOfWeights());
            hZj1b_slopeDown->Scale(hZj1b->GetSumOfWeights() / hZj1b_slopeDown->GetSumOfWeights());
            hZj2b_slopeDown->Scale(hZj2b->GetSumOfWeights() / hZj2b_slopeDown->GetSumOfWeights());
            
            hWj0b_slopeUp->Write(hWj0b_slopeUp->GetName());
            hWj1b_slopeUp->Write(hWj1b_slopeUp->GetName());
            hWj2b_slopeUp->Write(hWj2b_slopeUp->GetName());
            hZj0b_slopeUp->Write(hZj0b_slopeUp->GetName());
            hZj1b_slopeUp->Write(hZj1b_slopeUp->GetName());
            hZj2b_slopeUp->Write(hZj2b_slopeUp->GetName());
            hWj0b_slopeDown->Write(hWj0b_slopeDown->GetName());
            hWj1b_slopeDown->Write(hWj1b_slopeDown->GetName());
            hWj2b_slopeDown->Write(hWj2b_slopeDown->GetName());
            hZj0b_slopeDown->Write(hZj0b_slopeDown->GetName());
            hZj1b_slopeDown->Write(hZj1b_slopeDown->GetName());
            hZj2b_slopeDown->Write(hZj2b_slopeDown->GetName());
            
            delete hWj0b_slopeUp_0; delete hWj0b_slopeUp;
            delete hWj1b_slopeUp_0; delete hWj1b_slopeUp;
            delete hWj2b_slopeUp_0; delete hWj2b_slopeUp;
            delete hZj0b_slopeUp_0; delete hZj0b_slopeUp;
            delete hZj1b_slopeUp_0; delete hZj1b_slopeUp;
            delete hZj2b_slopeUp_0; delete hZj2b_slopeUp;
            delete hWj0b_slopeDown_0; delete hWj0b_slopeDown;
            delete hWj1b_slopeDown_0; delete hWj1b_slopeDown;
            delete hWj2b_slopeDown_0; delete hWj2b_slopeDown;
            delete hZj0b_slopeDown_0; delete hZj0b_slopeDown;
            delete hZj1b_slopeDown_0; delete hZj1b_slopeDown;
            delete hZj2b_slopeDown_0; delete hZj2b_slopeDown;
            
            TH1F * hZH_vhqcdUp_0      = new TH1F("ZH_vhqcdUp_0"     , "", nbins, xlow, xup);
            TH1F * hWH_vhqcdUp_0      = new TH1F("WH_vhqcdUp_0"     , "", nbins, xlow, xup);
            TH1F * hZH_SM_vhqcdUp_0   = new TH1F("ZH_SM_vhqcdUp_0"  , "", nbins, xlow, xup);
            TH1F * hWH_SM_vhqcdUp_0   = new TH1F("WH_SM_vhqcdUp_0"  , "", nbins, xlow, xup);
            TH1F * hZH_vhqcdDown_0    = new TH1F("ZH_vhqcdDown_0"   , "", nbins, xlow, xup);
            TH1F * hWH_vhqcdDown_0    = new TH1F("WH_vhqcdDown_0"   , "", nbins, xlow, xup);
            TH1F * hZH_SM_vhqcdDown_0 = new TH1F("ZH_SM_vhqcdDown_0", "", nbins, xlow, xup);
            TH1F * hWH_SM_vhqcdDown_0 = new TH1F("WH_SM_vhqcdDown_0", "", nbins, xlow, xup);
            
            TCut cutvhqcdUp   = Form("(1-0.5*(1-%s))", cutvhqcd.GetTitle());
            TCut cutvhqcdDown = Form("(1-1.5*(1-%s))", cutvhqcd.GetTitle());
            ev->ZH->Project("ZH_vhqcdUp_0", var, cutmc * cutvhqcdUp);
            ev->WH->Project("WH_vhqcdUp_0", var, cutmc * cutvhqcdUp);
            ev->ZH_SM->Project("ZH_SM_vhqcdUp_0", var, cutmc * cutvhqcdUp);
            ev->WH_SM->Project("WH_SM_vhqcdUp_0", var, cutmc * cutvhqcdUp);
            ev->ZH->Project("ZH_vhqcdDown_0", var, cutmc * cutvhqcdDown);
            ev->WH->Project("WH_vhqcdDown_0", var, cutmc * cutvhqcdDown);
            ev->ZH_SM->Project("ZH_SM_vhqcdDown_0", var, cutmc * cutvhqcdDown);
            ev->WH_SM->Project("WH_SM_vhqcdDown_0", var, cutmc * cutvhqcdDown);
            
            TH1F * hZH_vhqcdUp      = rebinner->rebin(hZH_vhqcdUp_0     , newnbins, "ZH_CMS_vhbb_boost_QCD_8TeVUp");
            TH1F * hWH_vhqcdUp      = rebinner->rebin(hWH_vhqcdUp_0     , newnbins, "WH_CMS_vhbb_boost_QCD_8TeVUp");
            TH1F * hZH_SM_vhqcdUp   = rebinner->rebin(hZH_SM_vhqcdUp_0  , newnbins, "ZH_SM_CMS_vhbb_boost_QCD_8TeVUp");
            TH1F * hWH_SM_vhqcdUp   = rebinner->rebin(hWH_SM_vhqcdUp_0  , newnbins, "WH_SM_CMS_vhbb_boost_QCD_8TeVUp");
            TH1F * hZH_vhqcdDown    = rebinner->rebin(hZH_vhqcdDown_0   , newnbins, "ZH_CMS_vhbb_boost_QCD_8TeVDown");
            TH1F * hWH_vhqcdDown    = rebinner->rebin(hWH_vhqcdDown_0   , newnbins, "WH_CMS_vhbb_boost_QCD_8TeVDown");
            TH1F * hZH_SM_vhqcdDown = rebinner->rebin(hZH_SM_vhqcdDown_0, newnbins, "ZH_SM_CMS_vhbb_boost_QCD_8TeVDown");
            TH1F * hWH_SM_vhqcdDown = rebinner->rebin(hWH_SM_vhqcdDown_0, newnbins, "WH_SM_CMS_vhbb_boost_QCD_8TeVDown");
            
            hZH_vhqcdUp->Scale(hZH->GetSumOfWeights() / hZH_vhqcdUp->GetSumOfWeights());
            hWH_vhqcdUp->Scale(hWH->GetSumOfWeights() / hWH_vhqcdUp->GetSumOfWeights());
            hZH_SM_vhqcdUp->Scale(hZH_SM->GetSumOfWeights() / hZH_SM_vhqcdUp->GetSumOfWeights());
            hWH_SM_vhqcdUp->Scale(hWH_SM->GetSumOfWeights() / hWH_SM_vhqcdUp->GetSumOfWeights());
            hZH_vhqcdDown->Scale(hZH->GetSumOfWeights() / hZH_vhqcdDown->GetSumOfWeights());
            hWH_vhqcdDown->Scale(hWH->GetSumOfWeights() / hWH_vhqcdDown->GetSumOfWeights());
            hZH_SM_vhqcdDown->Scale(hZH_SM->GetSumOfWeights() / hZH_SM_vhqcdDown->GetSumOfWeights());
            hWH_SM_vhqcdDown->Scale(hWH_SM->GetSumOfWeights() / hWH_SM_vhqcdDown->GetSumOfWeights());
            
            hZH_vhqcdUp->Write(hZH_vhqcdUp->GetName());
            hWH_vhqcdUp->Write(hWH_vhqcdUp->GetName());
            hZH_SM_vhqcdUp->Write(hZH_SM_vhqcdUp->GetName());
            hWH_SM_vhqcdUp->Write(hWH_SM_vhqcdUp->GetName());
            hZH_vhqcdDown->Write(hZH_vhqcdDown->GetName());
            hWH_vhqcdDown->Write(hWH_vhqcdDown->GetName());
            hZH_SM_vhqcdDown->Write(hZH_SM_vhqcdDown->GetName());
            hWH_SM_vhqcdDown->Write(hWH_SM_vhqcdDown->GetName());
            
            delete hZH_vhqcdUp_0; delete hZH_vhqcdUp;
            delete hWH_vhqcdUp_0; delete hWH_vhqcdUp;
            delete hZH_SM_vhqcdUp_0; delete hZH_SM_vhqcdUp;
            delete hWH_SM_vhqcdUp_0; delete hWH_SM_vhqcdUp;
            delete hZH_vhqcdDown_0; delete hZH_vhqcdDown;
            delete hWH_vhqcdDown_0; delete hWH_vhqcdDown;
            delete hZH_SM_vhqcdDown_0; delete hZH_SM_vhqcdDown;
            delete hWH_SM_vhqcdDown_0; delete hWH_SM_vhqcdDown;
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


EventsJ12 * Read(bool isControlRegion, const TCut cutallmc, const TCut cutalldata)
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
#ifdef QCDSHAPE
    QCD.Add(indir + prefix + "QCD" + suffix);
#endif
    
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

#ifdef MJJANALYSIS
    /// Get more MC statistics
    if (treename.EndsWith("test")) {
        const TString traintreename = Form("tree_%s_%s", channel.Data(), "train");
    
        ZH.AddFile(indir + prefix + Form("ZH%i", massH) + suffix, TChain::kBigNumber, traintreename);
        WH.AddFile(indir + prefix + Form("WH%i", massH) + suffix, TChain::kBigNumber, traintreename);
        Wj.AddFile(indir + prefix + "Wj" + suffix, TChain::kBigNumber, traintreename);
        Zj.AddFile(indir + prefix + "Zj" + suffix, TChain::kBigNumber, traintreename);
        TT.AddFile(indir + prefix + "TT" + suffix, TChain::kBigNumber, traintreename);
        s_Top.AddFile(indir + prefix + "s_Top" + suffix, TChain::kBigNumber, traintreename);
        VV.AddFile(indir + prefix + "VV" + suffix, TChain::kBigNumber, traintreename);
#ifdef QCDSHAPE
        QCD.AddFile(indir + prefix + "QCD" + suffix, TChain::kBigNumber, traintreename);
#endif
        Wj_syst.AddFile(indir + prefix + "WjHW" + suffix, TChain::kBigNumber, traintreename);
        Zj_syst.AddFile(indir + prefix + "ZjHW" + suffix, TChain::kBigNumber, traintreename);
        TT_syst.AddFile(indir + prefix + "TTPowheg" + suffix, TChain::kBigNumber, traintreename);
        ZH_SM.AddFile(indir + prefix + "ZH125" + suffix, TChain::kBigNumber, traintreename);
        WH_SM.AddFile(indir + prefix + "WH125" + suffix, TChain::kBigNumber, traintreename);
    }
#endif  // MJJANALYSIS

    // Start copying trees
    const TCut cutHF = "eventFlav==5";
    const TCut cutLF = "eventFlav!=5";
    
    const TCut cutVVHF = "eventFlav==5 && processname!=\"WW\"";
    const TCut cutVVLF = "eventFlav!=5 || processname==\"WW\"";

    const TCut cut2b = "abs(hJet_flavour[0])==5 && abs(hJet_flavour[1])==5";
    const TCut cut1b = "(abs(hJet_flavour[0])==5 && abs(hJet_flavour[1])!=5) || (abs(hJet_flavour[0])!=5 && abs(hJet_flavour[1])==5)";
    const TCut cut0b = "abs(hJet_flavour[0])!=5 && abs(hJet_flavour[1])!=5";
    
    EventsJ12 * ev = new EventsJ12();
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
#ifndef MJJANALYSIS
    ev->Zj0b = (TTree *) Zj.CopyTree(cutallmc + cut0b);
    std::clog << "... DONE: Zj0b copy tree." << std::endl;
    ev->Zj1b = (TTree *) Zj.CopyTree(cutallmc + cut1b);
    std::clog << "... DONE: Zj1b copy tree." << std::endl;
    ev->Zj2b = (TTree *) Zj.CopyTree(cutallmc + cut2b);
    std::clog << "... DONE: Zj2b copy tree." << std::endl;
#else
    ev->Zj0b = (TTree *) Zj.CopyTree(cutallmc + cut0b * "processname!=\"ZJetsHT50\"");
    std::clog << "... DONE: Zj0b copy tree." << std::endl;
    ev->Zj1b = (TTree *) Zj.CopyTree(cutallmc + cut1b * "processname!=\"ZJetsHT50\"");
    std::clog << "... DONE: Zj1b copy tree." << std::endl;
    ev->Zj2b = (TTree *) Zj.CopyTree(cutallmc + cut2b * "processname!=\"ZJetsHT50\"");
    std::clog << "... DONE: Zj2b copy tree." << std::endl;
#endif
    ev->TT = (TTree *) TT.CopyTree(cutallmc);
    std::clog << "... DONE: TT copy tree." << std::endl;
    ev->s_Top = (TTree *) s_Top.CopyTree(cutallmc);
    std::clog << "... DONE: s_Top copy tree." << std::endl;
    ev->VVLF = (TTree *) VV.CopyTree(cutallmc + cutVVLF);
    std::clog << "... DONE: VVLF copy tree." << std::endl;
    ev->VVHF = (TTree *) VV.CopyTree(cutallmc + cutVVHF);
    std::clog << "... DONE: VVHF copy tree." << std::endl;
#if defined(QCDSHAPE) && !defined(MJJANALYSIS)
    ev->QCD = (TTree *) QCD.CopyTree(cutallmc);
    std::clog << "... DONE: QCD copy tree." << std::endl;
#elif defined(QCDSHAPE)
    ev->QCD = (TTree *) QCD.CopyTree(cutallmc * "mindPhiMETJet_dPhi>0.5");
    std::clog << "... DONE: QCD copy tree." << std::endl;
#endif
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
    ev->ZH_SM = (TTree *) ZH_SM.CopyTree(cutallmc);
    std::clog << "... DONE: ZH_SM copy tree." << std::endl;
    ev->WH_SM = (TTree *) WH_SM.CopyTree(cutallmc);
    std::clog << "... DONE: WH_SM copy tree." << std::endl;
    
    ev->data_obs = (TTree *) data_obs.CopyTree(cutalldata);
    std::clog << "... DONE: data_obs copy tree." << std::endl;

    return ev;
}

////////////////////////////////////////////////////////////////////////////////
/// Main                                                                     ///
////////////////////////////////////////////////////////////////////////////////

void BDTShapeJ12(int nbins=500, long long newnbins=22, double rebinerrorf=0.25, bool freerebin=false, bool isControlRegion=false)  
// 2000,10101010, 0.25 for multi-BDT
// 500,16 for VV BDT
// 1500,101010, 0.25 for VV multi-BDT
{
    /// BINNING CHOICES
#ifndef VVANALYSIS
    double binning_BDT [4*3] = {22, 0.25, 0.25, 14, 0.25, 0.25, 11, 0.25, 0.25};
    double binning_mBDT[4*3] = {10101010, 0.15, 0.15, 8080808, 0.15, 0.15, 6060606, 0.15, 0.15};
    double binning_nBDT[4*3] = {30, 0.045, 0.25, 30, 0.025, 0.25, 30, 0.024, 0.22};
#else
    double binning_BDT [4*3] = {20, 0.25, 0.25, 14, 0.25, 0.25, 11, 0.25, 0.25};
    double binning_mBDT[4*3] = {101010, 0.15, 0.15, 80808, 0.15, 0.15, 60606, 0.15, 0.15};
    double binning_nBDT[4*3] = {25, 0.032, 0.35, 25, 0.032, 0.35, 25, 0.025, 0.35};
#endif

    const bool is_mBDT = (newnbins >= 100);
    const bool is_nBDT = (newnbins < 0);
    if (!is_mBDT)  assert(nbins > newnbins);
    if (is_nBDT)   newnbins = -1 * newnbins;
    
    gROOT->LoadMacro("tdrstyle.C");
    gROOT->ProcessLine("setTDRStyle()");
    gROOT->LoadMacro("HelperFunctions.h");
    
    TH1::SetDefaultSumw2(1);
    gROOT->SetBatch(1);
    
    //if (!TString(gROOT->GetVersion()).Contains("5.34")) {
    //    std::cout << "INCORRECT ROOT VERSION! Please use 5.34:" << std::endl;
    //    std::cout << "source /uscmst1/prod/sw/cms/slc5_amd64_gcc472/lcg/root/5.34.03-cms4/bin/thisroot.csh" << std::endl;
    //    std::cout << "Return without doing anything." << std::endl;
    //    return;
    //}

    if (gSystem->AccessPathName(g_plotdir))
        gSystem->mkdir(g_plotdir);

    const int begin = (isControlRegion) ?   1 : 0;
    const int end   = (isControlRegion) ? 1+5 : 1;
    
    
    ///-- ZnunuHighPt ----------------------------------------------------------
    channel = "ZnunuHighPt";
    std::clog << "\n*** " << channel << " ***\n" << std::endl;
    
    EventsJ12 * ev = Read(isControlRegion, g_cutallmc, g_cutalldata);
    ev->check();
    double errorffirst = rebinerrorf;
    double errorflast = rebinerrorf;
    if (!freerebin) {
        nbins       = 500;
        newnbins    = binning_BDT[0];
        errorffirst = binning_BDT[1];
        errorflast  = binning_BDT[2];
        if (is_mBDT) {
            nbins       = 2000;
#ifdef VVANALYSIS
            nbins       = 1500;
#endif
            newnbins    = binning_mBDT[0];
            errorffirst = binning_mBDT[1];
            errorflast  = binning_mBDT[2];
        } else if (is_nBDT) {
            nbins       = 750;
            newnbins    = binning_nBDT[0];
            errorffirst = binning_nBDT[1];
            errorflast  = binning_nBDT[2];
        }
    }

    /// Select scale factors
    const double * scalefactors     = g_scalefactors_ZnunuHighPt;
    const double * scalefactors_lnN = g_scalefactors_lnN_ZnunuHighPt;

    for (int ireg=begin; ireg<end; ireg++) {
        if (ireg==0) {  // signal region
            /// Reset rebinner
            delete rebinner; rebinner = 0;
            
            /// Make BDT plots
            /// Loop over systematics
            for (int isyst = 0; isyst < nsyst; isyst++) {
                ev->set_scalefactors(&scalefactors[0]);
                std::clog << "isyst " << isyst << ": scale factors: " << ev->sf_Wj0b << " (Wj0b), " << ev->sf_Wj1b << " (Wj1b), " << ev->sf_Wj2b << " (Wj2b), " << ev->sf_Zj0b << " (Zj0b), " << ev->sf_Zj1b << " (Zj1b), " << ev->sf_Zj2b << " (Zj2b), " << ev->sf_TT << " (TT)." << std::endl;
                
                /// Replace "ISYST" with the index of systematic
                TString var = (!is_mBDT) ? g_var : g_varslice;
                TString ssyst=Form("%i",isyst);
                var.ReplaceAll("ISYST",ssyst);
                
                TCut cutmc_test   = Form("2.0 * weightsMC[%i] * (selectFlags[0][%i]) * 19040/19624", isyst, isyst); // FIXME
                TCut cutdata_test = Form("selectFlags[0][%i]", isyst);
                
                MakePlots(ev, var, cutmc_test, cutdata_test, g_systematics[isyst], g_regions[ireg], "BDT", nbins, g_xlow, g_xup, newnbins, errorffirst, errorflast, scalefactors_lnN, "printStat:printCard:writeRoot:!plotData:plotLog:plotSig");
                
                //TCut lastpartition = "(BDTtt_125[0]>-0.5 && BDTvjlf_125[0]>-0.5 && BDTzz_125[0]>-0.3)";
                //MakePlots(ev, var+"+0", cutmc_test*lastpartition, cutdata_test*lastpartition, g_systematics[isyst], g_regions[ireg], "BDT", nbins, g_xlow, g_xup, newnbins, errorffirst, errorflast, scalefactors_lnN, "printStat:printCard:writeRoot:!plotData:plotLog:plotSig");
            }
        }
        
        if (g_manyplots) {
            /// Make plots of 13 input variables
            ev->set_scalefactors(&scalefactors[0]);
            TCut cutmc_ctrl = Form("weightsMC[0] * (selectFlags[%i][0]) * 19040/19624", ireg);  // FIXME
            TCut cutdata_ctrl = Form("(selectFlags[%i][0])", ireg);
            TCut cutmass_ctrl = "";
            if (ireg==0) {  // signal region
                cutmc_ctrl = Form("2.0 * weightsMC[0] * (selectFlags[%i][0]) * 19040/19624", ireg);  // FIXME
                cutmass_ctrl = "(HmassReg<100 || 140<HmassReg)";
            }
            
            if (ireg!=0) {  // control region
                MakePlots(ev, Form("BDTregular_%i[0]", massH), cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "BDT", nbins, g_xlow, g_xup, newnbins, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:plotLog:plotSig");
#ifdef VVANALYSIS
                MakePlots(ev, Form("BDTregular__sigzzhf_%i[0]", massH), cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "VV BDT", nbins, g_xlow, g_xup, newnbins, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:plotLog:plotSig");
#endif
            }
            //MakePlots(ev, "HmassReg", cutmc_ctrl*cutmass_ctrl, cutdata_ctrl*cutmass_ctrl, g_systematics[0], g_regions[ireg], "H mass [GeV]", 17, 0, 255, 17, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");  // apply mass veto
            MakePlots(ev, "HmassReg", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H mass [GeV]", 17, 0, 255, 17, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");  // don't apply mass veto
            MakePlots(ev, "HptReg", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H p_{T} [GeV]", 23, 115, 460, 23, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "max(hJet_ptReg[0],hJet_ptReg[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 p_{T} [GeV]", 20, 50, 350, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "min(hJet_ptReg[0],hJet_ptReg[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j2 p_{T} [GeV]", 13, 15, 210, 13, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "#Delta R(j1,j2)", 14, 0, 3.5, 14, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "abs(deltaPhi(hJet_phi[0], hJet_phi[1]))", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "#Delta #phi(j1,j2)", 16, 0, 3.2, 16, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "abs(hJet_eta[0]-hJet_eta[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "|#Delta #eta|(j1,j2)", 14, 0, 3.5, 14, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "min(abs(deltaPullAngle),pi)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "|#Delta #theta_{pull}|(j1,j2)", 16, 0, 3.2, 16, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "METtype1corr.et", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "Type-1 corr. pfMET [GeV]", 20, 155, 455, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "HMETdPhi", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "|#Delta #varphi|(H, pfMET)", 16, 1.6, 3.2, 16, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "mindPhiMETJet_dPhi", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "min |#Delta #varphi|(pfMET, j25)", 17, -0.1, 3.3, 17, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "max(hJet_csv_nominal[0],hJet_csv_nominal[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "CSV_{max}(j1,j2)", 30, 0.0, 1.08, 30, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "min(hJet_csv_nominal[0],hJet_csv_nominal[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "CSV_{min}(j1,j2)", 30, 0.0, 1.08, 30, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "naJets_Znn", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "# of add. jets", 8, 0, 8, 8, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
            MakePlots(ev, "MaxIf$(max(aJet_csv_nominal,0), aJet_pt>20 && abs(aJet_eta)<2.5)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "CSV_{max}(add. cj20)", 30, 0.0, 1.08, 30, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "min(mindRHJet_dR,5.653)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "min #DeltaR(H, add. j25)", 23, 0, 5.75, 23, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            
            // Bkg specific BDT
            MakePlots(ev, Form("BDTtt_%i[0]", massH), cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "BDTtt", 15, -1, 1, 15, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:plotLog:plotSig");
            MakePlots(ev, Form("BDTvjlf_%i[0]", massH), cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "BDTvjlf", 15, -1, 1, 15, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:plotLog:plotSig");
            MakePlots(ev, Form("BDTzz_%i[0]", massH), cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "BDTzz", 15, -1, 1, 15, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:plotLog:plotSig");
            //MakePlots(ev, Form("BDTtt__sigzzhf_%i[0]", massH), cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "BDTtt", 15, -1, 1, 15, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:plotLog:plotSig");
            //MakePlots(ev, Form("BDTvjlf__sigzzhf_%i[0]", massH), cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "BDTvjlf", 15, -1, 1, 15, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:plotLog:plotSig");
            //MakePlots(ev, Form("BDTzz__sigzzhf_%i[0]", massH), cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "BDTzz", 15, -1, 1, 15, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:plotLog:plotSig");
            
            // QCD check
            MakePlots(ev, "METtype1corr.et/sqrt(METtype1corr.sumet)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "Type-1 corr. pfMET/sqrt(sumET) [GeV]", 20, 0, 18, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "METtype1corr.sig", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "Type-1 corr. pfMET significance", 18, 0, 360, 18, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "nPVs", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "# of primary vertices", 20, 0, 40, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "mindPhiMETCtrJet_dPhi", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "min |#Delta #varphi|(pfMET, cj20)", 17, -0.1, 3.3, 17, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            
            // Angle (need the definition of the functions)
            //MakePlots(ev, "abs(evalCosThetaHbb(hJet_ptReg[0], hJet_pt[0], hJet_eta[0], hJet_phi[0], hJet_e[0], hJet_ptReg[1], hJet_pt[1], hJet_eta[1], hJet_phi[1], hJet_e[1]))", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "cos #theta* (H,j)", 20, 0, 1, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            //MakePlots(ev, "evalHMETMt(HptReg, H.phi, METtype1corr.et, METtype1corr.phi)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "mass_{T}(H,pfMET)", 25, 250, 1000, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            //MakePlots(ev, "evalHMETMassiveMt(HptReg, H.phi, METtype1corr.et, METtype1corr.phi)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "massive mass_{T}(H,pfMET)", 25, 250, 1000, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            
            // FatH.filteredmass
            MakePlots(ev, "FatH.filteredmass * (nfathFilterJets>0 && FatH.FatHiggsFlag==1)", cutmc_ctrl*cutmass_ctrl, cutdata_ctrl*cutmass_ctrl, g_systematics[0], g_regions[ireg], "H_{fj} mass [GeV]", 25, 0, 250, 25, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");  // apply mass veto
            //MakePlots(ev, "FatH.filteredpt * (nfathFilterJets>0 && FatH.FatHiggsFlag==1)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H_{fj} p_{T} [GeV]", 20, 50, 450, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            
            // FatHmassReg
            MakePlots(ev, "FatHmassReg * (nfathFilterJets>0 && FatH.FatHiggsFlag==1)", cutmc_ctrl*cutmass_ctrl, cutdata_ctrl*cutmass_ctrl, g_systematics[0], g_regions[ireg], "H_{fj} mass [GeV]", 25, 0, 250, 25, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");  // apply mass veto
            //MakePlots(ev, "FatHptReg * (nfathFilterJets>0 && FatH.FatHiggsFlag==1)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H_{fj} p_{T} [GeV]", 20, 50, 450, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            //MakePlots(ev, "fathFilterJets_ptReg[0] * (nfathFilterJets>0 && FatH.FatHiggsFlag==1)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H fj1 p_{T} [GeV]", 20, 50, 350, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            //MakePlots(ev, "fathFilterJets_ptReg[1] * (nfathFilterJets>0 && FatH.FatHiggsFlag==1)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H fj2 p_{T} [GeV]", 10, 20, 170, 10, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            
            // Out of the box
            //MakePlots(ev, "H.mass", cutmc_ctrl*cutmass_ctrl, cutdata_ctrl*cutmass_ctrl, g_systematics[0], g_regions[ireg], "H mass [GeV]", 17, 0, 255, 17, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");  // apply mass veto
            MakePlots(ev, "H.mass", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H mass [GeV]", 17, 0, 255, 17, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");  // don't apply mass veto
            MakePlots(ev, "H.pt", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H p_{T} [GeV] (no regression)", 23, 115, 460, 23, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "hJet_pt[0]", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 p_{T} [GeV] (no regression)", 20, 50, 350, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "hJet_pt[1]", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j2 p_{T} [GeV] (no regression)", 13, 15, 210, 13, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "MET.et", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "Type-1 corr. pfMET [GeV] (no JEC update)", 20, 155, 455, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "V.pt", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "V p_{T} [GeV] (no JEC update)", 20, 155, 455, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "Sum$(aJet_pt>30 && abs(aJet_eta)<4.5 && aJet_id==1 && aJet_puJetIdL>0)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "# of add. jets (p_{T} > 30 GeV)", 8, 0, 8, 8, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");

/*        
            if (g_regions[ireg]=="TT") {
                // Regression inputs
                MakePlots(ev, "nPVs", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "# of primary vertices (UNUSED)", 20, 0, 40, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
                MakePlots(ev, "rho25", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "energy density in |#eta|<2.5 [GeV] (UNUSED)", 20, 0, 30, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
                MakePlots(ev, "hJet_ptRaw[0]", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 raw p_{T} [GeV] (UNUSED)", 20, 50, 350, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
                MakePlots(ev, "smear_pt_res(hJet_ptRaw[0], hJet_genPt[0], hJet_eta[0])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 raw (JER corr.) p_{T} [GeV]", 20, 50, 350, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
                MakePlots(ev, "hJet_pt[0]", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 p_{T} [GeV]", 20, 50, 350, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
                MakePlots(ev, "evalEt(hJet_pt[0], hJet_eta[0], hJet_phi[0], hJet_e[0])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 E_{T} [GeV]", 20, 50, 350, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
                MakePlots(ev, "evalMt(hJet_pt[0], hJet_eta[0], hJet_phi[0], hJet_e[0])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 m_{T} [GeV]", 20, 50, 350, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
                MakePlots(ev, "hJet_eta[0]", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 #eta (UNUSED)", 24, -3.0, 3.0, 24, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
                MakePlots(ev, "hJet_ptLeadTrack[0]", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 leading track p_{T} [GeV]", 20, 0, 150, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
                MakePlots(ev, "max(0,hJet_vtx3dL[0])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 sec vtx 3-d flight length [cm]", 15, 0, 6, 15, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
                MakePlots(ev, "max(0,hJet_vtx3deL[0])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 sec vtx 3-d flight length error [cm]", 15, 0, 0.6, 15, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
                MakePlots(ev, "max(0,hJet_vtxMass[0])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 sec vtx mass [GeV]", 20, 0, 6.5, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
                MakePlots(ev, "max(0,hJet_vtxPt[0])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 sec vtx p_{T} [GeV]", 24, 0, 120, 24, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
                MakePlots(ev, "hJet_chf[0]", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 charged hadronic frac (UNUSED)", 21, 0, 1.05, 21, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
                MakePlots(ev, "hJet_cef[0]", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 charged electromagnetic frac", 21, 0, 1.05, 21, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
                MakePlots(ev, "hJet_nhf[0]", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 neutral hadronic frac (UNUSED)", 21, 0, 1.05, 21, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
                MakePlots(ev, "hJet_nef[0]", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 neutral electromagnetic frac(UNUSED)", 21, 0, 1.05, 21, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
                MakePlots(ev, "hJet_nch[0]", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 # charged particles (UNUSED)", 20, 0, 100, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
                MakePlots(ev, "hJet_nconstituents[0]", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 # constituents [GeV]", 24, 0, 120, 24, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
                MakePlots(ev, "hJet_JECUnc[0]", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 JEC uncertainty [GeV]", 25, 0.005, 0.025, 25, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
                MakePlots(ev, "max(0,hJet_SoftLeptptRel[0]*(hJet_SoftLeptIdlooseMu[0]==1 || hJet_SoftLeptId95[0]==1))", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 soft lepton p_{T, rel} [GeV]", 20, 0, 40, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
                MakePlots(ev, "max(0,hJet_SoftLeptPt[0]*(hJet_SoftLeptIdlooseMu[0]==1 || hJet_SoftLeptId95[0]==1))", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 soft lepton p_{T} [GeV]", 20, 0, 40, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
                MakePlots(ev, "max(0,hJet_SoftLeptdR[0]*(hJet_SoftLeptIdlooseMu[0]==1 || hJet_SoftLeptId95[0]==1))", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 soft lepton #DeltaR", 12, 0, 0.6, 12, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
                MakePlots(ev, "hJet_csv_nominal[0]", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 CSV", 30, 0.0, 1.08, 30, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
            }
*/

        }
    }  // end loop over regions
    delete ev;


    ///-- ZnunuMedPt -----------------------------------------------------------
    channel = "ZnunuMedPt";
    std::clog << "\n*** " << channel << " ***\n" << std::endl;

    ev = Read(isControlRegion, g_cutallmc, g_cutalldata);
    ev->check();
    if (!freerebin) {
        newnbins    = binning_BDT[3];
        errorffirst = binning_BDT[4];
        errorflast  = binning_BDT[5];
        if (is_mBDT) {
            newnbins    = binning_mBDT[3];
            errorffirst = binning_mBDT[4];
            errorflast  = binning_mBDT[5];
        } else if (is_nBDT) {
            newnbins    = binning_nBDT[3];
            errorffirst = binning_nBDT[4];
            errorflast  = binning_nBDT[5];
        }
    }

    /// Select scale factors
    scalefactors     = g_scalefactors_ZnunuMedPt;
    scalefactors_lnN = g_scalefactors_lnN_ZnunuMedPt;

    for (int ireg=begin; ireg<end; ireg++) {
        if (ireg==0) {  // signal region
            /// Reset rebinner
            delete rebinner; rebinner = 0;
            
            /// Make BDT plots
            /// Loop over systematics
            for (int isyst = 0; isyst < nsyst; isyst++) {
                ev->set_scalefactors(&scalefactors[0]);
                std::clog << "isyst " << isyst << ": scale factors: " << ev->sf_Wj0b << " (Wj0b), " << ev->sf_Wj1b << " (Wj1b), " << ev->sf_Wj2b << " (Wj2b), " << ev->sf_Zj0b << " (Zj0b), " << ev->sf_Zj1b << " (Zj1b), " << ev->sf_Zj2b << " (Zj2b), " << ev->sf_TT << " (TT)." << std::endl;
                
                /// Replace "ISYST" with the index of systematic
                TString var = (!is_mBDT) ? g_var : g_varslice;
                TString ssyst=Form("%i",isyst);
                var.ReplaceAll("ISYST",ssyst);
                
                TCut cutmc_test   = Form("2.0 * weightsMC[%i] * (selectFlags[0][%i]) * 19040/19624", isyst, isyst); // FIXME
                TCut cutdata_test = Form("selectFlags[0][%i]", isyst);
                
                MakePlots(ev, var, cutmc_test, cutdata_test, g_systematics[isyst], g_regions[ireg], "BDT", nbins, g_xlow, g_xup, newnbins, errorffirst, errorflast, scalefactors_lnN, "printStat:printCard:writeRoot:!plotData:plotLog:plotSig");
            }
        }
        
        if (g_manyplots) {
            /// Make plots of 13 input variables
            ev->set_scalefactors(&scalefactors[0]);
            TCut cutmc_ctrl = Form("weightsMC[0] * (selectFlags[%i][0]) * 19040/19624", ireg);  // FIXME
            TCut cutdata_ctrl = Form("(selectFlags[%i][0])", ireg);
            TCut cutmass_ctrl = "";
            if (ireg==0) {  // signal region
                cutmc_ctrl = Form("2.0 * weightsMC[0] * (selectFlags[%i][0]) * 19040/19624", ireg);  // FIXME
                cutmass_ctrl = "(HmassReg<100 || 140<HmassReg)";
            }
            
            if (ireg!=0) {  // control region
                MakePlots(ev, Form("BDTregular_%i[0]", massH), cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "BDT", nbins, g_xlow, g_xup, newnbins, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:plotLog:plotSig");
#ifdef VVANALYSIS
                MakePlots(ev, Form("BDTregular__sigzzhf_%i[0]", massH), cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "VV BDT", nbins, g_xlow, g_xup, newnbins, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:plotLog:plotSig");
#endif
            }
            //MakePlots(ev, "HmassReg", cutmc_ctrl*cutmass_ctrl, cutdata_ctrl*cutmass_ctrl, g_systematics[0], g_regions[ireg], "H mass [GeV]", 17, 0, 255, 17, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");  // apply mass veto
            MakePlots(ev, "HmassReg", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H mass [GeV]", 17, 0, 255, 17, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");  // don't apply mass veto
            MakePlots(ev, "HptReg", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H p_{T} [GeV]", 23, 115, 460, 23, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "max(hJet_ptReg[0],hJet_ptReg[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 p_{T} [GeV]", 20, 50, 350, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "min(hJet_ptReg[0],hJet_ptReg[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j2 p_{T} [GeV]", 13, 15, 210, 13, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "#Delta R(j1,j2)", 14, 0, 3.5, 14, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "abs(deltaPhi(hJet_phi[0], hJet_phi[1]))", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "#Delta #phi(j1,j2)", 16, 0, 3.2, 16, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "abs(hJet_eta[0]-hJet_eta[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "|#Delta #eta|(j1,j2)", 14, 0, 3.5, 14, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "min(abs(deltaPullAngle),pi)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "|#Delta #theta_{pull}|(j1,j2)", 16, 0, 3.2, 16, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "METtype1corr.et", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "Type-1 corr. pfMET [GeV]", 12, 120, 180, 12, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "HMETdPhi", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "|#Delta #varphi|(H, pfMET)", 16, 1.6, 3.2, 16, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "mindPhiMETJet_dPhi", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "min |#Delta #varphi|(pfMET, j25)", 17, -0.1, 3.3, 17, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "max(hJet_csv_nominal[0],hJet_csv_nominal[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "CSV_{max}(j1,j2)", 30, 0.0, 1.08, 30, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "min(hJet_csv_nominal[0],hJet_csv_nominal[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "CSV_{min}(j1,j2)", 30, 0.0, 1.08, 30, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "naJets_Znn", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "# of add. jets", 8, 0, 8, 8, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
            MakePlots(ev, "MaxIf$(max(aJet_csv_nominal,0), aJet_pt>20 && abs(aJet_eta)<2.5)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "CSV_{max}(add. cj20)", 30, 0.0, 1.08, 30, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "min(mindRHJet_dR,5.653)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "min #DeltaR(H, add. j25)", 23, 0, 5.75, 23, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            
            // Bkg specific BDT
            MakePlots(ev, Form("BDTtt_%i[0]", massH), cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "BDTtt", 15, -1, 1, 15, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:plotLog:plotSig");
            MakePlots(ev, Form("BDTvjlf_%i[0]", massH), cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "BDTvjlf", 15, -1, 1, 15, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:plotLog:plotSig");
            MakePlots(ev, Form("BDTzz_%i[0]", massH), cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "BDTzz", 15, -1, 1, 15, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:plotLog:plotSig");
            
            // QCD check
            MakePlots(ev, "METtype1corr.et/sqrt(METtype1corr.sumet)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "Type-1 corr. pfMET/sqrt(sumET) [GeV]", 20, 0, 12, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "METtype1corr.sig", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "Type-1 corr. pfMET significance", 15, 0, 270, 15, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "nPVs", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "# of primary vertices", 20, 0, 40, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "mindPhiMETCtrJet_dPhi", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "min |#Delta #varphi|(pfMET, cj20)", 17, -0.1, 3.3, 17, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            
            // Out of the box
            //MakePlots(ev, "H.mass", cutmc_ctrl*cutmass_ctrl, cutdata_ctrl*cutmass_ctrl, g_systematics[0], g_regions[ireg], "H mass [GeV]", 17, 0, 255, 17, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");  // apply mass veto
            MakePlots(ev, "H.mass", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H mass [GeV]", 17, 0, 255, 17, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");  // don't apply mass veto
            MakePlots(ev, "H.pt", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H p_{T} [GeV] (no regression)", 23, 115, 460, 23, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "hJet_pt[0]", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 p_{T} [GeV] (no regression)", 20, 50, 350, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "hJet_pt[1]", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j2 p_{T} [GeV] (no regression)", 13, 15, 210, 13, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "MET.et", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "Type-1 corr. pfMET [GeV] (no JEC update)", 12, 120, 180, 12, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "V.pt", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "V p_{T} [GeV] (no JEC update)", 20, 120, 220, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "Sum$(aJet_pt>30 && abs(aJet_eta)<4.5 && aJet_id==1 && aJet_puJetIdL>0)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "# of add. jets (p_{T} > 30 GeV)", 8, 0, 8, 8, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
        }
    }  // end loop over regions
    delete ev;


    ///-- ZnunuLowPt -----------------------------------------------------------
    channel = "ZnunuLowPt";
    std::clog << "\n*** " << channel << " ***\n" << std::endl;

    ev = Read(isControlRegion, g_cutallmc, g_cutalldata);
    ev->check();
    if (!freerebin) {
        newnbins    = binning_BDT[6];
        errorffirst = binning_BDT[7];
        errorflast  = binning_BDT[8];
        if (is_mBDT) {
            newnbins    = binning_mBDT[6];
            errorffirst = binning_mBDT[7];
            errorflast  = binning_mBDT[8];
        } else if (is_nBDT) {
            newnbins    = binning_nBDT[6];
            errorffirst = binning_nBDT[7];
            errorflast  = binning_nBDT[8];
        }
    }

    /// Select scale factors
    scalefactors     = g_scalefactors_ZnunuLowPt;
    scalefactors_lnN = g_scalefactors_lnN_ZnunuLowPt;
    
    for (int ireg=begin; ireg<end; ireg++) {
        if (ireg==0) {  // signal region
            /// Reset rebinner
            delete rebinner; rebinner = 0;
            
            /// Make BDT plots
            /// Loop over systematics
            for (int isyst = 0; isyst < nsyst; isyst++) {
                ev->set_scalefactors(&scalefactors[0]);
                std::clog << "isyst " << isyst << ": scale factors: " << ev->sf_Wj0b << " (Wj0b), " << ev->sf_Wj1b << " (Wj1b), " << ev->sf_Wj2b << " (Wj2b), " << ev->sf_Zj0b << " (Zj0b), " << ev->sf_Zj1b << " (Zj1b), " << ev->sf_Zj2b << " (Zj2b), " << ev->sf_TT << " (TT)." << std::endl;
                
                /// Replace "ISYST" with the index of systematic
                TString var = (!is_mBDT) ? g_var : g_varslice;
                TString ssyst=Form("%i",isyst);
                var.ReplaceAll("ISYST",ssyst);
                
                TCut cutmc_test   = Form("2.0 * weightsMC[%i] * (selectFlags[0][%i]) * 19040/19624", isyst, isyst); // FIXME
                TCut cutdata_test = Form("selectFlags[0][%i]", isyst);
                
                MakePlots(ev, var, cutmc_test, cutdata_test, g_systematics[isyst], g_regions[ireg], "BDT", nbins, g_xlow, g_xup, newnbins, errorffirst, errorflast, scalefactors_lnN, "printStat:printCard:writeRoot:!plotData:plotLog:plotSig");
            }
        }
        
        if (g_manyplots) {
            /// Make plots of 13 input variables
            ev->set_scalefactors(&scalefactors[0]);
            TCut cutmc_ctrl = Form("weightsMC[0] * (selectFlags[%i][0]) * 19040/19624", ireg);  // FIXME
            TCut cutdata_ctrl = Form("(selectFlags[%i][0])", ireg);
            TCut cutmass_ctrl = "";
            if (ireg==0) {  // signal region
                cutmc_ctrl = Form("2.0 * weightsMC[0] * (selectFlags[%i][0]) * 19040/19624", ireg);  // FIXME
                cutmass_ctrl = "(HmassReg<100 || 140<HmassReg)";
            }
            
            if (ireg!=0) {  // control region
                MakePlots(ev, Form("BDTregular_%i[0]", massH), cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "BDT", nbins, g_xlow, g_xup, newnbins, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:plotLog:plotSig");
#ifdef VVANALYSIS
                MakePlots(ev, Form("BDTregular__sigzzhf_%i[0]", massH), cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "VV BDT", nbins, g_xlow, g_xup, newnbins, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:plotLog:plotSig");
#endif
            }
            //MakePlots(ev, "HmassReg", cutmc_ctrl*cutmass_ctrl, cutdata_ctrl*cutmass_ctrl, g_systematics[0], g_regions[ireg], "H mass [GeV]", 17, 0, 255, 17, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");  // apply mass veto
            MakePlots(ev, "HmassReg", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H mass [GeV]", 17, 0, 255, 17, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");  // don't apply mass veto
            MakePlots(ev, "HptReg", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H p_{T} [GeV]", 24, 45, 405, 24, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "max(hJet_ptReg[0],hJet_ptReg[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 p_{T} [GeV]", 20, 15, 315, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "min(hJet_ptReg[0],hJet_ptReg[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j2 p_{T} [GeV]", 13, 15, 210, 13, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "#Delta R(j1,j2)", 14, 0, 3.5, 14, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "abs(deltaPhi(hJet_phi[0], hJet_phi[1]))", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "#Delta #phi(j1,j2)", 16, 0, 3.2, 16, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "abs(hJet_eta[0]-hJet_eta[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "|#Delta #eta|(j1,j2)", 14, 0, 3.5, 14, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "min(abs(deltaPullAngle),pi)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "|#Delta #theta_{pull}|(j1,j2)", 16, 0, 3.2, 16, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "METtype1corr.et", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "Type-1 corr. pfMET [GeV]", 10, 90, 140, 10, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "HMETdPhi", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "|#Delta #varphi|(H, pfMET)", 16, 1.6, 3.2, 16, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "mindPhiMETJet_dPhi", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "min |#Delta #varphi|(pfMET, j25)", 17, -0.1, 3.3, 17, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "max(hJet_csv_nominal[0],hJet_csv_nominal[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "CSV_{max}(j1,j2)", 30, 0.0, 1.08, 30, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "min(hJet_csv_nominal[0],hJet_csv_nominal[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "CSV_{min}(j1,j2)", 30, 0.0, 1.08, 30, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "naJets_Znn", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "# of add. jets", 8, 0, 8, 8, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
            MakePlots(ev, "MaxIf$(max(aJet_csv_nominal,0), aJet_pt>20 && abs(aJet_eta)<2.5)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "CSV_{max}(add. cj20)", 30, 0.0, 1.08, 30, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "min(mindRHJet_dR,5.653)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "min #DeltaR(H, add. j25)", 23, 0, 5.75, 23, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            
            // Bkg specific BDT
            MakePlots(ev, Form("BDTtt_%i[0]", massH), cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "BDTtt", 15, -1, 1, 15, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:plotLog:plotSig");
            MakePlots(ev, Form("BDTvjlf_%i[0]", massH), cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "BDTvjlf", 15, -1, 1, 15, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:plotLog:plotSig");
            MakePlots(ev, Form("BDTzz_%i[0]", massH), cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "BDTzz", 15, -1, 1, 15, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:plotLog:plotSig");
            
            // QCD check
            MakePlots(ev, "METtype1corr.et/sqrt(METtype1corr.sumet)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "Type-1 corr. pfMET/sqrt(sumET) [GeV]", 20, 0, 10, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "METtype1corr.sig", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "Type-1 corr. pfMET significance", 12, 0, 180, 12, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "nPVs", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "# of primary vertices", 20, 0, 40, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "mindPhiMETCtrJet_dPhi", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "min |#Delta #varphi|(pfMET, cj20)", 17, -0.1, 3.3, 17, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            
            // Out of the box
            //MakePlots(ev, "H.mass", cutmc_ctrl*cutmass_ctrl, cutdata_ctrl*cutmass_ctrl, g_systematics[0], g_regions[ireg], "H mass [GeV]", 17, 0, 255, 17, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");  // apply mass veto
            MakePlots(ev, "H.mass", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H mass [GeV]", 17, 0, 255, 17, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");  // don't apply mass veto
            MakePlots(ev, "H.pt", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H p_{T} [GeV] (no regression)", 24, 45, 405, 24, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "hJet_pt[0]", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 p_{T} [GeV] (no regression)", 20, 15, 315, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "hJet_pt[1]", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j2 p_{T} [GeV] (no regression)", 13, 15, 210, 13, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "MET.et", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "Type-1 corr. pfMET [GeV] (no JEC update)", 10, 90, 140, 10, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "V.pt", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "V p_{T} [GeV] (no JEC update)", 14, 90, 160, 14, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "Sum$(aJet_pt>30 && abs(aJet_eta)<4.5 && aJet_id==1 && aJet_puJetIdL>0)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "# of add. jets (p_{T} > 30 GeV)", 8, 0, 8, 8, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
        }
    }  // end loop over regions
    delete ev;

/*
    ///-- ZnunuLowCSV ----------------------------------------------------------
    /// FIXME: names given by g_regions[ireg] are wrong
    channel = "ZnunuLowCSV";
    std::clog << "\n*** " << channel << " ***\n" << std::endl;

    ev = Read(isControlRegion, g_cutallmc, g_cutalldata);
    ev->check();
    if (!freerebin) {  // use hardcoded binning choice
        newnbins = (!is_mBDT) ? 11 : 10101010;
    }

    /// Select scale factors
    scalefactors     = g_scalefactors_ZnunuLowCSV;
    scalefactors_lnN = g_scalefactors_lnN_ZnunuLowCSV;
    
    for (int ireg=begin; ireg<TMath::Min(end,3); ireg++) {
        if (ireg==0) {  // signal region
            /// Reset rebinner
            delete rebinner; rebinner = 0;
            
            /// Make BDT plots
            /// Loop over systematics
            for (int isyst = 0; isyst < nsyst; isyst++) {
                ev->set_scalefactors(&scalefactors[0]);
                std::clog << "isyst " << isyst << ": scale factors: " << ev->sf_Wj0b << " (Wj0b), " << ev->sf_Wj1b << " (Wj1b), " << ev->sf_Wj2b << " (Wj2b), " << ev->sf_Zj0b << " (Zj0b), " << ev->sf_Zj1b << " (Zj1b), " << ev->sf_Zj2b << " (Zj2b), " << ev->sf_TT << " (TT)." << std::endl;
                
                /// Replace "ISYST" with the index of systematic
                TString var = (!is_mBDT) ? g_var : g_varslice;
                TString ssyst=Form("%i",isyst);
                var.ReplaceAll("ISYST",ssyst);
                
                TCut cutmc_test   = Form("2.0 * weightsMC[%i] * (selectFlags[0][%i])", isyst, isyst);
                TCut cutdata_test = Form("selectFlags[0][%i]", isyst);

                MakePlots(ev, var, cutmc_test, cutdata_test, g_systematics[isyst], g_regions[ireg], "BDT", nbins, g_xlow, g_xup, newnbins, errorffirst, errorflast, scalefactors_lnN, "printStat:printCard:writeRoot:!plotData:plotLog:plotSig");
            }
        }
        
        if (g_manyplots) {
            /// Make plots of 13 input variables
            ev->set_scalefactors(&scalefactors[0]);
            TCut cutmc_ctrl = Form("weightsMC[0] * (selectFlags[%i][0])", ireg);
            TCut cutdata_ctrl = Form("(selectFlags[%i][0])", ireg);
            TCut cutmass_ctrl = "";
            if (ireg==0) {  // signal region
                cutmc_ctrl = Form("2.0 * weightsMC[0] * (selectFlags[%i][0])", ireg);
                cutmass_ctrl = "(HmassReg<100 || 140<HmassReg)";
            }
            
            if (ireg!=0) {  // control region
                MakePlots(ev, Form("BDTregular_%i[0]", massH), cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "BDT", nbins, g_xlow, g_xup, newnbins, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:plotLog:plotSig");
            }
            MakePlots(ev, "HmassReg", cutmc_ctrl*cutmass_ctrl, cutdata_ctrl*cutmass_ctrl, g_systematics[0], g_regions[ireg], "H mass [GeV]", 17, 0, 255, 17, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");  // apply mass veto
            MakePlots(ev, "HptReg", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H p_{T} [GeV]", 23, 115, 460, 23, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "max(hJet_ptReg[0],hJet_ptReg[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j1 p_{T} [GeV]", 20, 50, 350, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "min(hJet_ptReg[0],hJet_ptReg[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "H j2 p_{T} [GeV]", 13, 15, 210, 13, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "#Delta R(j1,j2)", 14, 0, 3.5, 14, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "abs(deltaPhi(hJet_phi[0], hJet_phi[1]))", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "#Delta #phi(j1,j2)", 16, 0, 3.2, 16, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "abs(hJet_eta[0]-hJet_eta[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "|#Delta #eta|(j1,j2)", 14, 0, 3.5, 14, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "min(abs(deltaPullAngle),pi)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "|#Delta #theta_{pull}|(j1,j2)", 16, 0, 3.2, 16, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "METtype1corr.et", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "Type-1 corr. pfMET [GeV]", 20, 155, 455, 20, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "HMETdPhi", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "|#Delta #varphi|(H, pfMET)", 16, 1.6, 3.2, 16, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "mindPhiMETJet_dPhi", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "min |#Delta #varphi|(pfMET, j25)", 17, -0.1, 3.3, 17, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "max(hJet_csv_nominal[0],hJet_csv_nominal[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "CSV_{max}(j1,j2)", 30, 0.0, 1.08, 30, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "min(hJet_csv_nominal[0],hJet_csv_nominal[1])", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "CSV_{min}(j1,j2)", 30, 0.0, 1.08, 30, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "naJets_Znn", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "# of add. jets", 8, 0, 8, 8, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:plotData:!plotLog:plotSig");
            MakePlots(ev, "MaxIf$(max(aJet_csv_nominal,0), aJet_pt>20 && abs(aJet_eta)<2.5)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "CSV_{max}(add. cj20)", 30, 0.0, 1.08, 30, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
            MakePlots(ev, "min(mindRHJet_dR,5.653)", cutmc_ctrl, cutdata_ctrl, g_systematics[0], g_regions[ireg], "min #DeltaR(H, add. j25)", 23, 0, 5.75, 23, errorffirst, errorflast, scalefactors_lnN, "!printStat:!printCard:!writeRoot:plotData:!plotLog:plotSig");
        }
    }  // end loop over regions
    delete ev;
*/
}
