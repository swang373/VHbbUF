////////////////////////////////////////////////////////////////////////////////
/// This macro runs on Step 4 with jmax = 9 in the following order
///   ZH         WH         WjLF       WjHF       ZjLF       ZjHF       TT         s_Top      VV         QCD
///   -1         0          1          2          3          4          5          6          7          8
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

//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130214/reload_20130216/";
const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130302/stitch/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130221/stitch/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130214/stitch/";
//const TString indir     = "root://eoscms//eos/cms/store/caf/user/degrutto/2013/Step4_20130207/";
const TString prefix    = "Step4_";
const TString suffix    = ".root";

const TCut    g_cutallmc    = "";          // used in CopyTree()
const TCut    g_cutalldata  = g_cutallmc;  // used in CopyTree()
const double  g_xlow        = -1.;
const double  g_xup         = 1.;

const int jmax              = 9;  // as in datacard
const int nsyst             = 13;
const int massH             = 125;
//const TString g_var         = Form("BDTregular_%i[ISYST]", massH);  // ISYST is replaced in the main function.
//const TString g_varslice    = Form("slice(BDTtt_%i[ISYST], BDTvjlf_%i[ISYST], BDTzz_%i[ISYST], BDTregular_%i[ISYST], -0.6, -0.7, -0.5)", massH, massH, massH, massH);
////const TString g_varslice    = Form("slice(BDTttbar_%i[ISYST], BDTvjl_%i[ISYST], BDTzz_%i[ISYST], BDTregular_%i[ISYST], -0.5, -0.4, -0.4)", massH, massH, massH, massH);
const TString g_dcname      = "vhbb_Znn_SF_$CHANNEL_8TeV.txt";
const TString g_wsname      = "vhbb_Znn_SF_$CHANNEL_8TeV.root";
const TString g_rootname    = "vhbb_Znn_SF_$CHANNEL_TH1.root";
const TString g_plotdir     = "plots/";

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

/// Scale factors of WjLF, WjHF, ZjLF, ZjHF, TT (as in datacard)
const double g_scalefactors_prefit[5] = {
    1.00, 1.30, 1.00, 1.30, 1.00
};
const double g_scalefactors_lnN_prefit[5] = {
    1.10, 1.40, 1.10, 1.20, 1.08
};
////////////////////////////////////////////////////////////////////////////////
/// END Configuration                                                        ///
////////////////////////////////////////////////////////////////////////////////

///_____________________________________________________________________________
/// Classes

class EventsJ9 {  // jmax = 9 in datacard, with 5 scale factors.
public:
    EventsJ9()
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
        ZjLF_syst(0),
        ZjHF_syst(0),
        TT_syst(0),
        ZH_SM(0),
        WH_SM(0),
        data_obs(0),
        sf_WjLF(1.), sf_WjHF(1.), sf_ZjLF(1.), sf_ZjHF(1.), sf_TT(1.) {}

    ~EventsJ9()
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
        delete ZjLF_syst;
        delete ZjHF_syst;
        delete TT_syst;
        delete ZH_SM;
        delete WH_SM;
        delete data_obs;
    }

    void check() const {
        assert(ZH        != 0 && ZH       ->GetEntriesFast() > 0);
        assert(WH        != 0 && WH       ->GetEntriesFast() > 0);
        assert(WjLF      != 0 && WjLF     ->GetEntriesFast() > 0);
        assert(WjHF      != 0 && WjHF     ->GetEntriesFast() > 0);
        assert(ZjLF      != 0 && ZjLF     ->GetEntriesFast() > 0);
        assert(ZjHF      != 0 && ZjHF     ->GetEntriesFast() > 0);
        assert(TT        != 0 && TT       ->GetEntriesFast() > 0);
        assert(s_Top     != 0 && s_Top    ->GetEntriesFast() > 0);
        assert(VVLF      != 0 && VVLF     ->GetEntriesFast() > 0);
        assert(VVHF      != 0 && VVHF     ->GetEntriesFast() > 0);
        //assert(QCD       != 0 && QCD      ->GetEntriesFast() > 0);
        assert(ZjLF_syst != 0 && ZjLF_syst->GetEntriesFast() > 0);
        assert(ZjHF_syst != 0 && ZjHF_syst->GetEntriesFast() > 0);
        assert(TT_syst   != 0 && TT_syst  ->GetEntriesFast() > 0);
        //assert(ZH_SM     != 0 && ZH_SM    ->GetEntriesFast() > 0);
        //assert(WH_SM     != 0 && WH_SM    ->GetEntriesFast() > 0);
        assert(data_obs  != 0 && data_obs ->GetEntriesFast() > 0);
        return;
    }
    
    void set_scale_factors(const double sf[5]) {
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
    TTree * ZjLF_syst;  // for ZJModel
    TTree * ZjHF_syst;  // for ZJModel
    TTree * TT_syst;    // for TTModel
    TTree * ZH_SM;      // for signal injection
    TTree * WH_SM;      // for signal injection
    TTree * data_obs;
    Double_t sf_WjLF;
    Double_t sf_WjHF;
    Double_t sf_ZjLF;
    Double_t sf_ZjHF;
    Double_t sf_TT;
};  // EventsJ9


///_____________________________________________________________________________
/// Functions

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

/// Check consistency between two histograms
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

// N.B. other channels might take absolute difference
TH1F * VaryModelErrors(const TH1F * h, const TH1F * hmodel, Double_t c, const char* newname="", bool absdiff=false, bool doUOflow=false)
{
    TH1F * hdummy = (TH1F *) h->Clone("dummy");
    hdummy->Sumw2();
    Int_t nbins = h->GetNbinsX();
    //if (TMath::Abs(c*c - 1.0) > 1e-10){
    //    hdummy->Error("Vary", "c must be +1 or -1!");  // allow c to be any constant?
    //    return hdummy;
    //}
    if (h == hmodel){
        hdummy->Warning("Vary", "The input histograms are the same. Return a clone. h=%s", h->GetName());
        hdummy->SetName(newname);
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
        hdummy->SetName(newname);
        return hnew;
    }
    if (hmodel->GetSumOfWeights() < 1e-10){
        hdummy->Warning("Vary", "Empty model histogram. Return a clone. hmodel=%s", hmodel->GetName());
        hdummy->SetName(newname);
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
        Double_t       diff           = bincontentm - bincontent;
        if (absdiff)   diff           = fabs(diff);
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

void MakePlots(const EventsJ9 * ev, const TString var, const TString var2, 
               const TCut cutmc, const TCut cutdata, const TString syst, const TString region, 
               const char* axtitle, int nbins, double xlow, double xup, 
               long long newnbins, double rebinerrorf, 
               const double scalefactors_lnN[5], 
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
    TH1F * hWjLF_0      = new TH1F("WjLF_0"     , "", nbins, xlow, xup);
    TH1F * hWjHF_0      = new TH1F("WjHF_0"     , "", nbins, xlow, xup);
    TH1F * hZjLF_0      = new TH1F("ZjLF_0"     , "", nbins, xlow, xup);
    TH1F * hZjHF_0      = new TH1F("ZjHF_0"     , "", nbins, xlow, xup);
    TH1F * hTT_0        = new TH1F("TT_0"       , "", nbins, xlow, xup);
    TH1F * hs_Top_0     = new TH1F("s_Top_0"    , "", nbins, xlow, xup);
    TH1F * hVVLF_0      = new TH1F("VVLF_0"     , "", nbins, xlow, xup);
    TH1F * hVVHF_0      = new TH1F("VVHF_0"     , "", nbins, xlow, xup);
    TH1F * hQCD_0       = new TH1F("QCD_0"      , "", nbins, xlow, xup);
    TH1F * hZjLF_syst_0 = new TH1F("ZjLF_syst_0", "", nbins, xlow, xup);
    TH1F * hZjHF_syst_0 = new TH1F("ZjHF_syst_0", "", nbins, xlow, xup);
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
    histos_0.push_back(hWjLF_0);
    histos_0.push_back(hWjHF_0);
    histos_0.push_back(hZjLF_0);
    histos_0.push_back(hZjHF_0);
    histos_0.push_back(hTT_0);
    histos_0.push_back(hs_Top_0);
    histos_0.push_back(hVVLF_0);
    histos_0.push_back(hVVHF_0);
    histos_0.push_back(hQCD_0);
    histos_0.push_back(hZjLF_syst_0);
    histos_0.push_back(hZjHF_syst_0);
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
    
    ev->ZjLF_syst->Project("ZjLF_syst_0", var, cutmc);
    std::clog << "... DONE: project ZjLF_syst_0." << std::endl;
    ev->ZjHF_syst->Project("ZjHF_syst_0", var, cutmc);
    std::clog << "... DONE: project ZjHF_syst_0." << std::endl;
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
    cutsf = Form("%f", ev->sf_WjLF);
    ev->WjLF->Project("+WjLF_0", var2, cutmc * cutsf);
    std::clog << "... DONE: project WjLF_0, scaled by " << cutsf << "." << std::endl;
    cutsf = Form("%f", ev->sf_WjHF);
    ev->WjHF->Project("+WjHF_0", var2, cutmc * cutsf);
    std::clog << "... DONE: project WjHF_0, scaled by " << cutsf << "." << std::endl;
    cutsf = Form("%f", ev->sf_ZjLF);
    ev->ZjLF->Project("+ZjLF_0", var2, cutmc * cutsf);
    std::clog << "... DONE: project ZjLF_0, scaled by " << cutsf << "." << std::endl;
    cutsf = Form("%f", ev->sf_ZjHF);
    ev->ZjHF->Project("+ZjHF_0", var2, cutmc * cutsf);
    std::clog << "... DONE: project ZjHF_0, scaled by " << cutsf << "." << std::endl;
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
    
    ev->ZjLF_syst->Project("+ZjLF_syst_0", var2, cutmc);
    std::clog << "... DONE: project ZjLF_syst_0." << std::endl;
    ev->ZjHF_syst->Project("+ZjHF_syst_0", var2, cutmc);
    std::clog << "... DONE: project ZjHF_syst_0." << std::endl;
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

    hmc_exp_0->Add(hWjLF_0);
    hmc_exp_0->Add(hWjHF_0);
    hmc_exp_0->Add(hZjLF_0);
    hmc_exp_0->Add(hZjHF_0);
    hmc_exp_0->Add(hTT_0);
    hmc_exp_0->Add(hs_Top_0);
    hmc_exp_0->Add(hVV_0);
    //hmc_exp_0->Add(hQCD_0);
    std::clog << "... DONE: add MC backgrounds to mc_exp_0." << std::endl;

    // Do the rebinning
    TH1F * hZH, * hWH, * hWjLF, * hWjHF, * hZjLF, * hZjHF, * hTT, * hs_Top, * hVVLF, * hVVHF, * hQCD,
         * hZjLF_syst, * hZjHF_syst, * hTT_syst, * hZH_SM, * hWH_SM, * hVH, * hVV, * hmc_exp, * hdata_obs;

    if (!(var.Contains("slice"))) {
        hmc_exp    = RebinByStatErrors(hmc_exp_0   , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "mc_exp");
        hZH        = RebinByStatErrors(hZH_0       , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "ZH");
        hWH        = RebinByStatErrors(hWH_0       , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "WH");
        hWjLF      = RebinByStatErrors(hWjLF_0     , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "WjLF");
        hWjHF      = RebinByStatErrors(hWjHF_0     , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "WjHF");
        hZjLF      = RebinByStatErrors(hZjLF_0     , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "ZjLF");
        hZjHF      = RebinByStatErrors(hZjHF_0     , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "ZjHF");
        hTT        = RebinByStatErrors(hTT_0       , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "TT");
        hs_Top     = RebinByStatErrors(hs_Top_0    , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "s_Top");
        hVVLF      = RebinByStatErrors(hVVLF_0     , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "VVLF");
        hVVHF      = RebinByStatErrors(hVVHF_0     , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "VVHF");
        hQCD       = RebinByStatErrors(hQCD_0      , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "QCD");
        hZjLF_syst = RebinByStatErrors(hZjLF_syst_0, newnbins, firstN, lastN, rebinerrorf, xlow, xup, "ZjLF_syst");
        hZjHF_syst = RebinByStatErrors(hZjHF_syst_0, newnbins, firstN, lastN, rebinerrorf, xlow, xup, "ZjHF_syst");
        hTT_syst   = RebinByStatErrors(hTT_syst_0  , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "TT_syst");
        //hZH_SM     = RebinByStatErrors(hZH_SM_0    , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "ZH_SM");
        //hWH_SM     = RebinByStatErrors(hWH_SM_0    , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "WH_SM");
        hVH        = RebinByStatErrors(hVH_0       , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "VH");
        hVV        = RebinByStatErrors(hVV_0       , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "VV");
        hdata_obs  = RebinByStatErrors(hdata_obs_0 , newnbins, firstN, lastN, rebinerrorf, xlow, xup, "data_obs");
    } else {
        // newnbins decides the binning in each slice (as opposed to the total).
        hmc_exp    = RebinByStatErrors(hmc_exp_0   , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "mc_exp");
        hZH        = RebinByStatErrors(hZH_0       , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "ZH");
        hWH        = RebinByStatErrors(hWH_0       , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "WH");
        hWjLF      = RebinByStatErrors(hWjLF_0     , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "WjLF");
        hWjHF      = RebinByStatErrors(hWjHF_0     , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "WjHF");
        hZjLF      = RebinByStatErrors(hZjLF_0     , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "ZjLF");
        hZjHF      = RebinByStatErrors(hZjHF_0     , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "ZjHF");
        hTT        = RebinByStatErrors(hTT_0       , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "TT");
        hs_Top     = RebinByStatErrors(hs_Top_0    , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "s_Top");
        hVVLF      = RebinByStatErrors(hVVLF_0     , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "VVLF");
        hVVHF      = RebinByStatErrors(hVVHF_0     , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "VVHF");
        hQCD       = RebinByStatErrors(hQCD_0      , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "QCD");
        hZjLF_syst = RebinByStatErrors(hZjLF_syst_0, newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "ZjLF_syst");
        hZjHF_syst = RebinByStatErrors(hZjHF_syst_0, newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "ZjHF_syst");
        hTT_syst   = RebinByStatErrors(hTT_syst_0  , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "TT_syst");
        //hZH_SM     = RebinByStatErrors(hZH_SM_0    , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "ZH_SM");
        //hWH_SM     = RebinByStatErrors(hWH_SM_0    , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "WH_SM");
        hVH        = RebinByStatErrors(hVH_0       , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "VH");
        hVV        = RebinByStatErrors(hVV_0       , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "VV");
        hdata_obs  = RebinByStatErrors(hdata_obs_0 , newnbins, firstN1, firstN2, firstN3, firstN4, lastN1, lastN2, lastN3, lastN4, rebinerrorf, xlow, xup, "data_obs");
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
    histos.push_back(hZjLF_syst);
    histos.push_back(hZjHF_syst);
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
        dc << "bin         "; for (Int_t j=0; j!=jmax+1; j++)  dc << channel_8TeV << " "; dc << std::endl;
        dc << "process     ZH         WH         WjLF       WjHF       ZjLF       ZjHF       TT         s_Top      VV         QCD        " << std::endl;
        dc << "process     -1         0          1          2          3          4          5          6          7          8          " << std::endl;
        dc << "rate        " << setw(10) << right << hZH->Integral() << " " << setw(10) << right << hWH->Integral() << " " << setw(10) << right << hWjLF->Integral() << " " << setw(10) << right << hWjHF->Integral() << " " << setw(10) << right << hZjLF->Integral() << " " << setw(10) << right << hZjHF->Integral() << " " << setw(10) << right << hTT->Integral() << " " << setw(10) << right << hs_Top->Integral() << " " << setw(10) << right << hVV->Integral() << " " << setw(10) << right << hQCD->Integral() << std::endl;
        dc.precision(2);
        dc << "-----------------------------------" << std::endl;
        dc << "" << std::endl;
        dc << "### Flat ########################## #####  ZH    WH    WjLF  WjHF  ZjLF  ZjHF  TT    s_Top VV    QCD  " << std::endl;
        dc << "CMS_vhbb_ZH_SF_" << channel_8TeV_1 << "     lnN    1.30  -     -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_WH_SF_" << channel_8TeV_1 << "     lnN    -     1.30  -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_WjLF_SF_" << channel_8TeV_1 << "   lnN    -     -     " << scalefactors_lnN[0] << "  -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_WjHF_SF_" << channel_8TeV_1 << "   lnN    -     -     -     " << scalefactors_lnN[1] << "  -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_ZjLF_SF_" << channel_8TeV_1 << "   lnN    -     -     -     -     " << scalefactors_lnN[2] << "  -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_ZjHF_SF_" << channel_8TeV_1 << "   lnN    -     -     -     -     -     " << scalefactors_lnN[3] << "  -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_TT_SF_" << channel_8TeV_1 << "     lnN    -     -     -     -     -     -     " << scalefactors_lnN[4] << "  -     -     -    " << std::endl;
        dc << "CMS_vhbb_s_Top_SF_" << channel_8TeV_1 << "  lnN    -     -     -     -     -     -     -     1.25  -     -    " << std::endl;
        dc << "CMS_vhbb_VV_SF_" << channel_8TeV_1 << "     lnN    -     -     -     -     -     -     -     -     1.30  -    " << std::endl;
        dc << "#CMS_vhbb_QCD_SF_" << channel_8TeV_1 << "   lnN    -     -     -     -     -     -     -     -     -     1.50 " << std::endl;
        dc << "### Shape ######################### #####  ZH    WH    WjLF  WjHF  ZjLF  ZjHF  TT    s_Top VV    QCD  " << std::endl;
        dc << "UEPS                                shape  1.00  1.00  -     -     -     -     -     1.00  1.00  -    " << std::endl;
        dc << "CMS_vhbb_eff_b_SIG                  shape  1.00  1.00  -     -     -     -     -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_eff_b                      shape  -     -     1.00  1.00  1.00  1.00  1.00  1.00  1.00  -    " << std::endl;
        dc << "CMS_vhbb_fake_b_8TeV                shape  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  -    " << std::endl;
        dc << "CMS_vhbb_res_j                      shape  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  -    " << std::endl;
        dc << "CMS_vhbb_scale_j                    shape  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  -    " << std::endl;
        dc << "CMS_vhbb_trigger_MET_Znunu_8TeV     shape  1.00  1.00  -     -     -     -     -     1.00  1.00  -    " << std::endl;
        dc << "CMS_vhbb_ZJModel_Znunu_8TeV         shape  -     -     -     -     1.00  1.00  -     -     -     -    " << std::endl;
        dc << "CMS_vhbb_TTModel_Znunu_8TeV         shape  -     -     -     -     -     -     1.00  -     -     -    " << std::endl;
        dc << "################################### #####  ZH    WH    WjLF  WjHF  ZjLF  ZjHF  TT    s_Top VV    QCD  " << std::endl;

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
            TH1F * hZjLF_modelUp = VaryModelErrors((TH1F *) hZjLF, (TH1F *) hZjLF_syst, 2.0, name);
            //TH1F * hZjLF_modelUp = VaryModelErrors((TH1F *) hZjLF, (TH1F *) hZjLF_syst, 1.0, name);
            name = "ZjHF_CMS_vhbb_ZJModel_Znunu_8TeVUp";
            TH1F * hZjHF_modelUp = VaryModelErrors((TH1F *) hZjHF, (TH1F *) hZjHF_syst, 2.0, name);
            //TH1F * hZjHF_modelUp = VaryModelErrors((TH1F *) hZjHF, (TH1F *) hZjHF_syst, 1.0, name);
            name = "ZjLF_CMS_vhbb_ZJModel_Znunu_8TeVDown";
            TH1F * hZjLF_modelDown = VaryModelErrors((TH1F *) hZjLF, (TH1F *) hZjLF_syst, 0.0, name);
            //TH1F * hZjLF_modelDown = VaryModelErrors((TH1F *) hZjLF, (TH1F *) hZjLF_syst, -1.0, name);
            name = "ZjHF_CMS_vhbb_ZJModel_Znunu_8TeVDown";
            TH1F * hZjHF_modelDown = VaryModelErrors((TH1F *) hZjHF, (TH1F *) hZjHF_syst, 0.0, name);
            //TH1F * hZjHF_modelDown = VaryModelErrors((TH1F *) hZjHF, (TH1F *) hZjHF_syst, -1.0, name);
            hZjLF_modelUp->Write(hZjLF_modelUp->GetName());
            hZjHF_modelUp->Write(hZjHF_modelUp->GetName());
            hZjLF_modelDown->Write(hZjLF_modelDown->GetName());
            hZjHF_modelDown->Write(hZjHF_modelDown->GetName());
            
            delete hZjLF_modelUp;
            delete hZjHF_modelUp;
            delete hZjLF_modelDown;
            delete hZjHF_modelDown;
            
            name = "TT_CMS_vhbb_TTModel_Znunu_8TeVUp";
            TH1F * hTT_modelUp = VaryModelErrors((TH1F *) hTT, (TH1F *) hTT_syst, 2.0, name);
            //TH1F * hTT_modelUp = VaryModelErrors((TH1F *) hTT, (TH1F *) hTT_syst, 1.0, name);
            name = "TT_CMS_vhbb_TTModel_Znunu_8TeVDown";
            TH1F * hTT_modelDown = VaryModelErrors((TH1F *) hTT, (TH1F *) hTT_syst, 0.0, name);
            //TH1F * hTT_modelDown = VaryModelErrors((TH1F *) hTT, (TH1F *) hTT_syst, -1.0, name);
            
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


EventsJ9 * Read(bool isControlRegion, const TCut cutallmc, const TCut cutalldata)
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
    
    EventsJ9 * ev = new EventsJ9();
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
    ev->ZjLF_syst = (TTree *) Zj_syst.CopyTree(cutallmc + cutLF);
    std::clog << "... DONE: ZjLF_syst copy tree." << std::endl;
    ev->ZjHF_syst = (TTree *) Zj_syst.CopyTree(cutallmc + cutHF);
    std::clog << "... DONE: ZjHF_syst copy tree." << std::endl;
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

void ScaleFactorJ9(Int_t nbins=500, Int_t newnbins=16, Double_t rebinerrorf=0.35, bool pqrs=false, bool isControlRegion=true)
{
    gROOT->LoadMacro("tdrstyle.C");
    gROOT->ProcessLine("setTDRStyle()");
    
    TH1::SetDefaultSumw2(1);
    gROOT->SetBatch(1);
    
    assert(newnbins <= nbins);  // nbins is used in plots before rebinning
    
    //if (!TString(gROOT->GetVersion()).Contains("5.34")) {
    //    std::cout << "INCORRECT ROOT VERSION! Please use 5.34:" << std::endl;
    //    std::cout << "source /uscmst1/prod/sw/cms/slc5_amd64_gcc472/lcg/root/5.34.03-cms4/bin/thisroot.csh" << std::endl;
    //    std::cout << "Return without doing anything." << std::endl;
    //    return;
    //}
    
    double fitresults_ZnunuHighPt[]     = {
         0.000,  0.000,  0.000,  0.000,  0.000,
        -1.035, -1.298,  1.174, -0.379, -1.111, 
        -0.084, -0.433, -0.096, -0.261,  0.437, 
        -0.064, -0.278, -0.072, -0.143,  0.197, 
    };
    double fitresults_lnN_ZnunuHighPt[] = {
         1.000,  1.000,  1.000,  1.000,  1.000,
         0.470,  0.685,  0.403,  0.492,  0.410,
         0.726,  0.859,  0.688,  0.770,  0.754,
         0.855,  0.901,  0.767,  0.831,  0.856,
    };
    
    double fitresults_ZnunuMedPt[]      = {
         0.000,  0.000,  0.000,  0.000,  0.000,
        -0.719, -1.819,  3.966,  2.075, -0.635, 
        -0.238, -0.690,  1.082,  1.132, -0.027, 
        -0.072, -0.340,  0.472,  0.716, -0.013,
    };
    double fitresults_lnN_ZnunuMedPt[]  = {
         1.000,  1.000,  1.000,  1.000,  1.000,
         0.343,  0.704,  0.300,  0.355,  0.324,
         0.735,  0.876,  0.705,  0.701,  0.660, 
         0.717,  0.935,  0.828,  0.803,  0.797, 
    };
    
    double fitresults_ZnunuLowPt[]      = {
         0.000,  0.000,  0.000,  0.000,  0.000,
         1.941, -1.024,  6.366, -0.037, -0.797, 
         0.505, -0.346,  2.600,  0.122,  0.105, 
         0.147, -0.155,  1.611,  0.014,  0.031,
    };
    double fitresults_lnN_ZnunuLowPt[]  = {
         1.000,  1.000,  1.000,  1.000,  1.000,
         0.453,  0.687,  0.439,  0.420,  0.277,
         0.667,  0.811,  0.722,  0.690,  0.742,
         0.819,  0.868,  0.783,  0.766,  0.737,
    };
    
    double fitresults_ZnunuLowCSV[]     = {
         0.000,  0.000,  0.000,  0.000,  0.000,
    };
    double fitresults_lnN_ZnunuLowCSV[] = { 
         1.000,  1.000,  1.000,  1.000,  1.000,
    };
    
    double scalefactors_ZnunuHighPt[5], scalefactors_lnN_ZnunuHighPt[5];
    double scalefactors_ZnunuMedPt [5], scalefactors_lnN_ZnunuMedPt [5];
    double scalefactors_ZnunuLowPt [5], scalefactors_lnN_ZnunuLowPt [5];
    double scalefactors_ZnunuLowCSV[5], scalefactors_lnN_ZnunuLowCSV[5];
    // 0 is the first iteration with prefit scale factors and uncertainties
    calc_scalefactors(5, 2, scalefactors_ZnunuHighPt, scalefactors_lnN_ZnunuHighPt, fitresults_ZnunuHighPt, fitresults_lnN_ZnunuHighPt);
    calc_scalefactors(5, 2, scalefactors_ZnunuMedPt , scalefactors_lnN_ZnunuMedPt , fitresults_ZnunuMedPt , fitresults_lnN_ZnunuMedPt );
    calc_scalefactors(5, 2, scalefactors_ZnunuLowPt , scalefactors_lnN_ZnunuLowPt , fitresults_ZnunuLowPt , fitresults_lnN_ZnunuLowPt );
    calc_scalefactors(5, 0, scalefactors_ZnunuLowCSV, scalefactors_lnN_ZnunuLowCSV, fitresults_ZnunuLowCSV, fitresults_lnN_ZnunuLowCSV);
    

    if (gSystem->AccessPathName(g_plotdir))
        gSystem->mkdir(g_plotdir);

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

    EventsJ9 * ev = Read(isControlRegion, g_cutallmc, g_cutalldata);
    ev->check();

    /// Select scale factors
    const double * scalefactors     = scalefactors_ZnunuHighPt;
    const double * scalefactors_lnN = scalefactors_lnN_ZnunuHighPt;

    for (Int_t ireg=begin; ireg<end; ireg++) {
        channel_CR = channel + "_" + g_regions[ireg];

        /// Initialize rebinning "bookmarks" to negative values
        /// These global variables that need to be reset when we switch channel.
        firstN = -1, lastN = -1;  // for single BDT
        firstN1 = -1, firstN2 = -1, firstN3 = -1, firstN4 = -1; // for PQRS BDT
        lastN1  = -1, lastN2  = -1, lastN3  = -1, lastN4  = -1; // for PQRS BDT

        /// Loop over systematics
        for (Int_t isyst = 0; isyst < nsyst; isyst++) {
            ev->set_scale_factors(scalefactors);
            std::clog << "isyst " << isyst << ": scale factors: " << ev->sf_WjLF << " (WjLF), " << ev->sf_WjHF << " (WjHF), " << ev->sf_ZjLF << " (ZjLF), " << ev->sf_ZjHF << " (ZjHF), " << ev->sf_TT << " (TT)." << std::endl;
            
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


    ///-- ZnunuMedPt -----------------------------------------------------------
    channel = "ZnunuMedPt";
    std::clog << "\n*** " << channel << " ***\n" << std::endl;

    ev = Read(isControlRegion, g_cutallmc, g_cutalldata);
    ev->check();

    /// Select scale factors
    scalefactors     = scalefactors_ZnunuMedPt;
    scalefactors_lnN = scalefactors_lnN_ZnunuMedPt;
    
    for (Int_t ireg=begin; ireg<end; ireg++) {
        channel_CR = channel + "_" + g_regions[ireg];

        /// Initialize rebinning "bookmarks" to negative values
        /// These global variables that need to be reset when we switch channel.
        firstN = -1, lastN = -1;  // for single BDT
        firstN1 = -1, firstN2 = -1, firstN3 = -1, firstN4 = -1; // for PQRS BDT
        lastN1  = -1, lastN2  = -1, lastN3  = -1, lastN4  = -1; // for PQRS BDT

        /// Loop over systematics
        for (Int_t isyst = 0; isyst < nsyst; isyst++) {
            ev->set_scale_factors(scalefactors);
            std::clog << "isyst " << isyst << ": scale factors: " << ev->sf_WjLF << " (WjLF), " << ev->sf_WjHF << " (WjHF), " << ev->sf_ZjLF << " (ZjLF), " << ev->sf_ZjHF << " (ZjHF), " << ev->sf_TT << " (TT)." << std::endl;
            
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


    ///-- ZnunuLowPt -----------------------------------------------------------
    channel = "ZnunuLowPt";
    std::clog << "\n*** " << channel << " ***\n" << std::endl;

    ev = Read(isControlRegion, g_cutallmc, g_cutalldata);
    ev->check();

    /// Select scale factors
    scalefactors     = scalefactors_ZnunuLowPt;
    scalefactors_lnN = scalefactors_lnN_ZnunuLowPt;

    for (Int_t ireg=begin; ireg<end; ireg++) {
        channel_CR = channel + "_" + g_regions[ireg];

        /// Initialize rebinning "bookmarks" to negative values
        /// These global variables that need to be reset when we switch channel.
        firstN = -1, lastN = -1;  // for single BDT
        firstN1 = -1, firstN2 = -1, firstN3 = -1, firstN4 = -1; // for PQRS BDT
        lastN1  = -1, lastN2  = -1, lastN3  = -1, lastN4  = -1; // for PQRS BDT

        /// Loop over systematics
        for (Int_t isyst = 0; isyst < nsyst; isyst++) {
            ev->set_scale_factors(scalefactors);
            std::clog << "isyst " << isyst << ": scale factors: " << ev->sf_WjLF << " (WjLF), " << ev->sf_WjHF << " (WjHF), " << ev->sf_ZjLF << " (ZjLF), " << ev->sf_ZjHF << " (ZjHF), " << ev->sf_TT << " (TT)." << std::endl;
            
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
            
            //TCut cutmc_ctrl = Form("weightsHCP[%i] * (selectFlags[%i][%i])", isyst, ireg, isyst);
            TCut cutmc_ctrl = Form("efflumi * PUweight * triggerFlags[41] * 11299./12101. * (selectFlags[%i][%i])", ireg, isyst);  // FIXME
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

    for (Int_t ireg=begin; ireg<end; ireg++) {
        channel_CR = channel + "_" + g_regions[ireg];
        if (ireg == 1)  // only two control regions for ZnunuLowCSV
            channel_CR = channel + "_" + g_regions[2];
        else if (ireg == 2)
            channel_CR = channel + "_" + g_regions[4];
        else
            continue;

        /// Initialize rebinning "bookmarks" to negative values
        /// These global variables that need to be reset when we switch channel.
        firstN = -1, lastN = -1;  // for single BDT
        firstN1 = -1, firstN2 = -1, firstN3 = -1, firstN4 = -1; // for PQRS BDT
        lastN1  = -1, lastN2  = -1, lastN3  = -1, lastN4  = -1; // for PQRS BDT

        /// Loop over systematics
        for (Int_t isyst = 0; isyst < nsyst; isyst++) {
            ev->set_scale_factors(scalefactors);
            std::clog << "isyst " << isyst << ": scale factors: " << ev->sf_WjLF << " (WjLF), " << ev->sf_WjHF << " (WjHF), " << ev->sf_ZjLF << " (ZjLF), " << ev->sf_ZjHF << " (ZjHF), " << ev->sf_TT << " (TT)." << std::endl;
            
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
    std::cout << scalefactors_ZnunuHighPt[0] << ", " << scalefactors_ZnunuHighPt[1] << ", " << scalefactors_ZnunuHighPt[2] << ", " << scalefactors_ZnunuHighPt[3] << ", " << scalefactors_ZnunuHighPt[4] << std::endl;
    std::cout << fixed << setprecision(2) << "ZnunuHighPt lnN: " << std::endl;
    std::cout << scalefactors_lnN_ZnunuHighPt[0] << ", " << scalefactors_lnN_ZnunuHighPt[1] << ", " << scalefactors_lnN_ZnunuHighPt[2] << ", " << scalefactors_lnN_ZnunuHighPt[3] << ", " << scalefactors_lnN_ZnunuHighPt[4] << std::endl;
    std::cout << fixed << setprecision(3) << "ZnunuMedPt: " << std::endl;
    std::cout << scalefactors_ZnunuMedPt[0] << ", " << scalefactors_ZnunuMedPt[1] << ", " << scalefactors_ZnunuMedPt[2] << ", " << scalefactors_ZnunuMedPt[3] << ", " << scalefactors_ZnunuMedPt[4] << std::endl;
    std::cout << fixed << setprecision(2) << "ZnunuMedPt lnN: " << std::endl;
    std::cout << scalefactors_lnN_ZnunuMedPt[0] << ", " << scalefactors_lnN_ZnunuMedPt[1] << ", " << scalefactors_lnN_ZnunuMedPt[2] << ", " << scalefactors_lnN_ZnunuMedPt[3] << ", " << scalefactors_lnN_ZnunuMedPt[4] << std::endl;
    std::cout << fixed << setprecision(3) << "ZnunuLowPt: " << std::endl;
    std::cout << scalefactors_ZnunuLowPt[0] << ", " << scalefactors_ZnunuLowPt[1] << ", " << scalefactors_ZnunuLowPt[2] << ", " << scalefactors_ZnunuLowPt[3] << ", " << scalefactors_ZnunuLowPt[4] << std::endl;
    std::cout << fixed << setprecision(2) << "ZnunuLowPt lnN: " << std::endl;
    std::cout << scalefactors_lnN_ZnunuLowPt[0] << ", " << scalefactors_lnN_ZnunuLowPt[1] << ", " << scalefactors_lnN_ZnunuLowPt[2] << ", " << scalefactors_lnN_ZnunuLowPt[3] << ", " << scalefactors_lnN_ZnunuLowPt[4] << std::endl;
    std::cout << fixed << setprecision(3) << "ZnunuLowCSV: " << std::endl;
    std::cout << scalefactors_ZnunuLowCSV[0] << ", " << scalefactors_ZnunuLowCSV[1] << ", " << scalefactors_ZnunuLowCSV[2] << ", " << scalefactors_ZnunuLowCSV[3] << ", " << scalefactors_ZnunuLowCSV[4] << std::endl;
    std::cout << fixed << setprecision(2) << "ZnunuLowCSV lnN: " << std::endl;
    std::cout << scalefactors_lnN_ZnunuLowCSV[0] << ", " << scalefactors_lnN_ZnunuLowCSV[1] << ", " << scalefactors_lnN_ZnunuLowCSV[2] << ", " << scalefactors_lnN_ZnunuLowCSV[3] << ", " << scalefactors_lnN_ZnunuLowCSV[4] << std::endl;
}
