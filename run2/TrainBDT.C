#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>

#include "TCut.h"
#include "TFile.h"
#include "TFormula.h"
#include "TH1.h"
#include "TTree.h"
#include "TString.h"
#include "TPRegexp.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Config.h"
#include "TMVA/VariableInfo.h"
#include "TMVA/MethodBDT.h"
//#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

#include "HelperFunctions.h"


////////////////////////////////////////////////////////////////////////////////
/// Configuration                                                            ///
////////////////////////////////////////////////////////////////////////////////

const TString indir     = "/afs/cern.ch/work/d/degrutto/public/MiniAOD/ZnnHbb_Phys14_PU20bx25/skimV11/step3/skim_ZnnH_classification/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130326/stitch/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130314/stitch/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130302/stitch/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130221/stitch/";
//const TString indir     = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130214/stitch/";
//const TString indir     = "root://eoscms//eos/cms/store/caf/user/degrutto/2013/Step4_20130207/";
const TString prefix    = "Step3_";
const TString suffix    = ".root";
const bool g_splitsig   = false;
const bool g_fullbkg    = true;
const bool g_hybrid     = false;

#if 1
/// Default BDT parameters
const int pars_ZnunuHighPt[4] = {300, 5, 3, 10};  // NTrees, MinNodeSize*100  , MaxDepth, AdaBoostBeta*100 <-- NTrees from 600 to 300, use MinNodeSize instead of MinNodeSize, 
const int pars_ZnunuMedPt [4] = {600, 200, 3, 10};  // NTrees, MinNodeSize, MaxDepth, AdaBoostBeta*100
const int pars_ZnunuLowPt [4] = {500, 250, 3, 10};  // NTrees, MinNodeSize, MaxDepth, AdaBoostBeta*100
const int pars_ZnunuLowCSV[4] = {500, 250, 3, 10};  // NTrees, MinNodeSize, MaxDepth, AdaBoostBeta*100
//***************************************************************************************************************
//*    Row      *    NTrees * nEventsMi *  MaxDepth * COMB_limi * COMB_sign * TMVA_kolS * TMVA_kolB * USER_psig *
//***************************************************************************************************************
//* ZnunuHighPt *       600 *       250 *         3 * 2.3828001 * 0.8805900 * 0.9938134 * 0.9105942 * 0.3636768 *
//* ZnunuMedPt  *       600 *       200 *         3 * 4.6406002 * 0.4617460 * 0.0757615 * 0.2361052 * 0.1640445 *
//* ZnunuLowPt  *       500 *       250 *         3 * 6.6561999 * 0.3087579 * 0.9411714 * 0.8741615 * 0.0935148 *
//* combo       *           *           *           *           *           *           *           *           *
//***************************************************************************************************************
#endif

#if 0  // V4a Step4_20130302
/// Default BDT parameters
const int pars_ZnunuHighPt[4] = {500, 300, 3, 10};  // NTrees, MinNodeSize, MaxDepth, AdaBoostBeta*100
const int pars_ZnunuMedPt [4] = {500, 300, 3, 10};  // NTrees, MinNodeSize, MaxDepth, AdaBoostBeta*100
const int pars_ZnunuLowPt [4] = {500, 200, 3, 10};  // NTrees, MinNodeSize, MaxDepth, AdaBoostBeta*100
const int pars_ZnunuLowCSV[4] = {500, 200, 3, 10};  // NTrees, MinNodeSize, MaxDepth, AdaBoostBeta*100
//***************************************************************************************************
//*    Row      *    NTrees * nEventsMi *  MaxDepth * COMB_limi * COMB_sign * TMVA_kolS * TMVA_kolB *
//***************************************************************************************************
//* ZnunuHighPt *       500 *       300 *         3 *     3.125 * 0.7305319 * 0.6633805 * 0.2561788 *
//* ZnunuMedPt  *       500 *       300 *         3 * 6.9061999 * 0.3598870 * 0.2118913 * 0.0062891 *
//* ZnunuLowPt  *       500 *       200 *         3 * 8.4531002 * 0.2538639 * 0.9942439 * 0.6261183 *
//* ZnunuLowCSV *       500 *       200 *         3 *           *           * 0.0821493 * 0.3214474 *
//* combo       *           *           *           *    2.7734 * 0.8081610 *           *           *
//***************************************************************************************************
#endif


/// Variables
UInt_t varcounter = 0;
const TMVA::VariableInfo g_variables[] = {
    TMVA::VariableInfo("HmassReg", 
                       "H mass", "GeV", varcounter++, 'F'),
    TMVA::VariableInfo("HptReg", 
                       "H p_{T}", "GeV", varcounter++, 'F'),
    TMVA::VariableInfo("max(hJet_ptReg[0],hJet_ptReg[1])", 
                       "H j1 p_{T}", "GeV", varcounter++, 'F'),
    TMVA::VariableInfo("min(hJet_ptReg[0],hJet_ptReg[1])", 
                       "H j2 p_{T}", "GeV", varcounter++, 'F'),
    TMVA::VariableInfo("deltaR := deltaR(Jet_eta[hJCidx[0]], Jet_phi[hJCidx[0]], Jet_eta[hJCidx[1]], Jet_phi[hJCidx[1]])", 
                       "#Delta R(j1,j2)", "", varcounter++, 'F'),
    TMVA::VariableInfo("deltaEta := abs(Jet_eta[hJCidx[0]]-Jet_eta[hJCidx[1]])", 
                       "|#Delta #eta|(j1,j2)", "", varcounter++, 'F'),
    //    TMVA::VariableInfo("deltaPullAngle := min(abs(deltaPullAngle),pi)", 
    //                   "|#Delta #theta_{pull}|(j1,j2)", "", varcounter++, 'F'),
    TMVA::VariableInfo("pfMET := met_pt", 
                       "Type-1 corr. pfMET", "GeV", varcounter++, 'F'),
    TMVA::VariableInfo("HMETdPhi := abs(deltaPhi(HCSV_phi , met_phi ))", 
                       "|#Delta #varphi|(H, pfMET)", "", varcounter++, 'F'),
    TMVA::VariableInfo("mindPhiMETJet_dPhi :=  MinIf$(abs(deltaPhi(met_phi,Jet_phi)) , Jet_pt>25 && abs(Jet_eta)<4.5 && Jet_puId==1)", 
                       "min |#Delta #varphi|(pfMET, j25)", "", varcounter++, 'F'),
    TMVA::VariableInfo("maxCSV := max(0,max(Jet_btagCSV[hJCidx[0]],Jet_btagCSV[hJCidx[1]]))", 
                       "CSV_{max}(j1,j2)", "", varcounter++, 'F'),
    TMVA::VariableInfo("minCSV := max(0,min(Jet_btagCSV[hJCidx[0]],Jet_btagCSV[hJCidx[1]]))", 
                       "CSV_{min}(j1,j2)", "", varcounter++, 'F'),
    TMVA::VariableInfo("naJets_Znn : = max(Sum$(Jet_pt>25 && abs(Jet_eta)<4.5 && Jet_puId==1) - 2, 0)", 
                       "# add. jets {p_{T} > 25}", "", varcounter++, 'I'),
    // new since HCP
    TMVA::VariableInfo("maxAddCSV := MaxIf$(max(Jet_btagCSV[aJCidx],0.), Jet_pt[aJCidx]>20 && abs(Jet_eta[aJCidx])<2.5 && Jet_puId[aJCidx]==1 )", 
                       "CSV_{max}(add. cj20)", "", varcounter++, 'F'),
    TMVA::VariableInfo("mindRAddJetH := Min$(deltaR(Jet_eta[hJCidx], Jet_phi[hJCidx], Jet_eta[aJCidx], Jet_phi[aJCidx]))", 
                       "min #DeltaR(H, add. j25)", "", varcounter++, 'F'),
};

//varcounter = 0;
const TMVA::VariableInfo g_variables_NoMjj[] = {
    //TMVA::VariableInfo("HmassReg", 
    //                   "H mass", "GeV", varcounter++, 'F'),
    TMVA::VariableInfo("HptReg", 
                       "H p_{T}", "GeV", varcounter++, 'F'),
    TMVA::VariableInfo("max(hJet_ptReg[0],hJet_ptReg[1])", 
                       "H j1 p_{T}", "GeV", varcounter++, 'F'),
    TMVA::VariableInfo("min(hJet_ptReg[0],hJet_ptReg[1])", 
                       "H j2 p_{T}", "GeV", varcounter++, 'F'),
    //TMVA::VariableInfo("deltaR := deltaR(hJet_eta[0], hJet_phi[0], hJet_eta[1], hJet_phi[1])", 
    //                   "#Delta R(j1,j2)", "", varcounter++, 'F'),
    //TMVA::VariableInfo("deltaEta := abs(hJet_eta[0]-hJet_eta[1])", 
    //                   "|#Delta #eta|(j1,j2)", "", varcounter++, 'F'),
    TMVA::VariableInfo("deltaPullAngle := min(abs(deltaPullAngle),pi)", 
                       "|#Delta #theta_{pull}|(j1,j2)", "", varcounter++, 'F'),
    TMVA::VariableInfo("pfMET := METtype1corr.et", 
                       "Type-1 corr. pfMET", "GeV", varcounter++, 'F'),
    TMVA::VariableInfo("HMETdPhi", 
                       "|#Delta #varphi|(H, pfMET)", "", varcounter++, 'F'),
    TMVA::VariableInfo("mindPhiMETJet_dPhi", 
                       "min |#Delta #varphi|(pfMET, j25)", "", varcounter++, 'F'),
    TMVA::VariableInfo("maxCSV := max(hJet_csv_nominal[0],hJet_csv_nominal[1])", 
                       "CSV_{max}(j1,j2)", "", varcounter++, 'F'),
    TMVA::VariableInfo("minCSV := min(hJet_csv_nominal[0],hJet_csv_nominal[1])", 
                       "CSV_{min}(j1,j2)", "", varcounter++, 'F'),
    TMVA::VariableInfo("naJets_Znn", 
                       "# add. jets {p_{T} > 25}", "", varcounter++, 'I'),
    // new since HCP
    TMVA::VariableInfo("maxAddCSV := MaxIf$(max(aJet_csv_nominal,0), aJet_pt>20 && abs(aJet_eta)<2.5)", 
                       "CSV_{max}(add. cj20)", "", varcounter++, 'F'),
    TMVA::VariableInfo("mindRAddJetH := min(mindRHJet_dR,5.653)", 
                       "min #DeltaR(H, add. j25)", "", varcounter++, 'F'),
};

//varcounter = 0;
// FIXME: use DJ values as default values?
const TMVA::VariableInfo g_variables_FJ[] = {    
    //TMVA::VariableInfo("FatHmass := FatH.filteredmass * (nfathFilterJets>0 && FatH.FatHiggsFlag==1)", 
    //                   "H_{FJ} mass", "GeV", varcounter++, 'F'),  // before regression
    TMVA::VariableInfo("FatHmass := (nfathFilterJets>0 && FatH.FatHiggsFlag==1) ? min(FatH.filteredmass,250) : 250", 
                       "H_{FJ} mass", "GeV", varcounter++, 'F'),  // before regression
    //TMVA::VariableInfo("FatHpt := FatH.filteredpt * (nfathFilterJets>0 && FatH.FatHiggsFlag==1)", 
    //                   "H_{FJ} p_{T}", "GeV", varcounter++, 'F'),  // before regression
    //TMVA::VariableInfo("FatHdeltaM := HmassReg - FatH.filteredmass * (nfathFilterJets>0 && FatH.FatHiggsFlag==1)", 
    //                   "H_{DJ} mass - H_{FJ} mass", "GeV", varcounter++, 'F'),  // before regression
    //TMVA::VariableInfo("FatHdeltaR := min(deltaR(Alt$(fathFilterJets_eta[0],0.), Alt$(fathFilterJets_phi[0],0.), Alt$(fathFilterJets_eta[1],1.06), Alt$(fathFilterJets_phi[1],1.06)) * (nfathFilterJets>0 && FatH.FatHiggsFlag==1) + 999*(!(nfathFilterJets>0 && FatH.FatHiggsFlag==1)),deltaR(0.,0.,1.06,1.06))", 
    //                   "#Delta R(fj1,fj2)", "", varcounter++, 'F'),
};

//varcounter = 0;
const TMVA::VariableInfo g_variables_FJReg[] = {
    //TMVA::VariableInfo("FatHmassReg := FatHmassReg * (nfathFilterJets>0 && FatH.FatHiggsFlag==1)", 
    //                   "H_{FJ} mass", "GeV", varcounter++, 'F'),
    TMVA::VariableInfo("FatHmassReg := (nfathFilterJets>0 && FatH.FatHiggsFlag==1) ? min(FatHmassReg,250) : 250", 
                       "H_{FJ} mass", "GeV", varcounter++, 'F'),
    //TMVA::VariableInfo("FatHptReg := FatHptReg * (nfathFilterJets>0 && FatH.FatHiggsFlag==1)", 
    //                   "H_{FJ} p_{T}", "GeV", varcounter++, 'F'),
    //TMVA::VariableInfo("FatHdeltaM := HmassReg - FatHmassReg * (nfathFilterJets>0 && FatH.FatHiggsFlag==1)", 
    //                   "H_{DJ} mass - H_{FJ} mass", "GeV", varcounter++, 'F'),
};

//varcounter = 0;
const TMVA::VariableInfo g_variables_Angle[] = {
    TMVA::VariableInfo("cosThetaHbb := abs(evalCosThetaHbb(hJet_ptReg[0], hJet_pt[0], hJet_eta[0], hJet_phi[0], hJet_e[0], hJet_ptReg[1], hJet_pt[1], hJet_eta[1], hJet_phi[1], hJet_e[1]))",
                       "cos #theta* (H,j)", "", varcounter++, 'F'),
    TMVA::VariableInfo("HMETMassiveMt := evalHMETMassiveMt(HmassReg, HptReg, H.phi, 91.2, METtype1corr.et, METtype1corr.phi)",
                       "massive mass_{T}(H,pfMET)", "GeV", varcounter++, 'F'),
    //TMVA::VariableInfo("HMETMt := evalHMETMt(HptReg, H.phi, METtype1corr.et, METtype1corr.phi)", 
    //                   "mass_{T}(H,pfMET)", "GeV", varcounter++, 'F'),
    //TMVA::VariableInfo("Mtjj := evalHMETMt(hJet_ptReg[0], hJet_phi[0], hJet_ptReg[1], hJet_phi[1])",
    //                   "mass_{T}(j1,j2)", "GeV", varcounter++, 'F'),
};

//varcounter = 0;
const TMVA::VariableInfo g_variables_KillTop[] = {
    TMVA::VariableInfo("TopHad_massReg",
                       "Top_{had} mass", "GeV", varcounter++, 'F'),
    TMVA::VariableInfo("j3ptfraction := Alt$(aJet_pt[0]/hJet_pt[0],0)",
                       "p_{T}(j3)/p_{T}(j1)", "", varcounter++, 'F')
};

//varcounter = 0;
const TMVA::VariableInfo g_variables_MET[] = {
    TMVA::VariableInfo("METtype1corr.et/sqrt(METtype1corr.sumet)",
                       "MET over sqrt(sumET)", "", varcounter++, 'F'),
    TMVA::VariableInfo("METtype1corr.sig",
                       "fancy MET significance", "", varcounter++, 'F'),
    TMVA::VariableInfo("METtype1corr.et/sqrt(HTeta45)",
                       "MET over sqrt(HT)", "", varcounter++, 'F')
};

//varcounter = 0;
const TMVA::VariableInfo g_supervariables[] = {
    //TMVA::VariableInfo("BDTregular_125[0]", "BDTregular", "", varcounter++, 'F'),
    TMVA::VariableInfo("BDTtt_125[0]", "BDTtt", "", varcounter++, 'F'),
    //TMVA::VariableInfo("BDTvj_125[0]", "BDTvj", "", varcounter++, 'F'),
    TMVA::VariableInfo("BDTvjlf_125[0]", "BDTvjlf", "", varcounter++, 'F'),
    //TMVA::VariableInfo("BDTvjhf_125[0]", "BDTvjhf", "", varcounter++, 'F'),
    //TMVA::VariableInfo("BDTvv_125[0]", "BDTvv", "", varcounter++, 'F'),
    TMVA::VariableInfo("BDTzz_125[0]", "BDTzz", "", varcounter++, 'F'),
    TMVA::VariableInfo("BDTzzhf_125[0]", "BDTzzhf", "", varcounter++, 'F'),
};

//const TMVA::VariableInfo g_supervariables[] = {
//    TMVA::VariableInfo("BDTregular_125[0]", "BDTregular", "", varcounter++, 'F'),
//    TMVA::VariableInfo("BDTregular_fj_125[0]", "BDTregular_fj", "", varcounter++, 'F'),
//    TMVA::VariableInfo("BDTregularhybrid_125[0]", "BDTregularhybrid", "", varcounter++, 'F'),
//    TMVA::VariableInfo("BDTregularhybrid_fj_125[0]", "BDTregularhybrid_fj", "", varcounter++, 'F'),
//};
////////////////////////////////////////////////////////////////////////////////
/// END Configuration                                                        ///
////////////////////////////////////////////////////////////////////////////////


///_____________________________________________________________________________
/// Functions

// Copied from tmvaglob.C
void NormalizeHists( TH1* sig, TH1* bkg )
{
    if (sig->GetSumw2N() == 0) sig->Sumw2();
    if (bkg->GetSumw2N() == 0) bkg->Sumw2();

    if (sig->GetSumOfWeights()!=0) {
        Float_t dx = (sig->GetXaxis()->GetXmax() - sig->GetXaxis()->GetXmin())/sig->GetNbinsX();
        sig->Scale( 1.0/sig->GetSumOfWeights()/dx );
    }
    if (bkg->GetSumOfWeights()!=0) {
        Float_t dx = (bkg->GetXaxis()->GetXmax() - bkg->GetXaxis()->GetXmin())/bkg->GetNbinsX();
        bkg->Scale( 1.0/bkg->GetSumOfWeights()/dx );
    }
}

// Copied from mvaeffs.C
TString GetFormula(TString formula)
{
   formula.ReplaceAll("S","x");
   formula.ReplaceAll("B","y");
   return formula;
}

///_____________________________________________________________________________
/// Core Function

void TrainChannelBDT(const TString& channel, int massH, TString whatbkg, TString whatsig, TString whatvar, const TString& outdir, int level, int NTrees, double MinNodeSize100 , int MaxDepth, int AdaBoostBeta100, int i_parlevel2, const TString& method)
{
    /// Parse options
    whatsig.ReplaceAll(" ", "");
    bool sigZH    = (whatsig == "") || (whatsig == "ZH");
    bool sigVH    = (whatsig == "VH");
    bool sigZZ    = (whatsig == "ZZ");
    bool sigZZHF  = (whatsig == "ZZHF");
    std::clog << "OPTION: sigZH=" << sigZH << ", sigVH=" << sigVH << ", sigZZ=" << sigZZ << ", sigZZHF=" << sigZZHF << std::endl;
    
    whatbkg.ReplaceAll(" ", "");
    bool bkgAll   = (whatbkg == "") || (whatbkg == "all") || (whatbkg == "ALL");
    bool bkgTT    = (whatbkg == "TT");
    bool bkgVj    = (whatbkg == "Vj");
    bool bkgVjLF  = (whatbkg == "VjLF");
    bool bkgZjLF  = (whatbkg == "ZjLF");
    bool bkgVjHF  = (whatbkg == "VjHF");
    bool bkgZjHF  = (whatbkg == "ZjHF");
    bool bkgST    = (whatbkg == "ST");
    bool bkgVV    = (whatbkg == "VV");
    bool bkgZZ    = (whatbkg == "ZZ");
    bool bkgZZHF  = (whatbkg == "ZZHF");
    std::clog << "OPTION: bkgAll=" << bkgAll << ", bkgTT=" << bkgTT << ", bkgVj=" << bkgVj << ", bkgVjLF=" << bkgVjLF << ", bkgZjLF=" << bkgZjLF << ", bkgVjHF=" << bkgVjHF << ", bkgZjHF=" << bkgZjHF << ", bkgST=" << bkgST << ", bkgVV=" << bkgVV << ", bkgZZ=" << bkgZZ << ", bkgZZHF=" << bkgZZHF << std::endl;
    
    whatvar.ReplaceAll(" ", "");
    bool useFJ         = (whatvar == "FJ");
    bool useFJReg      = (whatvar == "FJReg");
    bool useAngle      = (whatvar == "Angle");
    bool useKillTop    = (whatvar == "KillTop");
    bool useMET        = (whatvar == "MET");
    bool useNoMjj      = (whatvar == "NoMjj");
    bool useSuper      = (whatvar == "Super");
    std::clog << "OPTION: useFJ=" << useFJ << ", useFJReg=" << useFJReg << ", useAngle=" << useAngle << ", useKillTop=" << useKillTop << ", useMET=" << useMET << ", useNoMjj=" << useNoMjj << ", useSuper=" << useSuper << std::endl;
    
    /// Create a ROOT output file where TMVA will store ntuples, histograms, etc.
    TFile* out(0);
    if (level == 0) {
        out = TFile::Open( "TMVA_"+channel+".root", "RECREATE" );  // not kept
    } else {
        out = TFile::Open( outdir+"TMVA_"+channel+".root", "RECREATE" );
    }
    TString option = "";
    
    /// Setup TMVA Factory
    option = "!V:!Silent:Color:!DrawProgressBar:Transformations=I:AnalysisType=Classification";
    //option = "!V:!Silent:Color:!DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification";
    TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", out, option );
    
    /// Change directory and xml names
    (TMVA::gConfig().GetIONames()).fWeightFileDir       = outdir;
    (TMVA::gConfig().GetIONames()).fWeightFileExtension = "weights_" + channel;

    /// Add variables
    if (useNoMjj) {
        int nvars = sizeof(g_variables_NoMjj)/sizeof(g_variables_NoMjj[0]);
        for (int i=0; i<nvars; i++) {
            const TMVA::VariableInfo& var = g_variables_NoMjj[i];
            const TString labelexpression = (var.GetExpression() == var.GetLabel()) ? 
                                            var.GetExpression() : var.GetLabel() + ":=" + var.GetExpression();
            factory->AddVariable(labelexpression, var.GetTitle(), var.GetUnit(), var.GetVarType());
        }
    } else if (useSuper) {
        int nvars = sizeof(g_supervariables)/sizeof(g_supervariables[0]);
        for (int i=0; i<nvars; i++) {
            const TMVA::VariableInfo& var = g_supervariables[i];
            const TString labelexpression = (var.GetExpression() == var.GetLabel()) ? 
                                            var.GetExpression() : var.GetLabel() + ":=" + var.GetExpression();
            factory->AddVariable(labelexpression, var.GetTitle(), var.GetUnit(), var.GetVarType());
        }
    } else {
        int nvars = sizeof(g_variables)/sizeof(g_variables[0]);
        for (int i=0; i<nvars; i++) {
            const TMVA::VariableInfo& var = g_variables[i];
            const TString labelexpression = (var.GetExpression() == var.GetLabel()) ? 
                                            var.GetExpression() : var.GetLabel() + ":=" + var.GetExpression();
            factory->AddVariable(labelexpression, var.GetTitle(), var.GetUnit(), var.GetVarType());
        }
        
        if (useFJ) {
            nvars = sizeof(g_variables_FJ)/sizeof(g_variables_FJ[0]);
            for (int i=0; i<nvars; i++) {
                const TMVA::VariableInfo& var = g_variables_FJ[i];
                const TString labelexpression = (var.GetExpression() == var.GetLabel()) ? 
                                                var.GetExpression() : var.GetLabel() + ":=" + var.GetExpression();
                factory->AddVariable(labelexpression, var.GetTitle(), var.GetUnit(), var.GetVarType());
            }
        } else if (useFJReg) {
            nvars = sizeof(g_variables_FJReg)/sizeof(g_variables_FJReg[0]);
            for (int i=0; i<nvars; i++) {
                const TMVA::VariableInfo& var = g_variables_FJReg[i];
                const TString labelexpression = (var.GetExpression() == var.GetLabel()) ? 
                                                var.GetExpression() : var.GetLabel() + ":=" + var.GetExpression();
                factory->AddVariable(labelexpression, var.GetTitle(), var.GetUnit(), var.GetVarType());
            }
        } else if (useAngle) {
            nvars = sizeof(g_variables_Angle)/sizeof(g_variables_Angle[0]);
            for (int i=0; i<nvars; i++) {
                const TMVA::VariableInfo& var = g_variables_Angle[i];
                const TString labelexpression = (var.GetExpression() == var.GetLabel()) ? 
                                                var.GetExpression() : var.GetLabel() + ":=" + var.GetExpression();
                factory->AddVariable(labelexpression, var.GetTitle(), var.GetUnit(), var.GetVarType());
            }
        } else if (useKillTop) {
            nvars = sizeof(g_variables_KillTop)/sizeof(g_variables_KillTop[0]);
            for (int i=0; i<nvars; i++) {
                const TMVA::VariableInfo& var = g_variables_KillTop[i];
                const TString labelexpression = (var.GetExpression() == var.GetLabel()) ? 
                                                var.GetExpression() : var.GetLabel() + ":=" + var.GetExpression();
                factory->AddVariable(labelexpression, var.GetTitle(), var.GetUnit(), var.GetVarType());
            }
        } else if (useMET) {
            nvars = sizeof(g_variables_MET)/sizeof(g_variables_MET[0]);
            for (int i=0; i<nvars; i++) {
                const TMVA::VariableInfo& var = g_variables_MET[i];
                const TString labelexpression = (var.GetExpression() == var.GetLabel()) ? 
                                                var.GetExpression() : var.GetLabel() + ":=" + var.GetExpression();
                factory->AddVariable(labelexpression, var.GetTitle(), var.GetUnit(), var.GetVarType());
            }
        }
    }
    
    
    /// Get the trees from Step 4's
    TFile * ZH = TFile::Open(indir + prefix + Form("ZnnH%i", massH) + suffix);
    TFile * WH = TFile::Open(indir + prefix + Form("WlnH%i", massH) + suffix);
    TFile * Wj = TFile::Open(indir + prefix + "WJets" + suffix);
    TFile * Zj = TFile::Open(indir + prefix + "ZJets" + suffix);
    TFile * TT = TFile::Open(indir + prefix + "TTPythia8" + suffix);
    TFile * ST = TFile::Open(indir + prefix + "s_Top" + suffix);
    //    TFile * VV = TFile::Open(indir + prefix + "VV" + suffix); <-- missing
   
 //    TFile * ZZ = TFile::Open(indir + prefix + "ZZ" + suffix); <-- missing
    
    const TString treename_test = Form("tree_%s", channel.Data());
    const TString treename_train = Form("tree_%s_train", channel.Data());
    TTree * treeZH_test = (TTree *) ZH->Get(treename_test);
    TTree * treeWH_test = (TTree *) WH->Get(treename_test);
    TTree * treeWj_test = (TTree *) Wj->Get(treename_test);
    TTree * treeZj_test = (TTree *) Zj->Get(treename_test);
    TTree * treeTT_test = (TTree *) TT->Get(treename_test);
    TTree * treeST_test = (TTree *) ST->Get(treename_test);
    //    TTree * treeVV_test = (TTree *) VV->Get(treename_test); <--missing
    //    TTree * treeZZ_test = (TTree *) ZZ->Get(treename_test); <-- missing
    TTree * treeZH_train = (TTree *) ZH->Get(treename_train);
    TTree * treeWH_train = (TTree *) WH->Get(treename_train);
    TTree * treeWj_train = (TTree *) Wj->Get(treename_train);
    TTree * treeZj_train = (TTree *) Zj->Get(treename_train);
    TTree * treeTT_train = (TTree *) TT->Get(treename_train);
    TTree * treeST_train = (TTree *) ST->Get(treename_train);
    //    TTree * treeVV_train = (TTree *) VV->Get(treename_train); <-- missing
    //   TTree * treeZZ_train = (TTree *) ZZ->Get(treename_train);  <--missing
    
    /// Add signal trees
    if (sigZH) {
      std::cout << "before adding TRainging tree" << std::endl;
      factory->AddSignalTree(treeZH_train, 2.0, "Training");
      std::cout << "before adding test tree" << std::endl;
        factory->AddSignalTree(treeZH_test , 2.0, "Test"    );
        assert(treeZH_train->GetEntries() > 3);
    } else if (sigVH) {
      std::cout << "before adding TRainging tree" << std::endl;
        factory->AddSignalTree(treeZH_train, 2.0, "Training");
      std::cout << "before adding test tree" << std::endl;
        factory->AddSignalTree(treeZH_test , 2.0, "Test"    );
        factory->AddSignalTree(treeWH_train, 2.0, "Training");
        factory->AddSignalTree(treeWH_test , 2.0, "Test"    );
        assert(treeZH_train->GetEntries() > 3);
    } else if (sigZZ || sigZZHF) {
      //        factory->AddSignalTree(treeZZ_train, 2.0, "Training"); <-- missing
      //   factory->AddSignalTree(treeZZ_test , 2.0, "Test"    ); <--missing
      //   assert(treeZZ_train->GetEntries() > 3); <-- missing
    }
    
    /// Add background trees
    if (bkgAll) {
        factory->AddBackgroundTree(treeWj_train, 2.0, "Training");
      std::cout << "before adding test tree for Wj" << std::endl;
        factory->AddBackgroundTree(treeWj_test , 2.0, "Test"    );
        factory->AddBackgroundTree(treeZj_train, 2.0, "Training");
      std::cout << "before adding test tree for Zj" << std::endl;
        factory->AddBackgroundTree(treeZj_test , 2.0, "Test"    );
        factory->AddBackgroundTree(treeTT_train, 2.0, "Training");
      std::cout << "before adding test tree for TT" << std::endl;
        factory->AddBackgroundTree(treeTT_test , 2.0, "Test"    );
        assert(treeWj_train->GetEntries() > 3);
        assert(treeZj_train->GetEntries() > 3);
        assert(treeTT_train->GetEntries() > 3);
        
        if (g_fullbkg) {
            factory->AddBackgroundTree(treeST_train, 2.0, "Training");
	    std::cout << "before adding test tree for ST" << std::endl;
            factory->AddBackgroundTree(treeST_test , 2.0, "Test"    );            
            assert(treeST_train->GetEntries() > 3);
            if (!sigZZ && !sigZZHF) {
	      //                factory->AddBackgroundTree(treeVV_train, 2.0, "Training"); <-- missing
              //  factory->AddBackgroundTree(treeVV_test , 2.0, "Test"    ); <-- missing
	      //   assert(treeVV_train->GetEntries() > 3); <-- missing
            }
        }
    } else if (bkgTT) {
        factory->AddBackgroundTree(treeTT_train, 2.0, "Training");
        factory->AddBackgroundTree(treeTT_test , 2.0, "Test"    );
        assert(treeTT_train->GetEntries() > 3);
    } else if (bkgVj || bkgVjLF || bkgVjHF) {
        factory->AddBackgroundTree(treeWj_train, 2.0, "Training");
        factory->AddBackgroundTree(treeWj_test , 2.0, "Test"    );
        factory->AddBackgroundTree(treeZj_train, 2.0, "Training");
        factory->AddBackgroundTree(treeZj_test , 2.0, "Test"    );
        assert(treeWj_train->GetEntries() > 3);
        assert(treeZj_train->GetEntries() > 3);
    } else if (bkgZjLF || bkgZjHF) {
        factory->AddBackgroundTree(treeZj_train, 2.0, "Training");
        factory->AddBackgroundTree(treeZj_test , 2.0, "Test"    );
        assert(treeZj_train->GetEntries() > 3);
    } else if (bkgST) {
        factory->AddBackgroundTree(treeST_train, 2.0, "Training");
        factory->AddBackgroundTree(treeST_test , 2.0, "Test"    );
        assert(treeST_train->GetEntries() > 3);
    } else if (bkgVV) {
      //        factory->AddBackgroundTree(treeVV_train, 2.0, "Training"); <-- missing
      //  factory->AddBackgroundTree(treeVV_test , 2.0, "Test"    );<-- missing
      //  assert(treeVV_train->GetEntries() > 3);<-- missing
    } else if (bkgZZ || bkgZZHF) {
      //        factory->AddBackgroundTree(treeZZ_train, 2.0, "Training"); <-- missing
      //  factory->AddBackgroundTree(treeZZ_test , 2.0, "Test"    ); <-- missing
      //  assert(treeZZ_train->GetEntries() > 3); <-- missing
    }
    
    /// Set individual event weights (the variables must exist in the original TTree)
    //  factory->SetSignalWeightExpression    ("weightsMC[0]*(selectFlags[0][0])");
  factory->SetSignalWeightExpression    ("efflumi");
    // factory->SetBackgroundWeightExpression("weightsMC[0]*(selectFlags[0][0])");
factory->SetBackgroundWeightExpression("efflumi");
    //if (bkgAll || bkgVj || bkgVjLF || bkgVjHF || bkgZjLF || bkgZjHF) {
    //    factory->SetBackgroundWeightExpression("weightsMC[0]*(selectFlags[0][0])*max(0.5, (1.0 + (-0.0025)*(max(genZ.pt,130) - 130)) )");
    //} else {
    //    factory->SetBackgroundWeightExpression("weightsMC[0]*(selectFlags[0][0])");
    //}

    /// Apply additional cuts on the signal and background samples (can be different)
    TCut cuts = "";
    TCut cutb = "";
    if (sigZZHF)
        cuts = "(abs(Jet_mcFlavour[hJCidx[0]])==5 || abs(Jet_mcFlavour[hJCidx[1]])==5)";
    if (bkgVjLF || bkgZjLF)
        cutb = "!(abs(Jet_mcFlavour[hJCidx[0]])==5 || abs(Jet_mcFlavour[hJCidx[1]])==5)";
    if (bkgVjHF || bkgZjHF)
        cutb = "(abs(Jet_mcFlavour[hJCidx[0]])==5 || abs(Jet_mcFlavour[hJCidx[1]])==5)";
    if (bkgZZHF)
        cutb = "(abs(Jet_mcFlavour[hJCidx[0]])==5 || abs(Jet_mcFlavour[hJCidx[1]])==5)";
    if (bkgVV && (sigZZ || sigZZHF))
        cutb = "processname!=\"ZZ\"";

    if (g_splitsig) {
        /// instead of i % 2: (0,1) = (train, test),
        /// do i % 4: (0,1,2,3) = (reg train, reg test, spec train, spec test)
        if (sigZZ || sigZZHF) {
            std::cerr << "WARNING: you are training with ZZ or ZZHF as signal, and you want to split signal?" << std::endl;
        }
        if (bkgAll) {
            cuts += "(EVENT.event %4 == 0) || (EVENT.event %4 == 1)";
            factory->SetSignalWeightExpression("2.0*weightsMC[0]*(selectFlags[0][0])");
        } else if (bkgTT || bkgVj || bkgVjLF || bkgZjLF || bkgVjHF || bkgZjHF || bkgST || 
                   bkgVV || bkgZZ || bkgZZHF) {
            cuts += "(EVENT.event %4 == 2) || (EVENT.event %4 == 3)";
            factory->SetSignalWeightExpression("2.0*weightsMC[0]*(selectFlags[0][0])");
        }
    }
    if (g_hybrid) {
        /// split events either with FJ info or not
        if (useFJ || useFJReg) {
            cuts += "(nfathFilterJets>0 && FatH.FatHiggsFlag==1)";
            cutb += "(nfathFilterJets>0 && FatH.FatHiggsFlag==1)";
        } else {
            /// FIXME: remove the cuts to gain more statistical power?
            cuts += "(!(nfathFilterJets>0 && FatH.FatHiggsFlag==1))";
            cutb += "(!(nfathFilterJets>0 && FatH.FatHiggsFlag==1))";
        }
    }


    /// Tell the factory to use all remaining events in the trees after training for testing
    option = "!V:nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=None";
    factory->PrepareTrainingAndTestTree( cuts, cutb, option );

    /// Book MVA methods
    std::vector<TString> parameters;
    if (level == 0) {  // run only 1 trainig
      option = Form("!H:V:NTrees=%i:MinNodeSize=%.3f:MaxDepth=%i:BoostType=AdaBoost:AdaBoostBeta=%.2f:SeparationType=GiniIndex:nCuts=35:PruneMethod=NoPruning:!CreateMVAPdfs:!DoBoostMonitor", NTrees, float(MinNodeSize100)/100., MaxDepth, float(AdaBoostBeta100)/100.);
        //option = "!H:!V:NTrees=300:MinNodeSize=375:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.1:SeparationType=GiniIndex:nCuts=35:PruneMethod=NoPruning" );  // !superBDT
        //option =  "!H:V:NTrees=100:MinNodeSize=100:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.1:SeparationType=GiniIndex:nCuts=35:PruneMethod=NoPruning" );  // superBDT
        factory->BookMethod( TMVA::Types::kBDT, method, option );
        parameters.push_back(option);
        std::cout << "OPTION: " << option << std::endl;
    
    /// Optimization
    } else if (level == 1) {  // optimize the parameters for a given NTrees (5x2x1 steps)
        for (int j=0; j<5; j++) {
            int MinNodeSize_ = 150 + j*50;  // 150 - 350
            for (int k=0; k<2; k++) {
                int MaxDepth_ = 3 + k;  // 3 - 4
                for (int l=0; l<1; l++) {
                    int AdaBoostBeta100_ = 10 + l*10;  // 10 - 10
                    option = Form("!H:!V:NTrees=%i:MinNodeSize=%i:MaxDepth=%i:BoostType=AdaBoost:AdaBoostBeta=%.2f:SeparationType=GiniIndex:nCuts=35:NNodesMax=100000:PruneMethod=NoPruning", NTrees, MinNodeSize_, MaxDepth_, float(AdaBoostBeta100_)/100.);
                    factory->BookMethod( TMVA::Types::kBDT, method+Form("%lu",parameters.size()), option );
                    parameters.push_back(option);
                    std::cout << "OPTION: " << option << std::endl;
                }
            }
        }
    
    } else if (level == 2) {  // optimize the parameters for given NTrees and another parameter (5x2x1 steps)
        const char * pars_AdaBoostBeta[] = {  // 5
            "0.01", "0.1", "0.2", "0.3", "0.5",
        };
        const char * pars_SeparationType[] = {  // 3 + 2 dummies
            "GiniIndex", "GiniIndexWithLaplace", "MisClassificationError", "GiniIndex", "GiniIndex",
        };
        const char * pars_nCuts[] = {  // 5
            "-1", "20", "30", "35", "40",
        };
        const char * pars_NNodesMax[] = {  // 5
            "7", "9", "11", "13", "15",
        };
        assert(i_parlevel2 < 5);
        
        for (int j=0; j<5; j++) {
            int MinNodeSize_ = 150 + j*50;  // 150 - 350
            for (int k=0; k<2; k++) {
                int MaxDepth_ = 3 + k;  // 3 - 4
                for (int l=0; l<1; l++) {
                    int AdaBoostBeta100_ = 10 + l*10;  // 10 - 10
                    
                    /// For AdaBoostBeta (note: AdaBoostBeta100_ is ignored)
                    option = Form("!H:!V:NTrees=%i:MinNodeSize=%i:MaxDepth=%i:BoostType=AdaBoost:AdaBoostBeta=%s:SeparationType=GiniIndex:nCuts=35:NNodesMax=100000:PruneMethod=NoPruning", NTrees, MinNodeSize_, MaxDepth_, pars_AdaBoostBeta[i_parlevel2]);
                    
                    /// For SeparationType
                    //option = Form("!H:!V:NTrees=%i:MinNodeSize=%i:MaxDepth=%i:BoostType=AdaBoost:AdaBoostBeta=%.2f:SeparationType=%s:nCuts=35:NNodesMax=100000:PruneMethod=NoPruning", NTrees, MinNodeSize_, MaxDepth_, float(AdaBoostBeta100_)/100., pars_SeparationType[i_parlevel2]);
                
                    /// For nCuts
                    //option = Form("!H:!V:NTrees=%i:MinNodeSize=%i:MaxDepth=%i:BoostType=AdaBoost:AdaBoostBeta=%.2f:SeparationType=GiniIndex:nCuts=%s:NNodesMax=100000:PruneMethod=NoPruning", NTrees, MinNodeSize_, MaxDepth_, float(AdaBoostBeta100_)/100., pars_nCuts[i_parlevel2]);
                    
                    /// For NNodexMax (note: set MaxDepth to be 5)
                    //option = Form("!H:!V:NTrees=%i:MinNodeSize=%i:MaxDepth=%i:BoostType=AdaBoost:AdaBoostBeta=%.2f:SeparationType=GiniIndex:nCuts=35:NNodesMax=%s:PruneMethod=NoPruning", NTrees, MinNodeSize_, 5, float(AdaBoostBeta100_)/100., pars_NNodesMax[i_parlevel2]);
                
                    factory->BookMethod( TMVA::Types::kBDT, method+Form("%lu",parameters.size()), option );
                    parameters.push_back(option);
                    std::cout << "OPTION: " << option << std::endl;
                }
            }
        }
    }  // end if, else if (level)

    /// Train MVAs using the set of training events
    factory->TrainAllMethods();

    /// Evaluate all MVAs using the set of test events
    factory->TestAllMethods();

    /// Evaluate and compare performance of all configured MVAs
    factory->EvaluateAllMethods();

    /// Write FOM (not kept for level==0
    if (level != 0)  {
        Double_t sig, sep, roc, eff01, eff10, eff30, effArea;
        Double_t maxsignif=-1., maxpsignif=-1.;
        Double_t maxsignif_d=-1., maxpsignif_d=-1.;
        Double_t kolS, kolB;
        
        out->cd();
        TTree * fomtree = new TTree("fomtree", "fomtree");
        fomtree->Branch("NTrees", &NTrees);
        fomtree->Branch("MinNodeSize100", &MinNodeSize100);
        fomtree->Branch("MaxDepth", &MaxDepth);
        fomtree->Branch("AdaBoostBeta100", &AdaBoostBeta100);
        fomtree->Branch("i_parlevel2", &i_parlevel2);
        fomtree->Branch("TMVA_sig", &sig);
        fomtree->Branch("TMVA_sep", &sep);
        fomtree->Branch("TMVA_roc", &roc);
        fomtree->Branch("TMVA_eff01", &eff01);
        fomtree->Branch("TMVA_eff10", &eff10);
        fomtree->Branch("TMVA_eff30", &eff30);
        fomtree->Branch("TMVA_effArea", &effArea);
        fomtree->Branch("USER_sig", &maxsignif);
        fomtree->Branch("USER_psig", &maxpsignif);
        fomtree->Branch("USER_sig_d", &maxsignif_d);
        fomtree->Branch("USER_psig_d", &maxpsignif_d);
        fomtree->Branch("TMVA_kolS", &kolS);
        fomtree->Branch("TMVA_kolB", &kolB);
        
        for (UInt_t i=0; i<parameters.size(); i++) {
            const TString& par = parameters.at(i);
            TPRegexp regexp("NTrees=([0-9]+):MinNodeSize=([0-9]+):MaxDepth=([0-9]+):BoostType=AdaBoost:AdaBoostBeta=(0.[0-9]+):");
            TObjArray * matches = regexp.MatchS(par);
            assert(matches->GetLast() == 4);
            NTrees = ((TObjString *) matches->At(1))->GetString().Atoi();
            MinNodeSize100 = ((TObjString *) matches->At(2))->GetString().Atoi();
            MaxDepth = ((TObjString *) matches->At(3))->GetString().Atoi();
            Double_t AdaBoostBeta = ((TObjString *) matches->At(4))->GetString().Atof();
            AdaBoostBeta100 = int(AdaBoostBeta * 100.);
            std::cout << "NTrees=" << NTrees << ", MinNodeSize=" << MinNodeSize100 << ", MaxDepth=" << MaxDepth << ", AdaBoostBeta100=" << AdaBoostBeta100 << std::endl;
            
            const TString methodTitle = (level != 0) ? method+Form("%u",i) : method;
            TMVA::MethodBDT * theMethod = dynamic_cast<TMVA::MethodBDT *>(factory->GetMethod(methodTitle));
            
            Double_t err;
            sig = theMethod->GetSignificance();
            sep = theMethod->GetSeparation();
            roc = theMethod->GetROCIntegral();
            eff01 = theMethod->GetEfficiency("Efficiency:0.01", TMVA::Types::kTesting, err);
            //eff01err = err;
            eff10 = theMethod->GetEfficiency("Efficiency:0.10", TMVA::Types::kTesting, err);
            //eff10err = err;
            eff30 = theMethod->GetEfficiency("Efficiency:0.30", TMVA::Types::kTesting, err);
            //eff30err = err;
            effArea = theMethod->GetEfficiency("",              TMVA::Types::kTesting, err);
            
            std::cout << "sig=" << sig << ", sep=" << sep << ", roc=" << roc << std::endl;
            std::cout << "eff01=" << eff01 << ", eff10=" << eff10 << ", eff30=" << eff30 << ", effArea=" <<effArea << std::endl;
            
            TString hname = "Method_BDT/" + methodTitle + "/MVA_" + methodTitle;
            TH1* hsig   = dynamic_cast<TH1*>(out->Get( hname + "_S" ));  // nbins = 40
            TH1* hbgd   = dynamic_cast<TH1*>(out->Get( hname + "_B" ));  // nbins = 40
            TH1* hsigOv = dynamic_cast<TH1*>(out->Get( hname + "_Train_S" ));  // nbins = 40
            TH1* hbgdOv = dynamic_cast<TH1*>(out->Get( hname + "_Train_B" ));  // nbins = 40
            TH1* hsigE  = dynamic_cast<TH1*>(out->Get( hname + "_effS" ));  // nbins = 10000
            TH1* hbgdE  = dynamic_cast<TH1*>(out->Get( hname + "_effB" ));  // nbins = 10000
            
            TString hnameNm = "Method_BDT/" + methodTitle + "/" + g_variables_NoMjj[0].GetLabel() + "_";  // chose NoMjj as it's safer
            TH1* hsigNm = dynamic_cast<TH1*>(out->Get( hnameNm + "_Signal" ));  // from test trees
            TH1* hbkgNm = dynamic_cast<TH1*>(out->Get( hnameNm + "_Background" ));  // from test trees
            const Double_t nsig=hsigNm->Integral(), nbgd=hbkgNm->Integral();
            
            // Significance
            TFormula fsignif("fsignif", GetFormula("S/sqrt(S+B)"));
            TFormula fpsignif("fpsignif", GetFormula("S/(1.5+sqrt(B)+0.2*B)"));  // for 3 sigmas
            
            for (Int_t ibin=1; ibin<=hsigE->GetNbinsX(); ibin++) {
                Float_t S = hsigE->GetBinContent( ibin ) * nsig;
                Float_t B = hbgdE->GetBinContent( ibin ) * nbgd;
                Double_t signif = fsignif.Eval(S,B);
                Double_t psignif = fpsignif.Eval(S,B);
                if (maxsignif < signif) {
                    maxsignif = signif;
                    maxsignif_d = hsigE->GetXaxis()->GetBinCenter(ibin);
                }
                if (maxpsignif < psignif) {
                    maxpsignif = psignif;
                    maxpsignif_d = hsigE->GetXaxis()->GetBinCenter(ibin);
                }
            }
            
            // K-S test
            NormalizeHists(hsig, hbgd);
            NormalizeHists(hsigOv, hbgdOv);
            kolS = hsig->KolmogorovTest( hsigOv );
            kolB = hbgd->KolmogorovTest( hbgdOv );
            
            std::cout << "maxsignif=" << maxsignif << ", maxpsignif=" << maxpsignif << ", kolS=" << kolS << ", kolB=" << kolB << std::endl;
            std::cout << std::endl;
            
            fomtree->Fill();
        }
        fomtree->Write();
    }  // end if level != 0
    

    /// Clean up
    out->Close();
    ZH->Close();
    WH->Close();
    Wj->Close();
    Zj->Close();
    TT->Close();
    ST->Close();
    //    VV->Close(); <--missing
    //    ZZ->Close(); <-- missing

    delete factory;
    delete out;
    delete ZH;
    delete WH;
    delete Wj;
    delete Zj;
    delete TT;
    delete ST;
    //    delete VV; <-- missing
    //    delete ZZ; <-- missing
    
    return;
}

////////////////////////////////////////////////////////////////////////////////
/// Main                                                                     ///
////////////////////////////////////////////////////////////////////////////////
void TrainBDT(int massH=125, TString whatbkg="", TString whatsig="", TString whatvar="", TString outdir="weightsTEST/", int level=0, int i_NTrees=100000, int i_parlevel2=0, const TString method="BDT")
{
    // level==0 takes 0 additional parameter, use preset BDT parameters
    // level==1 takes 1 additional parameter: NTrees
    // level==2 takes 2 additional parameter: NTrees and another parameter depending what is enabled under Optimization

    gROOT->SetBatch(1);
    gROOT->LoadMacro("HelperFunctions.h" );  // make functions visible to TTreeFormula

    /// Make the output directory
    if (!outdir.EndsWith("/"))
        outdir += "/";
    if (gSystem->AccessPathName(outdir))
        gSystem->mkdir(outdir);
    outdir += Form("%i/", massH);
    if (gSystem->AccessPathName(outdir))
        gSystem->mkdir(outdir);
    
    assert(g_splitsig+g_fullbkg+g_hybrid <= 1);  // only 1 can be true at a time

    if (!TString(gROOT->GetVersion()).Contains("5.34")) {
        std::cout << "INCORRECT ROOT VERSION! Please use 5.34:" << std::endl;
        std::cout << "source /uscmst1/prod/sw/cms/slc5_amd64_gcc462/lcg/root/5.34.02-cms/bin/thisroot.csh" << std::endl;
        std::cout << "Return without doing anything." << std::endl;
        return;
    }    

    /// Load the library
    TMVA::Tools::Instance();
    int NTrees = 400, MinNodeSize100 = 0.05, MaxDepth = 3, AdaBoostBeta100 = 50;  // TMVA default
    TString channel = "ZnunuHighPt";
    
    if (level == 0) {
        ///-- ZnunuHighPt ------------------------------------------------------
        channel         = "ZnunuHighPt";
        NTrees          = pars_ZnunuHighPt[0];
        MinNodeSize100      = pars_ZnunuHighPt[1];
        MaxDepth        = pars_ZnunuHighPt[2];
        AdaBoostBeta100 = pars_ZnunuHighPt[3];
        TrainChannelBDT(channel, massH, whatbkg, whatsig, whatvar, outdir, level, NTrees, MinNodeSize100, MaxDepth, AdaBoostBeta100, 0, method);
        
        ///-- ZnunuMedPt -------------------------------------------------------
        channel         = "ZnunuMedPt";
        NTrees          = pars_ZnunuMedPt[0];
        MinNodeSize100      = pars_ZnunuMedPt[1];
        MaxDepth        = pars_ZnunuMedPt[2];
        AdaBoostBeta100 = pars_ZnunuMedPt[3];
        TrainChannelBDT(channel, massH, whatbkg, whatsig, whatvar, outdir, level, NTrees, MinNodeSize100, MaxDepth, AdaBoostBeta100, 0, method);

        ///-- ZnunuLowPt -------------------------------------------------------
        channel         = "ZnunuLowPt";
        NTrees          = pars_ZnunuLowPt[0];
        MinNodeSize100      = pars_ZnunuLowPt[1];
        MaxDepth        = pars_ZnunuLowPt[2];
        AdaBoostBeta100 = pars_ZnunuLowPt[3];
        TrainChannelBDT(channel, massH, whatbkg, whatsig, whatvar, outdir, level, NTrees, MinNodeSize100, MaxDepth, AdaBoostBeta100, 0, method);
/*
        ///-- ZnunuLowCSV ------------------------------------------------------
        channel         = "ZnunuLowCSV";
        NTrees          = pars_ZnunuLowCSV[0];
        MinNodeSize      = pars_ZnunuLowCSV[1];
        MaxDepth        = pars_ZnunuLowCSV[2];
        AdaBoostBeta100 = pars_ZnunuLowCSV[3];
        TrainChannelBDT(channel, massH, whatbkg, whatsig, whatvar, outdir, level, NTrees, MinNodeSize, MaxDepth, AdaBoostBeta100, 0, method);
*/
    } else if (level == 1 || level == 2) {
        ///-- ZnunuHighPt ------------------------------------------------------
        channel     = "ZnunuHighPt";
        TrainChannelBDT(channel, massH, whatbkg, whatsig, whatvar, outdir, level, i_NTrees, MinNodeSize100, MaxDepth, AdaBoostBeta100, i_parlevel2, method);
        
        ///-- ZnunuMedPt -------------------------------------------------------
        channel     = "ZnunuMedPt";
        TrainChannelBDT(channel, massH, whatbkg, whatsig, whatvar, outdir, level, i_NTrees, MinNodeSize100, MaxDepth, AdaBoostBeta100, i_parlevel2, method);
        
        ///-- ZnunuLowPt -------------------------------------------------------
        channel     = "ZnunuLowPt";
        TrainChannelBDT(channel, massH, whatbkg, whatsig, whatvar, outdir, level, i_NTrees, MinNodeSize100, MaxDepth, AdaBoostBeta100, i_parlevel2, method);
/*
        ///-- ZnunuLowCSV ------------------------------------------------------
        channel     = "ZnunuLowCSV";
        TrainChannelBDT(channel, massH, whatbkg, whatsig, whatvar, outdir, level, i_NTrees, MinNodeSize, MaxDepth, AdaBoostBeta100, i_parlevel2, method);
*/
    }
}

#if 0
    // Launch the GUI for the root macros
    TString TMVATEST = "$ROOTSYS/tmva/test/";
    gROOT->SetMacroPath( TMVATEST );
    gROOT->Macro( TMVATEST+"TMVAlogon.C" );
    gROOT->LoadMacro( TMVATEST+"TMVAGui.C" );
    TMVAGui( "TMVA_ZnunuHighPt.root" );
#endif
