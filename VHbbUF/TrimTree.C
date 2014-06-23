#include <cassert>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TFile.h"
#include "TString.h"
#include "TBranch.h"
#include "TBranchElement.h"
#include "TCut.h"
#include "TMath.h"
#include "TPRegexp.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TTree.h"
#include "TTreeFormula.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
//#include "TMVA/TXMLEngine.h"
#endif

#include "HelperTMVA.h"
#include "HelperFunctions.h"
//#include "HelperNtuples.h"
//#include "HelperVHbbDataFormats.h"

#define MAXSYST 20
#define MAXSELECT 20
#define MAXVAR 100
#define MAXWEIGHT 240

//#define STITCH

#define SPECIALQCD


////////////////////////////////////////////////////////////////////////////////
/// Configuration                                                            ///
////////////////////////////////////////////////////////////////////////////////

//const TString indir       = "dcache:/pnfs/cms/WAX/resilient/jiafu//ZnunuHbb/Step3_20130314/";
const TString indir       = "/eos/uscms/store/user/jiafu/ZnunuHbb/ZnunuHbb/Step3_20130314/";
const TString outdir      = "skim/";
const TString prefix      = "Step3_";
const TString suffix      = ".root";
const TString indirStep4  = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130404/stitch/";
//const TString indirStep4  = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130326/stitch/";
//const TString indirStep4  = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130314/stitch/";
//const TString indirStep4  = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130302/stitch/";
//const TString indirStep4  = "/uscms_data/d3/lpchbb/jiafu/ZnnH_postHCP/Step4_20130221/stitch/";
//TString weightdir         = "./";
TString weightdir         = "weights_20130401/";
//TString weightdir         = "weights_20130327_V4a/";
//TString weightdir         = "weights_20130323/";
//TString weightdir         = "weights_20130302/";

/// Format is ($WEIGHTDIR)/($MASS)/TMVAClassification_($TMVAMETHOD).weights_($CHANNEL).xml
/// Expect, for each dir, there are (i mass points) x (j channels)
/// e.g. weightsTT/125/TMVAClassification_BDT.weights_ZnunuHighPt.xml

const char* g_weights[] = {
    "weightsRegular",
    "weightsTT",
    //"weightsVj",
    "weightsVjLF",
    //"weightsZjLF",
    "weightsVjHF",
    //"weightsZjHF",
    //"weightsST",
    //"weightsVV",
    "weightsZZ",
    //"weightsZZHF",

    //"weightsRegular_FJ",
    //"weightsTT_FJ",
    //"weightsVjLF_FJ",
    //"weightsZZ_FJ",
    //"weightsZZHF_FJ",

    "weightsRegular_FJReg",
    //"weightsTT_FJReg",
    //"weightsVjLF_FJReg",
    //"weightsZZ_FJReg",
    //"weightsZZHF_FJReg",

    //"weightsRegular_Angle",

    //"weightsRegular_KillTop",

    //"weightsRegular_MET",

    "weightsRegular_NoMjj",

    "weightsRegular__sigZZHF",
    "weightsTT__sigZZHF",
    "weightsVjLF__sigZZHF",
    "weightsVV__sigZZHF",

    //"weightsRegularSplitSig",
    //"weightsTTSplitSig",
    //"weightsVjLFSplitSig",
    //"weightsZZSplitSig",
    //"weightsZZHFSplitSig",

    //"weightsRegularFullBkg",

    //"weightsRegularHybrid",
    //"weightsRegularHybrid_FJ",
    //"weightsRegularHybrid_FJReg",
};

const int g_masses[] = {
    //110,
    //115,
    //120,
    125,
    //130,
    //135,
    //140,
    //145,
    //150,
};

/// Search bins
const char* g_channels[] = {
    "ZnunuHighPt",
    "ZnunuMedPt",
    "ZnunuLowPt",
    //"ZnunuLowCSV",
};

/// Variables to keep
const char* g_variables[] = {
    "Vtype", "EVENT*", "event*", "lhe*", "nPVs",
    "lumi", "efflumi*", "PUweight*",
    "H*", "Hmass*", "Hpt*", "V*", "MET*",
    "nhJets", "naJets", "nvlep", "nalep",
    "hJet_pt*", "hJet_eta*", "hJet_phi*", "hJet_e*", "hJet_csv*", "hJet_id", "hJet_puJetIdL", "hJet_genPt", "hJet_flavour", "hJet_nhf", "hJet_nef",
    "hJet_ptRaw", "hJet_JECUnc", "hJet_vtxMass", "hJet_vtxPt", "hJet_vtx3dL", "hJet_vtx3deL", "hJet_chf", "hJet_cef", "hJet_nconstituents", "hJet_nch",
    "hJet_SoftLeptPt", "hJet_SoftLeptdR", "hJet_SoftLeptptRel", "hJet_SoftLeptIdlooseMu", "hJet_SoftLeptId95",
    "aJet_pt*", "aJet_eta*", "aJet_phi*", "aJet_e*", "aJet_csv*", "aJet_id", "aJet_puJetIdL", "aJet_genPt", "aJet_flavour",
    "FatH", "FatHmass*", "FatHpt*",
    "nfathFilterJets", "naJetsFat",
    "fathFilterJets_pt*", "fathFilterJets_eta*", "fathFilterJets_phi*", "fathFilterJets_e*", "fathFilterJets_csv*", "fathFilterJets_genPt", "fathFilterJets_flavour", "fathFilterJets_vtxMass",
    "aJetFat_pt*", "aJetFat_eta*", "aJetFat_phi*", "aJetFat_e*", "aJetFat_csv*",
    "mindPhiMET*", "mindRH*", "naJets_Znn*", "nalep_Znn", "nalep_pt5_Znn",
    "deltaPullAngle", "rho*",
    "triggerFlags", "hbhe", "ecalFlag", "cschaloFlag", "hcallaserFlag","trackingfailureFlag", "eebadscFlag", "isBadHcalEvent",
    "TopLep*", "TopHad*", "FatTopHad*",
    "nPdf", "PDFweight", "weightMCProd", "weightSignalEWK", "weightSignalQCD", "VtypeWithTau", "lheV_pt", "genZ", "genW", "genH",

    //"weightTrigMu", "weightTrigDLP*", "weightTrig2012*", "weightTrig",
    //"vLepton_pt", "vLepton_eta", "vLepton_phi", "vLepton_pfCombRelIso", "vLepton_vbtf", "vLepton_id95",
    //"aLepton_pt", "aLepton_eta", "aLepton_phi", "aLepton_pfCombRelIso", "aLepton_vbtf", "aLepton_id95",
};
////////////////////////////////////////////////////////////////////////////////
/// END Configuration                                                        ///
////////////////////////////////////////////////////////////////////////////////

///_____________________________________________________________________________
/// Classes & Functions

// Copied from TMVA::VariableInfo
std::vector<TString> ParseExpression(const TString& expression)
{
    std::vector<TString> expressions;
    expressions.push_back("");
    expressions.push_back(expression);
    if (expression.Contains(":=")) {
        Ssiz_t index = expression.Index(":=");
        expressions[0] = expression(0, index);
        expressions[1] = expression(index + 2, expression.Sizeof() - index - 2);
        expressions[0] = expressions[0].Strip(TString::kBoth);
        expressions[1] = expressions[1].Strip(TString::kBoth);
    }
    return expressions;
}

class TTreeFormulaWithSyst {
public:
    TTreeFormulaWithSyst(const char* name, const char* formula, TTree* tree,
        const std::vector<std::vector<std::pair<std::string, std::string> > > & systematics, bool b)
      : formulas(), isData(b)
    {
        // 0 is reserved for NONE
        const UInt_t nsyst = systematics.size();
        assert(nsyst > 0);
        assert(systematics.front().size() == 1);
        assert(systematics.front().front().first == systematics.front().front().second);

        TTreeFormula* ttf = 0;
        TString expression(formula);
        if (expression.Contains(":="))
            expression = ParseExpression(expression)[1];  // strip the label
        Bool_t verbose = false;
        for (UInt_t isyst = 0; isyst < nsyst; isyst++) {
            TString systexpression = expression;
            if (!isData) {
                const std::vector<std::pair<std::string, std::string> >& si = systematics.at(isyst);
                for (UInt_t jsyst = 0; jsyst < si.size(); jsyst++) {
                    const TString& str = si.at(jsyst).first;
                    const TString& repl = si.at(jsyst).second;
                    if (expression.Contains(str)){
                        TPRegexp regexp("(\\A|\\W)"+str+"(\\Z|\\W)");
                        regexp.Substitute(systexpression, "$1"+repl+"$2", "g");  // "g" for global
                    }
                }
            }
            ttf = new TTreeFormula(Form("%s_%i", name, isyst), systexpression, tree);
            ttf->SetQuickLoad(1);
            if (!ttf || !ttf->GetNdim()) {
                std::cerr << "ERROR: Failed to find any TTree variable: " << expression << std::endl;
            }
            formulas.push_back(ttf);
            changed.push_back(Bool_t(expression != systexpression));
            if (verbose && (expression != systexpression))
                std::cout << "VERBOSE: isyst " << isyst << " " << expression << " --> " << systexpression << std::endl;
        }  // end loop over systematics
    }

    ~TTreeFormulaWithSyst()
    {
        for (UInt_t i = 0; i < size(); i++)
            delete formulas.at(i);
    }

    void UpdateFormulaLeaves() {
        for (UInt_t i = 0; i < size(); i++)
            formulas.at(i)->UpdateFormulaLeaves();
        return;
    }

    void GetNdata() {
        for (UInt_t i = 0; i < size(); i++)
            formulas.at(i)->GetNdata();
        return;
    }

    Double_t EvalInstance(UInt_t i) const {
        return formulas.at(i)->EvalInstance();
    }

    Bool_t EvalBoolInstance(UInt_t i) const {
        // from TTreePlayer::CopyTree(...)
        const Int_t ndata = formulas.at(i)->GetNdata();
        Bool_t keep = kFALSE;
        for (Int_t current = 0; current<ndata && !keep; current++) {
            keep |= (formulas.at(i)->EvalInstance(current) != 0);
        }
        return keep;
    }

    void Is_Changed(std::vector<Bool_t>& in) {
        assert(in.size() == changed.size());
        for (UInt_t i=0; i < in.size(); i++) {
            in.at(i) = in.at(i) || changed.at(i);
        }
        return;
    }

    UInt_t size() const {return formulas.size();}

private:
    std::vector<TTreeFormula *> formulas;
    std::vector<Bool_t> changed;
    Bool_t isData;
};


std::vector<TString> GetVariablesFromXML(const TString& xml)
{
    void* doc = TMVA::gTools().xmlengine().ParseFile(xml, 1000000);  // the default buffer size in TXMLEngine::ParseFile is 100k. changed to 1000k.
    void* rootnode = TMVA::gTools().xmlengine().DocGetRootElement(doc);
    std::vector<TString> expressions;

    bool verbose = true;
    TString nodeName = "";
    UInt_t readNVar;
    void* ch = TMVA::gTools().GetChild(rootnode);
    while (ch!=0) {
        nodeName = TString( TMVA::gTools().GetName(ch) );
        if (nodeName == "Variables") {
            TMVA::gTools().ReadAttr(ch, "NVar", readNVar);
            TMVA::VariableInfo readVarInfo;
            int varIdx = 0;
            void* varch = TMVA::gTools().GetChild(ch);
            while (varch) {
                TMVA::gTools().ReadAttr(varch, "VarIndex", varIdx);
                readVarInfo.ReadFromXML(varch);
                varch = TMVA::gTools().GetNextChild(varch);
                expressions.push_back(readVarInfo.GetExpression());
                if (verbose)
                    std::cout << "VERBOSE: " << readVarInfo.GetLabel() << " := " << readVarInfo.GetExpression() << std::endl;
            }
        }
        ch = TMVA::gTools().GetNextChild(ch);
    }
    TMVA::gTools().xmlengine().FreeDoc(doc);
    assert(readNVar == expressions.size());
    return expressions;
}

void ResetDeleteBranches(TTree* tree)
{
    TObjArray* branches = tree->GetListOfBranches();
    Int_t nb = branches->GetEntriesFast();
    for (Int_t i = 0; i < nb; ++i) {
        TBranch* br = (TBranch*) branches->UncheckedAt(i);
        if (br->InheritsFrom(TBranchElement::Class())) {
            ((TBranchElement*) br)->ResetDeleteObject();
        }
    }
}

bool isZH(const TString& process) {
    return (process.BeginsWith("ZnnH") || process=="ZH");
}

bool isWH(const TString& process) {
    return (process.BeginsWith("WlnH") || process=="WH");
}

bool isZJ(const TString& process) {
    return (process.BeginsWith("ZJetsHT") || process.BeginsWith("ZJetsPtZ") || process=="Zj");
}

bool isWJ(const TString& process) {
    return (process.BeginsWith("WJetsPtW") || process=="Wj");
}

bool isTT(const TString& process) {
    return (process=="TTFullLeptMG" || process=="TTSemiLeptMG" || process=="TTHadronicMG" || process=="TT");
}


//TH1F* h_delta_ew_hw8 = new TH1F("h_delta_ew_hw8", "HW differential NLO EWK correction", 95, 25., 500.);
//TH1F* h_delta_ew_hll8 = new TH1F("h_delta_ew_hll8", "Hll differential NLO EWK correction", 95, 25., 500.);
//TH1F* h_delta_ew_hnn8 = new TH1F("h_delta_ew_hnn8", "H#nu#nu differential NLO EWK correction", 95, 25., 500.);
Double_t hw8_contents[95] = {0.00200484, -9.13928e-05, -0.00219974, -0.00202884, -0.00219318, -0.00299564, -0.00336771, -0.00409143, -0.0059936, -0.00750034, -0.00881597, -0.0111589, -0.0135328, -0.0151878, -0.0178145, -0.0189559, -0.0215716, -0.0233429, -0.0250507, -0.0278028, -0.0291708, -0.0319617, -0.0341056, -0.0365647, -0.0381327, -0.0399857, -0.0405974, -0.0444237, -0.0454833, -0.0482596, -0.0486239, -0.0534461, -0.0539187, -0.0559509, -0.0596689, -0.0610085, -0.0633554, -0.0639406, -0.0672306, -0.0699078, -0.0719664, -0.0740249, -0.0760835, -0.078142, -0.0802005, -0.0822591, -0.0843176, -0.0863761, -0.0884347, -0.0904932, -0.0925517, -0.0946103, -0.0966688, -0.0987273, -0.100786, -0.102844, -0.104903, -0.106961, -0.10902, -0.111079, -0.113137, -0.115196, -0.117254, -0.119313, -0.121371, -0.12343, -0.125488, -0.127547, -0.129605, -0.131664, -0.133722, -0.135781, -0.13784, -0.139898, -0.141957, -0.144015, -0.146074, -0.148132, -0.150191, -0.152249, -0.154308, -0.156366, -0.158425, -0.160483, -0.162542, -0.1646, -0.166659, -0.168718, -0.170776, -0.172835, -0.174893, -0.176952, -0.17901, -0.181069, -0.183127};
Double_t hll8_contents[95] = {0.000664024, -0.00357095, -0.00767076, -0.00967366, -0.0134844, -0.0157148, -0.0181885, -0.0209647, -0.0232788, -0.0252373, -0.0265634, -0.0275069, -0.0285776, -0.0281683, -0.0294206, -0.0299975, -0.0308047, -0.0311716, -0.030913, -0.0324821, -0.0323192, -0.0324639, -0.0319356, -0.0322621, -0.0331146, -0.0338905, -0.0345189, -0.0358591, -0.0358407, -0.040018, -0.0396389, -0.0407177, -0.0445103, -0.0441406, -0.0471215, -0.0463301, -0.0513777, -0.0536773, -0.0546446, -0.0568508, -0.0590333, -0.0612157, -0.0633981, -0.0655805, -0.067763, -0.0699454, -0.0721278, -0.0743103, -0.0764927, -0.0786751, -0.0808575, -0.08304, -0.0852224, -0.0874048, -0.0895872, -0.0917697, -0.0939521, -0.0961345, -0.098317, -0.100499, -0.102682, -0.104864, -0.107047, -0.109229, -0.111412, -0.113594, -0.115776, -0.117959, -0.120141, -0.122324, -0.124506, -0.126689, -0.128871, -0.131053, -0.133236, -0.135418, -0.137601, -0.139783, -0.141965, -0.144148, -0.14633, -0.148513, -0.150695, -0.152878, -0.15506, -0.157242, -0.159425, -0.161607, -0.16379, -0.165972, -0.168155, -0.170337, -0.172519, -0.174702, -0.176884};
Double_t hnn8_contents[95] = {0.0146846, 0.0136521, 0.0125801, 0.0117771, 0.010976, 0.00989665, 0.00929942, 0.00836484, 0.00781992, 0.00733247, 0.00688885, 0.00666833, 0.0063354, 0.00637412, 0.00662595, 0.0069015, 0.00716689, 0.00760953, 0.00823267, 0.00914484, 0.00960494, 0.0110894, 0.0122241, 0.0127155, 0.0126892, 0.0125873, 0.01278, 0.0128243, 0.0118519, 0.0116125, 0.0102697, 0.00960959, 0.00929141, 0.00807739, 0.00588976, 0.00522135, 0.00365527, 0.00214147, 0.000569382, 0.000322672, -0.0015679, -0.00345846, -0.00534903, -0.0072396, -0.00913017, -0.0110207, -0.0129113, -0.0148019, -0.0166924, -0.018583, -0.0204736, -0.0223641, -0.0242547, -0.0261453, -0.0280358, -0.0299264, -0.031817, -0.0337076, -0.0355981, -0.0374887, -0.0393793, -0.0412698, -0.0431604, -0.045051, -0.0469415, -0.0488321, -0.0507227, -0.0526132, -0.0545038, -0.0563944, -0.0582849, -0.0601755, -0.0620661, -0.0639566, -0.0658472, -0.0677378, -0.0696283, -0.0715189, -0.0734095, -0.0753001, -0.0771906, -0.0790812, -0.0809718, -0.0828623, -0.0847529, -0.0866435, -0.088534, -0.0904246, -0.0923152, -0.0942057, -0.0960963, -0.0979869, -0.0998774, -0.101768, -0.103659};

double ewkWlnH(double vpt){
    if (vpt<=25.) return  0.00200484;
    if (vpt>500.) return -0.183127;
    int corrBin = int((vpt - 25.) / 5.);
    return (hw8_contents[corrBin]);
}

double ewkZllH(double vpt){
    if (vpt<=25.) return  0.000664024;
    if (vpt>500.) return -0.176884;
    int corrBin = int((vpt - 25.) / 5.);
    return (hll8_contents[corrBin]);
}

double ewkZnnH(double vpt){
    if (vpt<=25.) return  0.0146846;
    if (vpt>500.) return -0.103659;
    int corrBin = int((vpt - 25.) / 5.);
    return (hnn8_contents[corrBin]);
}

Double_t nptZj0_contents[35] = {1.02398, 0.991482, 0.994731, 0.997151, 1.00393, 1.01583, 1.02323, 1.03775, 1.05809, 1.08722, 1.11672, 1.16524, 1.20547, 1.2458, 1.26416, 1.31411, 1.32892, 1.34192, 1.3495, 1.30883, 1.2895, 1.34906, 1.28428, 1.33623, 1.24839, 1.37714, 1.35714, 1.23531, 1.31833, 1.30503, 1.3387, 1.31301, 1.34943, 1.25387, 1.35959};  // overflow = 1.19855

Double_t nptZj1_contents[35] = {0.994659, 0.97887, 0.980088, 0.982158, 0.986126, 0.991133, 0.997379, 1.00386, 1.01533, 1.03179, 1.04813, 1.08557, 1.12985, 1.1638, 1.17632, 1.2076, 1.16545, 1.16063, 1.16332, 1.14679, 1.13625, 1.15521, 1.15901, 1.13647, 1.1199, 1.11558, 1.07588, 1.08819, 1.10949, 1.0724, 1.0907, 1.10413, 1.04055, 1.07357, 1.03401};  // overflow = 1.03687

double qcdZnnH(double vpt, unsigned int njets) {
    int corrBin = int(vpt / 10.);
    if (njets == 0) {
        if (vpt<0   ) return 1.0;
        if (vpt>350.) return 1.19855;
        return (nptZj0_contents[corrBin]);
    } else {
        if (vpt<0   ) return 1.0;
        if (vpt>350.) return 1.03687;
        return (nptZj1_contents[corrBin]);
    }
}

////////////////////////////////////////////////////////////////////////////////
/// Main                                                                     ///
////////////////////////////////////////////////////////////////////////////////



void TrimTree(TString process="ZnnH125", TString mvaMethod="BDT", Long64_t beginEntry=0, Long64_t endEntry=-1)
{
    // If process starts with "Step4_", then it runs from ready-made Step 4's;
    // else, it runs from Step 3's (default)
    // If mvaMethod contains exactly two ":", then it runs from ready-made Step 4's
    // and does the MVA tuning.

    // FIXME: store the timestamp of xml (Creator, Date under GeneralInfo)
    // ZnnH125 signal 19723 even 9882 odd 9841
    // ZnnH125 ZjLF 8705

    gROOT->SetBatch(1);
    gROOT->LoadMacro("HelperFunctions.h");  //< make functions visible to TTreeFormula

    if (!TString(gROOT->GetVersion()).Contains("5.34")) {
        std::cout << "INCORRECT ROOT VERSION! Please use 5.34:" << std::endl;
        std::cout << "source /uscmst1/prod/sw/cms/slc5_amd64_gcc462/lcg/root/5.34.02-cms/bin/thisroot.csh" << std::endl;
        std::cout << "Return without doing anything." << std::endl;
        return;
    }

    const bool isData = (process.BeginsWith("Data"));
    const bool reloadStep4 = process.Contains("Step4_");
    bool tuneMVA = (mvaMethod.CountChar(':') == 2);

    TString weightdirTuneMVA = "";
    TString outdirTuneMVA = outdir;
    int ntunes = 1;
    if (tuneMVA) {
        /// parse mvaMethod
        Ssiz_t pos = mvaMethod.Index(':');
        TString substr1 = mvaMethod(0, pos);
        TString substr2 = mvaMethod(pos+1, mvaMethod.Sizeof());
        pos = substr2.Index(':');
        TString substr3 = substr2(0, pos);
        TString substr4 = substr2(pos+1, substr2.Sizeof());

        weightdir         = "./";
        mvaMethod         = substr1;
        weightdirTuneMVA += substr3;
        outdirTuneMVA    += "skim" + substr3.ReplaceAll("weights","") + "/";
        ntunes            = substr4.Atoi();
        assert(ntunes >= 1);
    }
    if (tuneMVA && !reloadStep4) {
        std::cout << "MVA tuning must load from Step 4's. MVA tuning is turned off." << std::endl;
        tuneMVA = false;
        ntunes = 1;
    }

    TFile *input(0);
    TTree *inTree(0);
    if (!reloadStep4) {
        input = TFile::Open(indir + prefix + process + suffix);
    } else {
        input = TFile::Open(indirStep4 + process + suffix);
    }
    if (!input) {
        std::cout << "ERROR: Could not open input file." << std::endl;
        exit(1);
    }
    /// Make output directory if it doesn't exist
    if (gSystem->AccessPathName(outdir))
        gSystem->mkdir(outdir);
    if (tuneMVA && gSystem->AccessPathName(outdirTuneMVA))
        gSystem->mkdir(outdirTuneMVA);

    std::cout << "--- TrimTree                 : Using input file: " << input->GetName() << std::endl;

    TString outputname = "";
    if (!reloadStep4) {
        inTree = (TTree *) input->Get("tree");
        if (beginEntry == 0 && endEntry == -1)
            outputname = outdir + "Step4_" + process + suffix;
        else
            outputname = outdir + "Step4_" + process + Form("_%Li_%Li", beginEntry, endEntry) + suffix;
    } else {
        // Not supported
        if (!(beginEntry == 0 && endEntry == -1)) {
            std::cout << "WARNING: Reload of Step 4 doesn't support the use of beginEntry and endEntry. They are reset." << std::endl;
            beginEntry = 0;
            endEntry = -1;
        }
        if (!tuneMVA) {
            outputname = outdir + process + suffix;
        } else {
            outputname = outdirTuneMVA + process + suffix;
        }
        //outputname.ReplaceAll("Step4_", "Step4_reload_");
    }
    TFile *output = TFile::Open(outputname, "RECREATE");


    TString stitchcut = "";
#ifdef STITCH
    if (process == "WJetsPtW100")
        stitchcut = " && lheV_pt<180";  // append to the selection
    if (stitchcut != "")
        std::cout << "--- TrimTree                 : Stitching cut will be applied: " << stitchcut << std::endl;
#endif


    ///-- Read from global variables -------------------------------------------
    /// Set the variables to keep
    if (!reloadStep4) {
        const std::vector<TString> keepvariables(g_variables, g_variables + sizeof(g_variables)/sizeof(g_variables[0]));
        inTree->SetBranchStatus("*", 0);
        for (UInt_t i = 0; i < keepvariables.size(); i++) {
            inTree->SetBranchStatus(keepvariables.at(i), 1);
        }
    }

    /// Set the masses
    const std::vector<int> masses(g_masses, g_masses + sizeof(g_masses)/sizeof(g_masses[0]));

    /// Set the channels
    const std::vector<TString> channels(g_channels, g_channels + sizeof(g_channels)/sizeof(g_channels[0]));

    /// Set the weight files
    std::vector<TString> weights(g_weights, g_weights + sizeof(g_weights)/sizeof(g_weights[0]));
    if (tuneMVA) {
        weights.clear();
        weights.push_back(weightdirTuneMVA);
    }

    /// For TMVA
    TMVA::Tools::Instance(); //< This loads the library
    Float_t readerVars[MAXVAR];

    /// Get the systematic variations
    const std::vector<std::vector<std::pair<std::string, std::string> > >& systExpressions = GetSystExpressions();
    std::cout << "--- TrimTree                 : There are " << systExpressions.size()-1 << " systematic sources." << std::endl;


    /// Loop over channels
    for (UInt_t ichan = 0; ichan < channels.size(); ichan++) {
        const TString& channel = channels.at(ichan);
        std::cout << "--- TrimTree                 : Entering channel " << ichan+1 <<  "/" << channels.size() << ": " << channel << std::endl;

        /// Get the selections
        const std::vector<std::string>& selExpressions = GetSelExpressions(channel.Data());
        std::cout << "--- TrimTree                 : There are " << selExpressions.size() << " selections." << std::endl;

        const std::vector<std::string>& selMjjExpressions = GetSelMjjExpressions(channel.Data());
        assert(selMjjExpressions.size() == 1);

        /// Get the MC event weight
        const std::string weigExpression = GetWeightExpressions()[channel.Data()];

        /// Get the trigger
        const std::string trigExpression = GetTriggerExpressions()[channel.Data()];

        if (reloadStep4) {
            inTree = (TTree *) input->Get(Form("tree_%s_test", channel.Data()));  // only the test tree
            inTree->SetBranchStatus("*", 1);
            inTree->SetBranchStatus("BDT*", 0);  // remove all BDT from the previous Step 4
        }
        assert(inTree != 0);

        /// For input variables;
        std::vector<TString> inputExpressions;

        ///-- Setup TMVA Reader ----------------------------------------------------
        std::vector<TMVA::Reader *> readers;
        std::vector<std::pair<TString,TString> > readerNames;
        std::vector<std::vector<UInt_t> > readerVarIndices;
        for (UInt_t imass = 0; imass < masses.size(); imass++) {
            const int mass = masses.at(imass);
            for (UInt_t iread = 0; iread < weights.size(); iread++) {

                /// Loop over BDT[0-99] in the directory
                for (Int_t jread = 0; jread < ntunes; jread++) {
                    TMVA::Reader * reader = new TMVA::Reader("!Color:!Silent");
                    const TString weightdir1 = weightdir + weights.at(iread);
                    TString mvaMethod1 = (ntunes==1) ? mvaMethod : mvaMethod+Form("%i",jread);
                    const TString weightfile = Form("%s/%i/TMVAClassification_%s.weights_%s.xml", weightdir1.Data(), mass, mvaMethod1.Data(), channel.Data());
                    std::cout << "--- TrimTree                 : Reading " << weightfile << std::endl;

                    /// Get input variables from XML
                    const std::vector<TString>& variables = GetVariablesFromXML(weightfile);
                    assert(variables.size() < MAXVAR);
                    std::vector<UInt_t> varIndices;
                    for (UInt_t iexpr = 0; iexpr < variables.size(); iexpr++) {
                        /// Add this variable to reader
                        const TString expr = variables.at(iexpr);
                        reader->AddVariable(expr, &readerVars[iexpr]);

                        /// Check if this variable has been used previously
                        bool found = false;
                        for (UInt_t jexpr = 0; jexpr < inputExpressions.size(); jexpr++) {
                            if (expr == inputExpressions.at(jexpr)) {
                                varIndices.push_back(jexpr);
                                found = true;
                                break;
                            }
                        }
                        if (!found) {
                            inputExpressions.push_back(expr);
                            varIndices.push_back(inputExpressions.size()-1);
                        }
                    }
                    assert(varIndices.size() == variables.size());

                    /// Load TMVA weights
                    reader->BookMVA(mvaMethod1 + " method", weightfile);
                    readers.push_back(reader);
                    TString readerName = Form("%s_%i", weights.at(iread).Data(), mass);
                    // FIXME: only remove one layer
                    Ssiz_t pos = readerName.Index('/');
                    readerName = readerName(pos+1, readerName.Sizeof());
                    readerName.ReplaceAll("weights","").ToLower();
                    readerNames.push_back(std::pair<TString,TString>(mvaMethod1,readerName));
                    readerVarIndices.push_back(varIndices);
                }
            }
        }
        assert(readers.size() == masses.size() * weights.size() * ntunes);
        assert(readers.size() < MAXWEIGHT);


        ///-- Setup outTree branches -------------------------------------------
        TTree *outTree = inTree->CloneTree(0); // Do no copy the data yet
        outTree->SetAutoFlush(-120000000);  // default was 30 MB
        if (!reloadStep4)
            outTree->SetName(Form("%s_%s", inTree->GetName(), channel.Data()));
        /// The clone should not delete any shared i/o buffers.
        ResetDeleteBranches(outTree);

        char* processname = const_cast<char *>((process+"\0").Data());
        int nselect = selExpressions.size();
        int nsyst = systExpressions.size();
        bool isSignal=false, isControl=false;  // an event can be both signal and control due to systematics
        float scalefactor=1., slopefactor=0.;
        float weightSignalEWKNew=1., weightSignalEWKNewP=1., weightSignalEWKNewM=1.;  // FIXME: move to Step 3
        float weightSignalQCDNew=1., weightSignalQCDNewP=1., weightSignalQCDNewM=1.;  // FIXME: move to Step 3
        bool selectFlags[nselect][nsyst];
        bool selectMjj[nsyst];
        float weightsMC[MAXSYST];
        float discriminants[MAXWEIGHT][MAXSYST];
        assert(nselect < MAXSELECT);
        assert(nsyst < MAXSYST);

        if (!reloadStep4) {
            outTree->Branch("processname", processname, "processname/C");
            outTree->Branch("nselect", &nselect, "nselect/I");
            outTree->Branch("nsyst", &nsyst, "nsyst/I");
            outTree->Branch("isSignal", &isSignal, "isSignal/b");
            outTree->Branch("isControl", &isControl, "isControl/b");
            outTree->Branch("scalefactor", &scalefactor, "scalefactor/F");  // FIXME: not implemented
            outTree->Branch("slopefactor", &slopefactor, "slopefactor/F");
            outTree->Branch("weightSignalEWKNew", &weightSignalEWKNew, "weightSignalEWKNew/F");
            outTree->Branch("weightSignalEWKNewP", &weightSignalEWKNewP, "weightSignalEWKNewP/F");
            outTree->Branch("weightSignalEWKNewM", &weightSignalEWKNewM, "weightSignalEWKNewM/F");
            outTree->Branch("weightSignalQCDNew", &weightSignalQCDNew, "weightSignalQCDNew/F");
            outTree->Branch("weightSignalQCDNewP", &weightSignalQCDNewP, "weightSignalQCDNewP/F");
            outTree->Branch("weightSignalQCDNewM", &weightSignalQCDNewM, "weightSignalQCDNewM/F");
            outTree->Branch("selectFlags", selectFlags, Form("selectFlags[%i][%i]/b", nselect, nsyst));  // 2-d arrays must be of fixed size.
            outTree->Branch("selectMjj", selectMjj, "selectMjj[nsyst]/b");
            outTree->Branch("weightsMC", weightsMC, "weightsMC[nsyst]/F");
        }

        for (UInt_t iread = 0; iread < readers.size(); iread++) {
            TString name = readerNames.at(iread).first + readerNames.at(iread).second;
            outTree->Branch(name, &discriminants[iread][0], name+"[nsyst]/F");
        }

        if (beginEntry == 0 && endEntry == -1)
            std::cout << "--- Processing: " << inTree->GetEntriesFast() << " events" << std::endl;
        else
            std::cout << "--- Processing: " << beginEntry << " - " << endEntry << " from " << inTree->GetEntriesFast() << " events" << std::endl;
        TStopwatch sw;
        sw.Start();


        /// Create TTreeFormulasWithSyst
        TTreeFormulaWithSyst * ttfs = 0;
        std::vector<TTreeFormulaWithSyst *>::const_iterator formIt, formItEnd;
        std::vector<TTreeFormulaWithSyst *> inputFormulas;
        for (UInt_t iexpr = 0; iexpr < inputExpressions.size(); iexpr++) {
            TString expr = inputExpressions.at(iexpr);
            ttfs = new TTreeFormulaWithSyst("ttfin", expr, inTree, systExpressions, isData);
            inputFormulas.push_back(ttfs);
        }
        std::vector<TTreeFormulaWithSyst *> selFormulas;
        for (UInt_t iexpr = 0; !reloadStep4 && (iexpr < selExpressions.size()); iexpr++) {  // not needed when remaking Step 4
            TString expr = selExpressions.at(iexpr);
            expr += stitchcut;
#ifdef SPECIALQCD
            if (process.BeginsWith("QCD")) {
                expr.ReplaceAll("mindPhiMETJet_dPhi>0.5 && ", "");
                expr.ReplaceAll("mindPhiMETJet_dPhi>0.7 && ", "");
                expr.ReplaceAll("abs(deltaPhi(METtype1corr.phi,METnoPUCh.phi))<0.5 && ", "");  // FIXME: don't remove this?
                expr.ReplaceAll("METtype1corr.et/sqrt(METtype1corr.sumet)>3 && ", "");
            }
#endif
            ttfs = new TTreeFormulaWithSyst("ttfsel", expr, inTree, systExpressions, isData);
            selFormulas.push_back(ttfs);
        }
        TTreeFormulaWithSyst * weigFormula = new TTreeFormulaWithSyst("ttfweig", weigExpression.c_str(), inTree, systExpressions, isData);
        TTreeFormulaWithSyst * trigFormula = new TTreeFormulaWithSyst("ttftrig", trigExpression.c_str(), inTree, systExpressions, isData);
        TTreeFormulaWithSyst * mjjFormula = new TTreeFormulaWithSyst("ttfmjj", selMjjExpressions.front().c_str(), inTree, systExpressions, isData);

        /// Create TTreeFormulas
        double WJSlopeErr = 0.0020;
        double ZJSlopeErr = 0.0025;
        double WJSlope = -0.50;
        double ZJSlope = -1.00;
        TTreeFormula * genwptFormula = new TTreeFormula("ttfgenwpt", "genW.pt", inTree);
        TTreeFormula * genzptFormula = new TTreeFormula("ttfgenzpt", "genZ.pt", inTree);
        TTreeFormula * genhptFormula = new TTreeFormula("ttfgenhpt", "genH.pt", inTree);
        TTreeFormula * scaleFormula = new TTreeFormula("ttfscale", "(Sum$(abs(hJet_flavour)==5)==1) ? 2.0 : 1.0", inTree);
        TTreeFormula * naJetsFormula = new TTreeFormula("ttfnaJets", "naJets_Znn", inTree);
        genwptFormula->SetQuickLoad(1);
        genzptFormula->SetQuickLoad(1);
        genhptFormula->SetQuickLoad(1);
        scaleFormula->SetQuickLoad(1);
        naJetsFormula->SetQuickLoad(1);

        /// Speed up if a systematic does nothing to input variables
        std::vector<Bool_t> changed;
        changed.resize(nsyst, kFALSE);
        for (UInt_t i = 0; i < inputFormulas.size(); i++)
            inputFormulas.at(i)->Is_Changed(changed);
        //for (UInt_t isyst = 0; isyst < changed.size(); isyst++)
        //    if (!changed.at(isyst))
        //        std::cout << "VERBOSE: systematic " << isyst << " doesn't affect input variables." << std::endl;


        ///-- Loop over events -------------------------------------------------
        /// If isData and not passTrigger, the event is skipped
        /// If not selected as signal or control, the event is skipped
        Int_t curTree = inTree->GetTreeNumber();
        const Long64_t nentries = inTree->GetEntries();
        if (endEntry < 0)  endEntry = nentries;

        Long64_t ievt = 0;
        for (ievt=TMath::Max(ievt, beginEntry); ievt<TMath::Min(nentries, endEntry); ievt++) {
            if (ievt % 2000 == 0)
                std::cout << "--- ... Processing event: " << ievt << std::endl;

            const Long64_t local_entry = inTree->LoadTree(ievt);  // faster, but only for TTreeFormula
            if (local_entry < 0)  break;
            //inTree->GetEntry(ievt);  // same event as received by LoadTree()

            if (inTree->GetTreeNumber() != curTree) {
                curTree = inTree->GetTreeNumber();
                for (formIt=inputFormulas.begin(), formItEnd=inputFormulas.end(); formIt!=formItEnd; formIt++)
                    (*formIt)->UpdateFormulaLeaves();  // if using TChain
                for (formIt=selFormulas.begin(), formItEnd=selFormulas.end(); formIt!=formItEnd; formIt++)
                    (*formIt)->UpdateFormulaLeaves();  // if using TChain
                weigFormula->UpdateFormulaLeaves();    // if using TChain
                trigFormula->UpdateFormulaLeaves();    // if using TChain
                mjjFormula->UpdateFormulaLeaves();     // if using TChain
                genwptFormula->UpdateFormulaLeaves();  // if using TChain
                genzptFormula->UpdateFormulaLeaves();  // if using TChain
                genhptFormula->UpdateFormulaLeaves();  // if using TChain
                scaleFormula->UpdateFormulaLeaves();   // if using TChain
                naJetsFormula->UpdateFormulaLeaves();  // if using TChain
            }

            /// These need to be called when arrays of variable size are used in TTree.
            for (formIt=inputFormulas.begin(), formItEnd=inputFormulas.end(); formIt!=formItEnd; formIt++)
                (*formIt)->GetNdata();
            for (formIt=selFormulas.begin(), formItEnd=selFormulas.end(); formIt!=formItEnd; formIt++)
                (*formIt)->GetNdata();
            weigFormula->GetNdata();
            trigFormula->GetNdata();
            mjjFormula->GetNdata();
            genwptFormula->GetNdata();
            genzptFormula->GetNdata();
            genhptFormula->GetNdata();
            scaleFormula->GetNdata();
            naJetsFormula->GetNdata();

            if (isData) {
                bool passTrigger = trigFormula->EvalBoolInstance(0);
                if (!passTrigger)  continue;
            }

            /// Loop over systematics
            isSignal = false;
            isControl = false;
            for (UInt_t isyst = 0; isyst < systExpressions.size(); isyst++) {
                /// Loop over selections
                for (UInt_t isel = 0; isel < selFormulas.size(); isel++) {
                    selectFlags[isel][isyst] = selFormulas.at(isel)->EvalBoolInstance(isyst);
                    if(selectFlags[isel][isyst]) {
                        if (isel == 0)  // signal selection
                            isSignal = true;
                        else  // control selection
                            isControl = true;
                    }
                }
            }  // end loop over systematics

            if (!reloadStep4 && !isSignal && !isControl)  continue;  // not needed when remaking Step 4

            /// Loop over systematics again, evaluate BDT responses
            for (UInt_t isyst = 0; isyst < systExpressions.size(); isyst++) {
                if (isyst == 0 || changed.at(isyst)) {
                    /// Loop over readers
                    for (UInt_t iread = 0; iread < readers.size(); iread++) {
                        const std::vector<UInt_t>& varIndices = readerVarIndices.at(iread);
                        for (UInt_t iexpr = 0; iexpr < varIndices.size(); iexpr++) {
                            ttfs = inputFormulas.at(varIndices.at(iexpr));
                            readerVars[iexpr] = ttfs->EvalInstance(isyst);
                            //std::cout << readerVars[iexpr] << std::endl;
                        }
                        discriminants[iread][isyst] = readers.at(iread)->EvaluateMVA(readerNames.at(iread).first + " method");
                        //std::cout << "-->" << discriminants[iread][isyst] << std::endl;
                    }
                /// Speed up if a systematic does nothing to input variables
                } else {
                    /// Loop over readers
                    for (UInt_t iread = 0; iread < readers.size(); iread++) {
                        discriminants[iread][isyst] = discriminants[iread][0];
                    }
                }

                /// Evaluate event weight
                if (isData)
                    weightsMC[isyst] = 1.0;
                else
                    weightsMC[isyst] = weigFormula->EvalInstance(isyst);

                /// Evaluate Mjj selection
                selectMjj[isyst] = mjjFormula->EvalInstance(isyst);

                /// Evaluate scale factors and slopes
                if (isZJ(process)) {
                    double genzpt = genzptFormula->EvalInstance();
                    slopefactor = apply_pt_slope(genzpt, ZJSlope * ZJSlopeErr, 130);
                    scalefactor = scaleFormula->EvalInstance();
                } else if (isWJ(process)) {
                    double genwpt = genwptFormula->EvalInstance();
                    slopefactor = apply_pt_slope(genwpt, WJSlope * WJSlopeErr, 150);
                    scalefactor = scaleFormula->EvalInstance();
                } else if (isZH(process)) {
                    double genhpt = genhptFormula->EvalInstance();
                    double ewk = ewkZnnH(genhpt);
                    weightSignalEWKNew  = 1.0 + ewk;
                    weightSignalEWKNewP = 1.0 + ewk * 0.0;
                    weightSignalEWKNewM = 1.0 + ewk * 2.0;
                    double genzpt = genzptFormula->EvalInstance();
                    unsigned int naJets = naJetsFormula->EvalInstance();
                    double qcd = qcdZnnH(genzpt, naJets);
                    weightSignalQCDNew  = qcd;
                    weightSignalQCDNewP = qcd + qcd - 1.0;
                    weightSignalQCDNewP = qcd - qcd + 1.0;
                } else if (isWH(process)) {
                    double genhpt = genhptFormula->EvalInstance();
                    double ewk = ewkWlnH(genhpt);
                    weightSignalEWKNew  = 1.0 / (1.0 + ewk);
                    weightSignalEWKNewP = 1.0 / (1.0 + ewk * 0.0);
                    weightSignalEWKNewM = 1.0 / (1.0 + ewk * 2.0);
                } else if (process == "ZbbHinv") {
                    double genhpt = genhptFormula->EvalInstance();
                    double ewk = ewkZnnH(genhpt);
                    weightSignalEWKNew  = 1.0 + ewk;
                    weightSignalEWKNewP = 1.0 + ewk * 0.0;
                    weightSignalEWKNewM = 1.0 + ewk * 2.0;
                }

            }  // end loop over systematics again

            inTree->GetEntry(ievt);  // same event as received by LoadTree()
            outTree->Fill();  // fill it!
        }  // end loop over events
        if (endEntry == nentries)  endEntry = -1;

        /// Get elapsed time
        sw.Stop();
        std::cout << "--- End of event loop: ";
        sw.Print();

        std::cout << "--- TrimTree                 : Keep " << outTree->GetEntriesFast() << " out of " << nentries << " events." << std::endl;

        output->cd();
        outTree->Write();

        delete outTree;
        for (UInt_t iread = 0; iread < readers.size(); iread++)
            delete readers.at(iread);
        for (formIt=inputFormulas.begin(), formItEnd=inputFormulas.end(); formIt!=formItEnd; formIt++)
            delete *formIt;
        for (formIt=selFormulas.begin(), formItEnd=selFormulas.end(); formIt!=formItEnd; formIt++)
            delete *formIt;
        delete weigFormula;
        delete trigFormula;
        delete mjjFormula;
        delete genwptFormula;
        delete genzptFormula;
        delete genhptFormula;
        delete scaleFormula;
        delete naJetsFormula;
    }  // end loop over channels

    /// Split
    if (!reloadStep4) {
        TCut trainCut = "isSignal && (EVENT.event %2 == 0)";
        TCut testCut  = "isSignal && (EVENT.event %2 == 1)";
        TCut ctrlCut  = "isControl";
        TString outputname_tmp = outputname;
        outputname_tmp.ReplaceAll(".root","_tmp.root");
        TFile *output_tmp = TFile::Open(outputname_tmp, "RECREATE");
        systFlags_id.Write("systFlags_id");
        selectFlags_id.Write("selectFlags_id");

        for (UInt_t ichan = 0; ichan < channels.size(); ichan++) {
            const TString& channel = channels.at(ichan);
            TTree *outTree = (TTree *) output->Get(Form("%s_%s", inTree->GetName(), channel.Data()));

            TTree *testTree(0), *trainTree(0), *ctrlTree(0);
            output_tmp->cd();
            if (isData) {
                testTree = (TTree*) outTree->CopyTree("isSignal");
                testTree->SetName(Form("%s_%s_test", inTree->GetName(), channel.Data()));
                testTree->Write();
            } else {
                trainTree = (TTree*) outTree->CopyTree(trainCut);
                trainTree->SetName(Form("%s_%s_train", inTree->GetName(), channel.Data()));
                trainTree->Write();
                testTree = (TTree*) outTree->CopyTree(testCut);
                testTree->SetName(Form("%s_%s_test", inTree->GetName(), channel.Data()));
                testTree->Write();
            }
            ctrlTree = (TTree*) outTree->CopyTree(ctrlCut);
            ctrlTree->SetName(Form("%s_%s_ctrl", inTree->GetName(), channel.Data()));
            ctrlTree->Write();

            if (isData)
                std::cout << "--- TrimTree                 : Split " << channel << " into " << testTree->GetEntriesFast() << "(test), " << ctrlTree->GetEntriesFast() << "(ctrl)." << std::endl;
            else
                std::cout << "--- TrimTree                 : Split " << channel << " into " << trainTree->GetEntriesFast() << "(train), " << testTree->GetEntriesFast() << "(test), " << ctrlTree->GetEntriesFast() << "(ctrl)." << std::endl;

            delete trainTree;
            delete testTree;
            delete ctrlTree;
        }  // end loop over channels

        output_tmp->Close();
        output->Close();
        input->Close();
        delete input;
        delete output;
        delete output_tmp;

        /// Rename the file
        gSystem->Exec("mv " + outputname_tmp + " " + outputname);

    } else {
        output->Close();
        input->Close();
        delete input;
        delete output;
    }

    std::cout << "==> TrimTree is done!" << std::endl << std::endl;
    return;
}
