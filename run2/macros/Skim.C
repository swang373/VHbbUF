#include <iostream>

#include "TFile.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TString.h"
#include "TCut.h"
#include "TH1.h"
#include "TSystem.h"
#include "TROOT.h"

#include "XSec_8TeV19invfb.h"

const std::string tagMC        = "Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC";
const std::string tagData      = "Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_DATA";

void Skim(TString process="ZnnH125")
{
    gROOT->SetBatch(1);

    TString dirMC   = "dcache:/pnfs/cms/WAX/11/store/user/lpchbb/apana/" + tagMC + "/";
    TString dirData = "dcache:/pnfs/cms/WAX/11/store/user/lpchbb/apana/" + tagData + "/";

    TString fname   = "";
    TString outdir  = "skim/";
    TString prefix  = "Step2_";
    TString suffix  = ".root";
    TString dijet   = "DiJetPt_";

    TChain * chain  = new TChain("tree");

    // ZnnH
    if (process == "ZnnH110") {
        fname = dirMC + dijet + "ZH_ZToNuNu_HToBB_M-110_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);
    } else if (process == "ZnnH115") {
        fname = dirMC + dijet + "ZH_ZToNuNu_HToBB_M-115_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);
    } else if (process == "ZnnH120") {
        fname = dirMC + dijet + "ZH_ZToNuNu_HToBB_M-120_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);
    } else if (process == "ZnnH125") {
        fname = dirMC + dijet + "ZH_ZToNuNu_HToBB_M-125_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);
    } else if (process == "ZnnH130") {
        fname = dirMC + dijet + "ZH_ZToNuNu_HToBB_M-130_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);
    } else if (process == "ZnnH135") {
        fname = dirMC + dijet + "ZH_ZToNuNu_HToBB_M-135_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);
    } else if (process == "ZnnH140") {
        fname = dirMC + dijet + "ZH_ZToNuNu_HToBB_M-140_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);
    } else if (process == "ZnnH145") {
        fname = dirMC + dijet + "ZH_ZToNuNu_HToBB_M-145_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);
    } else if (process == "ZnnH150") {
        fname = dirMC + dijet + "ZH_ZToNuNu_HToBB_M-150_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);

    // WlnH
    } else if (process == "WlnH110") {
        fname = dirMC + dijet + "WH_WToLNu_HToBB_M-110_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);
    } else if (process == "WlnH115") {
        fname = dirMC + dijet + "WH_WToLNu_HToBB_M-115_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);
    } else if (process == "WlnH120") {
        fname = dirMC + dijet + "WH_WToLNu_HToBB_M-120_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);
    } else if (process == "WlnH125") {
        fname = dirMC + dijet + "WH_WToLNu_HToBB_M-125_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);
    } else if (process == "WlnH130") {
        fname = dirMC + dijet + "WH_WToLNu_HToBB_M-130_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);
    } else if (process == "WlnH135") {
        fname = dirMC + dijet + "WH_WToLNu_HToBB_M-135_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);
    } else if (process == "WlnH140") {
        fname = dirMC + dijet + "WH_WToLNu_HToBB_M-140_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);
    } else if (process == "WlnH145") {
        fname = dirMC + dijet + "WH_WToLNu_HToBB_M-145_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);
    } else if (process == "WlnH150") {
        fname = dirMC + dijet + "WH_WToLNu_HToBB_M-150_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);

    // ZllH
    } else if (process == "ZllH110") {
        fname = dirMC + dijet + "ZH_ZToLL_HToBB_M-110_8TeV-powheg-herwigpp3" + suffix;  // NOTE: '3' in name
        chain->Add(fname);
    } else if (process == "ZllH115") {
        fname = dirMC + dijet + "ZH_ZToLL_HToBB_M-115_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);
    } else if (process == "ZllH120") {
        fname = dirMC + dijet + "ZH_ZToLL_HToBB_M-120_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);
    } else if (process == "ZllH125") {
        fname = dirMC + dijet + "ZH_ZToLL_HToBB_M-125_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);
    } else if (process == "ZllH130") {
        fname = dirMC + dijet + "ZH_ZToLL_HToBB_M-130_8TeV-powheg-herwigpp3" + suffix;  // NOTE: '3' in name
        chain->Add(fname);
    } else if (process == "ZllH135") {
        fname = dirMC + dijet + "ZH_ZToLL_HToBB_M-135_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);
    } else if (process == "ZllH140") {
        fname = dirMC + dijet + "ZH_ZToLL_HToBB_M-140_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);
    } else if (process == "ZllH145") {
        fname = dirMC + dijet + "ZH_ZToLL_HToBB_M-145_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);
    } else if (process == "ZllH150") {
        fname = dirMC + dijet + "ZH_ZToLL_HToBB_M-150_8TeV-powheg-herwigpp" + suffix;
        chain->Add(fname);

    // WJets
    } else if (process == "WJetsPtW100") {
        fname = dirMC + dijet + "WJetsToLNu_PtW-100_TuneZ2star_8TeV-madgraph" + suffix;
        chain->Add(fname);

    // ZJets
    } else if (process == "ZJetsPtZ100") {
        fname = dirMC + dijet + "ZJetsToNuNu_PtZ-100_8TeV-madgraph" + suffix;
        chain->Add(fname);

    // DYJets
    } else if (process == "DYJetsM50"   ) {
        fname = dirMC + dijet + "DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball" + suffix;
        chain->Add(fname);
    } else if (process == "DYJetsPtZ70" ) {
        fname = dirMC + dijet + "DYJetsToLL_PtZ-70To100_TuneZ2star_8TeV-madgraph-tarball" + suffix;
        chain->Add(fname);
    } else if (process == "DYJetsPtZ100") {
        fname = dirMC + dijet + "DYJetsToLL_PtZ-100_TuneZ2star_8TeV-madgraph" + suffix;
        chain->Add(fname);

    // VV
    } else if (process == "WW") {
        fname = dirMC + dijet + "WW_TuneZ2star_8TeV_pythia6_tauola" + suffix;
        chain->Add(fname);
    } else if (process == "WZ") {
        fname = dirMC + dijet + "WZ_TuneZ2star_8TeV_pythia6_tauola" + suffix;
        chain->Add(fname);
    } else if (process == "ZZ") {
        fname = dirMC + dijet + "ZZ_TuneZ2star_8TeV_pythia6_tauola" + suffix;
        chain->Add(fname);

    // TT
    } else if (process == "TTPowheg") {
        fname = dirMC + dijet + "TT_CT10_TuneZ2star_8TeV-powheg-tauola" + suffix;
        chain->Add(fname);
        fname = dirMC + dijet + "TT_CT10_TuneZ2star_8TeV-powheg-tauola-7M" + suffix;
        chain->Add(fname);

    } else if (process == "TTFullLeptMG") {
        fname = dirMC + dijet + "TTJets_FullLeptMGDecays_8TeV-madgraph" + suffix;
        chain->Add(fname);

    } else if (process == "TTSemiLeptMG") {
        fname = dirMC + dijet + "TTJets_SemiLeptMGDecays_8TeV-madgraph" + suffix;
        chain->Add(fname);

    } else if (process == "TTHadronicMG") {
        fname = dirMC + dijet + "TTJets_HadronicMGDecays_8TeV-madgraph" + suffix;
        chain->Add(fname);

    // T/Tbar
    } else if (process == "T_tW"   ) {
        fname = dirMC + dijet + "T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola" + suffix;
        chain->Add(fname);
    } else if (process == "Tbar_tW") {
        fname = dirMC + dijet + "Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola" + suffix;
        chain->Add(fname);
    } else if (process == "T_s"    ) {
        fname = dirMC + dijet + "T_s-channel_TuneZ2star_8TeV-powheg-tauola" + suffix;
        chain->Add(fname);
    } else if (process == "Tbar_s" ) {
        fname = dirMC + dijet + "Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola" + suffix;
        chain->Add(fname);
    } else if (process == "T_t"    ) {
        fname = dirMC + dijet + "T_t-channel_TuneZ2star_8TeV-powheg-tauola" + suffix;
        chain->Add(fname);
    } else if (process == "Tbar_t" ) {
        fname = dirMC + dijet + "Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola" + suffix;
        chain->Add(fname);

    // Data
    } else if (process == "MET_ReReco") {
        fname = dirData + dijet + "METRun2012A-13Jul2012-v1_V42_v4" + suffix;
        chain->Add(fname);
        fname = dirData + dijet + "METRun2012A-recover-06Aug2012-v1_V42_v4" + suffix;
        chain->Add(fname);
        fname = dirData + dijet + "METRun2012B-13Jul2012-v1_v1_V42_v4" + suffix;
        chain->Add(fname);
        fname = dirData + dijet + "METRun2012C-24Aug2012-v1_V42_v4" + suffix;
        chain->Add(fname);
        fname = dirData + dijet + "METRun2012C-EcalRecover_11Dec2012-v1_v2" + suffix;
        chain->Add(fname);

    } else if (process == "MET_Prompt") {
        fname = dirData + dijet + "METRun2012Cv2-PromptReco-v2EdmV42_v3" + suffix;  // Run 201191 has to be excluded
        chain->Add(fname);
        fname = dirData + dijet + "METRun2012D" + suffix;
        chain->Add(fname);

    } else if (process == "SingleMu_ReReco") {
        fname = dirData + dijet + "SingleMuRun2012AJul13EdmV42" + suffix;
        chain->Add(fname);
        fname = dirData + dijet + "SingleMuRun2012AAug06EdmV42" + suffix;
        chain->Add(fname);
        fname = dirData + dijet + "SingleMuRun2012BJul13EdmV42" + suffix;
        chain->Add(fname);
        fname = dirData + dijet + "SingleMuRun2012CAug24RerecoEdmV42" + suffix;
        chain->Add(fname);
        fname = dirData + dijet + "SingleMuRun2012C-EcalRecover_11Dec2012-v1_v2" + suffix;
        chain->Add(fname);

    } else if (process == "SingleMu_Prompt") {
        fname = dirData + dijet + "SingleMuRun2012CPromptv2EdmV42" + suffix;  // Run 201191 has to be excluded
        chain->Add(fname);
        fname = dirData + dijet + "SingleMuRun2012CPromptV2TopUpEdmV42" + suffix;
        chain->Add(fname);
        fname = dirData + dijet + "SingleMuRun2012D-PromptReco-v1" + suffix;
        chain->Add(fname);

    } else if (process == "SingleEl_ReReco") {
        fname = dirData + dijet + "SingleElectronRun2012AJul13EdmV42b" + suffix;
        chain->Add(fname);
        fname = dirData + dijet + "SingleElectronRun2012AAug06EdmV42" + suffix;
        chain->Add(fname);
        fname = dirData + dijet + "SingleElectronRun2012BJul13EdmV42" + suffix;
        chain->Add(fname);
        fname = dirData + dijet + "SingleElectronRun2012CAug24RerecoEdmV42" + suffix;
        chain->Add(fname);
        fname = dirData + dijet + "SingleElectronRun2012C-EcalRecover_11Dec2012-v1_v2" + suffix;
        chain->Add(fname);

    } else if (process == "SingleEl_Prompt") {
        fname = dirData + dijet + "SingleElectronRun2012CPromptv2EdmV42" + suffix;  // Run 201191 has to be excluded
        chain->Add(fname);
        fname = dirData + dijet + "SingleElectronRun2012CPromptV2TopUpEdmV42" + suffix;
        chain->Add(fname);
        fname = dirData + dijet + "SingleElectronRun2012D-PromptReco-v1_v3" + suffix;
        chain->Add(fname);

    } else if (process == "DoubleEl_ReReco") {
        fname = dirData + dijet + "DoubleElectron_Run2012A-13Jul2012-v1_ProcFIXED" + suffix;
        chain->Add(fname);
        fname = dirData + dijet + "DoubleElectron_Run2012A-recover-06Aug2012-v1_ProcV2" + suffix;
        chain->Add(fname);
        fname = dirData + dijet + "DoubleElectron_Run2012B-13Jul2012-v1_ProcFIXED" + suffix;
        chain->Add(fname);
        fname = dirData + dijet + "DoubleElectronRun2012CAug24RerecoEdmV42" + suffix;
        chain->Add(fname);
        fname = dirData + dijet + "DoubleElectronRun2012C-EcalRecover_11Dec2012-v1_v2" + suffix;
        chain->Add(fname);

    } else if (process == "DoubleEl_Prompt") {
        fname = dirData + dijet + "DoubleElectron_Run2012C-PromptReco-v2_HBB_EDMNtupleV42_ProcV1_noRun201191" + suffix;  // Run 201191 has to be excluded
        chain->Add(fname);
        fname = dirData + dijet + "DoubleElectron_Run2012C-PromptReco-v2_HBB_EDMNtupleV42_ProcV2" + suffix;
        chain->Add(fname);
        fname = dirData + dijet + "DoubleElectronRun2012D" + suffix;
        chain->Add(fname);

    } else {
        std::cout << "Error: Unrecognized process. I quit!" << std::endl;
        return;
    }

    // Make selection
    TCut selection = (skimZll || skimWln || skimZnn);
    selection += skimHbb;

    if (process.Contains("ReReco") || process.Contains("Prompt")) {
        if (process.Contains("MET")) {
            selection += (trigZnn);
        } else if (process.Contains("SingleMu")) {
            selection += (trigZmm || trigWmn);
        } else if (process.Contains("SingleEl")) {
            selection += (trigWen);
        } else if (process.Contains("DoubleEl")) {
            selection += (trigZee);
        } else {
            std::cout << "Error: Unrecognized dataset. I quit!" << std::endl;
            return;
        }
    }

    if (process.Contains("Prompt")) {
        selection += exclRuns;
    }

    selection.Print();


    // Sum Count, CountWithPU, CountWithPUP, CountWithPUM, ...
    TObjArray * files = chain->GetListOfFiles();
    TIter next(files);
    TChainElement * chainElem = 0;
    TFile * f2 = 0;

    TH1F * count = new TH1F("Count", "Count", 1, 0, 2);
    TH1F * countWithPU = new TH1F("CountWithPU", "CountWithPU", 1, 0, 2);
    TH1F * countWithPUP = new TH1F("CountWithPUP", "CountWithPU plus one sigma", 1, 0, 2);
    TH1F * countWithPUM = new TH1F("CountWithPUM", "CountWithPU minus one sigma", 1, 0, 2);
    TH1F * countWithMCProd = new TH1F("CountWithMCProd", "CountWithMCProd", 1, 0, 2);
    TH1F * countWithPUMCProd = new TH1F("CountWithPUMCProd", "CountWithPUMCProd", 1, 0, 2);
    TH1F * countWithSignalQCDcorrections = new TH1F("countWithSignalQCDcorrections", "countWithSignalQCDcorrections", 1, 0, 2);

    TH1F * htemp = 0;
    while ((chainElem = (TChainElement*) next())) {
        f2 = TFile::Open("dcache:" + TString(chainElem->GetTitle()));
        htemp = (TH1F*) f2->Get("Count");
        count->SetBinContent(1, count->GetBinContent(1)+htemp->GetBinContent(1));
        htemp = (TH1F*) f2->Get("CountWithPU");
        countWithPU->SetBinContent(2, countWithPU->GetBinContent(2)+htemp->GetBinContent(1));
        htemp = (TH1F*) f2->Get("CountWithPUP");
        countWithPUP->SetBinContent(3, countWithPUP->GetBinContent(3)+htemp->GetBinContent(1));
        htemp = (TH1F*) f2->Get("CountWithPUM");
        countWithPUM->SetBinContent(4, countWithPUM->GetBinContent(4)+htemp->GetBinContent(1));
        htemp = (TH1F*) f2->Get("CountWithMCProd");
        countWithMCProd->SetBinContent(5, countWithMCProd->GetBinContent(5)+htemp->GetBinContent(1));
        htemp = (TH1F*) f2->Get("CountWithPUMCProd");
        countWithPUMCProd->SetBinContent(6, countWithPUMCProd->GetBinContent(6)+htemp->GetBinContent(1));
        htemp = (TH1F*) f2->Get("countWithSignalQCDcorrections");
        countWithSignalQCDcorrections->SetBinContent(7, countWithSignalQCDcorrections->GetBinContent(7)+htemp->GetBinContent(1));

        std::clog << process << ": skimmed from " << chainElem->GetTitle() << std::endl;
    }

    // Make output directory if it doesn't exist
    if (gSystem->AccessPathName(outdir))
        gSystem->mkdir(outdir);
    TString outname = outdir + prefix + Form("%s.root", process.Data());

    TFile* f1 = TFile::Open(outname, "RECREATE");
    TTree* t1 = (TTree*) chain->CopyTree(selection);
    std::clog << process << ": skimmed from " << chain->GetEntriesFast() << " to " << t1->GetEntriesFast() << " entries." << std::endl;

    t1->Write();
    count->Write();
    countWithPU->Write();
    countWithPUP->Write();
    countWithPUM->Write();
    countWithMCProd->Write();
    countWithPUMCProd->Write();
    countWithSignalQCDcorrections->Write();
    f1->Close();

    return;
}

// To run:
//root -l -b -q Skim.C+\(\"ZnnH125\"\)

