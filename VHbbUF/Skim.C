#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TString.h"
#include "TCut.h"
#include "TH1.h"
#include "TSystem.h"
#include "TROOT.h"

#include "HelperNtuples.h"
//#include "HelperFunctions.h"

//#define XROOTD

void Skim(TString process="ZnnH125")
{
    gROOT->LoadMacro("HelperFunctions.h" );  // make functions visible to TTreeFormula
    gROOT->SetBatch(1);

    TChain * chain  = new TChain("tree");
    TString fname   = "";
    TString dijet   = "DiJetPt_";
#ifndef XROOTD
    //TString dirMC   = "dcache:/pnfs/cms/WAX/resilient/jiafu/ZnunuHbb/" + tagMC + "/";
    //TString dirData = "dcache:/pnfs/cms/WAX/resilient/jiafu/ZnunuHbb/" + tagData + "/";
    TString dirMC   = "dcache:/pnfs/cms/WAX/11/store/user/lpchbb/apana/" + tagMC + "/";
    TString dirData = "dcache:/pnfs/cms/WAX/11/store/user/lpchbb/apana/" + tagData + "/";
#else
    TString dirMC   = "root://xrootd.ba.infn.it//store/user/arizzi/" + tagMC + "/";
    TString dirData = "root://xrootd.ba.infn.it//store/user/arizzi/" + tagData + "/";
#endif
    TString outdir  = "skim/";
    TString prefix  = "skim_";
    TString suffix  = ".root";

    // ZbbHinv
    if (process == "ZbbHinv105") {
        //fname = "dcache:/pnfs/cms/WAX/11/store/user/lpchbb/degrutto/ZH_ZToBB_HToInv_M-105_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/degrutto/ZH_ZToBB_HToInv_M-105_8TeV_pythia6/ZH_ZToBB_HToInv_M-105_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1.root";
        fname = "/uscms_data/d3/degrutto/CMSSW_5_3_3_patch2/src/VHbbAnalysis/VHbbDataFormats/bin/ZH_ZToBB_HToInv_M-105_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1.root";
        chain->Add(fname);
    } else if (process == "ZbbHinv115") {
        //fname = "dcache:/pnfs/cms/WAX/11/store/user/lpchbb/degrutto/ZH_ZToBB_HToInv_M-115_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/degrutto/ZH_ZToBB_HToInv_M-115_8TeV_pythia6/ZH_ZToBB_HToInv_M-115_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1.root";
        fname = "/uscms_data/d3/degrutto/CMSSW_5_3_3_patch2/src/VHbbAnalysis/VHbbDataFormats/bin/ZH_ZToBB_HToInv_M-115_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1.root";
        chain->Add(fname);
    } else if (process == "ZbbHinv125") {
        //fname = "dcache:/pnfs/cms/WAX/11/store/user/lpchbb/degrutto/ZH_ZToBB_HToInv_M-125_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/degrutto/ZH_ZToBB_HToInv_M-125_8TeV_pythia6/ZH_ZToBB_HToInv_M-125_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1.root";
        fname = "/uscms_data/d3/degrutto/CMSSW_5_3_3_patch2/src/VHbbAnalysis/VHbbDataFormats/bin/ZH_ZToBB_HToInv_M-125_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1.root";
        chain->Add(fname);
    } else if (process == "ZbbHinv135") {
        //fname = "dcache:/pnfs/cms/WAX/11/store/user/lpchbb/degrutto/ZH_ZToBB_HToInv_M-135_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/degrutto/ZH_ZToBB_HToInv_M-135_8TeV_pythia6/ZH_ZToBB_HToInv_M-135_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1.root";
        fname = "/uscms_data/d3/degrutto/CMSSW_5_3_3_patch2/src/VHbbAnalysis/VHbbDataFormats/bin/ZH_ZToBB_HToInv_M-135_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1.root";
        chain->Add(fname);
    } else if (process == "ZbbHinv145") {
        //fname = "dcache:/pnfs/cms/WAX/11/store/user/lpchbb/degrutto/ZH_ZToBB_HToInv_M-145_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/degrutto/ZH_ZToBB_HToInv_M-145_8TeV_pythia6/ZH_ZToBB_HToInv_M-145_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1.root";
        fname = "/uscms_data/d3/degrutto/CMSSW_5_3_3_patch2/src/VHbbAnalysis/VHbbDataFormats/bin/ZH_ZToBB_HToInv_M-145_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1.root";
        chain->Add(fname);
    } else if (process == "ZbbHinv150") {
        //fname = "dcache:/pnfs/cms/WAX/11/store/user/lpchbb/degrutto/ZH_ZToBB_HToInv_M-150_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/degrutto/ZH_ZToBB_HToInv_M-150_8TeV_pythia6/ZH_ZToBB_HToInv_M-150_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1.root";
        fname = "/uscms_data/d3/degrutto/CMSSW_5_3_3_patch2/src/VHbbAnalysis/VHbbDataFormats/bin/ZH_ZToBB_HToInv_M-150_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1.root";
        chain->Add(fname);

    // ZnnH
    } else if (process == "ZnnH110") {
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
    } else if (process == "WJetsPtW70" ) {
        fname = dirMC + dijet + "WJetsToLNu_PtW-70To100_TuneZ2star_8TeV-madgraph" + suffix;
        chain->Add(fname);
    } else if (process == "WJetsPtW100") {
        fname = dirMC + dijet + "WJetsToLNu_PtW-100_TuneZ2star_8TeV-madgraph" + suffix;
        chain->Add(fname);
        fname = dirMC + dijet + "WJetsToLNu_PtW-100_TuneZ2star_8TeV_ext-madgraph-tarball" + suffix;
        chain->Add(fname);
    } else if (process == "WJetsPtW180") {
        fname = dirMC + dijet + "WJetsToLNu_PtW-180_TuneZ2star_8TeV-madgraph-tarball" + suffix;
        chain->Add(fname);
    } else if (process == "WJetsHW"    ) {
        fname = dirMC + dijet + "WJetsToLNu_PtW-100_8TeV-herwigpp" + suffix;
        chain->Add(fname);
    
    // ZJets
    } else if (process == "ZJetsPtZ70") {
        fname = dirMC + dijet + "ZJetsToNuNu_PtZ-70To100_8TeV" + suffix;
        chain->Add(fname);

    } else if (process == "ZJetsPtZ100") {
        //fname = dirMC + dijet + "ZJetsToNuNu_PtZ-100_8TeV-madgraph" + suffix;  // bugged
        fname = dirMC + dijet + "ZJetsToNuNu_PtZ-100_8TeV-madgraph" + suffix;
        chain->Add(fname);
        
    } else if (process == "ZJetsHT50"  ) {
        fname = dirMC + dijet + "ZJetsToNuNu_50_HT_100_TuneZ2Star_8TeV_madgraph" + suffix;
        chain->Add(fname);
    } else if (process == "ZJetsHT100" ) {
        fname = dirMC + dijet + "ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph" + suffix;
        chain->Add(fname);
    } else if (process == "ZJetsHT200" ) {
        fname = dirMC + dijet + "ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph" + suffix;
        chain->Add(fname);
    } else if (process == "ZJetsHT400" ) {
        fname = dirMC + dijet + "ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph" + suffix;
        chain->Add(fname);
    } else if (process == "ZJetsHW"    ) {
        fname = dirMC + dijet + "ZJetsToNuNu_Pt-100_8TeV-herwigpp" + suffix;
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
    } else if (process == "DYJetsZmmToZbb") {
        fname = "/uscms_data/d2/jiafu/VHbbAnalysis/NtupleV42_FOR_CVS/src/VHbbAnalysis/VHbbDataFormats/bin/TestZmmToZbb.root";
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
    
    } else if (process == "TTMCatNLO") {
        //fname = dirMC + dijet + "TT_8TeV-mcatnlo-withGENWEIGHT" + suffix;
        fname = dirMC + dijet + "TT_8TeV-mcatnlo" + suffix;
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
    
    // QCD
    } else if (process == "QCDPt50"  ) {
        fname = dirMC + dijet + "QCD_Pt-50to80_TuneZ2star_8TeV_pythia6" + suffix;
        chain->Add(fname);
    } else if (process == "QCDPt80"  ) {
        fname = dirMC + dijet + "QCD_Pt-80to120_TuneZ2star_8TeV_pythia6" + suffix;
        chain->Add(fname);
    } else if (process == "QCDPt120" ) {
        fname = dirMC + dijet + "QCD_Pt-120to170_TuneZ2star_8TeV_pythia6" + suffix;
        chain->Add(fname);
    } else if (process == "QCDPt170" ) {
        fname = dirMC + dijet + "QCD_Pt-170to300_TuneZ2star_8TeV_pythia6" + suffix;
        chain->Add(fname);
    } else if (process == "QCDPt300" ) {
        fname = dirMC + dijet + "QCD_Pt-300to470_TuneZ2star_8TeV_pythia6" + suffix;
        chain->Add(fname);
    } else if (process == "QCDPt470" ) {
        fname = dirMC + dijet + "QCD_Pt-470to600_TuneZ2star_8TeV_pythia6" + suffix;
        chain->Add(fname);
    } else if (process == "QCDPt600" ) {
        fname = dirMC + dijet + "QCD_Pt-600to800_TuneZ2star_8TeV_pythia6" + suffix;
        chain->Add(fname);
    } else if (process == "QCDPt800" ) {
        fname = dirMC + dijet + "QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6" + suffix;
        chain->Add(fname);
    } else if (process == "QCDPt1000") {
        fname = dirMC + dijet + "QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6" + suffix;
        chain->Add(fname);
    } else if (process == "QCDPt1400") {
        fname = dirMC + dijet + "QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6" + suffix;
        chain->Add(fname);
    } else if (process == "QCDPt1800") {
        fname = dirMC + dijet + "QCD_Pt-1800_TuneZ2star_8TeV_pythia6" + suffix;
        chain->Add(fname);
    } else if (process == "QCDHT100" ) {
        fname = dirMC + dijet + "QCD_HT-100To250_TuneZ2star_8TeV-madgraph-pythia" + suffix;
        chain->Add(fname);
    } else if (process == "QCDHT250" ) {
        fname = dirMC + dijet + "QCD_HT-250To500_TuneZ2star_8TeV-madgraph-pythia6" + suffix;
        chain->Add(fname);

    // Data
    } else if (process == "Data_R") {
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
    
    } else if (process == "Data_P") {
        fname = dirData + dijet + "METRun2012Cv2-PromptReco-v2EdmV42_v3" + suffix;  // Run 201191 has to be excluded
        chain->Add(fname);
        fname = dirData + dijet + "METRun2012D" + suffix;
        chain->Add(fname);
    
    } else if (process == "Data_METBTag_R") {
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
    
    } else if (process == "Data_METBTag_P") {
        fname = dirData + dijet + "METRun2012Cv2-PromptReco-v2EdmV42_v3" + suffix;  // Run 201191 has to be excluded
        chain->Add(fname);
        fname = dirData + dijet + "METRun2012D" + suffix;
        chain->Add(fname);
    
    } else if (process == "Data_SingleMu_R") {
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
    
    } else if (process == "Data_SingleMu_P") {
        fname = dirData + dijet + "SingleMuRun2012CPromptv2EdmV42" + suffix;  // Run 201191 has to be excluded
        chain->Add(fname);
        fname = dirData + dijet + "SingleMuRun2012CPromptV2TopUpEdmV42" + suffix;  // Run 201191 has to be excluded
        chain->Add(fname);
        fname = dirData + dijet + "SingleMuRun2012D-PromptReco-v1" + suffix;
        chain->Add(fname);
    
    } else if (process == "Data_SingleEl_R") {
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
    
    } else if (process == "Data_SingleEl_P") {
        fname = dirData + dijet + "SingleElectronRun2012CPromptv2EdmV42" + suffix;  // Run 201191 has to be excluded
        chain->Add(fname);
        fname = dirData + dijet + "SingleElectronRun2012CPromptV2TopUpEdmV42" + suffix;  // Run 201191 has to be excluded
        chain->Add(fname);
        fname = dirData + dijet + "SingleElectronRun2012D-PromptReco-v1_v3" + suffix;
        chain->Add(fname);

    } else {
        std::cout << "Error: Unrecognized process. I quit!" << std::endl;
        return;
    }

    TCut selection = baseline.c_str();
    // Different baseline for dimuons
    if (process == "DYJetsZmmToZbb") {
        selection = baselineZmmToZbb.c_str();
    }

    // MET filters
    selection += metfilter.c_str();

    // JSON & trigger
    if (process == "Data_R" || process == "Data_P") {
        TCut trigger = mettrigger.c_str();
        selection += trigger;
        //selection.Print();
    } else if (process == "Data_METBTag_R" || process == "Data_METBTag_P") {
        TCut trigger = metbtagtrigger.c_str();
        selection += trigger;
        //selection.Print();
    } else if (process == "Data_SingleMu_R" || process == "Data_SingleMu_P") {
        TCut trigger = singlemutrigger.c_str();
        selection += trigger;
        //selection.Print();
    } else if (process == "Data_SingleEl_R" || process == "Data_SingleEl_P") {
        TCut trigger = singleeltrigger.c_str();
        selection += trigger;
        //selection.Print();
    }
    
    // Exclude runs in Prompt data that have been reprocessed
    if (process == "Data_P" || process == "Data_METBTag_P" ||
        process == "Data_SingleMu_P" || process == "Data_SingleEl_P") {
        TCut run = "(EVENT.run != 201191)";
        //TCut run = "(EVENT.run != 201191) && (!(EVENT.run <= 207883 && EVENT.run <= 208307))";
        selection += run;
        //selection.Print();
    }
    
    // MC bugs
    TCut mcbug = "";
    if (process == "WJetsPtW180") {
        mcbug += "(!(180<=lheV_pt && lheV_pt<220 && lheNj==1))";
    }
    selection += mcbug;

    // Sum Count, CountWithPU, CountWithPUP, CountWithPUM
    TObjArray * files = chain->GetListOfFiles();
    TIter next(files);
    TChainElement * chainElem = 0;
    TFile * f2 = 0;
    TH1D * h1 = new TH1D("Count", ";Counts", 16, 0, 16);
    TH1F * htemp = 0;
    while ((chainElem = (TChainElement*) next())) {
#ifndef XROOTD
        f2 = TFile::Open("dcache:" + TString(chainElem->GetTitle()));
#else
        f2 = TFile::Open(TString(chainElem->GetTitle()));
#endif
        htemp = (TH1F*) f2->Get("Count");
        h1->SetBinContent(1, h1->GetBinContent(1)+htemp->GetBinContent(1));
        htemp = (TH1F*) f2->Get("CountWithPU");
        h1->SetBinContent(2, h1->GetBinContent(2)+htemp->GetBinContent(1));
        htemp = (TH1F*) f2->Get("CountWithPUP");
        h1->SetBinContent(3, h1->GetBinContent(3)+htemp->GetBinContent(1));
        htemp = (TH1F*) f2->Get("CountWithPUM");
        h1->SetBinContent(4, h1->GetBinContent(4)+htemp->GetBinContent(1));
        htemp = (TH1F*) f2->Get("CountWithMCProd");
        h1->SetBinContent(5, h1->GetBinContent(5)+htemp->GetBinContent(1));
        htemp = (TH1F*) f2->Get("CountWithPUMCProd");
        h1->SetBinContent(6, h1->GetBinContent(6)+htemp->GetBinContent(1));
        htemp = (TH1F*) f2->Get("countWithSignalQCDcorrections");
        h1->SetBinContent(7, h1->GetBinContent(7)+htemp->GetBinContent(1));
        
        std::clog << process << ": skimmed from " << chainElem->GetTitle() << std::endl;
    }
    
    // LHE Count
    TH1D * h2 = new TH1D("LHECount", ";LHE Counts", 16, 0, 16);
    TString process_lhe = process;
    if (process_lhe.BeginsWith("WJets"))
        process_lhe = "WJets";
    else if (process_lhe.BeginsWith("ZJets"))
        process_lhe = "ZJets";
    else 
        process_lhe = "";
    const std::vector<std::string>& lhecuts = GetLHECuts(process_lhe.Data());
    for (unsigned int i=0; i < lhecuts.size(); i++) {
        TCut cut2 = lhecuts.at(i).c_str();
        cut2 += mcbug;
        h2->SetBinContent(i+1, chain->GetEntries(cut2));
    }
    
    
    // Make output directory if it doesn't exist
    if (gSystem->AccessPathName(outdir))
        gSystem->mkdir(outdir);
    TString outname = outdir + prefix + Form("%s.root", process.Data());

    TFile* f1 = TFile::Open(outname, "RECREATE");
    TTree* t1 = (TTree*) chain->CopyTree(selection);
    std::clog << process << ": skimmed from " << chain->GetEntriesFast() << " to " << t1->GetEntriesFast() << " entries." << std::endl;
    
    t1->Write();
    h1->Write();
    h2->Write();
    f1->Close();

    return;
}

// To run:
//root -l -b -q Skim.C+\(\"ZnnH125\"\)

