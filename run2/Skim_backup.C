#include <iostream>
#include "TFile.h"
#include "TFileInfo.h"
#include "TTree.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TString.h"
#include "TCut.h"
#include "TH1.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFileCollection.h"

#include "HelperNtuples.h"


//#include "HelperFunctions.h"

//#define XROOTD

void Skim(TString process="ZnnH125")
{
    gROOT->LoadMacro("HelperFunctions.h" );  // make functions visible to TTreeFormula
    gROOT->SetBatch(1);

    TChain * chain  = new TChain("tree");
    TString fname   = "";
    TString dijet   = "";
#ifndef XROOTD
    //TString dirMC   = "dcache:/pnfs/cms/WAX/resilient/jiafu/ZnunuHbb/" + tagMC + "/";
    //TString dirData = "dcache:/pnfs/cms/WAX/resilient/jiafu/ZnunuHbb/" + tagData + "/";
    //    TString dirMC   = "dcache:/pnfs/cms/WAX/11/store/user/lpchbb/apana/" + tagMC + "/";
    //    TString dirData = "dcache:/pnfs/cms/WAX/11/store/user/lpchbb/apana/" + tagData + "/";
    //tagMC         = "root://eoscms//eos/cms//store/user/capalmer/VHBBHeppyNtuples/V7/"

    //    TString dirMC   = "root://eoscms//eos/cms/store/user/capalmer/VHBBHeppyNtuples/V10/VHBBHeppyV10_";
    TString dirMC   = "root://xrootd.unl.edu//store/user/arizzi/VHBBHeppyV11/";
    TString dirData = "dcache:/pnfs/cms/WAX/11/store/user/lpchbb/apana/" + tagData + "/";
#else
    TString dirMC   = "root://xrootd.ba.infn.it//store/user/arizzi/" + tagMC + "/";
    TString dirData = "root://xrootd.ba.infn.it//store/user/arizzi/" + tagData + "/";
#endif
    TString outdir  = " /afs/cern.ch/work/d/degrutto/public/MiniAOD/ZnnHbb_Phys14_PU20bx25/skimV11/";
    TString prefix  = "skim_";
    TString suffix  = "tree*.root";
    

    // ZnnH
    if (process == "ZnnH125") {
  ///store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/150127_174533/0000/
      fname = dirMC + dijet + "ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp_VHBB_HEPPY_V11_ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp__Phys14DR-PU40bx25_PHYS14_25_V1-v1_150408_085844.root";
//ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp_VHBB_HEPPY_V10_ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1.root";

        chain->Add(fname);
    // WlnH
    } else if (process == "WlnH125") {
      fname = dirMC + dijet + "WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp_VHBB_HEPPY_V11_WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_150408_085753.root"; 
	 //"WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp_VHBB_HEPPY_V10_WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1.root";
        chain->Add(fname);
    // monoHiggs
    } else if (process == "monoH") {
        fname = "root://eoscms//eos/cms/store/user/degrutto/VHBBHeppyNtuplesV10/monohiggsbb.root";
        chain->Add(fname);
    // WJets
    } else if (process == "WJetsIncl" ) {
      fname = dirMC + dijet +"WJetsToLNu_13TeV-madgraph-pythia8-tauola_VHBB_HEPPY_V11_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1_150408_090322.root"; 
	  //"WJetsToLNu_13TeV-madgraph-pythia8-tauola_VHBB_HEPPY_V10_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1.root";
        chain->Add(fname);
    } else if (process == "WJetsHT100" ) {
      fname = dirMC + dijet + "WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola_VHBB_HEPPY_V11_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1_150408_090336.root"; 
	  //"WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola_VHBB_HEPPY_V10_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1.root";
        chain->Add(fname);
    } else if (process == "WJetsHT200") {
      fname = dirMC + dijet + "WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola_VHBB_HEPPY_V11_WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1_150408_090358.root";
//"WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola_VHBB_HEPPY_V10_WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1.root";
        chain->Add(fname);  

    } else if (process == "WJetsHT400") {
      fname = dirMC + dijet + "WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola_VHBB_HEPPY_V11_WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1_150408_090418.root";
//"WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola_VHBB_HEPPY_V10_WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1.root";        
	chain->Add(fname);
  
   } else if (process == "WJetsHT600") {
  //... missing
      fname = dirMC + dijet + "WJetsToLNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola_VHBB_HEPPY_V11_WJetsToLNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1_150408_090435.root"; 
	//"WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174711/0000/" + suffix;        chain->Add(fname);
  	chain->Add(fname);
     
    // ZJets
    } else if (process == "ZJetsHT100" ) {
      fname = dirMC + dijet + "WJetsToLNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola_VHBB_HEPPY_V11_WJetsToLNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1_150408_090435.root";
//"ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola_VHBB_HEPPY_V10_ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1.root";
        chain->Add(fname);
    } else if (process == "ZJetsHT200" ) {
      fname = dirMC + dijet + "ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola_VHBB_HEPPY_V11_ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1_150408_092839.root";
//"ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola_VHBB_HEPPY_V10_ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1.root";
        chain->Add(fname);
 } else if (process == "ZJetsHT400" ) {  // missing
      fname = dirMC + dijet + "ZJetsToNuNu_HT-400to600_Tune4C_13TeV-madgraph-tauola_VHBB_HEPPY_V11_ZJetsToNuNu_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v2_150409_121410.root";
"ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph" + suffix;
        chain->Add(fname);
 } else if (process == "ZJetsHT600" ) {  
      fname = dirMC + dijet + "ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola_VHBB_HEPPY_V11_ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1_150410_064455.root";
//"ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola_VHBB_HEPPY_V10_ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1.root";
        chain->Add(fname);
    
     // VV  // missing :(
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
    } else if (process == "TTPythia8") {
        fname = dirMC + dijet + "TT_Tune4C_13TeV-pythia8-tauola_VHBB_HEPPY_V11_TT_Tune4C_13TeV-pythia8-tauola__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_150409_121154.root";
        chain->Add(fname);

    } else if (process == "TTMad") {

      fname = dirMC + dijet + "TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola_VHBB_HEPPY_V11_TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1_150408_102615.root";
//"TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola_VHBB_HEPPY_V10_TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1.root";
      chain->Add(fname);
      //      TFileCollection fc("dum", "", "myfilelist.txt"); // The name is irrelevant
      //      fc.Add(fname);
      // chain->AddFileInfoList((TCollection*) fc.GetList() );


    // T/Tbar
    
    } else if (process == "TTMadv2") {


      //      fname = dirMC + dijet + "TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola_VHBB_HEPPY_V11_TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1_150408_102615.root";
//"TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola_VHBB_HEPPY_V10_TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1.root";
      //  chain->Add(fname);
	   TFileCollection* fc = new TFileCollection("mylist", "mylist", "myTTMadfilelist.txt");
	   chain->AddFileInfoList((TCollection*)fc->GetList());

	   //  fc.Add(fname);



    // T/Tbar
    }


 else if (process == "T_tW"   ) {
      fname = dirMC + dijet + "T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola_VHBB_HEPPY_V11_T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1_150408_102551.root";
//"T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola_VHBB_HEPPY_V10_T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1.root";
        chain->Add(fname);
    } else if (process == "Tbar_tW") {
      fname = dirMC + dijet + "Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola_VHBB_HEPPY_V11_Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1_150408_102531.root";
//"Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola_VHBB_HEPPY_V10_Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1.root";
        chain->Add(fname);
    } else if (process == "T_s"    ) {
      fname = dirMC + dijet + "TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola_VHBB_HEPPY_V11_TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1_150409_121235.root";//"TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola_VHBB_HEPPY_V10_TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1.root";
        chain->Add(fname);
    } else if (process == "Tbar_s" ) {
      fname = dirMC + dijet + "TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola_VHBB_HEPPY_V11_TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1_150409_121122.root";
//TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola_VHBB_HEPPY_V10_TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1.root";
        chain->Add(fname);
    } else if (process == "T_t"    ) {
      fname = dirMC + dijet + "TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola_VHBB_HEPPY_V11_TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1_150409_121347.root";
//"TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola_VHBB_HEPPY_V10_TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1.root";
        chain->Add(fname);
    } else if (process == "Tbar_t" ) {
      fname = dirMC + dijet + "TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola_VHBB_HEPPY_V11_TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1_150409_121140.root";
//"TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola_VHBB_HEPPY_V10_TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1.root";
        chain->Add(fname);
    
    // QCD
    } else if (process == "QCDHT100" ) {
      fname = dirMC + dijet + "QCD_HT-100To250_13TeV-madgraph_VHBB_HEPPY_V11_QCD_HT-100To250_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1_150408_090013.root";
//"QCD_HT-100To250_13TeV-madgraph_VHBB_HEPPY_V10_QCD_HT-100To250_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1.root";
        chain->Add(fname);
    } else if (process == "QCDHT250" ) {

      fname = dirMC + dijet + "QCD_HT_250To500_13TeV-madgraph_VHBB_HEPPY_V11_QCD_HT_250To500_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1_150408_090046.root";
//;"QCD_HT_250To500_13TeV-madgraph_VHBB_HEPPY_V10_QCD_HT_250To500_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v2.root";
        chain->Add(fname);
	fname = dirMC + dijet + "QCD_HT_250To500_13TeV-madgraph_VHBB_HEPPY_V11_QCD_HT_250To500_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v2_150408_090120.root";
//QCD_HT_250To500_13TeV-madgraph_VHBB_HEPPY_V10_QCD_HT_250To500_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1.root";
	 chain->Add(fname);
    } else if (process == "QCDHT500" ) {
      fname = dirMC + dijet + "QCD_HT-500To1000_13TeV-madgraph_VHBB_HEPPY_V11_QCD_HT-500To1000_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v1_150408_090150.root";
//"QCD_HT-500To1000_13TeV-madgraph_VHBB_HEPPY_V10_QCD_HT-500To1000_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1.root";
    chain->Add(fname);
    fname = dirMC + dijet + "QCD_HT-500To1000_13TeV-madgraph_VHBB_HEPPY_V11_QCD_HT-500To1000_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1_150408_090136.root";
//QCD_HT-500To1000_13TeV-madgraph_VHBB_HEPPY_V10_QCD_HT-500To1000_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v1.root";
    chain->Add(fname);
    }
 else if (process == "QCDHT1000" ) {
   fname = dirMC + dijet + "QCD_HT_1000ToInf_13TeV-madgraph_VHBB_HEPPY_V11_QCD_HT_1000ToInf_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1_150408_090301.root";
//"QCD_HT_1000ToInf_13TeV-madgraph_VHBB_HEPPY_V10_QCD_HT_1000ToInf_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v1.root";
    chain->Add(fname);
    //    fname = dirMC + dijet + "";
//"QCD_HT_1000ToInf_13TeV-madgraph_VHBB_HEPPY_V10_QCD_HT_1000ToInf_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1.root";
//    chain->Add(fname);
    // Data
    }
 else if (process == "Data_R") {
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

    // Sum Count, CountWithPU, CountWithPUP, CountWithPUM
    TObjArray * files = chain->GetListOfFiles();
    TIter next(files);
    TChainElement * chainElem = 0;
    TFile * f2 = 0;
    TH1D * h1 = new TH1D("Count", ";Counts", 16, 0, 16);
    TH1F * htemp = 0;

    while ((chainElem = (TChainElement*) next())) {
      //#ifndef XROOTD
      //      f2 = TFile::Open("dcache:" + TString(chainElem->GetTitle()));
      //#else
      std::cout << "chainElem->GetTitle() " << chainElem->GetTitle() << std::endl; 
        f2 = TFile::Open( TString(chainElem->GetTitle()));
	//#endif
        htemp = (TH1F*) f2->Get("Count");
        h1->SetBinContent(1, h1->GetBinContent(1)+htemp->GetBinContent(1));
	/*
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
        */
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
/*
cmsLs -R /store/user/capalmer/VHBBHeppyNtuples/V7/ >> Skim.C
dr-x(019)            1 2015-01-30 16:15:10 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-100To250_13TeV-madgraph
dr-x(019)            2 2015-01-30 16:15:26 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-500To1000_13TeV-madgraph
dr-x(019)            2 2015-01-30 16:15:08 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_1000ToInf_13TeV-madgraph
dr-x(019)            2 2015-01-30 16:15:22 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_250To500_13TeV-madgraph
dr-x(019)            2 2015-01-30 16:15:08 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8
dr-x(019)            2 2015-01-30 16:15:23 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8
dr-x(019)            3 2015-01-30 16:14:37 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8
dr-x(019)            5 2015-01-30 16:14:52 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8
dr-x(019)            3 2015-01-30 16:15:14 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8
dr-x(019)            2 2015-01-30 16:14:46 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8
dr-x(019)            1 2015-01-30 16:15:55 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8
dr-x(019)            2 2015-01-30 16:15:10 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8
dr-x(019)            1 2015-01-30 16:15:55 /store/user/capalmer/VHBBHeppyNtuples/V7/TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola
dr-x(019)            1 2015-01-30 16:15:28 /store/user/capalmer/VHBBHeppyNtuples/V7/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola
dr-x(019)            1 2015-01-30 16:16:02 /store/user/capalmer/VHBBHeppyNtuples/V7/TTWJets_Tune4C_13TeV-madgraph-tauola
dr-x(019)            1 2015-01-30 16:14:35 /store/user/capalmer/VHBBHeppyNtuples/V7/TTZJets_Tune4C_13TeV-madgraph-tauola
dr-x(019)            1 2015-01-30 16:15:56 /store/user/capalmer/VHBBHeppyNtuples/V7/TT_Tune4C_13TeV-pythia8-tauola
dr-x(019)            1 2015-01-30 16:14:47 /store/user/capalmer/VHBBHeppyNtuples/V7/TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola
dr-x(019)            1 2015-01-30 16:15:59 /store/user/capalmer/VHBBHeppyNtuples/V7/TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola
dr-x(019)            1 2015-01-30 16:14:45 /store/user/capalmer/VHBBHeppyNtuples/V7/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola
dr-x(019)            1 2015-01-30 16:14:42 /store/user/capalmer/VHBBHeppyNtuples/V7/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola
dr-x(019)            3 2015-01-30 16:15:11 /store/user/capalmer/VHBBHeppyNtuples/V7/WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp
dr-x(019)            1 2015-01-30 16:14:42 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola
dr-x(019)            1 2015-01-30 16:14:33 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola
dr-x(019)            1 2015-01-30 16:15:45 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola
dr-x(019)            1 2015-01-30 16:15:49 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola
dr-x(019)            1 2015-01-30 16:15:53 /store/user/capalmer/VHBBHeppyNtuples/V7/WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola
dr-x(019)            3 2015-01-30 21:18:20 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp
dr-x(019)            3 2015-01-30 16:15:09 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp
dr-x(019)            1 2015-01-30 16:15:35 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola
dr-x(019)            1 2015-01-30 16:15:51 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola
dr-x(019)            1 2015-01-30 16:15:13 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola
dr-x(019)            1 2015-01-30 16:14:39 /store/user/capalmer/VHBBHeppyNtuples/V7/ZZTo4L_Tune4C_13TeV-powheg-pythia8
dr-x(019)            1 2015-01-30 16:14:42 /store/user/capalmer/VHBBHeppyNtuples/V7/ZZTo4L_Tune4C_13TeV-powheg-pythia8/VHBB_HEPPY_V7_ZZTo4L_Tune4C_13TeV-powheg-pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:14:43 /store/user/capalmer/VHBBHeppyNtuples/V7/ZZTo4L_Tune4C_13TeV-powheg-pythia8/VHBB_HEPPY_V7_ZZTo4L_Tune4C_13TeV-powheg-pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175140
dr-x(019)            0 2015-01-30 22:10:53 /store/user/capalmer/VHBBHeppyNtuples/V7/ZZTo4L_Tune4C_13TeV-powheg-pythia8/VHBB_HEPPY_V7_ZZTo4L_Tune4C_13TeV-powheg-pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175140/0000
-r--(016)     34505330 2015-01-30 22:10:53 /store/user/capalmer/VHBBHeppyNtuples/V7/ZZTo4L_Tune4C_13TeV-powheg-pythia8/VHBB_HEPPY_V7_ZZTo4L_Tune4C_13TeV-powheg-pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175140/0000/tree_1.root
-r--(016)     30880909 2015-01-30 22:10:50 /store/user/capalmer/VHBBHeppyNtuples/V7/ZZTo4L_Tune4C_13TeV-powheg-pythia8/VHBB_HEPPY_V7_ZZTo4L_Tune4C_13TeV-powheg-pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175140/0000/tree_2.root
-r--(016)     34035877 2015-01-30 22:10:51 /store/user/capalmer/VHBBHeppyNtuples/V7/ZZTo4L_Tune4C_13TeV-powheg-pythia8/VHBB_HEPPY_V7_ZZTo4L_Tune4C_13TeV-powheg-pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175140/0000/tree_3.root
-r--(016)     34853919 2015-01-30 22:10:51 /store/user/capalmer/VHBBHeppyNtuples/V7/ZZTo4L_Tune4C_13TeV-powheg-pythia8/VHBB_HEPPY_V7_ZZTo4L_Tune4C_13TeV-powheg-pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175140/0000/tree_4.root
-r--(016)      3547605 2015-01-30 22:10:47 /store/user/capalmer/VHBBHeppyNtuples/V7/ZZTo4L_Tune4C_13TeV-powheg-pythia8/VHBB_HEPPY_V7_ZZTo4L_Tune4C_13TeV-powheg-pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175140/0000/tree_5.root
dr-x(019)            1 2015-01-30 16:15:15 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:17 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174929
dr-x(019)            0 2015-01-30 22:13:10 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174929/0000
-r--(016)    556000765 2015-01-30 22:11:43 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174929/0000/tree_1.root
-r--(016)    412459605 2015-01-30 22:11:26 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174929/0000/tree_10.root
-r--(016)    130136116 2015-01-30 22:11:17 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174929/0000/tree_11.root
-r--(016)    548175680 2015-01-30 22:13:10 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174929/0000/tree_2.root
-r--(016)    526936433 2015-01-30 22:12:03 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174929/0000/tree_3.root
-r--(016)    419606040 2015-01-30 22:11:20 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174929/0000/tree_4.root
-r--(016)    578243291 2015-01-30 22:11:30 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174929/0000/tree_5.root
-r--(016)    488384758 2015-01-30 22:11:19 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174929/0000/tree_6.root
-r--(016)    485183110 2015-01-30 22:11:46 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174929/0000/tree_7.root
-r--(016)    541960517 2015-01-30 22:12:40 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174929/0000/tree_8.root
-r--(016)    519222114 2015-01-30 22:11:23 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174929/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:15:53 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:54 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174920
dr-x(019)            0 2015-01-30 22:13:24 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174920/0000
-r--(016)    274178565 2015-01-30 22:12:29 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174920/0000/tree_1.root
-r--(016)    276792027 2015-01-30 22:11:46 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174920/0000/tree_10.root
-r--(016)    335110248 2015-01-30 22:13:24 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174920/0000/tree_2.root
-r--(016)    335927690 2015-01-30 22:12:19 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174920/0000/tree_3.root
-r--(016)    335765086 2015-01-30 22:13:06 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174920/0000/tree_4.root
-r--(016)    335591379 2015-01-30 22:13:00 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174920/0000/tree_5.root
-r--(016)    303611478 2015-01-30 22:13:04 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174920/0000/tree_6.root
-r--(016)    291970112 2015-01-30 22:12:11 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174920/0000/tree_7.root
-r--(016)    315523711 2015-01-30 22:13:15 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174920/0000/tree_8.root
-r--(016)    254027286 2015-01-30 22:12:44 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174920/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:15:36 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:37 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174910
dr-x(019)            0 2015-01-30 22:12:20 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174910/0000
-r--(016)    125433904 2015-01-30 22:12:00 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174910/0000/tree_1.root
-r--(016)    137267969 2015-01-30 22:11:37 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174910/0000/tree_10.root
-r--(016)     55726047 2015-01-30 22:11:48 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174910/0000/tree_11.root
-r--(016)    137141694 2015-01-30 22:11:48 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174910/0000/tree_2.root
-r--(016)    137573090 2015-01-30 22:12:00 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174910/0000/tree_3.root
-r--(016)    129848634 2015-01-30 22:11:36 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174910/0000/tree_4.root
-r--(016)    137151126 2015-01-30 22:12:06 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174910/0000/tree_5.root
-r--(016)    138261103 2015-01-30 22:11:37 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174910/0000/tree_6.root
-r--(016)    131054529 2015-01-30 22:12:12 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174910/0000/tree_7.root
-r--(016)    124909738 2015-01-30 22:12:20 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174910/0000/tree_8.root
-r--(016)    125531909 2015-01-30 22:11:43 /store/user/capalmer/VHBBHeppyNtuples/V7/ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174910/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:15:09 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1
dr-x(019)            1 2015-01-30 16:15:11 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:11 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp__Phys14DR-PU40bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:13 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp__Phys14DR-PU40bx25_PHYS14_25_V1-v1/150127_174523
dr-x(019)            0 2015-02-02 15:26:29 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp__Phys14DR-PU40bx25_PHYS14_25_V1-v1/150127_174523/0000
-r--(016)     52179534 2015-01-30 22:11:11 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp__Phys14DR-PU40bx25_PHYS14_25_V1-v1/150127_174523/0000/tree_1.root
dr-x(019)            1 2015-01-30 16:15:13 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/150127_174533
dr-x(019)            0 2015-02-04 15:41:08 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/150127_174533/0000
-r--(016)     47670462 2015-01-30 22:11:12 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/150127_174533/0000/tree_1.root
dr-x(019)            1 2015-01-30 16:15:11 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_174514
dr-x(019)            0 2015-01-30 22:11:17 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_174514/0000
-r--(016)     88126688 2015-01-30 22:11:17 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_174514/0000/tree_1.root
dr-x(019)            1 2015-01-30 16:14:39 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 21:18:21 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp__Phys14DR-PU40bx25_tsg_PHYS14_25_V1-v2
dr-x(019)            1 2015-01-30 21:18:21 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp__Spring14miniaod-141029_PU40bx50_PLS170_V6AN2-v1
dr-x(019)            1 2015-01-30 21:18:22 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp__Spring14miniaod-141029_PU40bx50_PLS170_V6AN2-v1/150127_174442
dr-x(019)            0 2015-01-30 22:11:07 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp__Spring14miniaod-141029_PU40bx50_PLS170_V6AN2-v1/150127_174442/0000
-r--(016)    174683982 2015-01-30 22:11:07 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp__Spring14miniaod-141029_PU40bx50_PLS170_V6AN2-v1/150127_174442/0000/tree_1.root
dr-x(019)            1 2015-01-30 21:18:22 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp__Phys14DR-PU40bx25_tsg_PHYS14_25_V1-v2/150130_155900
dr-x(019)            0 2015-01-30 22:11:11 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp__Phys14DR-PU40bx25_tsg_PHYS14_25_V1-v2/150130_155900/0000
-r--(016)    173942668 2015-01-30 22:11:11 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp__Phys14DR-PU40bx25_tsg_PHYS14_25_V1-v2/150130_155900/0000/tree_1.root
dr-x(019)            1 2015-01-30 16:14:41 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/150127_174425
dr-x(019)            0 2015-01-30 22:11:26 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/150127_174425/0000
-r--(016)    142448039 2015-01-30 22:11:26 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/150127_174425/0000/tree_1.root
-r--(016)     16244614 2015-01-30 22:10:47 /store/user/capalmer/VHBBHeppyNtuples/V7/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/150127_174425/0000/tree_2.root
dr-x(019)            1 2015-01-30 16:15:55 /store/user/capalmer/VHBBHeppyNtuples/V7/WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:56 /store/user/capalmer/VHBBHeppyNtuples/V7/WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175116
dr-x(019)            0 2015-01-30 22:11:47 /store/user/capalmer/VHBBHeppyNtuples/V7/WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175116/0000
-r--(016)     87442906 2015-01-30 22:11:47 /store/user/capalmer/VHBBHeppyNtuples/V7/WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175116/0000/tree_1.root
dr-x(019)            1 2015-01-30 16:15:51 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:53 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174711
dr-x(019)            0 2015-01-30 22:13:30 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174711/0000
-r--(016)    611611096 2015-01-30 22:13:30 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174711/0000/tree_1.root
-r--(016)    426728293 2015-01-30 22:12:26 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174711/0000/tree_10.root
-r--(016)    553839861 2015-01-30 22:11:52 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174711/0000/tree_2.root
-r--(016)    574775934 2015-01-30 22:11:52 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174711/0000/tree_3.root
-r--(016)    568976497 2015-01-30 22:13:25 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174711/0000/tree_4.root
-r--(016)    610315322 2015-01-30 22:13:30 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174711/0000/tree_5.root
-r--(016)    551558967 2015-01-30 22:11:53 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174711/0000/tree_6.root
-r--(016)    608993992 2015-01-30 22:13:29 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174711/0000/tree_7.root
-r--(016)    602221981 2015-01-30 22:12:16 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174711/0000/tree_8.root
-r--(016)    586517087 2015-01-30 22:13:30 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174711/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:15:48 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:49 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174703
dr-x(019)            0 2015-01-30 23:59:39 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174703/0000
-r--(016)    397722983 2015-01-30 22:12:49 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174703/0000/tree_1.root
-r--(016)    394078951 2015-01-30 22:12:15 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174703/0000/tree_10.root
-r--(016)    309309750 2015-01-30 22:12:53 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174703/0000/tree_11.root
-r--(016)    429659934 2015-01-30 22:14:29 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174703/0000/tree_2.root
-r--(016)    424661554 2015-01-30 22:12:31 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174703/0000/tree_3.root
-r--(016)    369080550 2015-01-30 22:13:06 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174703/0000/tree_4.root
-r--(016)    421403556 2015-01-30 23:59:39 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174703/0000/tree_5.root
-r--(016)    397818920 2015-01-30 22:12:51 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174703/0000/tree_6.root
-r--(016)    438162912 2015-01-30 22:12:42 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174703/0000/tree_7.root
-r--(016)    346379745 2015-01-30 22:11:48 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174703/0000/tree_8.root
-r--(016)    440032669 2015-01-30 22:12:18 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174703/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:14:35 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:14:39 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174655
dr-x(019)            0 2015-01-30 22:11:22 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174655/0000
-r--(016)    196656966 2015-01-30 22:11:00 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174655/0000/tree_1.root
-r--(016)    135337956 2015-01-30 22:11:01 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174655/0000/tree_10.root
-r--(016)    118796465 2015-01-30 22:10:52 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174655/0000/tree_11.root
-r--(016)    127850645 2015-01-30 22:10:56 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174655/0000/tree_12.root
-r--(016)     59215573 2015-01-30 22:10:49 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174655/0000/tree_13.root
-r--(016)    101700528 2015-01-30 22:10:49 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174655/0000/tree_14.root
-r--(016)     81887735 2015-01-30 22:11:10 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174655/0000/tree_15.root
-r--(016)     95648125 2015-01-30 22:11:00 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174655/0000/tree_16.root
-r--(016)    161617127 2015-01-30 22:11:22 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174655/0000/tree_17.root
-r--(016)     33980732 2015-01-30 22:10:54 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174655/0000/tree_18.root
-r--(016)    118272721 2015-01-30 22:11:19 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174655/0000/tree_19.root
-r--(016)     97695564 2015-01-30 22:11:08 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174655/0000/tree_2.root
-r--(016)     90209253 2015-01-30 22:10:55 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174655/0000/tree_20.root
-r--(016)    114458168 2015-01-30 22:10:50 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174655/0000/tree_21.root
-r--(016)     25425423 2015-01-30 22:10:53 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174655/0000/tree_22.root
-r--(016)    124285362 2015-01-30 22:11:10 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174655/0000/tree_3.root
-r--(016)    121225633 2015-01-30 22:10:49 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174655/0000/tree_4.root
-r--(016)    119111054 2015-01-30 22:10:53 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174655/0000/tree_5.root
-r--(016)     89701985 2015-01-30 22:11:07 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174655/0000/tree_6.root
-r--(016)     85756993 2015-01-30 22:10:50 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174655/0000/tree_7.root
-r--(016)    112398766 2015-01-30 22:11:07 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174655/0000/tree_8.root
-r--(016)    124608929 2015-01-30 22:11:01 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174655/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:14:43 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:14:46 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646
dr-x(019)            0 2015-02-03 17:17:07 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646/0000
-r--(016)     21135645 2015-01-30 22:10:53 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646/0000/tree_1.root
-r--(016)     18563620 2015-01-30 22:10:57 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646/0000/tree_10.root
-r--(016)     20820186 2015-01-30 22:10:54 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646/0000/tree_11.root
-r--(016)     20808732 2015-01-30 22:10:54 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646/0000/tree_12.root
-r--(016)     22745001 2015-01-30 22:10:55 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646/0000/tree_13.root
-r--(016)     18860752 2015-01-30 22:10:53 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646/0000/tree_14.root
-r--(016)     22744455 2015-01-30 22:10:55 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646/0000/tree_15.root
-r--(016)     22615223 2015-01-30 22:10:57 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646/0000/tree_16.root
-r--(016)     19142233 2015-01-30 22:10:53 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646/0000/tree_17.root
-r--(016)     22726235 2015-01-30 22:10:58 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646/0000/tree_18.root
-r--(016)     20516463 2015-01-30 22:10:53 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646/0000/tree_19.root
-r--(016)     18267052 2015-01-30 22:10:54 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646/0000/tree_2.root
-r--(016)     18945574 2015-01-30 22:10:53 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646/0000/tree_20.root
-r--(016)     19670268 2015-01-30 22:10:55 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646/0000/tree_21.root
-r--(016)     19562092 2015-01-30 22:10:58 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646/0000/tree_22.root
-r--(016)      4648169 2015-01-30 22:10:52 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646/0000/tree_23.root
-r--(016)     22763126 2015-01-30 22:10:54 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646/0000/tree_3.root
-r--(016)     18429794 2015-01-30 22:10:55 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646/0000/tree_4.root
-r--(016)     20332183 2015-01-30 22:10:55 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646/0000/tree_5.root
-r--(016)     21561374 2015-01-30 22:10:59 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646/0000/tree_6.root
-r--(016)     17120669 2015-01-30 22:10:53 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646/0000/tree_7.root
-r--(016)     22524063 2015-01-30 22:10:54 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646/0000/tree_8.root
-r--(016)     18194641 2015-01-30 22:10:53 /store/user/capalmer/VHBBHeppyNtuples/V7/WJetsToLNu_13TeV-madgraph-pythia8-tauola/VHBB_HEPPY_V7_WJetsToLNu_13TeV-madgraph-pythia8-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174646/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:15:11 /store/user/capalmer/VHBBHeppyNtuples/V7/WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v2
dr-x(019)            1 2015-01-30 16:15:13 /store/user/capalmer/VHBBHeppyNtuples/V7/WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:11 /store/user/capalmer/VHBBHeppyNtuples/V7/WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp__Phys14DR-PU40bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:12 /store/user/capalmer/VHBBHeppyNtuples/V7/WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp__Phys14DR-PU40bx25_PHYS14_25_V1-v1/150127_174506
dr-x(019)            0 2015-02-02 15:22:34 /store/user/capalmer/VHBBHeppyNtuples/V7/WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp__Phys14DR-PU40bx25_PHYS14_25_V1-v1/150127_174506/0000
-r--(016)     73479053 2015-01-30 22:11:13 /store/user/capalmer/VHBBHeppyNtuples/V7/WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp__Phys14DR-PU40bx25_PHYS14_25_V1-v1/150127_174506/0000/tree_1.root
dr-x(019)            1 2015-01-30 16:15:15 /store/user/capalmer/VHBBHeppyNtuples/V7/WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/150127_174458
dr-x(019)            0 2015-02-04 15:43:01 /store/user/capalmer/VHBBHeppyNtuples/V7/WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/150127_174458/0000
-r--(016)     67245486 2015-01-30 22:11:09 /store/user/capalmer/VHBBHeppyNtuples/V7/WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/150127_174458/0000/tree_1.root
dr-x(019)            1 2015-01-30 16:15:13 /store/user/capalmer/VHBBHeppyNtuples/V7/WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v2/150127_174448
dr-x(019)            0 2015-01-30 22:11:08 /store/user/capalmer/VHBBHeppyNtuples/V7/WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v2/150127_174448/0000
-r--(016)    127574687 2015-01-30 22:11:08 /store/user/capalmer/VHBBHeppyNtuples/V7/WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp/VHBB_HEPPY_V7_WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v2/150127_174448/0000/tree_1.root
dr-x(019)            1 2015-01-30 16:14:43 /store/user/capalmer/VHBBHeppyNtuples/V7/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/VHBB_HEPPY_V7_Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:14:44 /store/user/capalmer/VHBBHeppyNtuples/V7/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/VHBB_HEPPY_V7_Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175022
dr-x(019)            0 2015-01-30 22:11:45 /store/user/capalmer/VHBBHeppyNtuples/V7/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/VHBB_HEPPY_V7_Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175022/0000
-r--(016)    365897176 2015-01-30 22:11:02 /store/user/capalmer/VHBBHeppyNtuples/V7/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/VHBB_HEPPY_V7_Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175022/0000/tree_1.root
-r--(016)    359420426 2015-01-30 22:11:45 /store/user/capalmer/VHBBHeppyNtuples/V7/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/VHBB_HEPPY_V7_Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175022/0000/tree_2.root
-r--(016)     40350126 2015-01-30 22:10:56 /store/user/capalmer/VHBBHeppyNtuples/V7/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/VHBB_HEPPY_V7_Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175022/0000/tree_3.root
dr-x(019)            1 2015-01-30 16:14:47 /store/user/capalmer/VHBBHeppyNtuples/V7/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/VHBB_HEPPY_V7_T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:14:49 /store/user/capalmer/VHBBHeppyNtuples/V7/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/VHBB_HEPPY_V7_T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175031
dr-x(019)            0 2015-01-30 22:11:09 /store/user/capalmer/VHBBHeppyNtuples/V7/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/VHBB_HEPPY_V7_T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175031/0000
-r--(016)    345097451 2015-01-30 22:11:09 /store/user/capalmer/VHBBHeppyNtuples/V7/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/VHBB_HEPPY_V7_T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175031/0000/tree_1.root
-r--(016)    387840526 2015-01-30 22:11:01 /store/user/capalmer/VHBBHeppyNtuples/V7/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/VHBB_HEPPY_V7_T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175031/0000/tree_2.root
-r--(016)     40434238 2015-01-30 22:11:00 /store/user/capalmer/VHBBHeppyNtuples/V7/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/VHBB_HEPPY_V7_T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175031/0000/tree_3.root
dr-x(019)            1 2015-01-30 16:16:00 /store/user/capalmer/VHBBHeppyNtuples/V7/TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:16:02 /store/user/capalmer/VHBBHeppyNtuples/V7/TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174946
dr-x(019)            0 2015-02-04 13:50:48 /store/user/capalmer/VHBBHeppyNtuples/V7/TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174946/0000
-r--(016)    230129336 2015-01-30 22:12:04 /store/user/capalmer/VHBBHeppyNtuples/V7/TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174946/0000/tree_1.root
-r--(016)    274925017 2015-01-30 22:12:22 /store/user/capalmer/VHBBHeppyNtuples/V7/TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174946/0000/tree_2.root
-r--(016)    244611954 2015-02-04 13:50:48 /store/user/capalmer/VHBBHeppyNtuples/V7/TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174946/0000/tree_3.root
-r--(016)    245256313 2015-01-30 22:12:42 /store/user/capalmer/VHBBHeppyNtuples/V7/TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174946/0000/tree_4.root
-r--(016)    316239219 2015-01-30 22:13:19 /store/user/capalmer/VHBBHeppyNtuples/V7/TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174946/0000/tree_5.root
-r--(016)    315531606 2015-01-30 22:12:47 /store/user/capalmer/VHBBHeppyNtuples/V7/TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174946/0000/tree_6.root
-r--(016)    294638231 2015-01-30 22:12:37 /store/user/capalmer/VHBBHeppyNtuples/V7/TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174946/0000/tree_7.root
-r--(016)    292758000 2015-01-30 22:13:18 /store/user/capalmer/VHBBHeppyNtuples/V7/TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174946/0000/tree_8.root
-r--(016)    315802211 2015-01-30 22:12:16 /store/user/capalmer/VHBBHeppyNtuples/V7/TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174946/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:14:49 /store/user/capalmer/VHBBHeppyNtuples/V7/TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:14:50 /store/user/capalmer/VHBBHeppyNtuples/V7/TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174954
dr-x(019)            0 2015-01-30 22:11:00 /store/user/capalmer/VHBBHeppyNtuples/V7/TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174954/0000
-r--(016)    301502888 2015-01-30 22:11:00 /store/user/capalmer/VHBBHeppyNtuples/V7/TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174954/0000/tree_1.root
-r--(016)     34547658 2015-01-30 22:10:56 /store/user/capalmer/VHBBHeppyNtuples/V7/TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174954/0000/tree_2.root
dr-x(019)            1 2015-01-30 16:15:57 /store/user/capalmer/VHBBHeppyNtuples/V7/TT_Tune4C_13TeV-pythia8-tauola/VHBB_HEPPY_V7_TT_Tune4C_13TeV-pythia8-tauola__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:59 /store/user/capalmer/VHBBHeppyNtuples/V7/TT_Tune4C_13TeV-pythia8-tauola/VHBB_HEPPY_V7_TT_Tune4C_13TeV-pythia8-tauola__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/150127_174939
dr-x(019)            0 2015-02-03 09:56:35 /store/user/capalmer/VHBBHeppyNtuples/V7/TT_Tune4C_13TeV-pythia8-tauola/VHBB_HEPPY_V7_TT_Tune4C_13TeV-pythia8-tauola__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/150127_174939/0000
-r--(016)    529317069 2015-01-30 22:12:46 /store/user/capalmer/VHBBHeppyNtuples/V7/TT_Tune4C_13TeV-pythia8-tauola/VHBB_HEPPY_V7_TT_Tune4C_13TeV-pythia8-tauola__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/150127_174939/0000/tree_1.root
-r--(016)    454240966 2015-01-30 22:12:47 /store/user/capalmer/VHBBHeppyNtuples/V7/TT_Tune4C_13TeV-pythia8-tauola/VHBB_HEPPY_V7_TT_Tune4C_13TeV-pythia8-tauola__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/150127_174939/0000/tree_2.root
-r--(016)    509485639 2015-01-30 22:14:14 /store/user/capalmer/VHBBHeppyNtuples/V7/TT_Tune4C_13TeV-pythia8-tauola/VHBB_HEPPY_V7_TT_Tune4C_13TeV-pythia8-tauola__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/150127_174939/0000/tree_3.root
-r--(016)    551752076 2015-01-30 22:12:48 /store/user/capalmer/VHBBHeppyNtuples/V7/TT_Tune4C_13TeV-pythia8-tauola/VHBB_HEPPY_V7_TT_Tune4C_13TeV-pythia8-tauola__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/150127_174939/0000/tree_4.root
-r--(016)    460563400 2015-01-30 22:13:29 /store/user/capalmer/VHBBHeppyNtuples/V7/TT_Tune4C_13TeV-pythia8-tauola/VHBB_HEPPY_V7_TT_Tune4C_13TeV-pythia8-tauola__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/150127_174939/0000/tree_5.root
-r--(016)    449322478 2015-01-30 22:12:01 /store/user/capalmer/VHBBHeppyNtuples/V7/TT_Tune4C_13TeV-pythia8-tauola/VHBB_HEPPY_V7_TT_Tune4C_13TeV-pythia8-tauola__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/150127_174939/0000/tree_6.root
-r--(016)    368161407 2015-01-30 22:13:23 /store/user/capalmer/VHBBHeppyNtuples/V7/TT_Tune4C_13TeV-pythia8-tauola/VHBB_HEPPY_V7_TT_Tune4C_13TeV-pythia8-tauola__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/150127_174939/0000/tree_7.root
dr-x(019)            1 2015-01-30 16:14:39 /store/user/capalmer/VHBBHeppyNtuples/V7/TTZJets_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_TTZJets_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:14:41 /store/user/capalmer/VHBBHeppyNtuples/V7/TTZJets_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_TTZJets_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175107
dr-x(019)            0 2015-01-30 22:11:54 /store/user/capalmer/VHBBHeppyNtuples/V7/TTZJets_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_TTZJets_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175107/0000
-r--(016)    439469078 2015-01-30 22:11:54 /store/user/capalmer/VHBBHeppyNtuples/V7/TTZJets_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_TTZJets_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175107/0000/tree_1.root
-r--(016)      7400699 2015-01-30 22:10:47 /store/user/capalmer/VHBBHeppyNtuples/V7/TTZJets_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_TTZJets_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175107/0000/tree_2.root
dr-x(019)            1 2015-01-30 16:16:03 /store/user/capalmer/VHBBHeppyNtuples/V7/TTWJets_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_TTWJets_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:16:04 /store/user/capalmer/VHBBHeppyNtuples/V7/TTWJets_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_TTWJets_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175100
dr-x(019)            0 2015-01-30 22:14:11 /store/user/capalmer/VHBBHeppyNtuples/V7/TTWJets_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_TTWJets_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175100/0000
-r--(016)    423542366 2015-01-30 22:14:10 /store/user/capalmer/VHBBHeppyNtuples/V7/TTWJets_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_TTWJets_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175100/0000/tree_1.root
dr-x(019)            1 2015-01-30 16:15:30 /store/user/capalmer/VHBBHeppyNtuples/V7/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:32 /store/user/capalmer/VHBBHeppyNtuples/V7/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175004
dr-x(019)            0 2015-01-30 23:57:21 /store/user/capalmer/VHBBHeppyNtuples/V7/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175004/0000
-r--(016)    321956689 2015-01-30 22:11:40 /store/user/capalmer/VHBBHeppyNtuples/V7/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175004/0000/tree_1.root
-r--(016)    247928118 2015-01-30 22:12:00 /store/user/capalmer/VHBBHeppyNtuples/V7/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175004/0000/tree_2.root
-r--(016)    283123886 2015-01-30 22:11:39 /store/user/capalmer/VHBBHeppyNtuples/V7/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175004/0000/tree_3.root
-r--(016)    324399317 2015-01-30 23:57:21 /store/user/capalmer/VHBBHeppyNtuples/V7/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175004/0000/tree_4.root
-r--(016)    130635906 2015-01-30 22:11:40 /store/user/capalmer/VHBBHeppyNtuples/V7/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175004/0000/tree_5.root
dr-x(019)            1 2015-01-30 16:15:57 /store/user/capalmer/VHBBHeppyNtuples/V7/TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:58 /store/user/capalmer/VHBBHeppyNtuples/V7/TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175012
dr-x(019)            0 2015-01-30 22:12:39 /store/user/capalmer/VHBBHeppyNtuples/V7/TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175012/0000
-r--(016)    170028475 2015-01-30 22:12:39 /store/user/capalmer/VHBBHeppyNtuples/V7/TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/VHBB_HEPPY_V7_TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_175012/0000/tree_1.root
dr-x(019)            1 2015-01-30 16:15:12 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:11 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:12 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175428
dr-x(019)            0 2015-01-30 22:11:15 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175428/0000
-r--(016)     18927976 2015-01-30 22:11:11 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175428/0000/tree_1.root
-r--(016)     15924217 2015-01-30 22:11:13 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175428/0000/tree_2.root
-r--(016)     17601532 2015-01-30 22:11:13 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175428/0000/tree_3.root
-r--(016)     19996969 2015-01-30 22:11:13 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175428/0000/tree_4.root
-r--(016)     17819663 2015-01-30 22:11:15 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175428/0000/tree_5.root
dr-x(019)            1 2015-01-30 16:15:13 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175419
dr-x(019)            0 2015-01-30 22:11:13 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175419/0000
-r--(016)     59788590 2015-01-30 22:11:11 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175419/0000/tree_1.root
-r--(016)     60900797 2015-01-30 22:11:09 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175419/0000/tree_2.root
-r--(016)     44132417 2015-01-30 22:11:11 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175419/0000/tree_3.root
-r--(016)     67416376 2015-01-30 22:11:09 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175419/0000/tree_4.root
-r--(016)     34119536 2015-01-30 22:11:13 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175419/0000/tree_5.root
dr-x(019)            1 2015-01-30 16:15:57 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1
dr-x(019)            1 2015-01-30 16:15:58 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175305
dr-x(019)            0 2015-01-30 22:12:32 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175305/0000
-r--(016)     75570999 2015-01-30 22:11:48 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175305/0000/tree_1.root
-r--(016)     93292386 2015-01-30 22:11:56 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175305/0000/tree_2.root
-r--(016)     93692910 2015-01-30 22:12:32 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175305/0000/tree_3.root
-r--(016)     68625332 2015-01-30 22:11:55 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175305/0000/tree_4.root
-r--(016)     73787760 2015-01-30 22:12:09 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175305/0000/tree_5.root
-r--(016)     75145758 2015-01-30 22:12:20 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175305/0000/tree_6.root
-r--(016)     83595091 2015-01-30 22:11:58 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175305/0000/tree_7.root
-r--(016)     84972149 2015-01-30 22:12:25 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175305/0000/tree_8.root
-r--(016)     74978843 2015-01-30 22:11:50 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175305/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:14:48 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:14:48 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:14:49 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175410
dr-x(019)            0 2015-01-30 22:10:58 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175410/0000
-r--(016)        82481 2015-01-30 22:10:55 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175410/0000/tree_1.root
-r--(016)        85463 2015-01-30 22:10:58 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175410/0000/tree_2.root
-r--(016)        81734 2015-01-30 22:10:55 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175410/0000/tree_3.root
-r--(016)        82335 2015-01-30 22:10:58 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175410/0000/tree_4.root
-r--(016)        72735 2015-01-30 22:10:55 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175410/0000/tree_5.root
dr-x(019)            1 2015-01-30 16:14:49 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175403
dr-x(019)            0 2015-01-30 22:11:12 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175403/0000
-r--(016)       433420 2015-01-30 22:10:55 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175403/0000/tree_1.root
-r--(016)       450528 2015-01-30 22:10:56 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175403/0000/tree_2.root
-r--(016)       507577 2015-01-30 22:10:55 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175403/0000/tree_3.root
-r--(016)       481286 2015-01-30 22:10:55 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175403/0000/tree_4.root
-r--(016)       359674 2015-01-30 22:11:12 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-5to10_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175403/0000/tree_5.root
dr-x(019)            1 2015-01-30 16:15:15 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v1
dr-x(019)            1 2015-01-30 16:15:13 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1
dr-x(019)            1 2015-01-30 16:15:15 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v1
dr-x(019)            1 2015-01-30 16:15:17 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v1/150127_175258
dr-x(019)            0 2015-02-04 13:50:24 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v1/150127_175258/0000
-r--(016)     54074064 2015-01-30 22:11:15 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v1/150127_175258/0000/tree_1.root
-r--(016)     61436301 2015-01-30 22:11:18 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v1/150127_175258/0000/tree_2.root
-r--(016)     61143236 2015-01-30 22:11:30 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v1/150127_175258/0000/tree_3.root
-r--(016)     56704267 2015-01-30 22:11:24 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v1/150127_175258/0000/tree_4.root
-r--(016)     55695341 2015-02-04 13:50:24 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v1/150127_175258/0000/tree_5.root
-r--(016)     57980504 2015-01-30 22:11:29 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v1/150127_175258/0000/tree_6.root
-r--(016)     58593851 2015-01-30 22:11:16 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v1/150127_175258/0000/tree_7.root
-r--(016)     55306786 2015-01-30 22:11:24 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v1/150127_175258/0000/tree_8.root
-r--(016)     18745531 2015-01-30 22:11:21 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v1/150127_175258/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:15:15 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175251
dr-x(019)            0 2015-01-30 22:11:27 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175251/0000
-r--(016)     39263029 2015-01-30 22:11:17 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175251/0000/tree_1.root
-r--(016)     39833724 2015-01-30 22:11:14 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175251/0000/tree_2.root
-r--(016)     46724516 2015-01-30 22:11:14 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175251/0000/tree_3.root
-r--(016)     42898529 2015-01-30 22:11:27 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175251/0000/tree_4.root
-r--(016)     42349623 2015-01-30 22:11:24 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175251/0000/tree_5.root
-r--(016)     39261045 2015-01-30 22:11:14 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175251/0000/tree_6.root
-r--(016)     43614920 2015-01-30 22:11:17 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175251/0000/tree_7.root
-r--(016)     44781231 2015-01-30 22:11:19 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175251/0000/tree_8.root
-r--(016)     28252149 2015-01-30 22:11:19 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175251/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:15:18 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v1/150127_175243
dr-x(019)            0 2015-01-30 22:11:24 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v1/150127_175243/0000
-r--(016)     39177964 2015-01-30 22:11:19 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v1/150127_175243/0000/tree_1.root
-r--(016)     37644849 2015-01-30 22:11:12 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v1/150127_175243/0000/tree_2.root
-r--(016)     36546509 2015-01-30 22:11:14 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v1/150127_175243/0000/tree_3.root
-r--(016)     37157435 2015-01-30 22:11:24 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v1/150127_175243/0000/tree_4.root
-r--(016)     35265942 2015-01-30 22:11:21 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v1/150127_175243/0000/tree_5.root
-r--(016)     39264466 2015-01-30 22:11:20 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v1/150127_175243/0000/tree_6.root
-r--(016)     28739164 2015-01-30 22:11:17 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v1/150127_175243/0000/tree_7.root
-r--(016)     37135243 2015-01-30 22:11:15 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v1/150127_175243/0000/tree_8.root
-r--(016)     16413331 2015-01-30 22:11:18 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v1/150127_175243/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:14:52 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1
dr-x(019)            1 2015-01-30 16:14:52 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1
dr-x(019)            1 2015-01-30 16:14:46 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_castor_PHYS14_ST_V1-v1
dr-x(019)            1 2015-01-30 16:14:46 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:14:46 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:14:48 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175356
dr-x(019)            0 2015-01-30 22:12:03 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175356/0000
-r--(016)      1049992 2015-01-30 22:10:55 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175356/0000/tree_1.root
-r--(016)      1025120 2015-01-30 22:10:58 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175356/0000/tree_2.root
-r--(016)       989331 2015-01-30 22:12:03 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175356/0000/tree_3.root
-r--(016)      1074117 2015-01-30 22:10:58 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175356/0000/tree_4.root
-r--(016)       294554 2015-01-30 22:10:55 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175356/0000/tree_5.root
dr-x(019)            1 2015-01-30 16:14:47 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175346
dr-x(019)            0 2015-01-30 22:11:02 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175346/0000
-r--(016)     14427490 2015-01-30 22:11:01 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175346/0000/tree_1.root
-r--(016)     10337796 2015-01-30 22:11:02 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175346/0000/tree_2.root
-r--(016)     12833182 2015-01-30 22:10:59 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175346/0000/tree_3.root
-r--(016)     10053772 2015-01-30 22:10:59 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175346/0000/tree_4.root
-r--(016)     13770809 2015-01-30 22:11:00 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175346/0000/tree_5.root
-r--(016)      3375111 2015-01-30 22:10:58 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175346/0000/tree_6.root
dr-x(019)            1 2015-01-30 16:14:46 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175234
dr-x(019)            0 2015-01-30 22:11:06 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175234/0000
-r--(016)     28757569 2015-01-30 22:10:57 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175234/0000/tree_1.root
-r--(016)     28213597 2015-01-30 22:11:01 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175234/0000/tree_10.root
-r--(016)     14247344 2015-01-30 22:10:59 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175234/0000/tree_11.root
-r--(016)     28651632 2015-01-30 22:11:06 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175234/0000/tree_2.root
-r--(016)     28642387 2015-01-30 22:10:57 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175234/0000/tree_3.root
-r--(016)     28367880 2015-01-30 22:11:02 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175234/0000/tree_4.root
-r--(016)     28718637 2015-01-30 22:10:56 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175234/0000/tree_5.root
-r--(016)     28703221 2015-01-30 22:11:03 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175234/0000/tree_6.root
-r--(016)     26722435 2015-01-30 22:11:04 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175234/0000/tree_7.root
-r--(016)     28702837 2015-01-30 22:11:06 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175234/0000/tree_8.root
-r--(016)     27021334 2015-01-30 22:11:04 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175234/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:14:53 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175226
dr-x(019)            0 2015-01-30 22:11:12 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175226/0000
-r--(016)     19970398 2015-01-30 22:11:01 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175226/0000/tree_1.root
-r--(016)     20237629 2015-01-30 22:11:03 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175226/0000/tree_10.root
-r--(016)     10189320 2015-01-30 22:11:02 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175226/0000/tree_11.root
-r--(016)     20226585 2015-01-30 22:11:11 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175226/0000/tree_2.root
-r--(016)     17869693 2015-01-30 22:11:07 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175226/0000/tree_3.root
-r--(016)     20200265 2015-01-30 22:11:08 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175226/0000/tree_4.root
-r--(016)     20044908 2015-01-30 22:11:03 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175226/0000/tree_5.root
-r--(016)     19842907 2015-01-30 22:11:12 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175226/0000/tree_6.root
-r--(016)     19971560 2015-01-30 22:11:08 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175226/0000/tree_7.root
-r--(016)     18968824 2015-01-30 22:11:08 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175226/0000/tree_8.root
-r--(016)     20072619 2015-01-30 22:11:06 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/150127_175226/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:14:53 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/150127_175217
dr-x(019)            0 2015-01-30 22:11:08 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/150127_175217/0000
-r--(016)     15996909 2015-01-30 22:11:07 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/150127_175217/0000/tree_1.root
-r--(016)     15125107 2015-01-30 22:11:08 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/150127_175217/0000/tree_10.root
-r--(016)     11272133 2015-01-30 22:11:06 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/150127_175217/0000/tree_11.root
-r--(016)     16119692 2015-01-30 22:11:05 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/150127_175217/0000/tree_2.root
-r--(016)     13804803 2015-01-30 22:11:07 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/150127_175217/0000/tree_3.root
-r--(016)     16066998 2015-01-30 22:11:07 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/150127_175217/0000/tree_4.root
-r--(016)     15902576 2015-01-30 22:11:07 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/150127_175217/0000/tree_5.root
-r--(016)     16078947 2015-01-30 22:11:06 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/150127_175217/0000/tree_6.root
-r--(016)     14562264 2015-01-30 22:11:06 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/150127_175217/0000/tree_7.root
-r--(016)     14761433 2015-01-30 22:11:05 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/150127_175217/0000/tree_8.root
-r--(016)     15787856 2015-01-30 22:11:07 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/150127_175217/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:14:39 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v2
dr-x(019)            1 2015-01-30 16:14:39 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1
dr-x(019)            1 2015-01-30 16:14:39 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v2
dr-x(019)            1 2015-01-30 16:14:42 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v2/150127_175208
dr-x(019)            0 2015-01-30 22:11:01 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v2/150127_175208/0000
-r--(016)     16023045 2015-01-30 22:11:01 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v2/150127_175208/0000/tree_1.root
-r--(016)     21813268 2015-01-30 22:10:52 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v2/150127_175208/0000/tree_2.root
-r--(016)     20174149 2015-01-30 22:10:55 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v2/150127_175208/0000/tree_3.root
-r--(016)     20384294 2015-01-30 22:11:00 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v2/150127_175208/0000/tree_4.root
-r--(016)     18613566 2015-01-30 22:10:55 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v2/150127_175208/0000/tree_5.root
-r--(016)     22242009 2015-01-30 22:10:54 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v2/150127_175208/0000/tree_6.root
-r--(016)     20271557 2015-01-30 22:10:51 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v2/150127_175208/0000/tree_7.root
-r--(016)     22124054 2015-01-30 22:10:52 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v2/150127_175208/0000/tree_8.root
-r--(016)      2243936 2015-01-30 22:10:53 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v2/150127_175208/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:14:41 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175159
dr-x(019)            0 2015-01-30 22:26:45 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175159/0000
-r--(016)     13396137 2015-01-30 22:10:51 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175159/0000/tree_1.root
-r--(016)     10709067 2015-01-30 22:26:45 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175159/0000/tree_10.root
-r--(016)      8946325 2015-01-30 22:10:52 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175159/0000/tree_11.root
-r--(016)      4450868 2015-01-30 22:10:54 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175159/0000/tree_2.root
-r--(016)     12042668 2015-01-30 22:10:52 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175159/0000/tree_3.root
-r--(016)      9006352 2015-01-30 22:10:52 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175159/0000/tree_4.root
-r--(016)     11015201 2015-01-30 22:10:52 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175159/0000/tree_5.root
-r--(016)     11056847 2015-01-30 22:10:51 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175159/0000/tree_6.root
-r--(016)     11694547 2015-01-30 22:10:54 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175159/0000/tree_7.root
-r--(016)     11460404 2015-01-30 22:10:51 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175159/0000/tree_8.root
-r--(016)     10371498 2015-01-30 22:10:54 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/150127_175159/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:14:42 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v2/150127_175150
dr-x(019)            0 2015-01-30 22:24:47 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v2/150127_175150/0000
-r--(016)     10893159 2015-01-30 22:10:53 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v2/150127_175150/0000/tree_1.root
-r--(016)      7483366 2015-01-30 22:10:52 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v2/150127_175150/0000/tree_2.root
-r--(016)      9871679 2015-01-30 22:10:52 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v2/150127_175150/0000/tree_3.root
-r--(016)      9695750 2015-01-30 22:10:52 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v2/150127_175150/0000/tree_4.root
-r--(016)      8725207 2015-01-30 22:10:52 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v2/150127_175150/0000/tree_5.root
-r--(016)      9896843 2015-01-30 22:24:47 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v2/150127_175150/0000/tree_6.root
-r--(016)     11834656 2015-01-30 22:10:52 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v2/150127_175150/0000/tree_7.root
-r--(016)     11361535 2015-01-30 22:10:53 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v2/150127_175150/0000/tree_8.root
-r--(016)      9686802 2015-01-30 22:10:58 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8__Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v2/150127_175150/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:15:26 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v2
dr-x(019)            1 2015-01-30 16:15:24 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:26 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175337
dr-x(019)            0 2015-01-30 22:11:29 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175337/0000
-r--(016)       150148 2015-01-30 22:11:23 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175337/0000/tree_1.root
-r--(016)       137929 2015-01-30 22:11:24 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175337/0000/tree_2.root
-r--(016)       160482 2015-01-30 22:11:22 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175337/0000/tree_3.root
-r--(016)       139264 2015-01-30 22:11:29 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175337/0000/tree_4.root
-r--(016)       100556 2015-01-30 22:11:24 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175337/0000/tree_5.root
dr-x(019)            1 2015-01-30 16:15:28 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v2/150127_175329
dr-x(019)            0 2015-01-30 22:11:23 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v2/150127_175329/0000
-r--(016)      3613812 2015-01-30 22:11:23 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v2/150127_175329/0000/tree_1.root
-r--(016)      3567955 2015-01-30 22:11:23 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v2/150127_175329/0000/tree_2.root
-r--(016)      3785152 2015-01-30 22:11:23 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v2/150127_175329/0000/tree_3.root
-r--(016)      2705067 2015-01-30 22:11:23 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v2/150127_175329/0000/tree_4.root
-r--(016)      2080531 2015-01-30 22:11:23 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v2/150127_175329/0000/tree_5.root
dr-x(019)            1 2015-01-30 16:15:08 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:10 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:11 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175321
dr-x(019)            0 2015-01-30 22:11:07 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175321/0000
-r--(016)        82699 2015-01-30 22:11:07 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175321/0000/tree_1.root
-r--(016)        89255 2015-01-30 22:11:06 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175321/0000/tree_2.root
-r--(016)        84402 2015-01-30 22:11:06 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175321/0000/tree_3.root
-r--(016)        87395 2015-01-30 22:11:06 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175321/0000/tree_4.root
-r--(016)        76678 2015-01-30 22:11:06 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/150127_175321/0000/tree_5.root
dr-x(019)            1 2015-01-30 16:15:09 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175312
dr-x(019)            0 2015-01-30 22:11:07 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175312/0000
-r--(016)       508982 2015-01-30 22:11:06 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175312/0000/tree_1.root
-r--(016)       573015 2015-01-30 22:11:07 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175312/0000/tree_2.root
-r--(016)       603038 2015-01-30 22:11:07 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175312/0000/tree_3.root
-r--(016)       511312 2015-01-30 22:11:07 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175312/0000/tree_4.root
-r--(016)       487397 2015-01-30 22:11:07 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8/VHBB_HEPPY_V7_QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/150127_175312/0000/tree_5.root
dr-x(019)            1 2015-01-30 16:15:23 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_250To500_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_250To500_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:24 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_250To500_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_250To500_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v2
dr-x(019)            1 2015-01-30 16:15:26 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_250To500_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_250To500_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v2/150127_174559
dr-x(019)            0 2015-01-30 22:11:30 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_250To500_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_250To500_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v2/150127_174559/0000
-r--(016)     59938468 2015-01-30 22:11:21 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_250To500_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_250To500_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v2/150127_174559/0000/tree_1.root
-r--(016)     86669799 2015-01-30 22:11:30 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_250To500_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_250To500_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v2/150127_174559/0000/tree_2.root
-r--(016)     61977539 2015-01-30 22:11:21 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_250To500_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_250To500_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v2/150127_174559/0000/tree_3.root
-r--(016)     80455939 2015-01-30 22:11:24 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_250To500_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_250To500_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v2/150127_174559/0000/tree_4.root
-r--(016)     54974402 2015-01-30 22:11:25 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_250To500_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_250To500_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v2/150127_174559/0000/tree_5.root
-r--(016)     43499541 2015-01-30 22:11:21 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_250To500_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_250To500_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v2/150127_174559/0000/tree_6.root
dr-x(019)            1 2015-01-30 16:15:25 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_250To500_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_250To500_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174550
dr-x(019)            0 2015-01-30 22:11:33 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_250To500_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_250To500_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174550/0000
-r--(016)     80217896 2015-01-30 22:11:25 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_250To500_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_250To500_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174550/0000/tree_1.root
-r--(016)     48193634 2015-01-30 22:11:33 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_250To500_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_250To500_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174550/0000/tree_2.root
dr-x(019)            1 2015-01-30 16:15:09 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_1000ToInf_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_1000ToInf_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:08 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_1000ToInf_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_1000ToInf_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v1
dr-x(019)            1 2015-01-30 16:15:10 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_1000ToInf_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_1000ToInf_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v1/150127_174630
dr-x(019)            0 2015-01-30 22:11:23 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_1000ToInf_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_1000ToInf_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v1/150127_174630/0000
-r--(016)    313579620 2015-01-30 22:11:15 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_1000ToInf_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_1000ToInf_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v1/150127_174630/0000/tree_1.root
-r--(016)    314406362 2015-01-30 22:11:21 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_1000ToInf_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_1000ToInf_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v1/150127_174630/0000/tree_2.root
-r--(016)    133630284 2015-01-30 22:11:23 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_1000ToInf_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_1000ToInf_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v1/150127_174630/0000/tree_3.root
dr-x(019)            1 2015-01-30 16:15:11 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_1000ToInf_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_1000ToInf_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174638
dr-x(019)            0 2015-01-30 22:11:13 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_1000ToInf_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_1000ToInf_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174638/0000
-r--(016)    153735893 2015-01-30 22:11:13 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_1000ToInf_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_1000ToInf_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174638/0000/tree_1.root
-r--(016)     72623175 2015-01-30 22:11:08 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT_1000ToInf_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT_1000ToInf_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174638/0000/tree_2.root
dr-x(019)            1 2015-01-30 16:15:26 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-500To1000_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-500To1000_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:28 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-500To1000_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-500To1000_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v1
dr-x(019)            1 2015-01-30 16:15:29 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-500To1000_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-500To1000_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v1/150127_174618
dr-x(019)            0 2015-01-30 23:55:54 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-500To1000_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-500To1000_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v1/150127_174618/0000
-r--(016)    190481350 2015-01-30 22:11:42 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-500To1000_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-500To1000_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v1/150127_174618/0000/tree_1.root
-r--(016)    181531472 2015-01-30 22:11:28 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-500To1000_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-500To1000_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v1/150127_174618/0000/tree_2.root
-r--(016)    138566035 2015-01-30 22:11:23 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-500To1000_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-500To1000_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v1/150127_174618/0000/tree_3.root
-r--(016)    175341335 2015-01-30 22:11:27 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-500To1000_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-500To1000_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v1/150127_174618/0000/tree_4.root
-r--(016)    145684472 2015-01-30 22:11:40 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-500To1000_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-500To1000_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v1/150127_174618/0000/tree_5.root
-r--(016)    190195687 2015-01-30 22:11:26 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-500To1000_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-500To1000_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v1/150127_174618/0000/tree_6.root
-r--(016)    198622280 2015-01-30 22:11:29 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-500To1000_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-500To1000_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v1/150127_174618/0000/tree_7.root
-r--(016)     60328392 2015-01-30 23:55:54 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-500To1000_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-500To1000_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v1/150127_174618/0000/tree_8.root
dr-x(019)            1 2015-01-30 16:15:28 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-500To1000_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-500To1000_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174608
dr-x(019)            0 2015-01-30 22:11:26 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-500To1000_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-500To1000_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174608/0000
-r--(016)    191744720 2015-01-30 22:11:26 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-500To1000_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-500To1000_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174608/0000/tree_1.root
-r--(016)    144940231 2015-01-30 22:11:25 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-500To1000_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-500To1000_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174608/0000/tree_2.root
dr-x(019)            1 2015-01-30 16:15:12 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-100To250_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-100To250_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:13 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-100To250_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-100To250_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174543
dr-x(019)            0 2015-02-03 16:17:43 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-100To250_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-100To250_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174543/0000
-r--(016)     15538477 2015-01-30 22:11:11 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-100To250_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-100To250_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174543/0000/tree_1.root
-r--(016)     16922679 2015-01-30 22:11:13 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-100To250_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-100To250_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174543/0000/tree_2.root
-r--(016)     16475170 2015-01-30 22:11:11 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-100To250_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-100To250_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174543/0000/tree_3.root
-r--(016)     16665052 2015-01-30 22:11:16 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-100To250_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-100To250_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174543/0000/tree_4.root
-r--(016)     15466267 2015-01-30 22:41:06 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-100To250_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-100To250_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174543/0000/tree_5.root
-r--(016)     16811082 2015-01-30 22:11:18 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-100To250_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-100To250_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174543/0000/tree_6.root
-r--(016)     16643639 2015-01-30 22:11:13 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-100To250_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-100To250_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174543/0000/tree_7.root
-r--(016)     16977537 2015-01-30 22:11:13 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-100To250_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-100To250_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174543/0000/tree_8.root
-r--(016)      8495875 2015-01-30 22:11:13 /store/user/capalmer/VHBBHeppyNtuples/V7/QCD_HT-100To250_13TeV-madgraph/VHBB_HEPPY_V7_QCD_HT-100To250_13TeV-madgraph__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174543/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:15:05 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:07 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174903
dr-x(019)            0 2015-01-30 22:12:35 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174903/0000
-r--(016)    301307812 2015-01-30 22:12:35 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174903/0000/tree_1.root
-r--(016)    124933480 2015-01-30 22:11:08 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174903/0000/tree_10.root
-r--(016)    301078150 2015-01-30 22:11:51 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174903/0000/tree_2.root
-r--(016)    306091728 2015-01-30 22:11:42 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174903/0000/tree_3.root
-r--(016)    257656911 2015-01-30 22:11:09 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174903/0000/tree_4.root
-r--(016)    259077276 2015-01-30 22:11:31 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174903/0000/tree_5.root
-r--(016)    306884443 2015-01-30 22:11:44 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174903/0000/tree_6.root
-r--(016)    256290335 2015-01-30 22:11:29 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174903/0000/tree_7.root
-r--(016)    288896954 2015-01-30 22:11:12 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174903/0000/tree_8.root
-r--(016)    286700150 2015-01-30 22:11:28 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174903/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:15:50 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:52 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174853
dr-x(019)            0 2015-01-30 22:12:53 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174853/0000
-r--(016)    168616093 2015-01-30 22:12:14 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174853/0000/tree_1.root
-r--(016)    178928004 2015-01-30 22:12:53 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174853/0000/tree_10.root
-r--(016)    163874548 2015-01-30 22:12:31 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174853/0000/tree_11.root
-r--(016)    164312512 2015-01-30 22:12:53 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174853/0000/tree_12.root
-r--(016)    132372681 2015-01-30 22:12:05 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174853/0000/tree_2.root
-r--(016)    200232290 2015-01-30 22:12:10 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174853/0000/tree_3.root
-r--(016)    138622722 2015-01-30 22:12:30 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174853/0000/tree_4.root
-r--(016)    155611479 2015-01-30 22:12:07 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174853/0000/tree_5.root
-r--(016)    199450488 2015-01-30 22:12:19 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174853/0000/tree_6.root
-r--(016)    165898900 2015-01-30 22:11:51 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174853/0000/tree_7.root
-r--(016)    166449093 2015-01-30 22:11:44 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174853/0000/tree_8.root
-r--(016)    161735476 2015-01-30 22:12:37 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174853/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:15:57 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:58 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174841
dr-x(019)            0 2015-01-30 22:12:57 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174841/0000
-r--(016)    125230199 2015-01-30 22:12:53 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174841/0000/tree_1.root
-r--(016)     98520486 2015-01-30 22:12:06 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174841/0000/tree_10.root
-r--(016)    125099366 2015-01-30 22:12:47 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174841/0000/tree_2.root
-r--(016)    115156859 2015-01-30 22:11:59 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174841/0000/tree_3.root
-r--(016)    108227763 2015-01-30 22:11:56 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174841/0000/tree_4.root
-r--(016)    127718787 2015-01-30 22:12:57 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174841/0000/tree_5.root
-r--(016)    103462545 2015-01-30 22:12:49 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174841/0000/tree_6.root
-r--(016)    108359218 2015-01-30 22:11:52 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174841/0000/tree_7.root
-r--(016)    126833770 2015-01-30 22:11:58 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174841/0000/tree_8.root
-r--(016)    123627317 2015-01-30 22:11:53 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174841/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:15:35 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:37 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174832
dr-x(019)            0 2015-01-30 22:11:54 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174832/0000
-r--(016)     44448350 2015-01-30 22:11:41 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174832/0000/tree_1.root
-r--(016)     46419025 2015-01-30 22:11:41 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174832/0000/tree_10.root
-r--(016)     41609536 2015-01-30 22:11:40 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174832/0000/tree_11.root
-r--(016)     54041584 2015-01-30 22:11:36 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174832/0000/tree_2.root
-r--(016)     42299772 2015-01-30 22:11:36 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174832/0000/tree_3.root
-r--(016)     44533602 2015-01-30 22:11:47 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174832/0000/tree_4.root
-r--(016)     48593601 2015-01-30 22:11:36 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174832/0000/tree_5.root
-r--(016)     53565304 2015-01-30 22:11:40 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174832/0000/tree_6.root
-r--(016)     48184109 2015-01-30 22:11:54 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174832/0000/tree_7.root
-r--(016)     45959371 2015-01-30 22:11:36 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174832/0000/tree_8.root
-r--(016)     44504818 2015-01-30 22:11:36 /store/user/capalmer/VHBBHeppyNtuples/V7/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174832/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:15:11 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToMuMu_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToMuMu_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:12 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToMuMu_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToMuMu_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175131
dr-x(019)            0 2015-01-30 22:11:06 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToMuMu_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToMuMu_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175131/0000
-r--(016)      4853840 2015-01-30 22:11:06 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToMuMu_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToMuMu_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175131/0000/tree_1.root
dr-x(019)            1 2015-01-30 16:14:35 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:14:39 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124
dr-x(019)            0 2015-01-30 22:10:51 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124/0000
-r--(016)      4458112 2015-01-30 22:10:46 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124/0000/tree_1.root
-r--(016)     10096127 2015-01-30 22:10:48 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124/0000/tree_10.root
-r--(016)      8459057 2015-01-30 22:10:46 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124/0000/tree_11.root
-r--(016)     13176740 2015-01-30 22:10:49 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124/0000/tree_12.root
-r--(016)      4414704 2015-01-30 22:10:46 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124/0000/tree_13.root
-r--(016)      6362442 2015-01-30 22:10:49 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124/0000/tree_14.root
-r--(016)      2279358 2015-01-30 22:10:46 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124/0000/tree_15.root
-r--(016)      6780445 2015-01-30 22:10:49 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124/0000/tree_16.root
-r--(016)      8237363 2015-01-30 22:10:50 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124/0000/tree_17.root
-r--(016)      2446537 2015-01-30 22:10:46 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124/0000/tree_18.root
-r--(016)     16757496 2015-01-30 22:10:50 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124/0000/tree_19.root
-r--(016)      4368393 2015-01-30 22:10:46 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124/0000/tree_2.root
-r--(016)      2008608 2015-01-30 22:10:46 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124/0000/tree_20.root
-r--(016)      4148376 2015-01-30 22:10:48 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124/0000/tree_21.root
-r--(016)      4277184 2015-01-30 22:10:46 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124/0000/tree_22.root
-r--(016)      1270884 2015-01-30 22:10:46 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124/0000/tree_23.root
-r--(016)      8109510 2015-01-30 22:10:48 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124/0000/tree_3.root
-r--(016)      2882229 2015-01-30 22:10:47 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124/0000/tree_4.root
-r--(016)     11383682 2015-01-30 22:10:51 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124/0000/tree_5.root
-r--(016)      2680851 2015-01-30 22:10:46 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124/0000/tree_6.root
-r--(016)      4224727 2015-01-30 22:10:46 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124/0000/tree_7.root
-r--(016)      5287405 2015-01-30 22:10:46 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124/0000/tree_8.root
-r--(016)      7719877 2015-01-30 22:10:48 /store/user/capalmer/VHBBHeppyNtuples/V7/DYToEE_M-50_Tune4C_13TeV-pythia8/VHBB_HEPPY_V7_DYToEE_M-50_Tune4C_13TeV-pythia8__Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/150127_175124/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:15:32 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToMuMu_PtZ-180_M-50_13TeV-madgraph/VHBB_HEPPY_V7_DYJetsToMuMu_PtZ-180_M-50_13TeV-madgraph__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v3
dr-x(019)            1 2015-01-30 16:15:31 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToMuMu_PtZ-180_M-50_13TeV-madgraph/VHBB_HEPPY_V7_DYJetsToMuMu_PtZ-180_M-50_13TeV-madgraph__Phys14DR-PU40bx25_tsg_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:36 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToMuMu_PtZ-180_M-50_13TeV-madgraph/VHBB_HEPPY_V7_DYJetsToMuMu_PtZ-180_M-50_13TeV-madgraph__Phys14DR-PU40bx25_tsg_PHYS14_25_V1-v1/150127_174821
dr-x(019)            0 2015-01-30 22:11:36 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToMuMu_PtZ-180_M-50_13TeV-madgraph/VHBB_HEPPY_V7_DYJetsToMuMu_PtZ-180_M-50_13TeV-madgraph__Phys14DR-PU40bx25_tsg_PHYS14_25_V1-v1/150127_174821/0000
-r--(016)    160491860 2015-01-30 22:11:36 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToMuMu_PtZ-180_M-50_13TeV-madgraph/VHBB_HEPPY_V7_DYJetsToMuMu_PtZ-180_M-50_13TeV-madgraph__Phys14DR-PU40bx25_tsg_PHYS14_25_V1-v1/150127_174821/0000/tree_1.root
dr-x(019)            1 2015-01-30 16:15:36 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToMuMu_PtZ-180_M-50_13TeV-madgraph/VHBB_HEPPY_V7_DYJetsToMuMu_PtZ-180_M-50_13TeV-madgraph__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v3/150127_174813
dr-x(019)            0 2015-01-30 22:11:22 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToMuMu_PtZ-180_M-50_13TeV-madgraph/VHBB_HEPPY_V7_DYJetsToMuMu_PtZ-180_M-50_13TeV-madgraph__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v3/150127_174813/0000
-r--(016)    150133407 2015-01-30 22:11:22 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToMuMu_PtZ-180_M-50_13TeV-madgraph/VHBB_HEPPY_V7_DYJetsToMuMu_PtZ-180_M-50_13TeV-madgraph__Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v3/150127_174813/0000/tree_1.root
dr-x(019)            1 2015-01-30 16:15:52 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:54 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174806
dr-x(019)            0 2015-01-30 22:14:13 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174806/0000
-r--(016)    644535782 2015-01-30 22:13:50 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174806/0000/tree_1.root
-r--(016)    726119714 2015-01-30 22:13:32 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174806/0000/tree_10.root
-r--(016)    622440931 2015-01-30 22:12:58 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174806/0000/tree_11.root
-r--(016)    508837332 2015-01-30 22:11:55 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174806/0000/tree_2.root
-r--(016)    706429960 2015-01-30 22:13:53 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174806/0000/tree_3.root
-r--(016)    601293301 2015-01-30 22:11:54 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174806/0000/tree_4.root
-r--(016)    631208549 2015-01-30 22:12:21 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174806/0000/tree_5.root
-r--(016)    702610468 2015-01-30 22:12:01 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174806/0000/tree_6.root
-r--(016)    791839409 2015-01-30 22:14:13 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174806/0000/tree_7.root
-r--(016)    614251029 2015-01-30 22:12:47 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174806/0000/tree_8.root
-r--(016)    604021467 2015-01-30 22:13:30 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174806/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:15:26 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:28 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174755
dr-x(019)            0 2015-01-30 22:13:36 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174755/0000
-r--(016)    663903305 2015-01-30 22:13:06 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174755/0000/tree_1.root
-r--(016)    543200454 2015-01-30 22:13:03 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174755/0000/tree_10.root
-r--(016)    567130858 2015-01-30 22:12:34 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174755/0000/tree_11.root
-r--(016)    218217528 2015-01-30 22:12:00 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174755/0000/tree_12.root
-r--(016)    636657468 2015-01-30 22:13:16 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174755/0000/tree_2.root
-r--(016)    599260721 2015-01-30 22:13:36 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174755/0000/tree_3.root
-r--(016)    627002545 2015-01-30 22:13:22 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174755/0000/tree_4.root
-r--(016)    508782565 2015-01-30 22:12:13 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174755/0000/tree_5.root
-r--(016)    607441356 2015-01-30 22:11:48 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174755/0000/tree_6.root
-r--(016)    581139282 2015-01-30 22:11:33 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174755/0000/tree_7.root
-r--(016)    644435425 2015-01-30 22:12:58 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174755/0000/tree_8.root
-r--(016)    651677937 2015-01-30 22:13:25 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174755/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:15:41 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:43 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174748
dr-x(019)            0 2015-01-30 22:14:56 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174748/0000
-r--(016)    467195415 2015-01-30 22:11:51 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174748/0000/tree_1.root
-r--(016)    446400013 2015-01-30 22:11:46 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174748/0000/tree_10.root
-r--(016)    206557855 2015-01-30 22:11:37 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174748/0000/tree_11.root
-r--(016)    444842314 2015-01-30 22:11:44 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174748/0000/tree_2.root
-r--(016)    471563703 2015-01-30 22:11:58 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174748/0000/tree_3.root
-r--(016)    457185606 2015-01-30 22:13:33 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174748/0000/tree_4.root
-r--(016)    489364524 2015-01-30 22:12:27 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174748/0000/tree_5.root
-r--(016)    464224486 2015-01-30 22:14:56 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174748/0000/tree_6.root
-r--(016)    449579421 2015-01-30 22:14:53 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174748/0000/tree_7.root
-r--(016)    359407747 2015-01-30 22:11:37 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174748/0000/tree_8.root
-r--(016)    496435186 2015-01-30 22:12:18 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174748/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:15:32 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:35 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174730
dr-x(019)            0 2015-01-30 23:57:25 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174730/0000
-r--(016)    245708966 2015-01-30 22:12:31 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174730/0000/tree_1.root
-r--(016)    211307648 2015-01-30 22:12:39 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174730/0000/tree_10.root
-r--(016)    171624380 2015-01-30 22:12:17 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174730/0000/tree_2.root
-r--(016)    165940594 2015-01-30 22:12:05 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174730/0000/tree_3.root
-r--(016)    227376029 2015-01-30 22:11:49 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174730/0000/tree_4.root
-r--(016)    217951021 2015-01-30 22:12:57 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174730/0000/tree_5.root
-r--(016)    236341666 2015-01-30 22:12:30 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174730/0000/tree_6.root
-r--(016)    235660530 2015-01-30 23:57:25 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174730/0000/tree_7.root
-r--(016)    238694468 2015-01-30 22:11:27 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174730/0000/tree_8.root
-r--(016)    194966282 2015-01-30 23:57:16 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola/VHBB_HEPPY_V7_DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174730/0000/tree_9.root
dr-x(019)            1 2015-01-30 16:15:54 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_13TeV-madgraph-pythia8/VHBB_HEPPY_V7_DYJetsToLL_M-50_13TeV-madgraph-pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1
dr-x(019)            1 2015-01-30 16:15:56 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_13TeV-madgraph-pythia8/VHBB_HEPPY_V7_DYJetsToLL_M-50_13TeV-madgraph-pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174719
dr-x(019)            0 2015-01-30 22:12:29 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_13TeV-madgraph-pythia8/VHBB_HEPPY_V7_DYJetsToLL_M-50_13TeV-madgraph-pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174719/0000
-r--(016)     32021673 2015-01-30 22:11:52 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_13TeV-madgraph-pythia8/VHBB_HEPPY_V7_DYJetsToLL_M-50_13TeV-madgraph-pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174719/0000/tree_1.root
-r--(016)     32504962 2015-01-30 22:12:24 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_13TeV-madgraph-pythia8/VHBB_HEPPY_V7_DYJetsToLL_M-50_13TeV-madgraph-pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174719/0000/tree_2.root
-r--(016)     28219425 2015-01-30 22:11:57 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_13TeV-madgraph-pythia8/VHBB_HEPPY_V7_DYJetsToLL_M-50_13TeV-madgraph-pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174719/0000/tree_3.root
-r--(016)     32470954 2015-01-30 22:12:29 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_13TeV-madgraph-pythia8/VHBB_HEPPY_V7_DYJetsToLL_M-50_13TeV-madgraph-pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174719/0000/tree_4.root
-r--(016)     29545873 2015-01-30 22:11:51 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_13TeV-madgraph-pythia8/VHBB_HEPPY_V7_DYJetsToLL_M-50_13TeV-madgraph-pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174719/0000/tree_5.root
-r--(016)     29021221 2015-01-30 22:12:08 /store/user/capalmer/VHBBHeppyNtuples/V7/DYJetsToLL_M-50_13TeV-madgraph-pythia8/VHBB_HEPPY_V7_DYJetsToLL_M-50_13TeV-madgraph-pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1/150127_174719/0000/tree_6.root

*/
