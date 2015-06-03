#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
//#include "TChainElement.h"
#include "TString.h"
#include "TCut.h"
//#include "TH1.h"
#include "TSystem.h"
#include "TROOT.h"

#include "HelperNtuples.h"
//#include "HelperFunctions.h"

//#define XROOTD

void SkimRegression(TString process="ZnnH125", Long64_t skimentries=1000000000)
{
    //gROOT->LoadMacro("HelperFunctions.h" );  // make functions visible to TTreeFormula
    gROOT->SetBatch(1);

    TChain * chain  = new TChain("tree");
    TString fname   = "";
    TString dijet   = "";


    TString dirMC  = "/afs/cern.ch/work/d/degrutto/public/MiniAOD/ZnnHbb_Phys14_PU20bx25/skimV11/";
    TString prefix  = "skim_";
    TString suffix  = ".root";
    TString outdir  = "skim_regression/";                                                                                                                                                       


    /*#ifndef XROOTD
    //TString dirMC   = "dcache:/pnfs/cms/WAX/resilient/jiafu/ZnunuHbb/" + tagMC + "/";
    //TString dirData = "dcache:/pnfs/cms/WAX/resilient/jiafu/ZnunuHbb/" + tagData + "/";
    TString dirMC   = "dcache:/pnfs/cms/WAX/11/store/user/lpchbb/apana/" + tagMC + "/";
    TString dirData = "dcache:/pnfs/cms/WAX/11/store/user/lpchbb/apana/" + tagData + "/";
#else
    TString dirMC   = "root://xrootd.ba.infn.it//store/user/arizzi/" + tagMC + "/";
    TString dirData = "root://xrootd.ba.infn.it//store/user/arizzi/" + tagData + "/";
#endif
    TString outdir  = "skim_ZnnH_regression/";
    TString prefix  = "skim_";
    TString suffix  = ".root";

    */

    // ZbbHinv
    if (process == "ZbbHinv105") {
        fname = dirMC + dijet + "ZH_ZToBB_HToInv_M-105_8TeV_pythia6" + suffix;
        chain->Add(fname);
    } else if (process == "ZbbHinv115") {
        fname = dirMC + dijet + "ZH_ZToBB_HToInv_M-115_8TeV_pythia6" + suffix;
        chain->Add(fname);
    } else if (process == "ZbbHinv125") {
        fname = dirMC + dijet + "ZH_ZToBB_HToInv_M-125_8TeV_pythia6" + suffix;
        chain->Add(fname);
    } else if (process == "ZbbHinv135") {
        fname = dirMC + dijet + "ZH_ZToBB_HToInv_M-135_8TeV_pythia6" + suffix;
        chain->Add(fname);
    } else if (process == "ZbbHinv145") {
        fname = dirMC + dijet + "ZH_ZToBB_HToInv_M-145_8TeV_pythia6" + suffix;
        chain->Add(fname);
    } else if (process == "ZbbHinv150") {
        fname = dirMC + dijet + "ZH_ZToBB_HToInv_M-150_8TeV_pythia6" + suffix;
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

	//TTBar

    } else if (process == "TT") {
        fname = dirMC + dijet + prefix + "TTPythia8" + suffix;
        chain->Add(fname);
      fname = dirMC + dijet + prefix + "TTMadv2" + suffix;
        chain->Add(fname);


    } else {
        std::cout << "Error: Unrecognized process. I quit!" << std::endl;
        return;
    }

    TCut selection = regression.c_str();
    
    // Make output directory if it doesn't exist
    if (gSystem->AccessPathName(outdir))
        gSystem->mkdir(outdir);
    TString outname = outdir + prefix + Form("%s.root", process.Data());

    TFile* f1 = TFile::Open(outname, "RECREATE");
    TTree* t1 = (TTree*) chain->CopyTree(selection);

    if (t1->GetEntriesFast() > skimentries){
        TTree* t2 = t1->CopyTree("", "", skimentries);
        t2->SetName(TString(t1->GetName())+"_train");
        TTree* t3 = t1->CopyTree("", "", 1000000000, skimentries);
        t3->SetName(TString(t1->GetName())+"_test");
        t2->Write();
        t3->Write();
        assert(t1->GetEntriesFast() == t2->GetEntriesFast() + t3->GetEntriesFast());
        delete t1;
        std::clog << process << ": skimmed from " << chain->GetEntriesFast() << " to " << t2->GetEntriesFast() << " (training) and " << t3->GetEntriesFast() << " (test) entries." << std::endl;
    } else {
        t1->Write();
        std::clog << process << ": skimmed from " << chain->GetEntriesFast() << " to " << t1->GetEntriesFast() << " entries." << std::endl;
    }
    
    f1->Close();

    return;
}

// To run:
//root -l -b -q SkimRegression.C+\(\"ZnnH125\"\)

