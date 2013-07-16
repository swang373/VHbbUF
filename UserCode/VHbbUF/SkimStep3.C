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

//#include "HelperNtuples.h"
//#include "HelperFunctions.h"

//#define XROOTD

void SkimStep3(TString process="ZnnH125", Long64_t skimentries=1000000000)
{
    //gROOT->LoadMacro("HelperFunctions.h" );  // make functions visible to TTreeFormula
    gROOT->SetBatch(1);

    TChain * chain  = new TChain("tree");
    TString fname   = "";
    TString dijet   = "DiJetPt_";
#ifndef XROOTD
    TString tag     = "Step3_20130314";
    TString dirMC   = "dcache:/pnfs/cms/WAX/resilient/jiafu/ZnunuHbb/" + tag + "/";
    TString dirData = "dcache:/pnfs/cms/WAX/resilient/jiafu/ZnunuHbb/" + tag + "/";
    //TString dirMC   = "dcache:/pnfs/cms/WAX/11/store/user/lpchbb/apana/" + tagMC + "/";
    //TString dirData = "dcache:/pnfs/cms/WAX/11/store/user/lpchbb/apana/" + tagData + "/";
#else
    TString dirMC   = "root://xrootd.ba.infn.it//store/user/arizzi/" + tagMC + "/";
    TString dirData = "root://xrootd.ba.infn.it//store/user/arizzi/" + tagData + "/";
#endif
    TString outdir  = "skim_ZnnH_Step3/";
    TString prefix  = "Step3_";
    TString suffix  = ".root";

    if (process == "Wj") {
        fname = dirMC + prefix + "WJetsPtW70" + suffix;
        chain->Add(fname);
        fname = dirMC + prefix + "WJetsPtW100" + suffix;
        chain->Add(fname);
        fname = dirMC + prefix + "WJetsPtW180" + suffix;
        chain->Add(fname);
    
    } else if (process == "Zj") {
        fname = dirMC + prefix + "ZJetsHT50" + suffix;
        chain->Add(fname);
        fname = dirMC + prefix + "ZJetsHT100" + suffix;
        chain->Add(fname);
        fname = dirMC + prefix + "ZJetsHT200" + suffix;
        chain->Add(fname);
        fname = dirMC + prefix + "ZJetsHT400" + suffix;
        chain->Add(fname);
        fname = dirMC + prefix + "ZJetsPtZ100" + suffix;
        chain->Add(fname);
    
    } else if (process == "TT") {
        fname = dirMC + prefix + "TTFullLeptMG" + suffix;
        chain->Add(fname);
        fname = dirMC + prefix + "TTSemiLeptMG" + suffix;
        chain->Add(fname);
        fname = dirMC + prefix + "TTHadronicMG" + suffix;
        chain->Add(fname);
    
    } else if (process == "s_Top") {
        fname = dirMC + prefix + "T_t" + suffix;
        chain->Add(fname);
        fname = dirMC + prefix + "Tbar_t" + suffix;
        chain->Add(fname);
        fname = dirMC + prefix + "T_s" + suffix;
        chain->Add(fname);
        fname = dirMC + prefix + "Tbar_s" + suffix;
        chain->Add(fname);
        fname = dirMC + prefix + "T_tW" + suffix;
        chain->Add(fname);
        fname = dirMC + prefix + "Tbar_tW" + suffix;
        chain->Add(fname);
    
    } else if (process == "VV") {
        fname = dirMC + prefix + "WW" + suffix;
        chain->Add(fname);
        fname = dirMC + prefix + "WZ" + suffix;
        chain->Add(fname);
        fname = dirMC + prefix + "ZZ" + suffix;
        chain->Add(fname);
    
    } else if (process == "QCD") {
        //fname = dirMC + prefix + "QCDPt50" + suffix;
        //chain->Add(fname);
        fname = dirMC + prefix + "QCDPt80" + suffix;
        chain->Add(fname);
        fname = dirMC + prefix + "QCDPt120" + suffix;
        chain->Add(fname);
        fname = dirMC + prefix + "QCDPt170" + suffix;
        chain->Add(fname);
        fname = dirMC + prefix + "QCDPt300" + suffix;
        chain->Add(fname);
        fname = dirMC + prefix + "QCDPt470" + suffix;
        chain->Add(fname);
        fname = dirMC + prefix + "QCDPt600" + suffix;
        chain->Add(fname);
        fname = dirMC + prefix + "QCDPt800" + suffix;
        chain->Add(fname);
        fname = dirMC + prefix + "QCDPt1000" + suffix;
        chain->Add(fname);
        fname = dirMC + prefix + "QCDPt1400" + suffix;
        chain->Add(fname);
        fname = dirMC + prefix + "QCDPt1800" + suffix;
        chain->Add(fname);
    
    } else if (process == "data_obs") {
        fname = dirData + prefix + "Data_R" + suffix;
        chain->Add(fname);
        fname = dirData + prefix + "Data_P" + suffix;
        chain->Add(fname);
    
    } else if (process == "trig_METBTag") {
        fname = dirData + prefix + "Data_METBTag_R" + suffix;
        chain->Add(fname);
        fname = dirData + prefix + "Data_METBTag_P" + suffix;
        chain->Add(fname);
    
    } else if (process == "trig_SingleMu") {
        fname = dirData + prefix + "Data_SingleMu_R" + suffix;
        chain->Add(fname);
        fname = dirData + prefix + "Data_SingleMu_P" + suffix;
        chain->Add(fname);
    
    } else if (process == "trig_SingleEl") {
        fname = dirData + prefix + "Data_SingleEl_R" + suffix;
        chain->Add(fname);
        fname = dirData + prefix + "Data_SingleEl_P" + suffix;
        chain->Add(fname);
    
    } else if (process.BeginsWith("Data")) {
        fname = dirData + prefix + process + suffix;
        chain->Add(fname);
    
    } else {
        fname = dirMC + prefix + process + suffix;
        chain->Add(fname);
    }

    TCut selection = "(Vtype==2||Vtype==3) && METtype1corr.et>100 && H.pt>100 && hJet_pt[0]>60 && hJet_pt[1]>30 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.244 && H.mass<250 && mindPhiMETJet_dPhi>0.5 && naJets_Znn>=0 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && H.mass>0";  // modified TT+Wj
    
    // Make output directory if it doesn't exist
    if (gSystem->AccessPathName(outdir))
        gSystem->mkdir(outdir);
    TString outname = outdir + prefix + Form("%s.root", process.Data());

    TFile* f1 = TFile::Open(outname, "RECREATE");
    TTree* t1 = (TTree*) chain->CopyTree(selection);

    t1->Write();
    std::clog << process << ": skimmed from " << chain->GetEntriesFast() << " to " << t1->GetEntriesFast() << " entries." << std::endl;
        
    f1->Close();

    return;
}

// To run:
//root -l -b -q SkimStep3.C+\(\"ZnnH125\"\)
