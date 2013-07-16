#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TCut.h"
#include "TSystem.h"
#include "TROOT.h"


void Split(TString process="WJetsPtW100")
{
    gROOT->SetBatch(1);

    TString dir     = "skim_ZnnH_baseline/";
    TString prefix  = "skim_";
    TString suffix  = ".root";

    TString fname   = dir + prefix + process + suffix;
    TString fnameb  = fname;
    TString fnamec  = fname;
    TString fnamel  = fname;
    
    if (process.Contains("WJets")) {
        fnameb.ReplaceAll("WJets", "WbJets");
        fnamec.ReplaceAll("WJets", "WcJets");
        fnamel.ReplaceAll("WJets", "WlJets");
    } else if (process.Contains("ZJets")) {
        fnameb.ReplaceAll("ZJets", "ZbJets");
        fnamec.ReplaceAll("ZJets", "ZcJets");
        fnamel.ReplaceAll("ZJets", "ZlJets");
    } else {
        std::cout << "Error: Unrecognized process. I quit!" << std::endl;
        return;
    }

    TCut cutb = "eventFlav==5";
    TCut cutc = "eventFlav==4";
    TCut cutl = "eventFlav!=4 && eventFlav!=5";
    TFile* f0 = TFile::Open(fname, "READ");
    TTree* t0 = (TTree*) f0->Get("tree");
    Long64_t entries = t0->GetEntriesFast();
    

    TFile* f1 = TFile::Open(fnameb, "RECREATE");
    TTree* t1 = (TTree*) t0->CopyTree(cutb);
    Long64_t entriesb = t1->GetEntriesFast();
    t1->Write();
    f1->Close();
    
    TFile* f2 = TFile::Open(fnamec, "RECREATE");
    TTree* t2 = (TTree*) t0->CopyTree(cutc);
    Long64_t entriesc = t2->GetEntriesFast();
    t2->Write();
    f2->Close();
    
    TFile* f3 = TFile::Open(fnamel, "RECREATE");
    TTree* t3 = (TTree*) t0->CopyTree(cutl);
    Long64_t entriesl = t3->GetEntriesFast();
    t3->Write();
    f3->Close();
    
    assert(entries == entriesb + entriesc + entriesl);

    return;
}

//root -l -b -q Split.C+\(\"WJetsPtW100\"\)

