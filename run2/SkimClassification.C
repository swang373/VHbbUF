#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "TBranch.h"
#include "TBranchElement.h"
//#include "TChain.h"
//#include "TChainElement.h"
#include "TString.h"
#include "TCut.h"
//#include "TH1.h"
#include "TSystem.h"
#include "TROOT.h"

#include "HelperTMVA.h"
#include "HelperFunctions.h"


void ResetDeleteBranches(TTree* tree) {
    TObjArray* branches = tree->GetListOfBranches();
    Int_t nb = branches->GetEntriesFast();
    for (Int_t i = 0; i < nb; ++i) {
        TBranch* br = (TBranch*) branches->UncheckedAt(i);
        if (br->InheritsFrom(TBranchElement::Class())) {
            ((TBranchElement*) br)->ResetDeleteObject();
        }
    }
}


void SkimClassification(TString process="ZnnH125")
{
    gROOT->LoadMacro("HelperFunctions.h" );  // make functions visible to TTreeFormula
    gROOT->SetBatch(1);

    //TChain * chain  = new TChain("tree");
    //TString fname   = "";
    //TString dijet   = "DiJetPt_";
    //TString dirMC   = "dcache:/pnfs/cms/WAX/resilient/jiafu/ZnunuHbb/" + tagMC + "/";
    //TString dirData = "dcache:/pnfs/cms/WAX/resilient/jiafu/ZnunuHbb/" + tagData + "/";
    //TString indir   = "/afs/cern.ch/work/d/degrutto/public/MiniAOD/ZnnHbb_Phys14_PU20bx25/skimV11/step3/";
    TString indir   = "/afs/cern.ch/work/s/swang373/private/MiniAOD/ZnnHbb_Phys14_PU20bx25/skimV11/step3/"; // swang373, 30.06.2015
    //TString outdir  = "/afs/cern.ch/work/d/degrutto/public/MiniAOD/ZnnHbb_Phys14_PU20bx25/skimV11/step3/skim_ZnnH_classification/";
    TString outdir  = "/afs/cern.ch/work/s/swang373/private/MiniAOD/ZnnHbb_Phys14_PU20bx25/skimV11/step3/skim_ZnnH_classification/"; // swang373, 30.06.2015
    TString prefix  = "Step3_";
    TString suffix  = ".root";

    TFile *input = TFile::Open(indir + prefix + process + suffix);
    if (!input) {
        std::cout << "ERROR: Could not open input file." << std::endl;
        exit(1);
    }
    TTree *tree   = (TTree *) input->Get("tree");
    Long64_t entries = tree->GetEntriesFast();
    
    // Make output directory if it doesn't exist
    if (gSystem->AccessPathName(outdir))
        gSystem->mkdir(outdir);
    TString outname = outdir + prefix + Form("%s.root", process.Data());

    TFile* output = TFile::Open(outname, "RECREATE");

    // Get selections
    const std::vector < std::string > & selExpressions = GetSelExpressions("ZnunuHighPt") ;
    //    const UInt_t nsels = 3;  // ZnunuHighPt, ZnunuLowPt, ZnunuLowCSV
    const UInt_t nsels = 1;  // just one now...  ZnunuHighPt, ZnunuLowPt, ZnunuLowCSV

    //    assert(nsels == selExpressions.size()); <-- fixME
    // even-number events for training, odd-number events for testing
    TCut evenselection = "evt  %2 == 0";
    TCut oddselection  = "evt %2 == 1";

    TTreeFormula *ttf = 0;
    std::vector < TTreeFormula * >::const_iterator formIt, formItEnd;

    // Loop over selections
    std::vector<Long64_t> ventries;
    for (unsigned int i = 0; i < nsels; i++) {
        TString selname = "ZnunuHighPt";
        if (i == 1) {
            selname = "ZnunuMedPt";
	    //	    selExpressions = GetSelExpressions("ZnunuMedPt") ;
        } else if (i == 2) {
            selname = "ZnunuLowPt";
	    // selExpressions = GetSelExpressions("ZnunuLowPt") ;
        }

        TTree *t1 = (TTree*) tree->CloneTree(0);
        TTree *t2 = (TTree*) tree->CloneTree(0);
        
        t1->SetName(TString::Format("%s_%s_train", tree->GetName(), selname.Data()));
        t2->SetName(TString::Format("%s_%s", tree->GetName(), selname.Data()));
        
        // The clones should not delete any shared i/o buffers.
        ResetDeleteBranches(t1);
        ResetDeleteBranches(t2);
        
        ttf = new TTreeFormula(Form("ttfsel%i", i), selExpressions.at(i).c_str(), tree);
        ttf->SetQuickLoad(1);
        TTreeFormula *ttf1 = new TTreeFormula(Form("ttfeven%i", i), evenselection, tree);
        ttf1->SetQuickLoad(1);
        TTreeFormula *ttf2 = new TTreeFormula(Form("ttfodd%i", i), oddselection, tree);
        ttf2->SetQuickLoad(1);
        if (!ttf || !ttf->GetNdim()) {
            std::cerr << "ERROR: Failed to find any TTree variable from the selection: " << selExpressions.at(i) << std::endl;
            return;
        }

        
        /// Loop over events
        Int_t curTree = tree->GetTreeNumber();
        const Long64_t nentries = tree->GetEntries();
        for (Long64_t ievt = 0; ievt < nentries; ievt++) {
            Long64_t entryNumber = tree->GetEntryNumber(ievt);
            if (entryNumber < 0)  break;
            Long64_t localEntry = tree->LoadTree(entryNumber);
            if (localEntry < 0)  break;
            if (tree->GetTreeNumber() != curTree) {
                curTree = tree->GetTreeNumber();
                ttf ->UpdateFormulaLeaves();  // if using TChain
                ttf1->UpdateFormulaLeaves();  // if using TChain
                ttf2->UpdateFormulaLeaves();  // if using TChain
            }
            
            const Int_t ndata = ttf->GetNdata();
            Bool_t keep = kFALSE;
            for(Int_t current = 0; current<ndata && !keep; current++) {
                keep |= (bool(ttf->EvalInstance(current)) != 0);
            }
            if (!keep)  continue;


            bool even  = (bool) ttf1->EvalInstance();
            bool odd   = (bool) ttf2->EvalInstance();
            if (even && odd) {
                std::cerr << "ERROR: An event cannot be even and odd at the same time." << std::cout;
                return;
            }

            tree->GetEntry(entryNumber, 1);  // get all branches
            if (even) {
                t1->Fill();
            } else {
                t2->Fill();
            }
        }  // end loop over events

        t1->Write();
        t2->Write();
        
        ventries.push_back(t1->GetEntriesFast() + t2->GetEntriesFast());
        
        delete ttf;
        delete ttf1;
        delete ttf2;
    }
    
    std::clog << process << ": skimmed from " << entries << " to " << ventries[0] << " (ZnunuHighPt), " << ventries[1] << " (ZnunuLowPt), " << ventries[2] << " (ZnunuLowCSV) " << " entries." << std::endl;

    output->Close();
    input->Close();

    delete output;
    delete input;

    return;
}

// To run:
//root -l -b -q SkimClassification.C+\(\"ZnnH125\"\)
