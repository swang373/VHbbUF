#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
//#include "TTreeFormula.h"
#include "TString.h"
//#include "TObjString.h"
#include "TPRegexp.h"  // to use TStringToken
#include "TSystem.h"
#include "TROOT.h"

//#include "TMVA/Config.h"
//#include "TMVARegGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

#include "HelperTMVA.h"
#include "HelperFunctions.h"

//#define USE_WH


struct Label
{
    std::string xlabel;
    std::string unit;
    char type;
};

Label MakeLabel(TString str)
{
    if (str.CountChar(';') != 2)
        std::cout << "ERROR: Must contain 2 ';': " << str << std::endl;
    str.ReplaceAll(";;", "; ;");
    TStringToken token(str, ";");
    TString labeltype;
    Label label;

    if(token.NextToken())  label.xlabel = token;
    if(token.NextToken())  label.unit = token;
    if(token.NextToken())  labeltype = token;
    
    if (label.unit == " ")
        label.unit = "";
    if (labeltype != "F" && labeltype != "I")
        std::cout << "ERROR: Must be either 'F' or 'I': " << labeltype << std::endl;
    else
        label.type = (labeltype == "F") ? 'F' : 'I';
    
    return label;
}

////////////////////////////////////////////////////////////////////////////////
/// Main                                                                     ///
////////////////////////////////////////////////////////////////////////////////
void TrainRegression()
{
    gROOT->SetBatch(1);
    gROOT->LoadMacro( "HelperFunctions.h" );  // make functions visible to TTreeFormula

    // Much of the following code comes from the sample script TMVARegression.C. I've discarded excess code for the sake of clarity.
    // swang373, 07.07.2015
   
    // Load the library
    TMVA::Tools::Instance();
   
    // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
    TString outfileName( "TMVAReg.root" );
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

    // Create the factory object. 
    TMVA::Factory* factory = new TMVA::Factory( "TMVARegression", outputFile, 
                                                "!V:!Silent:!Color:!DrawProgressBar:Transformations=I:AnalysisType=Regression" );
 
    const std::vector<std::string> & inputExpressions      = GetInputExpressionsReg();
    const std::vector<std::string> & inputExpressionLabels = GetInputExpressionLabelsReg();
    assert(inputExpressions.size() == inputExpressionLabels.size());

    // Define the input variables that shall be used for the MVA training
    for (UInt_t iexpr=0; iexpr!=inputExpressions.size(); iexpr++){
        Label label = MakeLabel(inputExpressionLabels.at(iexpr));
        TString expr = inputExpressions.at(iexpr);
        factory->AddVariable(expr, label.xlabel, label.unit, label.type);
    }

    // Add the variable carrying the regression target
    factory->AddTarget( "Jet_mcPt[hJCidx]" );
 
    //--------------------------------------------------------------------------
    // Read training and test data
    TFile* input(0);
    TString dirname = "skim_regression/";
    TString prefix  = "skim_";
    TString suffix  = ".root";
    TTree *regTrainTree(0), *regTestTree(0);
    
    std::vector<std::string> processes;
    //processes.push_back("TT");
    processes.push_back("ZnnH125"); // swang373, 30.06.2015
    /*    processes.push_back("ZnnH110");
    processes.push_back("ZnnH115");
    processes.push_back("ZnnH120");
    processes.push_back("ZnnH125");
    processes.push_back("ZnnH130");
    processes.push_back("ZnnH135");
    processes.push_back("ZnnH140");
    processes.push_back("ZnnH145");
    processes.push_back("ZnnH150");
#ifdef USE_WH
    processes.push_back("WlnH110");
    processes.push_back("WlnH115");
    processes.push_back("WlnH120");
    processes.push_back("WlnH125");
    processes.push_back("WlnH130");
    processes.push_back("WlnH135");
    processes.push_back("WlnH140");
    processes.push_back("WlnH145");
    processes.push_back("WlnH150");
#endif
    */

    std::vector<TFile*> files;
    for (UInt_t i=0; i<processes.size(); i++){
        std::string process = processes.at(i);
        input = (TFile*) TFile::Open(dirname + prefix + process + suffix, "READ");
        if (!input) {
            std::cout << "ERROR: Could not open input file." << std::endl;
            exit(1);
        }
        std::cout << "--- TMVARegression           : Using input file: " << input->GetName() << std::endl;
        files.push_back(input);
        
        // --- Register the regression tree
        regTrainTree = (TTree*) input->Get("tree_train");
        regTestTree  = (TTree*) input->Get("tree_test");

        // Global event weights per tree (see below for setting event-wise weights)
        Double_t regWeight = 1.0;

        // You can add an arbitrary number of regression trees
        factory->AddRegressionTree(regTrainTree, regWeight, TMVA::Types::kTraining);
        factory->AddRegressionTree(regTestTree , regWeight, TMVA::Types::kTesting );
    }

    // Apply additional cuts on the signal and background samples (can be different)
    TCut mycut = "";
    //TCut mycut = "hJet_genPt[0] > 0. && hJet_genPt[1] > 0. && hJet_csv[0] > 0. && hJet_csv[1] > 0. && hJet_pt[0] > 20. && hJet_pt[1] > 20. && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5";

    // Tell the factory to use all remaining events in the trees after training for testing:
    factory->PrepareTrainingAndTestTree( mycut, "V:nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents" );

    // --- Book MVA methods
    
    factory->BookMethod( TMVA::Types::kBDT, "BDTG",
    "!H:V:NTrees=2000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.7:nCuts=200:MaxDepth=3:NNodesMax=15" );

    //--------------------------------------------------------------------------
    // Train MVAs using the set of training events
    factory->TrainAllMethods();

    // --- Evaluate all MVAs using the set of test events
    factory->TestAllMethods();

    // --- Evaluate and compare performance of all configured MVAs
    factory->EvaluateAllMethods();

    //--------------------------------------------------------------------------
    // Save the output
    outputFile->Close();

    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVARegression is done!" << std::endl;

    for (UInt_t i=0; i<files.size(); i++)
        files.at(i)->Close();

    delete outputFile;
    delete factory;

    // Launch the GUI for the root macros
    //gROOT->SetMacroPath( "$ROOTSYS/tmva/macros/" );
    //gROOT->Macro( "$ROOTSYS/tmva/macros/TMVAlogon.C" );
    //gROOT->LoadMacro( "$ROOTSYS/tmva/macros/TMVAGui.C" );
    //if (!gROOT->IsBatch()) TMVARegGui( outfileName );
}
