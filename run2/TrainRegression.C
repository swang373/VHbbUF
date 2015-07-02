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
void TrainRegression(TString myMethodList="BDTG")
{
    gROOT->SetBatch(1);
    gROOT->LoadMacro("HelperFunctions.h" );  // make functions visible to TTreeFormula

    /*  
  if (!TString(gROOT->GetVersion()).Contains("5.34")) {
        std::cout << "INCORRECT ROOT VERSION! Please use 5.34:" << std::endl;
        std::cout << "source /uscmst1/prod/sw/cms/slc5_amd64_gcc462/lcg/root/5.34.02-cms/bin/thisroot.csh" << std::endl;
        std::cout << "Return without doing anything." << std::endl;
        return;
    }
    */

    //TString curDynamicPath( gSystem->GetDynamicPath() );
    //gSystem->SetDynamicPath( "../lib:" + curDynamicPath );

    //TString curIncludePath(gSystem->GetIncludePath());
    //gSystem->SetIncludePath( " -I../include " + curIncludePath );

    // Load the library
    TMVA::Tools::Instance();


    //--------------------------------------------------------------------------
    // Default MVA methods to be trained + tested
    std::map<std::string, int> Use;

    // --- Mutidimensional likelihood and Nearest-Neighbour methods
    Use["PDERS"]           = 0;
    Use["PDEFoam"]         = 1;
    Use["KNN"]             = 1;
    //
    // --- Linear Discriminant Analysis
    Use["LD"]              = 1;
    //
    // --- Function Discriminant analysis
    Use["FDA_GA"]          = 1;
    Use["FDA_MC"]          = 0;
    Use["FDA_MT"]          = 0;
    Use["FDA_GAMT"]        = 0;
    //
    // --- Neural Network
    Use["MLP"]             = 1; 
    //
    // --- Support Vector Machine 
    Use["SVM"]             = 0;
    // 
    // --- Boosted Decision Tree
    Use["BDT"]             = 0;
    Use["BDT1"]            = 0;
    Use["BDTG"]            = 1;
    Use["BDTG1"]           = 0;

    //--------------------------------------------------------------------------
    std::cout << std::endl;
    std::cout << "==> Start TMVARegression" << std::endl;

    // Select methods (don't look at this code - not of interest)
    if (myMethodList != "") {
        for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

        std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
        for (UInt_t i=0; i<mlist.size(); i++) {
            std::string regMethod(mlist[i]);

            if (Use.find(regMethod) == Use.end()) {
                std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
                for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
                std::cout << std::endl;
                return;
            }
            Use[regMethod] = 1;
        }
    }

    //--------------------------------------------------------------------------
    // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
    TString outfileName( "TMVAReg.root" );
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

    // Create the factory object. Later you can choose the methods
    // whose performance you'd like to investigate. The factory will
    // then run the performance analysis for you.
    //
    // The first argument is the base of the name of all the
    // weightfiles in the directory weights/
    //
    // The second argument is the output file for the training results
    // All TMVA output can be suppressed by removing the "!" (not) in 
    // front of the "Silent" argument in the option string
    TMVA::Factory *factory = new TMVA::Factory( "TMVARegression", outputFile, 
                                                "!V:!Silent:!Color:!DrawProgressBar:Transformations=I:AnalysisType=Regression" );

    // If you wish to modify default settings
    // (please check "src/Config.h" to see all available global options)
    //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
    //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

    const std::vector<std::string> & inputExpressions      = GetInputExpressionsReg();
    const std::vector<std::string> & inputExpressionLabels = GetInputExpressionLabelsReg();
    assert(inputExpressions.size() == inputExpressionLabels.size());

    // Define the input variables that shall be used for the MVA training
    // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
    // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
    //factory->AddVariable( "var1", "Variable 1", "units", 'F' );
    //factory->AddVariable( "var2", "Variable 2", "units", 'F' );
    
    for (UInt_t iexpr=0; iexpr!=inputExpressions.size(); iexpr++){
        Label label = MakeLabel(inputExpressionLabels.at(iexpr));
        TString expr = inputExpressions.at(iexpr);
        factory->AddVariable(expr, label.xlabel, label.unit, label.type);
    }

    // You can add so-called "Spectator variables", which are not used in the MVA training,
    // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
    // input variables, the response values of all trained MVAs, and the spectator variables
    //factory->AddSpectator( "spec1 := var1*2",  "Spectator 1", "units", 'F' );
    //factory->AddSpectator( "spec2 := var1*3",  "Spectator 2", "units", 'F' );

    // Add the variable carrying the regression target
    //factory->AddTarget( "fvalue" );
    factory->AddTarget( "Jet_mcPt[hJCidx]" );

    // It is also possible to declare additional targets for multi-dimensional regression, ie:
    // -- factory->AddTarget( "fvalue2" );
    // BUT: this is currently ONLY implemented for MLP

    //--------------------------------------------------------------------------
    // Read training and test data
    TFile *input(0);
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

    std::vector<TFile *> files;
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

    // Set individual event weights (the variables must exist in the original TTree)
    //factory->SetWeightExpression( "var1", "Regression" );

    // Apply additional cuts on the signal and background samples (can be different)
    TCut mycut = ""; // for example: TCut mycut = "abs(var1)<0.5 && abs(var2-0.5)<1";
    //TCut mycut = "hJet_genPt[0] > 0. && hJet_genPt[1] > 0. && hJet_csv[0] > 0. && hJet_csv[1] > 0. && hJet_pt[0] > 20. && hJet_pt[1] > 20. && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5";

    // Tell the factory to use all remaining events in the trees after training for testing:
    factory->PrepareTrainingAndTestTree( mycut, "V:nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents" );

    // If no numbers of events are given, half of the events in the tree are used 
    // for training, and the other half for testing:
    //    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=Random:!V" );

    // --- Book MVA methods
    //
    // Please lookup the various method configuration options in the corresponding cxx files, eg:
    // src/MethodCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
    // it is possible to preset ranges in the option string in which the cut optimisation should be done:
    // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

    // PDE - RS method
    if (Use["PDERS"])
        factory->BookMethod( TMVA::Types::kPDERS, "PDERS",
                             "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=40:NEventsMax=60:VarTransform=None" );
    // And the options strings for the MinMax and RMS methods, respectively:
    //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );   
    //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );   

    if (Use["PDEFoam"])
        factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoam",
                             "!H:!V:MultiTargetRegression=F:TargetSelection=Mpv:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Compress=T:Kernel=None:Nmin=10:VarTransform=None" );

    // K-Nearest Neighbour classifier (KNN)
    if (Use["KNN"])
        factory->BookMethod( TMVA::Types::kKNN, "KNN",
                             "nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

    // Linear discriminant
    if (Use["LD"])
        factory->BookMethod( TMVA::Types::kLD, "LD", 
                             "!H:!V:VarTransform=None" );

    // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
    if (Use["FDA_MC"])
        factory->BookMethod( TMVA::Types::kFDA, "FDA_MC",
                             "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100):FitMethod=MC:SampleSize=100000:Sigma=0.1:VarTransform=D" );

    if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options) .. the formula of this example is good for parabolas
        factory->BookMethod( TMVA::Types::kFDA, "FDA_GA",
                             "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100):FitMethod=GA:PopSize=100:Cycles=3:Steps=30:Trim=True:SaveBestGen=1:VarTransform=Norm" );

    if (Use["FDA_MT"])
        factory->BookMethod( TMVA::Types::kFDA, "FDA_MT",
                             "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

    if (Use["FDA_GAMT"])
        factory->BookMethod( TMVA::Types::kFDA, "FDA_GAMT",
                             "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

    // Neural network (MLP)
    if (Use["MLP"])
        factory->BookMethod( TMVA::Types::kMLP, "MLP", 
                             "!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=20000:HiddenLayers=N+20:TestRate=6:TrainingMethod=BFGS:Sampling=0.3:SamplingEpoch=0.8:ConvergenceImprove=1e-6:ConvergenceTests=15:!UseRegulator" );

    // Support Vector Machine
    if (Use["SVM"])
        factory->BookMethod( TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

    // Boosted Decision Trees
    if (Use["BDT"])
        factory->BookMethod( TMVA::Types::kBDT, "BDT",
                             "!H:V:NTrees=100:nEventsMin=30:NodePurityLimit=0.5:BoostType=AdaBoostR2:SeparationType=RegressionVariance:nCuts=20:PruneMethod=CostComplexity:PruneStrength=30" );
//"!H:V:NTrees=60:nEventsMin=20:NodePurityLimit=0.5:BoostType=AdaBoostR2:SeparationType=RegressionVariance:nCuts=20:PruneMethod=CostComplexity:PruneStrength=30:DoBoostMonitor" );

    if (Use["BDT1"])
        factory->BookMethod( TMVA::Types::kBDT, "BDT1",
                             "!H:V:NTrees=100:nEventsMin=5:BoostType=AdaBoostR2:SeparationType=RegressionVariance:nCuts=20:PruneMethod=CostComplexity:PruneStrength=30" );

    if (Use["BDTG"])
        factory->BookMethod( TMVA::Types::kBDT, "BDTG",
                             "!H:V:NTrees=2000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.7:nCuts=200:MaxDepth=3:NNodesMax=15" );

    if (Use["BDTG1"])
        factory->BookMethod( TMVA::Types::kBDT, "BDTG1",
                             "!H:V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=20:MaxDepth=3:NNodesMax=15" );

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
