#include <cassert>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
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

#include "TMVA/Config.h"
//#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

#include "HelperTMVA.h"
#include "HelperFunctions.h"
//#include "HelperNtuples.h"


struct Label{
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
void TrainClassification(TString myMethodList="", TString analysis="highpt", int massH=125)
{
    //gROOT->SetBatch(1);
    gROOT->LoadMacro("HelperFunctions.h" );  // make functions visible to TTreeFormula

    if (!TString(gROOT->GetVersion()).Contains("5.34")) {
        std::cout << "INCORRECT ROOT VERSION! Please use 5.34:" << std::endl;
        std::cout << "source /uscmst1/prod/sw/cms/slc5_amd64_gcc462/lcg/root/5.34.02-cms/bin/thisroot.csh" << std::endl;
        std::cout << "Return without doing anything." << std::endl;
        return;
    }
    
    //TString curDynamicPath( gSystem->GetDynamicPath() );
    //gSystem->SetDynamicPath( "../lib:" + curDynamicPath );

    //TString curIncludePath(gSystem->GetIncludePath());
    //gSystem->SetIncludePath( " -I../include " + curIncludePath );

    // Load the library
    TMVA::Tools::Instance();
    
    // Welcome the user
    //TMVA::gTools().TMVAWelcomeMessage();

    //--------------------------------------------------------------------------
    // Default MVA methods to be trained + tested
    std::map<std::string, int> Use;

    // --- Cut optimisation
    Use["Cuts"]            = 1;
    Use["CutsD"]           = 1;
    Use["CutsPCA"]         = 0;
    Use["CutsGA"]          = 0;
    Use["CutsSA"]          = 0;
    // 
    // --- 1-dimensional likelihood ("naive Bayes estimator")
    Use["Likelihood"]      = 1;
    Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
    Use["LikelihoodPCA"]   = 1; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
    Use["LikelihoodKDE"]   = 0;
    Use["LikelihoodMIX"]   = 0;
    //
    // --- Mutidimensional likelihood and Nearest-Neighbour methods
    Use["PDERS"]           = 1;
    Use["PDERSD"]          = 0;
    Use["PDERSPCA"]        = 0;
    Use["PDEFoam"]         = 1;
    Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
    Use["KNN"]             = 1; // k-nearest neighbour method
    //
    // --- Linear Discriminant Analysis
    Use["LD"]              = 1; // Linear Discriminant identical to Fisher
    Use["Fisher"]          = 0;
    Use["FisherG"]         = 0;
    Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
    Use["HMatrix"]         = 0;
    //
    // --- Function Discriminant analysis
    Use["FDA_GA"]          = 1; // minimisation of user-defined function using Genetics Algorithm
    Use["FDA_SA"]          = 0;
    Use["FDA_MC"]          = 0;
    Use["FDA_MT"]          = 0;
    Use["FDA_GAMT"]        = 0;
    Use["FDA_MCMT"]        = 0;
    //
    // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
    Use["MLP"]             = 0; // Recommended ANN
    Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
    Use["MLPBNN"]          = 1; // Recommended ANN with BFGS training method and bayesian regulator
    Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
    Use["TMlpANN"]         = 0; // ROOT's own ANN
    //
    // --- Support Vector Machine 
    Use["SVM"]             = 1;
    // 
    // --- Boosted Decision Trees
    Use["BDT"]             = 1; // uses Adaptive Boost
    Use["BDT1"]            = 0; // uses Adaptive Boost
    Use["BDT2"]            = 0; // uses Adaptive Boost
    Use["BDT3"]            = 0; // uses Adaptive Boost
    Use["BDT4"]            = 0; // uses Adaptive Boost
    Use["BDT5"]            = 0; // uses Adaptive Boost
    Use["BDT6"]            = 0; // uses Adaptive Boost
    Use["BDT7"]            = 0; // uses Adaptive Boost
    Use["BDT8"]            = 0; // uses Adaptive Boost
    Use["BDT9"]            = 0; // uses Adaptive Boost
    Use["BDTG"]            = 0; // uses Gradient Boost
    Use["BDTG1"]           = 0; // uses Gradient Boost
    Use["BDTG2"]           = 0; // uses Gradient Boost
    Use["BDTG3"]           = 0; // uses Gradient Boost
    Use["BDTG4"]           = 0; // uses Gradient Boost
    Use["BDTG5"]           = 0; // uses Gradient Boost
    Use["BDTG6"]           = 0; // uses Gradient Boost
    Use["BDTB"]            = 0; // uses Bagging
    Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
    Use["BDTD1"]           = 0; // decorrelation + Adaptive Boost
    Use["BDTP"]            = 0; // pruning + Adaptive Boost
    Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting 
    Use["BDTF1"]           = 0; // allow usage of fisher discriminant for node splitting 
    // 
    // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
    Use["RuleFit"]         = 1;

    //--------------------------------------------------------------------------
    std::cout << std::endl;
    std::cout << "==> Start TMVAClassification" << std::endl;

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
    TString outfileName( "TMVA.root" );
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
    TMVA::Factory *factory = new TMVA::Factory( Form("TMVAClassification_mH%i", massH), outputFile, 
                                                "!V:!Silent:!Color:!DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

    UInt_t analysisbin = 0;
    analysis.ToLower();
    if (analysis == "highpt") {
        analysisbin = 0;
    } else if (analysis == "lowpt") {
        analysisbin = 1;
    } else if (analysis == "lowcsv") {
        analysisbin = 2;
    } else {
        std::cout << "Cannot recognize 'analysis = " << analysis << "'. 'highpt' is used." << std::endl;
    }

    // If you wish to modify default settings
    // (please check "src/Config.h" to see all available global options)
    //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
    //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";
    
    if (analysisbin == 1)
        (TMVA::gConfig().GetIONames()).fWeightFileDir = "weights_lowpt";
    else if (analysisbin == 2)
        (TMVA::gConfig().GetIONames()).fWeightFileDir = "weights_lowcsv";


    const std::vector<std::string> & inputExpressions      = GetInputExpressions();
    const std::vector<std::string> & inputExpressionLabels = GetInputExpressionLabels();
    assert(inputExpressions.size() == inputExpressionLabels.size());
    //const std::vector<std::string> & specExpressions       = GetSpecExpressions();
    //const std::vector<std::string> & specExpressionLabels  = GetSpecExpressionLabels();
    //assert(specExpressions.size() == specExpressionLabels.size());

    // Define the input variables that shall be used for the MVA training
    // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
    // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
    //factory->AddVariable( "myvar1 := var1+var2", 'F' );
    //factory->AddVariable( "myvar2 := var1-var2", "Expression 2", "", 'F' );
    //factory->AddVariable( "var3",                "Variable 3", "units", 'F' );
    //factory->AddVariable( "var4",                "Variable 4", "units", 'F' );
    
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

    //for (UInt_t iexpr=0; iexpr!=specExpressions.size(); iexpr++){
    //    Label label = MakeLabel(specExpressionLabels.at(iexpr));
    //    TString expr = specExpressions.at(iexpr);
    //    factory->AddSpectator(expr, label.xlabel, label.unit, label.type);
    //}

    //--------------------------------------------------------------------------
    // Read training and test data
    TFile *input(0);
    TString dirname = "skim_ZnnH_classification/";
    TString prefix  = "skim_";
    TString suffix  = ".root";
    TTree *trainTree(0), *testTree(0);
    
    std::vector<std::string> sigprocesses;
    if     (massH == 110)  sigprocesses.push_back("ZnnH110");
    else if(massH == 115)  sigprocesses.push_back("ZnnH115");
    else if(massH == 120)  sigprocesses.push_back("ZnnH120");
    else if(massH == 125)  sigprocesses.push_back("ZnnH125");
    else if(massH == 130)  sigprocesses.push_back("ZnnH130");
    else if(massH == 135)  sigprocesses.push_back("ZnnH135");

    std::vector<std::string> bkgprocesses;
    //bkgprocesses.push_back("ZJetsPtZ100");
    bkgprocesses.push_back("ZJetsHT100");
    bkgprocesses.push_back("ZJetsHT200");
    bkgprocesses.push_back("ZJetsHT400");
    //bkgprocesses.push_back("ZJetsHW");
    //bkgprocesses.push_back("WJetsPtW70");  // too few events
    bkgprocesses.push_back("WJetsPtW100");
    bkgprocesses.push_back("WW");
    bkgprocesses.push_back("WZ");
    bkgprocesses.push_back("ZZ");
    bkgprocesses.push_back("TTJets");
    bkgprocesses.push_back("T_s");
    bkgprocesses.push_back("T_t");
    bkgprocesses.push_back("T_tW");
    bkgprocesses.push_back("Tbar_s");
    bkgprocesses.push_back("Tbar_t");
    bkgprocesses.push_back("Tbar_tW");

    TString treename = "tree";
    if (analysisbin == 0) {
        treename  += "_ZnunuHighPt";
    } else if (analysisbin == 1) {
        treename  += "_ZnunuLowPt";
    } else if (analysisbin == 2) {
        treename  += "_ZnunuLowCSV";
    }

    std::vector<TFile *> files;
    for (UInt_t i=0; i<sigprocesses.size(); i++){
        const std::string & process = sigprocesses.at(i);
        input = (TFile*) TFile::Open(dirname + prefix + process + suffix, "READ");
        if (!input) {
            std::cout << "ERROR: Could not open input file." << std::endl;
            exit(1);
        }
        std::cout << "--- TMVAClassification       : Using input file: " << input->GetName() << std::endl;
        files.push_back(input);

        // --- Register the training and test trees
        trainTree       = (TTree*) input->Get(treename  + "_train");
        testTree        = (TTree*) input->Get(treename  + "_test");
        if (trainTree->GetEntries() < 5) {
            std::cout << "ERROR: Less than 5 events in a tree." << std::endl;
            exit(1);
        }
        
        // Global event weights per tree (see below for setting event-wise weights)
        //Double_t signalWeight     = 1.0;
        //Double_t backgroundWeight = 1.0;
        
        // To give different trees for training and testing, do as follows:
        factory->AddSignalTree(trainTree, 2.0, "Training");
        factory->AddSignalTree(testTree , 2.0, "Test"    );
    }
    
    for (UInt_t i=0; i<bkgprocesses.size(); i++){
        const std::string & process = bkgprocesses.at(i);
        input = (TFile*) TFile::Open(dirname + prefix + process + suffix, "READ");
        if (!input) {
            std::cout << "ERROR: Could not open input file." << std::endl;
            exit(1);
        }
        std::cout << "--- TMVAClassification       : Using input file: " << input->GetName() << std::endl;
        files.push_back(input);

        // --- Register the training and test trees
        trainTree       = (TTree*) input->Get(treename  + "_train");
        testTree        = (TTree*) input->Get(treename  + "_test");
        if (trainTree->GetEntries() < 5) {
            std::cout << "ERROR: Less than 5 events in a tree." << std::endl;
            exit(1);
        }
        
        // Global event weights per tree (see below for setting event-wise weights)
        //Double_t signalWeight     = 1.0;
        //Double_t backgroundWeight = 1.0;
        
        // To give different trees for training and testing, do as follows:
        factory->AddBackgroundTree(trainTree, 2.0, "Training");
        factory->AddBackgroundTree(testTree , 2.0, "Test"    );
    }

    // Set individual event weights (the variables must exist in the original TTree)
    //factory->SetSignalWeightExpression    ("weight1*weight2");
    //factory->SetBackgroundWeightExpression("weight1*weight2");

    const std::vector<std::string> weightExpressions = GetWeightExpressions();
    factory->SetSignalWeightExpression    ( weightExpressions.front() );
    factory->SetBackgroundWeightExpression( weightExpressions.front() );

    // Apply additional cuts on the signal and background samples (can be different)
    TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
    TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";


    // Tell the factory to use all remaining events in the trees after training for testing:
    factory->PrepareTrainingAndTestTree( mycuts, mycutb, "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:V" );

    // If no numbers of events are given, half of the events in the tree are used 
    // for training, and the other half for testing:
    //    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=Random:!V" );
    // To also specify the number of testing events, use:
    //    factory->PrepareTrainingAndTestTree( mycut,
    //                                         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );

    // --- Book MVA methods
    //
    // Please lookup the various method configuration options in the corresponding cxx files, eg:
    // src/MethodCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
    // it is possible to preset ranges in the option string in which the cut optimisation should be done:
    // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

    // Cut optimisation
    if (Use["Cuts"])
        factory->BookMethod( TMVA::Types::kCuts, "Cuts",
                             "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

    if (Use["CutsD"])
        factory->BookMethod( TMVA::Types::kCuts, "CutsD",
                             "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

    if (Use["CutsPCA"])
        factory->BookMethod( TMVA::Types::kCuts, "CutsPCA",
                             "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );

    if (Use["CutsGA"])
        factory->BookMethod( TMVA::Types::kCuts, "CutsGA",
                             "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

    if (Use["CutsSA"])
        factory->BookMethod( TMVA::Types::kCuts, "CutsSA",
                             "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

    // Likelihood ("naive Bayes estimator")
    if (Use["Likelihood"])
        factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood",
                             "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );

    // Decorrelated likelihood
    if (Use["LikelihoodD"])
        factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodD",
                             "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" );

    // PCA-transformed likelihood
    if (Use["LikelihoodPCA"])
        factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodPCA",
                             "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" ); 

    // Use a kernel density estimator to approximate the PDFs
    if (Use["LikelihoodKDE"])
        factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodKDE",
                             "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" ); 

    // Use a variable-dependent mix of splines and kernel density estimator
    if (Use["LikelihoodMIX"])
        factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodMIX",
                             "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ); 

    // Test the multi-dimensional probability density estimator
    // here are the options strings for the MinMax and RMS methods, respectively:
    //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
    //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
    if (Use["PDERS"])
        factory->BookMethod( TMVA::Types::kPDERS, "PDERS",
                             "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );

    if (Use["PDERSD"])
        factory->BookMethod( TMVA::Types::kPDERS, "PDERSD",
                             "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );

    if (Use["PDERSPCA"])
        factory->BookMethod( TMVA::Types::kPDERS, "PDERSPCA",
                             "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );

    // Multi-dimensional likelihood estimator using self-adapting phase-space binning
    if (Use["PDEFoam"])
        factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoam",
                             "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );

    if (Use["PDEFoamBoost"])
        factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoamBoost",
                             "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );

    // K-Nearest Neighbour classifier (KNN)
    if (Use["KNN"])
        factory->BookMethod( TMVA::Types::kKNN, "KNN",
                             "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

    // H-Matrix (chi2-squared) method
    if (Use["HMatrix"])
        factory->BookMethod( TMVA::Types::kHMatrix, "HMatrix", "!H:!V:VarTransform=None" );

    // Linear discriminant (same as Fisher discriminant)
    if (Use["LD"])
        factory->BookMethod( TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

    // Fisher discriminant (same as LD)
    if (Use["Fisher"])
        factory->BookMethod( TMVA::Types::kFisher, "Fisher", "!H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

    // Fisher with Gauss-transformed input variables
    if (Use["FisherG"])
        factory->BookMethod( TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );

    // Composite classifier: ensemble (tree) of boosted Fisher classifiers
    if (Use["BoostedFisher"])
        factory->BookMethod( TMVA::Types::kFisher, "BoostedFisher", 
                             "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring" );

    // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
    if (Use["FDA_MC"])
        factory->BookMethod( TMVA::Types::kFDA, "FDA_MC",
                             "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );

    if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
        factory->BookMethod( TMVA::Types::kFDA, "FDA_GA",
                             "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );

    if (Use["FDA_SA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
        factory->BookMethod( TMVA::Types::kFDA, "FDA_SA",
                             "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

    if (Use["FDA_MT"])
        factory->BookMethod( TMVA::Types::kFDA, "FDA_MT",
                             "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

    if (Use["FDA_GAMT"])
        factory->BookMethod( TMVA::Types::kFDA, "FDA_GAMT",
                             "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

    if (Use["FDA_MCMT"])
        factory->BookMethod( TMVA::Types::kFDA, "FDA_MCMT",
                             "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

    // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
    if (Use["MLP"])
        factory->BookMethod( TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

    if (Use["MLPBFGS"])
        factory->BookMethod( TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

    if (Use["MLPBNN"])
        factory->BookMethod( TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

    // CF(Clermont-Ferrand)ANN
    if (Use["CFMlpANN"])
        factory->BookMethod( TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...  

    // Tmlp(Root)ANN
    if (Use["TMlpANN"])
        factory->BookMethod( TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...

    // Support Vector Machine
    if (Use["SVM"])
        factory->BookMethod( TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

    // Boosted Decision Trees
    if (Use["BDTG"]) // Gradient Boost
        factory->BookMethod( TMVA::Types::kBDT, "BDTG",
                             "!H:!V:NTrees=850:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=20:NNodesMax=5" );

    if (Use["BDTG1"]) // Gradient Boost
        factory->BookMethod( TMVA::Types::kBDT, "BDTG1",
                             "!H:!V:NTrees=100:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad=True:GradBaggingFraction=0.5:nCuts=20:MaxDepth=3" );
//NTrees=150:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad=True:GradBaggingFraction=0.5:nCuts=20:NNodesMax=13
 
    if (Use["BDTG2"]) // Gradient Boost
        factory->BookMethod( TMVA::Types::kBDT, "BDTG2",
                             "!H:!V:NTrees=100:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad=True:GradBaggingFraction=0.5:nCuts=20:MaxDepth=3:nEventsMin=100" );
//NTrees=150:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad=True:GradBaggingFraction=0.5:nCuts=20:NNodesMax=11:nEventsMin=150

    if (Use["BDTG3"]) // Gradient Boost
        factory->BookMethod( TMVA::Types::kBDT, "BDTG3",
                             "!H:V:NTrees=100:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad=True:GradBaggingFraction=0.5:nCuts=20:MaxDepth=3:nEventsMin=150:DoBoostMonitor" );
//NTrees=150:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad=True:GradBaggingFraction=0.5:nCuts=20:MaxDepth=3:nEventsMin=150:DoBoostMonitor
//NTrees=150:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad=True:GradBaggingFraction=0.5:nCuts=20:NNodesMax=5:nEventsMin=150:DoBoostMonitor
//NTrees=150:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad=True:GradBaggingFraction=0.5:nCuts=20:NNodesMax=13:nEventsMin=60:DoBoostMonitor

    if (Use["BDTG4"]) // Gradient Boost
        factory->BookMethod( TMVA::Types::kBDT, "BDTG4",
                             "!H:!V:NTrees=100:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad=True:GradBaggingFraction=0.5:nCuts=20:NNodesMax=5:nEventsMin=100" );
//NTrees=150:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad=True:GradBaggingFraction=0.5:nCuts=20:NNodesMax=13:nEventsMin=110

    if (Use["BDTG5"]) // Gradient Boost
        factory->BookMethod( TMVA::Types::kBDT, "BDTG5",
                             "!H:!V:NTrees=150:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad=True:GradBaggingFraction=0.5:nCuts=20:NNodesMax=5:nEventsMin=150" );
//NTrees=150:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad=True:GradBaggingFraction=0.5:nCuts=20:NNodesMax=11:nEventsMin=100

    if (Use["BDTG6"]) // Gradient Boost
        factory->BookMethod( TMVA::Types::kBDT, "BDTG6",
                             "!H:!V:NTrees=150:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad=True:GradBaggingFraction=0.5:nCuts=20:NNodesMax=5:nEventsMin=50" );

    //! HOW TO TUNE
    //  First figure out the number of trees (50 - 300) based on available statistics.
    //  For each NNodesMax=5,7,9, find the nEventsMin (<200) after which the shape doesn't change.
    //  Then find the best nEventsMin (>40).

    if (Use["BDT"])  // Adaptive Boost
        factory->BookMethod( TMVA::Types::kBDT, "BDT",
                             "!H:V:NTrees=150:nEventsMin=100:NNodesMax=15:BoostType=AdaBoost:AdaBoostBeta=0.25:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:DoBoostMonitor" );
//NTrees=850:nEventsMin=150:MaxDepth=3

    if (Use["BDT1"])  // Adaptive Boost
        factory->BookMethod( TMVA::Types::kBDT, "BDT1",
                             "!H:!V:NTrees=100:nEventsMin=50:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
//NTrees=100:nEventsMin=200:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=35:PruneMethod=NoPruning
//NTrees=100:nEventsMin=120:MaxDepth=4:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=35:PruneMethod=NoPruning

    if (Use["BDT2"])  // Adaptive Boost
        factory->BookMethod( TMVA::Types::kBDT, "BDT2",
                             "!H:!V:NTrees=100:nEventsMin=90:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
//NTrees=100:nEventsMin=90:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning
//NTrees=150:nEventsMin=70:NNodesMax=9:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning

    if (Use["BDT3"])  // Adaptive Boost
        factory->BookMethod( TMVA::Types::kBDT, "BDT3",
                             "!H:!V:NTrees=100:nEventsMin=40:NNodesMax=15:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
//NTrees=150:nEventsMin=90:NNodesMax=9 0.982 0.547

    if (Use["BDT4"])  // Adaptive Boost
        factory->BookMethod( TMVA::Types::kBDT, "BDT4",
                             "!H:V:NTrees=100:nEventsMin=100:NNodesMax=15:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:DoBoostMonitor" );
//NTrees=100:nEventsMin=60:NNodesMax=15:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:DoBoostMonitor
//NTrees=150:nEventsMin=100:NNodesMax=9:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:DoBoostMonitor

    if (Use["BDT5"])  // Adaptive Boost
        factory->BookMethod( TMVA::Types::kBDT, "BDT5",
                             "!H:!V:NTrees=100:nEventsMin=80:NNodesMax=15:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
//NTrees=150:nEventsMin=70:NNodesMax=7:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning

    if (Use["BDT6"])  // Adaptive Boost
        factory->BookMethod( TMVA::Types::kBDT, "BDT6",
                             "!H:!V:NTrees=100:nEventsMin=60:NNodesMax=15:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:DoBoostMonitor" );
//NTrees=150:nEventsMin=120:NNodesMax=7:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:DoBoostMonitor

    if (Use["BDT7"])  // Adaptive Boost
        factory->BookMethod( TMVA::Types::kBDT, "BDT7",
                             "!H:!V:NTrees=100:nEventsMin=120:NNodesMax=15:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
//NTrees=150:nEventsMin=100:NNodesMax=5:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning

    if (Use["BDT8"])  // Adaptive Boost
        factory->BookMethod( TMVA::Types::kBDT, "BDT8",
                             "!H:!V:NTrees=100:nEventsMin=140:NNodesMax=15:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
//NTrees=150:nEventsMin=120:NNodesMax=11:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning

    if (Use["BDT9"])  // Adaptive Boost
        factory->BookMethod( TMVA::Types::kBDT, "BDT9",
                             "!H:!V:NTrees=100:nEventsMin=150:NNodesMax=15:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
//NTrees=100:nEventsMin=180:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning

    if (Use["BDTB"]) // Bagging
        factory->BookMethod( TMVA::Types::kBDT, "BDTB",
                             "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

    if (Use["BDTD"]) // Decorrelation + Adaptive Boost
        factory->BookMethod( TMVA::Types::kBDT, "BDTD",
                             "!H:!V:NTrees=400:nEventsMin=400:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );

    if (Use["BDTD1"]) // Decorrelation + Adaptive Boost
        factory->BookMethod( TMVA::Types::kBDT, "BDTD1",
                             "!H:!V:NTrees=100:nEventsMin=100:NNodesMax=5:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );

    if (Use["BDTP"])  // Prune + Adaptive Boost
        factory->BookMethod( TMVA::Types::kBDT, "BDTP",
                             "!H:!V:NTrees=850:nEventsMin=50:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:PruneStrength=1.5" );

    if (Use["BDTF"])  // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
        factory->BookMethod( TMVA::Types::kBDT, "BDTMitFisher",
                             "!H:!V:NTrees=50:nEventsMin=150:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

    if (Use["BDTF1"])  // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
        factory->BookMethod( TMVA::Types::kBDT, "BDTMitFisher1",
                             "!H:!V:NTrees=150:nEventsMin=100:NNodesMax=9:UseFisherCuts:MinLinCorrForFisher=0.8:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

    // RuleFit -- TMVA implementation of Friedman's method
    if (Use["RuleFit"])
        factory->BookMethod( TMVA::Types::kRuleFit, "RuleFit",
                             "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );

    // For an example of the category classifier usage, see: TMVAClassificationCategory


    //--------------------------------------------------------------------------
    // --- Now you can optimize the setting (configuration) of the MVAs using the set of training events
    // --- See $ROOTSYS/tmva/src/MethodBDT.cxx
    // fomType = ["Separation","ROCIntegral","SigEffAtBkgEff01","SigEffAtBkgEff001","BkgRejAtSigEff05","BkgEffAtSigEff05"]
    // fitType = ["Scan","FitGA","Minuit"]
    // tuneParameters for BDT = ["NTrees","MaxDepth","NodeMinEvents","AdaBoostBeta"]
    //factory->OptimizeAllMethods("SigEffAtBkgEff001","Scan");
    //factory->OptimizeAllMethods("ROCIntegral","GA");

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
    std::cout << "==> TMVAClassification is done!" << std::endl;

    for (UInt_t i=0; i<files.size(); i++)
        files.at(i)->Close();

    delete outputFile;
    delete factory;

    // Launch the GUI for the root macros
    //gROOT->SetMacroPath( "$ROOTSYS/tmva/macros/" );
    //gROOT->Macro( "$ROOTSYS/tmva/macros/TMVAlogon.C" );
    //gROOT->LoadMacro( "$ROOTSYS/tmva/macros/TMVAGui.C" );
    //if (!gROOT->IsBatch()) TMVAGui( outfileName );
}
