// @(#)root/tmva $Id: Class_ZH.C,v 1.2 2011/07/06 14:05:30 madfish Exp $
/**********************************************************************************
 * Project   : TMVA - a ROOT-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Root Macro: TMVAClassification                                                 *
 *                                                                                *
 * This macro provides examples for the training and testing of the               *
 * TMVA classifiers.                                                              *
 *                                                                                *
 * As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
 * and linearly correlated input variables.                                       *
 *                                                                                *
 * The methods to be used can be switched on and off by means of booleans, or     *
 * via the prompt command, for example:                                           *
 *                                                                                *
 *    root -l ./TMVAClassification.C\(\"Fisher,Likelihood\"\)                     *
 *                                                                                *
 * (note that the backslashes are mandatory)                                      *
 * If no method given, a default set of classifiers is used.                      *
 *                                                                                *
 * The output file "TMVA.root" can be analysed with the use of dedicated          *
 * macros (simply say: root -l <macro.C>), which can be conveniently              *
 * invoked through a GUI that will appear at the end of the run of this macro.    *
 * Launch the GUI via the command:                                                *
 *                                                                                *
 *    root -l ./TMVAGui.C                                                         *
 *                                                                                *
 **********************************************************************************/
//TString *cuts = new TString("H.mass>102 && H.mass< 126 && nHjj < 2 && Jet1.pt>20 && Jet2.pt>20 && Z.mass>75  && Z.mass<105 && ((Jet1.csv>0.9 && Jet2.csv>0.44) || (Jet2.csv>0.9 && Jet1.csv>0.44)) && Z.pt > 150 && jjdr < 1.5 && hjjzdPhi > 2.75");

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
//gSystem->Load( "TMVA/lib/libTMVA.1" )
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

void Class_ZH( TString myMethodList = "" )
{

double lumi = 0.186;
const double lzh = 7638.0.;
 
const double lwh = 1321.; 

const double lzjmad = 0.74456; 
const double lwj = 0.463; 
const double lwz = 42. ;
const double lww = 48.;
const double lzz = 357.; 

const double ltt = 5.05; 

const double lts = 499.; 
const double lTt = 23.; 
const double ltw= 46.; 

double lz0         = 0.586085; 
double lz10         = 30.858 * (2.03/33.3);  
double lz11         = 1080.613 * (0.78 /33.3); 
double lz13         = 60588.240 * (0.90 / 33.3);
double lz18         = 2879753.581 * ( 2.15 / 33.3);
double lz20         = 10.115 * (6.35/33.3); 
double lz21         = 104.075 * (4.54/33.3);
double lz23         = 4482.161 * (6.35/33.3);
double lz28         = 162135.973 * (8.71 / 33.3);
double lz30         =   (8.50/33.3) *       11.084;   
double lz31         =   (11.0/33.3) *       58.528;   
double lz33         =   (10.2/33.3) *     1729.905;   
double lz38         =   (23.1/33.3) *    77584.126;   
double lz40         =   (31.9/33.3) *       18.665;   
double lz41         =   (17.4/33.3) *       73.680;   
double lz43         =   (9.76/33.3) *      915.589;   
double lz48         =   (29.6/33.3) *    73920.903;   
double lz50         =   (15.0/33.3) *       21.582;   
double lz51         =   (13.1/33.3) *       50.731;   
double lz53         =   (19.4/33.3) *     1006.313;   
double lz58         =   (13.2/33.3) *   346530.969;  
double lz0b= 119.400  ; // 20.68. ;
double lz1b= 123.100 ; // 211000. ;
double lz2b= 17.400 ; // 30000.;
double lz3b= 39.900 ; // 30000.;
double lz0c= 150.100 ;//* 4./ 5. ; //256000. * 4./5.;
double lz1c= 113.200 ;
double lz2c= 17.160 ;
double lz3c= 36.200 ;

const double eelz0 = 0.586085;  //49530.  * 1.3;
const double eelz10 = 30.858 * (95.9/33.3);  //49530.  * 1.3;
const double eelz11 = 1080.613 * (98.4 /33.3); //1375000. ;
const double eelz13 = 60588.240 * (98.2 / 33.3); // 76945000.;
const double eelz18 = 2879753.581 * (95.7 / 33.3);// 3660.6800. ;
const double eelz20 = 10.115 * (87.3/33.3); //12850. ;
const double eelz21 = 104.075 * (90.9/33.3); //132000. ;
const double eelz23 = 4482.161 * (87.2/33.3); //5693000. ;
const double eelz28 = 162135.973 * (82.8 / 33.3); //05291000. ;
const double eelz30  =   (83.2/33.33) *     11.084;
const double eelz31  =   (77.9/33.33) *     58.528;
const double eelz33  =   (79.6/33.33) *   1729.905;
const double eelz38  =   (54.2/33.33) *  77584.126;
const double eelz40  =   (36.3/33.33) *     18.665;
const double eelz41  =   (65.2/33.33) *     73.680;
const double eelz43  =   (80.3/33.33) *    915.589;
const double eelz48  =   (41.3/33.33) *  73920.903;
const double eelz50  =   (69.5/33.33) *     21.582;
const double eelz51  =   (74.4/33.33) *     50.731;
const double eelz53  =   (61.3/33.33) *   1006.313;
const double eelz58  =   (73.4/33.33) * 346530.969;
const double eelz0b= 119.400.  ; // 20.68. ;
const double eelz1b= 123.100. ; // 211000. ;
const double eelz2b= 17.400. ; // 30000.;
const double eelz3b= 39.900. ; // 30000.;
const double eelz0c= 150.100. ;//* 4./ 5. ; //256000. * 4./5.;
const double eelz1c= 113.200. ;
const double eelz2c= 17.160. ;
const double eelz3c= 36.200. ;




   //ONLY OF THESE MAY BE TRUE!

   bool trainZJets = 0;
   bool trainZZ    = 0;
   bool trainTTbar = 0;
   bool trainALL   = 1;



   // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
   // if you use your private .rootrc, or run from a different directory, please copy the
   // corresponding lines from .rootrc

   // methods to be processed can be given as an argument; use format:
   //
   // mylinux~> root -l TMVAClassification.C\(\"myMethod1,myMethod2,myMethod3\"\)
   //
   // if you like to use a method via the plugin mechanism, we recommend using
   //
   // mylinux~> root -l TMVAClassification.C\(\"P_myMethod\"\)
   // (an example is given for using the BDT as plugin (see below),
   // but of course the real application is when you write your own
   // method based)

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

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
   Use["BDTG"]            = 0; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   // 
   // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 1;
   // ---------------------------------------------------------------

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

   // --------------------------------------------------------------------------------------------------

   // --- Here the preparation phase begins

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( "TMVA.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory is 
   // the only TMVA object you have to interact with
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in
   // front of the "Silent" argument in the option string


   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );


   // If you wish to modify default settings
   // (please check "src/Config.h" to see all available global options)
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
//   factory->AddVariable( "myvar1 := var1+var2", 'F' );
//   factory->AddVariable( "myvar2 := var1-var2", "Expression 2", "", 'F' );
//   factory->AddVariable( "var3",                "Variable 3", "units", 'F' );
   factory->AddVariable( "CosRest := abs(Jet1.cosTheta) + abs(Jet2.cosTheta)", 'F' );
   factory->AddVariable( "H.mass"    ,  "H Mass"    , "GeV"	, 'F' );
   factory->AddVariable( "H.pt"      ,  "H PT"      , "GeV" 	, 'F' );
   factory->AddVariable( "deltaPullAngle"      ,  "#Delta#theta_{Pull}"      , "rad" 	, 'F' );
   factory->AddVariable( "jjdr"       ,  "JJ #Delta R"     , ""       	, 'F' );
   factory->AddVariable( "Jet1.csv"      ,  "CSV 1"     , ""		, 'F' );
   factory->AddVariable( "Jet2.csv"      ,  "CSV 2"     , "" 	, 'F' );
   factory->AddVariable( "zmmMass"      ,  "Z #rightarrow #mu#mu Mass"    , "GeV"      , 'F' );
   factory->AddVariable( "zmmPt"        ,  "Z #rightarrow #mu#mu p_{T}"	     , "GeV"      , 'F' );
   factory->AddVariable( "hjjzmmdPhi"   ,  "H-Z #Delta #varphi", "rad"      , 'F' );
   factory->AddVariable( "Jet1.pt"       ,  "Jet 1 PT"     , ""		, 'F' );
//   factory->AddVariable( "Jet1.eta"      ,  "Jet 1 Eta"     , ""		, 'F' );
   factory->AddVariable( "Jet1.phi"      ,  "Jet 1 Phi"     , ""		, 'F' );
   factory->AddVariable( "Jet2.pt"      ,   "Jet 2 Pt"     , ""		, 'F' );
//   factory->AddVariable( "Jet2.eta"      ,  "Jet 2 Eta"     , ""		, 'F' );
   factory->AddVariable( "Jet2.phi"      ,  "Jet 2 Phi"     , ""		, 'F' );
   factory->AddVariable( "deltaEta := abs(Jet2.eta-Jet1.eta)"      ,  "J-J Delta Eta"     , ""		, 'F' );
//   factory->AddVariable( "nHjj"       ,  "Num H"   , ""		, 'F' );
// factory->AddVariable( "nofLeptons" ,  "Num Lep" , ""		, 'F' );

   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables
   //  factory->AddSpectator( "zeePt",  "Z pt", "GeV" , 'F' );
//   factory->AddSpectator( "spec1 := var1*2",  "Spectator 1", "units", 'F' );
 //  factory->AddSpectator( "spec2 := var1*3",  "Spectator 2", "units", 'F' );

   // Read training and test data
   // (it is also possible to use ASCII format as input -> see TMVA Users Guide)
   
//   if (gSystem->AccessPathName( fname ))  // file does not exist in local directory
//      gSystem->Exec("wget http://root.cern.ch/files/tmva_class_example.root");
 TFile *inputBtt  = TFile::Open("A-ZllHTTbar.root");                 
 TFile *inputBzz  = TFile::Open("A-ZllHZZ.root");
 TFile *inputts   = TFile::Open("A-ZllHT-s.root");
 TFile *inputTt   = TFile::Open("A-ZllHT-t.root");
 TFile *inputtw   = TFile::Open("A-ZllHT-tW.root");
 TFile *inputww   = TFile::Open("A-ZllHWW.root");
 TFile *inputwz   = TFile::Open("A-ZllHWZ.root");
 TFile *inputwj   = TFile::Open("A-ZllHWjets.root");
 TFile *inputz0   = TFile::Open("A-ZllHZ0Jets.root");
 TFile *input10   =  TFile::Open("A-ZllHZ1Jets0-100.root");
 TFile *input11   =  TFile::Open("A-ZllHZ1Jets100-300.root");
 TFile *input13   =  TFile::Open("A-ZllHZ1Jets300-800.root");
 TFile *input18   =  TFile::Open("A-ZllHZ1Jets800-1600.root");
 TFile *input20   =  TFile::Open("A-ZllHZ2Jets0-100.root");
 TFile *input21   =  TFile::Open("A-ZllHZ2Jets100-300.root");
 TFile *input23   =  TFile::Open("A-ZllHZ2Jets300-800.root");
 TFile *input28   =  TFile::Open("A-ZllHZ2Jets800-1600.root");
 TFile *input30   =  TFile::Open("A-ZllHZ3Jets0-100.root");
 TFile *input31   =  TFile::Open("A-ZllHZ3Jets100-300.root");
 TFile *input33   =  TFile::Open("A-ZllHZ3Jets300-800.root");
 TFile *input38   =  TFile::Open("A-ZllHZ3Jets800-1600.root");
 TFile *input40   =  TFile::Open("A-ZllHZ4Jets0-100.root");
 TFile *input41   =  TFile::Open("A-ZllHZ4Jets100-300.root");
 TFile *input43   =  TFile::Open("A-ZllHZ4Jets300-800.root");
 TFile *input48   =  TFile::Open("A-ZllHZ4Jets800-1600.root");
 TFile *input50   =  TFile::Open("A-ZllHZ5Jets0-100.root");
 TFile *input51   =  TFile::Open("A-ZllHZ5Jets100-300.root");
 TFile *input53   =  TFile::Open("A-ZllHZ5Jets300-800.root");
 TFile *input58   =  TFile::Open("A-ZllHZ5Jets800-1600.root");
 TFile *inputb0   =  TFile::Open("A-ZllHZbb0Jets.root");
 TFile *inputb1   =  TFile::Open("A-ZllHZbb1Jets.root");
 TFile *inputb2   =  TFile::Open("A-ZllHZbb2Jets.root");
 TFile *inputb3   =  TFile::Open("A-ZllHZbb3Jets.root");
 TFile *inputc0   =  TFile::Open("A-ZllHZcc0Jets.root");
 TFile *inputc1   =  TFile::Open("A-ZllHZcc1Jets.root");
 TFile *inputc2   =  TFile::Open("A-ZllHZcc2Jets.root");
 TFile *inputc3   =  TFile::Open("A-ZllHZcc3Jets.root");

 TFile *inputS = TFile::Open("A-ZllHZH115.root");

//   std::cout << "--- TMVAClassification       : Using input file: " << input->GetName() << std::endl;
   
   // --- Register the training and test trees


 TTree *signal        = (TTree*) inputS->Get("tree"); 
 TTree *backgroundtt  = (TTree*) inputBtt ->Get("tree");   
 TTree *backgroundzz  = (TTree*) inputBzz ->Get("tree"); 
//backgroundzj  = (TTree*) inputBzj ->Get("tree"); 
 TTree *backgroundTt = (TTree*) inputTt->Get("tree");
 TTree *backgroundts  = (TTree*) inputts ->Get("tree");
 TTree *backgroundtw  = (TTree*) inputtw ->Get("tree");
 TTree *backgroundww  = (TTree*) inputww ->Get("tree");
 TTree *backgroundwz  = (TTree*) inputwz ->Get("tree");
 TTree *backgroundwj  = (TTree*) inputwj ->Get("tree");
 TTree *backgroundz0  = (TTree*) inputz0 ->Get("tree");
 TTree *background10  = (TTree*) input10 ->Get("tree");
 TTree *background11  = (TTree*) input11 ->Get("tree");
 TTree *background13  = (TTree*) input13 ->Get("tree");
 TTree *background18  = (TTree*) input18 ->Get("tree");
 TTree *background20  = (TTree*) input20 ->Get("tree");
 TTree *background21  = (TTree*) input21 ->Get("tree");
 TTree *background23  = (TTree*) input23 ->Get("tree");
 TTree *background28  = (TTree*) input28 ->Get("tree");
 TTree *background30  = (TTree*) input30 ->Get("tree");
 TTree *background31  = (TTree*) input31 ->Get("tree");
 TTree *background33  = (TTree*) input33 ->Get("tree");
 TTree *background38  = (TTree*) input38 ->Get("tree");
 TTree *background40  = (TTree*) input40 ->Get("tree");
 TTree *background41  = (TTree*) input41 ->Get("tree");
 TTree *background43  = (TTree*) input43 ->Get("tree");
 TTree *background48  = (TTree*) input48 ->Get("tree");
 TTree *background50  = (TTree*) input50 ->Get("tree");
 TTree *background51  = (TTree*) input51 ->Get("tree");
 TTree *background53  = (TTree*) input53 ->Get("tree");
 TTree *background58  = (TTree*) input58 ->Get("tree");
 TTree *backgroundb0  = (TTree*) inputb0 ->Get("tree");
 TTree *backgroundb1  = (TTree*) inputb1 ->Get("tree");
 TTree *backgroundb2  = (TTree*) inputb2 ->Get("tree");
 TTree *backgroundb3  = (TTree*) inputb3 ->Get("tree");
 TTree *backgroundc0  = (TTree*) inputc0 ->Get("tree");
 TTree *backgroundc1  = (TTree*) inputc1 ->Get("tree");
 TTree *backgroundc2  = (TTree*) inputc2 ->Get("tree");
 TTree *backgroundc3  = (TTree*) inputc3 ->Get("tree");

// global event weights per tree (see below for setting event-wise weights)

  
   //Double_t backgroundWeight3 = lumi/lumiZJets;
   
   /*
   Double_t signalWeight1     = 1;
   Double_t backgroundWeight1 =1;
   Double_t backgroundWeight2 =1;
   Double_t backgroundWeight3 =1;
   */
   
   // You can add an arbitrary number of signal or background trees
  factory->AddSignalTree    ( signal	    , 0.186/7639.0  );


if(trainZJets && !trainTTbar && !trainZZ && !trainALL)  factory->AddBackgroundTree(backgroundzj   , 1  );
if(!trainZJets && !trainTTbar && trainZZ && !trainALL)  factory->AddBackgroundTree(backgroundzz   , 1  );
if(!trainZJets && trainTTbar && !trainZZ && !trainALL)  factory->AddBackgroundTree(backgroundtt   , 1  );

if(!trainZJets && !trainTTbar && !trainZZ && trainALL)
{   
factory->AddBackgroundTree(backgroundts,lumi/(lts/2.0));
factory->AddBackgroundTree(backgroundTt,lumi/(lTt/2.0));
factory->AddBackgroundTree(backgroundtw,lumi/(ltw/2.0));
factory->AddBackgroundTree(backgroundww,lumi/(lww/2.0));
factory->AddBackgroundTree(backgroundwz,lumi/(lwz/2.0));
factory->AddBackgroundTree(backgroundwj,lumi/(lwj/2.0));
factory->AddBackgroundTree(backgroundz0,lumi/(lz0 /2.0));
factory->AddBackgroundTree(background10,lumi/(lz10/2.0));
factory->AddBackgroundTree(background11,lumi/(lz11/2.0));
factory->AddBackgroundTree(background13,lumi/(lz13/2.0));
factory->AddBackgroundTree(background18,lumi/(lz18/2.0));
factory->AddBackgroundTree(background20,lumi/(lz20/2.0));
factory->AddBackgroundTree(background21,lumi/(lz21/2.0));
factory->AddBackgroundTree(background23,lumi/(lz23/2.0));
factory->AddBackgroundTree(background28,lumi/(lz28/2.0));
factory->AddBackgroundTree(background30,lumi/(lz30/2.0));
factory->AddBackgroundTree(background31,lumi/(lz31/2.0));
factory->AddBackgroundTree(background33,lumi/(lz33/2.0));
factory->AddBackgroundTree(background38,lumi/(lz38/2.0));
factory->AddBackgroundTree(background40,lumi/(lz40/2.0));
factory->AddBackgroundTree(background41,lumi/(lz41/2.0));
factory->AddBackgroundTree(background43,lumi/(lz43/2.0));
factory->AddBackgroundTree(background48,lumi/(lz48/2.0));
factory->AddBackgroundTree(background50,lumi/(lz50/2.0));
factory->AddBackgroundTree(background51,lumi/(lz51/2.0));
factory->AddBackgroundTree(background53,lumi/(lz53/2.0));
factory->AddBackgroundTree(background58,lumi/(lz58/2.0));
factory->AddBackgroundTree(backgroundb0,lumi/(lz0b/2.0));
factory->AddBackgroundTree(backgroundb1,lumi/(lz1b/2.0));
factory->AddBackgroundTree(backgroundb2,lumi/(lz2b/2.0));
factory->AddBackgroundTree(backgroundb3,lumi/(lz3b/2.0));
factory->AddBackgroundTree(backgroundc0,lumi/(lz0c/2.0));
factory->AddBackgroundTree(backgroundc1,lumi/(lz1c/2.0));
factory->AddBackgroundTree(backgroundc2,lumi/(lz2c/2.0));
factory->AddBackgroundTree(backgroundc3,lumi/(lz3c/2.0));

factory->AddBackgroundTree(backgroundzz   , lumi/(lzz/2.0));
  factory->AddBackgroundTree(backgroundtt   , lumi/(ltt/2.0));

}

//  factory->SetSignalWeightExpression("weight*");
//  factory->SetBackgroundWeightExpression( "weight*" );



   // To give different trees for training and testing, do as follows:
   //    factory->AddSignalTree( signalTrainingTree, signalTrainWeight, "Training" );
   //    factory->AddSignalTree( signalTestTree,     signalTestWeight,  "Test" );
   
   // Use the following code instead of the above two or four lines to add signal and background
   // training and test events "by hand"
   // NOTE that in this case one should not give expressions (such as "var1+var2") in the input
   //      variable definition, but simply compute the expression before adding the event
   //
   //     // --- begin ----------------------------------------------------------
   //     std::vector<Double_t> vars( 4 ); // vector has size of number of input variables
   //     Float_t  treevars[4], weight;
   //     
   //     // Signal
   //     for (UInt_t ivar=0; ivar<4; ivar++) signal->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
   //     for (UInt_t i=0; i<signal->GetEntries(); i++) {
   //        signal->GetEntry(i);
   //        for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
   //        // add training and test events; here: first half is training, second is testing
   //        // note that the weight can also be event-wise
   //        if (i < signal->GetEntries()/2.0) factory->AddSignalTrainingEvent( vars, signalWeight );
   //        else                              factory->AddSignalTestEvent    ( vars, signalWeight );
   //     }
   //   
   //     // Background (has event weights)
   //     background->SetBranchAddress( "weight", &weight );
   //     for (UInt_t ivar=0; ivar<4; ivar++) background->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
   //     for (UInt_t i=0; i<background->GetEntries(); i++) {
   //        background->GetEntry(i);
   //        for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
   //        // add training and test events; here: first half is training, second is testing
   //        // note that the weight can also be event-wise
   //        if (i < background->GetEntries()/2) factory->AddBackgroundTrainingEvent( vars, backgroundWeight*weight );
   //        else                                factory->AddBackgroundTestEvent    ( vars, backgroundWeight*weight );
   //     }
         // --- end ------------------------------------------------------------
   //
   // --- end of tree registration 

   // Set individual event weights (the variables must exist in the original TTree)
   //    for signal    : factory->SetSignalWeightExpression    ("weight1*weight2");
   //    for background: factory->SetBackgroundWeightExpression("weight1*weight2");
//   factory->SetBackgroundWeightExpression( "weight" );

/*
   factory->AddVariable( "H.mass"    ,  "H Mass"  , "GeV"	, 'F' );
   factory->AddVariable( "nHjj"       ,  "Num H"   , ""		, 'F' );
   //   factory->AddVariable( "H.pt"      ,  "H PT"    , "GeV" 	, 'F' );
   factory->AddVariable( "jjdr"       ,  "JJ DR"   , ""       	, 'F' );
   factory->AddVariable( "Jet1.csv"      ,  "CSV 1"   , ""		, 'F' );
   factory->AddVariable( "Jet2.csv"      ,  "CSV 2"   , "" 	, 'F' );
   factory->AddVariable( "Jet1.pt"       ,  "J1 PT"   , "GeV"	, 'F' );
   factory->AddVariable( "nofLeptons" ,  "Num Lep" , ""		, 'F' );
   factory->AddVariable( "zeeMass"    ,  "Z Mass"  , "GeV"      , 'F' );
   factory->AddVariable( "hjjzeedPhi" ,  "H-Z DPhi", "rad"      , 'F' );
   factory->AddVariable( "zeePt" ,  "Z PT", "GeV"      , 'F' );
*/

   // Apply additional cuts on the signal and background samples (can be different)
   //H.mass<126 && H.mass>102 &&
TCut mycuts,mycutb;
bool origCuts = false;

mycuts = "";//Jet1.pt > 60 && Jet2.pt > 20  && Jet1.csv > 0.68 && Jet2.csv >0.68 && BestCSVPair > 0 && zmmMass < 60";
mycutb = "";//Jet1.pt > 60 && Jet2.pt > 20  && Jet1.csv > 0.68 && Jet2.csv >0.68 && BestCSVPair > 0 && zmmMass < 60";
//mycuts = "Jet1.pt > 30 && Jet2.pt && zmmPt > 150 && H.pt > 150 && ((Jet1.csv > 0.9 && Jet2.csv >0.5) || (Jet2.csv > 0.9 && Jet1.csv >0.5)) && hjjzdPhi > 2.70 && zmmMass < 105 && zmmMass > 75";
//mycutb = "Jet1.pt > 30 && Jet2.pt && zmmPt > 150 && H.pt > 150 && ((Jet1.csv > 0.9 && Jet2.csv >0.5) || (Jet2.csv > 0.9 && Jet1.csv >0.5)) && hjjzdPhi > 2.70 && zmmMass < 105 && zmmMass > 75";


//mycuts = "H.pt > 150 && Jet1.csv > 0.5 && Jet2.csv > 0.5 && Z.pt > 150"; 
//mycutb = "H.pt > 150 && Jet1.csv > 0.5 && Jet2.csv > 0.5 && Z.pt > 150"; 



   // Tell the factory how to use the training and testing events
   //
   // If no numbers of events are given, half of the events in the tree are used 
   // for training, and the other half for testing:
   //    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
   // To also specify the number of testing events, use:
   // factory->PrepareTrainingAndTestTree( mycuts,"NSigTrain=300:NBkgTrain=300:NSigTest=300:NBkgTest=300:SplitMode=Random:!V" );


    //ZJETS
 if(trainZJets && !trainTTbar && !trainZZ && !trainALL ) factory->PrepareTrainingAndTestTree( mycuts, mycutb,"nTrain_Signal=100:nTrain_Background=100:SplitMode=Block:NormMode=None:!V" );

   //ZZ
 if(!trainZJets && !trainTTbar && trainZZ && !trainALL)   factory->PrepareTrainingAndTestTree( mycuts, mycutb,"nTrain_Signal=100:nTrain_Background=100:SplitMode=Block:NormMode=None:!V" );


   //TTbar
 if(!trainZJets && trainTTbar && !trainZZ && !trainALL)  factory->PrepareTrainingAndTestTree( mycuts, mycutb,"nTrain_Signal=5:nTrain_Background=5:SplitMode=Block:NormMode=None:!V" );
 
  //ALL
 if(!trainZJets && !trainTTbar && !trainZZ && trainALL)  factory->PrepareTrainingAndTestTree( mycuts, mycutb,"nTrain_Signal=2000:nTrain_Background=2000:SplitMode=Random:NormMode=None:!V" );

  
   // ---- Book MVA methods
   //
   // Please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   // Cut optimisation
   if (Use["Cuts"])
      factory->BookMethod( TMVA::Types::kCuts, "Cuts",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=20.68:VarProp=FSmart" );

   if (Use["CutsD"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsD",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=20.68:VarProp=FSmart:VarTransform=Decorrelate" );

   if (Use["CutsPCA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsPCA",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=20.68:VarProp=FSmart:VarTransform=PCA" );

   if (Use["CutsGA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsGA",
                           "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

   if (Use["CutsSA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsSA",
                           "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   // Likelihood ("naive Bayes estimator")
   if (Use["Likelihood"])
      factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood",
                           "H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );

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
                           "H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0333:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );

   if (Use["PDEFoamBoost"])
      factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoamBoost",
                           "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );

   // K-Nearest Neighbour classifier (KNN)
   if (Use["KNN"])
      factory->BookMethod( TMVA::Types::kKNN, "KNN",
                           "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

   // H-Matrix (chi2-squared) method
   if (Use["HMatrix"])
      factory->BookMethod( TMVA::Types::kHMatrix, "HMatrix", "!H:!V" );

   // Linear discriminant (same as Fisher discriminant)
   if (Use["LD"])
      factory->BookMethod( TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   // Fisher discriminant (same as LD)
   if (Use["Fisher"])
      factory->BookMethod( TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   // Fisher with Gauss-transformed input variables
   if (Use["FisherG"])
      factory->BookMethod( TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );

   // Composite classifier: ensemble (tree) of boosted Fisher classifiers
   if (Use["BoostedFisher"])
      factory->BookMethod( TMVA::Types::kFisher, "BoostedFisher", 
                           "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2" );

   // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
   if (Use["FDA_MC"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MC",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=10.68:Sigma=0.1" );

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
                           "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=20:NNodesMax=5" );

   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=250:nEventsMin=15:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

   if (Use["BDTB"]) // Bagging
      factory->BookMethod( TMVA::Types::kBDT, "BDTB",
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

   if (Use["BDTD"]) // Decorrelation + Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTD",
                           "!H:!V:NTrees=400:nEventsMin=400:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );

   // RuleFit -- TMVA implementation of Friedman's method
   if (Use["RuleFit"])
      factory->BookMethod( TMVA::Types::kRuleFit, "RuleFit",
                           "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );

   // For an example of the category classifier usage, see: TMVAClassificationCategory

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can optimize the setting (configuration) of the MVAs using the set of training events

   // factory->OptimizeAllMethods("SigEffAt001","Scan");
   // factory->OptimizeAllMethods("ROCIntegral","GA");

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can tell the factory to train, test, and evaluate the MVAs

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVAGui( outfileName );
}
