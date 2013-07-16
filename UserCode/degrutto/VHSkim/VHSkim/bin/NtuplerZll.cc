#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "DataFormats/FWLite/interface/Run.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

float ScaleCSV(float CSV)
{
if(CSV < 0.68) return 1.0;
if(CSV < 0.9) return  0.96;
else  return  0.94;

}



float ScaleIsoHLT(float pt1, float eta1)
{

float ptMin,ptMax,etaMin,etaMax,scale,error;
float s1 = 0;

TFile *scaleFile = new TFile ("ScaleFactor_muonEffsIsoToHLT_efficiency.root","read");
TTree *tscale = (TTree*) scaleFile->Get("tree");
int count = 0;
tscale->SetBranchAddress("ptMin",&ptMin);
tscale->SetBranchAddress("ptMax",&ptMax);
tscale->SetBranchAddress("etaMin",&etaMin);
tscale->SetBranchAddress("etaMax",&etaMax);
tscale->SetBranchAddress("scale",&scale);
tscale->SetBranchAddress("error",&error);

for(int jentry = 0; jentry < tscale->GetEntries(); jentry++)
  {
   tscale->GetEntry(jentry);
   if((pt1 > ptMin) && (pt1 < ptMax) && (eta1 > etaMin) && (eta1 < etaMax))
    {
    s1 = scale;
    count++;   
    }
  }

if(count == 0 || s1 == 0) 
{
 scaleFile->Close();
 return 1;
}


scaleFile->Close();
return (s1);
}


float ScaleID(float pt1, float eta1)
{

float ptMin,ptMax,etaMin,etaMax,scale,error;
float s1 = 0;

TFile *scaleFile = new TFile ("ScaleFactor_muonEffsRecoToIso_VBTF.root","read");
TTree *tscale = (TTree*) scaleFile->Get("tree");
int count = 0;
tscale->SetBranchAddress("ptMin",&ptMin);
tscale->SetBranchAddress("ptMax",&ptMax);
tscale->SetBranchAddress("etaMin",&etaMin);
tscale->SetBranchAddress("etaMax",&etaMax);
tscale->SetBranchAddress("scale",&scale);
tscale->SetBranchAddress("error",&error);

for(int jentry = 0; jentry < tscale->GetEntries(); jentry++)
  {

   tscale->GetEntry(jentry);
   if((pt1 > ptMin) && (pt1 < ptMax) && (eta1 > etaMin) && (eta1 < etaMax))
    {
    s1 = scale;
    count++;   
    }
   }

if(count == 0 || s1 == 0) 
{
 scaleFile->Close();
 return 1;
}

scaleFile->Close();
return (s1);

}


int main(int argc, char* argv[]) 
{
TTree *_outTree;
 float nHjj,jjdr,jjdPhi,jcvsmva,jbprob1,jbprob2,jTCHP1,jTCHP2,pfmet,jcvs,jjdistSVtx,zmmMass,zmmPt,hjjzmmdPhi,zeeMass,zeePt,hjjzeedPhi,wmnMt,wmnPt,hjjwmndPhi,wenMt,wenPt,hjjwendPhi,csvmva1,csvmva2,hjjzdPhi,hjjwdPhi,zmmEta,zeeEta,wmnPhi,wmnEta,wenEta,wenPhi,zeePhi,zmmPhi,deltaPullAngle,deltaPullAngleAK7,alpgendrcc,alpgendrbb,weight,eeweight,ThirdJetPt, minDeltaPhijetMET,SumJetsPt,cweightTrig,cweightID,cweightCSV; 
int nofLeptons,channel,numJets,BestCSVPair,bFlag,cFlag;

typedef struct 
{
  float mass;  //MT in case of W
  float pt;
  float eta;
  float phi;
} _TrackInfo;

typedef struct 
{
  float mass;  
  float pt;
  float eta;
  float phi;
  float weightID; 
  float weightTrig;  
} _MuInfo;

typedef struct 
{
  int sample;  
  int count;  
} _sampleInfo;

typedef struct 
{
  float et;  
  float sig;
  float phi;
} _METInfo;


struct 
{
  int run;
  int lumi;
  int event;
} _EventInfo;

typedef struct 
{
  float pt;
  float eta;
  float phi;
  float csv;
  float csvweight;
  float cosTheta;
  int numTracksSV;
} _JetInfo;

_METInfo MET;
_JetInfo Jet1,Jet2;
_TrackInfo H,Z,W,e1,e2;
_MuInfo mu1,mu2;
_sampleInfo sampleInfo;

  // define what muon you are using; this is necessary as FWLite is not 
  // capable of reading edm::Views
  //  using pat::Muon;

  // ----------------------------------------------------------------------
  // First Part: 
  //
  //  * enable the AutoLibraryLoader 
  //  * book the histograms of interest 
  //  * open the input file
  // ----------------------------------------------------------------------

  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  // parse arguments
  if ( argc < 2 ) {
    return 0;
  }

  // get the python configuration
  PythonProcessDesc builder(argv[1]);
  const edm::ParameterSet& in  = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput" );
  const edm::ParameterSet& out = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteOutput");
  const edm::ParameterSet& ana = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("Analyzer");

  // now get each parameter
  int maxEvents_( in.getParameter<int>("maxEvents") );
  unsigned int outputEvery_( in.getParameter<unsigned int>("outputEvery") );
  std::string outputFile_( out.getParameter<std::string>("fileName" ) );

std::string inputFile_( in.getParameter<std::string> ("fileName") );

  std::string sampleName(ana.getParameter<std::string>("sampleName"));
  float lumi = ana.getParameter<double>("lumi");
  bool isZinv_( ana.getParameter<bool>("isZinv"));  
  bool isWln_( ana.getParameter<bool>("isWln") );  
  bool isZll_ ( ana.getParameter<bool>("isZll") );


if(isZll_) channel = 1;
if(isWln_) channel = 2;
if(isZinv_) channel = 3;

sampleInfo.sample = 0;
weight = 1;
eeweight = 1;
///////////////////////////////////////////////////////////////////////
//From Michele's MC Studies:                                         // 
//See https://twiki.cern.ch/twiki/bin/viewauth/CMS/MCBackgroundList  // 
///////////////////////////////////////////////////////////////////////
//dataset                                          Tot_evts          Zee       |      Zmumu          Ztautau    sigLO(pb)   sigNLO(pb)  lumiNLO[\fb]
//Z0Jets_Tune//Z2_7TeV-alpgen-tauola                  1435909     478666 (33.3%) |478753 (33.3%)  478490 (33.3%)  1.929e+03   2450         0.586085    
//Z1Jets_pt//Z-0to100_Tune//Z2_7TeV-alpgen-tauola      14904390   14299454 (95.9%) |302244 (2.03%)  302692 (2.03%)  3.808e+02    483           30.858    
//Z1Jets_pt//Z-100to300_Tune//Z2_7TeV-alpgen-tauola    11994815   11806764 (98.4%) | 93929 (0.78%)   94122 (0.78%)  8.721e+00   11.1         1080.613    
//Z1Jets_pt//Z-300to800_Tune//Z2_7TeV-alpgen-tauola     5683177    5581364 (98.2%) | 50636 (0.89%)   51177 (0.90%)  7.386e-02  9.38e-02     60588.240    
//Z1Jets_pt//Z-800to1600_Tune//Z2_7TeV-alpgen-tauola     502517     480803 (95.7%) | 10806 (2.15%)   10908 (2.17%)  1.374e-04  1.745e-04  2879753.581    
//Z2Jets_pt//Z-0to100_Tune//Z2_7TeV-alpgen-tauola       1335287    1165914 (87.3%) | 84724 (6.35%)   84649 (6.34%)  1.035e+02    132           10.115    
//Z2Jets_pt//Z-100to300_Tune//Z2_7TeV-alpgen-tauola     1124019    1021498 (90.9%) | 51054 (4.54%)   51467 (4.58%)  8.534e+00   10.8          104.075    
//Z2Jets_pt//Z-300to800_Tune//Z2_7TeV-alpgen-tauola      655292     571713 (87.2%) | 41610 (6.35%)   41969 (6.40%)  1.151e-01  0.1462        4482.161    
//Z2Jets_pt//Z-800to1600_Tune//Z2_7TeV-alpgen-tauola      62244      51509 (82.8%) |  5422 (8.71%)    5313 (8.54%)  3.023e-04  3.839e-04   162135.973    
//Z3Jets_pt//Z-0to100_Tune//Z2_7TeV-alpgen-tauola        322554     268322 (83.2%) | 27416 (8.50%)   26816 (8.31%)  2.289e+01   29.1           11.084    
//Z3Jets_pt//Z-100to300_Tune//Z2_7TeV-alpgen-tauola      293815     228961 (77.9%) | 32446 (11.0%)   32408 (11.0%)  3.951e+00   5.02           58.528    
//Z3Jets_pt//Z-300to800_Tune//Z2_7TeV-alpgen-tauola      183370     145902 (79.6%) | 18658 (10.2%)   18810 (10.3%)  8.344e-02  0.106         1729.905    
//Z3Jets_pt//Z-800to1600_Tune//Z2_7TeV-alpgen-tauola      24439      13250 (54.2%) |  5657 (23.1%)    5532 (22.6%)  2.480e-04  3.15e-04     77584.126    
//Z4Jets_pt//Z-0to100_Tune//Z2_7TeV-alpgen-tauola        109008      39539 (36.3%) | 34794 (31.9%)   34675 (31.8%)  4.619e+00   5.84           18.665    
//Z4Jets_pt//Z-100to300_Tune//Z2_7TeV-alpgen-tauola      121573      79250 (65.2%) | 21178 (17.4%)   21145 (17.4%)  1.298e+00   1.65           73.680    
//Z4Jets_pt//Z-300to800_Tune//Z2_7TeV-alpgen-tauola       45752      36754 (80.3%) |  4467 (9.76%)    4531 (9.90%)  3.935e-02  4.997e-02      915.589    
//Z4Jets_pt//Z-800to1600_Tune//Z2_7TeV-alpgen-tauola      13084       5402 (41.3%) |  3875 (29.6%)    3807 (29.1%)  1.394e-04  1.77e-04     73920.903    
//Z5Jets_pt//Z-0to100_Tune//Z2_7TeV-alpgen-tauola         31100      21605 (69.5%) |  4666 (15.0%)    4829 (15.5%)  1.135e+00  1.441           21.582    
//Z5Jets_pt//Z-100to300_Tune//Z2_7TeV-alpgen-tauola       30657      22803 (74.4%) |  4012 (13.1%)    3842 (12.5%)  4.758e-01  0.6043          50.731    
//Z5Jets_pt//Z-300to800_Tune//Z2_7TeV-alpgen-tauola       24866      15254 (61.3%) |  4831 (19.4%)    4781 (19.2%)  1.946e-02  2.471e-02     1006.313    
//Z5Jets_pt//Z-800to1600_Tune//Z2_7TeV-alpgen-tauola      31666      23232 (73.4%) |  4169 (13.2%)    4265 (13.5%)  7.195e-05  9.138e-05   346530.969   


const float lumiZ1jets0_100 = 30.858 * (2.03/33.3);  //49530.  * 1.3;
const float lumiZ1jets100_300 = 1080.613 * (0.78 /33.3); //1375000. ;
const float lumiZ1jets300_800 = 60588.240 * (0.90 / 33.3); // 76945000.;
const float lumiZ1jets800_1600 = 2879753.581 * ( 2.15 / 33.3);// 3660000000. ;

const float lumiZ2jets0_100 = 10.115 * (6.35/33.3); //12850. ;
const float lumiZ2jets100_300 = 104.075 * (4.54/33.3); //132000. ;
const float lumiZ2jets300_800 = 4482.161 * (6.35/33.3); //5693000. ;
const float lumiZ2jets800_1600 = 162135.973 * (8.71 / 33.3); //05291000. ;

const float lumiZ3jets0_100     =   (8.50/33.3) *       11.084;   
const float lumiZ3jets100_300   =   (11.0/33.3) *       58.528;   
const float lumiZ3jets300_800   =   (10.2/33.3) *     1729.905;   
const float lumiZ3jets800_1600  =   (23.1/33.3) *    77584.126;   
const float lumiZ4jets0_100     =   (31.9/33.3) *       18.665;   
const float lumiZ4jets100_300   =   (17.4/33.3) *       73.680;   
const float lumiZ4jets300_800   =   (9.76/33.3) *      915.589;   
const float lumiZ4jets800_1600  =   (29.6/33.3) *    73920.903;   
const float lumiZ5jets0_100     =   (15.0/33.3) *       21.582;   
const float lumiZ5jets100_300   =   (13.1/33.3) *       50.731;   
const float lumiZ5jets300_800   =   (19.4/33.3) *     1006.313;   
const float lumiZ5jets800_1600  =   (13.2/33.3) *   346530.969;  



const float lumiZbb0jets = 119.400  ; // 204000. ;
const float lumiZbb1jets = 123.100 ; // 211000. ;
const float lumiZbb2jets = 17.400 ; // 30000.;
const float lumiZcc0jets = 150.100 ;//* 4./ 5. ; //256000. * 4./5.;
const float lumiZcc1jets = 113.200 ;
const float lumiZcc2jets = 17.160 ;
const float lumiZcc3jets = 36.200 ;
const float lumiZbb3jets = 39.900 ; // 30000.;

const float eelumiZ1jets0_100 = 30.858 * (95.9/33.3);  //49530.  * 1.3;
const float eelumiZ1jets100_300 = 1080.613 * (98.4 /33.3); //1375000. ;
const float eelumiZ1jets300_800 = 60588.240 * (98.2 / 33.3); // 76945000.;
const float eelumiZ1jets800_1600 = 2879753.581 * (95.7 / 33.3);// 3660000000. ;
const float eelumiZ2jets0_100 = 10.115 * (87.3/33.3); //12850. ;
const float eelumiZ2jets100_300 = 104.075 * (90.9/33.3); //132000. ;
const float eelumiZ2jets300_800 = 4482.161 * (87.2/33.3); //5693000. ;
const float eelumiZ2jets800_1600 = 162135.973 * (82.8 / 33.3); //05291000. ;
const float eelumiZbb0jets = 119.400  ; // 204000. ;
const float eelumiZbb1jets = 123.100 ; // 211000. ;
const float eelumiZbb2jets = 17.400 ; // 30000.;
const float eelumiZcc0jets = 150.100 ;//* 4./ 5. ; //256000. * 4./5.;
const float eelumiZcc1jets = 113.200 ;
const float eelumiZcc2jets = 17.160 ;
const float eelumiZcc3jets = 36.200 ;
const float eelumiZbb3jets = 39.900 ; // 30000.;

const float eelumiZ3jets0_100     =   (83.2/33.33) *     11.084;
const float eelumiZ3jets100_300   =   (77.9/33.33) *     58.528;
const float eelumiZ3jets300_800   =   (79.6/33.33) *   1729.905;
const float eelumiZ3jets800_1600  =   (54.2/33.33) *  77584.126;
const float eelumiZ4jets0_100     =   (36.3/33.33) *     18.665;
const float eelumiZ4jets100_300   =   (65.2/33.33) *     73.680;
const float eelumiZ4jets300_800   =   (80.3/33.33) *    915.589;
const float eelumiZ4jets800_1600  =   (41.3/33.33) *  73920.903;
const float eelumiZ5jets0_100     =   (69.5/33.33) *     21.582;
const float eelumiZ5jets100_300   =   (74.4/33.33) *     50.731;
const float eelumiZ5jets300_800   =   (61.3/33.33) *   1006.313;
const float eelumiZ5jets800_1600  =   (73.4/33.33) * 346530.969;




const float lumiZH = 870.0;
const float lumiZH2 = 7638.9;
 
const float lumiWH = 1321.; 

const float lumiZJ = 0.74456; 
const float lumiWJ = 0.463; 
const float lumiZ0Jets = 0.586;

const float lumiWZ = 42. ;
const float lumiWW = 48.;
const float lumiZZ =357.; 

const float lumiTT = 5.05; 

const float lumiTs = 499.; 
const float lumiTt = 23.; 
int counters[40] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

 
if(sampleName == "ZllHbb")            {sampleInfo.sample = 1; weight = lumi/lumiZH; eeweight = lumi/lumiZH;}          
if(sampleName == "ZH115")             {sampleInfo.sample = 1; weight = lumi/lumiZH2; eeweight = lumi/lumiZH2;}          
if(sampleName == "ZJetMad")           {sampleInfo.sample = 2; weight = lumi/lumiZJ; eeweight = lumi/lumiZJ;}
if(sampleName == "ZZ")                {sampleInfo.sample = 3; weight = lumi/lumiZZ; eeweight = lumi/lumiZZ;}
if(sampleName == "TTbar")             {sampleInfo.sample = 4; weight = lumi/lumiTT; eeweight = lumi/lumiTT;}

TFile *_outFile	= new TFile(outputFile_.c_str(), "recreate");	
_outTree = new TTree("tree", "myTree");

_outTree->Branch("eventInfo",      &_EventInfo	    	  ,  "run/I:lumi/I:event/I");
_outTree->Branch("sampleInfo"    , &sampleInfo      	  ,  "sample/I:count/I");
_outTree->Branch("H"		,  &H	            	  ,  "mass/F:pt/F:eta:phi/F");
_outTree->Branch("Jet1"		,  &Jet1	    	  ,  "pt/F:eta/F:phi/F:csv/F:csvweight/F:cosTheta/F:numTracksSV/I");
_outTree->Branch("Jet2"		,  &Jet2	    	  ,  "pt/F:eta/F:phi/F:csv/F:csvweight/F:cosTheta/F:numTracksSV/I");
_outTree->Branch("nHjj" 	,  &nHjj            	  ,  "nHjj/F"         );         	
_outTree->Branch("jjdr" 	,  &jjdr            	  ,  "jjdr/F"         );         	
_outTree->Branch("jjdPhi"  	,  &jjdPhi          	  ,  "jjdPhi/F"       );            	
_outTree->Branch("jcvsmva"  	,  &jcvsmva         	  ,  "jcvsmva/F"      );             	
//_outTree->Branch("jjdistSVtx"   ,  &jjdistSVtx      ,  "jjdistSVtx/F"    );                
_outTree->Branch("numJets"      ,  &numJets         	  ,  "numJets/I"       );                
_outTree->Branch("nofLeptons"   ,  &nofLeptons      	  ,  "nofLeptons/I"    );                
_outTree->Branch("channel"      ,  &channel         	  ,  "channel/I"       );                
_outTree->Branch("deltaPullAngle", &deltaPullAngle  	  ,  "deltaPullAngle/F");
_outTree->Branch("alpgendrcc"    , &alpgendrcc      	  ,  "alpgendrcc/F");
_outTree->Branch("alpgendrbb"    , &alpgendrbb      	  ,  "alpgendrbb/F");
_outTree->Branch("bFlag"        , &bFlag          	  ,  "bFlag/I");
_outTree->Branch("cFlag"        , &cFlag          	  ,  "cFlag/I");
_outTree->Branch("weight"        , &weight          	  ,  "weight/F");
_outTree->Branch("eeweight"        , &eeweight            ,  "eeweight/F");
_outTree->Branch("ThirdJetPt", &ThirdJetPt  	          ,  "ThirdJetPt/F");
_outTree->Branch("deltaPullAngleAK7", &deltaPullAngleAK7  ,  "deltaPullAngleAK7/F");
_outTree->Branch("BestCSVPair",       &BestCSVPair  	  ,  "BestCSVPair/I");
_outTree->Branch("cweightID",       &cweightID 		  ,  "cweightID/F");
_outTree->Branch("cweightTrig",       &cweightTrig 	  ,  "cweightTrig/F");
_outTree->Branch("cweightCSV",       &cweightCSV 	  ,  "cweightCSV/F");

//_outTree->Branch("ThirdJetPt",        &ThirdJetPt  ,  "ThirdJetPt/F");

if(isZll_)
{
_outTree->Branch("hjjzdPhi"     ,  &hjjzdPhi   ,   "hjjzdPhi/F" );                
_outTree->Branch("SumJetsPt",       &SumJetsPt  ,  "SumJetsPt/F");
_outTree->Branch("zeeMass"  	,  &zeeMass      ,   "zeeMass/F"    );             	
_outTree->Branch("zeeEta"  	,  &zeeEta      ,   "zeeEta/F"    );             	
_outTree->Branch("zeePt"  	,  &zeePt        ,   "zeePt/F"      );           	
_outTree->Branch("hjjzeedPhi"  ,   &hjjzeedPhi   ,   "hjjzeedPhi/F" );                
_outTree->Branch("zmmMass"  	,  &zmmMass      ,   "zmmMass/F"    );                                     	
_outTree->Branch("zmmPt"  	,  &zmmPt        ,   "zmmPt/F"      );           	
_outTree->Branch("zmmEta"  	,  &zmmEta       ,   "zmmEta/F"      );           	
_outTree->Branch("hjjzmmdPhi"  ,   &hjjzmmdPhi   ,   "hjjzmmdPhi/F" );                
_outTree->Branch("e1"		,  &e1	       ,   "mass/F:pt/F:eta:phi/F");
_outTree->Branch("e2"		,  &e2	       ,   "mass/F:pt/F:eta:phi/F");
_outTree->Branch("mu1"		,  &mu2	       ,   "mass/F:pt/F:eta:phi/F:weightID/F");
_outTree->Branch("mu2"		,  &mu2	       ,   "mass/F:pt/F:eta:phi/F:weightTrig/F");
}

if(isWln_)
{
_outTree->Branch("W"		,  &W	       ,   "mass/F:pt/F:eta:phi/F");
_outTree->Branch("hjjwdPhi"     ,  &hjjwdPhi   ,   "hjjwdPhi/F" );                
_outTree->Branch("wenMt"  	,  &wenMt      ,   "wenMt/F"      );             	
_outTree->Branch("wenPt"  	,  &wenPt      ,   "wenPt/F"      );           	
_outTree->Branch("hjjwendPhi" ,  &hjjwendPhi ,   "hjjwendPhi/F" );                
_outTree->Branch("wmnMt"  	,  &wmnMt      ,   "wmnMt/F"      );                                     	
_outTree->Branch("wmnPt"  	,  &wmnPt      ,   "wmnPt/F"      );           	
_outTree->Branch("hjjwmndPhi" ,  &hjjwmndPhi ,   "hjjwmndPhi/F" );                
}
                                                              
if(isZinv_)
{
_outTree->Branch("MET"		,  &MET	         ,   "et/F:sig/F:phi/F");
//_outTree->Branch("hjjmetdPhi"   ,  &hjjmetdPhi   ,   "hjjmetdPhi/F" );   
  _outTree->Branch("minDeltaPhijetMET"		,  &minDeltaPhijetMET	         ,   "minDeltaPhijetMET/F");   
}
  const int initValue = -99999;
  // loop the events
  int ievt=0;  

  int totalcount=0;
    std::cout << "reading file " << inputFile_  <<std::endl;
    TString *tempfn = new TString(inputFile_);

    if(!tempfn->Contains("root")) {std::cout << "File does not have .root extension!"<<std::endl;}
    else{
    TFile* inFile = TFile::Open(tempfn->Data());

    // open input file (can be located on castor)

        // ----------------------------------------------------------------------
        // Second Part: 
        //
        //  * loop the events in the input file 
        //  * receive the collections of interest via fwlite::Handle
        //  * fill the histograms
        //  * after the loop close the input file
        // ----------------------------------------------------------------------
        fwlite::Event ev(inFile);
        for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt)
        {
         nofLeptons=initValue; _EventInfo.event = 0; _EventInfo.run = 0; _EventInfo.lumi=0; numJets = initValue; SumJetsPt = 0; bFlag = initValue; cFlag = initValue;
 
         //assign initial (unphysical) values
 
	// break loop if maximal number of events is reached 
//	if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
	// simple event counter
//	if(outputEvery_!=0 ? (ievt>0 && ievt%outputEvery_==0) : false); 
	
        int countj =0; int counth=0;

    	fwlite::Handle<unsigned int > NofJets;
	NofJets.getByLabel(ev,"JetLepEdmNtuple", "NofJets");
        if (NofJets.isValid())  
        numJets =(int) * NofJets; 

    	fwlite::Handle<unsigned int > Events;
	Events.getByLabel(ev,"HjjEdmNtuple", "hjjEventNumber");
        if (Events.isValid())  
        _EventInfo.event = * Events; 

    	fwlite::Handle<unsigned int > Lumis;
	Lumis.getByLabel(ev,"HjjEdmNtuple", "hjjLumiBlock");
        if (Lumis.isValid())  
        _EventInfo.lumi = * Lumis; 

    	fwlite::Handle<unsigned int > Runs;
	Runs.getByLabel(ev,"HjjEdmNtuple", "hjjRunNumber");
        if (Runs.isValid())  
        _EventInfo.run = * Runs; 

//      std::cout << " RUN " << run << " LUMI " << lumi << " EVENT " << event <<std::endl;

        unsigned int nofMu = 0;
        unsigned int nofEle = 0;

    	fwlite::Handle< unsigned int > nofMuons;
	nofMuons.getByLabel(ev,"JetLepEdmNtuple", "NofMuons");
        if (nofMuons.isValid()) 
	  nofMu =  *nofMuons ;


    	  fwlite::Handle< unsigned int > nofElectrons;
	  nofElectrons.getByLabel(ev,"JetLepEdmNtuple", "NofElectrons");
          if (nofElectrons.isValid()) 
	  nofEle =  *nofElectrons;


          nofLeptons = nofMu + nofEle;

          nHjj = 0;


	  // MET info
	  MET.et = initValue, MET.sig = initValue, MET.phi = initValue;


	  fwlite::Handle<std::vector<float> > pfmets;
	  pfmets.getByLabel(ev,"MetEdmNtuple", "pfMetEt");
	  MET.et = (*pfmets)[0];
	  


	  fwlite::Handle<std::vector<float> > pfmetsigs;
	  pfmetsigs.getByLabel(ev,"MetEdmNtuple", "pfMetMEtSig");
	  MET.sig = (*pfmetsigs)[0];


	  float pfmetphi=initValue;
	  fwlite::Handle<std::vector<float> > pfmetphis;
	  pfmetphis.getByLabel(ev,"MetEdmNtuple", "pfMetPhi");
	  MET.phi = (*pfmetphis)[0];


	  // Handle to collection
	  fwlite::Handle<std::vector<float> > objs;

          objs.getByLabel(ev,"HjjEdmNtuple", "hjjMass");


          if (objs.isValid()) 
          {
            nHjj = objs->size();
          }
          counters[sampleInfo.sample]++;
          sampleInfo.count = counters[sampleInfo.sample];


          
	  fwlite::Handle<std::vector<float> > jetspt;
	  jetspt.getByLabel(ev,"JetEdmNtuple", "jetsPt");
	  for(unsigned int j = 0; j < jetspt->size(); j++)
             SumJetsPt += (*jetspt)[j];






          float csv1 =0; float csv2=0;
          float bestCSV = -99999;
          unsigned int BestPairIndex=0;

          minDeltaPhijetMET = 99999;
  	   float jpt1 = -1 ,jpt2 = -1, bpt1=-1, bpt2=-1;    
	  // looking for best CSV pair, and jet closer to the MET, and saving pt of the two b-jets

	  for(unsigned int i = 0; i < nHjj; i++) 
          {
	    // best CSV pair
	  fwlite::Handle<std::vector<float> > hjjJet1CSVs;
	  hjjJet1CSVs.getByLabel(ev,"HjjEdmNtuple", "hjjJet1CSV");
	  csv1 = (*hjjJet1CSVs)[i];

	  fwlite::Handle<std::vector<float> > hjjJet2CSVs;
	  hjjJet2CSVs.getByLabel(ev,"HjjEdmNtuple", "hjjJet2CSV");
	  csv2 = (*hjjJet2CSVs)[i];

	  fwlite::Handle<std::vector<float> > hjjJet1Pts;
	  hjjJet1Pts.getByLabel(ev,"HjjEdmNtuple", "hjjJet1Pt");
	  jpt1  = (*hjjJet1Pts)[i];
	  

	  fwlite::Handle<std::vector<float> > hjjJet2Pts;
	  hjjJet2Pts.getByLabel(ev,"HjjEdmNtuple", "hjjJet2Pt");
	  jpt2  = (*hjjJet2Pts)[i];

        
          if(csv1+csv2 > bestCSV) {BestPairIndex = i; bestCSV = csv1 + csv2;  bpt1 = jpt1;  bpt2 = jpt2; }
	  // closer jet to MET 
	  fwlite::Handle<std::vector<float> > hjjJet1Phis;
	  hjjJet1Phis.getByLabel(ev,"HjjEdmNtuple", "hjjJet1Phi");
	  float phi1 = (*hjjJet1Phis)[i];

	  if (  fabs(deltaPhi(phi1, MET.phi))< minDeltaPhijetMET) {
                minDeltaPhijetMET= fabs(deltaPhi(phi1, MET.phi))  ;
	  }
	  fwlite::Handle<std::vector<float> > hjjJet2Phis;
	  hjjJet2Phis.getByLabel(ev,"HjjEdmNtuple", "hjjJet2Phi");
	  float phi2 = (*hjjJet2Phis)[i];
         
	  if (  fabs(deltaPhi(phi2, MET.phi))< minDeltaPhijetMET) {
	    minDeltaPhijetMET = fabs(deltaPhi(phi2, MET.phi))  ;
	  }
	  
	  
          }

	  
          ThirdJetPt = -1;
         
          // now filling the third central jet pt, so for jet not in the best csv pair  
           for(unsigned int i = 0; i < nHjj; i++) 
	     {
	       fwlite::Handle<std::vector<float> > hjjJet1Pts;
	       hjjJet1Pts.getByLabel(ev,"HjjEdmNtuple", "hjjJet1Pt");
	       jpt1  = (*hjjJet1Pts)[i];
	       
               if ( (jpt1 ==  bpt1) || (jpt1 ==  bpt2)) { 
		 ;   } 
	       else {
		 if (jpt1>ThirdJetPt  )  ThirdJetPt = jpt1;
	       }
                  
	       fwlite::Handle<std::vector<float> > hjjJet2Pts;
	       hjjJet2Pts.getByLabel(ev,"HjjEdmNtuple", "hjjJet2Pt");
	       jpt2  = (*hjjJet2Pts)[i];
              
               if ( (jpt2 ==  bpt1) || (jpt2 ==  bpt2)) { 
		 ;   } 
	       else {
		 if (jpt2>ThirdJetPt  )  ThirdJetPt = jpt2;
	       }
                  
	        

	     }
          
	   

	  for(unsigned int i = 0; i < nHjj; i++) 
          {

          if(BestPairIndex == i) BestCSVPair = 1;
          else BestCSVPair = 0;
          // filling now only the Best CSV pair....
          if (BestCSVPair == 0) continue; 

          totalcount++;
         //assign initial (unphysical) values
         H.mass = initValue;  H.pt = initValue; jjdr = initValue; jjdPhi = initValue;  jcvsmva = initValue; jbprob1 = initValue; jbprob2 = initValue;
         jTCHP1= initValue;  jcvs = initValue;  jjdistSVtx = initValue; hjjzdPhi = initValue; W.mass = initValue;
         W.pt = initValue; hjjwdPhi = initValue; zmmMass = initValue; zmmPt = initValue; hjjzmmdPhi = initValue; jTCHP2= initValue;
         Jet2.pt = initValue; Jet1.pt = initValue;  zeeMass = initValue; zeePt = initValue; hjjzeedPhi = initValue; wmnMt = initValue; wmnPt = initValue; 
         hjjwmndPhi = initValue; wenMt = initValue; wenPt = initValue; hjjwendPhi = initValue; csvmva1 = initValue; csvmva2 = initValue; Jet1.csv = initValue; Jet2.csv=initValue;
         Jet1.numTracksSV = initValue; Jet2.numTracksSV = initValue;  Jet1.cosTheta = initValue; Jet2.cosTheta = initValue;  zeePhi = initValue;   zmmPhi = initValue; zmmEta = initValue; zeeEta = initValue;  wenPhi = initValue;	  wenEta = initValue; wmnEta = initValue; wenEta = initValue; deltaPullAngle = initValue; wmnPhi = initValue; deltaPullAngleAK7 = initValue; alpgendrcc = initValue; alpgendrbb = initValue;

 

	  H.mass = (*objs)[i];
	    
	  // looking to daughter jet jet with pt>ptcut  
	  fwlite::Handle<std::vector<float> > hjjpts;
	  hjjpts.getByLabel(ev,"HjjEdmNtuple", "hjjPt");
	  H.pt = (*hjjpts)[i];

	  fwlite::Handle<std::vector<float> > hjjPhis;
	  hjjPhis.getByLabel(ev,"HjjEdmNtuple", "hjjPhi");
	  H.phi = (*hjjPhis)[i];

	  fwlite::Handle<std::vector<float> > hjjEtas;
	  hjjEtas.getByLabel(ev,"HjjEdmNtuple", "hjjEta");
	  H.eta = (*hjjEtas)[i];
	    
	  fwlite::Handle<std::vector<float> > hjjJet1Pts;
	  hjjJet1Pts.getByLabel(ev,"HjjEdmNtuple", "hjjJet1Pt");
	  Jet1.pt = (*hjjJet1Pts)[i];

	  fwlite::Handle<std::vector<float> > hjjJet2Pts;
	  hjjJet2Pts.getByLabel(ev,"HjjEdmNtuple", "hjjJet2Pt");
	  Jet2.pt = (*hjjJet2Pts)[i];
 	    
	  fwlite::Handle<std::vector<float> > hjjJet1Etas;
	  hjjJet1Etas.getByLabel(ev,"HjjEdmNtuple", "hjjJet1Eta");
	  Jet1.eta = (*hjjJet1Etas)[i];
	    
	  fwlite::Handle<std::vector<float> > hjjJet2Etas;
	  hjjJet2Etas.getByLabel(ev,"HjjEdmNtuple", "hjjJet2Eta");
	  Jet2.eta = (*hjjJet2Etas)[i];
	    
	  fwlite::Handle<std::vector<float> > hjjJet1Phis;
	  hjjJet1Phis.getByLabel(ev,"HjjEdmNtuple", "hjjJet1Phi");
	  Jet1.phi = (*hjjJet1Phis)[i];
	    
	  fwlite::Handle<std::vector<float> > hjjJet2Phis;
	  hjjJet2Phis.getByLabel(ev,"HjjEdmNtuple", "hjjJet2Phi");
	  Jet2.phi = (*hjjJet2Phis)[i];

 /*
	  fwlite::Handle<std::vector<float> > hjjJet1CSVMVAs;
	  hjjJet1CSVMVAs.getByLabel(ev,"HjjEdmNtuple", "hjjJet1CSVMVA");
	  csvmva1 = (*hjjJet1CSVMVAs)[i];

	  fwlite::Handle<std::vector<float> > hjjJet2CSVMVAs;
	  hjjJet2CSVMVAs.getByLabel(ev,"HjjEdmNtuple", "hjjJet2CSVMVA");
	  csvmva2 = (*hjjJet2CSVMVAs)[i];

	  fwlite::Handle<std::vector<float> > hjjJet1JbProbs;
	  hjjJet1JbProbs.getByLabel(ev,"HjjEdmNtuple", "hjjJet1JbProb");
	  float bprob1 = (*hjjJet1JbProbs)[i];

	  fwlite::Handle<std::vector<float> > hjjJet2JbProbs;
	  hjjJet2JbProbs.getByLabel(ev,"HjjEdmNtuple", "hjjJet2JbProb");
	  float bprob2 = (*hjjJet2JbProbs)[i];


	  fwlite::Handle<std::vector<float> > hjjJet1TCHPs;
	  hjjJet1TCHPs.getByLabel(ev,"HjjEdmNtuple", "hjjJet1TKHP");
	  float jTCHP1 = (*hjjJet1TCHPs)[i];

	  fwlite::Handle<std::vector<float> > hjjJet2TCHPs;
	  hjjJet2TCHPs.getByLabel(ev,"HjjEdmNtuple", "hjjJet2TKHP");
	  float jTCHP2 = (*hjjJet2TCHPs)[i];

*/

	  fwlite::Handle<std::vector<float> > hjjJet1CSVs;
	  hjjJet1CSVs.getByLabel(ev,"HjjEdmNtuple", "hjjJet1CSV");
	  Jet1.csv = (*hjjJet1CSVs)[i];

	  fwlite::Handle<std::vector<float> > hjjJet2CSVs;
	  hjjJet2CSVs.getByLabel(ev,"HjjEdmNtuple", "hjjJet2CSV");
	  Jet2.csv = (*hjjJet2CSVs)[i];


          Jet1.csvweight = ScaleCSV(Jet1.csv); 
          Jet2.csvweight = ScaleCSV(Jet2.csv); 
          cweightCSV  = Jet1.csvweight * Jet2.csvweight;

/*	  float  jjdistSVtx =-9;
	  //	  fwlite::Handle<std::vector<float> > jjdistSVtxs;
	  fwlite::Handle< std::vector<float>  > jjdistSVtxs;
	  jjdistSVtxs.getByLabel(ev,"SVtxEdmNtuple", "hjjDistSVtx");
	  if (jjdistSVtxs->size() > 0) 
          {
	    jjdistSVtx= (*jjdistSVtxs)[i];
	  }
*/
	  fwlite::Handle< std::vector<float>  > hjjJet1DeltaThetas;
   	  hjjJet1DeltaThetas.getByLabel(ev,"SVtxEdmNtuple", "hjjJet1DeltaTheta");
	  if((hjjJet1DeltaThetas->size() > 0)&&(hjjJet1DeltaThetas.isValid())) 
          {
	    deltaPullAngle= (*hjjJet1DeltaThetas)[i];
            if(fabs(deltaPullAngle) > 7) deltaPullAngle = initValue; 
	  }

	  fwlite::Handle< std::vector<float>  > hjjJet1DeltaThetas7;
   	  hjjJet1DeltaThetas7.getByLabel(ev,"SVtxEdmNtuple", "hjjJet1DeltaThetaAK7");
	  if((hjjJet1DeltaThetas7->size() > 0)&&(hjjJet1DeltaThetas7.isValid())) 
          {
	    deltaPullAngleAK7= (*hjjJet1DeltaThetas7)[i];
            if(fabs(deltaPullAngleAK7) > 7) deltaPullAngleAK7 = initValue; 
	  }

/*
	  fwlite::Handle< std::vector<int>  > hjjJet1NTracks;
   	  hjjJet1NTracks.getByLabel(ev,"SVtxEdmNtuple", "hjjJet1SVNofTrk");
	  if ((hjjJet1NTracks->size() > 0) && (hjjJet1NTracks.isValid())) 
          {
	    Jet1.numTracksSV = (*hjjJet1NTracks)[i];
	  }

	  fwlite::Handle< std::vector<int>  > hjjJet2NTracks;
   	  hjjJet2NTracks.getByLabel(ev,"SVtxEdmNtuple", "hjjJet1SVNofTrk");
	  if ((hjjJet2NTracks->size() > 0) && (hjjJet2NTracks.isValid())) 
          {
	  Jet2.numTracksSV = (*hjjJet2NTracks)[i];
	  }
*/
	  fwlite::Handle< std::vector<float>  > hjjCosThetas1;
   	  hjjCosThetas1.getByLabel(ev,"SVtxEdmNtuple", "hjjJet1higgsHelicity");
	  if ((hjjCosThetas1->size() > 0) && (hjjCosThetas1.isValid())) 
          {
	  Jet1.cosTheta = (*hjjCosThetas1)[i];
          }

	  fwlite::Handle< std::vector<float>  > hjjCosThetas2;
   	  hjjCosThetas2.getByLabel(ev,"SVtxEdmNtuple", "hjjJet2higgsHelicity");
	  if ((hjjCosThetas2->size() > 0) && (hjjCosThetas2.isValid())) 
          {
	  Jet2.cosTheta = (*hjjCosThetas2)[i];
	  }



	  jjdr =  deltaR(Jet1.eta, Jet1.phi, Jet2.eta, Jet2.phi);
  
	  jjdPhi = fabs(deltaPhi(Jet1.phi, Jet2.phi));
		  



	  fwlite::Handle<std::vector<float> > zmmMasses;
	  zmmMasses.getByLabel(ev,"ZmmEdmNtuple", "zmmMass");

	  fwlite::Handle<unsigned int>  bFlags;
   	  bFlags.getByLabel(ev,"BMCTruthEdmNtuple", "alpgenFlavorFlagBB");
	  if(bFlags.isValid()) 
          {
	  bFlag = (*bFlags);
	  }

	  fwlite::Handle<unsigned int>  cFlags;
   	  cFlags.getByLabel(ev,"BMCTruthEdmNtuple", "alpgenFlavorFlagCC");
	  if(cFlags.isValid()) 
          {
	  cFlag = (*cFlags);
	  }


	  fwlite::Handle<float>  alpgendrccs;
   	  alpgendrccs.getByLabel(ev,"BMCTruthEdmNtuple", "alpgenDrcc");
	  if((alpgendrccs.isValid()) &&  ((*alpgendrccs) > 0)) 
          {
	  alpgendrcc = (*alpgendrccs);
	  }

	  fwlite::Handle<float>  alpgendrbbs;
   	  alpgendrbbs.getByLabel(ev,"BMCTruthEdmNtuple", "alpgenDrbb");
	  if((alpgendrbbs.isValid()) &&  ((*alpgendrbbs) > 0)) 
          {
	  alpgendrbb = (*alpgendrbbs);
	  }




          if ( zmmMasses.isValid() && isZll_) 
          {
	    
	    for  (unsigned int j = 0; j < zmmMasses->size(); ++j)
	    {
	      // choosing the Zmm closer to the zmass  
	      float zmmMassTmp = (*zmmMasses)[j];
	      
	      if (fabs(zmmMassTmp - 91) < fabs(zmmMass -91 ))
               {
		zmmMass = zmmMassTmp;
		// now get Pt and Phi 
		fwlite::Handle<std::vector<float> > zmmPts;
		zmmPts.getByLabel(ev,"ZmmEdmNtuple", "zmmPt");
		zmmPt = (*zmmPts)[j];
	        if(zmmPt == 0) zmmPt = -1;	
		fwlite::Handle<std::vector<float> > zmmPhis;
		zmmPhis.getByLabel(ev,"ZmmEdmNtuple", "zmmPhi");
		zmmPhi = (*zmmPhis)[j];
	
		fwlite::Handle<std::vector<float> > zmmEtas;
		zmmEtas.getByLabel(ev,"ZmmEdmNtuple", "zmmEta");
		zmmEta = (*zmmEtas)[j];
	
		fwlite::Handle<std::vector<float> > mu1eta;
		mu1eta.getByLabel(ev,"ZmmEdmNtuple", "zmmLeptDau1Eta");
		mu1.eta = (*mu1eta)[j];
	
		fwlite::Handle<std::vector<float> > mu1phi;
		mu1phi.getByLabel(ev,"ZmmEdmNtuple", "zmmLeptDau1Phi");
		mu1.phi = (*mu1phi)[j];
	

		fwlite::Handle<std::vector<float> > mu1pt;
		mu1pt.getByLabel(ev,"ZmmEdmNtuple", "zmmLeptDau1Pt");
		mu1.pt = (*mu1pt)[j];
	       
		fwlite::Handle<std::vector<float> > mu2eta;
		mu2eta.getByLabel(ev,"ZmmEdmNtuple", "zmmLeptDau2Eta");
		mu2.eta = (*mu2eta)[j];
	
		fwlite::Handle<std::vector<float> > mu2phi;
		mu2phi.getByLabel(ev,"ZmmEdmNtuple", "zmmLeptDau2Phi");
		mu2.phi = (*mu2phi)[j];
	

		fwlite::Handle<std::vector<float> > mu2pt;
		mu2pt.getByLabel(ev,"ZmmEdmNtuple", "zmmLeptDau2Pt");
		mu2.pt = (*mu2pt)[j];
               
                mu1.weightID = ScaleID(mu1.pt, mu1.eta);
                mu1.weightTrig = ScaleIsoHLT(mu1.pt, mu1.eta);
                
                mu2.weightID = ScaleID(mu2.pt, mu2.eta);
                mu2.weightTrig = ScaleIsoHLT(mu2.pt, mu2.eta);

                cweightID = mu1.weightID * mu2.weightID; 
                cweightTrig =  mu1.weightTrig + mu2.weightTrig  -  mu1.weightTrig * mu2.weightTrig; 


                mu1.mass = 0.105658; 
                mu2.mass = 0.105658; 
	
	      }
	    }
            if(zmmPhi != initValue)
	    hjjzmmdPhi = fabs(deltaPhi(H.phi, zmmPhi));
	    }
	  

	  // Zee
          fwlite::Handle<std::vector<float> > zeeMasses;
	  zeeMasses.getByLabel(ev,"ZeeEdmNtuple", "zeeMass");
          if ( zeeMasses.isValid() && isZll_ ) 
          { 
	    
	    for  (unsigned int j = 0; j < zeeMasses->size(); ++j){
	      // choosing the Zee closer to the zmass  
	      float zeeMassTmp = (*zeeMasses)[j];
	      
	      if (fabs(zeeMassTmp - 91) < fabs(zeeMass -91 )){
		zeeMass = zeeMassTmp;
		// now get Pt and Phi 
		fwlite::Handle<std::vector<float> > zeePts;
		zeePts.getByLabel(ev,"ZeeEdmNtuple", "zeePt");
		zeePt = (*zeePts)[j];
                if(zeePt == 0) zeePt = initValue;		
	
         	fwlite::Handle<std::vector<float> > zeePhis;
		zeePhis.getByLabel(ev,"ZeeEdmNtuple", "zeePhi");
		zeePhi = (*zeePhis)[j];
	
         	fwlite::Handle<std::vector<float> > zeeEtas;
		zeeEtas.getByLabel(ev,"ZeeEdmNtuple", "zeeEta");
		zeeEta = (*zeeEtas)[j];
		
		fwlite::Handle<std::vector<float> > e1eta;
		e1eta.getByLabel(ev,"ZeeEdmNtuple", "zeeLeptDau1Eta");
		e1.eta = (*e1eta)[j];
	
		fwlite::Handle<std::vector<float> > e1phi;
		e1phi.getByLabel(ev,"ZeeEdmNtuple", "zeeLeptDau1Phi");
		e1.phi = (*e1phi)[j];
	

		fwlite::Handle<std::vector<float> > e1pt;
		e1pt.getByLabel(ev,"ZeeEdmNtuple", "zeeLeptDau1Pt");
		e1.pt = (*e1pt)[j];
	
		fwlite::Handle<std::vector<float> > e2eta;
		e2eta.getByLabel(ev,"ZeeEdmNtuple", "zeeLeptDau2Eta");
		e2.eta = (*e2eta)[j];
	
		fwlite::Handle<std::vector<float> > e2phi;
		e2phi.getByLabel(ev,"ZeeEdmNtuple", "zeeLeptDau2Phi");
		e2.phi = (*e2phi)[j];
	

		fwlite::Handle<std::vector<float> > e2pt;
		e2pt.getByLabel(ev,"ZeeEdmNtuple", "zeeLeptDau2Pt");
		e2.pt = (*e2pt)[j];
	
	      }
	    }
	    
	   if(zeePhi != initValue)
	   hjjzeedPhi = fabs(deltaPhi(H.phi, zeePhi));
	  }


	  // Wen
	  
	  fwlite::Handle<std::vector<float> > wenMts;
	  wenMts.getByLabel(ev,"WenEdmNtuple", "wenMT");
	  
	  
          if ( wenMts.isValid() && isWln_ ) {

	    
	    //  std::cout << " wenMts->size() " <<wenMts->size() << std::endl;
	    for  (unsigned int j = 0; j < wenMts->size(); ++j){

	      // choosing the Wen closer to the zmass  
	      float wenMtTmp = (*wenMts)[j];
	      
	      if (fabs(wenMtTmp - 80) < fabs(wenMt -80 )){
		wenMt = wenMtTmp;
		// now get Pt and Phi 
		fwlite::Handle<std::vector<float> > wenPts;
		wenPts.getByLabel(ev,"WenEdmNtuple", "wenPt");
		wenPt = (*wenPts)[j];
	        if(wenPt == 0) wenPt = -1;

		fwlite::Handle<std::vector<float> > wenPhis;
		wenPhis.getByLabel(ev,"WenEdmNtuple", "wenPhi");
		wenPhi = (*wenPhis)[j];
	
		fwlite::Handle<std::vector<float> > wenEtas;
		wenEtas.getByLabel(ev,"WenEdmNtuple", "wenEta");
		wenEta = (*wenEtas)[j];
		
	      }
	    }
	    
	    wenMt = wenMt;
	    wenPt = wenPt;
            if(wenPhi != initValue)
	    hjjwendPhi =  fabs(deltaPhi(H.phi, wenPhi)); 
	  }
	  
	  // Wmn
	  fwlite::Handle<std::vector<float> > wmnMts;
	  
	  wmnMts.getByLabel(ev,"WmnEdmNtuple", "wmnMT");
	  
          if ( wmnMts.isValid() && isWln_) 
          {
	    
	    for  (unsigned int j = 0; j < wmnMts->size(); ++j){
	      // choosing the Wmn closer to the wmass  
	      float wmnMtTmp = (*wmnMts)[j];
	      
	      if (fabs(wmnMtTmp - 80) < fabs(wmnMt -80 )){
		wmnMt = wmnMtTmp;
		// now get Pt and Phi 
		fwlite::Handle<std::vector<float> >  wmnPhis;
  		wmnPhis.getByLabel(ev,"WmnEdmNtuple", "wmnPhi");
	        if (wmnPhis.isValid()) wmnPhi =  (*wmnPhis)[j];

		fwlite::Handle<std::vector<float> >  wmnEtas;
  		wmnEtas.getByLabel(ev,"WmnEdmNtuple", "wmnEta");
	        if (wmnEtas.isValid()) wmnEta =  (*wmnEtas)[j];


		fwlite::Handle<std::vector<float> > wmnPts;
		wmnPts.getByLabel(ev,"WmnEdmNtuple", "wmnPt");
	        if (wmnPts.isValid()) wmnPt = (*wmnPts)[j];
	        if(wmnPt == 0) wmnPt = -1;
	      }
	    }

             if(wmnPhi != initValue) 
             hjjwmndPhi = fabs(deltaPhi(H.phi, wmnPhi));
	  }

	  //if not isZll and we found a Z instead it discard the event (needed to suppress ttbar....):
           
	  // now filling hjj plots if at leat a Z/W  has been found
	  // discard bcases in which we found a W and Z at the same time, it helps for ttbar subtraction, but we need also to apply a smarter jet veto 
	  // if is ZinvH don't care... and fill anyway




          if(isWln_) 
          {
            if((wenMt > 0) || (wmnMt > 0))
            {
               if((fabs(wenMt - 80)) < fabs(wmnMt - 80)) {W.mass = wenMt; W.pt = wenPt; hjjwdPhi = hjjwendPhi;}
               else {W.mass = wmnMt; W.pt = wmnPt; hjjwdPhi = hjjwmndPhi;}
            }
          }


          weight = cweightID * cweightTrig * cweightCSV;
	  _outTree->Fill();

	
	  }
    }


     // close input file
    inFile->Close();
  
  // break loop if maximal number of events is reached:
  // this has to be done twice to stop the file loop as well
//  if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;

std::cout << "Events: " << ievt <<std::endl;
std::cout << "TotalCount: " << totalcount <<std::endl;
_outFile->cd();

_outTree->Write();
_outFile->Write();
_outFile->Close();
return 0;
}
}
