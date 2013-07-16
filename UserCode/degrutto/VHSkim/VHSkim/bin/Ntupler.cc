#include <TH1F.h>
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
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

TFile *scaleFile = new TFile ("IsoToHLT42.root","read");
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

TFile *scaleFile = new TFile ("ScaleEffs42.root","read");
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
  gROOT->Reset();

  TTree *_outTree;
 float nHjj,jjdr,jjdPhi,jcvsmva,jbprob1,jbprob2,jTCHP1,jTCHP2,jcvs,jjdistSVtx,zmmMass,zmmPt,hjjzmmdPhi,zeeMass,zeePt,hjjzeedPhi,wmnMt,wmnPt,hjjwmndPhi,wenMt,wenPt,hjjwendPhi,csvmva1,csvmva2,hjjzdPhi,hjjwdPhi,zmmEta,zeeEta,wmnPhi,wmnEta,wenEta,wenPhi,zeePhi,zmmPhi,deltaPullAngle,deltaPullAngleAK7,alpgendrcc,alpgendrbb, alpgenVpt,weightTrig,ThirdJetPt, minDeltaPhijetMET,  jetPt_minDeltaPhijetMET , PUweight;
   int nofLeptons15,nofLeptons20,  channel,numJets,numBJets, BestCSVPair;
   bool isMET80_CJ80, ispfMHT150, isMET80_2CJ20,isMET65_2CJ20, isJETID,isIsoMu17;
 
  typedef struct 
  {
    float mass;  //MT in case of W
    float pt;
    float eta;
    float phi;
  } _TrackInfo;
  
   typedef struct 
  {
    float mass;  //MT in case of W
    float pt;
    float eta;
    float phi;
    float aodCombRelIso;
    float pfCombRelIso;
    float photonIso;
    float neutralHadIso;
    float chargedHadIso;
    float particleIso;
    float dxy;
    float dz;
  } _MuInfo;
  
  typedef struct 
  {
    float et; 
    float sumet;   
    float sig;
    float phi;
  } _METInfo;
  
  typedef struct 
  {
    float mht;
    float ht;  
    float sig;
    float phi;
  } _MHTInfo;

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
    float cosTheta;
    int numTracksSV;
    float chf;
    float nhf;
    float cef;
    float nef;
    float nch;
    float nconstituents;
    float flavour;
    float genPt;
    float genEta;
    float genPhi;
    float JECUnc;

  } _JetInfo;
  
  _METInfo MET;
  _MHTInfo MHT;
  _JetInfo Jet1,Jet2, addJet1, addJet2;
  _TrackInfo H;
  _MuInfo mu1,mu2;
  
  
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
  
  
//  std::vector<std::string> inputFiles_( in.getParameter<std::vector<std::string> >("fileNames") );
  std::string inputFile( in.getParameter<std::string> ("fileName") );

  
  bool isZinv_( ana.getParameter<bool>("isZinv"));  
  bool isWln_( ana.getParameter<bool>("isWln") );  
  bool doJID_( ana.getParameter<bool>("doJID") );  
  bool isZll_( ana.getParameter<bool>("isZll") );
  std::string PUmcfileName_( in.getParameter<std::string> ("PUmcfileName") );
  std::string PUdatafileName_( in.getParameter<std::string> ("PUdatafileName") );
  bool isMC_( ana.getParameter<bool>("isMC") );  
  
    TFile *_outFile	= new TFile(outputFile_.c_str(), "recreate");	
  _outTree = new TTree("tree", "myTree");
  
  _outTree->Branch("H"		,  &H	            ,  "mass/F:pt/F:eta:phi/F");
  
  _outTree->Branch("Jet1"		,  &Jet1	    ,  "pt/F:eta/F:phi/F:csv/F:cosTheta/F:numTracksSV/I:chf/F:nhf:cef:nef:nch:nconstituents:flavour:genPt:genEta:genPhi:JECUnc/F");
  _outTree->Branch("Jet2"		,  &Jet2	    ,  "pt/F:eta/F:phi/F:csv/F:cosTheta/F:numTracksSV/I:chf/F:nhf:cef:nef:nch:nconstituents:flavour:genPt:genEta:genPhi:JECUnc/F");
  _outTree->Branch("addJet1"		,  &addJet1	    ,  "pt/F:eta/F:phi/F:csv/F:cosTheta/F:numTracksSV/I:chf/F:nhf:cef:nef:nch:nconstituents:flavour:genPt:genEta:genPhi:JECUnc/F");
  _outTree->Branch("addJet2"		,  &addJet2	    ,  "pt/F:eta/F:phi/F:csv/F:cosTheta/F:numTracksSV/I:chf/F:nhf:cef:nef:nch:nconstituents:flavour:genPt:genEta:genPhi:JECUnc/F");
  _outTree->Branch("nHjj" 	,  &nHjj            ,  "nHjj/F"         );         	
  _outTree->Branch("jjdr" 	,  &jjdr            ,  "jjdr/F"         );         	
  _outTree->Branch("jjdPhi"  	,  &jjdPhi          ,  "jjdPhi/F"       );            	
  _outTree->Branch("jcvsmva"  	,  &jcvsmva         ,  "jcvsmva/F"      );             	
  //_outTree->Branch("jjdistSVtx"   ,  &jjdistSVtx      ,  "jjdistSVtx/F"    );                
  _outTree->Branch("numJets"      ,  &numJets         ,  "numJets/I"       );                
  _outTree->Branch("numBJets"      ,  &numBJets         ,  "numBJets/I"       );                
  _outTree->Branch("nofLeptons15"   ,  &nofLeptons15      ,  "nofLeptons15/I"    );                
  _outTree->Branch("nofLeptons20"   ,  &nofLeptons20      ,  "nofLeptons20/I"    );                
  _outTree->Branch("deltaPullAngle", &deltaPullAngle  ,  "deltaPullAngle/F");
  _outTree->Branch("alpgendrcc"    , &alpgendrcc      ,  "alpgendrcc/F");
  _outTree->Branch("alpgendrbb"    , &alpgendrbb      ,  "alpgendrbb/F");
  _outTree->Branch("alpgenVpt"    , &alpgenVpt      ,  "alpgenVpt/F");
  _outTree->Branch("weightTrig"        , &weightTrig          ,  "weightTrig/F");
  _outTree->Branch("ThirdJetPt", &ThirdJetPt  ,  "ThirdJetPt/F");
  _outTree->Branch("deltaPullAngleAK7", &deltaPullAngleAK7  ,  "deltaPullAngleAK7/F");
  _outTree->Branch("BestCSVPair",       &BestCSVPair  ,  "BestCSVPair/I");
  _outTree->Branch("PUweight",       &PUweight  ,  "PUweight/F");
  //_outTree->Branch("ThirdJetPt",        &ThirdJetPt  ,  "ThirdJetPt/F");
  
  if(isZll_)
    {
      _outTree->Branch("hjjzdPhi"     ,  &hjjzdPhi   ,   "hjjzdPhi/F" );                
      _outTree->Branch("zeeMass"  	,  &zeeMass      ,   "zeeMass/F"    );             	
      _outTree->Branch("zeePt"  	,  &zeePt        ,   "zeePt/F"      );           	
      _outTree->Branch("hjjzeedPhi"  ,   &hjjzeedPhi   ,   "hjjzeedPhi/F" );                
      _outTree->Branch("zmmMass"  	,  &zmmMass      ,   "zmmMass/F"    );                                     	
      _outTree->Branch("zmmPt"  	,  &zmmPt        ,   "zmmPt/F"      );           	
      _outTree->Branch("hjjzmmdPhi"  ,   &hjjzmmdPhi   ,   "hjjzmmdPhi/F" );                
      mu1.mass = 0.105658;
      mu2.mass = 0.105658;
      _outTree->Branch("mu1"		,  &mu1	       ,   "mass/F:pt/F:eta:phi/F:aodCombRelIso/F:pfCombRelIso/F:photonIso/F:neutralHadIso/F:chargedHadIso/F:particleIso/F:dxy/F:dz/F");
      _outTree->Branch("mu2"		,  &mu2	       ,   "mass/F:pt/F:eta:phi/F:aodCombRelIso/F:pfCombRelIso/F:photonIso/F:neutralHadIso/F:chargedHadIso/F:particleIso/F:dxy/F:dz/F");
      _outTree->Branch("MET"		,  &MET	         ,   "et/F:sumet:sig/F:phi/F");
      _outTree->Branch("isIsoMu17", &isIsoMu17, "isIsoMu17/b"); 
    }
  
  if(isWln_)
    {
      _outTree->Branch("hjjwdPhi"     ,  &hjjwdPhi   ,   "hjjwdPhi/F" );                
      _outTree->Branch("wenMt"  	,  &wenMt      ,   "wenMt/F"      );             	
      _outTree->Branch("wenPt"  	,  &wenPt      ,   "wenPt/F"      );           	
      _outTree->Branch("hjjwendPhi" ,  &hjjwendPhi ,   "hjjwendPhi/F" );                
      _outTree->Branch("wmnMt"  	,  &wmnMt      ,   "wmnMt/F"      );                                     	
      _outTree->Branch("wmnPt"  	,  &wmnPt      ,   "wmnPt/F"      );           	
      _outTree->Branch("hjjwmndPhi" ,  &hjjwmndPhi ,   "hjjwmndPhi/F" );                
      _outTree->Branch("MET"		,  &MET	         ,   "et/F:sumet:sig/F:phi/F");
    }
  
  if(isZinv_)
    {
      _outTree->Branch("MET"		,  &MET	         ,   "et/F:sumet:sig/F:phi/F");
      _outTree->Branch("MHT"		,  &MHT	         ,   "mht/F:ht:sig/F:phi/F");
      //_outTree->Branch("hjjmetdPhi"   ,  &hjjmetdPhi   ,   "hjjmetdPhi/F" );   
      _outTree->Branch("minDeltaPhijetMET"		,  &minDeltaPhijetMET	         ,   "minDeltaPhijetMET/F");
      _outTree->Branch("jetPt_minDeltaPhijetMET"		,  &jetPt_minDeltaPhijetMET	         ,   "jetPt_minDeltaPhijetMET/F");
      _outTree->Branch("isMET80_CJ80", &isMET80_CJ80, "isMET80_CJ80/b"); 
      _outTree->Branch("isMET80_2CJ20", &isMET80_2CJ20, "isMET80_2CJ20/b"); 
      _outTree->Branch("isMET65_2CJ20", &isMET65_2CJ20, "isMET65_2CJ20/b"); 
      _outTree->Branch("ispfMHT150", &ispfMHT150, "ispfMHT150/b"); 
      _outTree->Branch("isJETID", &isJETID, "isJETID/b"); 
    }
  const int initValue = -99999;
  // loop the events
  int ievt=0;  
  int totalcount=0;
  TFile* inFile = new TFile(inputFile.c_str(), "read");
      
      
      
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

 	  if(isMC_){
 	  // PU weights
 
 	  edm::LumiReWeighting   LumiWeights_ = edm::LumiReWeighting(PUmcfileName_,PUdatafileName_ , "pileup", "pileup");
  
 	  double avg=0;
          fwlite::Handle< unsigned int >  PUintimeSizes;
          PUintimeSizes.getByLabel(ev,"BMCTruthEdmNtuple", "PUintimeSize");
 
 	  fwlite::Handle< unsigned int >  PUouttime1minusSizes;
 	  PUouttime1minusSizes.getByLabel(ev,"BMCTruthEdmNtuple", "PUouttime1minusSize");
 
 	  fwlite::Handle< unsigned int >  PUouttime1plusSizes;
 	  PUouttime1plusSizes.getByLabel(ev,"BMCTruthEdmNtuple", "PUouttime1plusSize");

 	   /*
 	   fwlite::Handle< unsigned int >  PUouttime1minusSizes;
 	   PUouttime1minusSizes.getByLabel(ev,"BMCTruthEdmNtuple", "PUouttime1minusSize");
 
 	   fwlite::Handle< unsigned int >  PUouttime1plusSizes;
 	   PUouttime1plusSizes.getByLabel(ev,"BMCTruthEdmNtuple", "PUouttime1plusSize");
 	   */
 
 	   if( PUintimeSizes.isValid() && PUouttime1minusSizes.isValid() && PUouttime1plusSizes.isValid()){
 	     avg = (double)( *PUintimeSizes );
 	   }
 	   PUweight = LumiWeights_.weight3BX( avg /3.);
 	  }
 




	  nofLeptons15=initValue; nofLeptons20=initValue; _EventInfo.event = 0; _EventInfo.run = 0; _EventInfo.lumi=0; numJets = initValue; 
	  
	  //assign initial (unphysical) values
	  
	  // break loop if maximal number of events is reached 
	  if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
	  // simple event counter
	  if(outputEvery_!=0 ? (ievt>0 && ievt%outputEvery_==0) : false); 
	  
	  
	  
	  fwlite::Handle<unsigned int > NofJets;
	  NofJets.getByLabel(ev,"JetLepEdmNtuple", "NofJets");
	  if (NofJets.isValid())  numJets =(int) * NofJets; 
	  
	  fwlite::Handle<unsigned int > NofBJets;
	  NofBJets.getByLabel(ev,"JetLepEdmNtuple", "NofBJets");
	  if (NofBJets.isValid())    numBJets =(int) * NofBJets; 
	  
	  fwlite::Handle<unsigned int > Events;
	  Events.getByLabel(ev,"HjjEdmNtuple", "hjjEventNumber");
	  if (Events.isValid())  _EventInfo.event = * Events; 

	  fwlite::Handle<unsigned int > Lumis;
	  Lumis.getByLabel(ev,"HjjEdmNtuple", "hjjLumiBlock");
	  if (Lumis.isValid()) _EventInfo.lumi = * Lumis; 

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
	  
	  
          nofLeptons15 = nofMu + nofEle;
	  
	  unsigned int nofMu20 =0;
          unsigned int nofEle20 =0;
	  
	  // now looking for muons and ele up to 20 GeV  
	  fwlite::Handle<std::vector<float> > muPts;
	  muPts.getByLabel(ev,"MuEdmNtuple", "muPt");
	  
	  if (muPts.isValid() ) {
	    for (unsigned int i =0; i< muPts->size(); i++   ){
	      double mpt = (*muPts)[i];  
	      if (mpt>20) nofMu20++; 
	    }
	  }

	  fwlite::Handle<std::vector<float> > elePts;
	  elePts.getByLabel(ev,"ElecEdmNtuple", "elecPt");
          if (elePts.isValid() ) {
	    for (unsigned int i =0; i< elePts->size(); i++   ){
	      double ept = (*elePts)[i];  
	      if (ept>20) nofEle20++; 
	    }
	  }
	  nofLeptons20 = nofMu20 + nofEle20;


	  isJETID = true;
	  
	  // check if any jet passes JET id
	  if (doJID_) {
	    fwlite::Handle<std::vector<float> > Jetchfs;
	    Jetchfs.getByLabel(ev,"JetEdmNtuple", "jetschf");
	    //Jet.chf = (*hjjJet1chfs)[i];
	    
	    fwlite::Handle<std::vector<float> > Jetnhfs;
	    Jetnhfs.getByLabel(ev,"JetEdmNtuple", "jetsnhf");
	    //Jet1.nhf = (*hjjJet1nhfs)[i];
	    
	    
	    fwlite::Handle<std::vector<float> > Jetcefs;
	    Jetcefs.getByLabel(ev,"JetEdmNtuple", "jetscef");
	    //	  Jet1.cef = (*hjjJet1cefs)[i];
	    
	    
	    fwlite::Handle<std::vector<float> > Jetnefs;
	    Jetnefs.getByLabel(ev,"JetEdmNtuple", "jetsnef");
	    //Jet1.nef = (*hjjJet1nefs)[i];
	    
	    fwlite::Handle<std::vector<float> > Jetnchs;
	    Jetnchs.getByLabel(ev,"JetEdmNtuple", "jetsnch");
	    //  Jet1.nch = (*hjjJet1nchs)[i];
	    
	    fwlite::Handle<std::vector<float> > Jetnconstituentss;
	    Jetnconstituentss.getByLabel(ev,"JetEdmNtuple", "jetsnconstituents");
	    //	  Jet.nconstituents = (*hjjJet1nconstituentss)[i];
	    
	    
	    fwlite::Handle<std::vector<float> > JetEtas;
	    JetEtas.getByLabel(ev,"JetEdmNtuple", "jetsEta");
	    
	    
	    if (JetEtas.isValid() &&  Jetchfs.isValid() && Jetnhfs.isValid()  && Jetcefs.isValid() &&  Jetnefs.isValid()  && Jetnchs.isValid() && Jetnconstituentss.isValid() );
	    
	    for ( unsigned int i =0; i<JetEtas->size(); i++){
	      
	      if ( (*Jetnhfs)[i]>0.99 || (*Jetnefs)[i]>0.99 || (*Jetnconstituentss)[i]==0 )  isJETID = false;
	      
	      if ( abs( (*JetEtas)[i])>2.4 ){
		if ( (*Jetchfs)[i]==0 || (*Jetnchs)[i]==0 || (*Jetcefs)[i]>0.99 )  isJETID = false;
	      }
	      if (isJETID == false) {
		

		//		std::cout << "found a jet with eta, nhf,nef,nconstituent,chf, nch, cef: " << (*JetEtas)[i]<< ", " << (*Jetnhfs)[i] << ", " << (*Jetnefs)[i] << ", " << (*Jetnconstituentss)[i] << ", " << (*Jetchfs)[i] << ", " << (*Jetnchs)[i] << ", " << (*Jetcefs)[i] << std::endl;
	      }
	    }
	  }

	  //	  std::cout << "isJETID? " << isJETID << std::endl;
	  
          nHjj = 0;
	  
	  ispfMHT150 = false;
	  isMET80_CJ80 = false;
	  isMET80_2CJ20 = false;
 	  isMET65_2CJ20 = false;

	  // hlt bits 
	  //hltPaths = cms.untracked.vstring("HLT_IsoMu17" , "HLT_DoubleMu7", "HLT_Mu13_Mu8_v4",  "HLT_Ele27_CaloIdVT_CaloIsoT_TrkId_TrkIsoT", "HLT_Ele27_WP80_PFMHT50",  "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL",  "HLT_MET120_v",  "HLT_CentralJet80_MET80", "HLT_PFMHT150_v", "HLT_DiCentralJet20_MET80", "HLT_DiCentralJet20_BTagIP_MET65"),

     	  fwlite::Handle< std::vector<unsigned int> > hltbits; //hardcoded..... 7-->MET80_CJ80, 8 -->pfMHT150 :-(
	  hltbits.getByLabel(ev,"HLTInfoMetJetEdmNtuple", "HLTBits" );
          if (hltbits.isValid()) {
	    ispfMHT150=  (*hltbits)[7];
	    isMET80_CJ80=  (*hltbits)[8];
            isIsoMu17 = (*hltbits)[0]; 
 	    isMET80_2CJ20=  (*hltbits)[9];
 	    isMET65_2CJ20=  (*hltbits)[10];
	  }
	  
	  
	  // MET/MHT info
	  MET.et = initValue; MET.sumet =initValue ; MET.sig = initValue; MET.phi = initValue;
          MHT.mht = initValue; MHT.ht =initValue ; MHT.sig = initValue; MHT.phi = initValue;
	  
	  fwlite::Handle<std::vector<float> > pfmets;
	  pfmets.getByLabel(ev,"MetEdmNtuple", "pfMetEt");
	  MET.et = (*pfmets)[0];

	  fwlite::Handle<std::vector<float> > pfmetsumets;
	  pfmetsumets.getByLabel(ev,"MetEdmNtuple", "pfMetSumEt");
	  if (pfmetsumets.isValid()) MET.sumet = (*pfmetsumets)[0];
	  


	  fwlite::Handle<std::vector<float> > pfmetsigs;
	  pfmetsigs.getByLabel(ev,"MetEdmNtuple", "pfMetMEtSig");
	  MET.sig = (*pfmetsigs)[0];


	 
	  fwlite::Handle<std::vector<float> > pfmetphis;
	  pfmetphis.getByLabel(ev,"MetEdmNtuple", "pfMetPhi");
	  MET.phi = (*pfmetphis)[0];



	  //vector<float>             "MhtEdmNtuple"             "pfMhtMht"        "PATv2"

	  /*
	  fwlite::Handle<std::vector<float> > pfmhts;
	  pfmhts.getByLabel(ev,"MhtEdmNtuple", "pfMhtMht");
	  if (pfmhts.isValid() )MHT.mht = (*pfmhts)[0];
	  
	  fwlite::Handle<std::vector<float> > pfmhthts;
	  pfmhthts.getByLabel(ev,"MhtEdmNtuple", "pfMhtHt");
	  if (pfmhthts.isValid()) MHT.ht = (*pfmhthts)[0];
	  


	  fwlite::Handle<std::vector<float> > pfmhtsigs;
	  pfmhtsigs.getByLabel(ev,"MhtEdmNtuple", "pfMhtMhtSig");
	  if (pfmhtsigs.isValid()) MHT.sig = (*pfmhtsigs)[0];


	
	  fwlite::Handle<std::vector<float> > pfmhtphis;
	  pfmhtphis.getByLabel(ev,"MhtEdmNtuple", "pfMhtPhi");
	  if (pfmhtphis.isValid()) MHT.phi = (*pfmhtphis)[0];

	  */


	  // Handle to collection
	  fwlite::Handle<std::vector<float> > objs;

          objs.getByLabel(ev,"HjjEdmNtuple", "hjjMass");


          if (objs.isValid()) 
          {
            nHjj = objs->size();
          }

          float csv1 =0; float csv2=0;
          float bestCSV = -99999;
          unsigned int BestPairIndex=0;

          minDeltaPhijetMET = 99999;
          jetPt_minDeltaPhijetMET = 99999;
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
	  double phi1 = (*hjjJet1Phis)[i];

	  if (  fabs(deltaPhi(phi1, MET.phi))< minDeltaPhijetMET) {
                minDeltaPhijetMET= fabs(deltaPhi(phi1, MET.phi))  ;
                jetPt_minDeltaPhijetMET = jpt1;
	  }
	  fwlite::Handle<std::vector<float> > hjjJet2Phis;
	  hjjJet2Phis.getByLabel(ev,"HjjEdmNtuple", "hjjJet2Phi");
	  double phi2 = (*hjjJet2Phis)[i];
	  
	  if (  fabs(deltaPhi(phi2, MET.phi))< minDeltaPhijetMET) {
	    minDeltaPhijetMET = fabs(deltaPhi(phi2, MET.phi))  ;
            jetPt_minDeltaPhijetMET = jpt2;
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
	       jTCHP1= initValue;  jcvs = initValue;  jjdistSVtx = initValue;  hjjzdPhi = initValue; hjjwdPhi = initValue; zmmMass = initValue; zmmPt = initValue; hjjzmmdPhi = initValue; jTCHP2= initValue;
	       Jet2.pt = initValue; Jet1.pt = initValue;  zeeMass = initValue; zeePt = initValue; hjjzeedPhi = initValue; wmnMt = initValue; wmnPt = initValue; 
	       hjjwmndPhi = initValue; wenMt = initValue; wenPt = initValue; hjjwendPhi = initValue; csvmva1 = initValue; csvmva2 = initValue; Jet1.csv = initValue; Jet2.csv=initValue;
	       Jet1.numTracksSV = initValue; Jet2.numTracksSV = initValue;  Jet1.cosTheta = initValue; Jet2.cosTheta = initValue;  zeePhi = initValue;  zeeEta = initValue; zmmPhi = initValue; zmmEta = initValue;  wenPhi = initValue;	  wenEta = initValue; wmnEta = initValue; wenEta = initValue; deltaPullAngle = initValue; wmnPhi = initValue; deltaPullAngleAK7 = initValue; alpgendrcc = initValue; alpgendrbb = initValue; alpgenVpt = initValue;
	       Jet1.chf = initValue; Jet1.nhf = initValue; Jet1.cef = initValue; Jet1.nef = initValue; Jet1.nch = initValue; Jet1.nconstituents = initValue;
	       Jet2.chf = initValue; Jet2.nhf = initValue; Jet2.cef = initValue; Jet2.nef = initValue; Jet2.nch = initValue; Jet2.nconstituents = initValue;
	       addJet1.pt = initValue; addJet2.pt = initValue;
	       
	       std::vector<float> jetspt; 
	       // now we have the pt of the three central jets we are interested in, and we can loop on the jet ntuples...
	       fwlite::Handle<std::vector<float> > JetsPts;
	       JetsPts.getByLabel(ev,"JetEdmNtuple", "jetsPt");

	       fwlite::Handle<std::vector<float> > JetsEtas;
	       JetsEtas.getByLabel(ev,"JetEdmNtuple", "jetsEta");

	       fwlite::Handle<std::vector<float> > JetsPhis;
	       JetsPhis.getByLabel(ev,"JetEdmNtuple", "jetsPhi");

	       fwlite::Handle<std::vector<float> > JetsCsvs;
	       JetsCsvs.getByLabel(ev,"JetEdmNtuple", "jetsCSV");

	       fwlite::Handle<std::vector<float> > JetsGenPts;
	       JetsGenPts.getByLabel(ev,"JetLepEdmNtuple", "genJetPt");

	       fwlite::Handle<std::vector<float> > JetsGenEtas;
	       JetsGenEtas.getByLabel(ev,"JetLepEdmNtuple", "genJetEta");

	       fwlite::Handle<std::vector<float> > JetsGenPhis;
	       JetsGenPhis.getByLabel(ev,"JetLepEdmNtuple", "genJetPhi");

	       fwlite::Handle<std::vector<float> > JetsJECUncs;
	       JetsJECUncs.getByLabel(ev,"JetLepEdmNtuple", "JECUnc");

	       fwlite::Handle<std::vector<float> > Jetsflavours;
	       Jetsflavours.getByLabel(ev,"JetEdmNtuple", "jetsflavour");




	       if ( JetsPts.isValid() ) jetspt= *JetsPts;
	       
	       for (unsigned int j =0 ; j<jetspt.size() ; ++j){
		 float pt =  jetspt[j]; 
		 if ( pt == bpt1 ){
		  
		   if (JetsGenPts.isValid())  Jet1.genPt = (*JetsGenPts)[j];
		   //if (JetsGenEtas.isValid())  Jet1.genEta = (*JetsGenEtas)[j];
		   //if (JetsGenPhis.isValid())  Jet1.genPhi = (*JetsGenPhis)[j];

		   if (JetsJECUncs.isValid())  Jet1.JECUnc = (*JetsJECUncs)[j];
 
		 }
		 if ( pt == bpt2 ){
		   
		   if (JetsGenPts.isValid())  Jet2.genPt = (*JetsGenPts)[j];
		   //if (JetsGenEtas.isValid())  Jet2.genEta = (*JetsGenEtas)[j];
		   //if (JetsGenPhis.isValid())  Jet2.genPhi = (*JetsGenPhis)[j];
		   if (JetsJECUncs.isValid())  Jet2.JECUnc = (*JetsJECUncs)[j];
 
		 }
		 if (pt==ThirdJetPt){

		   if (JetsEtas.isValid())  addJet1.eta = (*JetsEtas)[j];
		   if (JetsPts.isValid())  addJet1.pt = (*JetsPts)[j];
		   if (JetsPhis.isValid())  addJet1.phi = (*JetsPhis)[j];
		   if (JetsGenPts.isValid())  addJet1.genPt = (*JetsGenPts)[j];
		   // if (JetsGenEtas.isValid())  addJet1.genEta = (*JetsGenEtas)[j];
		   //if (JetsGenPhis.isValid())  addJet1.genPhi = (*JetsGenPhis)[j];
		   if (JetsJECUncs.isValid())  addJet1.JECUnc = (*JetsJECUncs)[j];
		   if (JetsCsvs.isValid())  addJet1.csv = (*JetsCsvs)[j];
		   if (Jetsflavours.isValid())  addJet1.flavour = (*Jetsflavours)[j];
		 }
		 if (pt< ThirdJetPt && pt<bpt1 && pt<bpt2){

                   if (JetsEtas.isValid() && abs((*JetsEtas)[j])<2.4 ){
		     if (JetsEtas.isValid())  addJet2.eta = (*JetsEtas)[j];
		     if (JetsPts.isValid())  addJet2.pt = (*JetsPts)[j];
		     if (JetsPhis.isValid())  addJet2.phi = (*JetsPhis)[j];
		     if (JetsGenPts.isValid())  addJet2.genPt = (*JetsGenPts)[j];
		     //if (JetsGenEtas.isValid())  addJet2.genEta = (*JetsGenEtas)[j];
		     //if (JetsGenPhis.isValid())  addJet2.genPhi = (*JetsGenPhis)[j];
		     if (JetsJECUncs.isValid())  addJet2.JECUnc = (*JetsJECUncs)[j];
		     if (JetsCsvs.isValid())  addJet2.csv = (*JetsCsvs)[j];
		     if (Jetsflavours.isValid())  addJet2.flavour = (*Jetsflavours)[j];
		   }
		 }
	       }
	       


	       
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
	       
	       fwlite::Handle<std::vector<float> > hjjJet1flavours;
	       hjjJet1flavours.getByLabel(ev,"HjjEdmNtuple", "hjjJet1PartonFlavour");
	       if (hjjJet1flavours.isValid()) Jet1.flavour = (*hjjJet1flavours)[i];

	       fwlite::Handle<std::vector<float> > hjjJet2flavours;
	       hjjJet2flavours.getByLabel(ev,"HjjEdmNtuple", "hjjJet2PartonFlavour");
	       if (hjjJet2flavours.isValid()) Jet2.flavour = (*hjjJet2flavours)[i];
	       
               


	       // JET id
	       if (doJID_) {
		 /// fix in the future....
		 /*
		 fwlite::Handle<std::vector<float> > hjjJet1chfs;
		 hjjJet1chfs.getByLabel(ev,"HjjEdmNtuple", "hjjJet1chf");
		 Jet1.chf = (*hjjJet1chfs)[i];
		 
		 fwlite::Handle<std::vector<float> > hjjJet2chfs;
		 hjjJet2chfs.getByLabel(ev,"HjjEdmNtuple", "hjjJet2chf");
		 Jet2.chf = (*hjjJet2chfs)[i];
		 
		 fwlite::Handle<std::vector<float> > hjjJet1nhfs;
		 hjjJet1nhfs.getByLabel(ev,"HjjEdmNtuple", "hjjJet1nhf");
		 Jet1.nhf = (*hjjJet1nhfs)[i];
		 
		 fwlite::Handle<std::vector<float> > hjjJet2nhfs;
		 hjjJet2nhfs.getByLabel(ev,"HjjEdmNtuple", "hjjJet2nhf");
		 Jet2.nhf = (*hjjJet2nhfs)[i];
		 
		 fwlite::Handle<std::vector<float> > hjjJet1cefs;
		 hjjJet1cefs.getByLabel(ev,"HjjEdmNtuple", "hjjJet1cef");
		 Jet1.cef = (*hjjJet1cefs)[i];
		 
		 fwlite::Handle<std::vector<float> > hjjJet2cefs;
		 hjjJet2cefs.getByLabel(ev,"HjjEdmNtuple", "hjjJet2cef");
		 Jet2.cef = (*hjjJet2cefs)[i];
		 
		 fwlite::Handle<std::vector<float> > hjjJet1nefs;
		 hjjJet1nefs.getByLabel(ev,"HjjEdmNtuple", "hjjJet1nef");
		 Jet1.nef = (*hjjJet1nefs)[i];
		 
		 fwlite::Handle<std::vector<float> > hjjJet2nefs;
		 hjjJet2nefs.getByLabel(ev,"HjjEdmNtuple", "hjjJet2nef");
		 Jet2.nef = (*hjjJet2nefs)[i];
		 
		 fwlite::Handle<std::vector<float> > hjjJet1nchs;
		 hjjJet1nchs.getByLabel(ev,"HjjEdmNtuple", "hjjJet1nch");
		 Jet1.nch = (*hjjJet1nchs)[i];
		 
		 fwlite::Handle<std::vector<float> > hjjJet2nchs;
		 hjjJet2nchs.getByLabel(ev,"HjjEdmNtuple", "hjjJet2nch");
		 Jet2.nch = (*hjjJet2nchs)[i];
		 
		 fwlite::Handle<std::vector<float> > hjjJet1nconstituentss;
		 hjjJet1nconstituentss.getByLabel(ev,"HjjEdmNtuple", "hjjJet1nconstituents");
		 Jet1.nconstituents = (*hjjJet1nconstituentss)[i];
		 
		 fwlite::Handle<std::vector<float> > hjjJet2nconstituentss;
		 hjjJet2nconstituentss.getByLabel(ev,"HjjEdmNtuple", "hjjJet2nch");
		 Jet2.nconstituents = (*hjjJet2nchs)[i];
		 */
	       }
	       
	       
	       
	       
	       
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
	       

	       fwlite::Handle<float>  alpgenVpts;
	       alpgenVpts.getByLabel(ev,"BMCTruthEdmNtuple", "alpgenVpt");
	       if((alpgenVpts.isValid()) &&  ((*alpgenVpts) > 0)) 
		 {
		   alpgenVpt = (*alpgenVpts);
		 }
	       



	       if ( zmmMasses.isValid() ) 
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
	
			   fwlite::Handle<std::vector<float> > aodIso;
			   aodIso.getByLabel(ev,"ZmmEdmNtuple", "zmmLeptDau1aodCombRelIso");
			   mu1.aodCombRelIso = (*aodIso)[j];
	
			   fwlite::Handle<std::vector<float> > pfIso;
			   pfIso.getByLabel(ev,"ZmmEdmNtuple", "zmmLeptDau1pfCombRelIso");
			   mu1.pfCombRelIso = (*pfIso)[j];
	
		   
			   fwlite::Handle<std::vector<float> > mu1eta;
			   mu1eta.getByLabel(ev,"ZmmEdmNtuple", "zmmLeptDau1Eta");
			   mu1.eta = (*mu1eta)[j];
			   
			   fwlite::Handle<std::vector<float> > mu1phi;
			   mu1phi.getByLabel(ev,"ZmmEdmNtuple", "zmmLeptDau1Phi");
			   mu1.phi = (*mu1phi)[j];
			   
			   
			   fwlite::Handle<std::vector<float> > mu1pt;
			   mu1pt.getByLabel(ev,"ZmmEdmNtuple", "zmmLeptDau1Pt");
			   mu1.pt = (*mu1pt)[j];
	
			   fwlite::Handle<std::vector<float> > aodIso2;
			   aodIso2.getByLabel(ev,"ZmmEdmNtuple", "zmmLeptDau2aodCombrelIso");
			   mu2.aodCombRelIso = (*aodIso2)[j];
	
			   fwlite::Handle<std::vector<float> > pfIso2;
			   pfIso2.getByLabel(ev,"ZmmEdmNtuple", "zmmLeptDau2pfCombRelIso");
			   mu2.pfCombRelIso = (*pfIso2)[j];
			   
			   fwlite::Handle<std::vector<float> > mu2eta;
			   mu2eta.getByLabel(ev,"ZmmEdmNtuple", "zmmLeptDau2Eta");
			   mu2.eta = (*mu2eta)[j];
			   
			   fwlite::Handle<std::vector<float> > mu2phi;
			   mu2phi.getByLabel(ev,"ZmmEdmNtuple", "zmmLeptDau2Phi");
			   mu2.phi = (*mu2phi)[j];
			   

			   fwlite::Handle<std::vector<float> > mu2pt;
			   mu2pt.getByLabel(ev,"ZmmEdmNtuple", "zmmLeptDau2Pt");
			   mu2.pt = (*mu2pt)[j];
			   
			 }
		     }
		   if(zmmPhi != initValue)
		     hjjzmmdPhi = fabs(deltaPhi(H.phi, zmmPhi));
		 }
	       
	       
	       // Zee
	       fwlite::Handle<std::vector<float> > zeeMasses;
	       zeeMasses.getByLabel(ev,"ZeeEdmNtuple", "zeeMass");
	       if ( zeeMasses.isValid() ) 
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
               
              if(isZll_)
              {

              fwlite::Handle<std::vector<float> > muPhoton;
	      muPhoton.getByLabel(ev,"MuEdmNtuple", "muphotonIso");
	
              fwlite::Handle<std::vector<float> > muNeutral;
	      muNeutral.getByLabel(ev,"MuEdmNtuple", "muneutralHadronIso");
	
              fwlite::Handle<std::vector<float> > muCharged;
	      muCharged.getByLabel(ev,"MuEdmNtuple", "muchargeHadronIso");
	
              fwlite::Handle<std::vector<float> > muParticle;
	      muParticle.getByLabel(ev,"MuEdmNtuple", "muparticleIso");

              fwlite::Handle<std::vector<float> > dxy;
	      dxy.getByLabel(ev,"JetLepEdmNtuple", "MuDxyPV");

              fwlite::Handle<std::vector<float> > dz;
	      dz.getByLabel(ev,"JetLepEdmNtuple", "MuDzPV");


	      if (muPts.isValid() ) 
               {
	        for (unsigned int i =0; i< muPts->size(); i++   )
                {
	         double mpt = (*muPts)[i];  
	           if (mu1.pt == mpt) 
                   {
                      mu1.photonIso = (*muPhoton)[i];
                      mu1.neutralHadIso = (*muNeutral)[i];
                      mu1.particleIso = (*muParticle)[i];
                      mu1.chargedHadIso = (*muCharged)[i];
                      mu1.dxy = (*dxy)[i];
                      mu1.dz = (*dz)[i];
                   }

	           if (mu2.pt == mpt) 
                   {
                      mu2.photonIso = (*muPhoton)[i];
                      mu2.neutralHadIso = (*muNeutral)[i];
                      mu2.chargedHadIso = (*muCharged)[i];
                      mu2.particleIso = (*muParticle)[i];
                      mu2.dxy = (*dxy)[i];
                      mu2.dz = (*dz)[i];
                   }
	
	        }
	       }
                  float cweightID = ScaleID(mu1.pt,mu1.eta) * ScaleID(mu2.pt,mu2.eta) ;

                  float weightTrig1 = ScaleIsoHLT(mu1.pt,mu1.eta);
                  float weightTrig2 = ScaleIsoHLT(mu2.pt,mu2.eta);
                  float cweightTrig = weightTrig1 + weightTrig2 - weightTrig1*weightTrig2;

//                  float cweightCSV = ScaleCSV(Jet1.csv) * ScaleCSV(Jet2.csv); 

                  weightTrig = cweightID * cweightTrig;

              }

	       _outTree->Fill();
	       
	       
	     } // closed Higgs loop
	
	}// closed event loop

      std::cout << "closing the file: " << inputFile << std::endl;

	   // close input file
      inFile->Close();
     
  
    std::cout << "Events: " << ievt <<std::endl;
    std::cout << "TotalCount: " << totalcount <<std::endl;

    
    
    _outFile->cd();
    
    _outTree->Write();
    _outFile->Write();
    _outFile->Close();
    return 0;
}


