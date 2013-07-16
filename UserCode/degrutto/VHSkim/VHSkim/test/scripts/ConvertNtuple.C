#include <iostream>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "THStack.h"
#include "TLatex.h"
#include "TStyle.h"
#include <iomanip>
#include <math.h>
#include <fstream>




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
 

float numCentralJets(float ni)
{
float n = (float) ni;
//return ((n+1)*n/2);
return ((1+sqrt(1+8*n))/2) ;
}


void Process(TString file)
{

  _METInfo MET;
  _JetInfo Jet1,Jet2, addJet1, addJet2;
  _TrackInfo H;
  _MuInfo mu1,mu2;
 
  float nHjj,jjdr,jjdPhi,jcvsmva,jbprob1,jbprob2,jTCHP1,jTCHP2,jcvs,jjdistSVtx,zmmMass,zmmPt,hjjzmmdPhi,zeeMass,zeePt,hjjzeedPhi,wmnMt,wmnPt,hjjwmndPhi,wenMt,wenPt,hjjwendPhi,csvmva1,csvmva2,hjjzdPhi,hjjwdPhi,zmmEta,zeeEta,wmnPhi,wmnEta,wenEta,wenPhi,zeePhi,zmmPhi,deltaPullAngle,deltaPullAngleAK7,alpgendrcc,alpgendrbb, alpgenVpt,weightTrig,ThirdJetPt, minDeltaPhijetMET,jetPt_minDeltaPhijetMET,PUweight,weight;
   bool isMET80_CJ80, ispfMHT150, isMET80_2CJ20,isMET65_2CJ20, isJETID,isIsoMu17;
   int nofLeptons15,nofLeptons20,  channel,numJets,numBJets, BestCSVPair;
 

  TFile *f = new TFile(file.Data(),"read");
  TTree *inTree = (TTree*) f->Get("tree");
  
  
  inTree->SetBranchAddress("H"		,  &H	                    );
  
  inTree->SetBranchAddress("Jet1"		,  &Jet1	    );
  inTree->SetBranchAddress("Jet2"		,  &Jet2	    );
  inTree->SetBranchAddress("addJet1"		,  &addJet1	    );
  inTree->SetBranchAddress("addJet2"		,  &addJet2	    );
  inTree->SetBranchAddress("nHjj" 	,  &nHjj           	    );
  inTree->SetBranchAddress("jjdr" 	,  &jjdr          	    );
  inTree->SetBranchAddress("jjdPhi"  	,  &jjdPhi          	    );
  inTree->SetBranchAddress("jcvsmva"  	,  &jcvsmva        	    );
  //inTree->SetBranchAddress("jjdistSVtx"   ,  &jjdistSVtx     	    );
  inTree->SetBranchAddress("numJets"      ,  &numJets         	    );
  inTree->SetBranchAddress("numBJets"      ,  &numBJets        	    );
  inTree->SetBranchAddress("nofLeptons15"   ,  &nofLeptons15        );
  inTree->SetBranchAddress("nofLeptons20"   ,  &nofLeptons20        );
  inTree->SetBranchAddress("deltaPullAngle", &deltaPullAngle  	    );
  inTree->SetBranchAddress("alpgendrcc"    , &alpgendrcc     	    );
  inTree->SetBranchAddress("alpgendrbb"    , &alpgendrbb     	    );
  inTree->SetBranchAddress("alpgenVpt"    , &alpgenVpt     	    );
  inTree->SetBranchAddress("weightTrig"        , &weightTrig        );
  inTree->SetBranchAddress("ThirdJetPt", &ThirdJetPt 	            );
  inTree->SetBranchAddress("deltaPullAngleAK7", &deltaPullAngleAK7  );
  inTree->SetBranchAddress("BestCSVPair",       &BestCSVPair 	    );
  inTree->SetBranchAddress("PUweight",       &PUweight 	            );
  inTree->SetBranchAddress("hjjzdPhi"     ,  &hjjzdPhi 	            );
  inTree->SetBranchAddress("zeeMass"  	,  &zeeMass    		    );
  inTree->SetBranchAddress("zeePt"  	,  &zeePt    	            );
  inTree->SetBranchAddress("hjjzeedPhi"  ,   &hjjzeedPhi  	    );
  inTree->SetBranchAddress("zmmMass"  	,  &zmmMass   	            );
  inTree->SetBranchAddress("zmmPt"  	,  &zmmPt     		    );
  inTree->SetBranchAddress("hjjzmmdPhi"  ,   &hjjzmmdPhi 	    );
  mu1.mass = 0.105658;
  mu2.mass = 0.105658;
  inTree->SetBranchAddress("mu1"		,  &mu1	       	    );
  inTree->SetBranchAddress("mu2"		,  &mu2	       	    );
  inTree->SetBranchAddress("MET"		,  &MET	            );
  //inTree->SetBranchAddress("isIsoMu17", &isIsoMu17		    );
 


  TString FOut("Repro-");
  FOut.Append(file.Data());
  TFile *_outFile	= new TFile(FOut.Data(), "recreate");	
 

 
  TTree *_outTree = new TTree("MyTuple1", "MyTuple1");







float ElecFlag=0;
float NumberJets=0;
float isB = 0;
float dEta= 0;

weight = 1.0/34.85;

_outTree->Branch("DeltaPhiElec_dij", &hjjzeedPhi    ,  "DeltaPhiElec_dij/F");  
_outTree->Branch("bbPt_dij",    &H.pt    ,    "bbPt_dij/F");
_outTree->Branch("btag1_dij",    &Jet1.csv   ,    "btag1_dij/F");
_outTree->Branch("btag2_dij",    &Jet2.csv   ,    "btag2_dij/F");
_outTree->Branch("DeltaEtaJ1J2_dij",    &dEta   ,    "DeltaEtaJ1J2_dij/F");
_outTree->Branch("JetDau1Pt_dij ",    &Jet1.pt   ,    "JetDau1Pt_dij/F");//leading jet pT)
_outTree->Branch("JetDau2Pt_dij ",    &Jet2.pt   ,    "JetDau2Pt_dij/F");//second jet pT)
_outTree->Branch("DeltaRbb_dij ",    &jjdr   ,    "DeltaRbb_dij/F");//deltaR between b jets in H cand)
//_outTree->Branch(("hely1_dij ",    &XXXX   ,    "hely1_dij/F");//i think this is cosTheta*)
_outTree->Branch("MET_dij",    &MET.et   ,    "MET_dij/F");
//_outTree->Branch("DeltaThetaFlow_dij",    &XXXX   ,    "DeltaThetaFlow_dij/F");
_outTree->Branch("NAddJet1_dij ",    &NumberJets   ,    "NAddJet1_dij/F");//num of ADDITIONAL jets (aside from bb) that are greater than 20 GeV in pT)
_outTree->Branch("ZMassElec_dij",    &zeeMass    ,    "ZMassElec_dij/F");
_outTree->Branch("ElecFlag_dij ",    &ElecFlag   ,    "ElecFlag_dij/F");//is there a Z->ee candidate?  this is 0 if it was a Z->mumu instead, 1 if Z->ee.  has to be one or the other to be in ntuple!)
//_outTree->Branch(("DeltaRbbGen_dij ",    &hjjzeedPXXXX   ,    "DeltaRbbGen_dij/F");//truth-level deltaR between b jets in H cand... not important for non-alpgen MC)
_outTree->Branch("genbevt_dij ",    &isB   ,    "genbevt_dij/F");
_outTree->Branch("bbMass_dij",      &H.mass	    ,  "bbMass_dij/F");
_outTree->Branch("diElecPt_dij"  ,  &zeePt	    ,  "diElecPt_dij/F");
_outTree->Branch("weight"  ,  &weight	    ,  "weight/F");
//_outTree->Branch("isIsoMu17", &isIsoMu17, "isIsoMu17/b"); 

  for(int jentry=0;jentry<inTree->GetEntries() ;jentry++)
  { 
    ElecFlag=0;
    NumberJets=0;
    isB = 0;
    dEta = fabs(Jet1.eta - Jet2.eta);


   inTree->GetEntry(jentry);

   NumberJets = numCentralJets(nHjj) - 2;
   if(zeeMass > 0) ElecFlag =  1; 
   if( fabs(Jet1.flavour) == 5 || fabs(Jet2.flavour) ==5 ) isB = 1;
   _outTree->Fill();
  }


_outFile->Write();
_outFile->Close();
}


void ConvertNtuple()
{

TString fname("42-ZJets100Mad.root");
Process(fname);
}
