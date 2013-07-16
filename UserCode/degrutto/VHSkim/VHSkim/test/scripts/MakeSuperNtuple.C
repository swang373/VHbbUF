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
  TString file2 = file;
  _METInfo MET;
  _JetInfo Jet1,Jet2, addJet1, addJet2;
  _TrackInfo H;
  _MuInfo mu1,mu2;
 
  float nHjj,jjdr,jjdPhi,jcvsmva,jbprob1,jbprob2,jTCHP1,jTCHP2,jcvs,jjdistSVtx,zmmMass,zmmPt,hjjzmmdPhi,zeeMass,zeePt,hjjzeedPhi,wmnMt,wmnPt,hjjwmndPhi,wenMt,wenPt,hjjwendPhi,csvmva1,csvmva2,hjjzdPhi,hjjwdPhi,zmmEta,zeeEta,wmnPhi,wmnEta,wenEta,wenPhi,zeePhi,zmmPhi,deltaPullAngle,deltaPullAngleAK7,alpgendrcc,alpgendrbb, alpgenVpt,weightTrig,ThirdJetPt, minDeltaPhijetMET,jetPt_minDeltaPhijetMET,PUweight,weight,BDTTT,BDTZJl,BDTZZ,BDTZJb;
   bool isMET80_CJ80, ispfMHT150, isMET80_2CJ20,isMET65_2CJ20, isJETID,isIsoMu17;
   int nofLeptons15,nofLeptons20,  channel,numJets,numBJets, BestCSVPair;
 
  file.Append("-TRZJl.root");
  TFile *fZJl = new TFile(file.Data(),"read");

  file.ReplaceAll("TRZJl","TRZJb");
  TFile *fZJb = new TFile(file.Data(),"read");

  file.ReplaceAll("TRZJb","TRZZ");
  TFile *fZZ = new TFile(file.Data(),"read");

  file.ReplaceAll("TRZZ","TRTT");
  TFile *fTT = new TFile(file.Data(),"read");


  TTree *inTreeZJl = (TTree*) fZJl->Get("tree");
  inTreeZJl->SetBranchAddress("BDT"		,  &BDTZJl                  );
 

  TTree *inTreeZJb= (TTree*) fZJb->Get("tree");
  inTreeZJb->SetBranchAddress("BDT"		,  &BDTZJb                  );
 

  TTree *inTreeZZ = (TTree*) fZZ->Get("tree");
  inTreeZZ->SetBranchAddress("BDT"		,  &BDTZZ                   );
  inTreeZZ->SetBranchAddress("H"		,  &H	                    );

  TTree *inTreeTT = (TTree*) fTT->Get("tree");
  inTreeZJl->SetBranchAddress("BDT"		,  &BDTTT                   );



  TString FOut("TMVA-");
  FOut.Append(file2.Data());
  TFile *_outFile	= new TFile(FOut.Data(), "recreate");	
 

 
  TTree *_outTree = new TTree("tree", "tree");

_outTree->Branch("H"		,  &H	            ,  "mass/F:pt/F:eta:phi/F");
_outTree->Branch("BDTZJl"		,  &BDTZJl. "BDTZJl/F"       );
_outTree->Branch("BDTZJb"		,  &BDTZJb. "BDTZJb/F"       );
_outTree->Branch("BDTTT"		,  &BDTTT. "BDTTT/F"       );
_outTree->Branch("BDTZZ"		,  &BDTZZ. "BDTZZ/F"       );



  for(int jentry=0;jentry<inTreeZZ->GetEntries() ;jentry++)
  { 
   inTreeZZ->GetEntry(jentry);
   inTreeZJl->GetEntry(jentry);
   inTreeZJb->GetEntry(jentry);
   inTreeTT->GetEntry(jentry);

   _outTree->Fill();
  }





 
////////////////////////////////
///EXTRA BRANCHES IF NEEDED
/*////////////////////////////

  inTreeZZ->SetBranchAddress("Jet1"		,  &Jet1	    );
  inTreeZZ->SetBranchAddress("Jet2"		,  &Jet2	    );
  inTreeZZ->SetBranchAddress("addJet1"		,  &addJet1	    );
  inTreeZZ->SetBranchAddress("addJet2"		,  &addJet2	    );
  inTreeZZ->SetBranchAddress("nHjj" 	,  &nHjj           	    );
  inTreeZZ->SetBranchAddress("jjdr" 	,  &jjdr          	    );
  inTreeZZ->SetBranchAddress("jjdPhi"  	,  &jjdPhi          	    );
  inTreeZZ->SetBranchAddress("jcvsmva"  	,  &jcvsmva        	    );
  //inTreeZZ->SetBranchAddress("jjdistSVtx"   ,  &jjdistSVtx     	    );
  inTreeZZ->SetBranchAddress("numJets"      ,  &numJets         	    );
  inTreeZZ->SetBranchAddress("numBJets"      ,  &numBJets        	    );
  inTreeZZ->SetBranchAddress("nofLeptons15"   ,  &nofLeptons15        );
  inTreeZZ->SetBranchAddress("nofLeptons20"   ,  &nofLeptons20        );
  inTreeZZ->SetBranchAddress("deltaPullAngle", &deltaPullAngle  	    );
  inTreeZZ->SetBranchAddress("alpgendrcc"    , &alpgendrcc     	    );
  inTreeZZ->SetBranchAddress("alpgendrbb"    , &alpgendrbb     	    );
  inTreeZZ->SetBranchAddress("alpgenVpt"    , &alpgenVpt     	    );
  inTreeZZ->SetBranchAddress("weightTrig"        , &weightTrig        );
  inTreeZZ->SetBranchAddress("ThirdJetPt", &ThirdJetPt 	            );
  inTreeZZ->SetBranchAddress("deltaPullAngleAK7", &deltaPullAngleAK7  );
  inTreeZZ->SetBranchAddress("BestCSVPair",       &BestCSVPair 	    );
  inTreeZZ->SetBranchAddress("PUweight",       &PUweight 	            );
  inTreeZZ->SetBranchAddress("hjjzdPhi"     ,  &hjjzdPhi 	            );
  inTreeZZ->SetBranchAddress("zeeMass"  	,  &zeeMass    		    );
  inTreeZZ->SetBranchAddress("zeePt"  	,  &zeePt    	            );
  inTreeZZ->SetBranchAddress("hjjzeedPhi"  ,   &hjjzeedPhi  	    );
  inTreeZZ->SetBranchAddress("zmmMass"  	,  &zmmMass   	            );
  inTreeZZ->SetBranchAddress("zmmPt"  	,  &zmmPt     		    );
  inTreeZZ->SetBranchAddress("hjjzmmdPhi"  ,   &hjjzmmdPhi 	    );
  mu1.mass = 0.105658;
  mu2.mass = 0.105658;
  inTreeZZ->SetBranchAddress("mu1"		,  &mu1	       	    );
  inTreeZZ->SetBranchAddress("mu2"		,  &mu2	       	    );
  inTreeZZ->SetBranchAddress("MET"		,  &MET	            );
  //inTreeZZ->SetBranchAddress("isIsoMu17", &isIsoMu17		    );
 
*/
/*
 
  _outTree->Branch("Jet1"		,  &Jet1	    ,  "pt/F:eta/F:phi/F:csv/F:cosTheta/F:numTracksSV/I:chf/F:nhf:cef:nef:nch:nconstituents:flavour:genPt:genEta:genPhi:JECUnc/F");
  _outTree->Branch("Jet2"		,  &Jet2	    ,  "pt/F:eta/F:phi/F:csv/F:cosTheta/F:numTracksSV/I:chf/F:nhf:cef:nef:nch:nconstituents:flavour:genPt:genEta:genPhi:JECUnc/F");
  _outTree->Branch("addJet1"		,  &addJet1	    ,  "pt/F:eta/F:phi/F:csv/F:cosTheta/F:numTracksSV/I:chf/F:nhf:cef:nef:nch:nconstituents:flavour:genPt:genEta:genPhi:JECUnc/F");
  _outTree->Branch("addJet2"		,  &addJet2	    ,  "pt/F:eta/F:phi/F:csv/F:cosTheta/F:numTracksSV/I:chf/F:nhf:cef:nef:nch:nconstituents:flavour:genPt:genEta:genPhi:JECUnc/F");
  _outTree->Branch("nHjj" 	,  &nHjj            ,  "nHjj/F"         );         	
  _outTree->Branch("jjdr" 	,  &jjdr            ,  "jjdr/F"         );         	
  _outTree->Branch("jjdPhi"  	,  &jjdPhi          ,  "jjdPhi/F"       );            	
  _outTree->Branch("jcvsmva"  	,  &jcvsmva         ,  "jcvsmva/F"      );             	
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
  
*/








_outFile->Write();
_outFile->Close();
}


void MakeSuperNtuple()
{

///Filenames should be XX-TRYY  Where = sampleName and YY = training name
/// Names/Training = TT,ZJb,ZJl,ZZ


//TString SAMPLE("TT");
//TString SAMPLE("ZJl");
//TString SAMPLE("ZJb");
TString SAMPLE("ZZ");
Process(SAMPLE);
}
