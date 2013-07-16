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

int main(int argc, char* argv[]) 
{
TTree *_outTree;
float pfmet,wmnMt,wmnPt,wenMt,wenPt,wmnPhi,wmnEta,wenEta,wenPhi;

typedef struct 
{
  float mt;  
  float pt;
  float eta;
  float phi;
} _WInfo;

typedef struct 
{
  float mass;  
  float pt;
  float eta;
  float phi;
} _TrackInfo;


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

_TrackInfo MuSim,MuGen,MuRec;
_WInfo W;
_METInfo MET;

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
  bool isMC ( ana.getParameter<bool>("isMC") );
  std::string outputFile_( out.getParameter<std::string>("fileName" ) );
  std::vector<std::string> inputFiles_( in.getParameter<std::vector<std::string> >("fileNames") );
  double lumi = 1000* ana.getParameter<double>("lumi");

//weight = 1;
///////////////////////////////////////////////////////////////////////
//From Michele's MC Studies:                                         // 
//See https://twiki.cern.ch/twiki/bin/viewauth/CMS/MCBackgroundList  // 
///////////////////////////////////////////////////////////////////////

 
//if(sampleName == "ZllHbb115")            {sampleInfo.sample = 1; weight = lumi/5217391.0;}          
//if(sampleName == "Z1Jets_ptZ-0to100")    {sampleInfo.sample = 2; weight = lumi/lumiZ1jets0_100;}


TFile *_outFile	= new TFile(outputFile_.c_str(), "recreate");	
_outTree = new TTree("tree", "myTree");

_outTree->Branch("eventInfo",      &_EventInfo	    ,  "run/I:lumi/I:event/I");
_outTree->Branch("W"		,  &W	            ,   "mt/F:pt/F:eta:phi/F");
if(isMC) _outTree->Branch("MuGen"		,  &MuGen            ,   "mass/F:pt/F:eta:phi/F");
if(isMC) _outTree->Branch("MuSim"		,  &MuSim            ,   "mass/F:pt/F:eta:phi/F");
_outTree->Branch("MuRec"		,  &MuRec            ,   "mass/F:pt/F:eta:phi/F");
/*_outTree->Branch("wenMt"  	,  &wenMt      ,   "wenMt/F"      );             	
_outTree->Branch("wenPt"  	,  &wenPt      ,   "wenPt/F"      );           	
_outTree->Branch("wmnMt"  	,  &wmnMt      ,   "wmnMt/F"      );                                     	
_outTree->Branch("wmnPt"  	,  &wmnPt      ,   "wmnPt/F"      );           	*/
_outTree->Branch("MET"		,  &MET	         ,   "et/F:sig/F:phi/F");
  
  const int initValue = -99999;

  const float muMass = 0.10565837;

  // loop the events
  int ievt=0;  
  for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile){
    // open input file (can be located on castor)
    std::cout << "reading file " << iFile+1 << "/" << inputFiles_.size() << ": " << inputFiles_[iFile]  <<std::endl;
    TFile* inFile = TFile::Open(inputFiles_[iFile].c_str());
    if( inFile ){
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
 
 
         //assign initial (unphysical) values
         MET.et = initValue; MET.sig = initValue; W.eta = initValue; W.mt = initValue; W.pt = initValue; wmnMt = initValue; wmnPt = initValue;  MuSim.eta = initValue; MuSim.mass = muMass; MuSim.pt = initValue;  MuGen.eta = initValue; MuGen.mass = muMass; MuGen.pt = initValue; MuGen.phi = initValue; MuSim.phi = initValue; MuSim.eta = MuRec.eta = initValue; MuRec.mass = muMass; MuRec.pt = initValue;
         wenMt = initValue; wenPt = initValue;_EventInfo.event = 0; _EventInfo.run = 0; _EventInfo.lumi=0;  wenPhi = initValue;	  wenEta = initValue; wmnEta = initValue; wenEta = initValue;  wmnPhi = initValue;   

	// break loop if maximal number of events is reached 
	if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
	// simple event counter
	if(outputEvery_!=0 ? (ievt>0 && ievt%outputEvery_==0) : false); 
	


        if(isMC)
        {
	  fwlite::Handle<std::vector<float> > genMuEtas;
	  genMuEtas.getByLabel(ev,"GenSimMuons", "genMuEta");
	  if(genMuEtas->size() > 0) MuGen.eta = (*genMuEtas)[0];

	  fwlite::Handle<std::vector<float> > genMuPhis;
	  genMuPhis.getByLabel(ev,"GenSimMuons", "genMuPhi");
	  if(genMuPhis->size() > 0) MuGen.phi = (*genMuPhis)[0];

	  fwlite::Handle<std::vector<float> > genMuPts;
	  genMuPts.getByLabel(ev,"GenSimMuons", "genMuPt");
	 if(genMuPts->size() > 0) MuGen.pt = (*genMuPts)[0];

	  fwlite::Handle<std::vector<float> > simMuEtas;
	  simMuEtas.getByLabel(ev,"GenSimMuons", "simMuEta");
	  if(simMuEtas->size() > 0)MuSim.eta = (*simMuEtas)[0];

	  fwlite::Handle<std::vector<float> > simMuPhis;
	  simMuPhis.getByLabel(ev,"GenSimMuons", "simMuPhi");
	  if(simMuPhis->size() > 0)MuSim.phi = (*simMuPhis)[0];

	  fwlite::Handle<std::vector<float> > simMuPts;
	  simMuPts.getByLabel(ev,"GenSimMuons", "simMuPt");
	 if(simMuPts->size() > 0) MuSim.pt = (*simMuPts)[0]; 
        }

        int countj =0; int counth=0;

    	  fwlite::Handle< unsigned int > evts;
	  evts.getByLabel(ev,"WmnEdmNtuple", "wmnEventNumber");
          if (evts.isValid()) 
	    _EventInfo.event  =  *evts;
	   
 	  fwlite::Handle< unsigned int > lums;
	  lums.getByLabel(ev,"WmnEdmNtuple", "wmnLumiBlock");
          if (lums.isValid())  
            _EventInfo.lumi  =  *lums;

   	  fwlite::Handle< unsigned int > runs;
	  evts.getByLabel(ev,"WmnEdmNtuple", "wmnRunNumber");
          if (runs.isValid()) 
	    _EventInfo.run  =  *runs;
	 


	  fwlite::Handle<std::vector<float> > wenMts;
	  wenMts.getByLabel(ev,"WenEdmNtuple", "wenMT");
	  
          if ( wenMts.isValid() ) {

	    
	    //  std::cout << " wenMts->size() " <<wenMts->size() << std::endl;
	    for  (unsigned int j = 0; j < wenMts->size(); ++j){

	      // choosing the Wen closer to the mt  
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
	  }
	  
	  // Wmn
	  fwlite::Handle<std::vector<float> > wmnMts;
	  
	  wmnMts.getByLabel(ev,"WmnEdmNtuple", "wmnMT");
	  
          if ( wmnMts.isValid() ) 
          {
	    
	    for  (unsigned int j = 0; j < wmnMts->size(); ++j){
	      // choosing the Wmn closer to the wmt  
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
                

                fwlite::Handle<std::vector<float> >  wmnMEtEts;
  		wmnMEtEts.getByLabel(ev,"WmnEdmNtuple", "wmnMEtEt");
	        if (wmnMEtEts.isValid()) MET.et =  (*wmnMEtEts)[j];


                fwlite::Handle<std::vector<float> >  wmnMEtPhis;
  		wmnMEtPhis.getByLabel(ev,"WmnEdmNtuple", "wmnMEtPhi");
	        if (wmnMEtPhis.isValid()) MET.phi =  (*wmnMEtPhis)[j];


	        fwlite::Handle<std::vector<float> >  wmnMEtSigs;
  		wmnMEtSigs.getByLabel(ev,"WmnEdmNtuple", "wmnMEtSig");
	        if (wmnMEtSigs.isValid()) MET.sig =  (*wmnMEtSigs)[j];


		fwlite::Handle<std::vector<float> > wmnPts;
		wmnPts.getByLabel(ev,"WmnEdmNtuple", "wmnPt");
	        if (wmnPts.isValid()) wmnPt = (*wmnPts)[j];
	        if(wmnPt == 0) wmnPt = -1;

		fwlite::Handle<std::vector<float> > recmPts;
                recmPts.getByLabel(ev,"WmnEdmNtuple", "wmnLeptDau1Pt");
	        if (recmPts.isValid()) MuRec.pt = (*recmPts)[j];

		fwlite::Handle<std::vector<float> > recmPhis;
                recmPhis.getByLabel(ev,"WmnEdmNtuple", "wmnLeptDau1Phi");
	        if (recmPhis.isValid()) MuRec.phi = (*recmPhis)[j];

		fwlite::Handle<std::vector<float> > recmEtas;
                recmEtas.getByLabel(ev,"WmnEdmNtuple", "wmnLeptDau1Eta");
	        if (recmEtas.isValid()) MuRec.eta = (*recmEtas)[j];



	      }
	    }

	  }

	  //if not isZll and we found a Z instead it discard the event (needed to suppress ttbar....):
           
	  // now filling hjj plots if at leat a Z/W  has been found
	  // discard bcases in which we found a W and Z at the same time, it helps for ttbar subtraction, but we need also to apply a smarter jet veto 
/*
            if((wenMt > 0) || (wmnMt > 0))
            {
               if((fabs(wenMt - 80)) < fabs(wmnMt - 80)) {W.mass = wenMt; W.pt = wenPt; hjjwdPhi = hjjwendPhi;}
               else {W.mass = wmnMt; W.pt = wmnPt; hjjwdPhi = hjjwmndPhi;}
            }
*/
	  // if is ZinvH don't care... and fill anyway

            if(wmnMt > 0)
            {
               W.mt= wmnMt; W.pt = wmnPt; W.eta = wmnEta;            }

          _outTree->Fill();

	
      }


     // close input file
    inFile->Close();
 }
  // break loop if maximal number of events is reached:
  // this has to be done twice to stop the file loop as well
std::cout << "Events: " << ievt <<std::endl;
_outFile->cd();

}
_outTree->Write();
_outFile->Write();
_outFile->Close();
return 0;

}
