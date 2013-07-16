#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
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


int main(int argc, char* argv[]) 
{
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
    std::cout << "Usage : " << argv[0] << " [parameters.py]" << std::endl;
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
  std::vector<std::string> inputFiles_( in.getParameter<std::vector<std::string> >("fileNames") );
  // edm::InputTag muons_( ana.getParameter<edm::InputTag>("muons") );
   
  double j1pt_cut_( ana.getParameter<double>("j2pt_cut") );  
  double j2pt_cut_( ana.getParameter<double>("j2pt_cut") );  
  double jCSVMVA_cut_( ana.getParameter<double>("jCSVMVA_cut") );   
  double jbProb_cut_( ana.getParameter<double>("jbProb_cut") ); 
  double jjDeltaRmin_cut_( ana.getParameter<double>("jjDeltaRmin_cut") ); 
  double jjDeltaRmax_cut_( ana.getParameter<double>("jjDeltaRmax_cut") ); 
  double pfMet_cut_( ana.getParameter<double>("pfMet_cut") ); 
  double jjDistSVtx_cut_( ana.getParameter<double>("jjDistSVtx_cut") ); 



  // book a set of histograms
  fwlite::TFileService fs = fwlite::TFileService(outputFile_.c_str());
  TFileDirectory dir = fs.mkdir("analyzeBasicPat");
  TH1F* zjjMass_  = dir.make<TH1F>("zjjMass"  , "mass"  ,   500,   0.,  1000.);
  TH1F* jjdr_  = dir.make<TH1F>("jjdr"  , "deltaR"  ,   100,   0.,  5.);
  TH1F* jcsvmva_  = dir.make<TH1F>("jcvsmva"  , "jet csvmva"  ,   100,   0,  1);
  TH1F* jbprob_  = dir.make<TH1F>("jbprob"  , "jet bprob"  ,   100,   0,  10);
  TH1F* pfmet_  = dir.make<TH1F>("pfmet"  , "pfmet"  ,   100,   0,  100);
  TH1F* jjdistSVtx_  = dir.make<TH1F>("jjdistSVtx"  , "jjdistSVtx"  ,   100,   0,  2);


  //  TH1F* muonEta_ = dir.make<TH1F>("muonEta" , "eta" ,   100,  -3.,    3.);
  //TH1F* muonPhi_ = dir.make<TH1F>("muonPhi" , "phi" ,   100,  -5.,    5.);  
  //TH1F* mumuMass_= dir.make<TH1F>("mumuMass", "mass",    90,   30., 120.);
  
  // loop the events
  int ievt=0;  
  for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile){
    // open input file (can be located on castor)
    std::cout << "reading the file: " << inputFiles_[iFile] << std::endl;
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
      for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt){
	//edm::EventBase const & event = ev;
	// break loop if maximal number of events is reached 
	if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
	// simple event counter
	//	if(outputEvery_!=0 ? (ievt>0 && ievt%outputEvery_==0) : false); 
	// std::cout << "  processing event: " << ievt << std::endl;
	
	// Handle to collection
	fwlite::Handle<std::vector<float> > objs;
	objs.getByLabel(ev,"ZjjEdmNtuple", "zjjMass");

	for(unsigned int i = 0; i < objs->size(); ++i){
	  double mass = (*objs)[i];

	  // if(fMin < mass && mass < fMax)

	  // looking to daughter jet jet with pt>ptcut  


	  fwlite::Handle<std::vector<float> > zjjJet1Pts;
	  zjjJet1Pts.getByLabel(ev,"ZjjEdmNtuple", "zjjJet1Pt");
	  double pt1 = (*zjjJet1Pts)[i];


	  fwlite::Handle<std::vector<float> > zjjJet2Pts;
	  zjjJet2Pts.getByLabel(ev,"ZjjEdmNtuple", "zjjJet2Pt");
	  double pt2= (*zjjJet2Pts)[i];

	  fwlite::Handle<std::vector<float> > zjjJet1Etas;
	  zjjJet1Etas.getByLabel(ev,"ZjjEdmNtuple", "zjjJet1Eta");
	  double eta1 = (*zjjJet1Etas)[i];

	  fwlite::Handle<std::vector<float> > zjjJet2Etas;
	  zjjJet2Etas.getByLabel(ev,"ZjjEdmNtuple", "zjjJet2Eta");
	  double eta2= (*zjjJet2Etas)[i];

	  fwlite::Handle<std::vector<float> > zjjJet1Phis;
	  zjjJet1Phis.getByLabel(ev,"ZjjEdmNtuple", "zjjJet1Phi");
	  double phi1 = (*zjjJet1Phis)[i];

	  fwlite::Handle<std::vector<float> > zjjJet2Phis;
	  zjjJet2Phis.getByLabel(ev,"ZjjEdmNtuple", "zjjJet2Phi");
	  double phi2= (*zjjJet2Phis)[i];

	  fwlite::Handle<std::vector<float> > zjjJet1CSVMVAs;
	  zjjJet1CSVMVAs.getByLabel(ev,"ZjjEdmNtuple", "zjjJet1CSVMVA");
	  double csvmva1 = (*zjjJet1CSVMVAs)[i];

	  fwlite::Handle<std::vector<float> > zjjJet2CSVMVAs;
	  zjjJet2CSVMVAs.getByLabel(ev,"ZjjEdmNtuple", "zjjJet2CSVMVA");
	  double csvmva2 = (*zjjJet2CSVMVAs)[i];

	  fwlite::Handle<std::vector<float> > zjjJet1JbProbs;
	  zjjJet1JbProbs.getByLabel(ev,"ZjjEdmNtuple", "zjjJet1JbProb");
	  double bprob1 = (*zjjJet1JbProbs)[i];

	  fwlite::Handle<std::vector<float> > zjjJet2JbProbs;
	  zjjJet2JbProbs.getByLabel(ev,"ZjjEdmNtuple", "zjjJet2JbProb");
	  double bprob2 = (*zjjJet2JbProbs)[i];

	  fwlite::Handle<std::vector<float> > pfmets;
	  pfmets.getByLabel(ev,"MetEdmNtuple", "pfMetPt");
	  double pfmet = (*pfmets)[i];

	  //	  fwlite::Handle<std::vector<float> > jjdistSVtxs;
	  fwlite::Handle< float  > jjdistSVtxs;
	  pfmets.getByLabel(ev,"SVtxEdmNtuple", "zjjDistSVtx");

	  /*	
  //	  double  jjdistSVtx= (*jjdistSVtxs)[i];
	
  float  jjdistSVtx= *jjdistSVtxs;

	  */
	  double dR = deltaR( eta1, phi1, eta2, phi2);  
	  

	  //	  std::cout <<" size of zjj daugther1 pt"<<<<std::endl;

	  if (pt1> j1pt_cut_ && pt2> j2pt_cut_  &&  csvmva1 > jCSVMVA_cut_ &&  csvmva2>jCSVMVA_cut_ &&  bprob1>jbProb_cut_ && bprob2>jbProb_cut_ &&  dR > jjDeltaRmin_cut_ && dR < jjDeltaRmax_cut_ && pfmet< pfMet_cut_ /* && jjdistSVtx > jjDistSVtx_cut_ */ ) {
	    //  std::cout <<" size of zjj masses "<<objs->size()<<std::endl;
	    //std::cout << " mass = " << mass << std::endl;
	    //std::cout <<" size of zjj daugther0 pt "<<zjjJet1Pts->size()<<std::endl;
	    //std::cout <<" zjj daugther0 pt "<<pt1<<std::endl;
	    //std::cout <<" zjj daugther1 pt "<<pt2<<std::endl;
	    //	  jCSVMVA_cut


	    zjjMass_->Fill(mass);

	      jjdr_->Fill(dR);
	    jcsvmva_->Fill(csvmva1);
	    jcsvmva_->Fill(csvmva2);

	    jbprob_->Fill(bprob1);
	    jbprob_->Fill(bprob2);

	    pfmet_->Fill(pfmet);

	    //	    jjdistSVtx_ -> Fill (jjdistSVtx);


	  }
      }
    }


      /*
	// loop muon collection and fill histograms
	for(std::vector<Muon>::const_iterator mu1=muons->begin(); mu1!=muons->end(); ++mu1){
	  muonPt_ ->Fill( mu1->pt () );
	  muonEta_->Fill( mu1->eta() );
	  muonPhi_->Fill( mu1->phi() );	  
	  if( mu1->pt()>20 && fabs(mu1->eta())<2.1 ){
	    for(std::vector<Muon>::const_iterator mu2=muons->begin(); mu2!=muons->end(); ++mu2){
	      if(mu2>mu1){ // prevent double conting
		if( mu1->charge()*mu2->charge()<0 ){ // check only muon pairs of unequal charge 
		  if( mu2->pt()>20 && fabs(mu2->eta())<2.1 ){
		    mumuMass_->Fill( (mu1->p4()+mu2->p4()).mass() );
		  }
		}
	      }
	    }
	  }
	}
      } 
      */ 
      // close input file
      inFile->Close();
    }
    // break loop if maximal number of events is reached:
    // this has to be done twice to stop the file loop as well
    if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
  }
  return 0;
}
