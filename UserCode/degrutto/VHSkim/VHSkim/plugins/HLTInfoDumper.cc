#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include <algorithm>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>

#include <TH1.h>
#include <TProfile.h>

#include <Math/VectorUtil.h>
#include <Math/GenVector/PxPyPzE4D.h>
#include <Math/GenVector/PxPyPzM4D.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include <vector>

using namespace edm;
using namespace std;
using namespace reco;
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

class HLTInfoDumper : public edm::EDProducer {
  //class HLTInfoDumper : public edm::EDAnalyzer {
public:
  HLTInfoDumper( const edm::ParameterSet & );

   
private:
  void produce(  edm::Event &, const edm::EventSetup & );
  //void beginRun(const Run& r, const EventSetup&iSetup) ;
  std::vector <std::string>  hltPaths_;
  edm::InputTag  trigTag_ ,  primaryVertices_ ;


  HLTConfigProvider  hltConfigProvider_;


};

HLTInfoDumper::HLTInfoDumper( const ParameterSet & cfg ) :
  hltPaths_(cfg.getUntrackedParameter<std::vector <std::string> >("hltPaths")),           
  trigTag_(cfg.getUntrackedParameter<edm::InputTag> ("TrigTag")),
  primaryVertices_(cfg.getParameter<InputTag>("primaryVertices")) {
  produces<int>( "numPV" ).setBranchAlias( "numPV" );
  produces<int>( "nTrkPV" ).setBranchAlias( "nTrkPV" );
  produces<float>( "chi2PV" ).setBranchAlias( "chi2PV" );
  produces<float>( "ndofPV" ).setBranchAlias( "ndofPV" );
  produces<float>( "zPV" ).setBranchAlias( "zPV" );
  produces<float>( "rhoPV" ).setBranchAlias( "rhoPV" );
  produces<vector <unsigned int > > ( "HLTBits" ).setBranchAlias( "HLTBits" );
  produces<vector <unsigned int > > ( "HLTPrescales" ).setBranchAlias( "HLTPrescales" );
  // produces<unsigned int  > ("PUintimeSize").setBranchAlias( "PUintimeSize"); 
//   produces<unsigned int  > ("PUouttime1plusSize").setBranchAlias( "PUouttime1plusSize"); 
//   produces<unsigned int  > ("PUouttime1minusSize").setBranchAlias( "PUouttime1minusSize");
//   produces<unsigned int  > ("PUouttime2plusSize").setBranchAlias( "PUouttime2plusSize");
//   produces<unsigned int  > ("PUouttime2minusSize").setBranchAlias( "PUouttime2minusSize");
  produces<float > ("rho").setBranchAlias( "rho");
}


  

  //  produces<vector <unsigned int> > ( "muHLTBits" ).setBranchAlias( "muonHLTBits" );
  // produces<vector <unsigned int> > ( "eleHLTBits" ).setBranchAlias( "eleHLTBits" );


//void HLTInfoDumper::beginRun(const Run& r, const EventSetup&iSetup) {
//  // passed as parameter to HLTConfigProvider::init(), not yet used
// bool isConfigChanged = false;
//  // isValidHltConfig_ used to short-circuit analyze() in case of problems
//  bool isValidHltConfig_ = hltConfigProvider_.init( r, iSetup, trigTag_.process(), isConfigChanged );
//
//
//std::cout << "hlt config trigger is valid??" << isValidHltConfig_ << std::endl; 

//}


void HLTInfoDumper::produce(  Event & evt, const EventSetup &iSetup ) {
  auto_ptr< vector <unsigned int> >  hltbits(new vector <unsigned int> );
  auto_ptr< vector <unsigned int> >  hltprescales(new vector <unsigned int> );
  //  auto_ptr< vector <unsigned int> >  muhltbits(new vector <unsigned int> );
  //auto_ptr< vector <unsigned int> >  elehltbits(new vector <unsigned int> );
 //  auto_ptr<unsigned int> pUintimeSize (new unsigned int); 
//   auto_ptr<unsigned int> pUouttime1plusSize (new unsigned int); 
//   auto_ptr<unsigned int> pUouttime1minusSize (new unsigned int); 
//   auto_ptr<unsigned int> pUouttime2plusSize (new unsigned int); 
//   auto_ptr<unsigned int> pUouttime2minusSize  (new unsigned int); 
 
   auto_ptr<float> rhoK6 (new float); 
 

  Handle<reco::VertexCollection> primaryVertices;  // Collection of primary Vertices
  evt.getByLabel(primaryVertices_, primaryVertices);
  auto_ptr<int> nVtxs( new int );
  auto_ptr<int> nTrkVtx( new int );
  auto_ptr<float> chi2Vtx( new float );
  auto_ptr<float> ndofVtx( new float );
  auto_ptr<float> zVtx( new float );
  auto_ptr<float> rhoVtx( new float );
  const reco::Vertex &pv = (*primaryVertices)[0];

  *nVtxs = -1;
  *nTrkVtx = -1;
  *chi2Vtx = -1.0;
  *ndofVtx = -1.0;
  *zVtx = -1000;
  *rhoVtx = -1000;
  if( !(pv.isFake()) ) {
    *nVtxs = primaryVertices->size();
    *nTrkVtx = pv.tracksSize();
    *chi2Vtx = pv.chi2();
    *ndofVtx = pv.ndof();
    *zVtx = pv.z();
    *rhoVtx = pv.position().Rho();
  }


  evt.put( nVtxs, "numPV" );
  evt.put( nTrkVtx, "nTrkPV" );
  evt.put( chi2Vtx, "chi2PV" );
  evt.put( ndofVtx, "ndofPV" );
  evt.put( zVtx, "zPV" );
  evt.put( rhoVtx, "rhoPV" );




  Handle<TriggerResults> triggerResults;
  if (!evt.getByLabel(trigTag_, triggerResults)) {
    //LogWarning("") << ">>> TRIGGER collection does not exist !!!";
    return;
  }
  
    
  evt.getByLabel(trigTag_, triggerResults); 
  const edm::TriggerNames & trigNames = evt.triggerNames(*triggerResults);
  //  LogWarning("")<<"Loop over triggers";

      bool isConfigChanged = false;
      bool isValidHltConfig_ = false; 
     
      bool found = false;
      unsigned int prescaleValue =0;
      bool fired = false;

      for(unsigned int index=0; index< hltPaths_.size() ; index++) {
      	fired =false;
	found = false;
	for (unsigned int i=0; i<triggerResults->size(); i++)  {

	const std::string trigName = trigNames.triggerName(i);
	size_t trigPath = trigName.find( hltPaths_[index]); // 0 if found, pos if not
	if (trigPath==0) found=true;
	// bit and prescales in the same order as the cfg
	//	cout << "menu contains the trigger " << trigName << " which has been found? --> " << found <<endl;
	if(found) { 
	  //bool prescaled=false;    
	  
	  isValidHltConfig_ = hltConfigProvider_.init( evt.getRun(), iSetup, trigTag_.process(), isConfigChanged );
	  if (isValidHltConfig_) {
	    for (unsigned int ps= 0; ps<  hltConfigProvider_.prescaleSize(); ps++){
	      prescaleValue = hltConfigProvider_.prescaleValue(ps, trigName) ;
	    }
	  }
	  if(  triggerResults->accept(i) ) {
	    {
              fired = true;
	      //   cout << "event " <<  evt.id().event() <<", run " <<  evt.id().run() << ", lumisection " << evt.luminosityBlock() <<"  fired the trigger " << trigName << ", with prescale value " <<  prescaleValue << endl;
	      hltbits->push_back(1);
	      hltprescales->push_back(prescaleValue);
	    }
	  } 
	  
	}
	if (found) break;
	}
	
	
	if (!fired){
	  hltbits->push_back(0);
	  hltprescales->push_back(0);
	  
	}
      }

  









 evt.put( hltbits, "HLTBits" );
 evt.put( hltprescales, "HLTPrescales" );
 
 
 // *pUintimeSize =-1;
//  *pUouttime1plusSize =-1;
//  *pUouttime1minusSize =-1;
//  *pUouttime2plusSize =-1;
//  *pUouttime2minusSize =-1;
 

// PileUp info
//   Handle<std::vector< PileupSummaryInfo > >  PupInfo;
//   pile info only in MC....
//  if (!evt.getByLabel(edm::InputTag("addPileupInfo"), PupInfo)) {
//     LogWarning("") << ">>> hjj collection does not exist !!!";
//     continue ;
//   }
//   evt.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);

//   std::vector<PileupSummaryInfo>::const_iterator PVI;
//   for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
//   if(PVI->getBunchCrossing()==0) *pUintimeSize= PVI->getPU_NumInteractions();
//   if(PVI->getBunchCrossing()==1) *pUouttime1plusSize= PVI->getPU_NumInteractions();
//   if(PVI->getBunchCrossing()==-1) *pUouttime1minusSize= PVI->getPU_NumInteractions();
//   if(PVI->getBunchCrossing()==2) *pUouttime2plusSize= PVI->getPU_NumInteractions();
//   if(PVI->getBunchCrossing()==-2) *pUouttime2minusSize= PVI->getPU_NumInteractions();
//   }

  *rhoK6  =-1;


  // rho:

  edm::Handle<double> rhoHandle;
  evt.getByLabel(edm::InputTag("kt6PFJets", "rho"),rhoHandle); // configure srcRho = cms.InputTag('kt6PFJets")
  *rhoK6 = *rhoHandle; 


 // evt.put( pUintimeSize,       "PUintimeSize" );
//  evt.put( pUouttime1plusSize, "PUouttime1plusSize" ); 
//  evt.put( pUouttime1minusSize,"PUouttime1minusSize" ); 
//  evt.put( pUouttime2plusSize, "PUouttime2plusSize" ); 
//  evt.put( pUouttime2minusSize,"PUouttime2minusSize" );

 evt.put( rhoK6,       "rho" );


}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( HLTInfoDumper );

