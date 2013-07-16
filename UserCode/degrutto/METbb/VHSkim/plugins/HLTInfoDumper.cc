#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


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

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include <vector>

using namespace edm;
using namespace std;
using namespace reco;


class HLTInfoDumper : public edm::EDProducer {
public:
  HLTInfoDumper( const edm::ParameterSet & );
   
private:
  void produce( edm::Event &, const edm::EventSetup & );
  edm::InputTag hjj_, met_ ;
  std::string hltPath_; 

};

HLTInfoDumper::HLTInfoDumper( const ParameterSet & cfg ) : 
  hjj_(cfg.getParameter<edm::InputTag>("hjj")),
  met_(cfg.getParameter<edm::InputTag>("met")),
  hltPath_(cfg.getParameter<std::string >("hltPath")) 
{
  produces<unsigned int> ( "hjjJet1HLTBit" ).setBranchAlias( "hjjJet1HLTBit" );
  produces<unsigned int> ( "hjjJet2HLTBit" ).setBranchAlias( "hjjJet2HLTBit" );

  produces<float> ( "hjjJet1HLTPt" ).setBranchAlias( "hjjJet1HLTPt" );
  produces<float> ( "hjjJet2HLTPt" ).setBranchAlias( "hjjJet2HLTPt" );

  produces<float> ( "hjjJet1HLTEta" ).setBranchAlias( "hjjJet1HLTEta" );
  produces<float> ( "hjjJet2HLTEta" ).setBranchAlias( "hjjJet2HLTEta" );

  produces<float> ( "hjjJet1HLTPhi" ).setBranchAlias( "hjjJet1HLTPhi" );
  produces<float> ( "hjjJet2HLTPhi" ).setBranchAlias( "hjjJet2HLTPhi" );

  produces<unsigned int> ( "metHLTBit" ).setBranchAlias( "metHLTBit" );


  produces<float> ( "metHLTEt" ).setBranchAlias( "metHLTEt" );


  produces<float> ( "metHLTPhi" ).setBranchAlias( "metHLTPhi" );



 
}



void HLTInfoDumper::produce( Event & evt, const EventSetup & ) {
  auto_ptr< unsigned int>  hjjjet1hltbit(new unsigned int );
  auto_ptr< unsigned int>  hjjjet2hltbit(new unsigned int );

  auto_ptr< float>  hjjjet1hltpt (new float);
  auto_ptr< float>  hjjjet2hltpt(new float );

  auto_ptr< float>  hjjjet1hlteta(new float);
  auto_ptr< float>  hjjjet2hlteta(new float );

  auto_ptr< float>  hjjjet1hltphi(new float);
  auto_ptr< float>  hjjjet2hltphi(new float );

  auto_ptr< unsigned int>  methltbit(new unsigned int );
  auto_ptr< float>  methltet(new float );
  auto_ptr< float>  methltphi(new float );

  *  hjjjet1hltbit = 0;
  *hjjjet2hltbit = 0;
  *hjjjet1hltpt = -9999;
  *hjjjet2hltpt = -9999;
  *hjjjet1hlteta = -9999;
  *hjjjet2hlteta= -9999;
  *hjjjet1hltphi = -9999;
  *hjjjet2hltphi= -9999;
  *methltbit = 0;
  *methltet = 0;
  *methltphi =-9999;


  // MET
  Handle<View<pat::MET> > metCollection;
  if (!evt.getByLabel(met_, metCollection)) {
    // LogWarning("") << ">>> MET collection does not exist !!!";
    return;
  }
  const pat::MET * met = & metCollection->at(0);
  const pat::TriggerObjectStandAloneCollection mHLTMatches =  met->triggerObjectMatchesByPath( hltPath_, false);
    unsigned int mHLTSize = mHLTMatches.size();


    //std::cout <<  "found a met with hlt size = " << mHLTSize << std::endl; 
    mHLTSize>0 ? *methltbit  = 1 : *methltbit = 0;  

    for (unsigned int  i = 0; i < mHLTSize ; i++) {
      const  pat::TriggerObject mHLT = mHLTMatches[i];
      // float mHLTDeltaPhi = deltaPhi(mHLT.phi(), met->phi());      
      *methltet = mHLT.et();      
      *methltphi = mHLT.phi();      
      // std::cout << "found a matching object[" << i << "] with deltaPhi "<<  mHLTDeltaPhi << "and et " << *methltet << "to compare with met et = " << met->pt()<< std::endl; 
    }    
    


  Handle<CandidateView > hjjColl;
  if (!evt.getByLabel(hjj_,hjjColl)) {
    // LogWarning("") << ">>> hjj collection does not exist !!!";
    return;
  }
for(CandidateView::const_iterator i = hjjColl->begin(); i != hjjColl->end(); ++i) {


  const Candidate * dau1 = i->daughter(0);
    const Candidate * dau2 = i->daughter(1);
    if(dau1 == 0|| dau2 == 0) 
      throw Exception(errors::InvalidReference) <<
	"one of the two daughter does not exist\n";
     const Candidate * c1 = dau1->masterClone().get();
    const pat::Jet * jet1 = dynamic_cast<const pat::Jet*>(c1);


    if (jet1 !=0) { const pat::TriggerObjectStandAloneCollection j1HLTMatches =  jet1->triggerObjectMatchesByPath( hltPath_, false);
    unsigned int jHLTSize = j1HLTMatches.size();
    // std::cout <<  "------found a jet1 with hlt size = " << jHLTSize << std::endl; 
    jHLTSize>0 ? *hjjjet1hltbit  = 1 : *hjjjet1hltbit = 0;  

    for (unsigned int  i = 0; i < jHLTSize ; i++) {
      const  pat::TriggerObject jHLT = j1HLTMatches[i];
      //float j1HLTDeltaR = deltaR(jHLT.eta(),  jHLT.phi(), jet1->eta(), jet1->phi());      
      *hjjjet1hltpt = jHLT.pt();      
      *hjjjet1hlteta = jHLT.eta();      
      *hjjjet1hltphi = jHLT.phi();      
      //std::cout << "found a matching object[" << i << "] with deltaR "<<  j1HLTDeltaR << "and pt " << *hjjjet1hltpt << "to compare with jet1 pt = " << jet1->pt()<< std::endl; 
    }    
    }



    const Candidate * c2 = dau2->masterClone().get();
    const pat::Jet * jet2 = dynamic_cast<const pat::Jet*>(c2);
	
   
    if (jet2 !=0) { const pat::TriggerObjectStandAloneCollection j2HLTMatches =  jet2->triggerObjectMatchesByPath( hltPath_, false);


      unsigned int jHLTSize = j2HLTMatches.size();
      // std::cout <<  "+++++++found a jet2 with hlt size = " << jHLTSize << std::endl; 
    jHLTSize>0 ? *hjjjet2hltbit  = 1 : *hjjjet2hltbit = 0;  

    for (unsigned int  i = 0; i < jHLTSize ; i++) {
      const  pat::TriggerObject jHLT = j2HLTMatches[i];
      // float j2HLTDeltaR = deltaR(jHLT.eta(), jHLT.phi(), jet2->eta(),  jet2->phi());      
      *hjjjet2hltpt = jHLT.pt();      
      *hjjjet2hlteta = jHLT.eta();      
      *hjjjet2hltphi = jHLT.phi();      
      // std::cout << "found a matching object[" << i << "] with deltaR "<<  j2HLTDeltaR << "and pt " << *hjjjet2hltpt << "to compare with jet2 pt = " << jet2->pt()<< std::endl; 
    }    
    }
	    
 }


 evt.put( hjjjet1hltbit, "hjjJet1HLTBit" );
 evt.put( hjjjet2hltbit, "hjjJet2HLTBit" );

 evt.put( hjjjet1hltpt, "hjjJet1HLTPt" );
 evt.put( hjjjet2hltpt, "hjjJet2HLTPt" );

 evt.put( hjjjet1hlteta, "hjjJet1HLTEta" );
 evt.put( hjjjet2hlteta, "hjjJet2HLTEta" );

 evt.put( hjjjet1hltphi, "hjjJet1HLTPhi" );
 evt.put( hjjjet2hltphi, "hjjJet2HLTPhi" );

 evt.put( methltbit, "metHLTBit" );
 evt.put( methltet, "metHLTEt" );
 evt.put( methltphi, "metHLTPhi" );


 
 
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( HLTInfoDumper );

