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



#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/JetFloatAssociation.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/getRef.h"

#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "SimDataFormats/JetMatching/interface/MatchedPartons.h"
#include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <vector>

using namespace edm;
using namespace std;
using namespace reco;


class BInfoDumper : public edm::EDProducer {
public:
  BInfoDumper( const edm::ParameterSet & );
   
private:
  void produce( edm::Event &, const edm::EventSetup & );
  edm::InputTag map_;
  edm::InputTag jet_;
  edm::InputTag z_;
  edm::InputTag allJet_;
  edm::InputTag met_;


};

BInfoDumper::BInfoDumper( const ParameterSet & cfg ) : 
  map_(  cfg.getParameter<edm::InputTag> ("map")),
  jet_ (cfg.getParameter<edm::InputTag> ("jet")),
  z_ (cfg.getParameter<edm::InputTag> ("z")),
  allJet_ (cfg.getParameter<edm::InputTag> ("allJet")),
  met_(cfg.getParameter<edm::InputTag>("met"))
{

  produces< vector<float> >( "BPt" ).setBranchAlias( "BPt" );
  produces< vector<float> >( "JPt" ).setBranchAlias( "JPt" );
  produces< vector<float> >( "BEta").setBranchAlias( "BEta" );
  produces< vector<float> >( "JEta" ).setBranchAlias( "JEta" );
  produces< vector<float> >( "BPhi" ).setBranchAlias( "BPhi" );
  produces< vector<float> >( "JPhi" ).setBranchAlias( "JPhi" );
  produces< float >( "BBmass" ).setBranchAlias( "BBmass" );
  produces< float >( "Zmass" ).setBranchAlias( "Zmass" );
  produces< float >( "systemPt" ).setBranchAlias( "systemPt" );
  produces< float >( "systemPhi" ).setBranchAlias( "systemPhi" );
  produces< float >( "systemMass" ).setBranchAlias( "systemMass" );
  produces<unsigned int >( "NofJets" ).setBranchAlias( "NofJets" );
  produces<float >( "ThirdJetPt" ).setBranchAlias( "ThirdJetPt" );
  produces<float >( "ThirdJetEta" ).setBranchAlias( "ThirdJetEta" );
  produces<float >( "ThirdJetPhi" ).setBranchAlias( "ThirdJetPhi" );

}



void BInfoDumper::produce( Event & evt, const EventSetup & ) {

  auto_ptr< vector<float> > bPt( new vector<float> );  
  auto_ptr< vector<float> > bPhi (new vector<float> );  
  auto_ptr< vector<float> > bEta( new vector<float> );  


  auto_ptr< vector<float> > jPt( new vector<float> );  
  auto_ptr< vector<float> > jPhi (new vector<float> );  
  auto_ptr< vector<float> > jEta( new vector<float> );  
  auto_ptr< float > bbmass( new float );  
 auto_ptr< float > zmass( new float );  
 auto_ptr< float > systempt( new float );  
 auto_ptr< float > systemmass( new float );  
 auto_ptr< float > systemphi( new float );  
  auto_ptr< unsigned int > nofJets ( new unsigned int );
  auto_ptr< float > thirdJetPt( new float );  
  auto_ptr< float > thirdJetEta( new float );  
  auto_ptr< float > thirdJetPhi( new float );  



  std::vector< GenParticle >  bs;
Handle< GenParticleMatch >  map; 
 Handle< edm::View<pat::Jet > >  jet;
 Handle< edm::View<MET > >  met;
 Handle< edm::View<pat::Jet > >  allJet;
Handle<CandidateView>  z;

* nofJets = -1;
  *thirdJetPt = -1.;
  *thirdJetEta = -99999.;
  *thirdJetPhi = -99999.;
  *systempt = -1;
 *systemphi = -99999.;


  std::cout << "before try " << std::endl; 
 try {
        evt.getByLabel (map_ ,map );
        evt.getByLabel (jet_ ,jet );
	evt.getByLabel (allJet_ ,allJet );
	evt.getByLabel (z_ ,z );
	evt.getByLabel (met_ ,met );
    } catch(std::exception& ce) {
        cerr << "[printJetFlavour] caught std::exception " << ce.what() << endl;
        return;
    }

  std::cout << "after catch " << std::endl; 

* nofJets = allJet->size();
 std::cout << "after allJet size " << std::endl;
 std::cout << "jet size " << * nofJets << std::endl; 
 math::XYZTLorentzVector ThirdJetp4;
 
//  if (* nofJets >2) {

//   for ( unsigned int k = 0;  k < *nofJets; k++){
//     if (((jet->begin())+ k)->pt() == ((allJet->begin())+ k)->pt()) {
//       std::cout << "jet[" << k << "] with pt "<< ((allJet->begin())+ k)->pt() <<" is used for building a Higgs" << std::endl; 
//     } 
//   }


//  }

 
 for ( unsigned int  j = 0; j < jet->size(); ++j ) {
   * bbmass = -1; 

 
   const pat::Jet & aJet = (*jet)[j] ;      
   // edm::Ref<std::vector<pat::Jet> > jRef(jet, j);
   edm::RefToBase<pat::Jet> jRef = jet->refAt(j); 
   reco::GenParticleRef  aRefFlav = (*map)[jRef];
   std::cout << "is aRefFlav nonnull?? " << aRefFlav.isNonnull() << std::endl;
   // if the jet is not associated to a bjets fill the third jet entry in the event
   if  ( !aRefFlav.isNonnull()) {

*thirdJetEta = ((allJet->begin()) + j)->eta();
    *thirdJetPhi = ((allJet->begin()) + j)->phi();
 *thirdJetPt = ((allJet->begin()) + j)->pt();
    const pat::Jet & thirdJet = *((allJet->begin()) + j) ;
    ThirdJetp4 += thirdJet.p4();

 

   }
   if  ( !aRefFlav.isNonnull()) continue; 
   const reco::GenParticle &  aFlav = *aRefFlav ;

 
   if (abs(aFlav.pdgId())==5) {   printf("[printJetFlavour] (pt,eta,phi) jet = %7.2f %6.3f %6.3f | parton = %7.2f %6.3f %6.3f | %4d\n",
		               aJet.et(),
		               aJet.eta(),
		               aJet.phi(), 
		               aFlav.pt(),
		               aFlav.eta(),
		               aFlav.phi(), 
		               aFlav.pdgId()
		            );
	      bPt->push_back( aFlav.pt()) ;
	      bEta->push_back( aFlav.eta()) ;
	      bPhi->push_back( aFlav.phi()) ;
	      jPt->push_back( aJet.pt()) ;
	      jEta->push_back( aJet.eta()) ;
	      jPhi->push_back( aJet.phi()) ;
	      bs.push_back(aFlav);
	}
 }

 //take now the two B if they exist
 std::cout << " n of  b partons" << bs.size()<< std::endl;
 if (bs.size() ==2) {

   math::XYZTLorentzVector b1 ( bs[0].p4())  ;
   math::XYZTLorentzVector b2 ( bs[1].p4())  ;
   math::XYZTLorentzVector bb =  b1 + b2; 
   std::cout << "b1 mother pdgId" << bs[0].mother()->mother()->pdgId()<< std::endl;
   std::cout << "b2 mother pdgId" << bs[1].mother()->mother()->pdgId()<< std::endl;
   * bbmass =   bb.M();
 std::cout << " bbMass" << *bbmass << std::endl;
 

 // now looking at the Z and summing the H,Z system

 double trueZmass = -1;
 const Candidate * trueZ = 0  ;
  for (unsigned int i =0; i < z->size();  i++ ){
   const Candidate & zed = (*z)[i];
 
   if ( abs(zed.mass() - trueZmass) > (abs(zed.mass()  - 91)) ) {
     std::cout << "after abs" <<  std::endl;
	  trueZ = &(*z)[i];
          trueZmass = trueZ->mass();
	  std::cout << "z mass" << trueZmass << std::endl;
	}
 }

   
 if (trueZmass != -1){
   math::XYZTLorentzVector system ;
 * zmass = trueZmass;
 if ( *nofJets>2) { 
   system  = bb + trueZ->p4() + ThirdJetp4 ; //+ met->at(0).p4();
 } else {
   system  = bb + trueZ->p4() ; //+ met->at(0).p4() ;
 }
   * systempt = system.pt();
   * systemphi = system.phi();
   * systemmass = system.mass();

   std::cout << " z mass, system pt, system mass, system phi,  thirdJetPt, thirdJetEta, thirdJetphi = " << *zmass << ", "<< *systempt << ", "<< *systemmass << ", "<<*systemphi << ", "<< *thirdJetPt << ", " << *thirdJetEta << ", " <<*thirdJetPhi <<  std::endl;

 }
 }


 evt.put( bPt, "BPt" ); 
 evt.put( bEta, "BEta" ); 
 evt.put( bPhi, "BPhi" ); 
 evt.put( jPt, "JPt" ); 
 evt.put( jEta, "JEta" ); 
 evt.put( jPhi, "JPhi" ); 
 evt.put( bbmass, "BBmass" ); 
 evt.put( zmass, "Zmass" ); 
 evt.put( systempt, "systemPt" ); 
 evt.put( systemmass, "systemMass" ); 
 evt.put( systemphi, "systemPhi" ); 
 evt.put( nofJets, "NofJets" );
 evt.put( thirdJetPt, "ThirdJetPt" );
 evt.put( thirdJetEta, "ThirdJetEta" );
 evt.put( thirdJetPhi, "ThirdJetPhi" );
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( BInfoDumper );

