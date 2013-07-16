#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h" 

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

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include <vector>

using namespace edm;
using namespace std;
using namespace reco;


class JetLepInfoDumper : public edm::EDProducer {
public:
  JetLepInfoDumper( const edm::ParameterSet & );
   
private:
  void produce( edm::Event &, const edm::EventSetup & );
  edm::InputTag j_;
  double btagcut_;
  edm::InputTag m_;
  edm::InputTag e_, primaryVertices_ ;


};

JetLepInfoDumper::JetLepInfoDumper( const ParameterSet & cfg ) : 

  j_(cfg.getParameter<edm::InputTag>("j")),
  btagcut_( cfg.getParameter<double>("btagcut")),
  m_(cfg.getParameter<edm::InputTag>("m")),
  e_(cfg.getParameter<edm::InputTag>("e")),
  primaryVertices_(cfg.getParameter<InputTag>("primaryVertices"))
{
  produces<unsigned int >( "NofJets" ).setBranchAlias( "NofJets" );
  produces<unsigned int >( "NofBJets" ).setBranchAlias( "NofBJets" );
  produces<unsigned int >( "NofMuons" ).setBranchAlias( "NofMuons" );
  // produces<float >( "ThirdMuonPt" ).setBranchAlias( "ThirdMuonPt" );
  produces<unsigned int >( "NofElectrons" ).setBranchAlias( "NofElectrons" );

  produces<vector<float>  >( "MuDxyPV" ).setBranchAlias( "MuDxyPV" );
  produces<vector<float>  >( "MuDzPV" ).setBranchAlias( "MuDzPV" );
  produces<vector<float>  >( "EleDxyPV" ).setBranchAlias( "EleDxyPV" );
  produces<vector<float>  >( "EleDzPV" ).setBranchAlias( "EleDzPV" );
  
  produces< vector<float> >( "JECUnc" ).setBranchAlias( "JECUnc" );
  produces< vector<float> >( "JECVal" ).setBranchAlias( "JECVal" );
  produces< vector<float> >( "genJetPt" ).setBranchAlias( "genJetPt" );
  produces< vector<float> >( "genJetEta" ).setBranchAlias( "genJetEta" );
  produces< vector<float> >( "genJetPhi" ).setBranchAlias( "genJetPhi" );
}



void JetLepInfoDumper::produce( Event & evt, const EventSetup & es) {
  auto_ptr< unsigned int > nofJets ( new unsigned int );
  auto_ptr< unsigned int > nofBJets ( new unsigned int );
  auto_ptr< unsigned int > nofMuons ( new unsigned int );
  auto_ptr< unsigned int > nofElectrons ( new unsigned int );

  auto_ptr< vector<float>  > mudxypv ( new vector<float> );
  auto_ptr< vector<float>  > mudzpv ( new vector<float> );

  auto_ptr< vector<float>  > eledxypv ( new vector<float> );
  auto_ptr< vector<float>  > eledzpv ( new vector<float> );


  auto_ptr< vector<float> > jecunc ( new vector<float> );
  auto_ptr< vector<float> > jecval ( new vector<float> );

  auto_ptr< vector<float> > genjetpt ( new vector<float> );
  auto_ptr< vector<float> > genjeteta ( new vector<float> );
  auto_ptr< vector<float> > genjetphi ( new vector<float> );


  //jet
   * nofJets = -1;
   * nofBJets = 0;
   // *thirdJetPt = -1.;
   // *thirdJetEta = -99999.;

  Handle<pat::JetCollection > jColl;
  if (!evt.getByLabel(j_,jColl)) {
    // LogWarning("") << ">>> hjj collection does not exist !!!";
    return;
  }


  * nofJets = jColl->size();
  //jets should be ordered in pt, let's see....
  /*    
  for(pat::JetCollection::const_iterator i = jColl->begin(); i != jColl->end(); ++i) {

  std::cout << "jet.pt() =" <<  i->pt() << std::endl;
 }
  */
  
  //yes, verified!!!!
 
  /// JEC  uncertainty for jets
 
  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  es.get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl); 
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);



   unsigned int n = 0;
   * nofJets < 999 ? n = * nofJets : n = 999;  
   for (unsigned int i = 0; i < n ; i++){

     //jetsPt->push_back(((jColl->begin()) + i)->pt());
     //jetsEta->push_back(((jColl->begin()) + i)->eta());
     //jetsbtag->push_back(((jColl->begin()) + i)->bDiscriminator("combinedSecondaryVertexBJetTags"));
    // std::cout << "jet[2].pt() =" << *thirdJetPt  << std::endl;
      //  std::cout << "jet[2]->parton.pt() =" <<   << std::endl;
     if ( ((jColl->begin()) + i)->bDiscriminator("combinedSecondaryVertexBJetTags") > btagcut_ ) (*nofBJets)+=1 ;
     double eta = ((jColl->begin()) + i)->eta();
     double pt = ((jColl->begin()) + i)->pt();

     jecUnc->setJetEta(eta);
     jecUnc->setJetPt(pt); // here you must use the CORRECTED jet pt
     double unc = jecUnc->getUncertainty(true);
     jecunc->push_back( unc);

     double jecvalue  = ((jColl->begin()) + i)->jecFactor("Uncorrected");
     jecval->push_back( jecvalue);
     // std::cout << "JEC unc for jet with pt, eta, jecval " << pt << ", " <<  eta << ", "<< jecvalue <<"is unc " <<unc << std::endl;  

     const reco::GenJet * genJet =  ((jColl->begin()) + i)->genJet();
     if (genJet!=0){
       //  std::cout << "jet->genJetPt.pt() =" << ((jColl->begin()) + i)->genJet()->pt()  << std::endl;}
       genjetpt->push_back(genJet->pt());
       genjeteta->push_back(genJet->eta());
       genjetphi->push_back(genJet->phi());

     }else{
       genjetpt->push_back(-99999.);
       genjeteta->push_back(-99999.);
       genjetphi->push_back(-99999.);
     }
   }

  // cms.PSet(
//     tag = cms.untracked.string("genJetPt"),
//     quantity = cms.untracked.string("genJet().pt")
//     ),
//     cms.PSet(
//     tag = cms.untracked.string("genJetEta"),
//     quantity = cms.untracked.string("genJet().eta")
//     ),
//     cms.PSet(
//     tag = cms.untracked.string("genJetPhi"),
//     quantity = cms.untracked.string("genJet().phi")
//     ),


 // dZ wrt PV

  Handle<reco::VertexCollection> primaryVertices;  // Collection of primary Vertices
  evt.getByLabel(primaryVertices_, primaryVertices);


  //muon
  * nofMuons = -1;


  Handle<pat::MuonCollection > mColl;
  if (!evt.getByLabel(m_,mColl)) {
    // LogWarning("") << ">>> hjj collection does not exist !!!";
    return;
  }


  * nofMuons = mColl->size();
  //jets should be ordered in pt, let's see....
    
  for(pat::MuonCollection::const_iterator i = mColl->begin(); i != mColl->end(); ++i) {
    reco::TrackRef mutrack = i->track();
    if (mutrack.isNonnull()){
    mudxypv->push_back(mutrack->dxy(primaryVertices->begin()->position() ));
    mudzpv->push_back(mutrack->dz(primaryVertices->begin()->position() ));
    //  std::cout << "muon.pt() =" <<  i->pt() << std::endl;
    }else{
      mudxypv->push_back( 9999);
      mudzpv->push_back( 9999);
    }
 }
  
   
  //yes, verified!!!!
 


  //electrons
  * nofElectrons = -1;
 

  Handle<pat::ElectronCollection > eColl;
  if (!evt.getByLabel(e_,eColl)) {
    // LogWarning("") << ">>> hjj collection does not exist !!!";
    return;
  }


  * nofElectrons = eColl->size();
  //jets should be ordered in pt, let's see....
   
  for(pat::ElectronCollection::const_iterator i = eColl->begin(); i != eColl->end(); ++i) {
    reco::GsfTrackRef eletrack =  i->gsfTrack() ;
    if (eletrack.isNonnull()){
      eledxypv->push_back( eletrack->dxy(primaryVertices->begin()->position() ));
      eledzpv->push_back(eletrack->dz(primaryVertices->begin()->position() ));}
    else{
      eledxypv->push_back(9999 );
      eledzpv->push_back(9999);
    }
  }
  //std::cout << "electron.pt() =" <<  i->pt() << std::endl;


  
   
  //yes, verified!!!!
  

 



 evt.put( nofJets, "NofJets" );
 evt.put( nofBJets, "NofBJets" );
 // evt.put( jetsPt, "JetsPt" );
 // evt.put( jetsEta, "JetsEta" );

 evt.put( nofMuons, "NofMuons" );
 evt.put( nofElectrons, "NofElectrons" );

 evt.put( mudxypv, "MuDxyPV" );
 evt.put( mudzpv, "MuDzPV" );

 evt.put( eledxypv, "EleDxyPV" );
 evt.put( eledzpv, "EleDzPV" );


 evt.put( jecunc, "JECUnc" );
 evt.put( jecval, "JECVal" );

 evt.put( genjetpt, "genJetPt" );
 evt.put( genjeteta, "genJetEta" );
 evt.put( genjetphi, "genJetPhi" );
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( JetLepInfoDumper );

