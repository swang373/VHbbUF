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
#include <TVector2.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include "DataFormats/Math/interface/deltaPhi.h"
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

#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"

#include "PhysicsTools/CandUtils/interface/CenterOfMassBooster.h"
#include "PhysicsTools/CandUtils/interface/Booster.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <vector>

using namespace edm;
using namespace std;
using namespace reco;


class SVtxInfoDumper : public edm::EDProducer {
public:
  SVtxInfoDumper( const edm::ParameterSet & );
   
private:
  void produce( edm::Event &, const edm::EventSetup & );
  double getDeltaTheta(const pat::Jet* j1, const pat::Jet* j2);
  TVector2  getTvect( const pat::Jet* patJet );
  //  double getHelicity( const pat::Jet* jet , TVector3 boost );
  double getNofChTrks( const pat::Jet* patJet );
 edm::InputTag hjj_;


};

SVtxInfoDumper::SVtxInfoDumper( const ParameterSet & cfg ) : 
 hjj_(cfg.getParameter<edm::InputTag>("hjj"))
{
  produces<std::vector<int> >( "hjjJet1NofSVtx" ).setBranchAlias( "hjjJet1NofSVtx" );
  produces<std::vector<int> >( "hjjJet2NofSVtx" ).setBranchAlias( "hjjJet2NofSVtx" );


  produces<std::vector<int> >( "hjjJet1SVNofTrk" ).setBranchAlias( "hjjJet1SVNofTrk" );
  produces<std::vector<int> >( "hjjJet2SVNofTrk" ).setBranchAlias( "hjjJet2SVNofTrk" );

  produces<std::vector<int> >( "hjjJet1NofChTrk" ).setBranchAlias( "hjjJet1NofChTrk" );
  produces<std::vector<int> >( "hjjJet2NofChTrk" ).setBranchAlias( "hjjJet2NofChTrk" );

 produces<std::vector<float> >( "hjjJet1SVdist" ).setBranchAlias( "hjjJet1SVdist" );
 produces<std::vector<float> >( "hjjJet2SVdist" ).setBranchAlias( "hjjJet2SVdist" );

 produces<std::vector<float> >( "hjjJet1SVdistErr" ).setBranchAlias( "hjjJet1SVdistErr" );
 produces<std::vector<float> >( "hjjJet2SVdistErr" ).setBranchAlias( "hjjJet2SVdistErr" );

 produces<std::vector<float> >( "hjjJet1SVdistSig" ).setBranchAlias( "hjjJet1SVdistSig" );
 produces<std::vector<float> >( "hjjJet2SVdistSig" ).setBranchAlias( "hjjJet2SVdistSig" );

 produces<std::vector<float> >( "hjjDistSVtx" ).setBranchAlias( "hjjDistSVtx" );

 produces<std::vector<float> >( "hjjJet1DeltaTheta" ).setBranchAlias( "hjjJet1DeltaTheta" );
 produces<std::vector<float> >( "hjjJet1DeltaThetaAK7" ).setBranchAlias( "hjjJet1DeltaThetaAK7" );

 produces<std::vector<float> >( "hjjJet1higgsHelicity" ).setBranchAlias( "hjjJet1higgsHelicity" );
 produces<std::vector<float> >( "hjjJet2higgsHelicity" ).setBranchAlias( "hjjJet2higgsHelicity" );




 
}

double SVtxInfoDumper::getNofChTrks( const pat::Jet* patJet ){

 unsigned int nOfconst = 0;

//re-reconstruct the jet direction with the charged tracks
  std::vector<reco::PFCandidatePtr>
    patJetpfc = patJet->getPFConstituents();
  for(size_t idx = 0; idx < patJetpfc.size(); idx++){
    if( patJetpfc.at(idx)->charge() != 0 ){
      nOfconst++;
    }
  }

  return nOfconst;
}


TVector2 SVtxInfoDumper::getTvect( const pat::Jet* patJet ){

  TVector2 t_Vect(0,0);
  TVector2 null(0,0);
  TVector2 ci(0,0);
  TLorentzVector pi(0,0,0,0);
  TLorentzVector J(0,0,0,0);
  TVector2 r(0,0);
  double patJetpfcPt = 1e10;
  double r_mag = 1e10;
  unsigned int nOfconst = 0;

//re-reconstruct the jet direction with the charged tracks
  std::vector<reco::PFCandidatePtr>
    patJetpfc = patJet->getPFConstituents();
  for(size_t idx = 0; idx < patJetpfc.size(); idx++){
    if( patJetpfc.at(idx)->charge() != 0 ){
      pi.SetPtEtaPhiE( patJetpfc.at(idx)->pt(),
patJetpfc.at(idx)->eta(), patJetpfc.at(idx)->phi(),
patJetpfc.at(idx)->energy() );
      J += pi;
      nOfconst++;
    }
  }

// if there are less than two charged tracks do not calculate the pull
//(there is not enough info). It returns a null vector

  if( nOfconst < 2 )
    return null;

  TVector2 v_J( J.Rapidity(), J.Phi() );
//calculate TVector using only charged tracks
  for(size_t idx = 0; idx < patJetpfc.size(); idx++){
    if( patJetpfc.at(idx)->charge() != 0  ){
      patJetpfcPt = patJetpfc.at(idx)->pt();
      pi.SetPtEtaPhiE( patJetpfc.at(idx)->pt(),
patJetpfc.at(idx)->eta(), patJetpfc.at(idx)->phi(),
patJetpfc.at(idx)->energy() );
      r.Set( pi.Rapidity() - J.Rapidity(), deltaPhi(
patJetpfc.at(idx)->phi(), J.Phi() ) );
      r_mag = r.Mod();
      t_Vect += ( patJetpfcPt / J.Pt() ) * r_mag * r;
    }
  }

  return t_Vect;

}





double SVtxInfoDumper::getDeltaTheta( const pat::Jet* j1, const pat::Jet* j2 ){

  double deltaTheta = 1e10;
  TLorentzVector pi(0,0,0,0);
  TLorentzVector v_j1(0,0,0,0);
  TLorentzVector v_j2(0,0,0,0);

//re-reconstruct the jet direction with the charged tracks

  std::vector<reco::PFCandidatePtr>
    j1pfc = j1->getPFConstituents();
  for(size_t idx = 0; idx < j1pfc.size(); idx++){
    if( j1pfc.at(idx)->charge() != 0 ){
      pi.SetPtEtaPhiE( j1pfc.at(idx)->pt(), j1pfc.at(idx)->eta(),
j1pfc.at(idx)->phi(), j1pfc.at(idx)->energy() );
      v_j1 += pi;
    }
  }

//re-reconstruct the jet direction with the charged tracks

  std::vector<reco::PFCandidatePtr>
    j2pfc = j2->getPFConstituents();
  for(size_t idx = 0; idx < j2pfc.size(); idx++){
    if( j2pfc.at(idx)->charge() != 0 ){
      pi.SetPtEtaPhiE( j2pfc.at(idx)->pt(), j2pfc.at(idx)->eta(),
j2pfc.at(idx)->phi(), j2pfc.at(idx)->energy() );
      v_j2 += pi;
    }
  }

  if( v_j2.Mag() <=0 or v_j1.Mag() <=0 )
    return deltaTheta = 1e10;

  TVector2 v2_j1( v_j1.Rapidity(), v_j1.Phi());
  TVector2 v2_j2( v_j2.Rapidity(), v_j2.Phi());
//use j1 to calculate the pull vector

  TVector2 t = getTvect(j1);

  if( t.Mod() == 0 )
    return deltaTheta = 1e10;

  Double_t deltaphi = deltaPhi( v_j2.Phi(), v_j1.Phi() );
  Double_t deltaeta = v_j2.Rapidity() - v_j1.Rapidity();
  TVector2 BBdir( deltaeta, deltaphi );

  deltaTheta = t.DeltaPhi(BBdir);

  return deltaTheta;

}

/*
double SVtxInfoDumper::getHelicity( const pat::Jet* jet , TVector3 boost ){
  double hel = 1e10;
  TLorentzVector j;
  j.SetPtEtaPhiE( jet->pt(), jet->eta(), jet->phi(), jet->energy() );
  j.Boost( -boost );
  hel = TMath::Cos( j.Vect().Angle( boost ) );
  return hel;
}
*/


void SVtxInfoDumper::produce( Event & evt, const EventSetup & ) {
  auto_ptr<vector<int> > hjjjet1NofSVtx ( new vector<int> );
  auto_ptr<vector<int> > hjjjet2NofSVtx( new vector<int> );  

  auto_ptr<vector<int> > hjjjet1SVNofTrk ( new vector<int> );
  auto_ptr<vector<int> > hjjjet2SVNofTrk( new vector<int> );  

  auto_ptr<vector<int> > hjjjet1NofChTrk ( new vector<int> );
  auto_ptr<vector<int> > hjjjet2NofChTrk( new vector<int> );  

  auto_ptr<vector<float> > hjjjet1SVdist ( new vector<float> );
  auto_ptr<vector<float> > hjjjet2SVdist( new vector<float> );  

  auto_ptr<vector<float> > hjjjet1SVdistErr ( new vector<float> );
  auto_ptr<vector<float> > hjjjet2SVdistErr( new vector<float> );  

  auto_ptr<vector<float> > hjjjet1SVdistSig ( new vector<float> );
  auto_ptr<vector<float> > hjjjet2SVdistSig( new vector<float> );  

  auto_ptr<vector<float> > hjjdistSVtx ( new vector<float> );

  auto_ptr<vector<float> > hjjjet1DeltaTheta  ( new vector<float> );
 auto_ptr<vector<float> > hjjjet1DeltaThetaAK7  ( new vector<float> );

  auto_ptr<vector<float> > hjjjet1higgsHelicity  ( new vector<float> );
  auto_ptr<vector<float> > hjjjet2higgsHelicity ( new vector<float> );

 


  // MET

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
    //  const Candidate * c1 = dau1->masterClone().get();
    const pat::Jet* jet1 = dynamic_cast<const pat::Jet*>(dau1);

    //const Candidate * c2 = dau2->masterClone().get();
    const pat::Jet* jet2 = dynamic_cast<const pat::Jet*>(dau2);


    // match an ak7pfjet for eacj jet with deltaR<0.5...
    std::vector< pat::Jet >  v_ak7PFJets; 

    edm::Handle<reco::PFJetCollection> ak7pfJets;
     if (!evt.getByLabel("ak7PFJets", ak7pfJets)) {
       std::cout << "ak7PFJets not found in the event" << std::endl;
    return;
  }

     double DeltaR =0;
     double minDeltaR =0;
     int pos = -1;  
     for (unsigned int p = 0; p < ak7pfJets->size() ; p++){
       //pat::Jet* jet7 = ak7pfJets.at[p] ;
	     const  pat::Jet  ak7jet = *((ak7pfJets->begin()) + p);
             DeltaR = deltaR(jet1->eta(), jet1->phi(),ak7jet.eta(), ak7jet.phi() );
             if ( (DeltaR < 0.5) || (DeltaR<minDeltaR )) {
              minDeltaR = DeltaR;
              pos = (int) p;
	     }
          
     }

     if (pos>-1 ) v_ak7PFJets.push_back( *((ak7pfJets->begin()) + pos) ); 
     DeltaR =0;
     minDeltaR =0;
     pos = -1;  

   
     for (unsigned int pp = 0; pp < ak7pfJets->size() ; pp++){
       //pat::Jet* jet7 = ak7pfJets.at[p] ;
       const  pat::Jet  ak7jet = *((ak7pfJets->begin()) + pp);
       DeltaR = deltaR(jet2->eta(), jet2->phi(),ak7jet.eta(), ak7jet.phi() );
       if ( (DeltaR < 0.5) || (DeltaR<minDeltaR )) {
	 minDeltaR = DeltaR;
	 pos = (int) pp;
	     }
          
     }
     if (pos>-1 ) v_ak7PFJets.push_back( *((ak7pfJets->begin()) + pos) );



    hjjjet1NofChTrk-> push_back( getNofChTrks ( jet1));
    hjjjet2NofChTrk-> push_back( getNofChTrks ( jet2));



    // now adding pull and costheta* for color studies, from jet1 and jet2  
    hjjjet1DeltaTheta->push_back( TMath::Abs( getDeltaTheta( jet1 ,jet2 ) ) ); 
    // pull for ak7 if the two jets don;t match to the say one...
    if ( v_ak7PFJets.size()==2 ){
      if (v_ak7PFJets[0].pt() != v_ak7PFJets[1].pt() ) {
   hjjjet1DeltaThetaAK7->push_back( TMath::Abs( getDeltaTheta( &v_ak7PFJets[0] ,&v_ak7PFJets[1] ) ) ); 
      }
    }


    math::XYZVector higgsBoostXYZ = i->p4().BoostToCM();
    TVector3 higgsBoost =TVector3( higgsBoostXYZ.X() ,higgsBoostXYZ.Y(), higgsBoostXYZ.Z() ); 
    
    
    //  hjjjet1higgsHelicity->push_back( getHelicity( jet1, higgsBoost ));
    // hjjjet2higgsHelicity->push_back( getHelicity( jet2, higgsBoost ));

    Booster hFrameBoost( i->boostToCM() );
    Candidate * boostedJ1_HFrame = dau1->clone();
    Candidate * boostedJ2_HFrame = dau2->clone();

    hFrameBoost.set( *boostedJ1_HFrame );
    hFrameBoost.set( *boostedJ2_HFrame );

    hjjjet1higgsHelicity->push_back( TMath::Cos( ROOT::Math::VectorUtil::Angle(  i->boostToCM(), boostedJ1_HFrame->momentum())) );
    hjjjet2higgsHelicity->push_back( TMath::Cos( ROOT::Math::VectorUtil::Angle(  i->boostToCM(), boostedJ2_HFrame->momentum())) );

    
  
    math::XYZVector dir1(0., 0. ,0.);
    math::XYZVector dir2(0., 0. ,0.);

    /*  
    const reco::SecondaryVertexTagInfo *svTagInfo1 =
      jet1->tagInfoSecondaryVertex("secondaryVertex");


    if( svTagInfo1 != 0  ){
      hjjjet1NofSVtx -> push_back(svTagInfo1->nVertices()) ;

      if (svTagInfo1->nVertices() < 1 )
	continue;

      const reco::Vertex &sv1 = svTagInfo1->secondaryVertex(0);
      
      
      // discard cases in which the sv is fake...
      if (sv1.isFake()) continue;
      
      
      hjjjet1SVNofTrk->push_back( sv1.tracksSize());
      
      // the precomputed transverse distance to the primary vertex
      
      
      Measurement1D distance1 = svTagInfo1->flightDistance(0, true);
      
      hjjjet1SVdist->push_back( distance1.value());
      
      hjjjet1SVdistErr->push_back( distance1.error());
      
      hjjjet1SVdistErr->push_back( distance1.significance());
      
      // the precomputed direction with respect to the primary vertex
      GlobalVector dir = svTagInfo1->flightDirection(0);
      
      // unfortunately CMSSW hsa all kinds of vectors,
      // and sometimes we need to convert them *sigh*
      dir1 = math::XYZVector(dir.x(), dir.y(), dir.z());
      
      
    }
    
    
    const reco::SecondaryVertexTagInfo *svTagInfo2 =
      jet2->tagInfoSecondaryVertex("secondaryVertex");
    if( svTagInfo2 != 0  ){
      hjjjet2NofSVtx ->push_back( svTagInfo2->nVertices() );
      
      if (svTagInfo2->nVertices() < 1 )
	continue;
      const reco::Vertex &sv2 = svTagInfo2->secondaryVertex(0);
      
      // discard cases in which the sv is fake...
      if (sv2.isFake()) continue;
      
      hjjjet2SVNofTrk ->push_back( sv2.tracksSize());
      
      // the precomputed transverse distance to the primary vertex
      Measurement1D distance2 = svTagInfo2->flightDistance(0, true);
      
      hjjjet2SVdist->push_back( distance2.value());
      
      hjjjet2SVdistErr->push_back( distance2.error());
      
      hjjjet2SVdistErr->push_back(distance2.significance());

      // the precomputed direction with respect to the primary vertex
      GlobalVector dir = svTagInfo2->flightDirection(0);
      
      // unfortunately CMSSW hsa all kinds of vectors,
      // and sometimes we need to convert them *sigh*
      dir2 =  math::XYZVector(dir.x(), dir.y(), dir.z());
      //v_sv.push_back(dir2);
    }
*/
    

   
    double dist_twoVertices = ( dir1 - dir2 ).R();
    // std::cout << "diatance beewtween two secondary vertices associated to different jets in the event == " << dist_twoVertices << std::endl;
    hjjdistSVtx->push_back(  dist_twoVertices);
    

    

	
 }
    
/*
 evt.put( hjjjet1NofSVtx, "hjjJet1NofSVtx" );
 evt.put( hjjjet2NofSVtx, "hjjJet2NofSVtx" );

 evt.put( hjjjet1SVNofTrk, "hjjJet1SVNofTrk" );
 evt.put( hjjjet2SVNofTrk, "hjjJet2SVNofTrk" );

 evt.put( hjjjet1NofChTrk, "hjjJet1NofChTrk" );
 evt.put( hjjjet2NofChTrk, "hjjJet2NofChTrk" );

 evt.put( hjjjet1SVdist, "hjjJet1SVdist" );
 evt.put( hjjjet2SVdist, "hjjJet2SVdist" );

 evt.put( hjjjet1SVdistErr, "hjjJet1SVdistErr" );
 evt.put( hjjjet2SVdistErr, "hjjJet2SVdistErr" );

 evt.put( hjjjet1SVdistSig, "hjjJet1SVdistSig" );
 evt.put( hjjjet2SVdistSig, "hjjJet2SVdistSig" );



 evt.put( hjjdistSVtx, "hjjDistSVtx" );
 
*/
 evt.put( hjjjet1DeltaTheta, "hjjJet1DeltaTheta" );
 evt.put( hjjjet1DeltaThetaAK7, "hjjJet1DeltaThetaAK7" );
 evt.put( hjjjet1higgsHelicity, "hjjJet1higgsHelicity" );
 evt.put( hjjjet2higgsHelicity, "hjjJet2higgsHelicity" );

 //std::cout << "after sv instance" << std::endl;
 
 
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( SVtxInfoDumper );

