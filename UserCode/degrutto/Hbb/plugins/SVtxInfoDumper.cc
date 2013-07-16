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

#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"


#include <vector>

using namespace edm;
using namespace std;
using namespace reco;


class SVtxInfoDumper : public edm::EDProducer {
public:
  SVtxInfoDumper( const edm::ParameterSet & );
   
private:
  void produce( edm::Event &, const edm::EventSetup & );
  edm::InputTag zjj_;


};

SVtxInfoDumper::SVtxInfoDumper( const ParameterSet & cfg ) : 
 zjj_(cfg.getParameter<edm::InputTag>("zjj"))
{
  produces<std::vector<int> >( "zjjJet1NofSVtx" ).setBranchAlias( "zjjJet1NofSVtx" );
  produces<std::vector<int> >( "zjjJet2NofSVtx" ).setBranchAlias( "zjjJet2NofSVtx" );


  produces<std::vector<int> >( "zjjJet1SVNofTrk" ).setBranchAlias( "zjjJet1SVNofTrk" );
  produces<std::vector<int> >( "zjjJet2SVNofTrk" ).setBranchAlias( "zjjJet2SVNofTrk" );

 produces<std::vector<float> >( "zjjJet1SVdist" ).setBranchAlias( "zjjJet1SVdist" );
 produces<std::vector<float> >( "zjjJet2SVdist" ).setBranchAlias( "zjjJet2SVdist" );

 produces<std::vector<float> >( "zjjJet1SVdistErr" ).setBranchAlias( "zjjJet1SVdistErr" );
 produces<std::vector<float> >( "zjjJet2SVdistErr" ).setBranchAlias( "zjjJet2SVdistErr" );

 produces<std::vector<float> >( "zjjJet1SVdistSig" ).setBranchAlias( "zjjJet1SVdistSig" );
 produces<std::vector<float> >( "zjjJet2SVdistSig" ).setBranchAlias( "zjjJet2SVdistSig" );

 produces<std::vector<float> >( "zjjDistSVtx" ).setBranchAlias( "zjjDistSVtx" );

 
}



void SVtxInfoDumper::produce( Event & evt, const EventSetup & ) {
  auto_ptr<vector<int> > zjjjet1NofSVtx ( new vector<int> );
  auto_ptr<vector<int> > zjjjet2NofSVtx( new vector<int> );  

  auto_ptr<vector<int> > zjjjet1SVNofTrk ( new vector<int> );
  auto_ptr<vector<int> > zjjjet2SVNofTrk( new vector<int> );  

  auto_ptr<vector<float> > zjjjet1SVdist ( new vector<float> );
  auto_ptr<vector<float> > zjjjet2SVdist( new vector<float> );  

  auto_ptr<vector<float> > zjjjet1SVdistErr ( new vector<float> );
  auto_ptr<vector<float> > zjjjet2SVdistErr( new vector<float> );  

  auto_ptr<vector<float> > zjjjet1SVdistSig ( new vector<float> );
  auto_ptr<vector<float> > zjjjet2SVdistSig( new vector<float> );  

  auto_ptr<vector<float> > zjjdistSVtx ( new vector<float> );


  // MET

  Handle<CandidateView > zjjColl;
  if (!evt.getByLabel(zjj_,zjjColl)) {
    // LogWarning("") << ">>> zjj collection does not exist !!!";
    return;
  }
for(CandidateView::const_iterator i = zjjColl->begin(); i != zjjColl->end(); ++i) {
  const Candidate * dau1 = i->daughter(0);
    const Candidate * dau2 = i->daughter(1);
    if(dau1 == 0|| dau2 == 0) 
      throw Exception(errors::InvalidReference) <<
	"one of the two daughter does not exist\n";
    //  const Candidate * c1 = dau1->masterClone().get();
    const pat::Jet * jet1 = dynamic_cast<const pat::Jet*>(dau1);

    //const Candidate * c2 = dau2->masterClone().get();
    const pat::Jet * jet2 = dynamic_cast<const pat::Jet*>(dau2);
    /*
    zjjjet1NofSVtx->push_back(  -1);  
    zjjjet1SVNofTrk->push_back(  -1);  
    zjjjet1SVdist->push_back(  -1);  
    zjjjet1SVdistErr->push_back(  -1);  
    zjjjet1SVdistSig->push_back(  -1);  

    zjjjet2NofSVtx->push_back(  -1); 
    zjjjet2SVNofTrk->push_back(  -1);   
    zjjjet2SVdist->push_back(  -1);      
    zjjjet2SVdistErr->push_back(  -1);  
    zjjjet2SVdistSig->push_back(  -1);  
    
    zjjdistSVtx->push_back(  -1);  
    */
	// store secondary vertices( directions....)....

    math::XYZVector dir1(0., 0. ,0.);
    math::XYZVector dir2(0., 0. ,0.);


    const reco::SecondaryVertexTagInfo *svTagInfo1 =
      jet1->tagInfoSecondaryVertex("secondaryVertex");

    
    if( svTagInfo1 != 0  ){
      zjjjet1NofSVtx -> push_back(svTagInfo1->nVertices()) ;

	if (svTagInfo1->nVertices() < 1 )
			continue;

		const reco::Vertex &sv1 = svTagInfo1->secondaryVertex(0);


		// discard cases in which the sv is fake...
		if (sv1.isFake()) continue;
		

		zjjjet1SVNofTrk->push_back( sv1.tracksSize());
		
		// the precomputed transverse distance to the primary vertex

	 
		Measurement1D distance1 = svTagInfo1->flightDistance(0, true);
	 
		zjjjet1SVdist->push_back( distance1.value());

		zjjjet1SVdistErr->push_back( distance1.error());

		zjjjet1SVdistErr->push_back( distance1.significance());

		// the precomputed direction with respect to the primary vertex
		GlobalVector dir = svTagInfo1->flightDirection(0);

		// unfortunately CMSSW hsa all kinds of vectors,
		// and sometimes we need to convert them *sigh*
			dir1 = math::XYZVector(dir.x(), dir.y(), dir.z());


    }
	

    const reco::SecondaryVertexTagInfo *svTagInfo2 =
      jet2->tagInfoSecondaryVertex("secondaryVertex");
    if( svTagInfo2 != 0  ){
      zjjjet2NofSVtx ->push_back( svTagInfo2->nVertices() );

	if (svTagInfo2->nVertices() < 1 )
			continue;
		const reco::Vertex &sv2 = svTagInfo2->secondaryVertex(0);
	
		// discard cases in which the sv is fake...
		if (sv2.isFake()) continue;
		
		zjjjet2SVNofTrk ->push_back( sv2.tracksSize());
		
		// the precomputed transverse distance to the primary vertex
		Measurement1D distance2 = svTagInfo2->flightDistance(0, true);

		zjjjet2SVdist->push_back( distance2.value());

		zjjjet2SVdistErr->push_back( distance2.error());

		zjjjet2SVdistErr->push_back(distance2.significance());

		// the precomputed direction with respect to the primary vertex
		GlobalVector dir = svTagInfo2->flightDirection(0);

		// unfortunately CMSSW hsa all kinds of vectors,
		// and sometimes we need to convert them *sigh*
		dir2 =  math::XYZVector(dir.x(), dir.y(), dir.z());
		//v_sv.push_back(dir2);
    }



    //unsigned int n_sv = 	v_sv.size();
	//std::cout << "n of secondary vertices associated to different jets in the event == " << n_sv << std::endl;
	//if (n_sv>1){
	 	  double dist_twoVertices = ( dir1 - dir2 ).R();
		  // std::cout << "diatance beewtween two secondary vertices associated to different jets in the event == " << dist_twoVertices << std::endl;
		  zjjdistSVtx->push_back(  dist_twoVertices);


    
 

		  /*		   std::cout << "before sv instance" << std::endl;    

	 std::cout << "zjjJet1NofSVtx " <<zjjJet1NofSVtx[i] << std::endl;    
	 std::cout << "zjjJet2NofSVtx " <<zjjJet2NofSVtx[i] << std::endl;    

	 std::cout << "zjjjet1SVNofTrk" <<zzjjjet1SVNofTrk[i] << std::endl;    
	 std::cout << "zjjjet2SVNofTrk" <<zzjjjet2SVNofTrk[i] << std::endl;    

	 std::cout << "zjjjet1SVdist" <<zzjjjet1SVdist[i] << std::endl;    
	 std::cout << "zjjjet2SVdist" <<zzjjjet2SVdist[i] << std::endl;    

	 std::cout << "zzjjDistSVtx" <<zjjDistSVtx[i] << std::endl;    
		  */  
 }


 evt.put( zjjjet1NofSVtx, "zjjJet1NofSVtx" );
 evt.put( zjjjet2NofSVtx, "zjjJet2NofSVtx" );

 evt.put( zjjjet1SVNofTrk, "zjjJet1SVNofTrk" );
 evt.put( zjjjet2SVNofTrk, "zjjJet2SVNofTrk" );

 evt.put( zjjjet1SVdist, "zjjJet1SVdist" );
 evt.put( zjjjet2SVdist, "zjjJet2SVdist" );

 evt.put( zjjjet1SVdistErr, "zjjJet1SVdistErr" );
 evt.put( zjjjet2SVdistErr, "zjjJet2SVdistErr" );

 evt.put( zjjjet1SVdistSig, "zjjJet1SVdistSig" );
 evt.put( zjjjet2SVdistSig, "zjjJet2SVdistSig" );



 evt.put( zjjdistSVtx, "zjjDistSVtx" );
 //std::cout << "after sv instance" << std::endl;
 
 
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( SVtxInfoDumper );

