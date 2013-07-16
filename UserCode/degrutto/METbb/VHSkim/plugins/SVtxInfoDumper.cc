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
  edm::InputTag hjj_;


};

SVtxInfoDumper::SVtxInfoDumper( const ParameterSet & cfg ) : 
 hjj_(cfg.getParameter<edm::InputTag>("hjj"))
{
  produces<std::vector<int> >( "hjjJet1NofSVtx" ).setBranchAlias( "hjjJet1NofSVtx" );
  produces<std::vector<int> >( "hjjJet2NofSVtx" ).setBranchAlias( "hjjJet2NofSVtx" );


  produces<std::vector<int> >( "hjjJet1SVNofTrk" ).setBranchAlias( "hjjJet1SVNofTrk" );
  produces<std::vector<int> >( "hjjJet2SVNofTrk" ).setBranchAlias( "hjjJet2SVNofTrk" );

 produces<std::vector<float> >( "hjjJet1SVdist" ).setBranchAlias( "hjjJet1SVdist" );
 produces<std::vector<float> >( "hjjJet2SVdist" ).setBranchAlias( "hjjJet2SVdist" );

 produces<std::vector<float> >( "hjjJet1SVdistErr" ).setBranchAlias( "hjjJet1SVdistErr" );
 produces<std::vector<float> >( "hjjJet2SVdistErr" ).setBranchAlias( "hjjJet2SVdistErr" );

 produces<std::vector<float> >( "hjjJet1SVdistSig" ).setBranchAlias( "hjjJet1SVdistSig" );
 produces<std::vector<float> >( "hjjJet2SVdistSig" ).setBranchAlias( "hjjJet2SVdistSig" );

 produces<std::vector<float> >( "hjjDistSVtx" ).setBranchAlias( "hjjDistSVtx" );

 
}



void SVtxInfoDumper::produce( Event & evt, const EventSetup & ) {
  auto_ptr<vector<int> > hjjjet1NofSVtx ( new vector<int> );
  auto_ptr<vector<int> > hjjjet2NofSVtx( new vector<int> );  

  auto_ptr<vector<int> > hjjjet1SVNofTrk ( new vector<int> );
  auto_ptr<vector<int> > hjjjet2SVNofTrk( new vector<int> );  

  auto_ptr<vector<float> > hjjjet1SVdist ( new vector<float> );
  auto_ptr<vector<float> > hjjjet2SVdist( new vector<float> );  

  auto_ptr<vector<float> > hjjjet1SVdistErr ( new vector<float> );
  auto_ptr<vector<float> > hjjjet2SVdistErr( new vector<float> );  

  auto_ptr<vector<float> > hjjjet1SVdistSig ( new vector<float> );
  auto_ptr<vector<float> > hjjjet2SVdistSig( new vector<float> );  

  auto_ptr<vector<float> > hjjdistSVtx ( new vector<float> );


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
    const pat::Jet * jet1 = dynamic_cast<const pat::Jet*>(dau1);

    //const Candidate * c2 = dau2->masterClone().get();
    const pat::Jet * jet2 = dynamic_cast<const pat::Jet*>(dau2);
    /*
    hjjjet1NofSVtx->push_back(  -1);  
    hjjjet1SVNofTrk->push_back(  -1);  
    hjjjet1SVdist->push_back(  -1);  
    hjjjet1SVdistErr->push_back(  -1);  
    hjjjet1SVdistSig->push_back(  -1);  

    hjjjet2NofSVtx->push_back(  -1); 
    hjjjet2SVNofTrk->push_back(  -1);   
    hjjjet2SVdist->push_back(  -1);      
    hjjjet2SVdistErr->push_back(  -1);  
    hjjjet2SVdistSig->push_back(  -1);  
    
    hjjdistSVtx->push_back(  -1);  
    */
	// store secondary vertices( directions....)....

    math::XYZVector dir1(0., 0. ,0.);
    math::XYZVector dir2(0., 0. ,0.);


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



    //unsigned int n_sv = 	v_sv.size();
	//std::cout << "n of secondary vertices associated to different jets in the event == " << n_sv << std::endl;
	//if (n_sv>1){
	 	  double dist_twoVertices = ( dir1 - dir2 ).R();
		  // std::cout << "diatance beewtween two secondary vertices associated to different jets in the event == " << dist_twoVertices << std::endl;
		  hjjdistSVtx->push_back(  dist_twoVertices);


    
 

		  /*		   std::cout << "before sv instance" << std::endl;    

	 std::cout << "hjjJet1NofSVtx " <<hjjJet1NofSVtx[i] << std::endl;    
	 std::cout << "hjjJet2NofSVtx " <<hjjJet2NofSVtx[i] << std::endl;    

	 std::cout << "hjjjet1SVNofTrk" <<zhjjjet1SVNofTrk[i] << std::endl;    
	 std::cout << "hjjjet2SVNofTrk" <<zhjjjet2SVNofTrk[i] << std::endl;    

	 std::cout << "hjjjet1SVdist" <<zhjjjet1SVdist[i] << std::endl;    
	 std::cout << "hjjjet2SVdist" <<zhjjjet2SVdist[i] << std::endl;    

	 std::cout << "zhjjDistSVtx" <<hjjDistSVtx[i] << std::endl;    
		  */  
 }


 evt.put( hjjjet1NofSVtx, "hjjJet1NofSVtx" );
 evt.put( hjjjet2NofSVtx, "hjjJet2NofSVtx" );

 evt.put( hjjjet1SVNofTrk, "hjjJet1SVNofTrk" );
 evt.put( hjjjet2SVNofTrk, "hjjJet2SVNofTrk" );

 evt.put( hjjjet1SVdist, "hjjJet1SVdist" );
 evt.put( hjjjet2SVdist, "hjjJet2SVdist" );

 evt.put( hjjjet1SVdistErr, "hjjJet1SVdistErr" );
 evt.put( hjjjet2SVdistErr, "hjjJet2SVdistErr" );

 evt.put( hjjjet1SVdistSig, "hjjJet1SVdistSig" );
 evt.put( hjjjet2SVdistSig, "hjjJet2SVdistSig" );



 evt.put( hjjdistSVtx, "hjjDistSVtx" );
 //std::cout << "after sv instance" << std::endl;
 
 
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( SVtxInfoDumper );

