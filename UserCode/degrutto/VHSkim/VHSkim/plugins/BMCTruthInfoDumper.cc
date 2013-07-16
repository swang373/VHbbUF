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

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

using namespace edm;
using namespace std;
using namespace reco;


class BMCTruthDumper : public edm::EDProducer {
public:
  BMCTruthDumper( const edm::ParameterSet & );
   
private:
  void produce( edm::Event &, const edm::EventSetup & );
  const float getDistanceForPairs (const reco::GenParticleRefVector &quark,	const reco::GenParticleRefVector &antiquark);
  double  ptPairCut_;



};

BMCTruthDumper::BMCTruthDumper( const ParameterSet & cfg ) : 

  ptPairCut_( cfg.getUntrackedParameter<double>("ptPairCut"))


{
  produces<unsigned int>("alpgenFlavorFlagBB").setBranchAlias("alpgenFlavorFlagBB");
  produces<unsigned int>("alpgenFlavorFlagCC").setBranchAlias("alpgenFlavorFlagCC");
  produces<float>("alpgenDrcc").setBranchAlias("alpgendrcc");
  produces<float>("alpgenDrbb").setBranchAlias("alpgendrbb");
  produces<float>("alpgenVpt").setBranchAlias("alpgenVpt");
  produces<vector<float> >("alpgenBpt").setBranchAlias("alpgenBpt");
  produces<vector<float> >("alpgenBeta").setBranchAlias("alpgenBeta");
  produces<vector<float> >("alpgenBphi").setBranchAlias("alpgenBeta");
  produces<unsigned int  > ("PUintimeSize").setBranchAlias( "PUintimeSize"); 
  produces<unsigned int  > ("PUouttime1plusSize").setBranchAlias( "PUouttime1plusSize"); 
  produces<unsigned int  > ("PUouttime1minusSize").setBranchAlias( "PUouttime1minusSize");
  produces<unsigned int  > ("PUouttime2plusSize").setBranchAlias( "PUouttime2plusSize");
  produces<unsigned int  > ("PUouttime2minusSize").setBranchAlias( "PUouttime2minusSize");
  produces<unsigned int  > ("PUouttime3plusSize").setBranchAlias( "PUouttime3plusSize");
  produces<unsigned int  > ("PUouttime3minusSize").setBranchAlias( "PUouttime3plusSize");


}


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
const float BMCTruthDumper::getDistanceForPairs (const reco::GenParticleRefVector &quark, const reco::GenParticleRefVector &antiquark) 
// Gets the maximum distance for quark/antiquark pairs (matched with the minimum
// distance.
{
  float drmax=-1; 

  // We loop over the quark list and get the minimum distance to an antiquark:

 
  for (reco::GenParticleRefVector::const_iterator xq = quark.begin();
       xq!=quark.end();++xq) {
    float drmin=-1;

    bool valid=true;

    for (reco::GenParticleRefVector::const_iterator xaq = antiquark.begin();
	 xaq!=antiquark.end();++xaq) {

       float dr = deltaR2 ((*xq)->eta(),(*xq)->phi(),
			  (*xaq)->eta(),(*xaq)->phi());  

      if (dr<drmin || drmin<0) {
	valid=true;
	if ((*xq)->pt()<ptPairCut_ || (*xaq)->pt()<ptPairCut_) valid=false;
	drmin=dr;
      }
    }
    
    // Once that for each quark we have the minimum... we consider it for the
    // maximum if one of the two passes the cut on pt
    if (valid && drmin>drmax) drmax=drmin;
  }

  if (drmax<0) return drmax;
  return sqrt(drmax);   // We return the actual DR value.
}



void BMCTruthDumper::produce( Event & evt, const EventSetup & ) {
  auto_ptr< unsigned int > alpgenflavorflagbb ( new unsigned int );  
  auto_ptr< unsigned int > alpgenflavorflagcc ( new unsigned int );  
  auto_ptr< float > alpgendrbb( new float );  
  auto_ptr< float > alpgendrcc( new float );  
  auto_ptr< float > alpgenvpt( new float );  
  auto_ptr< vector<float> > alpgenbpt( new vector<float> );  
  auto_ptr< vector<float> > alpgenbphi (new vector<float> );  
  auto_ptr< vector<float> > alpgenbeta( new vector<float> );  

  auto_ptr<unsigned int> pUintimeSize (new unsigned int); 
  auto_ptr<unsigned int> pUouttime1plusSize (new unsigned int); 
  auto_ptr<unsigned int> pUouttime1minusSize (new unsigned int); 
  auto_ptr<unsigned int> pUouttime2plusSize (new unsigned int); 
  auto_ptr<unsigned int> pUouttime2minusSize  (new unsigned int); 
  auto_ptr<unsigned int> pUouttime3plusSize (new unsigned int); 
  auto_ptr<unsigned int> pUouttime3minusSize  (new unsigned int); 

  reco::GenParticleRefVector charm;
  reco::GenParticleRefVector anticharm;
  reco::GenParticleRefVector bottom;
  reco::GenParticleRefVector antibottom;

  float ptmax_b=-1;
 
  float ymin_b=-1;

 

 // Variables for storing in the Event
  * alpgenflavorflagbb = 0;
  * alpgenflavorflagcc = 0;
  * alpgendrbb = 0;
  * alpgendrcc = 0;
  * alpgenvpt =0;
 

 

  // For production in ALPGEN (status code==3)
  
  int ncharms=0;
  int nbottoms=0;


  // We process the information

  //cout<<"-------------------------------------------------------"<<endl;

  edm::Handle<reco::GenParticleCollection> genpartH;
  evt.getByLabel("genParticles",genpartH);

  //  GenParticleCollection genpart = *genpartH;
  for (reco::GenParticleCollection::const_iterator xpart = genpartH->begin();
       xpart!=genpartH->end();++xpart) {

     // Z/W
    if (xpart->status()==3) {
      int hepid = abs(xpart->pdgId());

      if (hepid==23 || hepid==22 ) {

	//cout<<"CHARM QUARK: "<<hepid<<" "<<xpart->mass()<<" "<<xpart->pt()<<" "<<xpart->rapidity()<<endl;
	*alpgenvpt = xpart->pt() ;  // Only massive quarks
	   }
    }


    // For charm and bottom with status==3 we do some counting to try to
    // discover which class of event it is.
    if (xpart->status()==3) {
      int hepid = abs(xpart->pdgId());

      if (hepid==4) {

	//cout<<"CHARM QUARK: "<<hepid<<" "<<xpart->mass()<<" "<<xpart->pt()<<" "<<xpart->rapidity()<<endl;
	if (xpart->mass()>1.0) ++ncharms;  // Only massive quarks
      }
      else if (hepid==5) {

	//cout<<"BOTTOM QUARK: "<<xpart->pdgId()<<" "<<xpart->mass()<<endl;
	if (xpart->mass()>3.0) ++nbottoms;
      }
    }



    // We need to identify the charm and the bottom at parton level. The simplest
    // is to look for the strings/clusters and scan their mothers.

    if (xpart->pdgId()==92 || xpart->pdgId()==91) {
      reco::GenParticleRefVector mothers = xpart->motherRefVector();
      for (reco::GenParticleRefVector::const_iterator xmot = mothers.begin();
	   xmot!=mothers.end();++xmot) {

	// We store the found b's and c's:
	if ((*xmot)->pdgId()==4) {
	  //cout<<"charm: "<<(*xmot)->pt()<<" "<<(*xmot)->rapidity()<<endl;
	  charm.push_back(*xmot);
	  //if ((*xmot)->pt()>(*pt_singlec)) (*pt_singlec)=(*xmot)->pt();
	  //if (ymin_c<0 || fabs((*xmot)->rapidity())<ymin_c) ymin_c=fabs((*xmot)->rapidity());
	}
	else if ((*xmot)->pdgId()==-4) {
	  //cout<<"anticharm: "<<(*xmot)->pt()<<" "<<(*xmot)->rapidity()<<endl;
	  anticharm.push_back(*xmot);
          //if ((*xmot)->pt()>(*pt_singlec)) (*pt_singlec)=(*xmot)->pt();
          //if (ymin_c<0 || fabs((*xmot)->rapidity())<ymin_c) ymin_c=fabs((*xmot)->rapidity());
	}
	else if ((*xmot)->pdgId()==5) {
	  bottom.push_back(*xmot);
	   if ((*xmot)->pt()>ptmax_b) ptmax_b=(*xmot)->pt();
          if (ymin_b<0 || fabs((*xmot)->rapidity())<ymin_b) ymin_b=fabs((*xmot)->rapidity());
	}
        else if ((*xmot)->pdgId()==-5) {
	  antibottom.push_back(*xmot);
	  if ((*xmot)->pt()>ptmax_b) ptmax_b=(*xmot)->pt();
          if (ymin_b<0 || fabs((*xmot)->rapidity())<ymin_b) ymin_b=fabs((*xmot)->rapidity());
	}
      }
    }
  }

  // check to identify HF samples
  if (nbottoms>2 || ncharms>2) {
    //    cerr<<"HEPTF-ERROR: Problems identifying the type of sample... matching may be wrong!"<<endl;
  }

  // Filling control histograms:

  if (nbottoms>1 ) * alpgenflavorflagbb =1;
  if (ncharms>1 ) * alpgenflavorflagcc =1;

  // Some checks:

  if ( (charm.size()!=anticharm.size()) || (bottom.size()!=antibottom.size())) {
    //    cerr<<"HEPTF-ERROR: Problems with the number of charms and bottoms: "
    //	<<charm.size()<<" "<<anticharm.size()<<" "<<bottom.size()<<" "
    //	<<antibottom.size()<<endl;
  }

  // We process the variables we use in the matching:


    // For the usual pair production the argument is for each quark to look
    // for the closest antiquark (in y-phi), and use the maximum of the distance.
    // For both, charm and bottom:

    *alpgendrcc = getDistanceForPairs(charm,anticharm);

    
    *alpgendrbb = getDistanceForPairs(bottom,antibottom);
    
    for ( unsigned int b = 0; b< bottom.size(); ++b   ){
      alpgenbpt->push_back( bottom[b]->pt() );
      alpgenbeta->push_back( bottom[b]->eta() );
      alpgenbphi->push_back( bottom[b]->phi() );
    }

    for ( unsigned int ab = 0; ab< antibottom.size(); ++ab   ){
      alpgenbpt->push_back( antibottom[ab]->pt() );
      alpgenbeta->push_back( antibottom[ab]->eta() );
      alpgenbphi->push_back( antibottom[ab]->phi() );
    }
    



 
 evt.put( alpgenflavorflagbb, "alpgenFlavorFlagBB" ); 
 evt.put( alpgenflavorflagcc, "alpgenFlavorFlagCC" ); 
 evt.put( alpgendrbb, "alpgenDrbb" ); 
 evt.put( alpgendrcc, "alpgenDrcc" );
 evt.put( alpgenvpt, "alpgenVpt" );  
 evt.put( alpgenbpt, "alpgenBpt" ); 
 evt.put( alpgenbeta, "alpgenBeta" ); 
 evt.put( alpgenbphi, "alpgenBphi" ); 
 
*pUintimeSize =-1;
 *pUouttime1plusSize =-1;
 *pUouttime1minusSize =-1;
 *pUouttime2plusSize =-1;
 *pUouttime2minusSize =-1;
 

// PileUp info
   Handle<std::vector< PileupSummaryInfo > >  PupInfo;
//   pile info only in MC....
  if (!evt.getByLabel(edm::InputTag("addPileupInfo"), PupInfo)) {
     LogWarning("") << ">>> pile up info does not exist !!!";
    
   }
   evt.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);

   std::vector<PileupSummaryInfo>::const_iterator PVI;
   for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
   if(PVI->getBunchCrossing()==0) *pUintimeSize= PVI->getPU_NumInteractions();
   if(PVI->getBunchCrossing()==1) *pUouttime1plusSize= PVI->getPU_NumInteractions();
   if(PVI->getBunchCrossing()==-1) *pUouttime1minusSize= PVI->getPU_NumInteractions();
   if(PVI->getBunchCrossing()==2) *pUouttime2plusSize= PVI->getPU_NumInteractions();
   if(PVI->getBunchCrossing()==-2) *pUouttime2minusSize= PVI->getPU_NumInteractions();
   if(PVI->getBunchCrossing()==3) *pUouttime3plusSize= PVI->getPU_NumInteractions();
   if(PVI->getBunchCrossing()==-3) *pUouttime3minusSize= PVI->getPU_NumInteractions();
   }

 evt.put( pUintimeSize,       "PUintimeSize" );
 evt.put( pUouttime1plusSize, "PUouttime1plusSize" ); 
 evt.put( pUouttime1minusSize,"PUouttime1minusSize" ); 
 evt.put( pUouttime2plusSize, "PUouttime2plusSize" ); 
 evt.put( pUouttime2minusSize,"PUouttime2minusSize" );
 evt.put( pUouttime3plusSize, "PUouttime3plusSize" ); 
 evt.put( pUouttime3minusSize,"PUouttime3minusSize" );

 
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( BMCTruthDumper );

