#include "DataFormats/Common/interface/AssociationVector.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/Utilities/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "TH1.h"
#include "TH2.h"
#include <vector>
#include <string>
#include <iostream>
#include <sstream>


using namespace edm;
using namespace std;
using namespace reco;

 

typedef edm::AssociationVector<reco::CandidateRefProd, std::vector<double> > IsolationCollection;

class CosmicLHCAnalyzer: public edm::EDAnalyzer {
public:
  CosmicLHCAnalyzer(const edm::ParameterSet& pset);

private:
  virtual void analyze(const edm::Event& event, const edm::EventSetup& setup);
  virtual void endJob();
   void histo1D(const TH1F* , char* , char*) const ;
   void histo2D(const TH2F* , char* , char*) const ;
  InputTag src_muons ; 
  double ptCut , etaCut, massCut, dzCut, dxyCut;
  TH1F * Charge,* dxy,  *dz, *mass,  * ptMu, * etaMu, *size;
//*massAfterPtEtaCut, * massAfterPtEtaVtxCut ,
  TH2F * dz_vs_dxy, * dz_vs_dxy_afterCut;

 };

 void CosmicLHCAnalyzer::histo1D(const TH1F* hist,  char* cx, char*cy)  const {

  hist->GetXaxis()->SetTitle(cx);
  hist->GetYaxis()->SetTitle(cy);
 }

void CosmicLHCAnalyzer::histo2D(const TH2F* hist, char* cx, char*cy) const {

  hist->GetXaxis()->SetTitle(cx);
  hist->GetYaxis()->SetTitle(cy);
 }



CosmicLHCAnalyzer::CosmicLHCAnalyzer(const ParameterSet& pset):
  src_muons(pset.getParameter<InputTag>("src_muons")),
  ptCut(pset.getParameter<double>("ptCut")),
  etaCut(pset.getParameter<double>("etaCut")),
  massCut(pset.getParameter<double>("massCut")),
  dzCut(pset.getParameter<double>("dzCut")),
  dxyCut(pset.getParameter<double>("dxyCut"))

{
  edm::Service<TFileService> fs;
  //  cout<<"debug "<<endl;
  Charge = fs->make<TH1F>("Charge","Dimuons charge",2,-2,2);
  dxy = fs->make<TH1F>("dxy","ImpactParameter dxy",400,-200,200);
  dz = fs->make<TH1F>("dz","ImpactParameter dz",400,-200,200);
  mass = fs->make<TH1F>("mass","Dimuons mass",1000,0.,1000.);
  // massAfterPtEtaCut = fs->make<TH1F>("mass_pt_eta_cut","Dimuons mass after pt eta cut",200,0.,200.);
  //massAfterPtEtaVtxCut = fs->make<TH1F>("mass_pt_eta_vertex_cut","Dimuons mass after pt eta and vertex cut",200,0.,200.);
  ptMu = fs->make<TH1F>("ptMu","ptmu",200,0.,200.);
  etaMu = fs->make<TH1F>("etaMu","etaMu", 200, -6,6); 
  size = fs->make<TH1F>("size","muons size",10,0,10);
  dz_vs_dxy = fs->make<TH2F>("dz_vs_dxy","dz vs dxy",400,-200,200, 400,-200,200);  
  dz_vs_dxy_afterCut = fs->make<TH2F>("dz_vs_dxy_afterCut","dxy vs dz",400,-200,200, 400,-200,200);
}

void CosmicLHCAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  
  // looping on muons 
  Handle<TrackCollection> muons;
  event.getByLabel(src_muons,muons);
  std::cout<< "muons->size() = " << muons->size() << std::endl;
  size->Fill(muons->size());
  
  //  std::vector<math::XYZVector> v_posMuMomentum , v_negMuMomentum;
  std::vector<reco::Track> v_posMu, v_negMu; 
  for(size_t  i=0; i< muons->size(); ++ i ) {
    reco::Track muon = (*muons)[i] ;
    int q = muon.charge();
    Charge -> Fill(q);
    histo1D( Charge,  "q" , "#" );
    std::cout<< "lhclike standalone charge: " << q << std::endl; 
    if (q == 1){
      //   v_posMuMomentum.push_back(p);
      v_posMu.push_back(muon);
    }else {
      //  v_negMuMomentum.push_back(p);
      v_negMu.push_back(muon);
    } 
  }


  // evaluating dimuon invariant mass and all the variables
  for (size_t j=0; j <  v_posMu.size(); j++ ) {
    reco::Track mup = v_posMu[j] ;
    math::XYZVector pp = mup.momentum(); 
    //     std::cout<< "lhclike standalone momentum magnitude: " << mom.R() << std::endl; 
    double mupDxy= mup.dxy();

    double mupDz = mup.dz();
    dxy->Fill( mupDxy );
 
    dz->Fill(mupDz);
  
    dz_vs_dxy-> Fill(mupDxy,mupDz);
     histo2D( dz_vs_dxy,  "dxy (cm)", "dz (cm)");

    double mupPt = mup.pt();
    ptMu->Fill(mupPt);
 
    double mupEta = mup.eta();
    etaMu->Fill(mupEta);
 
    for (size_t k=0; k<  v_negMu.size(); k++ ){
      reco::Track mun = v_negMu[k] ;
      math::XYZVector pn = mun.momentum(); 
      //     std::cout<< "lhclike standalone momentum magnitude: " << mom.R() << std::endl; 
      double munDxy= mun.dxy();
      double munDz = mun.dz();
      dxy->Fill( munDxy );
      histo1D( dxy,   "dxy (cm) " , "#" );
      dz->Fill(munDz);
  histo1D( dz,  "dz (cm)", "#" );
      dz_vs_dxy-> Fill(munDxy,munDz);
      double munPt = mun.pt();
      ptMu->Fill(munPt);
   histo1D( ptMu ,  "p_{T} (GeV/c)" , "" );
      double munEta = mun.eta();
      etaMu->Fill(munEta);
   histo1D( etaMu ,  "#eta " , "" );
      // applying all the cut here 
      if ( mupPt > ptCut  && mupEta< etaCut  && munPt > ptCut  && munEta< etaCut){
	dz_vs_dxy_afterCut-> Fill(mupDxy,mupDz);
	dz_vs_dxy_afterCut-> Fill(munDxy,munDz);
	histo2D( dz_vs_dxy_afterCut,  "dxy (cm)", "dz (cm)");
	
	if (  mupDxy < dxyCut && mupDz < dzCut  && munDxy < dxyCut && munDz < dzCut){
	  //   double m = sqrt( 2 * ( posMuMomentum.R()* negMuMomentum.R() - posMuMomentum.x()* negMuMomentum.x()  - posMuMomentum.y()* negMuMomentum.y() - posMuMomentum.z()* negMuMomentum.z() ));
	  //     std::cout<< "lhclike dimuon standalone invariant mass: " << m << std::endl;    
	  double m = sqrt( 2 * ( pp.R()* pn.R() - pp.Dot( pn))) ;
      std::cout<< "lhclike dimuon standalone invariant mass: " << m<< std::endl;
      
      mass->Fill(m);
      histo1D( mass ,  "mass (GeV/c^{2})" , "" );
	}
      }
    }
  }
}
  
void CosmicLHCAnalyzer::endJob() {

}
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(CosmicLHCAnalyzer);
