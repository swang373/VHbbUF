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
#include "DataFormats/Math/interface/deltaPhi.h"


int main(int argc, char* argv[]) 
{

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
   
  double thirdJetPt_cut_( ana.getParameter<double>("thirdJetPt_cut") );  
  unsigned int nofLeptons_cut_( ana.getParameter<unsigned int>("nofLeptons_cut") );  
  double j1pt_cut_( ana.getParameter<double>("j2pt_cut") );  
  double j2pt_cut_( ana.getParameter<double>("j2pt_cut") );  
  double jCSVMVA_cut_( ana.getParameter<double>("jCSVMVA_cut") );   
  //double jCSV_cut_( ana.getParameter<double>("jCSV_cut") );   
  double jCSVH_cut_( ana.getParameter<double>("jCSVH_cut") );   
  double jCSVL_cut_( ana.getParameter<double>("jCSVL_cut") ); 
  double jbProb_cut_( ana.getParameter<double>("jbProb_cut") ); 
  double jTCHP_cut_( ana.getParameter<double>("jTCHP_cut") ); 
  double jjDeltaRmin_cut_( ana.getParameter<double>("jjDeltaRmin_cut") ); 
  double jjDeltaPhimin_cut_( ana.getParameter<double>("jjDeltaPhimin_cut") ); 
  double hjjmetDeltaPhimin_cut_( ana.getParameter<double>("hjjmetDeltaPhimin_cut") ); 
  double jjDeltaRmax_cut_( ana.getParameter<double>("jjDeltaRmax_cut") ); 
  double pfMet_cut_( ana.getParameter<double>("pfMet_cut") ); 
  double pfMetSig_cut_( ana.getParameter<double>("pfMetSig_cut") ); 
  double hjjpt_cut_( ana.getParameter<double>("hjjpt_cut") );  
  unsigned int nHjj_cut_( ana.getParameter<unsigned int >("nHjj_cut") );  
  double zmmmassMin_cut_( ana.getParameter<double>("zmmmassMin_cut") );  
  double zmmmassMax_cut_( ana.getParameter<double>("zmmmassMax_cut") );  
  double zmmptMin_cut_( ana.getParameter<double>("zmmptMin_cut") );  
  double hjjzmmDeltaPhimin_cut_( ana.getParameter<double>("hjjzmmDeltaPhimin_cut") ); 

  double zeemassMin_cut_( ana.getParameter<double>("zeemassMin_cut") );   
  double zeemassMax_cut_( ana.getParameter<double>("zeemassMax_cut") );  
  double zeeptMin_cut_( ana.getParameter<double>("zeeptMin_cut") );  
  double hjjzeeDeltaPhimin_cut_( ana.getParameter<double>("hjjzeeDeltaPhimin_cut") ); 

  double wmnmtMin_cut_( ana.getParameter<double>("wmnmtMin_cut") );   
  double wmnptMin_cut_( ana.getParameter<double>("wmnptMin_cut") );  
  double hjjwmnDeltaPhimin_cut_( ana.getParameter<double>("hjjwmnDeltaPhimin_cut") ); 

  double wenmtMin_cut_( ana.getParameter<double>("wenmtMin_cut") );   
  double wenptMin_cut_( ana.getParameter<double>("wenptMin_cut") );  
  double hjjwenDeltaPhimin_cut_( ana.getParameter<double>("hjjwenDeltaPhimin_cut") ); 
  bool isZinv_( ana.getParameter<bool>("isZinv") );  
  bool isWln_( ana.getParameter<bool>("isWln") );  
   bool isZll_( ana.getParameter<bool>("isZll") );  
  //  bool isZinv_( ana.getParameter<bool>("isZinv") );  

  // book a set of histograms
  fwlite::TFileService fs = fwlite::TFileService(outputFile_.c_str());
  TFileDirectory dir = fs.mkdir("analyzeBasicPat");
  TH1F* hjjMass_  = dir.make<TH1F>("hjjMass"  , "mass"  ,   500,   0.,  500.);
  TH1F* hjjzllMass_  = dir.make<TH1F>("hjjzllMass"  , "jjll mass"  ,   500,   0.,  500.);
  TH1F* hjjzllPt_  = dir.make<TH1F>("hjjzllPt"  , "jjll pt"  ,   500,   0.,  500.);
  TH1F* hjjMetMt_  = dir.make<TH1F>("hjjmeatMt"  , "jjmet mt"  ,   1000,   0.,  500.);
  TH1F* hjjMetPt_  = dir.make<TH1F>("hjjmetPt"  , "jjmet pt"  ,   500,   0.,  500.);
  TH1F* diffhjjMetPt_  = dir.make<TH1F>("diffhjjzllPt"  , "#Delta(jj pt  - met)"  ,   1000,   -500.,  500.);
  TH1F* nHjjbef_  = dir.make<TH1F>("nHjjbef"  , "# of jj pair before"  ,   20,   -0.5,  19.5);
  TH1F* nHjjaft_  = dir.make<TH1F>("nHjjaft"  , "# of jj pair after"  ,   10,   -0.5,  9.5);
  TH1F* hjjPt_  = dir.make<TH1F>("hjjPt"  , "hjjPt"  ,   500,   0.,  1000.);
  TH1F* j1Pt_  = dir.make<TH1F>("j1Pt"  , "j1Pt"  ,   500,   0.,  1000.);
  TH1F* j2Pt_  = dir.make<TH1F>("j2Pt"  , "j2Pt"  ,   500,   0.,  1000.);
  TH1F* j1Eta_  = dir.make<TH1F>("j1Eta"  , "j1Eta"  ,   500,   -5.,  5.);
  TH1F* j2Eta_  = dir.make<TH1F>("j2Eta"  , "j2Eta"  ,   500,   -5.,  5.);

  TH1F* jjdr_  = dir.make<TH1F>("jjdr"  , "deltaR"  ,   100,   0.,  5.);
  TH1F* jjdPhi_  = dir.make<TH1F>("jjdPhi"  , "deltaPhi"  ,   100,   0, 5.);

  TH1F* pullAK5_  = dir.make<TH1F>("pullAK5"  , "pull AK5"  ,   100,   -5., 5.);
  //  TH1F* pullAK7_  = dir.make<TH1F>("pullAK7"  , "pull AK7"  ,   100,   -5., 5.);
  TH1F* hele1_  = dir.make<TH1F>("hele1"  , "helicityjet1"  ,   110,   0., 1.1);
  TH1F* hele2_  = dir.make<TH1F>("hele2"  , "helicityjet2"  ,   110,   0., 1.1);

  TH1F* jjdEta_  = dir.make<TH1F>("jjdEta"  , "deltaEta"  ,   100,   0.,5.);
  TH1F* hjjmetdPhi_  = dir.make<TH1F>("hjjmetdPhi"  , "hjjmetdeltaPhi"  ,   100,   0.,  5.);
  TH1F* j1metdPhi_  = dir.make<TH1F>("j1metdPhi"  , "j1metdeltaPhi"  ,   100,   0.,  5.);
  TH1F* j2metdPhi_  = dir.make<TH1F>("j2metdPhi"  , "j1metdeltaPhi"  ,   100,   0.,  5.);
  TH1F* jcsvmva_  = dir.make<TH1F>("jcvsmva"  , "jet csvmva"  ,   100,   0,  1);
  TH1F* jcsv_  = dir.make<TH1F>("jcvs"  , "jet csv"  ,   100,   0,  1);
  TH1F* jbprob_  = dir.make<TH1F>("jbprob"  , "jet bprob"  ,   100,   0,  10);
  TH1F* jTCHP_  = dir.make<TH1F>("jTCHP"  , "jet TCHP"  ,   1000,   0,  100);
  TH1F* pfmet_  = dir.make<TH1F>("pfmet"  , "pfmet"  ,   500,   0,  1000);
  TH1F* pfmetsig_  = dir.make<TH1F>("pfmetsig"  , "pfmetsig"  ,   100,   0,  100);
  TH1F* zmmMass_  = dir.make<TH1F>("zmmMass"  , "zmmmass"  ,   500,   0.,  2000.);
  TH1F* zmmPt_  = dir.make<TH1F>("zmmPt"  , "zmmpt"  ,   500,   0.,  2000.);
  TH1F* hjjzmmdPhi_  = dir.make<TH1F>("hjjzmmdPhi"  , "hjjzmmdeltaPhi"  ,   100,   -5.,  5.);

  TH1F* zeeMass_  = dir.make<TH1F>("zeeMass"  , "zeemass"  ,   500,   0.,  2000.);
  TH1F* zeePt_  = dir.make<TH1F>("zeePt"  , "zeept"  ,   500,   0.,  2000.);
  TH1F* hjjzeedPhi_  = dir.make<TH1F>("hjjzeedPhi"  , "hjjzeedeltaPhi"  ,   100,   -5.,  5.);

  TH1F* wmnMt_  = dir.make<TH1F>("wmnMt"  , "wmnmt"  ,   500,   0.,  2000.);
  TH1F* wmnPt_  = dir.make<TH1F>("wmnPt"  , "wmnpt"  ,   500,   0.,  2000.);
  TH1F* hjjwmndPhi_  = dir.make<TH1F>("hjjwmndPhi"  , "hjjwmndeltaPhi"  ,   100,   -5.,  5.);

  TH1F* wenMt_  = dir.make<TH1F>("wenMt"  , "wenmt"  ,   500,   0.,  2000.);
  TH1F* wenPt_  = dir.make<TH1F>("wenPt"  , "wenpt"  ,   500,   0.,  2000.);
  TH1F* hjjwendPhi_  = dir.make<TH1F>("hjjwendPhi"  , "hjjwendeltaPhi"  ,   100,   -5.,  5.);
  TH1F* thirdJetPtBef_  = dir.make<TH1F>("thirdJetPtBef"  , "thirdJetPt"  ,   1000,   0,  10000);
  TH1F* thirdJetPtAft_  = dir.make<TH1F>("thirdJetPtAft"  , "thirdJetPt"  ,   1000,   0,  10000);
  TH1F* thirdJetEta_  = dir.make<TH1F>("thirdJetEta"  , "thirdJetEta"  ,   1000,   -5,  5);
  TH1F* jetsPt_  = dir.make<TH1F>("jetsPt"  , "jetsPt"  ,   1000,   0,  10000);
  TH1F* nofLeptons_  = dir.make<TH1F>("nofLeptons"  , "# of leptons"  ,   10,   -0.5,  9.5);
  TH1F* nofJets_  = dir.make<TH1F>("nofJets"  , "# of jets"  ,   10,   -0.5,  9.5);
  TH1F* nofBJets_  = dir.make<TH1F>("nofBJets"  , "# of bjets"  ,   10,   -0.5,  9.5);

  
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
	if(outputEvery_!=0 ? (ievt>0 && ievt%outputEvery_==0) : false); 
	// std::cout << "  processing event: " << ievt << std::endl;


	// ######### looking first at jets multiplicity     
	fwlite::Handle< std::vector<float> > JetsPts;
	JetsPts.getByLabel(ev,"JetLepEdmNtuple", "JetsPt");
	if (!JetsPts.isValid()) continue;
        nofJets_->Fill(JetsPts->size());

	///##### stand by now for the b-tag... I need to add something in the ntuple....

  
	for ( unsigned int j = 0 ; j < JetsPts->size(); ++j  ){
	  // jets are orders in pt
	  jetsPt_->Fill((*JetsPts)[j]);
	}

	fwlite::Handle< std::vector<float> > JetsEtas;
	JetsEtas.getByLabel(ev,"JetLepEdmNtuple", "JetsEta");
	
	
	//applying immediately the  lepton veto
	unsigned int nofMu =0, nofEle =0;
    	fwlite::Handle< unsigned int > nofMuons;
	nofMuons.getByLabel(ev,"JetLepEdmNtuple", "NofMuons");
        if (!nofMuons.isValid())  nofMu =  *nofMuons ;

    	fwlite::Handle< unsigned int > nofElectrons;
	nofElectrons.getByLabel(ev,"JetLepEdmNtuple", "NofElectrons");
        if (!nofElectrons.isValid()) 	 nofEle =  *nofElectrons ;
        unsigned int nofLeptons = nofMu + nofEle;
	//	std::cout << "n of lept " << nofLeptons << std::endl;

	nofLeptons_->Fill(nofLeptons);  
	if ( nofLeptons > nofLeptons_cut_) continue;

	// ############# jet var needed for jet veto ########## 
	float thirdJetPt = -1; 
	float thirdJetEta = -9999; 
	if (JetsPts->size()>2 ) {
	   thirdJetPt = (*JetsPts)[2];
	   thirdJetEta = (*JetsEtas)[2]; 
	}
	  

	// ############# Higgs cand loop ###################
	// Handle to collection
	fwlite::Handle<std::vector<float> > objs;
	objs.getByLabel(ev,"HjjEdmNtuple", "hjjMass");
        if (!objs.isValid()) continue;
        unsigned int nHjj = objs->size();
	nHjjbef_->Fill( nHjj);
	  if (nHjj>= nHjj_cut_) continue;
	  // protection on too much jj pairs....
	  nHjj<1000?  1 : nHjj = 999 ; 

 
	for(unsigned int i = 0; i < nHjj; ++i){

	  // ########### the first thing we cut on is the b-tag!!!!!!!!!!!!!!! #################################
	  fwlite::Handle<std::vector<float> > hjjJet1CSVMVAs;
	  hjjJet1CSVMVAs.getByLabel(ev,"HjjEdmNtuple", "hjjJet1CSVMVA");
	  float csvmva1 = (*hjjJet1CSVMVAs)[i];

	  fwlite::Handle<std::vector<float> > hjjJet2CSVMVAs;
	  hjjJet2CSVMVAs.getByLabel(ev,"HjjEdmNtuple", "hjjJet2CSVMVA");
	  float csvmva2 = (*hjjJet2CSVMVAs)[i];

	  if ((csvmva1 > jCSVMVA_cut_ &&  csvmva2>jCSVMVA_cut_)==0) continue;


	  fwlite::Handle<std::vector<float> > hjjJet1CSVs;
	  hjjJet1CSVs.getByLabel(ev,"HjjEdmNtuple", "hjjJet1CSV");
	  float csv1 = (*hjjJet1CSVs)[i];

	  fwlite::Handle<std::vector<float> > hjjJet2CSVs;
	  hjjJet2CSVs.getByLabel(ev,"HjjEdmNtuple", "hjjJet2CSV");
	  float csv2 = (*hjjJet2CSVs)[i];

	  //	  if ((csv1 > jCSV_cut_ &&  csv2>jCSV_cut_)==0) continue;
	  if(((csv1 > jCSVH_cut_ &&  csv2>jCSVL_cut_)==0) && ((csv1 > jCSVL_cut_ && csv2>jCSVH_cut_) == 0)) continue;

	  fwlite::Handle<std::vector<float> > hjjJet1JbProbs;
	  hjjJet1JbProbs.getByLabel(ev,"HjjEdmNtuple", "hjjJet1JbProb");
	  float bprob1 = (*hjjJet1JbProbs)[i];

	  fwlite::Handle<std::vector<float> > hjjJet2JbProbs;
	  hjjJet2JbProbs.getByLabel(ev,"HjjEdmNtuple", "hjjJet2JbProb");
	  float bprob2 = (*hjjJet2JbProbs)[i];

	  if ((bprob1>jbProb_cut_ && bprob2>jbProb_cut_)==0) continue; 

	  fwlite::Handle<std::vector<float> > hjjJet1TCHPs;
	  hjjJet1TCHPs.getByLabel(ev,"HjjEdmNtuple", "hjjJet1TKHP");
	  float jTCHP1 = (*hjjJet1TCHPs)[i];

	  fwlite::Handle<std::vector<float> > hjjJet2TCHPs;
	  hjjJet2TCHPs.getByLabel(ev,"HjjEdmNtuple", "hjjJet2TKHP");
	  float jTCHP2 = (*hjjJet2TCHPs)[i];

	  if ((jTCHP1 > jTCHP_cut_ && jTCHP2 > jTCHP_cut_)==0) continue;


	  // ######################### end of b-tag ############################################

	  float mass = (*objs)[i];

	  // if(fMin < mass && mass < fMax)

	  // looking to daughter jet jet with pt>ptcut  
	  fwlite::Handle<std::vector<float> > hjjPts;
	  hjjPts.getByLabel(ev,"HjjEdmNtuple", "hjjPt");
	  float hjjpt = (*hjjPts)[i];



          if ((hjjpt> hjjpt_cut_ )==0) continue;



	  fwlite::Handle<std::vector<float> > hjjPhis;
	  hjjPhis.getByLabel(ev,"HjjEdmNtuple", "hjjPhi");
	  float hjjphi = (*hjjPhis)[i];



	  float hjjzllmass =99999;
	  fwlite::Handle<std::vector<float> > hjjzllMasses;
	  hjjzllMasses.getByLabel(ev,"HjjZllEdmNtuple", "hjjzllMass");
	  if (hjjzllMasses.isValid())  hjjzllmass = (*hjjzllMasses)[i];

	  float hjjzllpt = 99999;
	  fwlite::Handle<std::vector<float> > hjjzllPts;
	  hjjzllPts.getByLabel(ev,"HjjZllEdmNtuple", "hjjzllPt");
	  if (hjjzllPts.isValid())  hjjzllpt = (*hjjzllPts)[i];



	  fwlite::Handle<std::vector<float> > hjjJet1Pts;
	  hjjJet1Pts.getByLabel(ev,"HjjEdmNtuple", "hjjJet1Pt");
	  float pt1 = (*hjjJet1Pts)[i];



          fwlite::Handle<std::vector<float> > hjjJet2Pts;
	  hjjJet2Pts.getByLabel(ev,"HjjEdmNtuple", "hjjJet2Pt");
	  float pt2= (*hjjJet2Pts)[i];

	  if (((pt1> j1pt_cut_ && pt2> j2pt_cut_) ||(pt1> j2pt_cut_ && pt2> j1pt_cut_ ))==0) continue;


          //now applying jet veto.... looking at the jet with pt > j1 and/or j2
 
          bool found = false;
	  //    ###########################  third jet pt: either the thurd jet or the fiusrt or second if greater than higgs jets ##################




	    for (  unsigned int j = 0 ; j < JetsPts->size(); ++j  ){
	      if (found==true) continue; 
 	      // jets are orders in pt
		float thirdJetPtTemp =   (*JetsPts)[j];
		float thirdJetEtaTemp =   (*JetsEtas)[j];
		//  std::cout << " thirdJetPtTemp == " << thirdJetPtTemp << ", while pt1 and pt2 "<< pt1 << ", " << pt2 << std::endl;  
		if ( (thirdJetPtTemp!= pt1 && thirdJetPtTemp!= pt2 ) && ((thirdJetPtTemp > pt1) || (thirdJetPtTemp > pt2))) {
		  thirdJetPt = thirdJetPtTemp;
		  thirdJetEta = thirdJetEtaTemp;
		  found = true; 
		}
              
	    }

	  /*
	    std::cout << "3rd jet pt  " <<  thirdJetPt<< std::endl;
	  
	   if (thirdJetPt> -1 ) {
	   std::cout << "for Hjj[" << i << "] we "<< "found in the event a jet with pt, eta " <<  thirdJetPt <<", " <<   thirdJetEta << 
	  "harder than the j1, j2 pt " << pt1 << ", "<< pt2 << std:: endl;
	   }
	  
	  */
	  
	  

	  
	  fwlite::Handle<std::vector<float> > hjjJet1Etas;
	  hjjJet1Etas.getByLabel(ev,"HjjEdmNtuple", "hjjJet1Eta");
	  float eta1 = (*hjjJet1Etas)[i];

	  fwlite::Handle<std::vector<float> > hjjJet2Etas;
	  hjjJet2Etas.getByLabel(ev,"HjjEdmNtuple", "hjjJet2Eta");
	  float eta2= (*hjjJet2Etas)[i];

	  fwlite::Handle<std::vector<float> > hjjJet1Phis;
	  hjjJet1Phis.getByLabel(ev,"HjjEdmNtuple", "hjjJet1Phi");
	  float phi1 = (*hjjJet1Phis)[i];

	  fwlite::Handle<std::vector<float> > hjjJet2Phis;
	  hjjJet2Phis.getByLabel(ev,"HjjEdmNtuple", "hjjJet2Phi");
	  float phi2= (*hjjJet2Phis)[i];




	  float pfmet = 0;

	  fwlite::Handle<std::vector<float> > pfmets;
	  pfmets.getByLabel(ev,"MetEdmNtuple", "pfMetEt");
	  pfmet = (*pfmets)[0];
	  
          if ((pfmet > pfMet_cut_)==0) continue; 

	  float pfmetsig = 0;

	  fwlite::Handle<std::vector<float> > pfmetsigs;
	  pfmetsigs.getByLabel(ev,"MetEdmNtuple", "pfMetMEtSig");
	  pfmetsig = (*pfmetsigs)[0];

          if ((pfmetsig > pfMetSig_cut_)==0) continue; 

	  float pfmetphi = 0;
	  fwlite::Handle<std::vector<float> > pfmetphis;
	  pfmetphis.getByLabel(ev,"MetEdmNtuple", "pfMetPhi");
	  pfmetphi = (*pfmetphis)[0];

	  float jjpfmetMt = 0;
	  fwlite::Handle<std::vector<float> > jjpfmetmts;
	  jjpfmetmts.getByLabel(ev,"HjjMetEdmNtuple", "hjjmetMT");
          if (jjpfmetmts.isValid())   jjpfmetMt = (*jjpfmetmts)[i];

	  float jjpfmetEt = 0;
	  fwlite::Handle<std::vector<float> > jjpfmetets;
	  jjpfmetets.getByLabel(ev,"HjjMetEdmNtuple", "hjjmetEt");
	  if (jjpfmetets.isValid()) jjpfmetEt = (*jjpfmetets)[i];




  
	  float dR = deltaR( eta1, phi1, eta2, phi2);
	  if ((dR > jjDeltaRmin_cut_ && dR < jjDeltaRmax_cut_)==0) continue; 
  
	  float dPhi = fabs(deltaPhi(phi1, phi2));
	  if ((dPhi > jjDeltaPhimin_cut_)==0) continue;

	  float dPhiHjjMet = fabs(deltaPhi(hjjphi, pfmetphi));
	  if ((dPhiHjjMet > hjjmetDeltaPhimin_cut_)==0) continue;

	  // ############ looking at pool and boost angle

          float pullAK5 = 0;
	  fwlite::Handle<std::vector<float> > pullsak5;
	  pullsak5.getByLabel(ev,"SVtxEdmNtuple", "hjjJet1DeltaTheta");
	  if (pullsak5.isValid()) pullAK5 = (*pullsak5)[i];


	  /*	
          float pullAK7 = 0;
	  fwlite::Handle<std::vector<float> > pullsak7;
	  std::cout << "before  AK7 " << std::endl;
  	  pullsak7.getByLabel(ev,"SVtxEdmNtuple", "hjjJet1DeltaThetaAK7");
	  std::cout << "after  AK7 " << std::endl;
	  if (pullsak7.isValid()) pullAK7 = (*pullsak7)[i];
	  std::cout << "after  AK7 is valid" << std::endl;
	  */


          float hele1 = 0;
	  fwlite::Handle<std::vector<float> > heles1;
	  heles1.getByLabel(ev,"SVtxEdmNtuple", "hjjJet1higgsHelicity");
	  if (heles1.isValid()) hele1 = fabs((*heles1)[i]);


          float hele2 = 0;
	  fwlite::Handle<std::vector<float> > heles2;
	  heles2.getByLabel(ev,"SVtxEdmNtuple", "hjjJet2higgsHelicity");
	  if (heles2.isValid()) hele2 = fabs((*heles2)[i]);



	  

	  
	  /// now looking at Zee, Zmm, Wmv, Wme if they exists
	  // Zmm
          bool isZmm = false;
          bool isWmn = false;
          bool isZee = false;
          bool isWen = false;


          float zmmMass = 0;
	  float zmmPt = 0;
	  float zmmPhi = 0;
	  fwlite::Handle<std::vector<float> > zmmMasses;
	  zmmMasses.getByLabel(ev,"ZmmEdmNtuple", "zmmMass");
          if ( zmmMasses.isValid()&& isZll_  ) {
	    
	    for  (unsigned int j = 0; j < zmmMasses->size(); ++j){
	      // choosing the Zmm closer to the zmass  
	      float zmmMassTmp = (*zmmMasses)[j];
	      
	      if (fabs(zmmMassTmp - 91) < fabs(zmmMass -91 )){
		zmmMass = zmmMassTmp;
		// now get Pt and Phi 
		fwlite::Handle<std::vector<float> > zmmPts;
		zmmPts.getByLabel(ev,"ZmmEdmNtuple", "zmmPt");
		zmmPt = (*zmmPts)[j];
		
		fwlite::Handle<std::vector<float> > zmmPhis;
		zmmPhis.getByLabel(ev,"ZmmEdmNtuple", "zmmPhi");
		zmmPhi = (*zmmPhis)[j];
		
	      }
	    }

	    float dPhiHjjZmm = fabs(deltaPhi(hjjphi, zmmPhi));
	    if (( zmmMass > zmmmassMin_cut_ ) && (zmmMass < zmmmassMax_cut_)   &&   (zmmPt > zmmptMin_cut_ ) && (dPhiHjjZmm > hjjzmmDeltaPhimin_cut_) ) {
	      
	      zmmMass_->Fill(zmmMass);
	      zmmPt_->Fill(zmmPt);
	      hjjzmmdPhi_->Fill(dPhiHjjZmm);
              isZmm = true; 
	    }
	  }

	  // Zee
          float zeeMass = 0;
	  float zeePt = 0;
	  float zeePhi = 0;
	  fwlite::Handle<std::vector<float> > zeeMasses;
	  zeeMasses.getByLabel(ev,"ZeeEdmNtuple", "zeeMass");
          if ( zeeMasses.isValid()  && isZll_) { 
	    
	    for  (unsigned int j = 0; j < zeeMasses->size(); ++j){
	      // choosing the Zee closer to the zmass  
	      float zeeMassTmp = (*zeeMasses)[j];
	      
	      if (fabs(zeeMassTmp - 91) < fabs(zeeMass -91 )){
		zeeMass = zeeMassTmp;
		// now get Pt and Phi 
		fwlite::Handle<std::vector<float> > zeePts;
		zeePts.getByLabel(ev,"ZeeEdmNtuple", "zeePt");
		zeePt = (*zeePts)[j];
		
		fwlite::Handle<std::vector<float> > zeePhis;
		zeePhis.getByLabel(ev,"ZeeEdmNtuple", "zeePhi");
		zeePhi = (*zeePhis)[j];
		
	      }
	    }
	    
	    
	    float dPhiHjjZee = fabs(deltaPhi(hjjphi, zeePhi));
	    if ((zeeMass > zeemassMin_cut_ ) && (zeeMass < zeemassMax_cut_)   && (zeePt > zeeptMin_cut_ ) && (dPhiHjjZee > hjjzeeDeltaPhimin_cut_)){
	      
	      zeeMass_->Fill(zeeMass);
	      zeePt_->Fill(zeePt);
	      hjjzeedPhi_->Fill(dPhiHjjZee);
              isZee = true; 
	    }
	  }
	  // Wen
          float wenMt = 0;
	  float wenPt = 0;
	  float wenPhi = 0;
	  
	  fwlite::Handle<std::vector<float> > wenMts;
	  wenMts.getByLabel(ev,"WenEdmNtuple", "wenMT");
	  
	  
          if ( wenMts.isValid() && isWln_) {
	    
	    
	    //  std::cout << " wenMts->size() " <<wenMts->size() << std::endl;
	    for  (unsigned int j = 0; j < wenMts->size(); ++j){

	      // choosing the Wen closer to the zmass  
	      float wenMtTmp = (*wenMts)[j];
	      
	      if (fabs(wenMtTmp - 80) < fabs(wenMt -80 )){
		wenMt = wenMtTmp;
		// now get Pt and Phi 
		fwlite::Handle<std::vector<float> > wenPts;
		wenPts.getByLabel(ev,"WenEdmNtuple", "wenPt");
		wenPt = (*wenPts)[j];
		
		fwlite::Handle<std::vector<float> > wenPhis;
		wenPhis.getByLabel(ev,"WenEdmNtuple", "wenPhi");
		wenPhi = (*wenPhis)[j];
		
	      }
	    }
	    
	    
	    float dPhiHjjWen = fabs(deltaPhi(hjjphi, wenPhi));
	    if ((wenMt > wenmtMin_cut_ ) && (wenPt > wenptMin_cut_ ) && (dPhiHjjWen > hjjwenDeltaPhimin_cut_)) {
	      
	      wenMt_->Fill(wenMt);
	      wenPt_->Fill(wenPt);
	      hjjwendPhi_->Fill(dPhiHjjWen);
              isWmn = true; 
	    }
	  }
	  
	  // Wmn
          float wmnMt = 0;
	  float wmnPt = 0;
	  float wmnPhi = 0;
	  fwlite::Handle<std::vector<float> > wmnMts;
	  
	  wmnMts.getByLabel(ev,"WmnEdmNtuple", "wmnMT");
	  
          if ( wmnMts.isValid() && isWln_ ) {
	    // std::cout << " wmnMts->size() " <<wmnMts->size() << std::endl;
	    
	    for  (unsigned int j = 0; j < wmnMts->size(); ++j){
	      //	    std::cout << " wmnMts->size() " <<wmnMts->size() << std::endl;
	      // choosing the Wmn closer to the zmass  
	      float wmnMtTmp = (*wmnMts)[j];
	      
	      if (fabs(wmnMtTmp - 80) < fabs(wmnMt -80 )){
		wmnMt = wmnMtTmp;
		// now get Pt and Phi 
		fwlite::Handle<std::vector<float> > wmnPts;
		wmnPts.getByLabel(ev,"WmnEdmNtuple", "wmnPt");
		wmnPt = (*wmnPts)[j];
		
		fwlite::Handle<std::vector<float> > wmnPhis;
		wmnPhis.getByLabel(ev,"WmnEdmNtuple", "wmnPhi");
		wmnPhi = (*wmnPhis)[j];
		
	      }
	    }
     	  
	    float dPhiHjjWmn = fabs(deltaPhi(hjjphi, wmnPhi));
	    if ((wmnMt > wmnmtMin_cut_ ) && (wmnPt > wmnptMin_cut_ ) && (dPhiHjjWmn > hjjwmnDeltaPhimin_cut_)){
	      
	      wmnMt_->Fill(wmnMt);
	      wmnPt_->Fill(wmnPt);
	      hjjwmndPhi_->Fill(dPhiHjjWmn);
              isWen = true; 
	    }
	  }

	  //if not isZll and we found a Z instead it discard the event (needed to suppress ttbar....):
           
	  // now filling hjj plots if at leat a Z/W  has been found
	  // discard cases in which we found a W and Z at the same time, it helps for ttbar subtraction, but we need also to apply a smarter jet veto 
	  // if is ZinvH don't care... and fill anyway
          bool isW = (isWmn  || isWen) ;
          bool isZ = (isZmm  || isZee) ;


	  if (isZ) {
	    hjjzllMass_->Fill(hjjzllmass);
	    hjjzllPt_->Fill(hjjzllpt);
	  }


          if (isZinv_) {
	    hjjMetMt_->Fill(jjpfmetMt);  
	    hjjMetPt_ ->Fill(jjpfmetEt); 
	    diffhjjMetPt_ ->Fill( hjjpt - pfmet );
	    j1metdPhi_ -> Fill (deltaPhi(  phi1, pfmetphi )); 
	    j2metdPhi_ -> Fill (deltaPhi(  phi2, pfmetphi )); 
  
	  }

          if ( (isZinv_ || isZ || isW ))// && (isZinv_ || !( isZ && isW)) )
	    {
          nHjjaft_->Fill( 1);
          j1Pt_->Fill(pt1);
          j2Pt_->Fill(pt2);
	  
	  j1Eta_->Fill(eta1);
	  j2Eta_->Fill(eta2);


          pullAK5_->Fill(pullAK5);
	  //	  pullAK7_->Fill(pullAK7);
          hele1_->Fill(hele1);
          hele2_->Fill(hele2);

	//#### now it make senses only for nBJets>=2 and with simmetric cuts on b-tag.........
	  if ( i ==0 ) nofBJets_->Fill(2);
          if (i>0) nofBJets_->Fill(1);

	  /*
	    std::cout << "for Hjj[" << i << "] we "<< " have that the third jet pt and  eta are "  <<  thirdJetPt <<", " <<   thirdJetEta <<  std:: endl;

	    std::cout << "the evnts has jets with pt  "; 
	for ( unsigned int j = 0 ; j < JetsPts->size(); ++j  ){
	// jets are orders in pt
	std::cout << (*JetsPts)[j] << " " ;

	}
	std::cout << " and Hjjj1 and hjjj2 pt:  " << pt1 << ", " << pt2 << std::endl;
	  */

	  thirdJetPtBef_->Fill(thirdJetPt);
	  if (thirdJetPt>  thirdJetPt_cut_ ) continue;
	    thirdJetPtAft_->Fill(thirdJetPt);
	    

	  hjjMass_->Fill(mass);
	  hjjPt_->Fill(hjjpt);
         



	 	  
	  jjdr_->Fill(dR);
	  jjdPhi_->Fill(dPhi);
	  jjdEta_->Fill(eta1 - eta2);

	

	  hjjmetdPhi_->Fill(dPhiHjjMet);
	  jcsvmva_->Fill(csvmva1);
	  jcsvmva_->Fill(csvmva2);

	  jcsv_->Fill(csv1);
	  jcsv_->Fill(csv2);
	  
	  jbprob_->Fill(bprob1);
	  jbprob_->Fill(bprob2);
	  
	  jTCHP_->Fill(jTCHP1);
	  jTCHP_->Fill(jTCHP2);
	  
	  
	  pfmet_->Fill(pfmet);
	  pfmetsig_->Fill(pfmetsig);



	 
	  thirdJetEta_->Fill(thirdJetEta );



	  }
	
	
      }
    }


      
    inFile->Close();
  }
  // break loop if maximal number of events is reached:
    // this has to be done twice to stop the file loop as well
  if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
}
return 0;
}
