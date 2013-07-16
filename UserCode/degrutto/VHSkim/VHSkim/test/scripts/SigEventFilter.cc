// -*- C++ -*-
//
// Package:    SigEventFilter
// Class:      SigEventFilter
// 
/**\class SigEventFilter SigEventFilter.cc MyCode/SigEventFilter/src/SigEventFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Mattew Fisher, mfisher@phys.ufl.edu
//         Created:  Sun Jul 19 11:14:58 EDT 2009
// $Id: SigEventFilter.cc,v 1.1 2011/07/12 14:02:54 madfish Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <string>
#include <iostream>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <DataFormats/L1GlobalMuonTrigger/interface/L1MuRegionalCand.h>

#include "MyCode/SigEventFilter/interface/Events.h"
#include "MyCode/SigEventFilter/interface/Runs.h"
#include "MyCode/SigEventFilter/interface/Lumis.h"
//
// class declaration
//

class SigEventFilter : public edm::EDFilter {
   public:
      explicit SigEventFilter(const edm::ParameterSet&);
      ~SigEventFilter();
      int totalEvts;
      int evtsWTF;
      int m;
      unsigned int run1,luminum,lumicnt;
    private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SigEventFilter::SigEventFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

  totalEvts=0;
  evtsWTF=0;
}


SigEventFilter::~SigEventFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
SigEventFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   totalEvts++;
   bool returnval = false;
   m=n;
   //NUMBER OF SIGNIFICANT EVENTS



   for(int i=0;i<n;i++)
   {
    if(/*(sigEvt[i]==iEvent.id().event()) &&*/ (sigRun[i]==iEvent.id().run()) && (sigLumi[i]==iEvent.luminosityBlock()) )
     returnval = true;
   }  

   luminum = iEvent.luminosityBlock();

   if(returnval)
     {
     evtsWTF++;
     } 
     std::cout <<"Passed Filter: " << evtsWTF << " out of "<< totalEvts <<std::endl;
    return(returnval);

    

}

// ------------ method called once each job just before starting event loop  ------------
void 
SigEventFilter::beginJob()
{
   run1 = 0;
   luminum=0;
//   std::cout << "#Runs: " ;
   lumicnt=0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SigEventFilter::endJob()
{
//        std::cout << luminum << ") ";
 
    std::cout << std::endl;
    std::cout <<"Number Lumis: " << lumicnt <<std::endl;
    std::cout <<"Passed Filter: " << evtsWTF << " out of "<< totalEvts <<std::endl;
    std::cout <<"Found : " << evtsWTF << " out of "<< m << " Events" <<std::endl;
}


//define this as a plug-in
DEFINE_FWK_MODULE(SigEventFilter);
