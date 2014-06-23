#if defined(__CINT__) && !defined(__MAKECINT__)
class loadFWLite {
  public:
    loadFWLite() {
      gSystem->Load("libFWCoreFWLite");
      AutoLibraryLoader::enable();
    }
};

static loadFWLite lfw;
#endif

#include "DataFormats/FWLite/interface/Handle.h"

#include <string>
#include <vector>
#include <iostream>

void triggerResultsByName_cint()
{
  std::vector<std::string> files;

//std::string filename = "dcache:/pnfs/cms/WAX/11/store/user/lpchbb/degrutto/ZH_ZToNuNu_HToBB_M-125_8TeV-powheg-herwigppSummer12_DR53X-PU_S10_START53_V7A-v1/degrutto/ZH_ZToNuNu_HToBB_M-125_8TeV-powheg-herwigpp/HBB_EDMNtupleV42/9803889241b1fc304f795d3b3875632d/PAT.edm_10_1_Rb5.root";
//std::string filename = "dcache:/pnfs/cms/WAX/11/store/user/lpchbb/degrutto/QCD_HT_250To500_TuneZ2star_8TeV/degrutto/QCD_HT-250To500_TuneZ2star_8TeV-madgraph-pythia6/QCD_HT_250To500_TuneZ2star_8TeV/93898e4b21312a4d4ac8927ff8d4ee62/PAT.edm_20_1_0wp.root";
std::string filename = "dcache:/pnfs/cms/WAX/11/store/mc/Summer12_DR53X/QCD_Pt-80to120_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v2/0000/0434D718-0AF4-E111-AC11-003048C6931C.root";


  files.push_back(filename);
  fwlite::ChainEvent ev(files);

  fwlite::Handle<edm::TriggerResults> hTriggerResults;


  bool expectedValue[4] = { true, false, true, true };
  int iEvent = 0;
  for (ev.toBegin(); ! ev.atEnd(); ++ev) {

    bool accept = false;
    edm::TriggerResultsByName resultsByName = ev.triggerResultsByName("HLT");
    std::vector<std::string> triggerNames = resultsByName.triggerNames();
    for (int i=0; i< triggerNames.size(); i++)
        std::cout << i << " " << triggerNames.at(i) << std::endl;

    //edm::TriggerResultsByName resultsByName = ev.triggerResultsByName("TEST");
    //if (resultsByName.isValid()) {
    //  std::cout << "From TriggerResultsByName, accept = "
    //            << resultsByName.accept("p") << "\n";
    //  accept = resultsByName.accept("p");
    //}
    //if (iEvent < 4 && expectedValue[iEvent] != accept) {
    //  std::cerr << "triggerResultsByName_cint.C, trigger results do not match expected values" << std::endl;
    //  abort();
    //}

    // Try again, but this time test what happens when the process does not
    // exist, the object should not be valid
    //edm::TriggerResultsByName resultsByName = ev.triggerResultsByName("DOESNOTEXIST");
    //if (resultsByName.isValid()) {
    //  std::cout << "From TriggerResultsByName, accept = "
    //            << resultsByName.accept("p") << "\n";
    //}
    //++iEvent;
    break;
  }
}
