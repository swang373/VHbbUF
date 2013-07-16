#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include "TROOT.h"
#include <vector>

output(TString fname)
{
TString rootname,filename;
rootname = fname;
rootname.ReplaceAll("txt","root");
std::cout << fname.Data() << std::endl;
TFile *outfile = new TFile(rootname.Data(),"RECREATE");
TTree *tree = new TTree("tree","tree");
filename = fname;
int nlines = (int) tree->ReadFile(fname.Data(),"ptMin/F:ptMax/F:etaMin/F:etaMax/F:scale/F:error/F");

std::cout<<"Read "<<nlines<<" lines of data."<<std::endl;
tree->Write(); outfile->Close();
}

CSVNtupler()
{
TString a("ScaleFactor_muonEffsIsoToHLT_efficiency.txt");
output(a);
TString a("ScaleFactor_muonEffsRecoToIso_VBTF.txt");
output(a);
}
  
