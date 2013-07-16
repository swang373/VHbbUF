using namespace std;
#include <iostream>
#include "TChain.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TTree.h"
#include "TROOT.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TFile.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TLatex.h"
#include "TStyle.h"
#include <iomanip>
#include <math.h>
#include <fstream>



double deltaPhi(double phi1, double phi2) 
{ 
    double PI = 3.14159265;
    double result = phi1 - phi2;
    while (result > PI) result -= 2*PI;
    while (result <= -PI) result += 2*PI;
    return result;
}


bool first = true; 


void Process(TString fname,float weight)
{



TString bdtboost("hJet_csv[0] > 0.24 && hJet_csv[1] > 0.24 && 2+Sum$(aJet_pt >  20 && abs(aJet_eta) < 2.5) < 4 && hJet_pt[0]>20 && hJet_pt[1] >20 && Vtype == 0 && V.mass > 75 && V.mass < 105  && H.pt > 100 && V.pt > 100");

TString ttbar("MET.et > 50 && V.mass > 120");
TString zjl("Sum$(aJet_pt >  20 && abs(aJet_eta) < 2.5) < 2 && V.mass > 75 && V.mass < 105 && H.pt > 100 && V.pt > 100");
TString zjh("hJet_csv[0] > 0.90 && hJet_csv[1] > 0.90 && 2+Sum$(aJet_pt >  20 && abs(aJet_eta) < 2.5)  < 4 && (H.mass < 90 || H.mass > 145) && V.mass > 75 && V.mass < 105 && abs(deltaPhi(H.phi,V.phi)) > 2.9");



/*here*/ TString cut = bdtboost;

         cut.Append(" && Vtype == 0 && vLepton_pt[0] > 20 && vLepton_pt[1] > 20 && hJet_pt[0] > 20 && hJet_pt[1] > 20 && abs(vLepton_eta[0]) < 2.4 && abs(vLepton_eta[1]) < 2.4 && abs(hJet_eta[0]) < 2.5 && abs(hJet_eta[1]) < 2.5 && hJet_id[0]==1 && hJet_id[1]==1");

 
if(first) {std::cout << "using cuts : " << cut.Data() << std::endl<<std::endl<<std::endl; first = false;}




if(fname.Contains("Mu"))
{
TString lfname(fname.Data()); 
TString lcut(cut.Data()); 
lcut.Append("&& (triggerFlags[0] || triggerFlags[13] || triggerFlags[14]) && EVENT.json");

TFile *file = new TFile(fname.Data(),"read");
TTree *tree  = (TTree*) file->Get("tree");
lfname.ReplaceAll("Test","B-");

std::cout << lfname.Data() <<std::endl<<std::endl;
std::cout << "cuts: " << lcut.Data() <<std::endl <<std::endl; 

TFile *nfile = new TFile(lfname.Data(),"recreate");
TTree *ntree = tree->CopyTree(lcut.Data());
nfile->Write();
std::cout << "wrote "<< ntree->GetEntries() << " Events" << std::endl<<std::endl;

delete file;
delete nfile;
}

//ZJets100
else if(fname.Contains("Z-100"))
{
TString lfname(fname.Data()); 
TString lcut(cut.Data()); 
lcut.Append("&& genZpt > 120 && EVENT.event % 2 != 0");

TFile *file = new TFile(fname.Data(),"read");
TTree *tree  = (TTree*) file->Get("tree");
lfname.ReplaceAll("Test","A-");

std::cout << lfname.Data() <<std::endl<<std::endl;
std::cout << "cuts: " << lcut.Data() <<std::endl <<std::endl; 


TFile *nfile = new TFile(lfname.Data(),"recreate");
TTree *ntree = tree->CopyTree(lcut.Data());
nfile->Write();

std::cout << "wrote "<< ntree->GetEntries() << " Events" << std::endl<<std::endl;
lfname.ReplaceAll("A-","B-");
lcut.ReplaceAll("!=","==");

TFile *file2 = new TFile(fname.Data(),"read");
TTree *tree2  = (TTree*) file->Get("tree");
TFile *nfile2 = new TFile(lfname.Data(),"recreate");
TTree *ntree2 = tree->CopyTree(lcut.Data());
nfile2->Write();
}


else if(fname.Contains("DYJetsToLL"))
{
TString lfname(fname.Data()); 
TString lcut(cut.Data()); 
lcut.Append("&& genZpt < 120 && EVENT.event % 2 != 0");
TFile *file = new TFile(fname.Data(),"read");
TTree *tree  = (TTree*) file->Get("tree");

lfname.ReplaceAll("Test","A-");

std::cout << lfname.Data() <<std::endl<<std::endl;
std::cout << "cuts: " << lcut.Data() <<std::endl <<std::endl; 

TFile *nfile = new TFile(lfname.Data(),"recreate");
TTree *ntree = tree->CopyTree(lcut.Data());
nfile->Write();

std::cout << "wrote "<< ntree->GetEntries() << " Events" << std::endl<<std::endl;
lfname.ReplaceAll("A-","B-");
lcut.ReplaceAll("!=","==");

TFile *file2 = new TFile(fname.Data(),"read");
TTree *tree2  = (TTree*) file->Get("tree");
TFile *nfile2 = new TFile(lfname.Data(),"recreate");
TTree *ntree2 = tree->CopyTree(lcut.Data());
nfile2->Write();
}
else
{
TString lfname(fname.Data());
TString lcut(cut.Data()); 
lcut.Append("&& EVENT.event % 2 != 0");
TFile *file = new TFile(fname.Data(),"read");
TTree *tree  = (TTree*) file->Get("tree");


lfname.ReplaceAll("Test","A-");
std::cout << lfname.Data() <<std::endl<<std::endl;
std::cout << "cuts: " << lcut.Data() <<std::endl <<std::endl; 
TFile *nfile = new TFile(lfname.Data(),"recreate");
TTree *ntree = tree->CopyTree(lcut.Data());
nfile->Write();

std::cout << "wrote "<< ntree->GetEntries() << " Events" << std::endl<<std::endl;
lfname.ReplaceAll("A-","B-");
lcut.ReplaceAll("!=","==");
TFile *file2 = new TFile(fname.Data(),"read");
TTree *tree2  = (TTree*) file->Get("tree");
TFile *nfile2 = new TFile(lfname.Data(),"recreate");
TTree *ntree2 = tree->CopyTree(lcut.Data());
nfile2->Write();
}

}




void CutProcessorBOOST()
{
const int size =36;
string titles[size] = {
"TestZH_ZToLL_HToBB_M-115_7TeV-powheg-herwigpp.root",
"TestDYJetsToLL_PtZ-100_TuneZ2_7TeV-madgraph-tauola.root",
"TestDYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola.root",
"TestDoubleMu_Run2010_May10Rereco.root",
"TestDoubleMu_Run2011A_Aug05ReReco.root",
"TestDoubleMu_Run2011A_PromptRecoV4.root",
"TestDoubleMu_Run2011A_PromptRecoV6.root",
"TestDoubleMu_Run2011B_PromptRecoV1.root",
"TestSingleMu_Run2010_May10Rereco.root",
"TestSingleMu_Run2011A_Aug05ReReco.root",
"TestSingleMu_Run2011A_PromptRecoV4.root",
"TestSingleMu_Run2011A_PromptRecoV6.root",
"TestSingleMu_Run2011B_PromptRecoV1.root",
"TestTTJets_TuneZ2_7TeV-madgraph-tauola.root",
"TestT_TuneZ2_s-channel_7TeV-powheg-tauola.root",
"TestT_TuneZ2_t-channel_7TeV-powheg-tauola.root",
"TestT_TuneZ2_tW-channel-DR_7TeV-powheg-tauola.root",
"TestT_TuneZ2_tW-channel-DS_7TeV-powheg-tauola.root",
"TestTbar_TuneZ2_s-channel_7TeV-powheg-tauola.root",
"TestTbar_TuneZ2_t-channel_7TeV-powheg-tauola.root",
"TestTbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola.root",
"TestTbar_TuneZ2_tW-channel-DS_7TeV-powheg-tauola.root",
"TestWJetsToLNu_Pt-100_7TeV-herwigpp.root",
"TestWJetsToLNu_PtW-100_TuneZ2_7TeV-madgraph.root",
"TestWJetsToLNu_TuneZ2_7TeV-madgraph-tauola.root",
"TestWW_TuneZ2_7TeV_pythia6_tauola.root",
"TestWZ_TuneZ2_7TeV_pythia6_tauola.root",
"TestZH_ZToLL_HToBB_M-100_7TeV-powheg-herwigpp.root",
"TestZH_ZToLL_HToBB_M-105_7TeV-powheg-herwigpp.root",
"TestZH_ZToLL_HToBB_M-110_7TeV-powheg-herwigpp.root",
"TestZH_ZToLL_HToBB_M-120_7TeV-powheg-herwigpp.root",
"TestZH_ZToLL_HToBB_M-125_7TeV-powheg-herwigpp.root",
"TestZH_ZToLL_HToBB_M-130_7TeV-powheg-herwigpp.root",
"TestZH_ZToLL_HToBB_M-135_7TeV-powheg-herwigpp.root",
"TestZJetsToLL_Pt-100_7TeV-herwigpp.root",
"TestZZ_TuneZ2_7TeV_pythia6_tauola.root",
};


for(int i=0; i < size; i++)
{
TString fTT1(titles[i]);
Process(fTT1,1);

std::cout <<std::endl<<"====================================================================================================="<<std::endl;
}

}

