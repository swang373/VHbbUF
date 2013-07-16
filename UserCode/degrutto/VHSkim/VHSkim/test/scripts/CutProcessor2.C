void Process(TString *fname)
{


//Control 1 ttbar to bb
TString *cut = new TString("Jet1.pt > 20 && Jet2.pt > 20  && Jet1.csv > 0.4 && Jet2.csv >0.4 && BestCSVPair > 0 && zmmMass > 60");


if(fname->Contains("ZllHMuv2.root"))
{
TFile *file = new TFile(fname->Data(),"read");
TTree *tree  = (TTree*) file->Get("tree");
fname->ReplaceAll("Zll","Repro-Zll");
std::cout << fname->Data();
TFile *nfile = new TFile(fname->Data(),"recreate");
TTree *ntree = tree->CopyTree(cut->Data());
nfile->Write();
std::cout << ", wrote " << ntree->GetEntries() << " events" << std::endl;

delete file;
delete nfile;
}
else
{
TFile *file = new TFile(fname->Data(),"read");
TTree *tree  = (TTree*) file->Get("tree");

fname->ReplaceAll("Zll","Repro-");
fname->ReplaceAll("Fixed-","Repro-");
std::cout << fname->Data();
TFile *nfile = new TFile(fname->Data(),"recreate");
TTree *ntree = tree->CopyTree(cut->Data());
nfile->Write();
std::cout << ", wrote "<< ntree->GetEntries() << " Events" << std::endl;
delete nfile;
delete file;
}

}




void CutProcessor2()
{
const int size =41;
string titles[size] = {
"Fixed-HZ0Jets.root",
"Fixed-HZ1Jets0-100.root",
"Fixed-HZ1Jets100-300.root",
"Fixed-HZ1Jets300-800.root",
"Fixed-HZ1Jets800-1600.root",
"Fixed-HZ2Jets0-100.root",
"Fixed-HZ2Jets100-300.root",
"Fixed-HZ2Jets300-800.root",
"Fixed-HZ2Jets800-1600.root",
"Fixed-HZ3Jets0-100.root",
"Fixed-HZ3Jets100-300.root",
"Fixed-HZ3Jets300-800.root",
"Fixed-HZ3Jets800-1600.root",
"Fixed-HZ4Jets0-100.root",
"Fixed-HZ4Jets100-300.root",
"Fixed-HZ4Jets300-800.root",
"Fixed-HZ4Jets800-1600.root",
"Fixed-HZ5Jets0-100.root",
"Fixed-HZ5Jets100-300.root",
"Fixed-HZ5Jets300-800.root",
"Fixed-HZ5Jets800-1600.root",
"ZllHZH115.root",
"ZllHZJetMad.root",
"ZllHZZ.root",
"Fixed-HZbb0Jets.root",
"Fixed-HZbb1Jets.root",
"Fixed-HZbb2Jets.root",
"Fixed-HZbb3Jets.root",
"Fixed-HZcc0Jets.root",
"Fixed-HZcc1Jets.root",
"Fixed-HZcc2Jets.root",
"Fixed-HZcc3Jets.root",
"ZllHMuv2.root",
"ZllHTTbar.root",
"ZllHWlnHbb115.root",
"ZllHWW.root",
"ZllHWZ.root",
"ZllHWjets.root",
"ZllHT-s.root",
"ZllHT-t.root",
"ZllHT-tW.root"
//"ZllHElev11v2.root",
//"ZllHZJets.root",
//"ZllHZllHbb.root",
}


for(int i=0;i<size;i++)
{
TString *fTT  = new TString(titles[i]);
Process(fTT);
}


//TString *fTT  = new TString("ZllHZJetMad.root");
//Process(fTT);

}

