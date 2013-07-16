

void Process(TString *fname)
{



TString *cut = new TString("Jet1.pt>20 && Jet2.pt >20 && BestCSVPair>0 && zmmMass > 75 && zmmMass < 105 && Jet1.csv > 0.5 && Jet2.csv >0.5");


TFile *file = new TFile(fname->Data(),"read");
TTree *tree  = (TTree*) file->Get("tree");

fname->ReplaceAll("42","Fast");


std::cout << fname->Data();
TFile *nfile = new TFile(fname->Data(),"recreate");
TTree *ntree = tree->CopyTree(cut->Data());
nfile->Write();
std::cout << ", wrote "<< ntree->GetEntries() << " Events" << std::endl;
delete nfile;
delete file;

}




void CutProcessorSimple()
{
const int size =3;
string titles[size] = {
"42-ZJetsMad.root",
"42-ZZ.root",
"42-ZllHbb115.root"
}


for(int i=0;i<size;i++)
{
TString *fTT  = new TString(titles[i]);
Process(fTT);
}


//TString *fTT  = new TString("ZllHZJetMad.root");
//Process(fTT);

}

