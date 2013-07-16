void Process(TString *fname)
{
TString *cut = new TString("H.mass<160 && H.mass> 40 && nHjj < 2 && Jet1.pt>20 && Jet2.pt>20 && Z.mass>75  && Z.mass<105 && ((Jet1.csv>0.9 && Jet2.csv>0.44) || (Jet2.csv>0.9 && Jet1.csv>0.44)) && Z.pt > 100 && deltaPullAngle > 0 && zmmMass >0");

TFile *file = new TFile(fname->Data(),"read");
TTree *tree  = (TTree*) file->Get("tree");
fname->ReplaceAll(".root","-repro.root");
TFile *nfile = new TFile(fname->Data(),"recreate");

TTree *ntree = tree->CopyTree(cut->Data());
nfile->Write();
delete file;
delete nfile;
}


void Dataizer()
{
TString *f1  = new TString("ZllHTTbar.root");
Process(f1);

TString *f2  = new TString("ZJets.root");
Process(f2);

TString *f3  = new TString("ZllHZZ.root");
Process(f3);

TString *f4 = new TString("ZllHZllHbb115.root");
Process(f4);
}

