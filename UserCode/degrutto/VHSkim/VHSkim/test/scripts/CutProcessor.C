bool first = true; 
int background = 0;
int numCentralJets(int ni)
{
float n = (float) ni;
//return ((n+1)*n/2);
return ((1+sqrt(1+8*n))/2) ;
}



void Process(TString *fname,float weight)
{


//Control 1 zjl to bb
TString *bdt = new TString("numCentralJets(nHjj) < 4 && Jet1.pt>20 && Jet2.pt >20 && BestCSVPair>0 && zmmMass > 75 && zmmMass < 105 &&  Jet1.csv > 0.5 && Jet2.csv > 0.5  && hjjzmmdPhi > 2.4");
TString *bdta = new TString("numCentralJets(nHjj) < 4 && Jet1.pt>20 && Jet2.pt >20 && BestCSVPair>0 && zmmMass > 75 && zmmMass < 105 &&  Jet1.csv > 0.5 && Jet2.csv > 0.5  && hjjzmmdPhi > 2.4 && eventInfo.event % 2 == 0");
TString *bdtb = new TString("numCentralJets(nHjj) < 4 && Jet1.pt>20 && Jet2.pt >20 && BestCSVPair>0 && zmmMass > 75 && zmmMass < 105 &&  Jet1.csv > 0.5 && Jet2.csv > 0.5  && hjjzmmdPhi > 2.4 && eventInfo.event % 2 != 0");




TCut cut5("numCentralJets(nHjj) < 4 && zmmPt > 100 && H.pt > 100 && Jet1.pt>20 && Jet2.pt >20 && BestCSVPair>0 && zmmMass > 75 && zmmMass < 105 &&  ((Jet1.csv > 0.9 && Jet2.csv >0.5) || (Jet2.csv > 0.9 && Jet1.csv >0.5)) && hjjzmmdPhi > 2.9");
TCut cut5a("numCentralJets(nHjj) < 4 && zmmPt > 100 && H.pt > 100 && Jet1.pt>20 && Jet2.pt >20 && BestCSVPair>0 && zmmMass > 75 && zmmMass < 105 &&  ((Jet1.csv > 0.9 && Jet2.csv >0.5) || (Jet2.csv > 0.9 && Jet1.csv >0.5)) && hjjzmmdPhi > 2.9 &&  eventInfo.event % 2 == 0");
TCut cut5b("numCentralJets(nHjj) < 4 && zmmPt > 100 && H.pt > 100 && Jet1.pt>20 && Jet2.pt >20 && BestCSVPair>0 && zmmMass > 75 && zmmMass < 105 &&  ((Jet1.csv > 0.9 && Jet2.csv >0.5) || (Jet2.csv > 0.9 && Jet1.csv >0.5)) && hjjzmmdPhi > 2.9 && eventInfo.event % 2 != 0");



TString *zjh =new TString("numJets < 4 && (H.mass < 100 || H.mass > 140) && Jet1.pt > 20 && Jet2.pt > 20  && ((Jet1.csv > 0.9 && Jet2.csv >0.5) || (Jet2.csv > 0.9 && Jet1.csv >0.5)) && BestCSVPair > 0 &&  zmmMass > 75 && zmmMass < 105");
TString *zjha =new TString("numJets < 4 && (H.mass < 100 || H.mass > 140) && Jet1.pt > 20 && Jet2.pt > 20  && ((Jet1.csv > 0.9 && Jet2.csv >0.5) || (Jet2.csv > 0.9 && Jet1.csv >0.5)) && BestCSVPair > 0 &&  zmmMass > 75 && zmmMass < 105 && eventInfo.event % 2 == 0");
TString *zjhb =new TString("numJets < 4 && (H.mass < 100 || H.mass > 140) && Jet1.pt > 20 && Jet2.pt > 20  && ((Jet1.csv > 0.9 && Jet2.csv >0.5) || (Jet2.csv > 0.9 && Jet1.csv >0.5)) && BestCSVPair > 0 &&  zmmMass > 75 && zmmMass < 105 && eventInfo.event % 2 != 0 ");

//&& H.pt > 50 && zmmPt > 50 
TString *ttbar = new TString("numJets < 4 &&  H.pt > 50 && zmmPt > 50 && Jet1.pt > 20 && Jet2.pt > 20  && Jet1.csv > 0.24 && Jet2.csv >0.24 && BestCSVPair > 0 &&  zmmMass > 120");
TString *ttbara = new TString("numJets < 4 && H.pt > 50 && zmmPt > 50 && Jet1.pt > 20 && Jet2.pt > 20  && Jet1.csv > 0.24 && Jet2.csv >0.24 && BestCSVPair > 0  && eventInfo.event % 2 == 0 &&  zmmMass > 120");
TString *ttbarb = new TString("numJets < 4 && H.pt > 50 && zmmPt > 50 && Jet1.pt > 20 && Jet2.pt > 20  && Jet1.csv > 0.24 && Jet2.csv >0.24 && BestCSVPair > 0  && eventInfo.event % 2 != 0 &&  zmmMass > 120");



TString *zjl = new TString("numJets < 4 && H.pt > 50 && zmmPt > 50 && Jet1.pt > 20 && Jet2.pt > 20  && Jet1.csv < 0.24 && Jet2.csv < 0.24 && BestCSVPair > 0 &&  zmmMass > 75 && zmmMass < 105");
TString *zjla = new TString("numJets < 4 && H.pt > 50 && zmmPt > 50 && Jet1.pt > 20 && Jet2.pt > 20  && Jet1.csv < 0.24 && Jet2.csv < 0.24 && BestCSVPair > 0  && eventInfo.event % 2 == 0 &&  zmmMass > 75 && zmmMass < 105");
TString *zjlb = new TString("numJets < 4 && H.pt > 50 && zmmPt > 50 && Jet1.pt > 20 && Jet2.pt > 20  && Jet1.csv < 0.24 && Jet2.csv <0.24 && BestCSVPair > 0  && eventInfo.event % 2 != 0 &&  zmmMass > 75 && zmmMass < 105");

TString *cut4 = new TString("numJets < 4 && H.pt > 50 && zmmPt > 50 && Jet1.pt>20 && Jet2.pt >20 && BestCSVPair>0 && zmmMass > 75 && zmmMass < 999 &&  ((Jet1.csv > 0.9 && Jet2.csv >0.5) || (Jet2.csv > 0.9 && Jet1.csv >0.5)) && hjjzmmdPhi > 2.8");
TString *cut4a = new TString("numJets < 4 && H.pt > 50 && zmmPt > 50 && Jet1.pt>20 && Jet2.pt >20 && BestCSVPair>0 && zmmMass > 75 && zmmMass < 999 &&  ((Jet1.csv > 0.9 && Jet2.csv >0.5) || (Jet2.csv > 0.9 && Jet1.csv >0.5)) && hjjzmmdPhi > 2.8 && numJets < 4 && eventInfo.event % 2 == 0");
TString *cut4b = new TString("numJets < 4 && H.pt > 50 && zmmPt > 50 && Jet1.pt>20 && Jet2.pt >20 && BestCSVPair>0 && zmmMass > 75 && zmmMass < 999 &&  ((Jet1.csv > 0.9 && Jet2.csv >0.5) || (Jet2.csv > 0.9 && Jet1.csv >0.5)) && hjjzmmdPhi > 2.8 && numJets < 4 && eventInfo.event % 2 != 0");

TString *relax4 =  new TString("numJets < 4  && BestCSVPair>0 &&  Jet2.csv > 0.24 && Jet1.csv >0.24 && zmmMass > 75 && zmmMass < 105 && hjjzmmdPhi > 2.0");
TString *relax4a = new TString("numJets < 4  && BestCSVPair>0 &&  Jet2.csv > 0.24 && Jet1.csv >0.24 && zmmMass > 75 && zmmMass < 105 && hjjzmmdPhi > 2.0 && eventInfo.event % 2 != 0");
TString *relax4b = new TString("numJets < 4  && BestCSVPair>0 &&  Jet2.csv > 0.24 && Jet1.csv >0.24 && zmmMass > 75 && zmmMass < 105 && hjjzmmdPhi > 2.0 && eventInfo.event % 2 == 0");





/*here*/ if(first) {std::cout << "using cuts : bdt" << std::endl; first = false;}
if(fname->Contains("1fb"))
{
TFile *file = new TFile(fname->Data(),"read");
TTree *tree  = (TTree*) file->Get("tree");
fname->ReplaceAll("Zll","B-Zll");
fname->ReplaceAll("42-","B-42-");
std::cout << fname->Data();
TFile *nfile = new TFile(fname->Data(),"recreate");
/*here*/TTree *ntree = tree->CopyTree(bdt->Data());
nfile->Write();
std::cout << ", wrote "<< ntree->GetEntries() << " Events" << std::endl;

delete file;
delete nfile;
}
else
{
TFile *file = new TFile(fname->Data(),"read");
TTree *tree  = (TTree*) file->Get("tree");

fname->ReplaceAll("42-","A-42-");
fname->ReplaceAll("Zll","A-42-");
fname->ReplaceAll("Fixed-","A-Zll");
std::cout << fname->Data();
TFile *nfile = new TFile(fname->Data(),"recreate");
/*here*/TTree *ntree = tree->CopyTree(bdta->Data());
nfile->Write();
if(fname->Contains("H115"))
{
  ofstream myfile;
  myfile.open ("signal");
  myfile << ntree->GetEntries()-10;
  myfile.close();
}
 else  if(ntree->GetEntries() >  0) background += ntree->GetEntries();
std::cout << ", wrote "<< ntree->GetEntries() << " Events" << std::endl;
delete nfile;

fname->ReplaceAll("A-42-","B-42-");
std::cout << fname->Data();
TFile *nfile2 = new TFile(fname->Data(),"recreate");
/*here*/TTree *ntree2 = tree->CopyTree(bdtb->Data());
std::cout << ", wrote "<< ntree2->GetEntries() << " Events" << std::endl;
nfile2->Write();
delete nfile2;
delete file;
}

}




void CutProcessor()
{
const int size =1;
string titles[size] = {
//"ZllHZH115.root",
"42-SingleMu1fb.root"
//"42-ZJetsMad.root",
//"42-TTbar.root",
//"42-WW.root",
//"42-WZ.root",
//"42-WJetsMad.root",
//"ZllHT-s.root",
//"ZllHT-tW.root",
//"42-ZZ.root"
//"ZllHT-tW.root"
//"ZllHWlnHbb115.root",
//"ZllHElev11v2.root",
//"ZllHZJets.root",
//"ZllHZllHbb.root",
}


float weight[11] ={
7638,
357,
0.186,
11.51,
48.0,
42,
0.463,
499.0,
23.0,
46.0,
17.14
};


for(int i=0;i<size;i++)
{

TString *fTT  = new TString(titles[i]);
Process(fTT,0.186/weight[i]);
}
  ofstream myfile;
  myfile.open ("background");
  myfile << background-10;
  myfile.close();

//TString *fTT  = new TString("ZllHZJetMad.root");
//Process(fTT);

}

