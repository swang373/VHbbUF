#include <sstream>
#include <string>
#include "stack_micro.h"
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
#include <iomanip>

using namespace std;

//#include "stackfast.h"

int Option =0;
//Options:
//0 = Data and MC Comparison
//1 = Data Only
//2 = MC Comapare No Data
//3 = Signal MC Only
//4 = 2D Significance Reporting Only!

int numCentralJets(int ni)
{
float n = (float) ni;
//return ((n+1)*n/2);
return ((1+sqrt(1+8*n))/2) ;
}



float max(float csv1, float  csv2)
{
if(csv1 > csv2) return csv1;
else return csv2;
}

float min(float csv1, float  csv2)
{
if(csv1 > csv2) return csv2;
else return csv1;
}



bool alpgenPass(float alpgendrbb, float  alpgendrcc, bool heavyFlav)
{

if ((heavyFlav) && (alpgendrbb > 0.7 || alpgendrcc > 0.7)) return true;
if ((!heavyFlav) && (alpgendrbb < 0.7 && alpgendrcc < 0.7)) return true; 
return false;
}



double thirdJetPt(double pt1, double pt2, double pt3,  double J1pt, double J2pt) 
{
 if ( pt1 == 0) return 0;
 if (( pt1 > J1pt ) || ( pt1 > J2pt))   return pt1; 
 if (( pt2 > J1pt ) || ( pt2 > J2pt))   return pt2;
 if (pt3>0) return pt3; 
 return 0;
}

void makeBDTplots(TCut XXX, TString YYY)
{
 
  TString YYYLep = YYY;

  YYYLep.Append("-BDT-zmmMass");
  makePlots("zmmMass","#mu#mu Invariant Mass [GeV]",XXX,YYYLep.Data(),0.01,5,50,200,Option,1,0,0,9999);


  YYYLep.ReplaceAll("zmmMass","zmmPt");
  makePlots("zmmPt","#mu#mu p_{T} [GeV]",XXX,YYYLep.Data(),0.01,5,0,200,Option,1,0,0,9999);

  YYYLep.ReplaceAll("zmmPt","Jet1Pt");
  makePlots("Jet1.pt", "b Jet 1 P_{T} [GeV]",XXX ,YYYLep.Data(),  0.01,  10,0,350,Option,1,0,0,9999);

  YYYLep.ReplaceAll("Jet1Pt","Jet2Pt");
  makePlots("Jet2.pt", "b Jet 2 P_{T} [GeV]",XXX ,YYYLep.Data(),  0.01,  5,0,200,Option,1,0,0,9999);

  YYYLep.ReplaceAll("Jet2Pt","Jet1Eta");
  makePlots("Jet1.eta", "b Jet 1 #eta",XXX ,YYYLep.Data(),  0.01,  0.2,-3,3,Option,1,0,0,9999);

  YYYLep.ReplaceAll("Jet1Eta","Jet2Eta");
  makePlots("Jet2.eta", "b Jet 2 #eta",XXX ,YYYLep.Data(),  0.01,  0.2,-3,3,Option,1,0,0,9999);

  YYYLep.ReplaceAll("Jet2Eta","Jet1Phi");
  makePlots("Jet1.phi", "b Jet 1 #varphi [rad]",XXX ,YYYLep.Data(),  0.01,  0.2,-4,4,Option,1,0,0,9999);

  YYYLep.ReplaceAll("Jet1Phi","Jet2Phi");
  makePlots("Jet2.phi", "b Jet 2 #varphi [rad]",XXX ,YYYLep.Data(),  0.01,  0.2,-4,4,Option,1,0,0,9999);


  YYYLep.ReplaceAll("Jet2Phi","HPt");
  makePlots("H.pt", "bb p_{T} [GeV]",XXX ,YYYLep.Data(),  0.01,  20,100,740,Option,1,0,0,9999);

  YYYLep.ReplaceAll("HPt","HMass");
  makePlots("H.mass", "bb Mass [GeV]",XXX ,YYYLep.Data(),  0.01,  10,0,600,Option,1,0,105,125);

  YYYLep.ReplaceAll("HMass","Jet1CSV");
  makePlots("Jet1.csv","b Jet 1 CSV", XXX ,YYYLep.Data(),  0.01,  0.02,0.4,1,Option,1,0,0,9999);

  YYYLep.ReplaceAll("Jet1CSV","Jet2CSV");
  makePlots("Jet2.csv", "b Jet 2 CSV",XXX ,YYYLep.Data(),  0.01,  0.02,0.4,1,Option,1,0,0,9999);

  YYYLep.ReplaceAll("Jet2CSV","bbZdphi");
  makePlots("hjjzmmdPhi", "bb-Z |#Delta#varphi| [rad]",XXX ,YYYLep.Data(),  0.01,0.2,0,3.18,Option,1,0,0,9999);

}


void makeXplots(TCut XXX, TString YYY)
{
 
  TString YYYLep = YYY;

  YYYLep.Append("-1-Leptons-Eta1");
  makePlots("mu1.eta","Muon 1 #eta",XXX,YYYLep.Data(),0.1,0.2,-3,3,Option,1,0,0,9999);

  YYYLep.ReplaceAll("Eta1","Pt1");
  makePlots("mu1.pt","Muon 1 p_{T} [GeV]",XXX,YYYLep.Data(),0.1,5,10,200,Option,1,0,0,9999);

  YYYLep.ReplaceAll("Pt1","Phi1");
  makePlots("mu1.phi","Muon 1 #varphi [rad]",XXX,YYYLep.Data(),0.1,0.1,-4,4,Option,1,0,0,9999);

  YYYLep.ReplaceAll("Phi1","Eta2");
  makePlots("mu2.eta","Muon 2 #eta",XXX,YYYLep.Data(),0.1,0.2,-3,3,Option,1,0,0,9999);

  YYYLep.ReplaceAll("Eta2","Pt2");
  makePlots("mu2.pt","Muon 2 p_{T} [GeV]",XXX,YYYLep.Data(),0.1,5,10,200,Option,1,0,0,9999);

  YYYLep.ReplaceAll("Pt2","Phi2");
  makePlots("mu2.phi","Muon 2 #varphi [rad]",XXX,YYYLep.Data(),0.1,0.1,-4,4,Option,1,0,0,9999);

  YYYLep.ReplaceAll("1-Leptons","2-Leptons");
  YYYLep.ReplaceAll("Phi2","DeltaPhi");
  makePlots("abs(deltaPhi(mu1.Phi,mu2.phi))"," #mu-#mu |#Delta#varphi| [rad]",XXX,YYYLep.Data(),0.1,0.1,-0.1,3.2,Option,1,0,0,9999);


  YYYLep.ReplaceAll("DeltaPhi","DeltaPt");
  makePlots("mu1.pt - mu2.pt","#mu-#mu #Delta p_{T} [GeV]",XXX,YYYLep.Data(),0.01,2,-4,200,Option,1,0,0,9999);
  
  YYYLep.ReplaceAll("DeltaPt","DeltaEta");
  makePlots("abs(mu1.eta - mu2.eta)","#mu #mu |#Delta#eta|",XXX,YYYLep.Data(),0.01,0.1,0,5,Option,1,0,0,9999);


  YYYLep.ReplaceAll("DeltaEta","zmmMass");
//  makePlots("zmmMass","#mu#mu Invariant Mass [GeV]",XXX,YYYLep.Data(),0.01,5,0,200,Option,1,0,0,9999);
  makePlots("zmmMass","#mu#mu Invariant Mass [GeV]",XXX,YYYLep.Data(),0.01,1,75,105,Option,1,0,0,9999);


  YYYLep.ReplaceAll("zmmMass","zmmPt");
  makePlots("zmmPt","#mu#mu p_{T} [GeV]",XXX,YYYLep.Data(),0.01,5,0,200,Option,1,0,0,9999);

  YYYLep.ReplaceAll("2-Leptons","3-Jets");
  YYYLep.ReplaceAll("zmmPt","Jet1Pt");
  makePlots("Jet1.pt", "b Jet 1 P_{T} [GeV]",XXX ,YYYLep.Data(),  0.01,  10,0,350,Option,1,0,0,9999);

  YYYLep.ReplaceAll("Jet1Pt","Jet2Pt");
  makePlots("Jet2.pt", "b Jet 2 P_{T} [GeV]",XXX ,YYYLep.Data(),  0.01,  5,0,200,Option,1,0,0,9999);

  YYYLep.ReplaceAll("Jet2Pt","Jet1Eta");
  makePlots("Jet1.eta", "b Jet 1 #eta",XXX ,YYYLep.Data(),  0.01,  0.2,-3,3,Option,1,0,0,9999);

  YYYLep.ReplaceAll("Jet1Eta","Jet2Eta");
  makePlots("Jet2.eta", "b Jet 2 #eta",XXX ,YYYLep.Data(),  0.01,  0.2,-3,3,Option,1,0,0,9999);

  YYYLep.ReplaceAll("Jet2Eta","Jet1Phi");
  makePlots("Jet1.phi", "b Jet 1 #varphi [rad]",XXX ,YYYLep.Data(),  0.01,  0.2,-4,4,Option,1,0,0,9999);

  YYYLep.ReplaceAll("Jet1Phi","Jet2Phi");
  makePlots("Jet2.phi", "b Jet 2 #varphi [rad]",XXX ,YYYLep.Data(),  0.01,  0.2,-4,4,Option,1,0,0,9999);

  YYYLep.ReplaceAll("3-Jets","4-Jets");
  YYYLep.ReplaceAll("Jet2Phi","NumJets");
  makePlots("numJets","Number of Jets",XXX,YYYLep.Data(),0.01,1,0,20,Option,1,0,0,9999);
  

  YYYLep.ReplaceAll("NumJets","HPt");
  makePlots("H.pt", "bb p_{T} [GeV]",XXX ,YYYLep.Data(),  0.01,  20,100,740,Option,1,0,0,9999);

  YYYLep.ReplaceAll("HPt","HMass");
  makePlots("H.mass", "bb Mass [GeV]",XXX ,YYYLep.Data(),  0.01,  10,0,600,Option,1,0,105,125);

  YYYLep.ReplaceAll("HMass","bbDphi");
  makePlots("abs(deltaPhi(Jet1.phi,Jet2.phi))", "b-b |#Delta#varphi| [rad]",XXX ,YYYLep.Data(),  0.01,0.1,0,3.5,Option,1,0,0,9999);

  YYYLep.ReplaceAll("bbDphi","bbDeta");
  makePlots("abs(Jet1.eta-Jet2.eta)", "b-b |#Delta#eta|",XXX ,YYYLep.Data(),  0.01,0.1,0,5,Option,1,0,0,9999);

  YYYLep.ReplaceAll("bbDeta","bbDpt");
  makePlots("Jet1.pt-Jet2.pt", "b-b #Deltap_{T} [GeV]",XXX ,YYYLep.Data(),  0.01,5,-1,200,Option,1,0,0,9999);

  YYYLep.ReplaceAll("bbDpt","maxCSV");
//  makePlots("max(Jet1.csv,Jet2.csv)","b Jet Max CSV", XXX ,YYYLep.Data(),  0.01,  0.02,0,1,Option,1,0,0,9999);
  makePlots("max(Jet1.csv,Jet2.csv)","b Jet Max CSV", XXX ,YYYLep.Data(),  0.01,  0.01,0.88,1,Option,1,0,0,9999);

  YYYLep.ReplaceAll("maxCSV","minCSV");
//  makePlots("min(Jet1.csv,Jet2.csv)", "b Jet Min CSV",XXX ,YYYLep.Data(),  0.01,  0.02,0,1,Option,1,0,0,9999);
  makePlots("min(Jet1.csv,Jet2.csv)", "b Jet Min CSV",XXX ,YYYLep.Data(),  0.01,  0.02,0.3,1,Option,1,0,0,9999);

  YYYLep.ReplaceAll("4-Jets","5-Jets");
  YYYLep.ReplaceAll("minCSV","bbZdphi");
//  makePlots("hjjzmmdPhi", "bb-Z |#Delta#varphi| [rad]",XXX ,YYYLep.Data(),  0.01,0.1,0,3.5,Option,1,0,0,9999);
  makePlots("hjjzmmdPhi", "bb-Z |#Delta#varphi| [rad]",XXX ,YYYLep.Data(),  0.01,0.02,2.76,3.18,Option,1,0,0,9999);

  YYYLep.ReplaceAll("bbZdphi","bbzmmPt");
  makePlots("H.pt-Z.Pt", "bb-Z #Deltap_{T} [GeV]",XXX ,YYYLep.Data(),  0.01,10,-1,400,Option,1,0,0,9999);

  YYYLep.ReplaceAll("bbzmmPt","TotMass");
  makePlots("Mass(H.pt, H.eta,  H.phi, H.mass, Z.Pt,  Z.eta,  Z.Phi,  Z.Mass)", "bb+Z Mass_{inv} [GeV]",XXX ,YYYLep.Data(),  0.001,50,-1,1500,Option,1,0,0,9999);

  YYYLep.ReplaceAll("TotMass","TotPt");
  makePlots("PT(H.pt, H.eta,  H.phi, H.mass, Z.Pt,  Z.eta,  Z.Phi,  Z.Mass)", "bb+Z p_{T} [GeV]",XXX ,YYYLep.Data(),  0.01,10,-1,500,Option,1,0,0,9999);

  YYYLep.ReplaceAll("TotPt","JZB");
  makePlots("JZB","JZB",XXX ,YYYLep.Data(),  0.01,5,0,200,Option,1,0,0,9999);

  YYYLep.ReplaceAll("JZB","Met");
  makePlots("MET.et","MET Et [GeV]",XXX ,YYYLep.Data(),  0.01,5,0,100,Option,1,0,0,9999);

}



void micro_Zll() 
{
  
           //Usage
           //Variable  |  title     | cut     |  fileName | ymin | bin separation |  xmin |  xmax | doData |  doLog |  compareShape (Norm to Data area)



TCut ttbar("ttbar","numCentralJets(nHjj) > 1 &&  Jet1.pt > 20 && Jet2.pt > 20  && MET.et > 50 && BestCSVPair > 0 &&  zmmMass > 120");


TString s0("0-MCCOMP-");
TString s1("1-MCCOMP-");
TString s2("2-MCCOMP-");
TString s3("3-MCCOMP-");
TString s4("4-MCCOMP-");
TString s5("5-MCCOMP-");
TString stt("tt-MCCOMP-");

TString sC1("C1-MCCOMP-");


TCut cut0("cut0","Jet1.pt>20 && Jet2.pt >20 && BestCSVPair>0");
TCut cut1("cut1","Jet1.pt>20 && Jet2.pt >20 && BestCSVPair>0 && zmmMass > 75 && zmmMass < 105");
TCut cut2("cut2","Jet1.pt>20 && Jet2.pt >20 && BestCSVPair>0 && zmmMass > 75 && zmmMass < 105 &&  ((Jet1.csv > 0.9 && Jet2.csv >0.5) || (Jet2.csv > 0.9 && Jet1.csv >0.5))");
TCut cut3("cut3","Jet1.pt>20 && Jet2.pt >20 && BestCSVPair>0 && zmmMass > 75 && zmmMass < 105 &&  ((Jet1.csv > 0.9 && Jet2.csv >0.5) || (Jet2.csv > 0.9 && Jet1.csv >0.5)) && hjjzmmdPhi > 2.9");
TCut cut4("cut4","H.pt > 100 && zmmPt > 100 && Jet1.pt>20 && Jet2.pt >20 && BestCSVPair>0 && zmmMass > 75 && zmmMass < 105 &&  ((Jet1.csv > 0.9 && Jet2.csv >0.5) || (Jet2.csv > 0.9 && Jet1.csv >0.5)) && hjjzmmdPhi > 2.9");

TCut cut5("cut5","numCentralJets(nHjj) < 4 && zmmPt > 100 && H.pt > 100 && Jet1.pt>20 && Jet2.pt >20 && BestCSVPair>0 && zmmMass > 75 && zmmMass < 105 &&  ((Jet1.csv > 0.9 && Jet2.csv >0.5) || (Jet2.csv > 0.9 && Jet1.csv >0.5)) && hjjzmmdPhi > 2.9");

//TCut mass("H.mass > 88 && H.mass < 122");
//TCut amass("H.mass < 90");

//TCut relax4("numJets < 4  && BestCSVPair>0 && Jet1.csv > 0.24 && Jet2.csv > 0.24 && zmmMass > 75 && zmmMass < 105 && hjjzmmdPhi > 2.0");
//TCut relax3("numJets < 4  && BestCSVPair>0 &&  zmmMass > 75 && zmmMass < 105 && hjjzmmdPhi > 2.0 ");
//makeXplots(cut0,s0);
//makeXplots(cut1,s1);
//makeXplots(cut2,s2);
//makeXplots(cut3,s3);
//makeXplots(zjl,stt);

TCut zjl("zjl","zmmPt > 100 && Jet1.csv < 0.24 && Jet2.csv <0.24 && Jet1.pt > 20 && Jet2.pt > 20 && numCentralJets(nHjj) < 4 && BestCSVPair > 0 &&  zmmMass > 75 && zmmMass < 105 && H.pt > 100");
TCut zjh("zjh","numCentralJets(nHjj)  < 4 && (H.mass < 95 || H.mass > 145) && Jet1.pt > 20 && Jet2.pt > 20  && Jet1.csv > 0.898 && Jet2.csv >0.898 && BestCSVPair > 0 &&  zmmMass > 75 && zmmMass < 105 && hjjzmmdPhi > 2.9");
//makeXplots(cut4+amass,sC1);

makePlots("H.pt", "bb p_{T} [GeV]",zjh,"zjlpt",  0.01,25,0,300,Option,1,0,0,9999);

makePlots("H.mass", "bb Mass [GeV]",zjh,"zjlmass",  0.01,25,0,500,Option,1,0,0,9999);
//makePlots("H.mass", "bb Mass [GeV]",cut5+"H.mass > 95 && H.mass < 125","cutpumass",  0.01,1,95,125,Option,1,0,0,9999);
//makePlots("H.mass", "bb Mass [GeV]",zjh,"zjhpumass",  0.01,25,0,500,Option,1,0,0,9999);
//makePlots("numCentralJets(nHjj)", "Num. Central Jets",cut4,"njets",  0.01,1,0,12,Option,1,0,0,9999);
//makePlots("addJet1.pt", "Third Jet 1 p_{T}",cut4,"3jet",  0.01,10,0,200,Option,1,0,0,9999);
//makePlots("addJet2.pt", "Fourth Jet 1 p_{T}",cut4,"4jet",  0.01,10,0,200,Option,1,0,0,9999);
//makePlots("Jet1.pt", "Jet 1 p_{T} [GeV]",cut1,"cut1pt1",  0.01,25,0,300,Option,1,0,0,9999);
//makePlots("Jet2.pt", "Jet 2 p_{T} [GeV]",cut1,"cut1pt2",  0.01,25,0,300,Option,1,0,0,9999);
//makePlots("numJets", "Number of Jets",cut1,"cut1njets",  0.01,1,0,12,Option,1,0,0,9999);
//makePlots("Jet1.pt", "Jet 1 p_{T} (ee channel of Z)","Jet1.pt > 20 && Jet2.pt > 20 && BestCSVPair >0 && zmmMass > 60","Jet1Ptee",  0.01,10,0,300,Option,1,0,0,9999);

//makePlots("numJets","Number of Jets",cut3,"numJets",0.01,1,0,20,Option,1,0,0,9999);
//makePlots("ThirdJetPt", "Third Jet p_{T} [GeV]",cut1,"thirdjet",  0.01,10,-20,300,Option,1,0,0,9999);

//makePlots("zmmMass","#mu#mu Invariant Mass [GeV]","zmmMass > 50 && Jet1.pt > 20 && Jet2.pt > 20 && mu1.pt > 20 && mu2.pt > 20 && BestCSVPair > 0","zMass",0.01,5,0,200,Option,1,0,0,9999);
//makePlots("H.mass", "bb Mass [GeV]",zjh,"cut3bb",  0.01,10,0,450,Option,1,0,0,9999);
//makePlots("numJets","Number of Jets",zjh,"numJets",0.01,1,0,20,Option,1,0,0,9999);


/*
makePlots("H.mass", "bb Mass [GeV]",cut4,"cut4bb",  0.01,1,80,150,Option,1,0,0,9999);
makePlots("H.mass", "bb Mass [GeV]",cut4+"H.mass > 95 && H.mass < 125","cut4bb120",  0.01,1,80,150,Option,1,0,0,9999);
makePlots("H.mass", "bb Mass [GeV]",cut4+"H.mass > 100 && H.mass < 130","cut4bb130",  0.01,1,80,150,Option,1,0,0,9999);
makePlots("H.mass", "bb Mass [GeV]",cut4+"H.mass > 105 && H.mass < 135","cut4bb135",  0.01,1,80,150,Option,1,0,0,9999);
makePlots("H.mass", "bb Mass [GeV]",cut4+"H.mass > 110 && H.mass < 140","cut4bb140",  0.01,1,80,150,Option,1,0,0,9999);
*/

//  makePlots("Jet1.csv","b Jet Max CSV", relax3 ,"jet1csvlog",  0.01,  0.02,0,1,Option,1,0,0,9999);
//makePlots("hjjzmmdPhi", "H-Z |#Delta#varphi| [rad]",cut2,"deltaPhi",  0.01,0.1,0,3.5,Option,1,0,0,9999);
// makePlots("Jet1.csv","b Jet Max CSV", relax4 ,"jet1csv",  0.01,  0.02,0,1,Option,0,0,0,9999);
//TString *relax4 = new TString("numJets < 4  && H.pt > 50 && zmmPt > 50  && BestCSVPair>0 && Jet1.csv >0.5 &&  Jet2.csv > 0.5 && zmmMass > 50");


// /*here*/ makePlots("H.mass", "bb Mass [GeV]",cutnb,"cut4bbMassCorrected",  0.01,10,0,300,Option,1,0,90,120);
// /*here*/ makePlots("H.pt", "bb p_{T} [GeV]",cut4,"cut4bbPTCorrected",  0.01,10,10,300,Option,1,0,90,120);

// /*here*/ makePlots("H.pt", "bb p_{T} [GeV]",CCC,"CCCbbptlin",  0.01,10,40,200,Option,0,0,90,120);
// /*here*/ makePlots("numJets", "Number of Jets",CCC,"CCCnjetslin",  0.01,1,0,14,Option,0,0,0,9999);


//  makePlots("H.mass", "bb Mass [GeV]",boost70,"temp",  0.01,10,0,300,Option,1,0,0,9999);
/*
std::cout << "void SignalFOM60(TH2F *SIGEE) {" <<std::endl;


for(int i=0;i<20;i++)
 for(int j=0;j<20;j++)
 {
  int  ii=i*10;
  int  jj=j*10;

  TString cutsee = TString::Format("Jet1.pt>60 && Jet2.pt >20 && BestCSVPair>0 && zeeMass > 75 && zeeMass < 105 &&  ((Jet1.csv > 0.9 && Jet2.csv >0.5) || (Jet2.csv > 0.9 && Jet1.csv >0.5)) && hjjzeedPhi > 2.8 && numJets < 4 && H.mass > 90 && H.mass < 120 &&  zeePt > %i && H.pt > %i", ii, jj);

//  std::cout << "SIGMM->SetBinContent(" << i+1 << "," << j+1 << ",";
//  makePlots("H.mass", "bb Mass [GeV]",cutsmm.Data(),"temp",  0.01,10,0,300,Option,1,0,0,9999);
//  std::cout << ");" <<std::endl;

  std::cout << "SIGEE->SetBinContent(" << i+1 << "," << j+1 << ",";
  makePlots("H.mass", "bb Mass [GeV]",cutsee.Data(),"temp",  0.01,10,0,300,Option,1,0,0,9999);
  std::cout << ");" <<std::endl;



 }

std::cout << "}" <<std::endl;
//  */
//makePlots("H.mass", "HMASS",cut3+"numJets < 3","HMASS",  0.01,10,0,100,Option,1,0,0,9999);
//makePlots("hjjzmmdPhi", "bb-Z |#Delta#varphi| [rad]",cut2+"MET.et < 50","dPhi",  0.01,0.1,0,3.5,Option,1,0,0,9999);
//makePlots("MET.et", "MET E_{T} [GeV]",cut2+"hjjzmmdPhi > 2.8","MET",  0.01,0.1,0,3.5,Option,1,0,0,9999);
}
