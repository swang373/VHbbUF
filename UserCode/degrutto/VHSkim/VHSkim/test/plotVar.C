#include <iostream>
#include "plot.h"

double deltaPhi(double phi1, double phi2) 
{ 
    double PI = 3.14159265;
    double result = phi1 - phi2;
    while (result > PI) result -= 2*PI;
    while (result <= -PI) result += 2*PI;
    return result;
}

void plotVar()
{
TString *cuts = new TString("hjjJet1Pt > 30 && hjjJet2Pt > 30 && ((hjjJet1CSV>0.9 && hjjJet2CSV>0.44) || (hjjJet2CSV>0.9 && hjjJet1CSV>0.44)) && hjjPt > 0 && hjjPt < 150");	
//TString *cuts = new TString("hjjJet1Pt > 20 && hjjJet2Pt > 20 && ((hjjJet1CSV>0.9 && hjjJet2CSV>0.44) || (hjjJet2CSV>0.9 && hjjJet1CSV>0.44)) && hjjPt > 0 && hjjPt < 100 && zmmMass>75 && zmmMass<105");	
//TString cuts = new TString("hjjJet1Pt > 20 && hjjJet2Pt > 20 && (hjjPt > 100 && zmmMass>75 && zmmMass<105");	






int nbin = 20;
int xmin = 50;
int xmax = 150;

std::string var,name,axis,chan;
float max;
var = "?";


std::cout << "Which Channel? (ee or mm)" << std::endl;
std::cin >> chan;
std::cin.ignore();

std::cout << "Which Variable Would You Like To Plot? (type: ? to see a list)" << std::endl;
std::cin >> var;

if(var=="?")
 {outputVars(); std::cout << std::endl <<std::endl; std::cout << "Which Variable Would You Like To Plot?" << std::endl;std::cin.ignore();std::cin >> var;}

std::cout << "Name of Plot?" << std::endl;
std::cin.ignore();
std::getline(std::cin, name);
std::cout << "X Axis Title?" << std::endl;
std::getline(std::cin, axis);
std::cout << "X Axis Min" << std::endl;
std::cin >> xmin;
std::cout << "X Axis Max" << std::endl;
std::cin >> xmax;
std::cout << "X Axis Number of Bins" << std::endl;
std::cin >> nbin;



TString *V = new TString(var);
TString *N = new TString(name);
TString *A = new TString(axis);
TString *C = new TString(chan);

double lumi;

if(chan == "mm" || chan == "MM")
{
 cuts->ReplaceAll("zee","zmm");	
 std::cout << cuts->Data() <<std::endl;
 lumi = 0.124;
}
if(chan == "ee" || chan == "EE")
{
 cuts->ReplaceAll("zmm","zee");	
 std::cout << cuts->Data() <<std::endl;
 lumi =0.138;
}


gROOT->SetStyle("Plain");

const double lumiZJetsMad      = 0.744;
const double lumiZ1jets0_100   = 1.6;
const double lumiZ1jets100_300 = 18.5;
const double lumiZ1jets300_800 = 1619.0;
const double lumiZ1jets800_1600= 185930;
const double lumiZ2jets0_100   = 1.93;
const double lumiZ2jets100_300 = 2.3;
const double lumiZ2jets300_800 = 293;
const double lumiZ2jets800_1600= 42403;
const double lumiZbb0jets      = 29.8;
const double lumiZbb1jets      = 61;
const double lumiZbb2jets      = 17.4;
const double lumiZcc0jets      = 50;
const double lumiZcc1jets      = 57;
const double lumiZcc2jets      = 8;
const double lumiTTbar =  5.0    ; 
const double lumiZZ =     357.0  ; 
const double lumiZH =     870.0  ; 

TH1F *ZJets = new TH1F("ZJets","ZJets",nbin,xmin,xmax);

/*TH1F *ZJ2 = new TH1F("ZJ2","ZJ2",nbin,xmin,xmax);
TH1F *ZJ3 = new TH1F("ZJ3","ZJ3",nbin,xmin,xmax);
TH1F *ZJ4 = new TH1F("ZJ4","ZJ4",nbin,xmin,xmax);
TH1F *ZJ5 = new TH1F("ZJ5","ZJ5",nbin,xmin,xmax);
TH1F *ZJ6 = new TH1F("ZJ6","ZJ6",nbin,xmin,xmax);
TH1F *ZJ7 = new TH1F("ZJ7","ZJ7",nbin,xmin,xmax);
TH1F *ZJ8 = new TH1F("ZJ8","ZJ8",nbin,xmin,xmax);*/
TH1F *ZZ = new TH1F("ZZ","ZZ",nbin,xmin,xmax);
TH1F *ZH = new TH1F("ZH","ZH",nbin,xmin,xmax);
TH1F *TTbar = new TH1F("TTbar","TTbar",nbin,xmin,xmax);
TH1F *Muv1 = new TH1F("Muv1","Muv1",nbin,xmin,xmax);
TH1F *Muv2 = new TH1F("Muv2","Muv2",nbin,xmin,xmax);
TH1F *Ev1 = new TH1F("Ev1","Ev1",nbin,xmin,xmax);
TH1F *Ev2 = new TH1F("Ev2","Ev2",nbin,xmin,xmax);


//Define what Plot You Want


TChain *chainTT = new TChain("Events");
chainTT->Add("/data/1a/madfish/JuneData/Jun7/TTbar.root");
chainTT ->SetWeight(lumi/lumiTTbar,"global");
chainTT->Project("TTbar",V->Data(),cuts->Data());

//chainTT->Draw(V->Data(),"hjjMass<102 && hjjMass> 126 && hjjMass@.size() < 2 && hjjJet1Pt>20 && hjjJet2Pt>20 && zmmMass>75  && zmmMass<105 && ((hjjJet1CSV>0.9 && hjjJet2CSV>0.44) || (hjjJet2CSV>0.9 && hjjJet1CSV>0.44)) && zmmPt > 150 && deltaPhi(hjjPhi,zmmPhi) > 2.75 && hjjPt > 150");




TChain *chainZH = new TChain("Events");
chainZH ->SetWeight(lumi/lumiZH,"global");
chainZH->Add("/data/1a/degrutto/data/Spring11/ZllHbb/ZllHbb/EdmNtuples_Zll_2_1_8l0.root");
chainZH->Project("ZH",V->Data(),cuts->Data());


TChain *chainZZ = new TChain("Events");
chainZZ->SetWeight(lumi/lumiZZ,"global");
chainZZ->Add("/data/1a/madfish/JuneData/Jun7/ZZ.root");
chainZZ->Project("ZZ",V->Data(),cuts->Data());

TChain *chainZJM = new TChain("Events");
chainZJM->SetWeight(lumi/lumiZJetsMad ,"global");
chainZJM->Add("/data/1a/madfish/JuneData/Jun7/ZJetMad.root");
chainZJM->Project("ZJets",V->Data(),cuts->Data());


/*
TChain *chainZJ0 = new TChain("Events");
chainZJ0->SetWeight(lumi/lumiZ1jets0_100 ,"global");
chainZJ0->Add("/data/1a/madfish/JuneData/Z1Jets_ptZ_0to100.root");
chainZJ0->Project("ZJets",V->Data(),cuts->Data());

TChain *chainZJ100 = new TChain("Events");
chainZJ100->SetWeight(lumi/lumiZ1jets100_300 ,"global");
chainZJ100->Add("/data/1a/madfish/JuneData/Z1Jets_ptZ_100to300.root");
chainZJ100->Project("ZJ2",V->Data(),cuts->Data());


TChain *chainZJ300 = new TChain("Events");
chainZJ300->SetWeight(lumi/lumiZ1jets300_800 ,"global");
chainZJ300->Add("/data/1a/madfish/JuneData/Z1Jets_ptZ_300to800.root");
chainZJ300->Project("ZJ3",V->Data(),cuts->Data());

TChain *chainZJ800 = new TChain("Events");
chainZJ800->SetWeight(lumi/lumiZ1jets800_1600 ,"global");
chainZJ800->Add("/data/1a/madfish/JuneData/Z1Jets_ptZ_800to1600.root");
chainZJ800->Project("ZJ4",V->Data(),cuts->Data());

TChain *chainZ2J0 = new TChain("Events");
chainZ2J0->SetWeight(lumi/lumiZ2jets0_100 ,"global");
chainZ2J0->Add("/data/1a/madfish/JuneData/Z2Jets_ptZ_0to100.root");
chainZ2J0->Project("ZJ5",V->Data(),cuts->Data());

TChain *chainZ2J100 = new TChain("Events");
chainZ2J100->SetWeight(lumi/lumiZ2jets100_300 ,"global");
chainZ2J100->Add("/data/1a/madfish/JuneData/Z2Jets_ptZ_100to300.root");
chainZ2J100->Project("ZJ6",V->Data(),cuts->Data());

TChain *chainZ2J300 = new TChain("Events");
chainZ2J300->SetWeight(lumi/lumiZ2jets300_800 ,"global");
chainZ2J300->Add("/data/1a/madfish/JuneData/Z2Jets_ptZ_300to800.root");
chainZ2J300->Project("ZJ7",V->Data(),cuts->Data());

TChain *chainZ2J800 = new TChain("Events");
chainZ2J800->SetWeight(lumi/lumiZ2jets800_1600 ,"global");
chainZ2J800->Add("/data/1a/madfish/JuneData/Z2Jets_ptZ_800to1600.root");
chainZ2J800->Project("ZJ8",V->Data(),cuts->Data());
*/

TChain *CMuv1 = new TChain("Events");
CMuv1->Add("/data/1a/madfish/JuneData/Muv1.root");
CMuv1->Project("Muv1",V->Data(),cuts->Data());


TChain *CMuv2 = new TChain("Events");
CMuv2->Add("/data/1a/degrutto/data/data11/Muv2/Muv2.root");
CMuv2->Project("Muv2",V->Data(),cuts->Data());

TChain *CEv1 = new TChain("Events");
CEv1->Add("/data/1a/madfish/JuneData/Ev1.root");
CEv1->Project("Ev1",V->Data(),cuts->Data());


TChain *CEv2 = new TChain("Events");
CEv2->Add("/data/1a/degrutto/data/data11/Elev2/Elev2.root");
CEv2->Project("Ev2",V->Data(),cuts->Data());

Ev1->Add(Ev2);
Muv1->Add(Muv2);

Ev1->SetMarkerStyle(kFullCircle);
Ev1->SetMarkerSize(1.2);
Ev1->SetMarkerColor(kBlack);
Ev1->SetLineWidth(2);
Ev1->SetLineColor(kBlack);

Muv1->SetMarkerStyle(kFullCircle);
Muv1->SetMarkerSize(1.2);
Muv1->SetMarkerColor(kBlack);
Muv1->SetLineWidth(2);
Muv1->SetLineColor(kBlack);
//gStyle->SetErrorX(.5);
gStyle->SetEndErrorSize(2);

/*ZJets->Add(ZJ2);
ZJets->Add(ZJ3);
ZJets->Add(ZJ4);
ZJets->Add(ZJ5);
ZJets->Add(ZJ6);
ZJets->Add(ZJ7);
ZJets->Add(ZJ8);*/

THStack * hs = new THStack("hs","");

ZJets->SetFillColor(kAzure+2);
ZZ->SetFillColor(kViolet);
TTbar->SetFillColor(kRed+1);
ZH->SetLineColor(kRed);

float maximum=0;

hs->Add(TTbar);
hs->Add(ZJets);
hs->Add(ZZ);
//hs->Add(ZH);

std::cout << "TTbar: "<< TTbar->Integral() <<std::endl;
std::cout << "ZJets: "<< ZJets->Integral() <<std::endl;
std::cout << "ZZ: "<< ZZ->Integral() <<std::endl;
std::cout << "ZH: "<< ZH->Integral() <<std::endl;
std::cout << "Mu Data: "<< Muv1->Integral() <<std::endl;
std::cout << "Ele Data: "<< Ev1->Integral() <<std::endl;

if(chan == "mm" || chan == "MM")
{
if(hs->GetMaximum() > Muv1->GetMaximum()) maximum = hs->GetMaximum();
  else maximum = Muv1->GetMaximum();
}
else if(chan == "ee" || chan == "EE")
{
if(hs->GetMaximum() > Ev1->GetMaximum()) maximum = hs->GetMaximum();
  else maximum = Ev1->GetMaximum();
}

hs->SetMaximum(maximum);
Muv1->SetMaximum(maximum);
hs->Draw("HIST");
/*hs->GetXaxis()->SetNdivisions(506);
hs->GetYaxis()->SetNdivisions(506);*/
hs->GetXaxis()->SetTitle(A->Data());
//hs->SetTitle(axis);

if(chan == "mm" || chan == "MM")
{
Muv1->Draw("PE1SAME");
}
else
if(chan == "ee" || chan == "ee")
{
Ev1->Draw("PE1SAME");
}


TLegend* legend=new TLegend(0.65,0.5,0.9,0.70);
legend->SetLineColor(0);
legend->SetFillColor(0);



//leg = new TLegend(0.20,0.7,0.35,0.85);
legend->AddEntry(Muv1,"Data","f");
//legend->AddEntry(ZH,"ZH","f");
legend->AddEntry(ZJets,"Z+Jets","f");
legend->AddEntry(ZZ,"ZZ","f");
legend->AddEntry(TTbar,"tt#bar","f");
legend->SetShadowColor(kWhite);
legend->Draw("SAME");

TLatex latex;
latex.SetNDC();
latex.SetTextSize(0.04);

latex.SetTextAlign(31); // align right
latex.DrawLatex(0.70,0.93,"#sqrt{s} = 7 TeV");
latex.SetTextAlign(31); // align right
latex.DrawLatex(0.8,0.8,Form("#int #font[12]{L} dt = %.0f pb^{-1}",lumi*1000));

latex.SetTextAlign(11); // align left
latex.DrawLatex(0.12,0.93,"CMS simulation 2011");





//chainZ2J100->Draw(V->Data(),"hjjMass>102 && hjjMass< 126 && hjjMass@.size() < 2 && hjjJet1Pt>20 && hjjJet2Pt>20 && zmmMass>75  && zmmMass<105 && ((hjjJet1CSV>0.9 && hjjJet2CSV>0.44) || (hjjJet2CSV>0.9 && hjjJet1CSV>0.44)) && zmmPt > 150 &&  hjjPt > 150 ");
//





/*

TFile *fzj = new TFile("ZllHZbb0Jets.root","read");
TFile *fzh = new TFile("ZllHZllHbb115.root");

TTree *treezh= (TTree*) fzh->Get("tree");
TTree *treezj= (TTree*) fzj->Get("tree");
std::string var,name,axis;
float max;

std::cout << "Which Variable Would You Like To Plot?" << std::endl;
std::cin >> var;
std::cout << "Name of Plot?" << std::endl;
std::cin.ignore();
std::getline(std::cin, name);
std::cout << "X Axis Title?" << std::endl;
std::getline(std::cin, axis);


TString *V = new TString(var);
TString *N = new TString(name);
TString *A = new TString(axis);

if(max != 0)
{
if (treezh->GetMaximum(V->Data()) > treezj->GetMaximum(V->Data())) 
  {max = treezh->GetMaximum(V->Data());}
else max = treezj->GetMaximum(V->Data());
}
else {std::cout << "Maximum X Axis? "<< std::endl; std::cin >>max;}

std::cout << "max "<< max << std::endl;
TH1F *thistzh = new TH1F("thistzh",N->Data(),40,0,max+fabs(max/5.0));
TH1F *thistzj = new TH1F("thistzj",N->Data(),40,0,max+fabs(max/5.0));

TString *cuts = new TString("");

//TString *cuts = new TString("hjjJet1CSV > 0.6 && hjjJet2CSV > 0.6 && hjjJet1Pt > 20 && hjjJet2Pt > 20 && deltaPullAngle > -10 && H.pt > 150 && Z.pt > 150 && hjjzdPhi > 2.7");

//TString *cuts = new TString("hjjJet1CSV > 0.44 && hjjJet2CSV > 0.44 && Z.mass > 60 && Z.mass < 120 && hjjMass > 80 && hjjMass < 140 && hjjJet1Pt > 20 && hjjJet2Pt > 20 && deltaPullAngle > -10 && H.pt > 150");


treezh->Project("thistzh",V->Data(),cuts->Data());
treezj->Project("thistzj",V->Data(),cuts->Data());

thistzh->GetXaxis()->SetTitle(A->Data());
thistzh->SetLineColor(kBlue);
thistzh->SetLineWidth(3);
thistzj->SetLineColor(kRed);
thistzj->SetLineWidth(3);

thistzh->SetStats(0);
thistzj->SetStats(0);

float zhmax = thistzh->GetMaximum()/thistzh->Integral();
float zjmax = thistzh->GetMaximum()/thistzh->Integral();

if(zhmax > zjmax)
{
thistzh->DrawNormalized();
thistzj->DrawNormalized("SAME");
}
else
{
thistzj->DrawNormalized();
thistzh->DrawNormalized("SAME");
}



TLegend *tl = new TLegend(0.53,0.75,0.63,0.85);
tl->SetTextSize(0.04);
tl->SetTextFont(42);
tl->SetBorderSize(0);
tl->SetFillColor(0);
tl->AddEntry(thistzh,"ZH-115");
tl->AddEntry(thistzj,"Z+Jets");
tl->SetFillColor(0);
tl->Draw("SAME");

*/

}
