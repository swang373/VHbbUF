#include <iostream>

void plotVar()
{
gROOT->SetStyle("Plain");
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

//TString *cuts = new TString("Jet1.csv > 0.6 && Jet2.csv > 0.6 && Jet1.pt > 20 && Jet2.pt > 20 && deltaPullAngle > -10 && H.pt > 150 && Z.pt > 150 && hjjzdPhi > 2.7");

//TString *cuts = new TString("Jet1.csv > 0.44 && Jet2.csv > 0.44 && Z.mass > 60 && Z.mass < 120 && H.mass > 80 && H.mass < 140 && Jet1.pt > 20 && Jet2.pt > 20 && deltaPullAngle > -10 && H.pt > 150");


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



}
