#ifndef STACK_COMMON_ZLL_H
#define STACK_COMMON_ZLL_H


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

bool CorrectIt  = 0;

const int canvasSizeX = 500;
const int canvasSizeY = 500;

const Color_t hLineColor = kBlack;
const Color_t hFillColor = 0;


const Color_t vjLineColor = kAzure+2;
const Color_t vjFillColor = kAzure+6;

const Color_t vvLineColor = kViolet+3;
const Color_t vvFillColor = kViolet-5;


const Color_t ttLineColor =  kRed+3;
const Color_t ttFillColor = kRed+1;
double ratioToDataIntMG;


double lumi;

bool option2 = 1; 
const double lumiZH = 7652.4;
 
const double lumiWH = 1321.; 

const double lumiZJ = 11.902; 
const double lumiWJ = 2.6; 


const double lumiWZ = 233.073;
const double lumiWW = 98.506;
const double lumiZZ = 659.0; 
//const double lumiZZ = 357; 

const double lumiTT = 15.84; 

const double lumiTs = 499.; 
const double lumiTt = 23.; 
const double lumiTtW = 46.; 

TCanvas *c1 = new TCanvas("c1","Stack plot", 300,300,479,510);




double Mass(double pt, double eta, double phi, double mass, double pt2, double eta2, double phi2, double mass2)
{
  TLorentzVector m1, m2, msum;
  m1.SetPtEtaPhiM(pt, eta, phi, mass);
  m2.SetPtEtaPhiM(pt2, eta2, phi2, mass2);
  msum = m1 + m2;
  return(msum.M());
}

double PT(double pt, double eta, double phi, double mass, double pt2, double eta2, double phi2, double mass2)
{
  TLorentzVector m1, m2, msum;
  m1.SetPtEtaPhiM(pt, eta, phi, mass);
  m2.SetPtEtaPhiM(pt2, eta2, phi2, mass2);
  msum = m1 + m2;
  return(msum.Pt());
}



double deltaPhi(double phi1, double phi2) 
{ 
    double PI = 3.14159265;
    double result = phi1 - phi2;
    while (result > PI) result -= 2*PI;
    while (result <= -PI) result += 2*PI;
    return result;
}

double deltaR(double eta1, double phi1, double eta2, double phi2) 
{
    double deta = eta1 - eta2;
    double dphi = deltaPhi(phi1, phi2);
    return sqrt(deta*deta + dphi*dphi);
}



void setHisto(TH1 * h, Color_t fill, Color_t line, double scale) {
  h->SetFillColor(fill);
  h->SetLineColor(line);
  h->Scale(scale);
}



void stat(TH1 * hZH, TH1 * hWH, TH1 * hZJl, TH1 * hZJb,  TH1 * hWJ,  TH1 * hZZ, TH1 * hWZ, TH1 * hWW,   TH1 * hTT, TH1 *hTs , TH1 *hTt , TH1 *hTtW ,  TH1* hdata, bool scaleToDataInt,int option, float rebin, float mMin, float mMax) 
{


   double a = mMin/rebin +1, b = mMax/rebin;

   if(option2){   std::cout << "Integrating From " << mMin << "-" <<mMax <<std::endl;
   std::cout << "or Bin " << a << "-" << b <<std::endl;}
   
   double iZH = hZH->Integral(a,b);
   double errZH = sqrt(lumi/lumiZH *  iZH);

   double iWH = hWH->Integral(a,b);
   double errWH = sqrt(lumi/lumiWH *  iWH);

   double iZJl = hZJl->Integral(a,b);
   double iZJb = hZJb->Integral(a,b);

   double iWJ = hWJ->Integral(a,b);
   double errWJ = sqrt(lumi/lumiWJ *  iWJ);

   double iZZ = hZZ->Integral(a,b);
   double errZZ = sqrt(lumi/lumiZZ *  iZZ);

   double iWW = hWW->Integral(a,b);
   double errWW = sqrt(lumi/lumiWW *  iWW);

   double iWZ = hWZ->Integral(a,b);
   double errWZ = sqrt(lumi/lumiWZ *  iWZ);

   double iTT = hTT->Integral(a,b);
   double errTT = sqrt(lumi/lumiTT *  iTT);

   double iTs = hTs->Integral(a,b);
   double errTs = sqrt(lumi/lumiTs *  iTs);

   double iTt = hTt->Integral(a,b);
   double errTt = sqrt(lumi/lumiTt *  iTt);

   double iTtW = hTtW->Integral(a,b);
   double errTtW = sqrt(lumi/lumiTtW *  iTtW);

   double idata = hdata->Integral(a,b); 
   //  double errData = sqrt(idata);
   if(option2)
   {
   std::cout.setf(ios::fixed,ios::floatfield);
   std::cout.precision(2);
   std::cout <<"ZH  = ";
   std::cout.precision(2);
   std::cout << iZH << " +/- " << errZH <<std::endl;
   std::cout.precision(2);
   std::cout <<"WH  = ";
   std::cout.precision(2);
   std::cout << iWH << " +/- " << errWH <<std::endl;
   std::cout.precision(2);
   std::cout <<"ZJl  = ";
   std::cout.precision(2);
   std::cout << iZJl << std::endl;
   std::cout.precision(2);

   std::cout <<"ZJb  = ";
   std::cout.precision(2);
   std::cout << iZJb << std::endl;


   std::cout.precision(2);
   std::cout <<"WJ = "; 
   std::cout.precision(2);
   std::cout << iWJ << " +/- " << errWJ <<std::endl; 
   std::cout.precision(2);

   std::cout <<"WW = "; 
   std::cout.precision(2);
   std::cout << iWW << " +/- " << errWW <<std::endl; 
   std::cout.precision(2);

   std::cout <<"WZ = "; 
   std::cout.precision(2);
   std::cout << iWZ << " +/- " << errWZ <<std::endl; 
   std::cout.precision(2);

   std::cout <<"ZZ = "; 
   std::cout.precision(2);
   std::cout << iZZ << " +/- " << errZZ <<std::endl; 
   std::cout.precision(2);

   std::cout <<"TT  = "; 
   std::cout.precision(2);
   std::cout << iTT << " +/- " << errTT <<std::endl; 
   std::cout.precision(2);
  std::cout <<"Ts = "; 
   std::cout.precision(2);
   std::cout << iTs << " +/- " << errTs <<std::endl; 
   std::cout.precision(2);
  std::cout <<"Tt = "; 
   std::cout.precision(2);
   std::cout << iTt << " +/- " << errTt <<std::endl; 
   std::cout.precision(2);
  std::cout <<"TtW = "; 
   std::cout.precision(2);
   std::cout << iTtW << " +/- " << errTtW<<std::endl; 
   std::cout.precision(2);
   float B =  iZJl+iZJb+iWJ+iZZ+iWW+iWZ+iTT+iTs+iTt+iTtW;
  std::cout << "Total B = " << B  << std::endl;
if(option < 2) {   std::cout <<"Data  = ";
   std::cout <<idata  << std::endl; }
  std::cout << "Data/MC = " << idata/B  << std::endl;

   float S = iZH;
   float siggie = S/(1.5+sqrt(B)+0.1*B);


   std::cout.precision(3);
  std::cout << "Significance: " << siggie << std::endl;
}
//else std::cout <<siggie;


 //if one wnats one can scale any MC to the data events 
   if (scaleToDataInt && (option < 2)) {
     double ratioToDataInt = idata / ( iZZ + iWZ + iWW + iZJl +  iZJb + iWJ + iTT + iTs + iTt + iTtW ); 
     ratioToDataIntMG = idata / ( iZZ + iWZ + iWW +  iWJ + iTT + iTs + iTt + iTtW +iZJl +  iZJb); 
     hZH->Scale(ratioToDataInt);
     hWH->Scale(ratioToDataInt);
     hZJl->Scale(ratioToDataInt);
     hZJb->Scale(ratioToDataInt);
     hWJ->Scale(ratioToDataInt);
     hZZ->Scale(ratioToDataInt);
     hWW->Scale(ratioToDataInt);
     hWZ->Scale(ratioToDataInt);
     hTT->Scale(ratioToDataInt);
     hTs->Scale(ratioToDataInt);
     hTt->Scale(ratioToDataInt);
     hTtW->Scale(ratioToDataInt);
   }


 }

 void makeStack(TH1 * hZH, TH1 * hWH, TH1 * hZJl, TH1 * hZJb, TH1 * hWJ, TH1 * hZZ, TH1 * hWZ, TH1 * hWW, TH1 * hTT, TH1 * hTs, TH1* hTt, TH1 * hTtW,  TH1 * hdata, double min, int option, bool logScale, bool scaleToDataInt, TString *title, float epsilonBin, float a, float b) {

  
   setHisto(hZH, 0, kBlack, lumi/lumiZH);
   hZH->SetLineWidth(3);
   
   setHisto(hZJl, kGreen + 3, kGreen +4,lumi/lumiZJ);
   setHisto(hZJb, kSpring, kSpring, lumi/lumiZJ);
   setHisto(hWH, hFillColor, hLineColor, lumi/lumiWH);

   setHisto(hWJ, vjFillColor, vjLineColor, lumi/lumiWJ);

   setHisto(hZZ, kRed -2, kRed -2, lumi/lumiZZ);
    setHisto(hWZ, vvFillColor, vvLineColor, lumi/lumiWZ);
   setHisto(hWW, vvFillColor, vvLineColor, lumi/lumiWW);

   setHisto(hTT, kBlue, kBlue, lumi/lumiTT);
   setHisto(hTs, kBlue, kBlue, lumi/lumiTs);
   setHisto(hTt, kBlue, kBlue, lumi/lumiTt);
   setHisto(hTtW, kBlue, kBlue, lumi/lumiTtW);






   if(CorrectIt) {hTT->Scale(1.01); hZJl->Scale(0.95); hZJb->Scale(1.06); }

   stat(hZH, hWH, hZJl,hZJb,hWJ, hZZ,hWZ, hWW, hTT, hTs, hTt, hTtW, hdata, scaleToDataInt,option,epsilonBin,a,b);

//   hZH->Add(hWH);

   
   hZJl->Add(hWJ);
   hZZ->Add(hWZ);
   hZZ->Add(hWW);

   hTT->Add(hTs);
   hTT->Add(hTt);
   hTT->Add(hTtW);

   THStack * hs = new THStack("hs","");

  
   hs->Add(hZZ);
   hs->Add(hTT);
   hs->Add(hZJl);
   hs->Add(hZJb);
   hs->Draw("HIST");






 
   TString yTag = TString::Format("Number of Events / %.2f", epsilonBin);

   hs->SetMinimum(min);
   hs->GetXaxis()->SetTitle(title->Data());
   hs->GetYaxis()->SetTitle(yTag.Data());
   hs->GetXaxis()->SetTitleOffset(1.1);
   hs->GetYaxis()->SetTitleOffset(1.2);

   hs->GetYaxis()->SetLabelSize(.04);
   hs->GetYaxis()->SetTitleSize(.05);
   hs->GetXaxis()->SetLabelSize(.04);
   hs->GetXaxis()->SetTitleSize(.06);

   hdata->SetMinimum(min);
   hdata->GetXaxis()->SetTitle(title->Data());
   hdata->GetYaxis()->SetTitle(yTag.Data());
   hdata->GetXaxis()->SetTitleOffset(1.1);
   hdata->GetYaxis()->SetTitleOffset(1.2);

   hdata->GetYaxis()->SetLabelOffset(0.0);
   hdata->GetXaxis()->SetLabelSize(.04);
   hdata->GetYaxis()->SetLabelSize(.04);
   hdata->GetYaxis()->SetTitleSize(.05);
   hdata->GetXaxis()->SetTitleSize(.06);
   hdata->SetTitle("");
   hdata->SetStats(0);

if(option < 2) 
{

  c1->SetLeftMargin(  87./479 );
  c1->SetRightMargin( 42./479 );
  c1->SetTopMargin(  30./510 );
  c1->SetBottomMargin( 80./510 ); 
  c1->SetFillColor(0);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetFrameFillStyle(0);
  c1->SetFrameLineWidth(2);
  c1->SetFrameBorderMode(0);
  int lineWidth(2);

if( logScale )
  {
    lineWidth = 1;
  }

 // titles and axis, marker size
  TString ytitle;
int ndivx(506);
  int ndivy(506);
  float markerSize(2.);

  // canvas name
  TString cname("c");
  TString ctitle;

  // legend position and scale;
  float xl_  = 0.;
  float yl_  = 0.6;
  float scalel_ = 0.05;
      ndivx = 120;

if( logScale )
	{
	  ndivx=506; 
	  ndivy = 506;
	}
      else
	{
	  ndivy = 506;
          ndivx=506;
	}

      if( logScale )
	{
	  markerSize = 1.1;
	}
      else
	{	
	  markerSize = 1.2;
	}


      if( logScale )
	{
	  xl_ = 0.60;
	  yl_ = 0.65;
	}
      else
	{
	  xl_ = 0.60;
	  yl_ = 0.60;
	  scalel_ = 0.06;
	}
      hs->GetXaxis()->SetNdivisions(ndivx);
      hs->GetYaxis()->SetNdivisions(ndivy);
      // log plots, so the maximum should be one order of magnitude more...
      hs->SetMaximum( TMath::Max( hdata->GetMaximum(), hs->GetMaximum()) + 0.5) ;
      //hs->SetMaximum( 29.8 )  ;
      if (logScale) 
      hs->SetMaximum( TMath::Max( hdata->GetMaximum(), hs->GetMaximum()) +250) ;


     hdata->SetMarkerStyle(kFullCircle);
     hdata->SetMarkerSize(markerSize);
     hdata->SetMarkerColor(kBlack);
     hdata->SetLineWidth(markerSize);
     hdata->SetLineColor(kBlack);
     gStyle->SetEndErrorSize(1);
     if(option == 0)   hdata->Draw("PE1SAME");
     if(option == 1){hdata->SetLineWidth(3);   hdata->Draw();}
     else
     {
     hdata->GetXaxis()->SetLabelSize(0);
     hdata->GetYaxis()->SetLabelSize(0);
     hdata->GetXaxis()->SetNdivisions(ndivx);
     hdata->GetYaxis()->SetNdivisions(ndivy);
     }
    
}




hZH->SetMinimum(min);
hZH->SetLineWidth(3);
hZH->GetXaxis()->SetTitle(title->Data());
hZH->GetYaxis()->SetTitle(yTag.Data());
hZH->GetXaxis()->SetTitleOffset(1.1);
hZH->GetYaxis()->SetTitleOffset(1.2);
hZH->GetYaxis()->SetLabelOffset(0.0);
hZH->GetXaxis()->SetLabelSize(.04);
hZH->GetYaxis()->SetLabelSize(.04);
hZH->GetYaxis()->SetTitleSize(.05);
hZH->GetXaxis()->SetTitleSize(.06);
hZH->SetTitle("");
hZH->SetStats(0);

if(option==0) {hZH->Draw("SAME");}
if(option==2) {hZH->Draw("SAME");}
if(option==3) {hZH->Draw();}


 // middle right
 TLegend* legend;   
 legend=new TLegend(0.65,0.70,0.9,0.85);
 //else legend = new TLegend(0.6,0.72,0.8,0.80);
   legend->SetLineColor(0);  legend->SetFillColor(0);





   //leg = new TLegend(0.20,0.7,0.35,0.85);
   if(option < 2)
     if (option <2) legend->AddEntry(hdata,"data", "pl");
   legend->AddEntry("hZH","ZH","pl");
   legend->AddEntry(hZZ,"ZZ+WW+WZ","f");
   legend->AddEntry(hZJb,"Z + b Jets","f");
   legend->AddEntry(hZJl,"W/Z + udcs Jets","f");
   legend->AddEntry(hTT,"t#bar{t} + st","f"); 




  legend->SetShadowColor(kWhite);
//   legend->Draw();

//*/

 TLatex latex;
   latex.SetNDC();
   latex.SetTextSize(0.03);

   latex.SetTextAlign(31); // align right
   latex.DrawLatex(0.90,0.96,"#sqrt{s} = 7 TeV");
   if (lumi > 0.) {
     latex.SetTextAlign(31); // align right
     latex.DrawLatex(0.89,0.89,Form("#int #font[12]{L} dt = %.2f fb^{-1}",lumi+0.005));
   }
   latex.SetTextAlign(11); // align left
   if(scaleToDataInt)   latex.DrawLatex(0.12,0.96,"CMS preliminary 2011: Shape");
   else   latex.DrawLatex(0.12,0.96,"CMS preliminary 2011");



 
 }






 void makePlots(const char * var, string xtitle,   TCut cut,  const char * plot, double min , float epsilonBin , double xMin, double xMax,  int option , bool logScale,bool scaleToDataInt, float a, float b) {

  std::cout << "Cut: " << cut.GetName() << std::endl; 
  c1->SetLeftMargin(  87./479 );
  c1->SetRightMargin( 42./479 );
  c1->SetTopMargin(  30./510 );
  c1->SetBottomMargin( 80./510 ); 
  c1->SetFillColor(0);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetFrameFillStyle(0);
  c1->SetFrameLineWidth(2);
  c1->SetFrameBorderMode(0);



 
 int lineWidth(2);

if( logScale )
  {
    lineWidth = 1;
  }

  // titles and axis, marker size
//  TString xtitle;i
  TString ytitle;
int ndivx(506);
  int ndivy(506);
  float markerSize(2.);

  // canvas name
  TString cname("c");
  TString ctitle;

  // legend position and scale;
  float xl_  = 0.;
  float yl_  = 0.6;
  float scalel_ = 0.05;
      ndivx = 120;

if( logScale )
	{
	  ndivx=506; 
	  ndivy = 506;
	}
      else
	{
	  ndivy = 506;
          ndivx=506;
	}

      if( logScale )
	{
	  markerSize = 1.1;
	}
      else
	{	
	  markerSize = 1.2;
	}


      if( logScale )
	{
	  xl_ = 0.60;
	  yl_ = 0.65;
	}
      else
	{
	  xl_ = 0.60;
	  yl_ = 0.60;
	  scalel_ = 0.06;
	}
 
   bool eeChan = false;
   TString *V = new TString(var);
   TString *C = new TString(cut.GetTitle());
   if(V->Contains("zee") || C->Contains("zee")) eeChan = true;
   else eeChan = false;
   
   if(eeChan) lumi = 1.085;
   else lumi = 1.085;

   if(option > 1) lumi = 1;
  
   TString *xTitle = new TString(xtitle);

   TChain * zhtree = new TChain("tree"); 
//   zhtree->Add("ZllHZllHbb/EdmNtuples_Zll_2_1_8l0.root");
   zhtree->Add("42-ZllHbb115.root");


   TChain * whtree = new TChain("tree"); 
   whtree->Add("ZllHWlnHbb115.root");

   TChain * zjtree = new TChain("tree"); 
   zjtree->Add("42-ZJetsMad.root");

   TChain * wjtree = new TChain("tree"); 
   wjtree->Add("42-WJetsMad.root");



   TChain * zztree = new TChain("tree"); 
   zztree->Add("42-ZZ.root");

    TChain * tttree = new TChain("tree"); 
   tttree->Add("42-TTbar.root");


   TChain * ststree = new TChain("tree"); 
   ststree->Add("ZllHT-s.root");

   TChain * stttree = new TChain("tree"); 
   stttree->Add("ZllHT-t.root");

   TChain * sttwtree = new TChain("tree"); 
   sttwtree->Add("ZllHT-tW.root");

  TChain * wwtree = new TChain("tree"); 
   wwtree->Add("42-WW.root");

   TChain * wztree = new TChain("tree"); 
   wztree->Add("42-WZ.root");







   TChain * datatree2011= new TChain("tree");

if(eeChan) datatree2011->Add("Elev11v2/Elev2.root"); 
else datatree2011->Add("42-SingleMu1fb.root");

   //##/data/1a/degrutto/data/data11/METv2/METv2.root");
  //  ###datatree2011->Add("");
  




  TH1F *hdata = new TH1F ("hdata", "hdata", (xMax-xMin)/epsilonBin, xMin, xMax);
  TH1F *hZH = new TH1F ("hZH", "hZH", (xMax-xMin)/epsilonBin, xMin, xMax);
  TH1F *hWH = new TH1F ("hWH", "hWH", (xMax-xMin)/epsilonBin, xMin, xMax);
  
  TH1F *hZJl = new TH1F ("hZJl", "hZJl", (xMax-xMin)/epsilonBin, xMin, xMax);
  TH1F *hZJb = new TH1F ("hZJb", "hZJb", (xMax-xMin)/epsilonBin, xMin, xMax);
  TH1F *hWJ = new TH1F ("hWJ", "hWJ", (xMax-xMin)/epsilonBin, xMin, xMax);


  TH1F *hZZ = new TH1F ("hZZ", "hZZ", (xMax-xMin)/epsilonBin, xMin, xMax);
  TH1F *hWZ = new TH1F ("hWZ", "hWZ", (xMax-xMin)/epsilonBin, xMin, xMax);
  TH1F *hWW = new TH1F ("hWW", "hWW", (xMax-xMin)/epsilonBin, xMin, xMax);

  TH1F *hTT = new TH1F ("hTT", "hTT", (xMax-xMin)/epsilonBin, xMin, xMax);
  TH1F *hTs = new TH1F ("hTs", "hTs", (xMax-xMin)/epsilonBin, xMin, xMax);
  TH1F *hTt = new TH1F ("hTt", "hTt", (xMax-xMin)/epsilonBin, xMin, xMax);
  TH1F *hTtW = new TH1F ("hTtW", "hTtW", (xMax-xMin)/epsilonBin, xMin, xMax);

if((option == 0) || (option == 2))
{



zjtree->Project("hZJl", var, cut+"!((abs(Jet1.flavour) == 5) && (abs(Jet2.flavour) == 5))");
zjtree->Project("hZJb", var, cut+"((abs(Jet1.flavour) == 5) && (abs(Jet2.flavour) == 5))");

TString *mccut = new TString("(");
mccut->Append(cut.GetTitle()); 
mccut->Append(")*PUweight"); 
std::cout << "mcc cut " << mccut ->Data() << std::endl;

zhtree->Project("hZH", var, cut);
  
zztree->Project("hZZ", var, mccut->Data());

tttree->Project("hTT", var, mccut->Data());

ststree->Project("hTs", var, cut);
stttree->Project("hTt", var, cut);
wztree->Project("hWZ", var, mccut->Data());
////whtree->Project("hWH", var, cut);
wwtree->Project("hWW", var, mccut->Data());
wjtree->Project("hWJ", var, mccut->Data());
sttwtree->Project("hTtW", var, cut);
}
  
if(option == 3) zhtree->Project("hZH", var, cut);
 

if(option < 2)  datatree2011->Project("hdata", var, cut) ;
    
  
  makeStack(hZH, hWH, hZJl, hZJb, hWJ, hZZ, hWZ, hWW, hTT, hTs, hTt, hTtW, hdata, min, option, logScale, scaleToDataInt, xTitle, epsilonBin,a,b);

float maxSig1=0; 
float maxSig2=0; 
float maxSig3=0; 
float maxSig4=0; 
float maxSig5=0;
 
float tempSigU1=0;float tempSigD1=0;float tempSigU2=0;float tempSigD2=0;float tempSigU3=0;float tempSigD3=0;float tempSigU4=0;float tempSigD4;float x = 0.1; 
int goodbin1=0;int goodbin2=0;int goodbin3=0;int goodbin4=0;int goodbin5; 
bool GT1=0;bool GT2=0;bool GT3=0;bool GT4=0;bool GT5; 
float mS1=0;float mS2=0;float mS3=0;float mB1=0;float mB2=0;float mB3=0;float mB4=0;float mS4=0;float mS5=0;float mB5; 

for(int i=1; i < hZH->GetNbinsX();i++)
{

   float Bu = hZZ->Integral(i,9999) + hTT->Integral(i,9999) + hZJl->Integral(i,9999) +  hZJb->Integral(i,9999);
   float Bd = hZZ->Integral(0,i) + hTT->Integral(0,i) + hZJl->Integral(0,i)+ hZJb->Integral(0,i);
   float Su = hZH ->Integral(i,9999);
   float Sd = hZH ->Integral(0,i);

   tempSigU1 =  Su/(1.5+sqrt(Bu + x*x*Bu*Bu));
   tempSigU1 =  Sd/(1.5+sqrt(Bd + x*x*Bd*Bd));
   
   tempSigU2 = Su / (1.5 + sqrt(Bu));
   tempSigD2 = Sd / (1.5 + sqrt(Bd));
   
   tempSigU3 = Su / (1.5 + sqrt(Bu) + x*Bu);
   tempSigD3 = Sd / (1.5 + sqrt(Bd) + x*Bd);

   tempSigU4 = Su / (sqrt(Bu + x*x*Bu*Bu));
   tempSigD4 = Sd / (sqrt(Bd + x*x*Bd*Bd));

float   tempSigU5 = Su / (sqrt(Su + Bu));
float   tempSigD5 = Sd / (sqrt(Bd + Bd));



   std::cout.precision(8);
/*   std::cout << "Sig5->SetBinContent(" << i << "," << tempSigU5 <<");" << std::endl;
   std::cout << "Sig4->SetBinContent(" << i << "," << tempSigU4 <<");" << std::endl;
   std::cout << "Sig3->SetBinContent(" << i << "," << tempSigU2 <<");" << std::endl;
   std::cout << "Sig2->SetBinContent(" << i << "," << tempSigU2 <<");" << std::endl;
   std::cout << "Sig1->SetBinContent(" << i << "," << tempSigU1 <<");" << std::endl;
*/



   if(tempSigU1 > maxSig1) {maxSig1 = tempSigU1; GT1 = true;  goodbin1 = i; mS1 = Su; mB1 = Bu;}
   if(tempSigD1 > maxSig1) {maxSig1 = tempSigD1; GT1 = false; goodbin1 = i; mS1 = Sd; mB1 = Bd;}

   if(tempSigU2 > maxSig2) {maxSig2 = tempSigU2; GT2 = true;  goodbin2 = i; mS2 = Su; mB2 = Bu;}
   if(tempSigD2 > maxSig2) {maxSig2 = tempSigD2; GT2 = false; goodbin2 = i; mS2 = Sd; mB2 = Bd;}

   if(tempSigU3 > maxSig3) {maxSig3 = tempSigU3; GT3 = true;  goodbin3 = i; mS3 = Su; mB3 = Bu;}
   if(tempSigD3 > maxSig3) {maxSig3 = tempSigD3; GT3 = false; goodbin3 = i; mS3 = Sd; mB3 = Bd;}

   if(tempSigU4 > maxSig4) {maxSig4 = tempSigU4; GT4 = true;  goodbin4 = i; mS4 = Su; mB4 = Bu;}
   if(tempSigD4 > maxSig4) {maxSig4 = tempSigD4; GT4 = false; goodbin4 = i; mS4 = Sd; mB4 = Bd;}

   if(tempSigU5 > maxSig5) {maxSig5 = tempSigU4; GT5 = true;  goodbin5 = i; mS5 = Su; mB5 = Bu;}
   if(tempSigD5 > maxSig5) {maxSig5 = tempSigD4; GT5 = false; goodbin5 = i; mS5 = Sd; mB5 = Bd;}

}


  std::cout << "Cut1 S/(1.5+sqrt(B + x*x*B*B)) At: " << hZH ->GetBinCenter(goodbin1) << " GT: " << GT1 << " S: " << mS1 << " B: " << mB1 <<  " " << maxSig1  <<std::endl;
  std::cout << "Cut2 S/(1.5 + sqrt(B))         At: " << hZH ->GetBinCenter(goodbin2) << " GT: " << GT2 << " S: " << mS2 << " B: " << mB2 <<  " " << maxSig2  <<std::endl;
  std::cout << "Cut3 S/(1.5 + sqrt(B) + x*B)   At: " << hZH ->GetBinCenter(goodbin3) << " GT: " << GT3 << " S: " << mS3 << " B: " << mB3 <<  " " << maxSig3  <<std::endl;
  std::cout << "Cut4 S/(sqrt(B+x*x*Bd*Bd))     At: " << hZH ->GetBinCenter(goodbin4) << " GT: " << GT4 << " S: " << mS4 << " B: " << mB4 <<  " " << maxSig4  <<std::endl;
  std::cout << std::endl << std::endl;


  if (logScale) c1->SetLogy();
if(option2)
{
  c1->SaveAs((std::string(plot)+".eps").c_str());
  c1->SaveAs((std::string(plot)+".png").c_str());
  c1->SaveAs((std::string(plot)+".pdf").c_str());
  
  TFile * out = new TFile("plot.root", "RECREATE");

  c1->Write();
  c1->SaveAs("hPlot.C");


  out->Close(); 
}

/*
  hZH->Delete();
  hWH->Delete();
  
  hWJ->Delete();
  hZZ ->Delete();
  hWZ ->Delete();
  hWW ->Delete();

  hTT ->Delete();
  hTs ->Delete();
  hTt ->Delete();
  hTtW ->Delete();

  hZJl->Delete();
  hZJb->Delete();
zhtree->Delete();
whtree ->Delete(); 
zjtree  ->Delete();
zztree->Delete();
tttree->Delete();
datatree2011->Delete();
hdata->Delete();
wztree->Delete();
wwtree->Delete();
ststree->Delete();
stttree->Delete();
sttwtree ->Delete();
*/

}


#endif
