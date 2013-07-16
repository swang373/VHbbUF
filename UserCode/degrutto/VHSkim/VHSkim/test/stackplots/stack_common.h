#ifndef STACK_COMMON_H
#define STACK_COMMON_H


#include <iostream>
using namespace std;
#include "TChain.h"

#include "TGraphAsymmErrors.h"

const int canvasSizeX = 500;
const int canvasSizeY = 500;

const Color_t hLineColor = kBlack;
const Color_t hFillColor = 10;


const Color_t vjLineColor = kAzure+2;
const Color_t vjFillColor = kAzure+6;

const Color_t vvLineColor = kViolet+3;
const Color_t vvFillColor = kViolet-5;


const Color_t ttLineColor =  kRed+3;
const Color_t ttFillColor = kRed+1;

 

cout << "loading stack_common.h" << endl;

const double lumi = .182;


const double lumiZH = 3839.; 
const double lumiWH = 1321.; 

const double lumiZJ = 0.375; 
const double lumiWJ = 0.463; 


const double lumiWZ = 42. ;
const double lumiWW = 48.;
const double lumiZZ =357.; 

const double lumiTT = 5.; 

const double lumiTs = 499.; 
const double lumiTt = 23.; 
const double lumiTtW = 46.; 


const double mMin = 70;
const double mMax = 110;





TCanvas *c1 = new TCanvas("c1","Stack plot", 300,300,479,510);

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

 // histogram limits, in linear and logarithmic
  int nbin_(100);
  float xmin_(0.), xmax_(100.); 
  float ymin_(0.), ymax_(40.); 
  float yminl_(0.1), ymaxl_(200.); 

  // titles and axis, marker size
  TString xtitle;
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




void setHisto(TH1 * h, Color_t fill, Color_t line, double scale, int rebin) {
  h->SetFillColor(fill);
  h->SetLineColor(line);
  h->Scale(scale);
  h->Rebin(rebin);  

}



void stat(TH1 * hZH, TH1 * hWH, TH1 * hZJ, TH1 * hWJ,  TH1 * hZZ, TH1 * hWZ, TH1 * hWW,   TH1 * hTT, TH1 *hTs , TH1 *hTt , TH1 *hTtW ,  TH1* hdata, int rebin, bool scaleToDataInt ) 
{
  double a = mMin/rebin +1, b = mMax/rebin;


  double iZH = hZH->Integral(a, b);
  double errZH = sqrt(lumi/lumiZH *  iZH);

  double iWH = hWH->Integral(a, b);
  double errWH = sqrt(lumi/lumiWH *  iWH);

  double iZJ = hZJ->Integral(a, b);
  double errZJ = sqrt(lumi/lumiZJ *  iZJ);

  double iWJ = hWJ->Integral(a, b);
  double errWJ = sqrt(lumi/lumiWJ *  iWJ);

  double iZZ = hZZ->Integral(a, b);
  double errZZ = sqrt(lumi/lumiZZ *  iZZ);

  double iWW = hWW->Integral(a, b);
  double errWW = sqrt(lumi/lumiWW *  iWW);

  double iWZ = hWZ->Integral(a, b);
  double errWZ = sqrt(lumi/lumiWZ *  iWZ);

  double iTT = hTT->Integral(a, b);
  double errTT = sqrt(lumi/lumiTT *  iTT);

  double iTs = hTs->Integral(a, b);
  double errTs = sqrt(lumi/lumiTs *  iTs);

  double iTt = hTt->Integral(a, b);
  double errTt = sqrt(lumi/lumiTt *  iTt);

  double iTtW = hTtW->Integral(a, b);
  double errTtW = sqrt(lumi/lumiTtW *  iTtW);



  double idata = hdata != 0 ? hdata->Integral(a, b) : 0;
  //  double errData = sqrt(idata);

  std::cout.setf(0,ios::floatfield);
  std::cout.setf(ios::fixed,ios::floatfield);
  std::cout.precision(1);
  std::cout <<"ZH (" << mMin << ", " << mMax << ") = ";
  std::cout.precision(8);
  std::cout << iZH << "+/- " << errZH <<std::endl;
  std::cout.precision(1);
  std::cout <<"WH (" << mMin << ", " << mMax << ") = ";
  std::cout.precision(8);
  std::cout << iWH << "+/- " << errWH <<std::endl;
  std::cout.precision(1);
  std::cout <<"ZJ (" << mMin << ", " << mMax << ") = ";
  std::cout.precision(8);
  std::cout << iZJ << "+/- " << errZJ <<std::endl;
  std::cout.precision(1);
  std::cout <<"WJ(" << mMin << ", " << mMax << ") = "; 
  std::cout.precision(8);
  std::cout << iWJ << "+/- " << errWJ <<std::endl; 
  std::cout.precision(1);

  std::cout <<"WW(" << mMin << ", " << mMax << ") = "; 
  std::cout.precision(8);
  std::cout << iWW << "+/- " << errWW <<std::endl; 
  std::cout.precision(1);

  std::cout <<"WZ(" << mMin << ", " << mMax << ") = "; 
  std::cout.precision(8);
  std::cout << iWZ << "+/- " << errWZ <<std::endl; 
  std::cout.precision(1);

  std::cout <<"ZZ(" << mMin << ", " << mMax << ") = "; 
  std::cout.precision(8);
  std::cout << iZZ << "+/- " << errZZ <<std::endl; 
  std::cout.precision(1);

  std::cout <<"TT (" << mMin << ", " << mMax << ") = "; 
  std::cout.precision(8);
  std::cout << iTT << "+/- " << errTT <<std::endl; 
  std::cout.precision(1);
 std::cout <<"Ts(" << mMin << ", " << mMax << ") = "; 
  std::cout.precision(8);
  std::cout << iTs << "+/- " << errTs <<std::endl; 
  std::cout.precision(1);
 std::cout <<"Tt(" << mMin << ", " << mMax << ") = "; 
  std::cout.precision(8);
  std::cout << iTt << "+/- " << errTt <<std::endl; 
  std::cout.precision(1);
 std::cout <<"TtW(" << mMin << ", " << mMax << ") = "; 
  std::cout.precision(8);
  std::cout << iTtW << "+/- " << errTtW<<std::endl; 
  std::cout.precision(1);


  std::cout <<"data (" << mMin << ", " << mMax << ") = ";
  

//if one wnats one can scale any MC to the data events 
  if (scaleToDataInt) {
    double ratioToDataInt = idata / ( iZZ + iZW + iWW + iZJ + iWJ + iTT + iTs + iTt + iTtW ) ; 
    hZH->Scale( ratioToDataInt);
    hWH->Scale( ratioToDataInt);
    hZJ->Scale( ratioToDataInt);
    hWJ->Scale( ratioToDataInt);
    hZZ->Scale( ratioToDataInt);
    hWW->Scale( ratioToDataInt);
    hWZ->Scale( ratioToDataInt);
    hTT->Scale( ratioToDataInt);
    hTs->Scale( ratioToDataInt);
    hTt->Scale( ratioToDataInt);
    hTtW->Scale( ratioToDataInt);
  }
}

void makeStack(TH1 * hZH, TH1 * hWH, TH1 * hZJ, TH1 * hWJ, TH1 * hZZ, TH1 * hWZ, TH1 * hWW, TH1 * hTT, TH1 * hTs, TH1* hTt, TH1 * hTtW,  TH1 * hdata, double min, int rebin , bool doData, bool logScale, bool scaleToDataInt) {
  

  setHisto(hZH, hFillColor, hLineColor, lumi/lumiZH, rebin);
  setHisto(hWH, hFillColor, hLineColor, lumi/lumiWH, rebin);
  
  setHisto(hZJ, vjFillColor, vjLineColor, lumi/lumiZJ, rebin);
  setHisto(hWJ, vjFillColor, vjLineColor, lumi/lumiWJ, rebin);
  
  setHisto(hZZ, vvFillColor, vvLineColor, lumi/lumiZZ, rebin);
  setHisto(hWZ, vvFillColor, vvLineColor, lumi/lumiWZ, rebin);
  setHisto(hWW, vvFillColor, vvLineColor, lumi/lumiWW, rebin);

  setHisto(hTT, ttFillColor, ttLineColor, lumi/lumiTT, rebin);
  setHisto(hTs, ttFillColor, ttLineColor, lumi/lumiTs, rebin);
  setHisto(hTt, ttFillColor, ttLineColor, lumi/lumiTt, rebin);
  setHisto(hTtW, ttFillColor, ttLineColor, lumi/lumiTtW, rebin);



  if (doData) hdata->Rebin(rebin); 

  stat(hZH, hWH, hZJ, hWJ, hZZ,hWZ, hWW, hTT, hTs, hTt, hTtW, hdata, rebin,scaleToDataInt);

  hZH->Add(hWH);

  hZJ->Add(hWJ);

  hZZ->Add(hWZ);
  hZZ->Add(hWW);
  
  hTT->Add(hTs);
  hTT->Add(hTt);
  hTT->Add(hTtW);







  THStack * hs = new THStack("hs","");


 
   
    
  
  //hs->Add(h5);
  hs->Add(hTT);
  hs->Add(hZZ);
  hs->Add(hZJ);
  hs->Add(hZH);
  
  hs->Draw("HIST");
  if(hdata != 0) {
    
    
    /* TGraphAsymmErrors* dataGraph = (TGraphAsymmErrors*)hdata;
    dataGraph->SetMarkerStyle(kFullCircle);
    dataGraph->SetMarkerColor(kBlack);
    dataGraph->SetMarkerSize(markerSize);
    // Remove the horizontal bars (at Michael's request)
    double x_(0), y_(0);
    for( int ii=0; ii<dataGraph->GetN(); ii++ )
      {
	dataGraph->SetPointEXlow(ii,0);
	dataGraph->SetPointEXhigh(ii,0);
	dataGraph->GetPoint(ii,x_,y_ );
	if( y_==0 )
	  {
	    dataGraph->RemovePoint( ii );
	    ii--;
	  }  
      }
dataGraph->Draw("pesame");
    */
      
  if(hdata != 0) {
   hdata->SetMarkerStyle(kFullCircle);
       // hdata->SetMarkerStyle(9);
    hdata->SetMarkerSize(markerSize);
    hdata->SetMarkerColor(kBlack);
    hdata->SetLineWidth(lineWidth);
    hdata->SetLineColor(kBlack);
    //gStyle->SetErrorX(.5);
    gStyle->SetEndErrorSize(1);

    if(hdata != 0)  hdata->Draw("PE1SAME");
    hdata->GetXaxis()->SetLabelSize(0);
    hdata->GetYaxis()->SetLabelSize(0);
    hdata->GetXaxis()->SetNdivisions(ndivx);
    hdata->GetYaxis()->SetNdivisions(ndivy);
  }  
  hs->GetXaxis()->SetNdivisions(ndivx);
  hs->GetYaxis()->SetNdivisions(ndivy);
    // log plots, so the maximum should be one order of magnitude more...
    

   
  hs->SetMaximum( TMath::Max( hdata->GetMaximum(), hs->GetMaximum()) +10) ;
  //hs->SetMaximum( 29.8 )  ;
  if (logScale) {
    // hs->SetMaximum( pow(10 , -0.3 + int(log( hdata->GetMaximum() )  )  ));
    // hs->SetMaximum( pow(10 , +1.5 + int(log( hdata->GetMaximum() )  )  ));
    hs->SetMaximum( TMath::Max( hdata->GetMaximum(), hs->GetMaximum()) +100) ;
  } 
  
  
  }
  hs->SetMinimum(min);
  
  hs->GetXaxis()->SetTitle("M(#mu^{+} #mu^{-}) [GeV]");
  

  
  std::string yTag = "";
  switch(rebin) {
  case 1: yTag = "number of events/ 1 GeV"; break;
  case 2: yTag = "number of events/ 2 GeV"; break;
  case 2.5: yTag = "number of events/ 2.5 GeV"; break;
  case 3: yTag = "number of events/ 3 GeV"; break;
  case 4: yTag = "number of events/ 4 GeV"; break;
  case 5: yTag = "number of events/ 5 GeV"; break;
  case 10: yTag = "number of events/ 10 GeV"; break;
  case 20: yTag = "number of events/ 20 GeV"; break;
  default:
    std::cerr << ">>> ERROR: set y tag for rebin = " << rebin << std::endl;
  };
  
  hs->GetYaxis()->SetTitle(yTag.c_str());
  


   hs->GetXaxis()->SetTitleOffset(1.0);
   hs->GetYaxis()->SetTitleOffset(1.2);
   
   hs->GetYaxis()->SetLabelOffset(0.0);
   hs->GetXaxis()->SetLabelSize(.05);
   hs->GetYaxis()->SetLabelSize(.04);
   hs->GetYaxis()->SetTitleSize(.05);

  

// middle right
TLegend* legend=new TLegend(0.6,0.72,0.8,0.80);   

  if (logScale) { TLegend* legend=new TLegend(0.7,0.68,0.9,0.82);   }
  legend->SetLineColor(0);
  legend->SetFillColor(0);





  //leg = new TLegend(0.20,0.7,0.35,0.85);
  if(hdata != 0)
    if (doData) legend->AddEntry(hdata,"data", "pl");
  legend->AddEntry(hZH,"ZH+WH","f");
  legend->AddEntry(hZZ,"ZZ+WW+WZ","f");
  legend->AddEntry(hZJ,"Z/W + jets","f");
  legend->AddEntry(hTT,"t#bar{t} + st","f"); 

  


   
 legend->SetShadowColor(kWhite);
  legend->Draw();
 


TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.04);

  latex.SetTextAlign(31); // align right
  latex.DrawLatex(0.90,0.96,"#sqrt{s} = 7 TeV");
  if (lumi > 0.) {
    latex.SetTextAlign(31); // align right
    latex.DrawLatex(0.88,0.84,Form("#int #font[12]{L} dt = %.3f fb^{-1}",lumi));
  }
  latex.SetTextAlign(11); // align left
  latex.DrawLatex(0.12,0.96,"CMS preliminary 2011");


}






void makePlots(const char * var,    TCut cut,  int rebin, const char * plot, double min = 0.001, unsigned int nbins, double xMin, double xMax,  bool doData = true, bool logScale=false, bool scaleToDataInt = true) {
  
  
 
  
  TChain * zhEvents = new TChain("Events"); 
  zhEvents->Add("/data/1a/degrutto/data/Spring11/ZinvHbb/ZinvHbb115/ZinvHbb.root");
  
  TChain * whEvents = new TChain("Events"); 
  whEvents->Add("/data/1a/degrutto/data/Spring11/ZinvHbb/WlnHbb115/WlnHbb.root");
  
  TChain * zjEvents = new TChain("Events"); 
  zjEvents->Add("/data/1a/degrutto/data/Spring11/ZinvHbb/Zinvjets/Zinvjets.root");
  
  
  TChain * wjEvents = new TChain("Events"); 
  wjEvents->Add("/data/1a/degrutto/data/Spring11/ZinvHbb/Wjets/Wjets.root");
  
  TChain * zzEvents = new TChain("Events"); 
  zzEvents->Add("/data/1a/degrutto/data/Spring11/ZinvHbb/ZZ/ZZ.root");
 
  TChain * wwEvents = new TChain("Events"); 
  wwEvents->Add("/data/1a/degrutto/data/Spring11/ZinvHbb/WW/WW.root");

  TChain * wzEvents = new TChain("Events"); 
  wzEvents->Add("/data/1a/degrutto/data/Spring11/ZinvHbb/WZ/WZ.root");
  
  TChain * wwEvents = new TChain("Events"); 
  wwEvents->Add("/data/1a/degrutto/data/Spring11/ZinvHbb/WW/WW.root");
  
  TChain * ttEvents = new TChain("Events"); 
  ttEvents->Add("/data/1a/degrutto/data/Spring11/ZinvHbb/TTbar/ttbar.root");
  
  TChain * stsEvents = new TChain("Events"); 
  stsEvents->Add("/data/1a/degrutto/data/Spring11/ZinvHbb/T-s/Ts.root");
  
  TChain * sttEvents = new TChain("Events"); 
  sttEvents->Add("/data/1a/degrutto/data/Spring11/ZinvHbb/T-t/Tt.root");
  
  TChain * sttwEvents = new TChain("Events"); 
  sttwEvents->Add("/data/1a/degrutto/data/Spring11/ZinvHbb/T-tW/TtW.root");
  




  TChain * dataEvents2011= new TChain("Events");

  dataEvents2011->Add("/data/1a/degrutto/data/data11/METv2/METv2.root");
  //  ###dataEvents2011->Add("");
  




  TH1F *hZH = new TH1F ("hZH", "hZH", nbins, xMin, xMax);
  TH1F *hWH = new TH1F ("hWH", "hWH", nbins, xMin, xMax);
  
  TH1F *hZJ = new TH1F ("hZJ", "hZJ", nbins, xMin, xMax);
  TH1F *hWJ = new TH1F ("hWJ", "hWJ", nbins, xMin, xMax);

  TH1F *hZZ = new TH1F ("hZZ", "hZZ", nbins, xMin, xMax);
  TH1F *hWZ = new TH1F ("hWZ", "hWZ", nbins, xMin, xMax);
  TH1F *hWW = new TH1F ("hWW", "hWW", nbins, xMin, xMax);

  TH1F *hTT = new TH1F ("hTT", "hTT", nbins, xMin, xMax);
  TH1F *hTs = new TH1F ("hTs", "hTs", nbins, xMin, xMax);
  TH1F *hTt = new TH1F ("hTt", "hTt", nbins, xMin, xMax);
  TH1F *hTtW = new TH1F ("hTtW", "hTtW", nbins, xMin, xMax);


  zhEvents->Project("hZH", var, cut);
  whEvents->Project("hWH", var, cut);
  
  zjEvents->Project("hZJ", var, cut);
  wjEvents->Project("hWJ", var, cut);
  
  zzEvents->Project("hZZ", var, cut);
  wzEvents->Project("hWZ", var, cut);
  wwEvents->Project("hWW", var, cut);

  ttEvents->Project("hTT", var, cut);
  stsEvents->Project("hTs", var, cut);
  sttEvents->Project("hTt", var, cut);
  sttwEvents->Project("hTtW", var, cut);
  

  
 
  if(!doData) TH1F *hdata =  0;
  if (doData) { 
    TH1F *hdata = new TH1F ("hdata", "hdata", nbins, xMin, xMax);
    
    dataEvents2011->Project("hdata", var, cut) ;
    
  }
  
  makeStack(hZH, hWH, hZJ, hWJ, hZZ, hWZ, hWW, hTT, hTs, hTt, hTtW, hdata, min, rebin, doData, logScale, scaleToDataInt);
 
  
  

 
  if (logScale) c1->SetLogy();

  c1->SaveAs((std::string(plot)+".eps").c_str());
  c1->SaveAs((std::string(plot)+".gif").c_str());
  c1->SaveAs((std::string(plot)+".pdf").c_str());
  
  TFile * out = new TFile("plot.root", "RECREATE");

  c1->Write();
  c1->SaveAs("hPlot.C");
}


#endif
