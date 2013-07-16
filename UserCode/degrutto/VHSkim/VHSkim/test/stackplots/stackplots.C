#include <sstream>
#include <string>
#include "stack_common.h"


double deltaPhi(double phi1, double phi2) { 
  double result = phi1 - phi2;
  while (result > TMath::Pi()) result -= 2* TMath::Pi();
  while (result <= -TMath::Pi()) result += 2* TMath::Pi();
  return result;
}

double deltaR(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = deltaPhi(phi1, phi2);
  return sqrt(deta*deta + dphi*dphi);
}



void stackplots() {

  
   makePlots("pfMetEt",  "hjjJet1CSV>0.5 && hjjJet2CSV>0.5",  20, "pfMET",  0.1, 500, 0 ,500, true, true, false);

   hs->GetXaxis()->SetTitle("pfMET [GeV]");
   //hs->GetYaxis()->SetTitle("Events / 20 ");
   // hs->SetMaximum( 1005);
   c1->SaveAs("met.png");




}


