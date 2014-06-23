/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Gaussian core + exponential tail on high side + exponential tail on low side
// Souvik Das
// 8/23/2013

#include "Riostream.h" 

#include "RooExpGaussExp.h" 
//#include "RooAbsReal.h" 
//#include "RooAbsCategory.h" 
//#include <math.h> 
//#include "TMath.h" 

ClassImp(RooExpGaussExp) 

RooExpGaussExp::RooExpGaussExp(const char *name, const char *title, 
                       RooAbsReal& _x,
                       RooAbsReal& _mean,
                       RooAbsReal& _sigma,
                       RooAbsReal& _lambda,
                       RooAbsReal& _kappa) :
  RooAbsPdf(name,title), 
  x("x","x",this,_x),
  mean("mean","mean",this,_mean),
  sigma("sigma","sigma",this,_sigma),
  lambda("lambda","lambda",this,_lambda),
  kappa("kappa","kappa",this,_kappa)
{ 
} 


RooExpGaussExp::RooExpGaussExp(const RooExpGaussExp& other, const char* name) :  
  RooAbsPdf(other,name), 
  x("x",this,other.x),
  mean("mean",this,other.mean),
  sigma("sigma",this,other.sigma),
  lambda("lambda",this,other.lambda),
  kappa("kappa",this,other.kappa)
{ 
} 

Double_t RooExpGaussExp::evaluate() const 
{
  Double_t std = (x-mean)/sigma;
  Double_t result = 0;
  assert (lambda >= kappa);
  
  if (std<kappa) {
    result = exp(kappa*kappa/2.-kappa*std);
  } else if (kappa<=std && std<lambda) {
    result = exp(-std*std/2.);
  } else {
    result = exp(lambda*lambda/2.-lambda*std);
  }
  
  return result; 
} 


////_____________________________________________________________________________
//Int_t RooExpGaussExp::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
//{
//  if (matchArgs(allVars,analVars,x)) return 1 ;
//  if (matchArgs(allVars,analVars,mean)) return 2 ;
//  return 0 ;
//}
//
//
////_____________________________________________________________________________
//Double_t RooExpGaussExp::analyticalIntegral(Int_t code, const char* rangeName) const
//{
//  assert(code==1 || code==2) ;
//
//  static Double_t root2 = sqrt(2.) ;
//  static Double_t rootPiBy2 = sqrt(atan2(0.0,-1.0)/2.0);
//  Double_t xscale = root2*sigma;
//  Double_t stdmax = (x.max(rangeName)-mean)/sigma;
//  Double_t stdmin = (x.min(rangeName)-mean)/sigma;
//  Double_t ret = 0;
//  assert (lambda >= kappa);
//
//  if(code==1){
//    if (stdmax<kappa) {  // stdmin < stdmax < kappa < lambda
//      if(kappa == 0.0) {
//        ret = (x.max(rangeName) - x.min(rangeName));
//      } else {
//        ret = ( exp( kappa*kappa/2.-kappa*stdmax ) - exp( kappa*kappa/2.-kappa*stdmin ) )*(-sigma/kappa);
//      }
//cout << "1" << kappa << " " << lambda << " " << mean << " " << sigma << " " << stdmin << " " << stdmax << " " << ret << endl;
//
//    } else if (stdmin>lambda) {  // kappa < lambda < stdmin < stdmax
//      if(lambda == 0.0) {
//        ret = (x.max(rangeName) - x.min(rangeName));
//      } else {
//        ret = ( exp( lambda*lambda/2.-lambda*stdmax ) - exp( lambda*lambda/2.-lambda*stdmin ) )*(-sigma/lambda);
//      }
//cout << "2" << kappa << " " << lambda << " " << mean << " " << sigma << " " << stdmin << " " << stdmax << " " << ret << endl;
//
//    } else if (kappa<stdmin && stdmax<lambda) {  // kappa < stdmin < stdmax < lambda
//      ret = rootPiBy2*sigma*(RooMath::erf(stdmax/root2)-RooMath::erf(stdmin/root2));
//cout << "3" << kappa << " " << lambda << " " << mean << " " << sigma << " " << stdmin << " " << stdmax << " " << ret << endl;
//
//    } else if (kappa<stdmin) {  // kappa < stdmin < lambda < stdmax
//      ret = rootPiBy2*sigma*(RooMath::erf(lambda/root2)-RooMath::erf(stdmin/root2));
//      if(lambda == 0.0) {
//        ret += (x.max(rangeName) - 0.0);
//      } else {
//        ret += ( exp( lambda*lambda/2.-lambda*stdmax ) - exp( -lambda*lambda/2. ) )/(-sigma/lambda);
//      }
//cout << "4" << kappa << " " << lambda << " " << mean << " " << sigma << " " << stdmin << " " << stdmax << " " << ret << endl;
//
//    } else if (stdmax<lambda) {  // stdmin < kappa < stdmax < lambda
//      ret = rootPiBy2*sigma*(RooMath::erf(stdmax/root2)-RooMath::erf(kappa/root2));
//      if(kappa == 0.0) {
//        ret += (0.0 - x.min(rangeName));
//      } else {
//        ret += ( exp( -kappa*kappa/2. ) - exp( kappa*kappa/2.-kappa*stdmin ) )*(-sigma/kappa);
//      }
//cout << "5" << kappa << " " << lambda << " " << mean << " " << sigma << " " << stdmin << " " << stdmax << " " << ret << endl;
//cout << "  " << rootPiBy2*sigma*(RooMath::erf(stdmax/root2)-RooMath::erf(kappa/root2)) << " " << ( exp( -kappa*kappa/2. ) - exp( kappa*kappa/2.-kappa*stdmin ) )*(-sigma/kappa) << endl;
//
//    } else {  // stdmin < kappa < lambda < stdmax
//      assert(stdmin <= kappa && lambda <= stdmax);
//      ret = rootPiBy2*sigma*(RooMath::erf(lambda/root2)-RooMath::erf(kappa/root2));
//      if(kappa == 0.0) {
//        ret += (0.0 - x.min(rangeName));
//      } else {
//        ret += ( exp( -kappa*kappa/2. ) - exp( kappa*kappa/2.-kappa*stdmin ) )*(-sigma/kappa);
//      }
//      if(lambda == 0.0) {
//        ret += (x.max(rangeName) - 0.0);
//      } else {
//        ret += ( exp( lambda*lambda/2.-lambda*stdmax ) - exp( -lambda*lambda/2. ) )/(-sigma/lambda);
//      }
//      //cout << "stdmin="  << stdmin << ", stdmax=" << stdmax << ", lambda=" << lambda << ", mean=" << mean << ", sigma=" << sigma << ", ret=" << ret <<  endl;
//cout << "6" << kappa << " " << lambda << " " << mean << " " << sigma << " " << stdmin << " " << stdmax << " " << ret << endl;
//cout << "  " << rootPiBy2*sigma*(RooMath::erf(lambda/root2)-RooMath::erf(kappa/root2)) << " " << ( exp( -kappa*kappa/2. ) - exp( kappa*kappa/2.-kappa*stdmin ) )*(-sigma/kappa) << " " << ( exp( lambda*lambda/2.-lambda*stdmax ) - exp( -lambda*lambda/2. ) )/(-sigma/lambda) << endl;
//
//    }
//
//  } else if(code==2) {
//    cout << "Error in RooExpGaussExp::analyticalIntegral: Not implemented " << endl;
//
//  } else {
//    cout << "Error in RooExpGaussExp::analyticalIntegral" << endl;
//  }
//  //cout << ret << endl;
//  return ret ;
//}

