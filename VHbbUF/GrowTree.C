#include <cassert>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TFile.h"
#include "TString.h"
#include "TBranch.h"
#include "TBranchElement.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include "TMath.h"
//#include "TPRegexp.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TTree.h"
#include "TTreeFormula.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#endif

#include "HelperTMVA.h"
#include "HelperFunctions.h"
#include "HelperNtuples.h"
#include "HelperVHbbDataFormats.h"

#include "EquationSolver.h"
#include "Math/GenVector/LorentzVector.h"
namespace math {
    typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > XYZTLorentzVector;
}

#define MAXJ 130
#define MAXL 110

//#define CSVSYST
#ifdef CSVSYST
#include "BTagReshaping.h"  // FIXME: local version is out of date with V6 Step 2
#endif

//#define JECFWLITE
#ifdef JECFWLITE
#include "JECFWLiteStandalone.h"
#endif

#define STITCH

// FIXME: also patch V.pt
#define PATCHMETTYPE1CORR


// FIXME: change it to output the delta, not met+delta
float shift_met_by_Nvtx(float met, float metphi, int Nvtx, int EVENT_run)
{
    // from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/JetMETCorrections/Type1MET/python/pfMETsysShiftCorrections_cfi.py?revision=1.6&view=markup
    float mex = met * cos(metphi);
    float mey = met * sin(metphi);
    float px = 0.0, py = 0.0;
    if (EVENT_run != 1) {
        //pfMEtSysShiftCorrParameters_2012runABCvsNvtx_data
        px = +0.2661 + 0.3217*Nvtx;
        py = -0.2251 - 0.1747*Nvtx;
    } else {
        //pfMEtSysShiftCorrParameters_2012runABCvsNvtx_mc
        px = +0.1166 + 0.0200*Nvtx;
        py = +0.2764 - 0.1280*Nvtx;
    }
    mex -= px;
    mey -= py;
    return TMath::Sqrt(mex * mex + mey * mey);
}

// FIXME: change it to output the delta, not met+delta
float shift_metphi_by_Nvtx(float met, float metphi, int Nvtx, int EVENT_run)
{
    // from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/JetMETCorrections/Type1MET/python/pfMETsysShiftCorrections_cfi.py?revision=1.6&view=markup
    float mex = met * cos(metphi);
    float mey = met * sin(metphi);
    float px = 0.0, py = 0.0;
    if (EVENT_run != 1) {
        //pfMEtSysShiftCorrParameters_2012runABCvsNvtx_data
        px = +0.2661 + 0.3217 * Nvtx;
        py = -0.2251 - 0.1747 * Nvtx;
    } else {
        //pfMEtSysShiftCorrParameters_2012runABCvsNvtx_mc
        px = +0.1166 + 0.0200 * Nvtx;
        py = +0.2764 - 0.1280 * Nvtx;
    }
    mex -= px;
    mey -= py;
    if (mex == 0.0 && mey == 0.0)
        return 0.0;
    float phi1 = TMath::ATan2(mey, mex);
    float phi2 = TMath::ATan2(mey, mex) - 2.0 * M_PI;
    if (TMath::Abs(phi1 - metphi) < TMath::Abs(phi2 - metphi) + 0.5 * M_PI)
        return phi1;
    else
        return phi2;
}

TLorentzVector makePtEtaPhiE(double ptCor, double pt, double eta, double phi, double e, 
                             bool rescaleEnergy=true)
{
    TLorentzVector j;
    j.SetPtEtaPhiE(ptCor, eta, phi, (rescaleEnergy ? (e * ptCor / pt) : e));
    return j;
}

// Copied from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/TopQuarkAnalysis/SingleTop/src/TopProducer.cc?revision=1.9&view=markup
TLorentzVector getNu4Momentum(const TLorentzVector& TLepton, const TLorentzVector& TMET)
{
  const math::XYZTLorentzVector Lepton(TLepton.Px(), TLepton.Py(), TLepton.Pz(), TLepton.E());
  const math::XYZTLorentzVector MET(TMET.Px(), TMET.Py(), 0., TMET.E());

  double  mW = 80.38;

  //std::vector<math::XYZTLorentzVector> result;
  std::vector<TLorentzVector> result;

  //  double Wmt = sqrt(pow(Lepton.et()+MET.pt(),2) - pow(Lepton.px()+MET.px(),2) - pow(Lepton.py()+MET.py(),2) );
    
  double MisET2 = (MET.px()*MET.px() + MET.py()*MET.py());
  double mu = (mW*mW)/2 + MET.px()*Lepton.px() + MET.py()*Lepton.py();
  double a  = (mu*Lepton.pz())/(Lepton.energy()*Lepton.energy() - Lepton.pz()*Lepton.pz());
  double a2 = TMath::Power(a,2);
  double b  = (TMath::Power(Lepton.energy(),2.)*(MisET2) - TMath::Power(mu,2.))/(TMath::Power(Lepton.energy(),2) - TMath::Power(Lepton.pz(),2));
  double pz1(0),pz2(0),pznu(0);
  int nNuSol(0);

  //math::XYZTLorentzVector p4nu_rec;
  TLorentzVector p4nu_rec;
  math::XYZTLorentzVector p4W_rec;
  math::XYZTLorentzVector p4b_rec;
  math::XYZTLorentzVector p4Top_rec;
  math::XYZTLorentzVector p4lep_rec;

  p4lep_rec.SetPxPyPzE(Lepton.px(),Lepton.py(),Lepton.pz(),Lepton.energy());
  
  //math::XYZTLorentzVector p40_rec(0,0,0,0);

  if(a2-b > 0 ){
    //if(!usePositiveDeltaSolutions_)
    //  {
    //    result.push_back(p40_rec);
    //    return result;
    //  }
    double root = sqrt(a2-b);
    pz1 = a + root;
    pz2 = a - root;
    nNuSol = 2;

    //if(usePzPlusSolutions_)pznu = pz1;    
    //if(usePzMinusSolutions_)pznu = pz2;
    //if(usePzAbsValMinimumSolutions_){
      pznu = pz1;
      if(fabs(pz1)>fabs(pz2)) pznu = pz2;
    //}

    double Enu = sqrt(MisET2 + pznu*pznu);

    p4nu_rec.SetPxPyPzE(MET.px(), MET.py(), pznu, Enu);

    result.push_back(p4nu_rec);

  }
  else{
    //if(!useNegativeDeltaSolutions_){
    //  result.push_back(p40_rec);
    //  return result;
    //}
    //    double xprime = sqrt(mW;

    double ptlep = Lepton.pt(),pxlep=Lepton.px(),pylep=Lepton.py(),metpx=MET.px(),metpy=MET.py();

    double EquationA = 1;
    double EquationB = -3*pylep*mW/(ptlep);
    double EquationC = mW*mW*(2*pylep*pylep)/(ptlep*ptlep)+mW*mW-4*pxlep*pxlep*pxlep*metpx/(ptlep*ptlep)-4*pxlep*pxlep*pylep*metpy/(ptlep*ptlep);
    double EquationD = 4*pxlep*pxlep*mW*metpy/(ptlep)-pylep*mW*mW*mW/ptlep;

    std::vector<long double> solutions = EquationSolve<long double>((long double)EquationA,(long double)EquationB,(long double)EquationC,(long double)EquationD);

    std::vector<long double> solutions2 = EquationSolve<long double>((long double)EquationA,-(long double)EquationB,(long double)EquationC,-(long double)EquationD);

    double deltaMin = 14000*14000;
    double zeroValue = -mW*mW/(4*pxlep); 
    double minPx=0;
    double minPy=0;

    //    std::cout<<"a "<<EquationA << " b " << EquationB  <<" c "<< EquationC <<" d "<< EquationD << std::endl; 
      
    //if(usePxMinusSolutions_){
      for( int i =0; i< (int)solutions.size();++i){
      if(solutions[i]<0 ) continue;
      double p_x = (solutions[i]*solutions[i]-mW*mW)/(4*pxlep); 
      double p_y = ( mW*mW*pylep + 2*pxlep*pylep*p_x -mW*ptlep*solutions[i])/(2*pxlep*pxlep);
      double Delta2 = (p_x-metpx)*(p_x-metpx)+(p_y-metpy)*(p_y-metpy); 

      //      std::cout<<"intermediate solution1 met x "<<metpx << " min px " << p_x  <<" met y "<<metpy <<" min py "<< p_y << std::endl; 

      if(Delta2< deltaMin && Delta2 > 0){deltaMin = Delta2;
      minPx=p_x;
      minPy=p_y;}
      //     std::cout<<"solution1 met x "<<metpx << " min px " << minPx  <<" met y "<<metpy <<" min py "<< minPy << std::endl; 
      }
    //}

    //if(usePxPlusSolutions_){
      for( int i =0; i< (int)solutions2.size();++i){
        if(solutions2[i]<0 ) continue;
        double p_x = (solutions2[i]*solutions2[i]-mW*mW)/(4*pxlep); 
        double p_y = ( mW*mW*pylep + 2*pxlep*pylep*p_x +mW*ptlep*solutions2[i])/(2*pxlep*pxlep);
        double Delta2 = (p_x-metpx)*(p_x-metpx)+(p_y-metpy)*(p_y-metpy); 
        //  std::cout<<"intermediate solution2 met x "<<metpx << " min px " << minPx  <<" met y "<<metpy <<" min py "<< minPy << std::endl; 
        if(Delta2< deltaMin && Delta2 > 0){deltaMin = Delta2;
          minPx=p_x;
          minPy=p_y;
        }
        //        std::cout<<"solution2 met x "<<metpx << " min px " << minPx  <<" met y "<<metpy <<" min py "<< minPy << std::endl; 
      }
    //}

    double pyZeroValue= ( mW*mW*pxlep + 2*pxlep*pylep*zeroValue);
    double delta2ZeroValue= (zeroValue-metpx)*(zeroValue-metpx) + (pyZeroValue-metpy)*(pyZeroValue-metpy);

    if(deltaMin==14000*14000) return TLorentzVector(0,0,0,0);
    //if(deltaMin==14000*14000) return result.front();
    //    else std::cout << " test " << std::endl;

    if(delta2ZeroValue < deltaMin){
      deltaMin = delta2ZeroValue;
      minPx=zeroValue;
      minPy=pyZeroValue;}

    //    std::cout<<" MtW2 from min py and min px "<< sqrt((minPy*minPy+minPx*minPx))*ptlep*2 -2*(pxlep*minPx + pylep*minPy)  <<std::endl;
    ///    ////Y part   

    double mu_Minimum = (mW*mW)/2 + minPx*pxlep + minPy*pylep;
    double a_Minimum  = (mu_Minimum*Lepton.pz())/(Lepton.energy()*Lepton.energy() - Lepton.pz()*Lepton.pz());
    pznu = a_Minimum;

    //if(!useMetForNegativeSolutions_){
      double Enu = sqrt(minPx*minPx+minPy*minPy + pznu*pznu);
      p4nu_rec.SetPxPyPzE(minPx, minPy, pznu , Enu);
    //}
    //else{
    //  pznu = a;
    //  double Enu = sqrt(metpx*metpx+metpy*metpy + pznu*pznu);
    //  p4nu_rec.SetPxPyPzE(metpx, metpy, pznu , Enu);
    //}
    result.push_back(p4nu_rec);
  }
  return result.front();
}

void ResetDeleteBranches(TTree* tree)
{
    TObjArray* branches = tree->GetListOfBranches();
    Int_t nb = branches->GetEntriesFast();
    for (Int_t i = 0; i < nb; ++i) {
        TBranch* br = (TBranch*) branches->UncheckedAt(i);
        if (br->InheritsFrom(TBranchElement::Class())) {
            ((TBranchElement*) br)->ResetDeleteObject();
        }
    }
}


////////////////////////////////////////////////////////////////////////////////
/// Main                                                                     ///
////////////////////////////////////////////////////////////////////////////////
void GrowTree(TString process, std::string regMethod="BDTG", Long64_t beginEntry=0, Long64_t endEntry=-1)
{
    gROOT->SetBatch(1);
    TH1::SetDefaultSumw2(1);
    gROOT->LoadMacro("HelperFunctions.h");  //< make functions visible to TTreeFormula

    if (!TString(gROOT->GetVersion()).Contains("5.34")) {
        std::cout << "INCORRECT ROOT VERSION! Please use 5.34:" << std::endl;
        std::cout << "source /uscmst1/prod/sw/cms/slc5_amd64_gcc462/lcg/root/5.34.02-cms/bin/thisroot.csh" << std::endl;
        std::cout << "Return without doing anything." << std::endl;
        return;
    }

    //const TString indir   = "skim_ZnnH_baseline/";
    const TString indir   = "dcache:/pnfs/cms/WAX/resilient/jiafu/ZnunuHbb/skim_ZnnH_baseline/";
    const TString outdir  = "skim/";
    const TString prefix  = "skim_";
    const TString suffix  = ".root";

    TFile *input = TFile::Open(indir + prefix + process + suffix);
    if (!input) {
        std::cout << "ERROR: Could not open input file." << std::endl;
        exit(1);
    }
    /// Make output directory if it doesn't exist
    if (gSystem->AccessPathName(outdir))
        gSystem->mkdir(outdir);

    std::cout << "--- GrowTree                 : Using input file: " << input->GetName() << std::endl;
    TTree *inTree = (TTree *) input->Get("tree");
    TH1F  *hcount = (TH1F *) input->Get("Count");
    TFile *output(0);
    if (beginEntry == 0 && endEntry == -1)
        output = TFile::Open(outdir + "Step3_" + process + suffix, "RECREATE");
    else
        output = TFile::Open(outdir + "Step3_" + process + TString::Format("_%Li_%Li", beginEntry, endEntry) + suffix, "RECREATE");
    TTree *outTree = inTree->CloneTree(0); // Do no copy the data yet
    /// The clone should not delete any shared i/o buffers.
    ResetDeleteBranches(outTree);


    ///-- Set branch addresses -------------------------------------------------
    EventInfo EVENT;
    HiggsInfo H;
    FatHiggsInfo FatH;
    VInfo V;
    METInfo METtype1corr;
    METInfo METtype1diff;
    int Vtype, nPVs;
    int nhJets, naJets;
    float hJet_pt[2], hJet_eta[2], hJet_phi[2], hJet_e[2], hJet_ptRaw[2], 
        hJet_csv[2], hJet_csv_nominal[2], hJet_genPt[2], hJet_flavour[2], hJet_JECUnc[2],
        hJet_puJetIdL[2];
    float aJet_pt[MAXJ], aJet_eta[MAXJ], aJet_phi[MAXJ], aJet_e[MAXJ], 
        aJet_csv[MAXJ], aJet_csv_nominal[MAXJ], aJet_genPt[MAXJ], aJet_flavour[MAXJ], aJet_JECUnc[MAXJ],
        aJet_puJetIdL[MAXJ];
    bool hJet_id[2];
    bool aJet_id[MAXJ]; 
    int nfathFilterJets, naJetsFat;
    float fathFilterJets_pt[3], fathFilterJets_eta[3], fathFilterJets_phi[3], 
        fathFilterJets_e[3], fathFilterJets_ptRaw[3], fathFilterJets_csv[3], 
        fathFilterJets_csvivf[3], fathFilterJets_genPt[3], fathFilterJets_flavour[3], 
        fathFilterJets_JECUnc[3];
    float aJetFat_pt[MAXJ], aJetFat_eta[MAXJ], aJetFat_phi[MAXJ], aJetFat_e[MAXJ];
    float metUnc_et[24], metUnc_phi[24];  // could be changed to 12 in the future
    int nvlep, nalep;
    float vLepton_mass[2], vLepton_pt[2], vLepton_eta[2], vLepton_phi[2], vLepton_pfCombRelIso[2], vLepton_id95[2], vLepton_vbtf[2];
    int vLepton_type[2];
    float aLepton_mass[MAXL], aLepton_pt[MAXL], aLepton_eta[MAXL], aLepton_phi[MAXL], aLepton_pfCombRelIso[MAXL], aLepton_id95[MAXL], aLepton_vbtf[MAXL];
    int aLepton_type[MAXL];
    
    inTree->SetBranchAddress("EVENT", &EVENT);
    inTree->SetBranchAddress("H", &H);
    inTree->SetBranchAddress("FatH", &FatH);
    inTree->SetBranchAddress("V", &V);
    inTree->SetBranchAddress("METtype1corr", &METtype1corr);
    inTree->SetBranchAddress("METtype1diff", &METtype1diff);
    inTree->SetBranchAddress("Vtype", &Vtype);
    inTree->SetBranchAddress("nPVs", &nPVs);
    inTree->SetBranchAddress("nhJets", &nhJets);
    inTree->SetBranchAddress("naJets", &naJets);
    inTree->SetBranchAddress("naJetsFat", &naJetsFat);
    inTree->SetBranchAddress("hJet_pt", &hJet_pt);
    inTree->SetBranchAddress("hJet_eta", &hJet_eta);
    inTree->SetBranchAddress("hJet_phi", &hJet_phi);
    inTree->SetBranchAddress("hJet_e", &hJet_e);
    inTree->SetBranchAddress("hJet_ptRaw", &hJet_ptRaw);
    inTree->SetBranchAddress("hJet_csv", &hJet_csv);
    inTree->SetBranchAddress("hJet_csv_nominal", &hJet_csv_nominal);
    inTree->SetBranchAddress("hJet_genPt", &hJet_genPt);
    inTree->SetBranchAddress("hJet_flavour", &hJet_flavour);
    inTree->SetBranchAddress("hJet_JECUnc", &hJet_JECUnc);
    inTree->SetBranchAddress("hJet_puJetIdL", &hJet_puJetIdL);
    inTree->SetBranchAddress("aJet_pt", &aJet_pt);
    inTree->SetBranchAddress("aJet_eta", &aJet_eta);
    inTree->SetBranchAddress("aJet_phi", &aJet_phi);
    inTree->SetBranchAddress("aJet_e", &aJet_e);
    inTree->SetBranchAddress("aJet_csv", &aJet_csv);
    inTree->SetBranchAddress("aJet_csv_nominal", &aJet_csv_nominal);
    inTree->SetBranchAddress("aJet_genPt", &aJet_genPt);
    inTree->SetBranchAddress("aJet_flavour", &aJet_flavour);
    inTree->SetBranchAddress("aJet_JECUnc", &aJet_JECUnc);
    inTree->SetBranchAddress("aJet_puJetIdL", &aJet_puJetIdL);
    inTree->SetBranchAddress("hJet_id", &hJet_id);
    inTree->SetBranchAddress("aJet_id", &aJet_id);
    inTree->SetBranchAddress("nfathFilterJets", &nfathFilterJets);
    inTree->SetBranchAddress("fathFilterJets_pt", &fathFilterJets_pt);
    inTree->SetBranchAddress("fathFilterJets_eta", &fathFilterJets_eta);
    inTree->SetBranchAddress("fathFilterJets_phi", &fathFilterJets_phi);
    inTree->SetBranchAddress("fathFilterJets_e", &fathFilterJets_e);
    inTree->SetBranchAddress("fathFilterJets_ptRaw", &fathFilterJets_ptRaw);
    inTree->SetBranchAddress("fathFilterJets_csv", &fathFilterJets_csv);
    inTree->SetBranchAddress("fathFilterJets_csvivf", &fathFilterJets_csvivf);
    inTree->SetBranchAddress("fathFilterJets_genPt", &fathFilterJets_genPt);
    inTree->SetBranchAddress("fathFilterJets_flavour", &fathFilterJets_flavour);
    inTree->SetBranchAddress("fathFilterJets_JECUnc", &fathFilterJets_JECUnc);
    inTree->SetBranchAddress("aJetFat_pt", &aJetFat_pt);
    inTree->SetBranchAddress("aJetFat_eta", &aJetFat_eta);
    inTree->SetBranchAddress("aJetFat_phi", &aJetFat_phi);
    inTree->SetBranchAddress("aJetFat_e", &aJetFat_e);
    inTree->SetBranchAddress("metUnc_et", &metUnc_et);
    inTree->SetBranchAddress("metUnc_phi", &metUnc_phi);
    inTree->SetBranchAddress("nvlep", &nvlep);
    inTree->SetBranchAddress("nalep", &nalep);
    inTree->SetBranchAddress("vLepton_mass", &vLepton_mass);
    inTree->SetBranchAddress("vLepton_pt", &vLepton_pt);
    inTree->SetBranchAddress("vLepton_eta", &vLepton_eta);
    inTree->SetBranchAddress("vLepton_phi", &vLepton_phi);
    inTree->SetBranchAddress("vLepton_pfCombRelIso", &vLepton_pfCombRelIso);
    inTree->SetBranchAddress("vLepton_id95", &vLepton_id95);
    inTree->SetBranchAddress("vLepton_vbtf", &vLepton_vbtf);
    inTree->SetBranchAddress("vLepton_type", &vLepton_type);
    inTree->SetBranchAddress("aLepton_mass", &aLepton_mass);
    inTree->SetBranchAddress("aLepton_pt", &aLepton_pt);
    inTree->SetBranchAddress("aLepton_eta", &aLepton_eta);
    inTree->SetBranchAddress("aLepton_phi", &aLepton_phi);
    inTree->SetBranchAddress("aLepton_pfCombRelIso", &aLepton_pfCombRelIso);
    inTree->SetBranchAddress("aLepton_id95", &aLepton_id95);
    inTree->SetBranchAddress("aLepton_vbtf", &aLepton_vbtf);
    inTree->SetBranchAddress("aLepton_type", &aLepton_type);

    /// Increase reading speed
    //inTree->SetBranchStatus("*", 0);
    //inTree->SetBranchStatus("EVENT", 1);
    //inTree->SetBranchStatus("H", 1);
    //inTree->SetBranchStatus("V", 1);
    //inTree->SetBranchStatus("METtype1corr", 1);
    //inTree->SetBranchStatus("nPVs", 1);
    //inTree->SetBranchStatus("naJets", 1);
    //inTree->SetBranchStatus("hJet_*", 1);
    //inTree->SetBranchStatus("aJet_*", 1);
    //inTree->SetBranchStatus("metUnc_et", 1);
    //inTree->SetBranchStatus("metUnc_phi", 1);
    //inTree->SetBranchStatus("rho25", 1);
    inTree->SetBranchStatus("SimBs_*", 0);
    inTree->SetBranchStatus("Sv_*", 0);
    inTree->SetBranchStatus("TkSharing", 0);
    inTree->SetBranchStatus("weightTrigMET*", 0);
    inTree->SetBranchStatus("weightTrig2CJet20", 0);

    ///-- Make new branches ----------------------------------------------------
    int EVENT_run, EVENT_event;  // set these as TTree index?
    float lumi_ = lumi, efflumi, efflumi_old, 
        efflumi_UEPS_up, efflumi_UEPS_down;
    float hJet_ptReg[2], 
        hJet_ptReg_res_j_up[2], hJet_ptReg_res_j_down[2], hJet_ptReg_scale_j_up[2], hJet_ptReg_scale_j_down[2];
#ifdef CSVSYST
    float hJet_csv_reshaped[2], 
        hJet_csv_eff_b_up[2], hJet_csv_eff_b_down[2], hJet_csv_fake_b_up[2], hJet_csv_fake_b_down[2];
#endif
    float HptNorm, HptGen, HptReg,
        HptReg_res_j_up, HptReg_res_j_down, HptReg_scale_j_up, HptReg_scale_j_down;
    float HmassNorm, HmassGen, HmassReg,
        HmassReg_res_j_up, HmassReg_res_j_down, HmassReg_scale_j_up, HmassReg_scale_j_down;
    float aJet_pt_res_j_up[MAXJ], aJet_pt_res_j_down[MAXJ], aJet_pt_scale_j_up[MAXJ], aJet_pt_scale_j_down[MAXJ];
    float fathFilterJets_ptReg[3], 
        fathFilterJets_ptReg_res_j_up[3], fathFilterJets_ptReg_res_j_down[3], fathFilterJets_ptReg_scale_j_up[3], fathFilterJets_ptReg_scale_j_down[3];
    float FatHptNorm, FatHptGen, FatHptReg,
        FatHptReg_res_j_up, FatHptReg_res_j_down, FatHptReg_scale_j_up, FatHptReg_scale_j_down;
    float FatHmassNorm, FatHmassGen, FatHmassReg,
        FatHmassReg_res_j_up, FatHmassReg_res_j_down, FatHmassReg_scale_j_up, FatHmassReg_scale_j_down;
#ifdef CSVSYST
    float fathFilterJets_csv_reshaped[3];
#endif
    float MET_xyshift,
        MET_res_j_up, MET_res_j_down, MET_scale_j_up, MET_scale_j_down;
    float METphi_xyshift,
        METphi_res_j_up, METphi_res_j_down, METphi_scale_j_up, METphi_scale_j_down;
    float HTeta25, HTeta45;
    
    int nalep_Znn, nalep_pt5_Znn; 
    int naJets_Znn, 
        naJets_Znn_res_j_up, naJets_Znn_res_j_down, naJets_Znn_scale_j_up, naJets_Znn_scale_j_down; 
    int naCtrJets_pt20to25_Znn, 
        naCtrJets_pt20to25_Znn_res_j_up, naCtrJets_pt20to25_Znn_res_j_down, naCtrJets_pt20to25_Znn_scale_j_up, naCtrJets_pt20to25_Znn_scale_j_down;
    float mindPhiMETJet_pt;
    float mindPhiMETJet_dPhi, // absolute values
        mindPhiMETJet_dPhi_res_j_up, mindPhiMETJet_dPhi_res_j_down, mindPhiMETJet_dPhi_scale_j_up, mindPhiMETJet_dPhi_scale_j_down;
    float mindPhiMETCtrJet_pt;
    float mindPhiMETCtrJet_dPhi, // absolute values
        mindPhiMETCtrJet_dPhi_res_j_up, mindPhiMETCtrJet_dPhi_res_j_down, mindPhiMETCtrJet_dPhi_scale_j_up, mindPhiMETCtrJet_dPhi_scale_j_down;
    float mindRHJet_pt, 
        mindRHJet_pt_res_j_up, mindRHJet_pt_res_j_down, mindRHJet_pt_scale_j_up, mindRHJet_pt_scale_j_down;
    float mindRHJet_dR, 
        mindRHJet_dR_res_j_up, mindRHJet_dR_res_j_down, mindRHJet_dR_scale_j_up, mindRHJet_dR_scale_j_down;
    float mindRHJet_csv_nominal;
    float weightTrigMu, weightTrigDLP_ElecMay10, weightTrigDLP_Elec;
    float TopLep_mass, TopLep_pt, TopLep_eta, TopLep_phi, TopLep_massReg, TopLep_ptReg, TopLep_jetPt, TopLep_csv_nominal;  // FIXME: add systematic vars?
    float TopLep_nuPt, TopLep_nuEta, TopLep_nuPhi;
    float TopHad_mass, TopHad_pt, TopHad_eta, TopHad_phi, TopHad_massReg, TopHad_ptReg, TopHad_j3Pt;  // FIXME: add systematic vars?
    float FatTopHad_mass, FatTopHad_pt, FatTopHad_eta, FatTopHad_phi, FatTopHad_massReg, FatTopHad_ptReg, FatTopHad_j3Pt;

    outTree->Branch("EVENT_run", &EVENT_run, "EVENT_run/I");
    outTree->Branch("EVENT_event", &EVENT_event, "EVENT_event/I");
    outTree->Branch("lumi", &lumi_, "lumi/F");
    outTree->Branch("efflumi", &efflumi, "efflumi/F");
    outTree->Branch("efflumi_old", &efflumi_old, "efflumi_old/F");
    outTree->Branch("efflumi_UEPS_up", &efflumi_UEPS_up, "efflumi_UEPS_up/F");
    outTree->Branch("efflumi_UEPS_down", &efflumi_UEPS_down, "efflumi_UEPS_down/F");
    outTree->Branch("hJet_ptReg", &hJet_ptReg, "hJet_ptReg[2]/F");
    outTree->Branch("hJet_ptReg_res_j_up", &hJet_ptReg_res_j_up, "hJet_ptReg_res_j_up[2]/F");
    outTree->Branch("hJet_ptReg_res_j_down", &hJet_ptReg_res_j_down, "hJet_ptReg_res_j_down[2]/F");
    outTree->Branch("hJet_ptReg_scale_j_up", &hJet_ptReg_scale_j_up, "hJet_ptReg_scale_j_up[2]/F");
    outTree->Branch("hJet_ptReg_scale_j_down", &hJet_ptReg_scale_j_down, "hJet_ptReg_scale_j_down[2]/F");
#ifdef CSVSYST
    outTree->Branch("hJet_csv_reshaped", &hJet_csv_reshaped, "hJet_csv_reshaped[2]/F");
    outTree->Branch("hJet_csv_eff_b_up", &hJet_csv_eff_b_up, "hJet_csv_eff_b_up[2]/F");
    outTree->Branch("hJet_csv_eff_b_down", &hJet_csv_eff_b_down, "hJet_csv_eff_b_down[2]/F");
    outTree->Branch("hJet_csv_fake_b_up", &hJet_csv_fake_b_up, "hJet_csv_fake_b_up[2]/F");
    outTree->Branch("hJet_csv_fake_b_down", &hJet_csv_fake_b_down, "hJet_csv_fake_b_down[2]/F");
#endif
    outTree->Branch("HptNorm", &HptNorm, "HptNorm/F");
    outTree->Branch("HptGen", &HptGen, "HptGen/F");
    outTree->Branch("HptReg", &HptReg, "HptReg/F");
    outTree->Branch("HptReg_res_j_up", &HptReg_res_j_up, "HptReg_res_j_up/F");
    outTree->Branch("HptReg_res_j_down", &HptReg_res_j_down, "HptReg_res_j_down/F");
    outTree->Branch("HptReg_scale_j_up", &HptReg_scale_j_up, "HptReg_scale_j_up/F");
    outTree->Branch("HptReg_scale_j_down", &HptReg_scale_j_down, "HptReg_scale_j_down/F");
    outTree->Branch("HmassNorm", &HmassNorm, "HmassNorm/F");
    outTree->Branch("HmassGen", &HmassGen, "HmassGen/F");
    outTree->Branch("HmassReg", &HmassReg, "HmassReg/F");
    outTree->Branch("HmassReg_res_j_up", &HmassReg_res_j_up, "HmassReg_res_j_up/F");
    outTree->Branch("HmassReg_res_j_down", &HmassReg_res_j_down, "HmassReg_res_j_down/F");
    outTree->Branch("HmassReg_scale_j_up", &HmassReg_scale_j_up, "HmassReg_scale_j_up/F");
    outTree->Branch("HmassReg_scale_j_down", &HmassReg_scale_j_down, "HmassReg_scale_j_down/F");
    outTree->Branch("aJet_pt_res_j_up", &aJet_pt_res_j_up, "aJet_pt_res_j_up[naJets]/F");
    outTree->Branch("aJet_pt_res_j_down", &aJet_pt_res_j_down, "aJet_pt_res_j_down[naJets]/F");
    outTree->Branch("aJet_pt_scale_j_up", &aJet_pt_scale_j_up, "aJet_pt_scale_j_up[naJets]/F");
    outTree->Branch("aJet_pt_scale_j_down", &aJet_pt_scale_j_down, "aJet_pt_scale_j_down[naJets]/F");
    outTree->Branch("fathFilterJets_ptReg", &fathFilterJets_ptReg, "fathFilterJets_ptReg[nfathFilterJets]/F");
    outTree->Branch("fathFilterJets_ptReg_res_j_up", &fathFilterJets_ptReg_res_j_up, "fathFilterJets_ptReg_res_j_up[nfathFilterJets]/F");
    outTree->Branch("fathFilterJets_ptReg_res_j_down", &fathFilterJets_ptReg_res_j_down, "fathFilterJets_ptReg_res_j_down[nfathFilterJets]/F");
    outTree->Branch("fathFilterJets_ptReg_scale_j_up", &fathFilterJets_ptReg_scale_j_up, "fathFilterJets_ptReg_scale_j_up[nfathFilterJets]/F");
    outTree->Branch("fathFilterJets_ptReg_scale_j_down", &fathFilterJets_ptReg_scale_j_down, "fathFilterJets_ptReg_scale_j_down[nfathFilterJets]/F");
#ifdef CSVSYST
    outTree->Branch("fathFilterJets_csv_reshaped", &fathFilterJets_csv_reshaped, "fathFilterJets_csv_reshaped[nfathFilterJets]/F");
#endif
    outTree->Branch("FatHptNorm", &FatHptNorm, "FatHptNorm/F");
    outTree->Branch("FatHptGen", &FatHptGen, "FatHptGen/F");
    outTree->Branch("FatHptReg", &FatHptReg, "FatHptReg/F");
    outTree->Branch("FatHptReg_res_j_up", &FatHptReg_res_j_up, "FatHptReg_res_j_up/F");
    outTree->Branch("FatHptReg_res_j_down", &FatHptReg_res_j_down, "FatHptReg_res_j_down/F");
    outTree->Branch("FatHptReg_scale_j_up", &FatHptReg_scale_j_up, "FatHptReg_scale_j_up/F");
    outTree->Branch("FatHptReg_scale_j_down", &FatHptReg_scale_j_down, "FatHptReg_scale_j_down/F");
    outTree->Branch("FatHmassNorm", &FatHmassNorm, "FatHmassNorm/F");
    outTree->Branch("FatHmassGen", &FatHmassGen, "FatHmassGen/F");
    outTree->Branch("FatHmassReg", &FatHmassReg, "FatHmassReg/F");
    outTree->Branch("FatHmassReg_res_j_up", &FatHmassReg_res_j_up, "FatHmassReg_res_j_up/F");
    outTree->Branch("FatHmassReg_res_j_down", &FatHmassReg_res_j_down, "FatHmassReg_res_j_down/F");
    outTree->Branch("FatHmassReg_scale_j_up", &FatHmassReg_scale_j_up, "FatHmassReg_scale_j_up/F");
    outTree->Branch("FatHmassReg_scale_j_down", &FatHmassReg_scale_j_down, "FatHmassReg_scale_j_down/F");
    outTree->Branch("MET_xyshift", &MET_xyshift, "MET_xyshift/F");
    outTree->Branch("MET_res_j_up", &MET_res_j_up, "MET_res_j_up/F");
    outTree->Branch("MET_res_j_down", &MET_res_j_down, "MET_res_j_down/F");
    outTree->Branch("MET_scale_j_up", &MET_scale_j_up, "MET_scale_j_up/F");
    outTree->Branch("MET_scale_j_down", &MET_scale_j_down, "MET_scale_j_down/F");
    outTree->Branch("METphi_xyshift", &METphi_xyshift, "METphi_xyshift/F");
    outTree->Branch("METphi_res_j_up", &METphi_res_j_up, "METphi_res_j_up/F");
    outTree->Branch("METphi_res_j_down", &METphi_res_j_down, "METphi_res_j_down/F");
    outTree->Branch("METphi_scale_j_up", &METphi_scale_j_up, "METphi_scale_j_up/F");
    outTree->Branch("METphi_scale_j_down", &METphi_scale_j_down, "METphi_scale_j_down/F");
    outTree->Branch("HTeta25", &HTeta25, "HTeta25/F");
    outTree->Branch("HTeta45", &HTeta45, "HTeta45/F");
    
    outTree->Branch("nalep_Znn", &nalep_Znn, "nalep_Znn/I");
    outTree->Branch("nalep_pt5_Znn", &nalep_pt5_Znn, "nalep_pt5_Znn/I");
    outTree->Branch("naJets_Znn", &naJets_Znn, "naJets_Znn/I");
    outTree->Branch("naJets_Znn_res_j_up", &naJets_Znn_res_j_up, "naJets_Znn_res_j_up/I");
    outTree->Branch("naJets_Znn_res_j_down", &naJets_Znn_res_j_down, "naJets_Znn_res_j_down/I");
    outTree->Branch("naJets_Znn_scale_j_up", &naJets_Znn_scale_j_up, "naJets_Znn_scale_j_up/I");
    outTree->Branch("naJets_Znn_scale_j_down", &naJets_Znn_scale_j_down, "naJets_Znn_scale_j_down/I");
    outTree->Branch("naCtrJets_pt20to25_Znn", &naCtrJets_pt20to25_Znn, "naCtrJets_pt20to25_Znn/I");
    outTree->Branch("naCtrJets_pt20to25_Znn_res_j_up", &naCtrJets_pt20to25_Znn_res_j_up, "naCtrJets_pt20to25_Znn_res_j_up/I");
    outTree->Branch("naCtrJets_pt20to25_Znn_res_j_down", &naCtrJets_pt20to25_Znn_res_j_down, "naCtrJets_pt20to25_Znn_res_j_down/I");
    outTree->Branch("naCtrJets_pt20to25_Znn_scale_j_up", &naCtrJets_pt20to25_Znn_scale_j_up, "naCtrJets_pt20to25_Znn_scale_j_up/I");
    outTree->Branch("naCtrJets_pt20to25_Znn_scale_j_down", &naCtrJets_pt20to25_Znn_scale_j_down, "naCtrJets_pt20to25_Znn_scale_j_down/I");
    
    outTree->Branch("mindPhiMETJet_pt", &mindPhiMETJet_pt, "mindPhiMETJet_pt/F");
    outTree->Branch("mindPhiMETJet_dPhi", &mindPhiMETJet_dPhi, "mindPhiMETJet_dPhi/F");
    outTree->Branch("mindPhiMETJet_dPhi_res_j_up", &mindPhiMETJet_dPhi_res_j_up, "mindPhiMETJet_dPhi_res_j_up/F");
    outTree->Branch("mindPhiMETJet_dPhi_res_j_down", &mindPhiMETJet_dPhi_res_j_down, "mindPhiMETJet_dPhi_res_j_down/F");
    outTree->Branch("mindPhiMETJet_dPhi_scale_j_up", &mindPhiMETJet_dPhi_scale_j_up, "mindPhiMETJet_dPhi_scale_j_up/F");
    outTree->Branch("mindPhiMETJet_dPhi_scale_j_down", &mindPhiMETJet_dPhi_scale_j_down, "mindPhiMETJet_dPhi_scale_j_down/F");
    
    outTree->Branch("mindPhiMETCtrJet_pt", &mindPhiMETCtrJet_pt, "mindPhiMETCtrJet_pt/F");
    outTree->Branch("mindPhiMETCtrJet_dPhi", &mindPhiMETCtrJet_dPhi, "mindPhiMETCtrJet_dPhi/F");
    outTree->Branch("mindPhiMETCtrJet_dPhi_res_j_up", &mindPhiMETCtrJet_dPhi_res_j_up, "mindPhiMETCtrJet_dPhi_res_j_up/F");
    outTree->Branch("mindPhiMETCtrJet_dPhi_res_j_down", &mindPhiMETCtrJet_dPhi_res_j_down, "mindPhiMETCtrJet_dPhi_res_j_down/F");
    outTree->Branch("mindPhiMETCtrJet_dPhi_scale_j_up", &mindPhiMETCtrJet_dPhi_scale_j_up, "mindPhiMETCtrJet_dPhi_scale_j_up/F");
    outTree->Branch("mindPhiMETCtrJet_dPhi_scale_j_down", &mindPhiMETCtrJet_dPhi_scale_j_down, "mindPhiMETCtrJet_dPhi_scale_j_down/F");
    
    outTree->Branch("mindRHJet_pt", &mindRHJet_pt, "mindRHJet_pt/F");
    outTree->Branch("mindRHJet_pt_res_j_up", &mindRHJet_pt_res_j_up, "mindRHJet_pt_res_j_up/F");
    outTree->Branch("mindRHJet_pt_res_j_down", &mindRHJet_pt_res_j_down, "mindRHJet_pt_res_j_down/F");
    outTree->Branch("mindRHJet_pt_scale_j_up", &mindRHJet_pt_scale_j_up, "mindRHJet_pt_scale_j_up/F");
    outTree->Branch("mindRHJet_pt_scale_j_down", &mindRHJet_pt_scale_j_down, "mindRHJet_pt_scale_j_down/F");
    outTree->Branch("mindRHJet_dR", &mindRHJet_dR, "mindRHJet_dR/F");
    outTree->Branch("mindRHJet_dR_res_j_up", &mindRHJet_dR_res_j_up, "mindRHJet_dR_res_j_up/F");
    outTree->Branch("mindRHJet_dR_res_j_down", &mindRHJet_dR_res_j_down, "mindRHJet_dR_res_j_down/F");
    outTree->Branch("mindRHJet_dR_scale_j_up", &mindRHJet_dR_scale_j_up, "mindRHJet_dR_scale_j_up/F");
    outTree->Branch("mindRHJet_dR_scale_j_down", &mindRHJet_dR_scale_j_down, "mindRHJet_dR_scale_j_down/F");
    outTree->Branch("mindRHJet_csv_nominal", &mindRHJet_csv_nominal, "mindRHJet_nominal/F");
    
    outTree->Branch("weightTrigMu", &weightTrigMu, "weightTrigMu/F");
    outTree->Branch("weightTrigDLP_ElecMay10", &weightTrigDLP_ElecMay10, "weightTrigDLP_ElecMay10/F");
    outTree->Branch("weightTrigDLP_Elec", &weightTrigDLP_Elec, "weightTrigDLP_Elec/F");
    
    outTree->Branch("TopLep_mass", &TopLep_mass, "TopLep_mass/F");
    outTree->Branch("TopLep_pt", &TopLep_pt, "TopLep_pt/F");
    outTree->Branch("TopLep_eta", &TopLep_eta, "TopLep_eta/F");
    outTree->Branch("TopLep_phi", &TopLep_phi, "TopLep_phi/F");
    outTree->Branch("TopLep_massReg", &TopLep_massReg, "TopLep_massReg/F");
    outTree->Branch("TopLep_ptReg", &TopLep_ptReg, "TopLep_ptReg/F");
    outTree->Branch("TopLep_jetPt", &TopLep_jetPt, "TopLep_jetPt/F");
    outTree->Branch("TopLep_csv_nominal", &TopLep_csv_nominal, "TopLep_csv_nominal/F");
    outTree->Branch("TopLep_nuPt", &TopLep_nuPt, "TopLep_nuPt/F");
    outTree->Branch("TopLep_nuEta", &TopLep_nuEta, "TopLep_nuEta/F");
    outTree->Branch("TopLep_nuPhi", &TopLep_nuPhi, "TopLep_nuPhi/F");
    outTree->Branch("TopHad_mass", &TopHad_mass, "TopHad_mass/F");
    outTree->Branch("TopHad_pt", &TopHad_pt, "TopHad_pt/F");
    outTree->Branch("TopHad_eta", &TopHad_eta, "TopHad_eta/F");
    outTree->Branch("TopHad_phi", &TopHad_phi, "TopHad_phi/F");
    outTree->Branch("TopHad_massReg", &TopHad_massReg, "TopHad_massReg/F");
    outTree->Branch("TopHad_ptReg", &TopHad_ptReg, "TopHad_ptReg/F");
    outTree->Branch("TopHad_j3Pt", &TopHad_j3Pt, "TopHad_j3Pt/F");
    outTree->Branch("FatTopHad_mass", &FatTopHad_mass, "FatTopHad_mass/F");
    outTree->Branch("FatTopHad_pt", &FatTopHad_pt, "FatTopHad_pt/F");
    outTree->Branch("FatTopHad_eta", &FatTopHad_eta, "FatTopHad_eta/F");
    outTree->Branch("FatTopHad_phi", &FatTopHad_phi, "FatTopHad_phi/F");
    outTree->Branch("FatTopHad_massReg", &FatTopHad_massReg, "FatTopHad_massReg/F");
    outTree->Branch("FatTopHad_ptReg", &FatTopHad_ptReg, "FatTopHad_ptReg/F");
    outTree->Branch("FatTopHad_j3Pt", &FatTopHad_j3Pt, "FatTopHad_j3Pt/F");

    /// Get effective lumis
    std::map < std::string, float > efflumis = GetLumis();
    efflumi = efflumis[process.Data()];
    assert(efflumi > 0);
    efflumi_old       = efflumi;
    efflumi_UEPS_up   = efflumi * hcount->GetBinContent(2) / hcount->GetBinContent(3);
    efflumi_UEPS_down = efflumi * hcount->GetBinContent(2) / hcount->GetBinContent(4);

    TTreeFormula* ttf_lheweight = new TTreeFormula("ttf_lheweight", Form("%f", efflumi), inTree);
#ifdef STITCH
    std::map < std::string, std::string > lheweights = GetLHEWeights();
    TString process_lhe = process;
    if (process_lhe.BeginsWith("WJets") && process_lhe != "WJetsHW")
        process_lhe = "WJets";
    else if (process_lhe.BeginsWith("ZJets") && process_lhe != "ZJetsHW")
        process_lhe = "ZJets";
    else 
        process_lhe = "";
    TString lheweight = lheweights[process_lhe.Data()];
    if (lheweight != "") {
        delete ttf_lheweight;
        
        // Bug fix for ZJetsPtZ100
        if (process == "ZJetsPtZ100")
            lheweight.ReplaceAll("lheV_pt", "999");
        std::cout << "BUGFIX: " << lheweight << std::endl;
        ttf_lheweight = new TTreeFormula("ttf_lheweight", lheweight, inTree);
    }
#endif
    ttf_lheweight->SetQuickLoad(1);

    ///-- Setup TMVA Reader ----------------------------------------------------
    TMVA::Tools::Instance();  //< This loads the library
    TMVA::Reader * reader = new TMVA::Reader("!Color:!Silent");

    /// Get the variables
    const std::vector < std::string > & inputExpressionsReg = GetInputExpressionsReg();
    const UInt_t nvars = inputExpressionsReg.size();
    Float_t readerVars[nvars];
    int idx_rawpt = -1, idx_pt = -1, idx_et = -1, idx_mt = -1;
    for (UInt_t iexpr = 0; iexpr < nvars; iexpr++) {
        const TString& expr = inputExpressionsReg.at(iexpr);
        reader->AddVariable(expr, &readerVars[iexpr]);
        if      (expr.BeginsWith("breg_rawptJER := "))  idx_rawpt = iexpr;
        else if (expr.BeginsWith("breg_pt := "))        idx_pt = iexpr;
        else if (expr.BeginsWith("breg_et := "))        idx_et = iexpr;
        else if (expr.BeginsWith("breg_mt := "))        idx_mt = iexpr;
    }
    assert(idx_rawpt!=-1 && idx_pt!=-1 && idx_et!=-1 && idx_mt!=-1);

    /// Setup TMVA regression inputs
    const std::vector < std::string > & inputExpressionsReg0 = GetInputExpressionsReg0();
    const std::vector < std::string > & inputExpressionsReg1 = GetInputExpressionsReg1();
    assert(inputExpressionsReg0.size() == nvars);
    assert(inputExpressionsReg1.size() == nvars);

    /// Load TMVA weights
    TString weightdir  = "weights/";
    TString weightfile = weightdir + "TMVARegression_" + regMethod + ".testweights.xml";
    reader->BookMVA(regMethod + " method", weightfile);

    ///-- Setup TMVA Reader (filter jets) ---------------------------------------
    TMVA::Reader * readerFJ = new TMVA::Reader("!Color:!Silent");

    /// Get the variables (filter jets)
    const std::vector < std::string > & inputExpressionsFJReg = GetInputExpressionsFJReg();
    const UInt_t nvarsFJ = inputExpressionsFJReg.size();
    Float_t readerVarsFJ[nvarsFJ];
    int idx_fjrawpt = -1, idx_fjpt = -1, idx_fjet = -1, idx_fjmt = -1;
    for (UInt_t iexpr = 0; iexpr < nvarsFJ; iexpr++) {
        const TString& expr = inputExpressionsFJReg.at(iexpr);
        readerFJ->AddVariable(expr, &readerVarsFJ[iexpr]);
        //if      (expr.BeginsWith("breg_fjrawpt := "))  idx_fjrawpt = iexpr;
        if      (expr.BeginsWith("breg_fjrawptJER := "))  idx_fjrawpt = iexpr;
        else if (expr.BeginsWith("breg_fjpt := "))     idx_fjpt = iexpr;
        else if (expr.BeginsWith("breg_fjet := "))     idx_fjet = iexpr;
        else if (expr.BeginsWith("breg_fjmt := "))     idx_fjmt = iexpr;
    }
    assert(idx_fjrawpt!=-1 && idx_fjpt!=-1 && idx_fjet!=-1 && idx_fjmt!=-1);

    /// Setup TMVA regression inputs (filter jets)
    const std::vector < std::string > & inputExpressionsFJReg0 = GetInputExpressionsFJReg0();
    const std::vector < std::string > & inputExpressionsFJReg1 = GetInputExpressionsFJReg1();
    const std::vector < std::string > & inputExpressionsFJReg2 = GetInputExpressionsFJReg2();
    assert(inputExpressionsFJReg0.size() == nvarsFJ);
    assert(inputExpressionsFJReg1.size() == nvarsFJ);
    assert(inputExpressionsFJReg2.size() == nvarsFJ);

    /// Load TMVA weights (filter jets)
    weightdir  = "weights/";
    weightfile = weightdir + "TMVARegressionFJ_" + regMethod + ".testweights.xml";
    readerFJ->BookMVA(regMethod + " method", weightfile);


    if (beginEntry == 0 && endEntry == -1)
        std::cout << "--- Processing: " << inTree->GetEntriesFast() << " events" << std::endl;
    else
        std::cout << "--- Processing: " << beginEntry << " - " << endEntry << " from " << inTree->GetEntriesFast() << " events" << std::endl;
    TStopwatch sw;
    sw.Start();

    /// Create TTreeFormulas
    TTreeFormula *ttf = 0;
    std::vector < TTreeFormula * >::const_iterator formIt, formItEnd;
    std::vector < TTreeFormula * > inputFormulasReg0;
    std::vector < TTreeFormula * > inputFormulasReg1;
    std::vector < TTreeFormula * > inputFormulasFJReg0;
    std::vector < TTreeFormula * > inputFormulasFJReg1;
    std::vector < TTreeFormula * > inputFormulasFJReg2;
    for (UInt_t iexpr = 0; iexpr < nvars; iexpr++) {
        ttf = new TTreeFormula(Form("ttfreg%i_0", iexpr), inputExpressionsReg0.at(iexpr).c_str(), inTree);
        ttf->SetQuickLoad(1);
        inputFormulasReg0.push_back(ttf);
        ttf = new TTreeFormula(Form("ttfreg%i_1", iexpr), inputExpressionsReg1.at(iexpr).c_str(), inTree);
        ttf->SetQuickLoad(1);
        inputFormulasReg1.push_back(ttf);
    }
    for (UInt_t iexpr = 0; iexpr < nvarsFJ; iexpr++) {
        ttf = new TTreeFormula(Form("ttffjreg%i_0", iexpr), inputExpressionsFJReg0.at(iexpr).c_str(), inTree);
        ttf->SetQuickLoad(1);
        inputFormulasFJReg0.push_back(ttf);
        ttf = new TTreeFormula(Form("ttffjreg%i_1", iexpr), inputExpressionsFJReg1.at(iexpr).c_str(), inTree);
        ttf->SetQuickLoad(1);
        inputFormulasFJReg1.push_back(ttf);
        ttf = new TTreeFormula(Form("ttffjreg%i_2", iexpr), inputExpressionsFJReg2.at(iexpr).c_str(), inTree);
        ttf->SetQuickLoad(1);
        inputFormulasFJReg2.push_back(ttf);
    }


#ifdef CSVSYST
    /// Setup b-tagging reshaping
    TString csvFile = "csvdiscr.root";
    BTagShapeInterface *csvShape = new BTagShapeInterface(csvFile, 0.0, 0.0, true);     // use 4 points!
    BTagShapeInterface *csvShape_eff_b_up = new BTagShapeInterface(csvFile, 1.5, 0.0, true);    // use 4 points!
    BTagShapeInterface *csvShape_eff_b_down = new BTagShapeInterface(csvFile, -1.5, 0.0, true); // use 4 points!
    BTagShapeInterface *csvShape_fake_b_up = new BTagShapeInterface(csvFile, 0.0, 1.0, true);   // use 4 points!
    BTagShapeInterface *csvShape_fake_b_down = new BTagShapeInterface(csvFile, 0.0, -1.0, true);        // use 4 points!
#endif

#ifdef JECFWLITE
    JECFWLiteStandalone jec("jec", "AK5PFchs");
#endif

    ///-- Loop over events -----------------------------------------------------
    Int_t curTree = inTree->GetTreeNumber();
    const Long64_t nentries = inTree->GetEntries();
    if (endEntry < 0)  endEntry = nentries;

    Long64_t ievt = 0;
    for (ievt=TMath::Max(ievt, beginEntry); ievt<TMath::Min(nentries, endEntry); ievt++) {
        if (ievt % 2000 == 0)
            std::cout << "--- ... Processing event: " << ievt << std::endl;

        const Long64_t local_entry = inTree->LoadTree(ievt);  // faster, but only for TTreeFormula
        if (local_entry < 0)  break;
        inTree->GetEntry(ievt);  // same event as received by LoadTree()

        if (inTree->GetTreeNumber() != curTree) {
            curTree = inTree->GetTreeNumber();
            for (formIt=inputFormulasReg0.begin(), formItEnd=inputFormulasReg0.end(); formIt!=formItEnd; formIt++)
                (*formIt)->UpdateFormulaLeaves();  // if using TChain
            for (formIt=inputFormulasReg1.begin(), formItEnd=inputFormulasReg1.end(); formIt!=formItEnd; formIt++)
                (*formIt)->UpdateFormulaLeaves();  // if using TChain
            for (formIt=inputFormulasFJReg0.begin(), formItEnd=inputFormulasFJReg0.end(); formIt!=formItEnd; formIt++)
                (*formIt)->UpdateFormulaLeaves();  // if using TChain
            for (formIt=inputFormulasFJReg1.begin(), formItEnd=inputFormulasFJReg1.end(); formIt!=formItEnd; formIt++)
                (*formIt)->UpdateFormulaLeaves();  // if using TChain
            for (formIt=inputFormulasFJReg2.begin(), formItEnd=inputFormulasFJReg2.end(); formIt!=formItEnd; formIt++)
                (*formIt)->UpdateFormulaLeaves();  // if using TChain
            ttf_lheweight->UpdateFormulaLeaves();
        }

        /// These need to be called when arrays of variable size are used in TTree.
        for (formIt=inputFormulasReg0.begin(), formItEnd=inputFormulasReg0.end(); formIt!=formItEnd; formIt++)
            (*formIt)->GetNdata();
        for (formIt=inputFormulasReg1.begin(), formItEnd=inputFormulasReg1.end(); formIt!=formItEnd; formIt++)
            (*formIt)->GetNdata();
        for (formIt=inputFormulasFJReg0.begin(), formItEnd=inputFormulasFJReg0.end(); formIt!=formItEnd; formIt++)
            (*formIt)->GetNdata();
        for (formIt=inputFormulasFJReg1.begin(), formItEnd=inputFormulasFJReg1.end(); formIt!=formItEnd; formIt++)
            (*formIt)->GetNdata();
        for (formIt=inputFormulasFJReg2.begin(), formItEnd=inputFormulasFJReg2.end(); formIt!=formItEnd; formIt++)
            (*formIt)->GetNdata();
        ttf_lheweight->GetNdata();

        /// Fill branches
        EVENT_run = EVENT.run;
        EVENT_event = EVENT.event;

#ifdef STITCH        
        efflumi           = ttf_lheweight->EvalInstance();
        efflumi_UEPS_up   = efflumi * hcount->GetBinContent(2) / hcount->GetBinContent(3);
        efflumi_UEPS_down = efflumi * hcount->GetBinContent(2) / hcount->GetBinContent(4);
#endif

        bool verbose = false;
        double rawpt_ = 0., pt_ = 0., e_ = 0.;
        for (Int_t ihj = 0; ihj < 2; ihj++) {

#ifdef JECFWLite
            hJet_JECUnc[ihj] = uncert(hJet_eta[ihj], hJet_pt[ihj], (EVENT_run == 1), false);
#endif

            /// Evaluate TMVA regression output
            for (UInt_t iexpr = 0; iexpr < nvars; iexpr++) {
                if (ihj==0) {
                    readerVars[iexpr] = inputFormulasReg0.at(iexpr)->EvalInstance();
                } else if (ihj==1) {
                    readerVars[iexpr] = inputFormulasReg1.at(iexpr)->EvalInstance();
                }
            }
            assert(TMath::AreEqualRel(hJet_pt[ihj], readerVars[idx_pt], 1.E-12));
            const double hJet_ptRawJER    = readerVars[idx_rawpt];  // hJet_ptRaw[ihj] is before smearing
            hJet_ptReg[ihj]               = (reader->EvaluateRegression(regMethod + " method"))[0];
            if (verbose)  std::cout << readerVars[idx_pt] << " " << readerVars[idx_rawpt] << " " << readerVars[idx_et] << " " << readerVars[idx_mt] << " " << hJet_pt[ihj] << " " << hJet_ptReg[ihj] << " " << hJet_genPt[ihj] << std::endl;

            // res_j_up
            pt_                           = smear_pt_resErr(hJet_pt[ihj] , hJet_genPt[ihj], hJet_eta[ihj], true);
            rawpt_                        = pt_ / hJet_pt[ihj] * hJet_ptRawJER;
            e_                            = pt_ / hJet_pt[ihj] * hJet_e[ihj];
            readerVars[idx_pt]            = pt_;
            readerVars[idx_rawpt]         = rawpt_;
            readerVars[idx_et]            = evalEt(pt_, hJet_eta[ihj], hJet_phi[ihj], e_);
            readerVars[idx_mt]            = evalMt(pt_, hJet_eta[ihj], hJet_phi[ihj], e_);
            hJet_ptReg_res_j_up[ihj]      = (reader->EvaluateRegression(regMethod + " method"))[0];
            if (verbose)  std::cout << readerVars[idx_pt] << " " << readerVars[idx_rawpt] << " " << readerVars[idx_et] << " " << readerVars[idx_mt] << " " << hJet_pt[ihj] << " " << hJet_ptReg[ihj] << " " << hJet_ptReg_res_j_up[ihj] << std::endl;

            // res_j_down
            pt_                           = smear_pt_resErr(hJet_pt[ihj] , hJet_genPt[ihj], hJet_eta[ihj], false);
            rawpt_                        = pt_ / hJet_pt[ihj] * hJet_ptRawJER;
            e_                            = pt_ / hJet_pt[ihj] * hJet_e[ihj];
            readerVars[idx_pt]            = pt_;
            readerVars[idx_rawpt]         = rawpt_;
            readerVars[idx_et]            = evalEt(pt_, hJet_eta[ihj], hJet_phi[ihj], e_);
            readerVars[idx_mt]            = evalMt(pt_, hJet_eta[ihj], hJet_phi[ihj], e_);
            hJet_ptReg_res_j_down[ihj]    = (reader->EvaluateRegression(regMethod + " method"))[0];
            if (verbose)  std::cout << readerVars[idx_pt] << " " << readerVars[idx_rawpt] << " " << readerVars[idx_et] << " " << readerVars[idx_mt] << " " << hJet_pt[ihj] << " " << hJet_ptReg[ihj] << " " << hJet_ptReg_res_j_down[ihj] << std::endl;

            // scale_j_up
            pt_                           = smear_pt_scaleErr(hJet_pt[ihj] , hJet_JECUnc[ihj], true);
            rawpt_                        = pt_ / hJet_pt[ihj] * hJet_ptRawJER;
            e_                            = pt_ / hJet_pt[ihj] * hJet_e[ihj];
            readerVars[idx_pt]            = pt_;
            readerVars[idx_rawpt]         = rawpt_;
            readerVars[idx_et]            = evalEt(pt_, hJet_eta[ihj], hJet_phi[ihj], e_);
            readerVars[idx_mt]            = evalMt(pt_, hJet_eta[ihj], hJet_phi[ihj], e_);
            hJet_ptReg_scale_j_up[ihj]    = (reader->EvaluateRegression(regMethod + " method"))[0];
            if (verbose)  std::cout << readerVars[idx_pt] << " " << readerVars[idx_rawpt] << " " << readerVars[idx_et] << " " << readerVars[idx_mt] << " " << hJet_pt[ihj] << " " << hJet_ptReg[ihj] << " " << hJet_ptReg_scale_j_up[ihj] << std::endl;

            // scale_j_down
            pt_                           = smear_pt_scaleErr(hJet_pt[ihj] , hJet_JECUnc[ihj], false);
            rawpt_                        = pt_ / hJet_pt[ihj] * hJet_ptRawJER;
            e_                            = pt_ / hJet_pt[ihj] * hJet_e[ihj];
            readerVars[idx_pt]            = pt_;
            readerVars[idx_rawpt]         = rawpt_;
            readerVars[idx_et]            = evalEt(pt_, hJet_eta[ihj], hJet_phi[ihj], e_);
            readerVars[idx_mt]            = evalMt(pt_, hJet_eta[ihj], hJet_phi[ihj], e_);
            hJet_ptReg_scale_j_down[ihj]  = (reader->EvaluateRegression(regMethod + " method"))[0];
            if (verbose)  std::cout << readerVars[idx_pt] << " " << readerVars[idx_rawpt] << " " << readerVars[idx_et] << " " << readerVars[idx_mt] << " " << hJet_pt[ihj] << " " << hJet_ptReg[ihj] << " " << hJet_ptReg_scale_j_down[ihj] << std::endl;

#ifdef CSVSYST
            hJet_csv_reshaped[ihj]    = csvShape            ->reshape(hJet_eta[ihj], hJet_ptReg[ihj], hJet_csv[ihj], hJet_flavour[ihj]);
            hJet_csv_eff_b_up[ihj]    = csvShape_eff_b_up   ->reshape(hJet_eta[ihj], hJet_ptReg[ihj], hJet_csv[ihj], hJet_flavour[ihj]);
            hJet_csv_eff_b_down[ihj]  = csvShape_eff_b_down ->reshape(hJet_eta[ihj], hJet_ptReg[ihj], hJet_csv[ihj], hJet_flavour[ihj]);
            hJet_csv_fake_b_up[ihj]   = csvShape_fake_b_up  ->reshape(hJet_eta[ihj], hJet_ptReg[ihj], hJet_csv[ihj], hJet_flavour[ihj]);
            hJet_csv_fake_b_down[ihj] = csvShape_fake_b_down->reshape(hJet_eta[ihj], hJet_ptReg[ihj], hJet_csv[ihj], hJet_flavour[ihj]);
#endif
        }

        const TLorentzVector p4Zero                     = TLorentzVector(0., 0., 0., 0.);

        const TLorentzVector& hJet_p4Norm_0             = makePtEtaPhiE(hJet_pt[0]                , hJet_pt[0], hJet_eta[0], hJet_phi[0], hJet_e[0]);
        const TLorentzVector& hJet_p4Norm_1             = makePtEtaPhiE(hJet_pt[1]                , hJet_pt[1], hJet_eta[1], hJet_phi[1], hJet_e[1]);
        const TLorentzVector& hJet_p4Gen_0              = hJet_genPt[0] > 0 ? 
                                                          makePtEtaPhiE(hJet_genPt[0]             , hJet_pt[0], hJet_eta[0], hJet_phi[0], hJet_e[0]) : p4Zero;
        const TLorentzVector& hJet_p4Gen_1              = hJet_genPt[1] > 0 ? 
                                                          makePtEtaPhiE(hJet_genPt[1]             , hJet_pt[1], hJet_eta[1], hJet_phi[1], hJet_e[1]) : p4Zero;
        const TLorentzVector& hJet_p4Reg_0              = makePtEtaPhiE(hJet_ptReg[0]             , hJet_pt[0], hJet_eta[0], hJet_phi[0], hJet_e[0]);
        const TLorentzVector& hJet_p4Reg_1              = makePtEtaPhiE(hJet_ptReg[1]             , hJet_pt[1], hJet_eta[1], hJet_phi[1], hJet_e[1]);
        const TLorentzVector& hJet_p4Reg_res_j_up_0     = makePtEtaPhiE(hJet_ptReg_res_j_up[0]    , hJet_pt[0], hJet_eta[0], hJet_phi[0], hJet_e[0]);
        const TLorentzVector& hJet_p4Reg_res_j_up_1     = makePtEtaPhiE(hJet_ptReg_res_j_up[1]    , hJet_pt[1], hJet_eta[1], hJet_phi[1], hJet_e[1]);
        const TLorentzVector& hJet_p4Reg_res_j_down_0   = makePtEtaPhiE(hJet_ptReg_res_j_down[0]  , hJet_pt[0], hJet_eta[0], hJet_phi[0], hJet_e[0]);
        const TLorentzVector& hJet_p4Reg_res_j_down_1   = makePtEtaPhiE(hJet_ptReg_res_j_down[1]  , hJet_pt[1], hJet_eta[1], hJet_phi[1], hJet_e[1]);
        const TLorentzVector& hJet_p4Reg_scale_j_up_0   = makePtEtaPhiE(hJet_ptReg_scale_j_up[0]  , hJet_pt[0], hJet_eta[0], hJet_phi[0], hJet_e[0]);
        const TLorentzVector& hJet_p4Reg_scale_j_up_1   = makePtEtaPhiE(hJet_ptReg_scale_j_up[1]  , hJet_pt[1], hJet_eta[1], hJet_phi[1], hJet_e[1]);
        const TLorentzVector& hJet_p4Reg_scale_j_down_0 = makePtEtaPhiE(hJet_ptReg_scale_j_down[0], hJet_pt[0], hJet_eta[0], hJet_phi[0], hJet_e[0]);
        const TLorentzVector& hJet_p4Reg_scale_j_down_1 = makePtEtaPhiE(hJet_ptReg_scale_j_down[1], hJet_pt[1], hJet_eta[1], hJet_phi[1], hJet_e[1]);

        HptNorm             = (hJet_p4Norm_0             + hJet_p4Norm_1            ).Pt();
        HptGen              = (hJet_p4Gen_0              + hJet_p4Gen_1             ).Pt();
        HptReg              = (hJet_p4Reg_0              + hJet_p4Reg_1             ).Pt();
        HptReg_res_j_up     = (hJet_p4Reg_res_j_up_0     + hJet_p4Reg_res_j_up_1    ).Pt();
        HptReg_res_j_down   = (hJet_p4Reg_res_j_down_0   + hJet_p4Reg_res_j_down_1  ).Pt();
        HptReg_scale_j_up   = (hJet_p4Reg_scale_j_up_0   + hJet_p4Reg_scale_j_up_1  ).Pt();
        HptReg_scale_j_down = (hJet_p4Reg_scale_j_down_0 + hJet_p4Reg_scale_j_down_1).Pt();

        HmassNorm             = (hJet_p4Norm_0             + hJet_p4Norm_1            ).M();
        HmassGen              = (hJet_p4Gen_0              + hJet_p4Gen_1             ).M();
        HmassReg              = (hJet_p4Reg_0              + hJet_p4Reg_1             ).M();
        HmassReg_res_j_up     = (hJet_p4Reg_res_j_up_0     + hJet_p4Reg_res_j_up_1    ).M();
        HmassReg_res_j_down   = (hJet_p4Reg_res_j_down_0   + hJet_p4Reg_res_j_down_1  ).M();
        HmassReg_scale_j_up   = (hJet_p4Reg_scale_j_up_0   + hJet_p4Reg_scale_j_up_1  ).M();
        HmassReg_scale_j_down = (hJet_p4Reg_scale_j_down_0 + hJet_p4Reg_scale_j_down_1).M();

        for (Int_t iaj = 0; iaj < naJets; iaj++) {

#ifdef JECFWLite
            aJet_JECUnc[iaj] = uncert(aJet_eta[iaj], aJet_pt[iaj], (EVENT_run == 1), false);
#endif

            aJet_pt_res_j_up[iaj]     = smear_pt_resErr(aJet_pt[iaj], aJet_genPt[iaj], aJet_eta[iaj], true);
            aJet_pt_res_j_down[iaj]   = smear_pt_resErr(aJet_pt[iaj], aJet_genPt[iaj], aJet_eta[iaj], false);
            aJet_pt_scale_j_up[iaj]   = smear_pt_scaleErr(aJet_pt[iaj], aJet_JECUnc[iaj], true);
            aJet_pt_scale_j_down[iaj] = smear_pt_scaleErr(aJet_pt[iaj], aJet_JECUnc[iaj], false);
        }


        for (Int_t ihj = 0; ihj < nfathFilterJets; ihj++) {

#ifdef JECFWLite
            fathFilterJets_JECUnc[ihj] = uncert(fathFilterJets_eta[ihj], fathFilterJets_pt[ihj], (EVENT_run == 1), false);
#endif

            /// Evaluate TMVA regression output (filter jets)
            /// Only for filter jets pT > 15, |eta| < 2.5
            if (fathFilterJets_pt[ihj] > 15. && fabs(fathFilterJets_eta[ihj]) < 2.5) {
                for (UInt_t iexpr = 0; iexpr < nvarsFJ; iexpr++) {
                    if (ihj==0) {
                        readerVarsFJ[iexpr] = inputFormulasFJReg0.at(iexpr)->EvalInstance();
                    } else if (ihj==1) {
                        readerVarsFJ[iexpr] = inputFormulasFJReg1.at(iexpr)->EvalInstance();
                    } else if (ihj==2) {
                        readerVarsFJ[iexpr] = inputFormulasFJReg2.at(iexpr)->EvalInstance();
                    }
                }
                assert(TMath::AreEqualRel(fathFilterJets_pt[ihj], readerVarsFJ[idx_pt], 1.E-12));
                const double fathFilterJets_ptRawJER    = readerVarsFJ[idx_fjrawpt];
                fathFilterJets_ptReg[ihj]               = (readerFJ->EvaluateRegression(regMethod + " method"))[0];
                // res_j_up
                pt_                                     = smear_pt_resErr(fathFilterJets_pt[ihj] , fathFilterJets_genPt[ihj], fathFilterJets_eta[ihj], true);
                rawpt_                                  = pt_ / fathFilterJets_pt[ihj] * fathFilterJets_ptRawJER;
                e_                                      = pt_ / fathFilterJets_pt[ihj] * fathFilterJets_e[ihj];
                readerVarsFJ[idx_pt]                    = pt_;
                readerVarsFJ[idx_rawpt]                 = rawpt_;
                readerVarsFJ[idx_et]                    = evalEt(pt_, fathFilterJets_eta[ihj], fathFilterJets_phi[ihj], e_);
                readerVarsFJ[idx_mt]                    = evalMt(pt_, fathFilterJets_eta[ihj], fathFilterJets_phi[ihj], e_);
                fathFilterJets_ptReg_res_j_up[ihj]      = (readerFJ->EvaluateRegression(regMethod + " method"))[0];  // FIXME: don't propagate through regression?

                // res_j_down
                pt_                                     = smear_pt_resErr(fathFilterJets_pt[ihj] , fathFilterJets_genPt[ihj], fathFilterJets_eta[ihj], false);
                rawpt_                                  = pt_ / fathFilterJets_pt[ihj] * fathFilterJets_ptRawJER;
                e_                                      = pt_ / fathFilterJets_pt[ihj] * fathFilterJets_e[ihj];
                readerVarsFJ[idx_pt]                    = pt_;
                readerVarsFJ[idx_rawpt]                 = rawpt_;
                readerVarsFJ[idx_et]                    = evalEt(pt_, fathFilterJets_eta[ihj], fathFilterJets_phi[ihj], e_);
                readerVarsFJ[idx_mt]                    = evalMt(pt_, fathFilterJets_eta[ihj], fathFilterJets_phi[ihj], e_);
                fathFilterJets_ptReg_res_j_down[ihj]    = (readerFJ->EvaluateRegression(regMethod + " method"))[0];  // FIXME: don't propagate through regression?

                // scale_j_up
                pt_                                     = smear_pt_scaleErr(fathFilterJets_pt[ihj] , fathFilterJets_JECUnc[ihj], true);
                rawpt_                                  = pt_ / fathFilterJets_pt[ihj] * fathFilterJets_ptRawJER;
                e_                                      = pt_ / fathFilterJets_pt[ihj] * fathFilterJets_e[ihj];
                readerVarsFJ[idx_pt]                    = pt_;
                readerVarsFJ[idx_rawpt]                 = rawpt_;
                readerVarsFJ[idx_et]                    = evalEt(pt_, fathFilterJets_eta[ihj], fathFilterJets_phi[ihj], e_);
                readerVarsFJ[idx_mt]                    = evalMt(pt_, fathFilterJets_eta[ihj], fathFilterJets_phi[ihj], e_);
                fathFilterJets_ptReg_scale_j_up[ihj]    = (readerFJ->EvaluateRegression(regMethod + " method"))[0];  // FIXME: don't propagate through regression?

                // scale_j_down
                pt_                                     = smear_pt_scaleErr(fathFilterJets_pt[ihj] , fathFilterJets_JECUnc[ihj], false);
                rawpt_                                  = pt_ / fathFilterJets_pt[ihj] * fathFilterJets_ptRawJER;
                e_                                      = pt_ / fathFilterJets_pt[ihj] * fathFilterJets_e[ihj];
                readerVarsFJ[idx_pt]                    = pt_;
                readerVarsFJ[idx_rawpt]                 = rawpt_;
                readerVarsFJ[idx_et]                    = evalEt(pt_, fathFilterJets_eta[ihj], fathFilterJets_phi[ihj], e_);
                readerVarsFJ[idx_mt]                    = evalMt(pt_, fathFilterJets_eta[ihj], fathFilterJets_phi[ihj], e_);
                fathFilterJets_ptReg_scale_j_down[ihj]  = (readerFJ->EvaluateRegression(regMethod + " method"))[0];  // FIXME: don't propagate through regression?
            } else if (fathFilterJets_pt[ihj] > 0.) {
                fathFilterJets_ptReg[ihj]               = fathFilterJets_pt[ihj];
                fathFilterJets_ptReg_res_j_up[ihj]      = fathFilterJets_pt[ihj];
                fathFilterJets_ptReg_res_j_down[ihj]    = fathFilterJets_pt[ihj];
                fathFilterJets_ptReg_scale_j_up[ihj]    = fathFilterJets_pt[ihj];
                fathFilterJets_ptReg_scale_j_down[ihj]  = fathFilterJets_pt[ihj];
            } else {
                fathFilterJets_ptReg[ihj]               = -99.;
                fathFilterJets_ptReg_res_j_up[ihj]      = -99.;
                fathFilterJets_ptReg_res_j_down[ihj]    = -99.;
                fathFilterJets_ptReg_scale_j_up[ihj]    = -99.;
                fathFilterJets_ptReg_scale_j_down[ihj]  = -99.;
            }
#ifdef CSVSYST
            fathFilterJets_csv_reshaped[ihj] = csvShape->reshape(fathFilterJets_eta[ihj], fathFilterJets_ptReg[ihj], fathFilterJets_csv[ihj], fathFilterJets_flavour[ihj]);
#endif
        }

        assert(nfathFilterJets != 1); // either 0, 2 or 3
        if (nfathFilterJets>0 && fathFilterJets_pt[0]>0. && fathFilterJets_pt[1]>0.){
            const TLorentzVector& fathFilterJets_p4Norm_0             = makePtEtaPhiE(fathFilterJets_pt[0]                , fathFilterJets_pt[0], fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_e[0]);
            const TLorentzVector& fathFilterJets_p4Norm_1             = makePtEtaPhiE(fathFilterJets_pt[1]                , fathFilterJets_pt[1], fathFilterJets_eta[1], fathFilterJets_phi[1], fathFilterJets_e[1]);
            const TLorentzVector& fathFilterJets_p4Norm_2             = (nfathFilterJets > 2 && fathFilterJets_pt[2] > 0) ? 
                                                                        makePtEtaPhiE(fathFilterJets_pt[2]                , fathFilterJets_pt[2], fathFilterJets_eta[2], fathFilterJets_phi[2], fathFilterJets_e[2]) : p4Zero;
            const TLorentzVector& fathFilterJets_p4Gen_0              = fathFilterJets_genPt[0] > 0 ? 
                                                                        makePtEtaPhiE(fathFilterJets_genPt[0]             , fathFilterJets_pt[0], fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_e[0]) : p4Zero;
            const TLorentzVector& fathFilterJets_p4Gen_1              = fathFilterJets_genPt[1] > 0 ? 
                                                                        makePtEtaPhiE(fathFilterJets_genPt[1]             , fathFilterJets_pt[1], fathFilterJets_eta[1], fathFilterJets_phi[1], fathFilterJets_e[1]) : p4Zero;
            const TLorentzVector& fathFilterJets_p4Gen_2              = (nfathFilterJets > 2 && fathFilterJets_pt[2] > 0) && fathFilterJets_genPt[2] > 0 ? 
                                                                        makePtEtaPhiE(fathFilterJets_genPt[2]             , fathFilterJets_pt[2], fathFilterJets_eta[2], fathFilterJets_phi[2], fathFilterJets_e[2]) : p4Zero;
            const TLorentzVector& fathFilterJets_p4Reg_0              = makePtEtaPhiE(fathFilterJets_ptReg[0]             , fathFilterJets_pt[0], fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_e[0]);
            const TLorentzVector& fathFilterJets_p4Reg_1              = makePtEtaPhiE(fathFilterJets_ptReg[1]             , fathFilterJets_pt[1], fathFilterJets_eta[1], fathFilterJets_phi[1], fathFilterJets_e[1]);
            const TLorentzVector& fathFilterJets_p4Reg_2              = (nfathFilterJets > 2 && fathFilterJets_pt[2] > 0) ? 
                                                                        makePtEtaPhiE(fathFilterJets_ptReg[2]             , fathFilterJets_pt[2], fathFilterJets_eta[2], fathFilterJets_phi[2], fathFilterJets_e[2]) : p4Zero;
            const TLorentzVector& fathFilterJets_p4Reg_res_j_up_0     = makePtEtaPhiE(fathFilterJets_ptReg_res_j_up[0]    , fathFilterJets_pt[0], fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_e[0]);
            const TLorentzVector& fathFilterJets_p4Reg_res_j_up_1     = makePtEtaPhiE(fathFilterJets_ptReg_res_j_up[1]    , fathFilterJets_pt[1], fathFilterJets_eta[1], fathFilterJets_phi[1], fathFilterJets_e[1]);
            const TLorentzVector& fathFilterJets_p4Reg_res_j_up_2     = (nfathFilterJets > 2 && fathFilterJets_pt[2] > 0) ? 
                                                                        makePtEtaPhiE(fathFilterJets_ptReg_res_j_up[2]    , fathFilterJets_pt[2], fathFilterJets_eta[2], fathFilterJets_phi[2], fathFilterJets_e[2]) : p4Zero;
            const TLorentzVector& fathFilterJets_p4Reg_res_j_down_0   = makePtEtaPhiE(fathFilterJets_ptReg_res_j_down[0]  , fathFilterJets_pt[0], fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_e[0]);
            const TLorentzVector& fathFilterJets_p4Reg_res_j_down_1   = makePtEtaPhiE(fathFilterJets_ptReg_res_j_down[1]  , fathFilterJets_pt[1], fathFilterJets_eta[1], fathFilterJets_phi[1], fathFilterJets_e[1]);
            const TLorentzVector& fathFilterJets_p4Reg_res_j_down_2   = (nfathFilterJets > 2 && fathFilterJets_pt[2] > 0) ? 
                                                                        makePtEtaPhiE(fathFilterJets_ptReg_res_j_down[2]  , fathFilterJets_pt[2], fathFilterJets_eta[2], fathFilterJets_phi[2], fathFilterJets_e[2]) : p4Zero;
            const TLorentzVector& fathFilterJets_p4Reg_scale_j_up_0   = makePtEtaPhiE(fathFilterJets_ptReg_scale_j_up[0]  , fathFilterJets_pt[0], fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_e[0]);
            const TLorentzVector& fathFilterJets_p4Reg_scale_j_up_1   = makePtEtaPhiE(fathFilterJets_ptReg_scale_j_up[1]  , fathFilterJets_pt[1], fathFilterJets_eta[1], fathFilterJets_phi[1], fathFilterJets_e[1]);
            const TLorentzVector& fathFilterJets_p4Reg_scale_j_up_2   = (nfathFilterJets > 2 && fathFilterJets_pt[2] > 0) ? 
                                                                        makePtEtaPhiE(fathFilterJets_ptReg_scale_j_up[2]  , fathFilterJets_pt[2], fathFilterJets_eta[2], fathFilterJets_phi[2], fathFilterJets_e[2]) : p4Zero;
            const TLorentzVector& fathFilterJets_p4Reg_scale_j_down_0 = makePtEtaPhiE(fathFilterJets_ptReg_scale_j_down[0], fathFilterJets_pt[0], fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_e[0]);
            const TLorentzVector& fathFilterJets_p4Reg_scale_j_down_1 = makePtEtaPhiE(fathFilterJets_ptReg_scale_j_down[1], fathFilterJets_pt[1], fathFilterJets_eta[1], fathFilterJets_phi[1], fathFilterJets_e[1]);
            const TLorentzVector& fathFilterJets_p4Reg_scale_j_down_2 = (nfathFilterJets > 2 && fathFilterJets_pt[2] > 0) ? 
                                                                        makePtEtaPhiE(fathFilterJets_ptReg_scale_j_down[2], fathFilterJets_pt[2], fathFilterJets_eta[2], fathFilterJets_phi[2], fathFilterJets_e[2]) : p4Zero;

            FatHptNorm             = (fathFilterJets_p4Norm_0             + fathFilterJets_p4Norm_1             + fathFilterJets_p4Norm_2            ).Pt();
            FatHptGen              = (fathFilterJets_p4Gen_0              + fathFilterJets_p4Gen_1              + fathFilterJets_p4Gen_2             ).Pt();
            FatHptReg              = (fathFilterJets_p4Reg_0              + fathFilterJets_p4Reg_1              + fathFilterJets_p4Reg_2             ).Pt();
            FatHptReg_res_j_up     = (fathFilterJets_p4Reg_res_j_up_0     + fathFilterJets_p4Reg_res_j_up_1     + fathFilterJets_p4Reg_res_j_up_2    ).Pt();
            FatHptReg_res_j_down   = (fathFilterJets_p4Reg_res_j_down_0   + fathFilterJets_p4Reg_res_j_down_1   + fathFilterJets_p4Reg_res_j_down_2  ).Pt();
            FatHptReg_scale_j_up   = (fathFilterJets_p4Reg_scale_j_up_0   + fathFilterJets_p4Reg_scale_j_up_1   + fathFilterJets_p4Reg_scale_j_up_2  ).Pt();
            FatHptReg_scale_j_down = (fathFilterJets_p4Reg_scale_j_down_0 + fathFilterJets_p4Reg_scale_j_down_1 + fathFilterJets_p4Reg_scale_j_down_2).Pt();

            FatHmassNorm             = (fathFilterJets_p4Norm_0             + fathFilterJets_p4Norm_1             + fathFilterJets_p4Norm_2            ).M();
            FatHmassGen              = (fathFilterJets_p4Gen_0              + fathFilterJets_p4Gen_1              + fathFilterJets_p4Gen_2             ).M();
            FatHmassReg              = (fathFilterJets_p4Reg_0              + fathFilterJets_p4Reg_1              + fathFilterJets_p4Reg_2             ).M();
            FatHmassReg_res_j_up     = (fathFilterJets_p4Reg_res_j_up_0     + fathFilterJets_p4Reg_res_j_up_1     + fathFilterJets_p4Reg_res_j_up_2    ).M();
            FatHmassReg_res_j_down   = (fathFilterJets_p4Reg_res_j_down_0   + fathFilterJets_p4Reg_res_j_down_1   + fathFilterJets_p4Reg_res_j_down_2  ).M();
            FatHmassReg_scale_j_up   = (fathFilterJets_p4Reg_scale_j_up_0   + fathFilterJets_p4Reg_scale_j_up_1   + fathFilterJets_p4Reg_scale_j_up_2  ).M();
            FatHmassReg_scale_j_down = (fathFilterJets_p4Reg_scale_j_down_0 + fathFilterJets_p4Reg_scale_j_down_1 + fathFilterJets_p4Reg_scale_j_down_2).M();
        } else {
            FatHptNorm                  = -99.;
            FatHptGen                   = -99.;
            FatHptReg                   = -99.;
            FatHptReg_res_j_up          = -99.;
            FatHptReg_res_j_down        = -99.;
            FatHptReg_scale_j_up        = -99.;
            FatHptReg_scale_j_down      = -99.;
            FatHmassNorm                = -99.;
            FatHmassGen                 = -99.;
            FatHmassReg                 = -99.;
            FatHmassReg_res_j_up        = -99.;
            FatHmassReg_res_j_down      = -99.;
            FatHmassReg_scale_j_up      = -99.;
            FatHmassReg_scale_j_down    = -99.;
        }


        MET_res_j_up          = metUnc_et[4];
        MET_res_j_down        = metUnc_et[5];
        MET_scale_j_up        = metUnc_et[2];
        MET_scale_j_down      = metUnc_et[3];
        METphi_res_j_up       = metUnc_phi[4];
        METphi_res_j_down     = metUnc_phi[5];
        METphi_scale_j_up     = metUnc_phi[2];
        METphi_scale_j_down   = metUnc_phi[3];

#ifdef PATCHMETTYPE1CORR
        TVector2 MET_original, MET_patched;
        MET_original.SetMagPhi(METtype1corr.et, METtype1corr.phi);
        MET_patched .SetMagPhi(METtype1diff.et, METtype1diff.phi);
        MET_patched         = MET_original + MET_patched;
        METtype1corr.et     = MET_patched.Mod();
        METtype1corr.phi    = MET_patched.Phi();
        METtype1corr.sumet  = METtype1corr.sumet + METtype1diff.sumet;
        
        MET_original.SetMagPhi(MET_res_j_up, METphi_res_j_up);
        MET_patched .SetMagPhi(METtype1diff.et, METtype1diff.phi);
        MET_patched         = MET_original + MET_patched;
        MET_res_j_up        = MET_patched.Mod();
        METphi_res_j_up     = MET_patched.Phi();
        
        MET_original.SetMagPhi(MET_res_j_down, METphi_res_j_down);
        MET_patched .SetMagPhi(METtype1diff.et, METtype1diff.phi);
        MET_patched         = MET_original + MET_patched;
        MET_res_j_down      = MET_patched.Mod();
        METphi_res_j_down   = MET_patched.Phi();
        
        MET_original.SetMagPhi(MET_scale_j_up, METphi_scale_j_up);
        MET_patched .SetMagPhi(METtype1diff.et, METtype1diff.phi);
        MET_patched         = MET_original + MET_patched;
        MET_scale_j_up      = MET_patched.Mod();
        METphi_scale_j_up   = MET_patched.Phi();
        
        MET_original.SetMagPhi(MET_scale_j_down, METphi_scale_j_down);
        MET_patched .SetMagPhi(METtype1diff.et, METtype1diff.phi);
        MET_patched         = MET_original + MET_patched;
        MET_scale_j_down    = MET_patched.Mod();
        METphi_scale_j_down = MET_patched.Phi();
#endif

        MET_xyshift           = shift_met_by_Nvtx(METtype1corr.et, METtype1corr.phi, nPVs, EVENT_run);
        METphi_xyshift        = shift_metphi_by_Nvtx(METtype1corr.et, METtype1corr.phi, nPVs, EVENT_run);

        HTeta25 = 0.;  // sum pt instead of sum et
        HTeta45 = 0.;  // sum pt instead of sum et
        for (Int_t ihj = 0; ihj < 2; ihj++) {
            if (hJet_ptReg[ihj] > 20. && fabs(hJet_eta[ihj]) < 2.5) {
                HTeta25 += hJet_ptReg[ihj];
            }
            if (hJet_ptReg[ihj] > 25. && fabs(hJet_eta[ihj]) < 4.5) {
                HTeta45 += hJet_ptReg[ihj];
            }
        }
        for (Int_t iaj = 0; iaj < naJets; iaj++) {
            if (aJet_pt[iaj] > 20. && fabs(aJet_eta[iaj]) < 2.5) {
                HTeta25 += aJet_pt[iaj];
            }
            if (aJet_pt[iaj] > 25. && fabs(aJet_eta[iaj]) < 4.5) {
                HTeta45 += aJet_pt[iaj];
            }
        }
        
        naJets_Znn = 0;
        naJets_Znn_res_j_up = 0;
        naJets_Znn_res_j_down = 0;
        naJets_Znn_scale_j_up = 0;
        naJets_Znn_scale_j_down = 0;
        naCtrJets_pt20to25_Znn = 0;
        naCtrJets_pt20to25_Znn_res_j_up = 0;
        naCtrJets_pt20to25_Znn_res_j_down = 0;
        naCtrJets_pt20to25_Znn_scale_j_up = 0;
        naCtrJets_pt20to25_Znn_scale_j_down = 0;
        mindPhiMETJet_dPhi = 999.;
        mindPhiMETJet_dPhi_res_j_up = 999.;
        mindPhiMETJet_dPhi_res_j_down = 999.;
        mindPhiMETJet_dPhi_scale_j_up = 999.;
        mindPhiMETJet_dPhi_scale_j_down = 999.;
        mindPhiMETCtrJet_dPhi = 999.;
        mindPhiMETCtrJet_dPhi_res_j_up = 999.;
        mindPhiMETCtrJet_dPhi_res_j_down = 999.;
        mindPhiMETCtrJet_dPhi_scale_j_up = 999.;
        mindPhiMETCtrJet_dPhi_scale_j_down = 999.;
        mindRHJet_dR = 999.;
        mindRHJet_dR_res_j_up = 999.;
        mindRHJet_dR_res_j_down = 999.;
        mindRHJet_dR_scale_j_up = 999.;
        mindRHJet_dR_scale_j_down = 999.;
        mindRHJet_csv_nominal = -10.;
        
        /// Loop on Higgs jets
        for (Int_t ihj = 0; ihj < 2; ihj++) {
            // min dPhi with central jets
            if (hJet_ptReg[ihj] > 20. && fabs(hJet_eta[ihj]) < 2.5) {
                float dPhi = fabs(deltaPhi(METtype1corr.phi, hJet_phi[ihj]));
                if (mindPhiMETCtrJet_dPhi > dPhi) {
                    mindPhiMETCtrJet_dPhi = dPhi;
                    mindPhiMETCtrJet_pt = hJet_ptReg[ihj];
                }
            }
            if (hJet_ptReg_res_j_up[ihj] > 20. && fabs(hJet_eta[ihj]) < 2.5) {
                float dPhi = fabs(deltaPhi(METphi_res_j_up, hJet_phi[ihj]));
                if (mindPhiMETCtrJet_dPhi_res_j_up > dPhi) {
                    mindPhiMETCtrJet_dPhi_res_j_up = dPhi;
                }
            }
            if (hJet_ptReg_res_j_down[ihj] > 20. && fabs(hJet_eta[ihj]) < 2.5) {
                float dPhi = fabs(deltaPhi(METphi_res_j_down, hJet_phi[ihj]));
                if (mindPhiMETCtrJet_dPhi_res_j_down > dPhi) {
                    mindPhiMETCtrJet_dPhi_res_j_down = dPhi;
                }
            }
            if (hJet_ptReg_scale_j_up[ihj] > 20. && fabs(hJet_eta[ihj]) < 2.5) {
                float dPhi = fabs(deltaPhi(METphi_scale_j_up, hJet_phi[ihj]));
                if (mindPhiMETCtrJet_dPhi_scale_j_up > dPhi) {
                    mindPhiMETCtrJet_dPhi_scale_j_up = dPhi;
                }
            }
            if (hJet_ptReg_scale_j_down[ihj] > 20. && fabs(hJet_eta[ihj]) < 2.5) {
                float dPhi = fabs(deltaPhi(METphi_scale_j_down, hJet_phi[ihj]));
                if (mindPhiMETCtrJet_dPhi_scale_j_down > dPhi) {
                    mindPhiMETCtrJet_dPhi_scale_j_down = dPhi;
                }
            }
            
            // min dPhi with all jets
            if (hJet_ptReg[ihj] > 25. && fabs(hJet_eta[ihj]) < 4.5 && hJet_id[ihj] && hJet_puJetIdL[ihj]>0.) {
                float dPhi = fabs(deltaPhi(METtype1corr.phi, hJet_phi[ihj]));
                if (mindPhiMETJet_dPhi > dPhi) {
                    mindPhiMETJet_dPhi = dPhi;
                    mindPhiMETJet_pt = hJet_ptReg[ihj];
                }
            }
            if (hJet_ptReg_res_j_up[ihj] > 25. && fabs(hJet_eta[ihj]) < 4.5 && hJet_id[ihj] && hJet_puJetIdL[ihj]>0.) {
                float dPhi = fabs(deltaPhi(METphi_res_j_up, hJet_phi[ihj]));
                if (mindPhiMETJet_dPhi_res_j_up > dPhi) {
                    mindPhiMETJet_dPhi_res_j_up = dPhi;
                }
            }
            if (hJet_ptReg_res_j_down[ihj] > 25. && fabs(hJet_eta[ihj]) < 4.5 && hJet_id[ihj] && hJet_puJetIdL[ihj]>0.) {
                float dPhi = fabs(deltaPhi(METphi_res_j_down, hJet_phi[ihj]));
                if (mindPhiMETJet_dPhi_res_j_down > dPhi) {
                    mindPhiMETJet_dPhi_res_j_down = dPhi;
                }
            }
            if (hJet_ptReg_scale_j_up[ihj] > 25. && fabs(hJet_eta[ihj]) < 4.5 && hJet_id[ihj] && hJet_puJetIdL[ihj]>0.) {
                float dPhi = fabs(deltaPhi(METphi_scale_j_up, hJet_phi[ihj]));
                if (mindPhiMETJet_dPhi_scale_j_up > dPhi) {
                    mindPhiMETJet_dPhi_scale_j_up = dPhi;
                }
            }
            if (hJet_ptReg_scale_j_down[ihj] > 25. && fabs(hJet_eta[ihj]) < 4.5 && hJet_id[ihj] && hJet_puJetIdL[ihj]>0.) {
                float dPhi = fabs(deltaPhi(METphi_scale_j_down, hJet_phi[ihj]));
                if (mindPhiMETJet_dPhi_scale_j_down > dPhi) {
                    mindPhiMETJet_dPhi_scale_j_down = dPhi;
                }
            }
        }

        /// Loop on additional jets
        for (Int_t iaj = 0; iaj < naJets; iaj++) {
            // min dPhi with central jets
            if (aJet_pt[iaj] > 20. && fabs(aJet_eta[iaj]) < 2.5) {
                if (aJet_pt[iaj] <= 25.)  naCtrJets_pt20to25_Znn++;
                float dPhi = fabs(deltaPhi(METtype1corr.phi, aJet_phi[iaj]));
                if (mindPhiMETCtrJet_dPhi > dPhi) {
                    mindPhiMETCtrJet_dPhi = dPhi;
                    mindPhiMETCtrJet_pt = aJet_pt[iaj];
                }
            }
            if (aJet_pt_res_j_up[iaj] > 20. && fabs(aJet_eta[iaj]) < 2.5) {
                if (aJet_pt_res_j_up[iaj] <= 25.)  naCtrJets_pt20to25_Znn_res_j_up++;
                float dPhi = fabs(deltaPhi(METphi_res_j_up, aJet_phi[iaj]));
                if (mindPhiMETCtrJet_dPhi_res_j_up > dPhi) {
                    mindPhiMETCtrJet_dPhi_res_j_up = dPhi;
                }
            }
            if (aJet_pt_res_j_down[iaj] > 20. && fabs(aJet_eta[iaj]) < 2.5) {
                if (aJet_pt_res_j_down[iaj] <= 25.)  naCtrJets_pt20to25_Znn_res_j_down++;
                float dPhi = fabs(deltaPhi(METphi_res_j_down, aJet_phi[iaj]));
                if (mindPhiMETCtrJet_dPhi_res_j_down > dPhi) {
                    mindPhiMETCtrJet_dPhi_res_j_down = dPhi;
                }
            }
            if (aJet_pt_scale_j_up[iaj] > 20. && fabs(aJet_eta[iaj]) < 2.5) {
                if (aJet_pt_scale_j_up[iaj] <= 25.)  naCtrJets_pt20to25_Znn_scale_j_up++;
                float dPhi = fabs(deltaPhi(METphi_scale_j_up, aJet_phi[iaj]));
                if (mindPhiMETCtrJet_dPhi_scale_j_up > dPhi) {
                    mindPhiMETCtrJet_dPhi_scale_j_up = dPhi;
                }
            }
            if (aJet_pt_scale_j_down[iaj] > 20. && fabs(aJet_eta[iaj]) < 2.5) {
                if (aJet_pt_scale_j_down[iaj] <= 25.)  naCtrJets_pt20to25_Znn_scale_j_down++;
                float dPhi = fabs(deltaPhi(METphi_scale_j_down, aJet_phi[iaj]));
                if (mindPhiMETCtrJet_dPhi_scale_j_down > dPhi) {
                    mindPhiMETCtrJet_dPhi_scale_j_down = dPhi;
                }
            }
            
            // min dPhi with all jets
            if (aJet_pt[iaj] > 25. && fabs(aJet_eta[iaj]) < 4.5 && aJet_id[iaj] && aJet_puJetIdL[iaj]>0.) {
                naJets_Znn++;
                float dPhi = fabs(deltaPhi(METtype1corr.phi, aJet_phi[iaj]));
                if (mindPhiMETJet_dPhi > dPhi) {
                    mindPhiMETJet_dPhi = dPhi;
                    mindPhiMETJet_pt = aJet_pt[iaj];
                }
            }
            if (aJet_pt_res_j_up[iaj] > 25. && fabs(aJet_eta[iaj]) < 4.5 && aJet_id[iaj] && aJet_puJetIdL[iaj]>0.) {
                naJets_Znn_res_j_up++;
                float dPhi = fabs(deltaPhi(METphi_res_j_up, aJet_phi[iaj]));
                if (mindPhiMETJet_dPhi_res_j_up > dPhi) {
                    mindPhiMETJet_dPhi_res_j_up = dPhi;
                }
            }
            if (aJet_pt_res_j_down[iaj] > 25. && fabs(aJet_eta[iaj]) < 4.5 && aJet_id[iaj] && aJet_puJetIdL[iaj]>0.) {
                naJets_Znn_res_j_down++;
                float dPhi = fabs(deltaPhi(METphi_res_j_down, aJet_phi[iaj]));
                if (mindPhiMETJet_dPhi_res_j_down > dPhi) {
                    mindPhiMETJet_dPhi_res_j_down = dPhi;
                }
            }
            if (aJet_pt_scale_j_up[iaj] > 25. && fabs(aJet_eta[iaj]) < 4.5 && aJet_id[iaj] && aJet_puJetIdL[iaj]>0.) {
                naJets_Znn_scale_j_up++;
                float dPhi = fabs(deltaPhi(METphi_scale_j_up, aJet_phi[iaj]));
                if (mindPhiMETJet_dPhi_scale_j_up > dPhi) {
                    mindPhiMETJet_dPhi_scale_j_up = dPhi;
                }
            }
            if (aJet_pt_scale_j_down[iaj] > 25. && fabs(aJet_eta[iaj]) < 4.5 && aJet_id[iaj] && aJet_puJetIdL[iaj]>0.) {
                naJets_Znn_scale_j_down++;
                float dPhi = fabs(deltaPhi(METphi_scale_j_down, aJet_phi[iaj]));
                if (mindPhiMETJet_dPhi_scale_j_down > dPhi) {
                    mindPhiMETJet_dPhi_scale_j_down = dPhi;
                }
            }

            // min dR with additional jets
            if (aJet_pt[iaj] > 25. && fabs(aJet_eta[iaj]) < 4.5 && aJet_id[iaj] && aJet_puJetIdL[iaj]>0.) {
                float dR = deltaR(aJet_eta[iaj], aJet_phi[iaj], H.eta, H.phi);
                if (mindRHJet_dR > dR) {
                    mindRHJet_dR = dR;
                    mindRHJet_pt = aJet_pt[iaj];
                    mindRHJet_csv_nominal = aJet_csv_nominal[iaj];
                }
            }
            if (aJet_pt_res_j_up[iaj] > 25. && fabs(aJet_eta[iaj]) < 4.5 && aJet_id[iaj] && aJet_puJetIdL[iaj]>0.) {
                float dR = deltaR(aJet_eta[iaj], aJet_phi[iaj], H.eta, H.phi);
                if (mindRHJet_dR_res_j_up > dR) {
                    mindRHJet_dR_res_j_up = dR;
                    mindRHJet_pt_res_j_up = aJet_pt_res_j_up[iaj];
                }
            }
            if (aJet_pt_res_j_down[iaj] > 25. && fabs(aJet_eta[iaj]) < 4.5 && aJet_id[iaj] && aJet_puJetIdL[iaj]>0.) {
                float dR = deltaR(aJet_eta[iaj], aJet_phi[iaj], H.eta, H.phi);
                if (mindRHJet_dR_res_j_down > dR) {
                    mindRHJet_dR_res_j_down = dR;
                    mindRHJet_pt_res_j_down = aJet_pt_res_j_down[iaj];
                }
            }
            if (aJet_pt_scale_j_up[iaj] > 25. && fabs(aJet_eta[iaj]) < 4.5 && aJet_id[iaj] && aJet_puJetIdL[iaj]>0.) {
                float dR = deltaR(aJet_eta[iaj], aJet_phi[iaj], H.eta, H.phi);
                if (mindRHJet_dR_scale_j_up > dR) {
                    mindRHJet_dR_scale_j_up = dR;
                    mindRHJet_pt_scale_j_up = aJet_pt_scale_j_up[iaj];
                }
            }
            if (aJet_pt_scale_j_down[iaj] > 25. && fabs(aJet_eta[iaj]) < 4.5 && aJet_id[iaj] && aJet_puJetIdL[iaj]>0.) {
                float dR = deltaR(aJet_eta[iaj], aJet_phi[iaj], H.eta, H.phi);
                if (mindRHJet_dR_scale_j_down > dR) {
                    mindRHJet_dR_scale_j_down = dR;
                    mindRHJet_pt_scale_j_down = aJet_pt_scale_j_down[iaj];
                }
            }
        }

        nalep_Znn = 0;
        nalep_pt5_Znn = 0;
        /// Loop on V leptons
        for (Int_t ivl = 0; ivl < nvlep; ivl++) {
            if (vLepton_pt[ivl] > 15. && vLepton_pfCombRelIso[ivl] < 0.15 && (vLepton_id95[ivl]==7 || vLepton_vbtf[ivl]==1)) {
                nalep_Znn++;
            }
            if (vLepton_pt[ivl] > 5.  && vLepton_pfCombRelIso[ivl] < 0.3  && (vLepton_id95[ivl]==7 || vLepton_vbtf[ivl]==1)) {
                if (deltaR(vLepton_eta[ivl], vLepton_phi[ivl], hJet_eta[0], hJet_phi[0]) > 0.3 && deltaR(vLepton_eta[ivl], vLepton_phi[ivl], hJet_eta[1], hJet_phi[1]) > 0.3) {
                    nalep_pt5_Znn++;
                }
            }
        }

        /// Loop on additional leptons
        for (Int_t ial = 0; ial < nalep; ial++) {
            if (aLepton_pt[ial] > 15. && aLepton_pfCombRelIso[ial] < 0.15 && (aLepton_id95[ial]==7 || aLepton_vbtf[ial]==1)) {
                nalep_Znn++;
            }
            if (aLepton_pt[ial] > 5.  && aLepton_pfCombRelIso[ial] < 0.3  && (aLepton_id95[ial]==7 || aLepton_vbtf[ial]==1)) {
                if (deltaR(aLepton_eta[ial], aLepton_phi[ial], hJet_eta[0], hJet_phi[0]) > 0.3 && deltaR(aLepton_eta[ial], aLepton_phi[ial], hJet_eta[1], hJet_phi[1]) > 0.3) {
                    nalep_pt5_Znn++;
                }
            }
        }

        weightTrigMu = 0.;
        weightTrigDLP_ElecMay10 = 0.;
        weightTrigDLP_Elec = 0.;
        double TnPMuTriggerRecoNew[70]={0.990, 0.995, 0.997, 1.007, 1.001, 1.005, 1.004, 1.002, 1.004, 
            0.996, 0.993, 0.997, 0.995, 1.001, 0.995, 0.997, 1.002, 1.000, 1.006, 0.997, 0.989, 0.995, 0.995, 
            0.998, 0.998, 1.001, 1.002, 1.002, 1.002, 0.999, 0.981, 0.994, 0.991, 0.999, 0.999, 1.000, 1.003, 
            1.000, 1.004, 0.995, 0.989, 0.996, 0.993, 0.999, 0.997, 1.000, 1.000, 0.998, 1.000, 0.992, 0.990, 
            0.997, 0.991, 0.998, 0.996, 0.999, 1.000, 0.998, 0.998, 0.989, 0.993, 0.996, 0.991, 0.999, 0.997, 
            1.001, 1.000, 1.000, 0.998, 0.994} ;
        double TnPMuTriggerOnlyIso22[70]={0.63053367,0.835756,0.873449,0.905668,0.919801,0.924853,0.896946,
            0.856479,0.837122,0.630819,0.634994,0.846767,0.870174,0.904408,0.920462,0.924441,0.898895,0.860216,
            0.843505,0.638490,0.628831,0.844661,0.868911,0.893475,0.915531,0.918174,0.893783,0.856110,0.841024,
            0.650392,0.645827,0.840136,0.861151,0.885003,0.909036,0.908888,0.889770,0.850571,0.851518,0.644835,
            0.639849,0.844856,0.864901,0.885730,0.909001,0.908483,0.889374,0.852705,0.845170,0.647222,0.642904,
            0.849120,0.864003,0.888698,0.908981,0.907638,0.885068,0.855384,0.850608,0.656233,0.617986,0.841628,
            0.854240,0.874969,0.896921,0.900516,0.874691,0.844891,0.840691,0.637464};
        double TnpElecRecoNewEps[20]={1.001,1.006, 1.002, 1.009, 1.004, 0.999, 0.993, 1.013, 1.000, 0.995, 0.994, 1.011, 
       1.020, 0.994, 0.997, 1.017, 1.029, 0.980, 0.999, 1.042}; 
        double TnpElecIDNewEps[20]={0.934,1.002,0.996,0.937,0.954,0.991,0.996,0.961,0.966,0.995,0.989,0.965,
       0.967,0.994,0.992,0.973,0.984,0.992,0.986,0.981};
        double TnpElecTriggerNewMay10[20]={0.9558,0.9675,0.9736,0.9335,0.9661,0.9746,0.9712,0.9704,0.9699,0.9781,
       0.9777,0.9735,0.9709,0.9769,0.9757,0.9743,0.9735,0.9781,0.9766,0.982};
        double TnpElecTriggerNewLP[20]={0.87424,0.97113,0.96896,0.85186,0.96298,0.97755,0.97572,0.96233,0.97202,
       0.97886,0.97622,0.96992,0.97202,0.97818,0.97640,0.97268,0.97470,0.98197,0.98075,0.97461};

        if (Vtype==2) {
            for(int i=0;i<7;i++) {
                double ptTnPL=20.+5.*i; 
                double ptTnPH=25.+5.*i;
                if(i==6) ptTnPH=10000.; 
                if(vLepton_pt[0]>ptTnPL && vLepton_pt[0]<ptTnPH && vLepton_eta[0]>=-2.4 && vLepton_eta[0]<-2.1) weightTrigMu=TnPMuTriggerOnlyIso22[10*i]*TnPMuTriggerRecoNew[10*i];
                if(vLepton_pt[0]>ptTnPL && vLepton_pt[0]<ptTnPH && vLepton_eta[0]>=-2.1 && vLepton_eta[0]<-1.5) weightTrigMu=TnPMuTriggerOnlyIso22[10*i+1]*TnPMuTriggerRecoNew[10*i+1];
                if(vLepton_pt[0]>ptTnPL && vLepton_pt[0]<ptTnPH && vLepton_eta[0]>=-1.5 && vLepton_eta[0]<-1.) weightTrigMu=TnPMuTriggerOnlyIso22[10*i+2]*TnPMuTriggerRecoNew[10*i+2];
                if(vLepton_pt[0]>ptTnPL && vLepton_pt[0]<ptTnPH && vLepton_eta[0]>=-1. && vLepton_eta[0]<-0.5) weightTrigMu=TnPMuTriggerOnlyIso22[10*i+3]*TnPMuTriggerRecoNew[10*i+3];
                if(vLepton_pt[0]>ptTnPL && vLepton_pt[0]<ptTnPH && vLepton_eta[0]>=-0.5 && vLepton_eta[0]<0.) weightTrigMu=TnPMuTriggerOnlyIso22[10*i+4]*TnPMuTriggerRecoNew[10*i+4];
                if(vLepton_pt[0]>ptTnPL && vLepton_pt[0]<ptTnPH && vLepton_eta[0]>=0. && vLepton_eta[0]<0.5) weightTrigMu=TnPMuTriggerOnlyIso22[10*i+5]*TnPMuTriggerRecoNew[10*i+5];
                if(vLepton_pt[0]>ptTnPL && vLepton_pt[0]<ptTnPH && vLepton_eta[0]>=0.5 && vLepton_eta[0]<1.) weightTrigMu=TnPMuTriggerOnlyIso22[10*i+6]*TnPMuTriggerRecoNew[10*i+6];
                if(vLepton_pt[0]>ptTnPL && vLepton_pt[0]<ptTnPH && vLepton_eta[0]>=1. && vLepton_eta[0]<1.5) weightTrigMu=TnPMuTriggerOnlyIso22[10*i+7]*TnPMuTriggerRecoNew[10*i+7];
                if(vLepton_pt[0]>ptTnPL && vLepton_pt[0]<ptTnPH && vLepton_eta[0]>=1.5 && vLepton_eta[0]<2.1) weightTrigMu=TnPMuTriggerOnlyIso22[10*i+8]*TnPMuTriggerRecoNew[10*i+8];
                if(vLepton_pt[0]>ptTnPL && vLepton_pt[0]<ptTnPH && vLepton_eta[0]>=2.1 && vLepton_eta[0]<2.4) weightTrigMu=TnPMuTriggerOnlyIso22[10*i+9]*TnPMuTriggerRecoNew[10*i+9];
            }
        } else if (Vtype==3) {
            for(int i=0;i<5;i++){
                double ptTnPL=30.+5.*i; 
                double ptTnPH=35.+5.*i;
                if(i==4) ptTnPH=10000.;
                if(vLepton_pt[0]>ptTnPL && vLepton_pt[0]<ptTnPH && vLepton_eta[0]>=-2.5 && vLepton_eta[0]<-1.5) weightTrigDLP_ElecMay10=TnpElecTriggerNewMay10[4*i]*TnpElecIDNewEps[4*i]*TnpElecRecoNewEps[4*i];
                if(vLepton_pt[0]>ptTnPL && vLepton_pt[0]<ptTnPH && vLepton_eta[0]>=-1.5 && vLepton_eta[0]<0.) weightTrigDLP_ElecMay10=TnpElecTriggerNewMay10[4*i+1]*TnpElecIDNewEps[4*i+1]*TnpElecRecoNewEps[4*i+1];
                if(vLepton_pt[0]>ptTnPL && vLepton_pt[0]<ptTnPH && vLepton_eta[0]>=0. && vLepton_eta[0]<1.5) weightTrigDLP_ElecMay10=TnpElecTriggerNewMay10[4*i+2]*TnpElecIDNewEps[4*i+2]*TnpElecRecoNewEps[4*i+2];
                if(vLepton_pt[0]>ptTnPL && vLepton_pt[0]<ptTnPH && vLepton_eta[0]>=1.5 && vLepton_eta[0]<2.5) weightTrigDLP_ElecMay10=TnpElecTriggerNewMay10[4*i+3]*TnpElecIDNewEps[4*i+3]*TnpElecRecoNewEps[4*i+3];

                if(vLepton_pt[0]>ptTnPL && vLepton_pt[0]<ptTnPH && vLepton_eta[0]>=-2.5 && vLepton_eta[0]<-1.5) weightTrigDLP_Elec=TnpElecTriggerNewLP[4*i]*TnpElecIDNewEps[4*i]*TnpElecRecoNewEps[4*i];
                if(vLepton_pt[0]>ptTnPL && vLepton_pt[0]<ptTnPH && vLepton_eta[0]>=-1.5 && vLepton_eta[0]<0.) weightTrigDLP_Elec=TnpElecTriggerNewLP[4*i+1]*TnpElecIDNewEps[4*i+1]*TnpElecRecoNewEps[4*i+1];
                if(vLepton_pt[0]>ptTnPL && vLepton_pt[0]<ptTnPH && vLepton_eta[0]>=0. && vLepton_eta[0]<1.5) weightTrigDLP_Elec=TnpElecTriggerNewLP[4*i+2]*TnpElecIDNewEps[4*i+2]*TnpElecRecoNewEps[4*i+2];
                if(vLepton_pt[0]>ptTnPL && vLepton_pt[0]<ptTnPH && vLepton_eta[0]>=1.5 && vLepton_eta[0]<2.5) weightTrigDLP_Elec=TnpElecTriggerNewLP[4*i+3]*TnpElecIDNewEps[4*i+3]*TnpElecRecoNewEps[4*i+3];
            }
        }

        /// SemiLeptonic Top
        TopLep_mass = -99.;
        TopLep_pt = -99.;
        TopLep_eta = -99.;
        TopLep_phi = -99.;
        TopLep_massReg = -99.;
        TopLep_ptReg = -99.;
        TopLep_jetPt = -99.;
        TopLep_csv_nominal = -10.;
        if ((Vtype==2 || Vtype==3) && vLepton_pt[0]>0.) {
            TLorentzVector Top_lep, Top_MET, Top_jet, Top_jetReg;
            Top_MET.SetPtEtaPhiE(METtype1corr.et, 0., METtype1corr.phi, METtype1corr.et);

            TVector3 Top_lep_temp;
            Top_lep_temp.SetPtEtaPhi(vLepton_pt[0], vLepton_eta[0], vLepton_phi[0]);
            float Top_lep_mass = (Vtype==2) ? 0.105 : 0.005;
            float Top_lep_e = TMath::Sqrt(Top_lep_temp.Mag2() + Top_lep_mass * Top_lep_mass);
            Top_lep.SetPtEtaPhiE(vLepton_pt[0], vLepton_eta[0], vLepton_phi[0], Top_lep_e);
            TLorentzVector Top_nu = getNu4Momentum(Top_lep, Top_MET);
            //TLorentzVector Top_W = Top_lep + Top_nu;

            float bestCSVhJet = -9999.;
            int bestCSVhJet_id = -99;
            // Pick the best b-tagged Higgs jet
            for (Int_t ihj = 0; ihj < 2; ihj++) {
                if (hJet_pt[ihj] > 30. && fabs(hJet_eta[ihj]) < 2.5 && hJet_csv_nominal[ihj]>0.) {
                    if (bestCSVhJet < hJet_csv_nominal[ihj]) {
                        bestCSVhJet = hJet_csv_nominal[ihj];
                        bestCSVhJet_id = ihj;
                    }
                }
            }
            float bestCSVaJet = -9999.;
            int bestCSVaJet_id = -99;
            // Pick the best b-tagged additional jet
            for (Int_t iaj = 0; iaj < naJets; iaj++) {
                if (aJet_pt[iaj] > 30. && fabs(aJet_eta[iaj]) < 2.5 && aJet_csv_nominal[iaj]>0.) {
                    if (bestCSVaJet < aJet_csv_nominal[iaj]) {
                        bestCSVaJet = aJet_csv_nominal[iaj];
                        bestCSVaJet_id = iaj;
                    }
                }
            }
            int use_hJet_or_aJet = -1;
            double diffCSVCut = 0.35;
            if (bestCSVhJet_id == -99 && bestCSVaJet_id == -99)  // should not happen...
                use_hJet_or_aJet = -1;
            else if (bestCSVaJet_id == -99)
                use_hJet_or_aJet = 0;
            else if (bestCSVhJet_id == -99)
                use_hJet_or_aJet = 1;
            else {
                double diffCSV = (bestCSVhJet - bestCSVaJet) / (bestCSVhJet + bestCSVaJet);
                if (fabs(diffCSV) > diffCSVCut) {
                    if (diffCSV > 0.)
                        use_hJet_or_aJet = 0;
                    else
                        use_hJet_or_aJet = 1;
                } else {  // just use hJet then...
                    use_hJet_or_aJet = 0;
                }
            }

            if (use_hJet_or_aJet != -1) {
                if (use_hJet_or_aJet == 1) {  // aJet
                    Top_jet.SetPtEtaPhiE(aJet_pt[bestCSVaJet_id], aJet_eta[bestCSVaJet_id], aJet_phi[bestCSVaJet_id], aJet_e[bestCSVaJet_id]);
                    Top_jetReg = Top_jet;
                    TopLep_jetPt = aJet_pt[bestCSVaJet_id];
                    TopLep_csv_nominal = aJet_csv_nominal[bestCSVaJet_id];
                } else {  // hJet
                    Top_jet.SetPtEtaPhiE(hJet_pt[bestCSVhJet_id], hJet_eta[bestCSVhJet_id], hJet_phi[bestCSVhJet_id], hJet_e[bestCSVhJet_id]);
                    Top_jetReg.SetPtEtaPhiE(hJet_ptReg[bestCSVhJet_id], hJet_eta[bestCSVhJet_id], hJet_phi[bestCSVhJet_id], hJet_e[bestCSVhJet_id] * hJet_ptReg[bestCSVhJet_id] / hJet_pt[bestCSVhJet_id]);
                    TopLep_jetPt = hJet_pt[bestCSVhJet_id];
                    TopLep_csv_nominal = hJet_csv_nominal[bestCSVhJet_id];
                }
                TLorentzVector Top = Top_lep + Top_nu + Top_jet;
                TopLep_mass = Top.M();
                TopLep_pt = Top.Pt();
                TopLep_eta = Top.Eta();
                TopLep_phi = Top.Phi();
                TLorentzVector TopReg = Top_lep + Top_nu + Top_jetReg;
                TopLep_massReg = TopReg.M();
                TopLep_ptReg = TopReg.Pt();
                TopLep_nuPt = Top_nu.Pt();
                TopLep_nuEta = (TopLep_nuPt > 0.) ? Top_nu.Eta() : -99.;
                TopLep_nuPhi = Top_nu.Phi();
            }

        } else if (Vtype==4 && nalep>0 && aLepton_pt[0]>0.) {
            // Take the additional lepton (perhaps also consider the closest jet to get tau?)
            TLorentzVector Top_lep, Top_MET, Top_jet, Top_jetReg;
            Top_MET.SetPtEtaPhiE(METtype1corr.et, 0., METtype1corr.phi, METtype1corr.et);

            TVector3 Top_lep_temp;
            Top_lep_temp.SetPtEtaPhi(aLepton_pt[0], aLepton_eta[0], aLepton_phi[0]);
            float Top_lep_mass = (abs(aLepton_type[0])==13) ? 0.105 : 0.005;
            float Top_lep_e = TMath::Sqrt(Top_lep_temp.Mag2() + Top_lep_mass * Top_lep_mass);
            Top_lep.SetPtEtaPhiE(aLepton_pt[0], aLepton_eta[0], aLepton_phi[0], Top_lep_e);
            TLorentzVector Top_nu = getNu4Momentum(Top_lep, Top_MET);
            //TLorentzVector Top_W = Top_lep + Top_nu;

            float bestCSVhJet = -9999.;
            int bestCSVhJet_id = -99;
            // Pick the best b-tagged Higgs jet
            for (Int_t ihj = 0; ihj < 2; ihj++) {
                if (hJet_pt[ihj] > 30. && fabs(hJet_eta[ihj]) < 2.5 && hJet_csv_nominal[ihj]>0.) {
                    if (bestCSVhJet < hJet_csv_nominal[ihj]) {
                        bestCSVhJet = hJet_csv_nominal[ihj];
                        bestCSVhJet_id = ihj;
                    }
                }
            }
            float bestCSVaJet = -9999.;
            int bestCSVaJet_id = -99;
            // Pick the best b-tagged additional jet
            for (Int_t iaj = 0; iaj < naJets; iaj++) {
                if (aJet_pt[iaj] > 30. && fabs(aJet_eta[iaj]) < 2.5 && aJet_csv_nominal[iaj]>0.) {
                    if (bestCSVaJet < aJet_csv_nominal[iaj]) {
                        bestCSVaJet = aJet_csv_nominal[iaj];
                        bestCSVaJet_id = iaj;
                    }
                }
            }
            int use_hJet_or_aJet = -1;
            double diffCSVCut = 0.35;
            if (bestCSVhJet_id == -99 && bestCSVaJet_id == -99)  // should not happen...
                use_hJet_or_aJet = -1;
            else if (bestCSVaJet_id == -99)
                use_hJet_or_aJet = 0;
            else if (bestCSVhJet_id == -99)
                use_hJet_or_aJet = 1;
            else {
                double diffCSV = (bestCSVhJet - bestCSVaJet) / (bestCSVhJet + bestCSVaJet);
                if (fabs(diffCSV) > diffCSVCut) {
                    if (diffCSV > 0.)
                        use_hJet_or_aJet = 0;
                    else
                        use_hJet_or_aJet = 1;
                } else {  // just use hJet then...
                    use_hJet_or_aJet = 0;
                }
            }

            if (use_hJet_or_aJet != -1) {
                if (use_hJet_or_aJet == 1) {  // aJet
                    Top_jet.SetPtEtaPhiE(aJet_pt[bestCSVaJet_id], aJet_eta[bestCSVaJet_id], aJet_phi[bestCSVaJet_id], aJet_e[bestCSVaJet_id]);
                    Top_jetReg = Top_jet;
                    TopLep_jetPt = aJet_pt[bestCSVaJet_id];
                    TopLep_csv_nominal = aJet_csv_nominal[bestCSVaJet_id];
                } else {  // hJet
                    Top_jet.SetPtEtaPhiE(hJet_pt[bestCSVhJet_id], hJet_eta[bestCSVhJet_id], hJet_phi[bestCSVhJet_id], hJet_e[bestCSVhJet_id]);
                    Top_jetReg.SetPtEtaPhiE(hJet_ptReg[bestCSVhJet_id], hJet_eta[bestCSVhJet_id], hJet_phi[bestCSVhJet_id], hJet_e[bestCSVhJet_id] * hJet_ptReg[bestCSVhJet_id] / hJet_pt[bestCSVhJet_id]);
                    TopLep_jetPt = hJet_pt[bestCSVhJet_id];
                    TopLep_csv_nominal = hJet_csv_nominal[bestCSVhJet_id];
                }
                TLorentzVector Top = Top_lep + Top_nu + Top_jet;
                TopLep_mass = Top.M();
                TopLep_pt = Top.Pt();
                TopLep_eta = Top.Eta();
                TopLep_phi = Top.Phi();
                TLorentzVector TopReg = Top_lep + Top_nu + Top_jetReg;
                TopLep_massReg = TopReg.M();
                TopLep_ptReg = TopReg.Pt();
                TopLep_nuPt = Top_nu.Pt();
                TopLep_nuEta = (TopLep_nuPt > 0.) ? Top_nu.Eta() : -99.;
                TopLep_nuPhi = Top_nu.Phi();
            }
        }

        /// Hadronic Top
        TopHad_mass = -99.;
        TopHad_pt = -99.;
        TopHad_eta = -99.;
        TopHad_phi = -99.;
        TopHad_massReg = -99.;
        TopHad_ptReg = -99.;
        TopHad_j3Pt = -99.;
        FatTopHad_mass = -99.;
        FatTopHad_pt = -99.;
        FatTopHad_eta = -99.;
        FatTopHad_phi = -99.;
        FatTopHad_massReg = -99.;
        FatTopHad_ptReg = -99.;
        FatTopHad_j3Pt = -99.;
        if (Vtype==2 || Vtype==3 || Vtype==4) {
            float mindRaJet = TMath::Min(2.0, (2.0 * 173. / H.pt));  // set initial dR to 2.0 or (2 * m_top / pT_top)
            int mindRaJet_id = -99;
            for (Int_t iaj = 0; iaj < naJets; iaj++) {
                if (aJet_pt[iaj] > 30. && fabs(aJet_eta[iaj]) < 4.5 && aJet_id[iaj] && aJet_puJetIdL[iaj]>0.) {
                    float dR = deltaR(aJet_eta[iaj], aJet_phi[iaj], H.eta, H.phi);
                    if (mindRaJet > dR) {
                        mindRaJet = dR;
                        mindRaJet_id = iaj;
                    }
                }
            }
            //const TLorentzVector& hJet_p4Norm_0 = makePtEtaPhiE(hJet_pt[0]   , hJet_pt[0], hJet_eta[0], hJet_phi[0], hJet_e[0]);
            //const TLorentzVector& hJet_p4Norm_1 = makePtEtaPhiE(hJet_pt[1]   , hJet_pt[1], hJet_eta[1], hJet_phi[1], hJet_e[1]);
            //const TLorentzVector& hJet_p4Reg_0  = makePtEtaPhiE(hJet_ptReg[0], hJet_pt[0], hJet_eta[0], hJet_phi[0], hJet_e[0]);
            //const TLorentzVector& hJet_p4Reg_1  = makePtEtaPhiE(hJet_ptReg[1], hJet_pt[1], hJet_eta[1], hJet_phi[1], hJet_e[1]);
            TLorentzVector Top_3rdjet;
            if (mindRaJet_id == -99)
                Top_3rdjet = p4Zero;
            Top_3rdjet.SetPtEtaPhiE(aJet_pt[mindRaJet_id], aJet_eta[mindRaJet_id], aJet_phi[mindRaJet_id], aJet_e[mindRaJet_id]);
            TLorentzVector Top = hJet_p4Norm_0 + hJet_p4Norm_1 + Top_3rdjet;
            TopHad_mass = Top.M();  // = H.mass if there is no 3rd jets
            TopHad_pt = Top.Pt();
            TopHad_eta = Top.Eta();
            TopHad_phi = Top.Phi();
            TLorentzVector TopReg = hJet_p4Reg_0 + hJet_p4Reg_1 + Top_3rdjet;  // = HmassReg if there is no 3rd jets
            TopHad_massReg = TopReg.M();
            TopHad_ptReg = TopReg.Pt();
            TopHad_j3Pt = Top_3rdjet.Pt();
            
            if (nfathFilterJets>0 && fathFilterJets_pt[0]>0. && fathFilterJets_pt[1]>0.){
                float mindRaJetFat = TMath::Min(2.0, (2.0 * 173. / FatH.pt));  // set initial dR to 2.0 or (2 * m_top / pT_top)
                int mindRaJetFat_id = -99.;
                for (Int_t iaj = 0; iaj < naJetsFat; iaj++) {
                    //if (aJetFat_pt[iaj] > 30. && fabs(aJetFat_eta[iaj]) < 4.5 && aJetFat_id[iaj] && aJetFat_puJetIdL[iaj]>0.) {
                    if (aJetFat_pt[iaj] > 30. && fabs(aJetFat_eta[iaj]) < 4.5) {
                        float dR = deltaR(aJetFat_eta[iaj], aJetFat_phi[iaj], FatH.eta, FatH.phi);
                        if (mindRaJetFat > dR) {
                            mindRaJetFat = dR;
                            mindRaJetFat_id = iaj;
                        }
                    }
                }
                const TLorentzVector& fathFilterJets_p4Norm_0 = makePtEtaPhiE(fathFilterJets_pt[0]   , fathFilterJets_pt[0], fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_e[0]);
                const TLorentzVector& fathFilterJets_p4Norm_1 = makePtEtaPhiE(fathFilterJets_pt[1]   , fathFilterJets_pt[1], fathFilterJets_eta[1], fathFilterJets_phi[1], fathFilterJets_e[1]);
                const TLorentzVector& fathFilterJets_p4Norm_2 = (nfathFilterJets > 2 && fathFilterJets_pt[2] > 0) ? 
                                                                makePtEtaPhiE(fathFilterJets_pt[2]   , fathFilterJets_pt[2], fathFilterJets_eta[2], fathFilterJets_phi[2], fathFilterJets_e[2]) : p4Zero;
                const TLorentzVector& fathFilterJets_p4Reg_0  = makePtEtaPhiE(fathFilterJets_ptReg[0], fathFilterJets_pt[0], fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_e[0]);
                const TLorentzVector& fathFilterJets_p4Reg_1  = makePtEtaPhiE(fathFilterJets_ptReg[1], fathFilterJets_pt[1], fathFilterJets_eta[1], fathFilterJets_phi[1], fathFilterJets_e[1]);
                const TLorentzVector& fathFilterJets_p4Reg_2  = (nfathFilterJets > 2 && fathFilterJets_pt[2] > 0) ? 
                                                                makePtEtaPhiE(fathFilterJets_ptReg[2], fathFilterJets_pt[2], fathFilterJets_eta[2], fathFilterJets_phi[2], fathFilterJets_e[2]) : p4Zero;
                TLorentzVector FatTop_3rdjet = fathFilterJets_p4Norm_2;
                if (fathFilterJets_pt[2] <= 15. && mindRaJetFat_id != -99 && aJetFat_pt[mindRaJetFat_id] > 0.)
                    FatTop_3rdjet.SetPtEtaPhiE(aJetFat_pt[mindRaJetFat_id], aJetFat_eta[mindRaJetFat_id], aJetFat_phi[mindRaJetFat_id], aJetFat_e[mindRaJetFat_id]);  // aJetFat_pt is > 20
                TLorentzVector FatTop_3rdjetReg = fathFilterJets_p4Reg_2;
                if (fathFilterJets_ptReg[2] <= 15. && mindRaJetFat_id != -99 && aJetFat_pt[mindRaJetFat_id] > 0.)
                    FatTop_3rdjetReg.SetPtEtaPhiE(aJetFat_pt[mindRaJetFat_id], aJetFat_eta[mindRaJetFat_id], aJetFat_phi[mindRaJetFat_id], aJetFat_e[mindRaJetFat_id]);  // aJetFat_pt is > 20
                TLorentzVector FatTop = fathFilterJets_p4Norm_0 + fathFilterJets_p4Norm_1 + FatTop_3rdjet;
                FatTopHad_mass = FatTop.M();
                FatTopHad_pt = FatTop.Pt();
                FatTopHad_eta = FatTop.Eta();
                FatTopHad_phi = FatTop.Phi();
                TLorentzVector FatTopReg = fathFilterJets_p4Reg_0 + fathFilterJets_p4Reg_1 + FatTop_3rdjetReg;
                FatTopHad_massReg = FatTopReg.M();
                FatTopHad_ptReg = FatTopReg.Pt();
                FatTopHad_j3Pt = FatTop_3rdjet.Pt();

                //if (FatTopHad_j3Pt > 10000.) std::cout << FatTopHad_j3Pt << " " << fathFilterJets_pt[0] << " " << fathFilterJets_pt[1] << " " << fathFilterJets_pt[2] << " " << aJetFat_pt[mindRaJetFat_id] << " " << hJet_pt[0] << " " << mindRaJetFat << std::endl;
            }
        }

        outTree->Fill();  // fill it!
    }  // end loop over TTree entries

    /// Get elapsed time
    sw.Stop();
    std::cout << "--- End of event loop: ";
    sw.Print();

    output->cd();
    outTree->Write();
    output->Close();
    input->Close();

    delete input;
    delete output;
    delete reader;
    delete readerFJ;

    for (formIt=inputFormulasReg0.begin(), formItEnd=inputFormulasReg0.end(); formIt!=formItEnd; formIt++)
        delete *formIt;
    for (formIt=inputFormulasReg1.begin(), formItEnd=inputFormulasReg1.end(); formIt!=formItEnd; formIt++)
        delete *formIt;
    for (formIt=inputFormulasFJReg0.begin(), formItEnd=inputFormulasFJReg0.end(); formIt!=formItEnd; formIt++)
        delete *formIt;
    for (formIt=inputFormulasFJReg1.begin(), formItEnd=inputFormulasFJReg1.end(); formIt!=formItEnd; formIt++)
        delete *formIt;
    for (formIt=inputFormulasFJReg2.begin(), formItEnd=inputFormulasFJReg2.end(); formIt!=formItEnd; formIt++)
        delete *formIt;
    delete ttf_lheweight;

#ifdef CSVSYST
    delete csvShape;
    delete csvShape_eff_b_up;
    delete csvShape_eff_b_down;
    delete csvShape_fake_b_up;
    delete csvShape_fake_b_down;
#endif

    std::cout << "==> GrowTree is done!" << std::endl << std::endl;
    return;
}
