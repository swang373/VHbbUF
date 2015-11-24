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
#include "Math/QuantFuncMathCore.h"

namespace math {
    typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > XYZTLorentzVector;
}

#define MAXJ 100
#define MAXL 110


// FIXME: change it to output the delta, not met+delta
double shift_met_by_Nvtx(double met, double metphi, int Nvtx, int EVENT_run)
{
    // from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/JetMETCorrections/Type1MET/python/pfMETsysShiftCorrections_cfi.py?revision=1.6&view=markup
    double mex = met * cos(metphi);
    double mey = met * sin(metphi);
    double px = 0.0, py = 0.0;
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
double shift_metphi_by_Nvtx(double met, double metphi, int Nvtx, int EVENT_run)
{
    // from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/JetMETCorrections/Type1MET/python/pfMETsysShiftCorrections_cfi.py?revision=1.6&view=markup
    double mex = met * cos(metphi);
    double mey = met * sin(metphi);
    double px = 0.0, py = 0.0;
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
    double phi1 = TMath::ATan2(mey, mex);
    double phi2 = TMath::ATan2(mey, mex) - 2.0 * M_PI;
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


TLorentzVector makePtEtaPhiM(double ptCor, double pt, double eta, double phi, double m, 
                             bool print=false)
{
    TLorentzVector j;
    if (print ) { 
      std::cout << "evaluating vector with ptCor, pt, eta, phi, m " << ptCor << " , " <<  pt << " , " << eta << " , " << phi << " , " << m << std::endl;
    }
    j.SetPtEtaPhiM(ptCor, eta, phi,  m );
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

    const TString indir   = "/afs/cern.ch/work/d/degrutto/public/MiniAOD/ZnnHbb_Phys14_PU20bx25/skimV11/";
    const TString outdir  = "/afs/cern.ch/work/d/degrutto/public/MiniAOD/ZnnHbb_Phys14_PU20bx25/skimV11/step3/";
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
    double hJet_pt[MAXJ], hJet_eta[MAXJ], hJet_phi[MAXJ], hJet_m[MAXJ], hJet_ptRaw[MAXJ], hJet_genPt[MAXJ];
    int hJCidx[2];

    inTree->SetBranchStatus("*", 1);

    inTree->SetBranchStatus("hJCidx",1);
    inTree->SetBranchStatus("Jet_*",1);
    inTree->SetBranchAddress("hJCidx", &hJCidx);
    inTree->SetBranchAddress("Jet_pt", &hJet_pt);

    inTree->SetBranchAddress("Jet_eta", &hJet_eta);
    inTree->SetBranchAddress("Jet_phi", &hJet_phi);
    inTree->SetBranchAddress("Jet_mass", &hJet_m);
    inTree->SetBranchAddress("Jet_rawPt", &hJet_ptRaw);

    inTree->SetBranchAddress("Jet_mcPt", &hJet_genPt);


    ///-- Make new branches ----------------------------------------------------
    int EVENT_run, EVENT_event;  // set these as TTree index?
    float lumi_ = lumi, efflumi, efflumi_old, 
        efflumi_UEPS_up, efflumi_UEPS_down;
    float hJet_ptReg[2];
    float HptNorm, HptGen, HptReg;
    float HmassNorm, HmassGen, HmassReg;

    outTree->Branch("EVENT_run", &EVENT_run, "EVENT_run/I");
    outTree->Branch("EVENT_event", &EVENT_event, "EVENT_event/I");
    outTree->Branch("lumi", &lumi_, "lumi/F");
    outTree->Branch("efflumi", &efflumi, "efflumi/F");
    outTree->Branch("efflumi_old", &efflumi_old, "efflumi_old/F");
    outTree->Branch("efflumi_UEPS_up", &efflumi_UEPS_up, "efflumi_UEPS_up/F");
    outTree->Branch("efflumi_UEPS_down", &efflumi_UEPS_down, "efflumi_UEPS_down/F");
    outTree->Branch("hJet_ptReg", &hJet_ptReg, "hJet_ptReg[2]/F");
    
    outTree->Branch("HptNorm", &HptNorm, "HptNorm/F");
    outTree->Branch("HptGen", &HptGen, "HptGen/F");
    outTree->Branch("HptReg", &HptReg, "HptReg/F");
    outTree->Branch("HmassNorm", &HmassNorm, "HmassNorm/F");
    outTree->Branch("HmassGen", &HmassGen, "HmassGen/F");
    outTree->Branch("HmassReg", &HmassReg, "HmassReg/F");

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

    // regression stuff here

    
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
    //    assert(idx_rawpt!=-1 && idx_pt!=-1 && idx_et!=-1 && idx_mt!=-1);
    assert(idx_rawpt!=-1 && idx_pt!=-1 );

    /// Setup TMVA regression inputs
    const std::vector < std::string > & inputExpressionsReg0 = GetInputExpressionsReg0();
    const std::vector < std::string > & inputExpressionsReg1 = GetInputExpressionsReg1();
    assert(inputExpressionsReg0.size() == nvars);
    assert(inputExpressionsReg1.size() == nvars);

    /// Load TMVA weights
    TString weightdir  = "weights/";
    TString weightfile = weightdir + "TMVARegression_" + regMethod + ".testweights.xml";
    reader->BookMVA(regMethod + " method", weightfile);
    
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
	//        efflumi_UEPS_up   = efflumi * hcount->GetBinContent(2) / hcount->GetBinContent(3);
        //efflumi_UEPS_down = efflumi * hcount->GetBinContent(2) / hcount->GetBinContent(4);
#endif
    
        bool verbose = false;
 
	for (Int_t ihj = 0; ihj < 2; ihj++) {

   
            /// Evaluate TMVA regression output
            for (UInt_t iexpr = 0; iexpr < nvars; iexpr++) {
                if (ihj==0) {
                    readerVars[iexpr] = inputFormulasReg0.at(iexpr)->EvalInstance();

                } else if (ihj==1) {
                    readerVars[iexpr] = inputFormulasReg1.at(iexpr)->EvalInstance();
                }
            }

	    hJet_ptReg[ihj]               = (reader->EvaluateRegression(regMethod + " method"))[0];
            if (verbose)  std::cout << readerVars[idx_pt] << " " << readerVars[idx_rawpt] <<  " " << hJet_pt[ihj] << " " << hJet_ptReg[ihj] << " " << hJet_genPt[ihj] << std::endl;
        const TLorentzVector p4Zero                     = TLorentzVector(0., 0., 0., 0.);
	//	int idx =  hJCidx[0] ;
	//	std::cout << "the regressed pt for jet 0 is " << hJet_ptReg[0] << "; the hJCidx is " << hJCidx[0] << ", hence the origianl pt is " <<  hJet_pt[idx] << std::endl;

	
       
        const TLorentzVector& hJet_p4Norm_0             = makePtEtaPhiM(hJet_pt[hJCidx[0]]                , hJet_pt[hJCidx[0]], hJet_eta[hJCidx[0]], hJet_phi[hJCidx[0]], hJet_m[hJCidx[0]]);
        const TLorentzVector& hJet_p4Norm_1             = makePtEtaPhiM(hJet_pt[hJCidx[1]]                , hJet_pt[hJCidx[1]], hJet_eta[hJCidx[1]], hJet_phi[hJCidx[1]], hJet_m[hJCidx[1]]);
        const TLorentzVector& hJet_p4Gen_0              = hJet_genPt[hJCidx[0]] > 0 ? 
                                                          makePtEtaPhiM(hJet_genPt[hJCidx[0]]             , hJet_pt[hJCidx[0]], hJet_eta[hJCidx[0]], hJet_phi[hJCidx[0]], hJet_m[hJCidx[0]]) : p4Zero;
        const TLorentzVector& hJet_p4Gen_1              = hJet_genPt[hJCidx[1]] > 0 ? 
                                                          makePtEtaPhiM(hJet_genPt[hJCidx[1]]             , hJet_pt[hJCidx[1]], hJet_eta[hJCidx[1]], hJet_phi[hJCidx[1]], hJet_m[hJCidx[1]]) : p4Zero;
        const TLorentzVector& hJet_p4Reg_0              = makePtEtaPhiM(hJet_ptReg[0]             , hJet_pt[hJCidx[0]], hJet_eta[hJCidx[0]], hJet_phi[hJCidx[0]], hJet_m[hJCidx[0]]);
        const TLorentzVector& hJet_p4Reg_1              = makePtEtaPhiM(hJet_ptReg[1]             , hJet_pt[hJCidx[1]], hJet_eta[hJCidx[1]], hJet_phi[hJCidx[1]], hJet_m[hJCidx[1]]);
        HptNorm             = (hJet_p4Norm_0             + hJet_p4Norm_1            ).Pt();
        HptGen              = (hJet_p4Gen_0              + hJet_p4Gen_1             ).Pt();
        HptReg              = (hJet_p4Reg_0              + hJet_p4Reg_1             ).Pt();
        HmassNorm             = (hJet_p4Norm_0             + hJet_p4Norm_1            ).M();
        HmassGen              = (hJet_p4Gen_0              + hJet_p4Gen_1             ).M();
        HmassReg              = (hJet_p4Reg_0              + hJet_p4Reg_1             ).M();
	//        std::cout << "HmassReg is " << HmassReg << std::endl; 
	
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

    std::cout << "==> GrowTree is done!" << std::endl << std::endl;
    return;
}
