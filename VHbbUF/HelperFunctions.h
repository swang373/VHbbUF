#ifndef ROOT_TMath
#include "TMath.h"
#endif

#ifndef ROOT_TRandom
#include "TRandom.h"
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


inline double quadsum(double a, double b)
{
    return TMath::Sqrt(a*a + b*b);
}

//______________________________________________________________________________
inline double deltaPhi(double phi1, double phi2)
{
  double result = phi1 - phi2;
  while (result > M_PI) result -= 2.0*M_PI;
  while (result <= -M_PI) result += 2.0*M_PI;
  return result;
}

inline float deltaPhi(float phi1, float phi2)
{
  float result = phi1 - phi2;
  while (result > float(M_PI)) result -= float(2.0*M_PI);
  while (result <= -float(M_PI)) result += float(2.0*M_PI);
  return result;
}

inline double deltaR(double eta1, double phi1, double eta2, double phi2)
{
  double deta = eta1 - eta2;
  double dphi = deltaPhi(phi1, phi2);
  return TMath::Sqrt(deta*deta + dphi*dphi);
}

inline float deltaR(float eta1, float phi1, float eta2, float phi2)
{
  float deta = eta1 - eta2;
  float dphi = deltaPhi(phi1, phi2);
  return TMath::Sqrt(deta*deta + dphi*dphi);
}

//______________________________________________________________________________
inline float deltaPhiMETjets(float metphi, float jetphi, float jetpt, float jeteta, float minpt=30, float maxeta=2.5)
{
    if(jetpt <= minpt || TMath::Abs(jeteta) >= maxeta)
        return 999.;
    else
        return deltaPhi(metphi, jetphi);
}

inline float projectionMETOntoJet(float met, float metphi, float jet, float jetphi, bool onlyPositive=true, float threshold=M_PI/4.0)
{
    float deltaphi = deltaPhi(metphi, jetphi);
    float met_dot_jet = met * jet * TMath::Cos(deltaphi);
    float jetsq = jet * jet;
    float projection = met_dot_jet / jetsq * jet;
    
    if (onlyPositive && TMath::Abs(deltaphi) >= threshold)
        return 0.0;
    else
        return projection;
}

inline bool closestToMET(float metphi, float jetphi, float jet1phi, float jet2phi)
{
    float deltaphi1 = deltaPhi(metphi, jet1phi);
    float deltaphi2 = deltaPhi(metphi, jet2phi);
    float mindeltaphi = TMath::Min(deltaphi1, deltaphi2);
    return (jetphi <= mindeltaphi+1e-6);
}

//______________________________________________________________________________
#include "TLorentzVector.h"
inline double evalEt(double pt, double eta, double phi, double e)
{
    TLorentzVector j;
    j.SetPtEtaPhiE(pt, eta, phi, e);
    return j.Et();
}

inline double evalMt(double pt, double eta, double phi, double e)
{
    TLorentzVector j;
    j.SetPtEtaPhiE(pt, eta, phi, e);
    return j.Mt();
}

inline double evalDijetPt(double ptCor0, double pt0, double eta0, double phi0, double e0, 
                          double ptCor1, double pt1, double eta1, double phi1, double e1,
                          bool rescaleEnergy=true)
{
    TLorentzVector j0, j1;
    j0.SetPtEtaPhiE(ptCor0, eta0, phi0, (rescaleEnergy ? e0 * ptCor0 / pt0 : e0));
    j1.SetPtEtaPhiE(ptCor1, eta1, phi1, (rescaleEnergy ? e1 * ptCor1 / pt1 : e1));
    return (j0+j1).Pt();
}

inline double evalDijetMass(double ptCor0, double pt0, double eta0, double phi0, double e0, 
                            double ptCor1, double pt1, double eta1, double phi1, double e1,
                            bool rescaleEnergy=true)
{
    TLorentzVector j0, j1;
    j0.SetPtEtaPhiE(ptCor0, eta0, phi0, (rescaleEnergy ? e0 * ptCor0 / pt0 : e0));
    j1.SetPtEtaPhiE(ptCor1, eta1, phi1, (rescaleEnergy ? e1 * ptCor1 / pt1 : e1));
    return (j0+j1).M();
}

inline double evalDijetMassDrop(double pt0, double eta0, double phi0, double e0, 
                                double pt1, double eta1, double phi1, double e1)
{
    TLorentzVector j0, j1;
    j0.SetPtEtaPhiE(pt0, eta0, phi0, e0);
    j1.SetPtEtaPhiE(pt1, eta1, phi1, e1);
    return TMath::Min(j0.M(), j1.M())/(j0+j1).M();
}

//______________________________________________________________________________
//inline float patchMETtype1corr(float met, float metphi, float met_diff, float metphi_diff)
//{
//    TVector2 j0, j1;
//    j0.SetMagPhi(met, metphi);
//    j1.SetMagPhi(met_diff, metphi_diff);
//    return (j0+j1).Mod();
//}

//inline float patchMETtype1corr_phi(float met, float metphi, float met_diff, float metphi_diff)
//{
//    TVector2 j0, j1;
//    j0.SetMagPhi(met, metphi);
//    j1.SetMagPhi(met_diff, metphi_diff);
//    return (j0+j1).Phi();
//}


//______________________________________________________________________________
inline double evalCosThetaHbb(double ptCor0, double pt0, double eta0, double phi0, double e0, 
                              double ptCor1, double pt1, double eta1, double phi1, double e1,
                              bool rescaleEnergy=true)
{
    TLorentzVector j0, j1;
    j0.SetPtEtaPhiE(ptCor0, eta0, phi0, (rescaleEnergy ? e0 * ptCor0 / pt0 : e0));
    j1.SetPtEtaPhiE(ptCor1, eta1, phi1, (rescaleEnergy ? e1 * ptCor1 / pt1 : e1));
    TLorentzVector jsum = j0+j1;
    TVector3 boostsum = jsum.BoostVector();
    j0.Boost(-boostsum);
    j1.Boost(-boostsum);
    TVector3 boost;
    if(int(pt0) % 2 == 0)
        boost = j0.BoostVector();
    else
        boost = j1.BoostVector();

    double cosTheta = boost.Dot(boostsum) / (boost.Mag() * boostsum.Mag());
    return cosTheta;
}

inline double evalHMETMt(double pt0, double phi0, double pt1, double phi1)
{
    return TMath::Sqrt(2.0 * pt0 * pt1 * (1.0 - TMath::Cos(deltaPhi(phi0, phi1)) ));
}

inline double evalHMETMassiveMt(double m0, double pt0, double phi0, double m1, double pt1, double phi1)
{
    return TMath::Sqrt(m0*m0 + m1*m1 + 2.0 * (TMath::Sqrt((m0*m0+pt0*pt0)*(m1*m1+pt1*pt1)) - pt0*pt1*TMath::Cos(deltaPhi(phi0, phi1)) ));
}

//______________________________________________________________________________
#define JERSUMMER11

inline float smear_pt_res(float pt, float genpt, float eta)
{
    eta = fabs(eta);
    if (genpt>15. && (fabs(pt - genpt) / pt)<0.5) {  // limit the effect to the core
        double res    = 1.0;
        double resErr = 0.0;
#ifdef JERSUMMER11
        // from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
        if (eta <= 0.5) {
            res    = 1.052;
            resErr = quadsum(0.012, 0.062);
        } else if (0.5 < eta && eta <= 1.1) {
            res    = 1.057;
            resErr = quadsum(0.012, 0.062);
        } else if (1.1 < eta && eta <= 1.7) {
            res    = 1.096;
            resErr = quadsum(0.017, 0.063);
        } else if (1.7 < eta && eta <= 2.3) {
            res    = 1.134;
            resErr = quadsum(0.035, 0.087);
        } else {
            res    = 1.288;
            resErr = quadsum(0.127, 0.155);
        }
#else
        // from VHbb analysis
        if (eta <= 1.1) {
            res    = 1.05;
            resErr = 0.05;
        } else if (1.1 < eta && eta <= 2.5) {
            res    = 1.10;
            resErr = 0.10;
        } else {
            res    = 1.30;
            resErr = 0.20;
        }
#endif
        float deltapt = (pt - genpt) * res;
        return TMath::Max(float(0.), genpt + deltapt);
    }
    return pt;
}

inline float smear_pt_resErr(float pt, float genpt, float eta, bool up)
{
    eta = fabs(eta);
    if (genpt>15. && (fabs(pt - genpt) / pt)<0.5) {  // limit the effect to the core
        double res    = 1.0;
        double resErr = 0.0;
#ifdef JERSUMMER11
        // from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
        if (eta <= 0.5) {
            res    = 1.052;
            resErr = quadsum(0.012, 0.062);
        } else if (0.5 < eta && eta <= 1.1) {
            res    = 1.057;
            resErr = quadsum(0.012, 0.062);
        } else if (1.1 < eta && eta <= 1.7) {
            res    = 1.096;
            resErr = quadsum(0.017, 0.063);
        } else if (1.7 < eta && eta <= 2.3) {
            res    = 1.134;
            resErr = quadsum(0.035, 0.087);
        } else {
            res    = 1.288;
            resErr = quadsum(0.127, 0.155);
        }
#else
        // from VHbb analysis
        if (eta <= 1.1) {
            res    = 1.05;
            resErr = 0.05;
        } else if (1.1 < eta && eta <= 2.5) {
            res    = 1.10;
            resErr = 0.10;
        } else {
            res    = 1.30;
            resErr = 0.20;
        }
#endif
        float deltapt = (pt - genpt) * resErr;
        if (up)
            return TMath::Max(float(0.), pt + deltapt);
        else
            return TMath::Max(float(0.), pt - deltapt);
    }
    return pt;
}

inline float smear_pt_scaleErr(float pt, float scaleErr, bool up)
{
    if (up)
        return pt * (1.0 + scaleErr);
    else
        return pt * (1.0 - scaleErr);
}

float apply_pt_slope(float pt, float slope, float start, float initweight=1.0, float deltaweight=0.5)
{
    if (slope<0.)
        return TMath::Max(initweight-deltaweight, (initweight + slope * (TMath::Max(pt,start) - start)) );
    else
        return TMath::Min(initweight+deltaweight, (initweight + slope * (TMath::Max(pt,start) - start)) );
}

//______________________________________________________________________________
inline float smear_gen_to_reco(float pt, float genpt, float eta)
{
    eta = fabs(eta);
    /// Values from http://arxiv.org/abs/1301.5023
    float N    = -5.090;
    float S    =  0.512;
    float C    =  0.033;
    float M    =  0.325;
    float sgnN = (N>0) ? 1.0 : -1.0;

    /// Values from AN_2010/371 Table 4
    if (eta <= 0.5) {
        N    =  3.96859;
        S    =  0.18348;
        C    =  0.00000;
        M    =  0.62627;
        sgnN = (N>0) ? 1.0 : -1.0;
    } else if (0.5 < eta && eta <= 1.0) {
        N    =  3.55226;
        S    =  0.24026;
        C    =  0.00000;
        M    =  0.52571;
        sgnN = (N>0) ? 1.0 : -1.0;
    } else if (1.0 < eta && eta <= 1.5) {
        N    =  4.54826;
        S    =  0.22652;
        C    =  0.00000;
        M    =  0.58963;
        sgnN = (N>0) ? 1.0 : -1.0;
    } else if (1.5 < eta && eta <= 2.0) {
        N    =  4.62622;
        S    =  0.23664;
        C    =  0.00000;
        M    =  0.48738;
        sgnN = (N>0) ? 1.0 : -1.0;
    } else if (2.0 < eta && eta <= 2.5) {
        N    =  2.53324;
        S    =  0.34306;
        C    =  0.00000;
        M    =  0.28662;
        sgnN = (N>0) ? 1.0 : -1.0;
    } else if (2.5 < eta && eta <= 3.0) {
        N    = -3.33814;
        S    =  0.73360;
        C    =  0.00000;
        M    =  0.08264;
        sgnN = (N>0) ? 1.0 : -1.0;
    } else {
        N    =  2.95397;
        S    =  0.11619;
        C    =  0.00000;
        M    =  0.96086;
        sgnN = (N>0) ? 1.0 : -1.0;
    }
    float res = genpt * TMath::Sqrt(sgnN * N*N/genpt/genpt + S*S*TMath::Power(genpt,M)/genpt + C*C);
    float deltapt = gRandom->Gaus(1.0, res);
    
    return TMath::Max(float(0.), genpt - deltapt);
}

inline float smear_gen_by_double_gaussian(float genpt, float p0, float p1, float p2, float p3, float p4)
{
    bool gaus1_or_gaus2     = (gRandom->Uniform(0,1)) < p4;  // 0: use narrow gaus1, 1: use wide gaus2
    float jetres            = gaus1_or_gaus2 ? 
                              gRandom->Gaus(p2, p1 * p3) :  // wide gaus2
                              gRandom->Gaus(p0, p1);  // narrow gaus1
    return TMath::Max(float(0.), genpt * jetres);
}

//______________________________________________________________________________
//For 2012A HLT_DiCentralPFJet30_PFMHT80, valid for pfMET > 100 GeV:
inline float scaleDiJet30MHT80_2012A(float x)
{
  if(x<100) return 0;
  return   (1e0 - exp(-0.04197*(x-75.73))) * 0.9721 ;
}
//For 2012B HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned, valid for pfMET > 100 GeV:
inline float scaleSumpT100MET100_2012B(float x)
{
  if(x<100) return 0;
  return   (1e0 - exp(-0.06704*(x-96.84))) * 0.9199 ;
}
//For 2012A+B HLT_PFMET150, valid for pfMET > 150 GeV:
inline float scalePFMET150_2012AB(float x)
{
  if(x<150) return 0;
  return  (1e0 - exp(-0.07135*(x-147.4))) * 0.9707;
}

//For 2012A HLT_PFMET150 OR HLT_DiCentralPFJet30_PFMHT80, valid for pfMET > 100 GeV:
inline float scalePFMET150orDiJetMET_2012A(float x)
{
  if(x<100) return 0;
  return  (1e0 - exp(-0.0412*(x-75.52))) * 0.9772;
}
//For 2012B HLT_PFMET150 OR HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned, valid for pfMET > 100 GeV:
inline float scalePFMET150orDiJetMET_2012B(float x)
{
  if(x<100) return 0;
  return  (1e0 - exp(-0.05482*(x-95.59))) * 0.9702;
}
//For 2012C HLT_PFMET150 OR HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned, valid for pfMET > 100 GeV:
inline float scalePFMET150orDiJetMET_2012C(float x)
{
  if(x<100) return 0;
  return  (1e0 - exp(-0.05627*(x-95.15))) * 0.9659;
}

//For 2012D HLT_PFMET150 OR HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned, valid for pfMET > 100 GeV:
inline float scalePFMET150orDiJetMET_2012D(float x)
{
  return  scalePFMET150orDiJetMET_2012C(x);
}

//______________________________________________________________________________
inline float triggercorrMET(float met)
{
    if (met < 110.)  return 0.940;
    if (met < 120.)  return 0.977;
    if (met < 130.)  return 0.984;
    if (met < 140.)  return 0.988;
    if (met < 150.)  return 1.000;
    return 1.;
}

inline float triggercorrMET_up(float met)
{
    if (met < 110.)  return 0.98;
    return 1.;
}

inline float triggercorrMET_down(float met)
{
    if (met < 110.) return 0.886;
    if (met < 120.) return 0.905;
    if (met < 130.) return 0.953;
    if (met < 140.) return 0.967;
    if (met < 150.) return 0.970;
    if (met < 200.) return 0.980;
    if (met < 250.) return 0.990;
    return 1.;
}

inline float triggercorrCSV(float csv)
{
    // naJets_Znn <= 1
    if (csv>0.244) return (1.375 + (-0.375) * csv);
    return (1.375 + (-0.375) * 0.244);
    
    // naJets_Znn <= 0
    //if (csv>0.244) return (1.389 + (-0.417) * csv);
    //return (1.389 + (-0.417) * 0.244);
}

inline float triggercorrCSV_up(float csv)
{
    if (csv>0.244) return (1.750 + (-0.750) * csv);
    return (1.750 + (-0.750) * 0.244);
}

inline float triggercorrCSV_down(float csv)
{
    return 1.;
}

//______________________________________________________________________________
inline float triggerweight2012AB(float pfmet)
{
    // lumi recorded by HLT_PFMET150_*
    // 2012A:  703.529924 pb
    // 2012B: 4388.772348 pb
    float weightPFMET150 = scalePFMET150_2012AB(pfmet);
    float weightDiJet30MHT80 = scaleDiJet30MHT80_2012A(pfmet);
    float weightSumpT100MET100 = scaleSumpT100MET100_2012B(pfmet);
    return TMath::Max(0., (703.5*(weightPFMET150+weightDiJet30MHT80 - (weightPFMET150*weightDiJet30MHT80)) + 4388.8*(weightPFMET150+weightSumpT100MET100 - (weightPFMET150*weightSumpT100MET100)) )/5092.3);
}

inline float triggerweight2012ABC(float pfmet)
{
    // lumi recorded by HLT_PFMET150_* 
    // 2012A        :   809.379300 pb
    // 2012A-recover:    82.135672 pb
    // 2012B        :  4403.763061 pb
    // 2012C-rereco :   495.002899 pb
    // 2012C        :  6311.311736 pb
    // TOTAL        : 12101.592668 pb
    float weightPFMET150orDiJetMET_2012A = scalePFMET150orDiJetMET_2012A(pfmet);
    float weightPFMET150orDiJetMET_2012B = scalePFMET150orDiJetMET_2012B(pfmet);
    float weightPFMET150orDiJetMET_2012C = scalePFMET150orDiJetMET_2012C(pfmet);
    return TMath::Max(0., (891.514972*weightPFMET150orDiJetMET_2012A + 4403.763061*weightPFMET150orDiJetMET_2012B + 6806.314635*weightPFMET150orDiJetMET_2012C )/12101.592668);
    
    // Prompt only
    //return TMath::Max(0., (698.190585*weightPFMET150orDiJetMET_2012A + 4427.198353*weightPFMET150orDiJetMET_2012B + 6806.235235*weightPFMET150orDiJetMET_2012C)/11931.624173);
}

inline float triggerweight2012ABC_up(float pfmet)
{
    float weight = triggerweight2012ABC(pfmet);
    if (pfmet < 100)
        return weight;
    else if (pfmet < 150) 
        weight *= (1.0 + 0.01);
    else if (pfmet < 200) 
        weight *= (1.0 - 0.02);
    return TMath::Min(1.0, weight * (1.0 + 0.03));
}

inline float triggerweight2012ABC_down(float pfmet)
{
    float weight = triggerweight2012ABC(pfmet);
    if (pfmet < 100)
        return weight;
    else if (pfmet < 150) 
        weight *= (1.0 + 0.01);
    else if (pfmet < 200) 
        weight *= (1.0 - 0.02);
    return TMath::Min(1.0, weight * (1.0 - 0.03));
}

inline float triggerweight2012ABCD(float pfmet)
{
    double bincontents[40] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.323, 0.331, 0.305, 0.415, 0.535, 0.640, 0.732, 0.795, 0.870, 0.907, 0.942, 0.955, 0.963, 0.969, 0.989, 0.989, 0.991, 0.993, 0.997, 0.997, 0.997, 1, 1, 1, 1, 1};  // already factored in triggercorrMET()
    if (pfmet>200.) return 1.0;
    int binid = int(pfmet/5.);
    return bincontents[binid];
}

inline float triggerweight2012ABCD_up(float pfmet)
{
    return (triggerweight2012ABCD(pfmet) * triggercorrMET_up(pfmet) / triggercorrMET(pfmet));
}

inline float triggerweight2012ABCD_down(float pfmet)
{
    return (triggerweight2012ABCD(pfmet) * triggercorrMET_down(pfmet) / triggercorrMET(pfmet));
}

//______________________________________________________________________________
inline float triggercorr2012ABCD(bool mettriggerbit, bool metcsvtriggerbit, float pfmet, float maxcsv)
{
    float corr = float(mettriggerbit || metcsvtriggerbit) * triggercorrMET(pfmet);
    if (!mettriggerbit && metcsvtriggerbit)  // if passed only MET+CSV trigger
        corr *= triggercorrCSV(maxcsv);
    return corr;
}

inline float triggercorr2012ABCD_MET_up(bool mettriggerbit, bool metcsvtriggerbit, float pfmet, float maxcsv)
{
    float corr = float(mettriggerbit || metcsvtriggerbit) * triggercorrMET_up(pfmet);
    if (!mettriggerbit && metcsvtriggerbit)  // if passed only MET+CSV trigger
        corr *= triggercorrCSV(maxcsv);
    return corr;
}

inline float triggercorr2012ABCD_MET_down(bool mettriggerbit, bool metcsvtriggerbit, float pfmet, float maxcsv)
{
    float corr = float(mettriggerbit || metcsvtriggerbit) * triggercorrMET_down(pfmet);
    if (!mettriggerbit && metcsvtriggerbit)  // if passed only MET+CSV trigger
        corr *= triggercorrCSV(maxcsv);
    return corr;
}

inline float triggercorr2012ABCD_CSV_up(bool mettriggerbit, bool metcsvtriggerbit, float pfmet, float maxcsv)
{
    float corr = float(mettriggerbit || metcsvtriggerbit) * triggercorrMET(pfmet);
    if (!mettriggerbit && metcsvtriggerbit)  // if passed only MET+CSV trigger
        corr *= triggercorrCSV_up(maxcsv);
    return corr;
}

inline float triggercorr2012ABCD_CSV_down(bool mettriggerbit, bool metcsvtriggerbit, float pfmet, float maxcsv)
{
    float corr = float(mettriggerbit || metcsvtriggerbit) * triggercorrMET(pfmet);
    if (!mettriggerbit && metcsvtriggerbit)  // if passed only MET+CSV trigger
        corr *= triggercorrCSV_down(maxcsv);
    return corr;
}
