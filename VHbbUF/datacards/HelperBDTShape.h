#include <vector>
#include <iostream>
#include "TH1.h"
#include "TMath.h"

//#define REBINLASTTHREE

#define REBINFIRSTTHREE


///_____________________________________________________________________________
/// Rebinner

/// partnbins       : The number of bins in the returned TH1F. To make m=3 
///                   partitions of n=80 bins, input 808080. The max nbins in 
///                   each partition is 99; the max npartitions is 9 (due to 
///                   limit of long long). The number of rightmost partition 
///                   is indicated by the leftmost 2 digits.
/// errorffirst     : thresholds of the error fractions. Both the leftmost and 
/// errorflast        rightmost bins are grouped until the error fractions are 
///                   lower than the threshold. If it's outside of [0,1], 
///                   call TH1F::Rebin().
/// xlow, xup       : the boundaries of the returned TH1F. If they are
///                   the same, the boundaries of the input TH1F are used.

void findFirstN(UInt_t& firstN, const UInt_t begin, const TH1F* bkg, const double errorf) {
    assert(bkg != 0 && bkg->Integral() > 0);
    double firstbincontent  = bkg->GetBinContent(begin);
    double firstbinerror    = bkg->GetBinError(begin);
    double firstbinerror2   = firstbinerror * firstbinerror;
    double firstbinerrorf   = (firstbincontent > 1e0) ? (TMath::Sqrt(firstbinerror2) / firstbincontent) : 1.0;
    // Find firstN such that its error fraction becomes less than the threshold
    while ((firstbinerrorf > errorf) && (firstN <= (UInt_t) bkg->GetNbinsX())) {
        //hdummy->Info("Rebin", "firstN=%2i int. binc=%f int. bine=%f frac=%f", firstN, firstbincontent, sqrt(firstbinerror2), firstbinerrorf);
        //std::cout << firstN << " " << firstbincontent << " " << TMath::Sqrt(firstbinerror2) << " " << firstbinerrorf << std::endl;
        
        const UInt_t ibin = begin + firstN;
        const Double_t bincontent    = bkg->GetBinContent(ibin);
        const Double_t binerror      = bkg->GetBinError(ibin);
        firstbincontent             += bincontent;
        firstbinerror2              += (binerror * binerror);
        firstbinerrorf               = (firstbincontent > 1e0) ? (TMath::Sqrt(firstbinerror2) / firstbincontent) : 1.0;
        firstN ++;
    }
}

void findLastN(UInt_t& lastN, const UInt_t end, const TH1F* bkg, const double errorf) {
    assert(bkg != 0 && bkg->Integral() > 0);
    double lastbincontent   = bkg->GetBinContent(end-1);
    double lastbinerror     = bkg->GetBinError(end-1);
    double lastbinerror2    = lastbinerror * lastbinerror;
    double lastbinerrorf    = (lastbincontent > 1e0) ? (TMath::Sqrt(lastbinerror2) / lastbincontent) : 1.0;
    // Find lastN such that its error fraction becomes less than the threshold
    while ((lastbinerrorf > errorf) && (lastN <= (UInt_t) bkg->GetNbinsX())) {
        //hdummy->Info("Rebin", "lastN=%2i int. binc=%f int. bine=%f frac=%f", lastN, lastbincontent, sqrt(lastbinerror2), lastbinerrorf);
                
        const UInt_t ibin = end-1 - lastN;
        const Double_t bincontent    = bkg->GetBinContent(ibin);
        const Double_t binerror      = bkg->GetBinError(ibin);
        lastbincontent              += bincontent;
        lastbinerror2               += (binerror * binerror);
        lastbinerrorf                = (lastbincontent > 1e0) ? (TMath::Sqrt(lastbinerror2) / lastbincontent) : 1.0;
        lastN ++;
    }
}

#ifdef REBINLASTTHREE
void findLastN2(UInt_t& lastN, const UInt_t end, const double bkgrej, const double bkgint, const TH1F* bkg, const double errorf) {
    assert(bkg != 0 && bkg->Integral() > 0);
    assert(bkgrej > 0.0 && bkgrej < 1.0);
    double lastbincontent   = bkg->GetBinContent(end-1);
    double lastbinerror     = bkg->GetBinError(end-1);
    double lastbinerror2    = lastbinerror * lastbinerror;
    double lastbinerrorf    = (lastbincontent > 1e0) ? (TMath::Sqrt(lastbinerror2) / lastbincontent) : 1.0;
    // Find lastN such that its error fraction becomes less than the threshold
    //while ((lastbinerrorf > errorf || sig->Integral(end - lastN, end-1)/sigint < sigeff) /*&& (lastN <= oldnbins)*/) {
    while ((lastbinerrorf > errorf || (1.0 - bkg->Integral(end-1 - lastN+1, end-1)/bkgint) > bkgrej) && (lastN <= (UInt_t) bkg->GetNbinsX())) {
        //hdummy->Info("Rebin", "lastN=%2i int. binc=%f int. bine=%f frac=%f", lastN, lastbincontent, sqrt(lastbinerror2), lastbinerrorf);
                
        const UInt_t ibin = end-1 - lastN;
        const Double_t bincontent    = bkg->GetBinContent(ibin);
        const Double_t binerror      = bkg->GetBinError(ibin);
        lastbincontent              += bincontent;
        lastbinerror2               += (binerror * binerror);
        lastbinerrorf                = (lastbincontent > 1e0) ? (TMath::Sqrt(lastbinerror2) / lastbincontent) : 1.0;
        lastN ++;
    }
}
#endif

class Rebinner {
public:
    Rebinner(long long nb, double eff, double efl, double xl, double xu)
      : partnbins(nb),
        errorffirst(eff),
        errorflast(efl),
        xlow(xl),
        xup(xu),
        recreate(true),
        totalnbins(0),
        sig(0),
        bkg(0),
        nbinsvec(),
        lowedges()
    {
        /// Parse number of partitions and number of bins
        npartitions = 1;
        while (true) {
            nb = nb / 100;
            if (nb == 0)
                break;
            npartitions++;
        }
        
        for (int ipart=0; ipart < npartitions; ipart++) {
            int npb = partnbins % (long long)(pow(100, ipart+1)) / (long long)(pow(100, ipart));
            nbinsvec.push_back(npb);
            totalnbins += npb;
        }
        if (npartitions == 1)
            assert(totalnbins == partnbins);
        
        //lowedges.resize(totalnbins+1);
    }
    
    Rebinner(long long nb)
    {
        Rebinner(nb, 0.25, 0.25, 0., 0.);
    }
    
    ~Rebinner() {}
    
    void set_signal_backgr(const TH1F * s, const TH1F * b) {
        sig = s;
        bkg = b;
    }

    TH1F * rebin(const TH1F * h, long long nb, const char * newname ="")
    {
        TH1F * hdummy = (TH1F *) h->Clone("dummy");
        hdummy->Sumw2();
        const int    in_nbins = h->GetNbinsX();
        const TAxis* in_xaxis = h->GetXaxis();
#if !defined(REBINFIRSTTHREE) && !defined(REBINLASTTHREE)
        const unsigned int minoldnbins = 3;
#else
        const unsigned int minoldnbins = 5;
#endif
        if (in_nbins == nb) {  // just return a clone with newname
            //hdummy->Warning("Rebin", "nbins is the same. Return a clone. h=%s", h->GetName());
            hdummy->SetName(newname);
            return hdummy;
        }
        if (in_xaxis->IsVariableBinSize()) {
            hdummy->Warning("Rebin", "h has variable bin sizes. Return a clone. h=%s", h->GetName());
            hdummy->SetName(newname);
            return hdummy;
        }
        if (errorffirst < 0.0 || errorffirst > 1.0 || errorflast < 0.0 || errorflast > 1.0) {  // do TH1F::Rebin instead
            hdummy->Warning("Rebin", "errorffirst or errorflast is not within [0,1]. Call TH1F::Rebin(). h=%s", h->GetName());
            hdummy->SetName(newname);
            hdummy->Rebin(float(in_nbins)/float(nb));
            return hdummy;
        }
        
        /// From here on, don't allow underflow or overflow
        if (h->GetBinContent(0) > 0. || h->GetBinContent(in_nbins+1) > 0.) {
            hdummy->Error("Rebin", "Underflow and/or overflow is not empty! h=%s", h->GetName());
            return hdummy;
        }
        if (in_nbins % npartitions != 0) {
            hdummy->Error("Rebin", "in_nbins is not an exact divider of npartitions. in_nbins=%i, npartitions=%i", in_nbins, npartitions);
            return hdummy;
        }
        if (in_nbins / npartitions < int(minoldnbins)) {
            hdummy->Error("Rebin", "in_nbins / npartitions is less than minimum. in_nbins=%i, npartitions=%i", in_nbins, npartitions);
            return hdummy;
        }
        
        if (!recreate) {
            for (int ipart=0; ipart < npartitions; ipart++) {
                UInt_t newnbins = nbinsvec.at(ipart);
                UInt_t oldnbins = in_nbins/npartitions;
                if (newnbins < minoldnbins || newnbins > oldnbins) {
                    hdummy->Error("Rebin", "newnbins must be within [%i,%i]. newnbins=%i", minoldnbins, oldnbins, newnbins);
                    return hdummy;
                }
            }
            /// Check lowedges
            for (int iedge=0; iedge < totalnbins+1; iedge++) {
                bool found = false;
                double edge1 = lowedges.at(iedge);
                for (int ibin=0; ibin < in_nbins+2; ibin++) {
                    double edge2 = h->GetBinLowEdge(ibin+1);
                    if (TMath::AreEqualRel(edge1, edge2, 1E-10)) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    hdummy->Error("Rebin", "a lowedge is not found. lowedge=%f", edge1);
                    return hdummy;
                }
            }

        } else {
            lowedges.clear();
            lowedges.resize(totalnbins+1);
            int iedge = 0;
            
            for (int ipart=0; ipart < npartitions; ipart++) {
                UInt_t newnbins = nbinsvec.at(ipart);
                UInt_t oldnbins = in_nbins/npartitions;
                if (newnbins < minoldnbins || newnbins > oldnbins) {
                    hdummy->Error("Rebin", "newnbins must be within [%i,%i]. newnbins=%i", minoldnbins, oldnbins, newnbins);
                    return hdummy;
                }

                const UInt_t begin = (ipart * oldnbins) + 1;
                const UInt_t end   = ((ipart+1) * oldnbins) + 1;
                
                double olderrorffirst = errorffirst;  // FIXME
                double olderrorflast  = errorflast;   // FIXME
                if (ipart != npartitions-1) {
                    errorffirst = 0.15;  // FIXME
                    errorflast  = errorffirst;  // FIXME
                }
                
                /// Group bins from the left
                UInt_t firstN = 1;  // starts with 1 bin
                findFirstN(firstN, begin, bkg, errorffirst);
                
                /// Group bins from the right
                UInt_t lastN = 1;  // starts with 1 bin
                findLastN(lastN, end, bkg, errorflast);
                
#if !defined(REBINFIRSTTHREE) && !defined(REBINLASTTHREE)
                /// Group bins in the middle
                int nmiddlebins = oldnbins - firstN - lastN;
                int ngroup = newnbins - 2;
                int nbg = nmiddlebins / ngroup;
                //hdummy->Info("Rebin", "firstN=%2i lastN=%2i nmiddlebins=%i before adjusting nmiddlebins", firstN, lastN, nmiddlebins);
                if (nmiddlebins < ngroup){
                    hdummy->Error("Rebin", "Not enough number of bins! Perhaps the errorffirst or errorflast is too low or newnbins is too many. errorffirst=%f, errorflast=%f, newnbins=%i", errorffirst, errorflast, newnbins);
                    return hdummy;
                } else {
                    // Group more bins from the left until nmiddlebins can be divided by ngroup
                    while ((nbg*ngroup != nmiddlebins)){
                        firstN++;  
                        nmiddlebins = oldnbins - firstN - lastN;
                        nbg = nmiddlebins / ngroup;
                    }
                }
                //hdummy->Info("Rebin", "firstN=%2i lastN=%2i nmiddlebins=%i  after adjusting nmiddlebins", firstN, lastN, nmiddlebins);

                if (ipart == 0) {  // the lowest edge
                    lowedges.at(iedge) = in_xaxis->GetBinLowEdge(begin);
                    iedge++;
                }
                for (UInt_t i=1; i<newnbins; i++){
                    UInt_t ibin = begin + firstN + (i-1)*nbg;
                    lowedges.at(iedge) = in_xaxis->GetBinLowEdge(ibin);
                    iedge++;
                    if (i == newnbins-1)
                        assert(ibin == (end - lastN));
                }
                lowedges.at(iedge) = in_xaxis->GetBinLowEdge(end);
                iedge++;

#elif defined(REBINLASTTHREE)
                /// Group bins from the right
                UInt_t lastN2 = 1;  // starts with 1 bin
                findLastN2(lastN2, end - lastN, (1.0 - bkg->Integral(end - lastN, end-1)/bkg->Integral(begin, end-1))*0.99, bkg->Integral(begin, end-1 - lastN), bkg, errorflast);
                //findLastN2(lastN2, end - lastN, sig->Integral(end - lastN, end-1)/sig->Integral(begin, end-1)*1.30, sig->Integral(begin, end-1 - lastN));
                
                UInt_t lastN3 = 1;  // starts with 1 bin
                findLastN2(lastN3, end - lastN - lastN2, (1.0 - bkg->Integral(end - lastN - lastN2, end-1 - lastN)/bkg->Integral(begin, end-1 - lastN))*0.99, bkg->Integral(begin, end-1 - lastN - lastN2));
                //findLastN2(lastN3, end - lastN - lastN2, sig->Integral(end - lastN - lastN2, end-1 - lastN)/sig->Integral(begin, end-1 - lastN)*1.50, sig->Integral(begin, end-1 - lastN - lastN2));
                
                std::cout << (1.0 - bkg->Integral(end - lastN, end-1)/bkg->Integral(begin, end-1)) << std::endl;
                std::cout << (1.0 - bkg->Integral(end - lastN - lastN2, end-1 - lastN)/bkg->Integral(begin, end-1 - lastN)) << std::endl;
                std::cout << (1.0 - bkg->Integral(end - lastN - lastN2 - lastN3, end-1 - lastN - lastN2)/bkg->Integral(begin, end-1 - lastN - lastN2)) << std::endl;
                
                /// Group bins in the middle
                int nmiddlebins = oldnbins - firstN - lastN - lastN2 - lastN3;
                int ngroup = newnbins-4;
                int nbg = nmiddlebins / ngroup;
                //hdummy->Info("Rebin", "firstN=%2i lastN=%2i nmiddlebins=%i before adjusting nmiddlebins", firstN, lastN, nmiddlebins);
                if (nmiddlebins < ngroup) {
                    hdummy->Error("Rebin", "Not enough number of bins! Perhaps the errorffirst or errorflast is too low or newnbins is too many. errorffirst=%f, errorflast=%f, newnbins=%i", errorffirst, errorflast, newnbins);
                    return hdummy;
                } else {
                    // Group more bins from the left until nmiddlebins can be divided by ngroup
                    while ((nbg*ngroup != nmiddlebins)){
                        firstN++;  
                        nmiddlebins = oldnbins - firstN - lastN - lastN2 - lastN3;
                        nbg = nmiddlebins / ngroup;
                    }
                }
                //hdummy->Info("Rebin", "firstN=%2i lastN=%2i nmiddlebins=%i  after adjusting nmiddlebins", firstN, lastN, nmiddlebins);

                if (ipart == 0) {  // the lowest edge
                    lowedges.at(iedge) = in_xaxis->GetBinLowEdge(begin);
                    iedge++;
                }
                for (UInt_t i=1; i<newnbins-2; i++){
                    UInt_t ibin = begin + firstN + (i-1)*nbg;
                    lowedges.at(iedge) = in_xaxis->GetBinLowEdge(ibin);
                    iedge++;
                    if (i == newnbins-2-1)
                        assert(ibin == (end - lastN - lastN2 - lastN3));
                }
                lowedges.at(iedge) = in_xaxis->GetBinLowEdge(end - lastN - lastN2);
                iedge++;
                lowedges.at(iedge) = in_xaxis->GetBinLowEdge(end - lastN);
                iedge++;
                lowedges.at(iedge) = in_xaxis->GetBinLowEdge(end);
                iedge++;

#elif defined(REBINFIRSTTHREE)
                /// Group bins from the left
                UInt_t firstN2 = 1;  // starts with 1 bin
                findFirstN(firstN2, begin + firstN, bkg, errorffirst*0.99);
                
                UInt_t firstN3 = 1;  // starts with 1 bin
                findFirstN(firstN3, begin + firstN + firstN2, bkg, errorffirst*0.99*0.99);
                
                /// Reset when errorffirst == errorflast
                if(TMath::AreEqualRel(errorffirst, errorflast,1.E-7)) {
                    firstN2 = 0;
                    firstN3 = 0;
                }
                
                double integral=0., error=0.;
                integral = bkg->IntegralAndError(begin, begin+firstN-1, error);
                std::cout << firstN << ": " << error/integral << std::endl;
                integral = bkg->IntegralAndError(begin+firstN, begin+firstN+firstN2-1, error);
                std::cout << firstN2 << ": " << error/integral << std::endl;
                integral = bkg->IntegralAndError(begin+firstN+firstN2, begin+firstN+firstN2+firstN3-1, error);
                std::cout << firstN3 << ": " << error/integral << std::endl;
                
                /// Group bins in the middle
                int nmiddlebins = oldnbins - firstN - firstN2 - firstN3 - lastN;
                int ngroup = newnbins-4;
                if(TMath::AreEqualRel(errorffirst, errorflast,1.E-7)) {
                    ngroup = newnbins-2;
                }
                int nbg = nmiddlebins / ngroup;
                //hdummy->Info("Rebin", "firstN=%2i lastN=%2i nmiddlebins=%i before adjusting nmiddlebins", firstN, lastN, nmiddlebins);
                if (nmiddlebins < ngroup) {
                    hdummy->Error("Rebin", "Not enough number of bins! Perhaps the errorffirst or errorflast is too low or newnbins is too many. errorffirst=%f, errorflast=%f, newnbins=%i", errorffirst, errorflast, newnbins);
                    return hdummy;
                } else {
                    // Group more bins from the left until nmiddlebins can be divided by ngroup
                    while ((nbg*ngroup != nmiddlebins)){
                        firstN++;
                        nmiddlebins = oldnbins - firstN - firstN2 - firstN3 - lastN;
                        nbg = nmiddlebins / ngroup;
                    }
                }
                //hdummy->Info("Rebin", "firstN=%2i lastN=%2i nmiddlebins=%i  after adjusting nmiddlebins", firstN, lastN, nmiddlebins);

                UInt_t firstI = 1;
                if (ipart == 0) {  // the lowest edges
                    lowedges.at(iedge) = in_xaxis->GetBinLowEdge(begin);
                    iedge++;
                    
                    if(!TMath::AreEqualRel(errorffirst, errorflast,1.E-7)) {
                        lowedges.at(iedge) = in_xaxis->GetBinLowEdge(begin + firstN);
                        iedge++;
                        firstI++;
                        lowedges.at(iedge) = in_xaxis->GetBinLowEdge(begin + firstN + firstN2);
                        iedge++;
                        firstI++;
                    }
                }
                for (UInt_t i=firstI; i<newnbins; i++){
                    UInt_t ibin = begin + firstN + firstN2 + firstN3 + (i-firstI)*nbg;
                    lowedges.at(iedge) = in_xaxis->GetBinLowEdge(ibin);
                    iedge++;
                    if (i == newnbins-1)
                        assert(ibin == (end - lastN));
                }
                lowedges.at(iedge) = in_xaxis->GetBinLowEdge(end);
                iedge++;
#endif
                hdummy->Info("Rebin", "Transform %i bins into %i bins {|%i|,|%i|x%i,|%i|}", oldnbins, newnbins, firstN, nbg, ngroup, lastN);
                hdummy->Info("Rebin", "Last 3 bins correspond to %.3f, %.3f, %.3f, %.3f", lowedges.at(iedge-4), lowedges.at(iedge-3), lowedges.at(iedge-2), lowedges.at(iedge-1));
                
                //for(UInt_t ie=0; ie<lowedges.size(); ie++)
                //    std::cout << lowedges.at(ie) << std::endl;
                
                errorffirst = olderrorffirst;  // FIXME
                errorflast  = olderrorflast;   // FIXME
            }
            
            assert(iedge == totalnbins+1);
            recreate = false;
        }  // end if else
    
    
        //for (int iedge=0; iedge < totalnbins+1; iedge++)
        //    std::cout << "iedge " << iedge << ": " << lowedges.at(iedge) << std::endl;

        TH1F * hnew = (TH1F *) h->Clone("hnew");  // temp histogram with the same binning
        hnew->Sumw2();
        TH1F * hnewrebin = (TH1F *) hnew->Rebin(totalnbins, "hnewrebin", &lowedges[0]);  // temp histogram with the new binning of variable size
        hnewrebin->Sumw2();
    
        // Sanity checks
        const TAxis* anewrebin = hnewrebin->GetXaxis();
        if (! TMath::AreEqualRel(in_xaxis->GetXmin(), anewrebin->GetXmin(),1.E-10) ||
            ! TMath::AreEqualRel(in_xaxis->GetXmax(), anewrebin->GetXmax(),1.E-10) ) {
            hdummy->Error("Rebin", "Axis limits have changed! Please debug! h=%s", h->GetName());
            return hdummy;
        }
        if (TMath::AreEqualRel(xlow, xup, 1E-10)) {  // if newxlow == newxup, keep the boundaries
            xlow = in_xaxis->GetXmin();
            xup = in_xaxis->GetXmax();
        }
        TString newtitle = Form("%s;%s;%s", h->GetTitle(), h->GetXaxis()->GetTitle(), h->GetYaxis()->GetTitle());
        TH1F * hnewrebinfixed = new TH1F(newname, newtitle, totalnbins, xlow, xup);  // final histogram with the new binning of fixed size
        hnewrebinfixed->Sumw2();
        for (Int_t i=0; i < totalnbins+2; i++){  // include UOflow
            hnewrebinfixed->SetBinContent(i, hnewrebin->GetBinContent(i));
            hnewrebinfixed->SetBinError  (i, hnewrebin->GetBinError(i));
        }
        // Sanity checks
        if (! TMath::AreEqualRel(h->GetSumOfWeights(), hnewrebinfixed->GetSumOfWeights(), 1E-7)){  // they can differ slightly
            hdummy->Error("Rebin", "Sum of weights has changed! Please debug! h=%s, %f --> %f", h->GetName(), h->GetSumOfWeights(), hnewrebinfixed->GetSumOfWeights());
            return hdummy;
        }
        Double_t oldsumoferrors2 = 0., newsumoferrors2 = 0.;
        for (Int_t i=0; i < in_nbins+2; i++){  // include UOflow
            oldsumoferrors2 += (h->GetBinError(i) * h->GetBinError(i));
            if (i < totalnbins+2)
                newsumoferrors2 += (hnewrebinfixed->GetBinError(i) * hnewrebinfixed->GetBinError(i));
        }
        if (! TMath::AreEqualRel(oldsumoferrors2, newsumoferrors2, 1E-10)){
            hdummy->Error("Rebin", "Sum of errors squared has changed! Please debug! h=%s, %f --> %f", h->GetName(), oldsumoferrors2, newsumoferrors2);
            return hdummy;
        }
        double firstbinerrorf = (hnewrebinfixed->GetBinError(1) / hnewrebinfixed->GetBinContent(1));
        double lastbinerrorf  = (hnewrebinfixed->GetBinError(totalnbins) / hnewrebinfixed->GetBinContent(totalnbins));
        //hdummy->Info("Rebin", "firstbinerrorf=%f, lastbinerrorf=%f", firstbinerrorf, lastbinerrorf);
        if (firstbinerrorf < 0 || lastbinerrorf < 0)
            hdummy->Error("Rebin", "Insensible error fractions! Please debug! h=%s, firstbinerrorf=%f, lastbinerrorf=%f", h->GetName(), firstbinerrorf, lastbinerrorf);
        delete hdummy;
        delete hnew;
        delete hnewrebin;
        return hnewrebinfixed;

    }

private:
    int npartitions;
    long long partnbins;
    double errorffirst;
    double errorflast;
    double xlow;
    double xup;
    bool recreate;
    int totalnbins;
    const TH1F * sig;
    const TH1F * bkg;
    std::vector<int> nbinsvec;
    std::vector<double> lowedges;
    
};


///_____________________________________________________________________________
/// Normalizer

class Normalizer {
public:
    Normalizer()
     : nsf_res_j(1.), nsf_scale_j(1.), nsf_eff_b(1.), nsf_fake_b(1.) {}
    
    ~Normalizer();
    
    void set_nuisancefactors(const double nsf[4]) {
        nsf_res_j   = nsf[0];
        nsf_scale_j = nsf[1];
        nsf_eff_b   = nsf[2];
        nsf_fake_b  = nsf[3];
        return;
    }
    
    double normalize(const TString& syst, std::vector<TH1*>& histos) {
        if (syst == "NONE") {
            norms.clear();
            for (UInt_t ih = 0; ih < histos.size(); ih++) {
                norms.push_back(histos.at(ih)->GetSumOfWeights());
            }
        } else {
            assert(norms.size() == histos.size());
            double nsf = 1.;
            if (syst == "CMS_vhbb_res_jUp" || syst == "CMS_vhbb_res_jDown")
                nsf = nsf_res_j;
            else if (syst == "CMS_vhbb_scale_jUp" || syst == "CMS_vhbb_scale_jDown")
                nsf = nsf_scale_j;
            else if (syst == "CMS_vhbb_eff_bUp" || syst == "CMS_vhbb_eff_bDown")
                nsf = nsf_eff_b;
            else if (syst == "CMS_vhbb_fake_b_8TeVUp" || syst == "CMS_vhbb_fake_b_8TeVDown")
                nsf = nsf_fake_b;
            else 
                return nsf;
            
            for (UInt_t ih = 0; ih < histos.size(); ih++) {
                double oldsumofweights = histos.at(ih)->GetSumOfWeights();
                double norm = norms.at(ih);
                double newsumofweights = norm * (1.0 + nsf * (oldsumofweights / norm - 1.0));
                histos.at(ih)->Scale(newsumofweights / oldsumofweights);
            }
        }
        return 1.;
    }

private:
    double nsf_res_j;
    double nsf_scale_j;
    double nsf_eff_b;
    double nsf_fake_b;
    std::vector<double> norms;
};


///_____________________________________________________________________________
/// Functions

void setHisto(TH1 * h, const TString& p)
{
    //h->Sumw2();
    if (p == "VH") {
        h->SetLineColor  (2);
        h->SetFillColor  (2);
        h->SetMarkerColor(2);
    } else if (p == "VH_1") {
        h->SetLineColor  (kPink - 4);
        h->SetFillColor  (kPink - 4);
        h->SetMarkerColor(kPink - 4);
    } else if (p == "data_obs" || p == "Data") {
        h->SetMarkerSize(0.8);
        h->SetMarkerStyle(20);
    } else {
        h->SetLineColor(kBlack);
        if (p == "WjLF") {
            h->SetFillColor  (814);
            h->SetMarkerColor(814);
        } else if (p == "WjHFc") {
            h->SetFillColor  (816);
            h->SetMarkerColor(816);
        } else if (p == "WjHFb") {
            h->SetFillColor  (820);
            h->SetMarkerColor(820);
        } else if (p == "ZjLF") {
            h->SetFillColor  (401);
            h->SetMarkerColor(401);
        } else if (p == "ZjHFc") {
            h->SetFillColor  (41);
            h->SetMarkerColor(41);
        } else if (p == "ZjHFb") {
            h->SetFillColor  (5);
            h->SetMarkerColor(5);
        } else if (p == "TT") {
            h->SetFillColor  (596);
            h->SetMarkerColor(596);
        } else if (p == "s_Top") {
            h->SetFillColor  (840);
            h->SetMarkerColor(840);
        } else if (p == "VV") {
            h->SetFillColor  (922);
            h->SetMarkerColor(922);
        } else if (p == "VVHF") {
            h->SetFillColor  (920);
            h->SetMarkerColor(920);
        } else if (p == "QCD") {
            h->SetFillColor  (616);
            h->SetMarkerColor(616);
        }
    }
}

void setHisto_Znn(TH1 * h, const TString& p)
{
    //h->Sumw2();
    if (p == "VH") {
        h->SetLineColor  (kRed);
        h->SetFillColor  (kRed);
        h->SetMarkerColor(kRed);
    } else if (p == "VH_1") {
        h->SetLineColor  (kPink - 4);
        h->SetFillColor  (kPink - 4);
        h->SetMarkerColor(kPink - 4);
    } else if (p == "data_obs" || p == "Data") {
        h->SetMarkerSize(0.8);
        h->SetMarkerStyle(20);
    } else {
        h->SetLineColor(kBlack);
        if (p == "WjLF") {
            h->SetFillColor  (kSpring - 6);
            h->SetMarkerColor(kSpring - 6);
        } else if (p == "WjHFc") {
            h->SetFillColor  (kGreen - 3);
            h->SetMarkerColor(kGreen - 3);
        } else if (p == "WjHFb") {
            h->SetFillColor  (kSpring);
            h->SetMarkerColor(kSpring);
        } else if (p == "ZjLF") {
            h->SetFillColor  (kYellow + 2);
            h->SetMarkerColor(kYellow + 2);
        } else if (p == "ZjHFc") {
            h->SetFillColor  (kYellow - 3);
            h->SetMarkerColor(kYellow - 3);
        } else if (p == "ZjHFb") {
            h->SetFillColor  (kYellow);
            h->SetMarkerColor(kYellow);
        } else if (p == "TT") {
            h->SetFillColor  (kBlue);
            h->SetMarkerColor(kBlue);
        } else if (p == "s_Top") {
            h->SetFillColor  (kTeal);
            h->SetMarkerColor(kTeal);
        } else if (p == "VV") {
            h->SetFillColor  (kGray + 2);
            h->SetMarkerColor(kGray + 2);
        } else if (p == "VVHF") {
            h->SetFillColor  (kGray);
            h->SetMarkerColor(kGray);
        } else if (p == "QCD") {
            h->SetFillColor  (kMagenta + 1);
            h->SetMarkerColor(kMagenta + 1);
        }
    }
}


void FormatFileName(TString& str)
{
    str.ReplaceAll(" ", "");
    const char *specials = ".[]{}()$*/?&";
    if (str.First(specials) != -1) {
        str.ReplaceAll(".","");
        str.ReplaceAll("[","_");
        str.ReplaceAll("]","_");
        str.ReplaceAll("{","_");
        str.ReplaceAll("}","_");
        str.ReplaceAll("(","_");
        str.ReplaceAll(")","_");
        str.ReplaceAll("$","_");
        str.ReplaceAll("*","_");
        str.ReplaceAll("/","_");
        str.ReplaceAll("?","_");
        str.ReplaceAll("&","_");
    }
    str.Remove(TString::kTrailing, '_');
    str.ReplaceAll("<=","leq");
    str.ReplaceAll(">=","geq");
    str.ReplaceAll("<","lt");
    str.ReplaceAll(">","gt");
    if (str.Sizeof() > 200)  str = str(0,200);
}

// Check consistency between two histograms
bool CheckConsistency(const TH1 * h1, const TH1 * h2)
{
    if (h1->GetNbinsX() != h2->GetNbinsX()){
        h1->Error("CheckConsistency", "Two input histograms have different number of bins!");
        return false;
    }
    const TAxis * a1 = h1->GetXaxis();
    const TAxis * a2 = h2->GetXaxis();
    if (! TMath::AreEqualRel(a1->GetXmin(), a2->GetXmin(),1.E-12) ||
        ! TMath::AreEqualRel(a1->GetXmax(), a2->GetXmax(),1.E-12) ) {
        h1->Error("CheckConsistency", "Two input histograms have different axis limits!");
        return false;
    }
    const TArrayD * h1Array = a1->GetXbins(); 
    const TArrayD * h2Array = a2->GetXbins(); 
    Int_t fN = h1Array->fN;
    if ( fN != 0 ) {
        if ( h2Array->fN != fN ) {
            h1->Error("CheckConsistency", "Two input histograms have different bin limits!");
            return false;
        }
        else {
            for ( int i = 0; i < fN; ++i ) {
                if ( ! TMath::AreEqualRel( h1Array->GetAt(i), h2Array->GetAt(i), 1E-10 ) ) {
                    h1->Error("CheckConsistency", "Two input histograms have different bin limits!");
                    return false;
                }
            }
        }
    }
    return true;
}

TH1F * VaryStatErrors(const TH1F * h, const Double_t c, const char * newname="", bool doUOflow=false)
{
    TH1F * hdummy = (TH1F *) h->Clone("dummy");
    hdummy->Sumw2();
    Int_t nbins = h->GetNbinsX();
    if (TMath::Abs(c*c - 1.0) > 1e-10){
        hdummy->Error("Vary", "c must be +1 or -1!");
        return hdummy;
    }
    TH1F * hnew = (TH1F *) h->Clone(newname);
    if (h->GetSumw2N() == 0){
        hdummy->Warning("Vary", "Sumw2 is not created prior to this.");
        hnew->Sumw2();
    }
    if (h->GetSumOfWeights() < 1e-10){
        hdummy->Warning("Vary", "Empty histogram. Return a clone. h=%s", h->GetName());
        return hnew;
    }
    UInt_t begin=1, end=nbins+1;
    if (doUOflow){  // do underflow and overflow
        begin=0; end=nbins+2;
    } else if ((h->GetBinContent(0)>1e-6) || (h->GetBinContent(nbins+1)>1e-6)){
        hdummy->Warning("Vary", "Underflow and/or overflow are not empty and not modified.");
    }
    Double_t sumoferrors2 = 0.;
    for (UInt_t i=begin; i<end; i++){
        const Double_t bincontent     = h->GetBinContent(i);
        const Double_t binerror       = h->GetBinError(i);
        const Double_t newbincontent  = TMath::Max(0., bincontent + (binerror * c));
        // delta error = sqrt(old error)
        // new error = old error**2 +/- delta error**2
        //           = old error**2 +/- old error
        const Double_t newbinerror2   = TMath::Max(0., (binerror * binerror) + (binerror * c));
        const Double_t newbinerror    = TMath::Sqrt(newbinerror2);
        sumoferrors2                 += (binerror * binerror);

        hnew->SetBinContent(i, newbincontent);
        hnew->SetBinError  (i, newbinerror);
    }
    // new normalization = old normalization +/- sqrt of sum of old stat errors
    Double_t newsumofweights = h->GetSumOfWeights() + (TMath::Sqrt(sumoferrors2) * c);
    if (newsumofweights < 0.){
        hdummy->Warning("Vary", "Negative new sum of weights!");
        delete hdummy;
        return hnew;
    }
    hnew->Scale(newsumofweights / (hnew->GetSumOfWeights()));
    delete hdummy;
    return hnew;
}

// N.B. other channels might take absolute difference
TH1F * VaryModelErrors(const TH1F * h, const TH1F * hmodel, Double_t c, const char* newname="", bool absdiff=false, bool doUOflow=false)
{
    TH1F * hdummy = (TH1F *) h->Clone("dummy");
    hdummy->Sumw2();
    Int_t nbins = h->GetNbinsX();
    //if (TMath::Abs(c*c - 1.0) > 1e-10){
    //    hdummy->Error("Vary", "c must be +1 or -1!");  // allow c to be any constant?
    //    return hdummy;
    //}
    if (h == hmodel){
        hdummy->Warning("Vary", "The input histograms are the same. Return a clone. h=%s", h->GetName());
        hdummy->SetName(newname);
        return hdummy;
    } else {
        if(!CheckConsistency(h, hmodel)) {
            hdummy->Error("Vary", "The input histograms are not consistent! h=%s", h->GetName());
            return hdummy;
        }
    }
    TH1F * hnew = (TH1F *) h->Clone(newname);
    if (h->GetSumw2N() == 0){
        hdummy->Warning("Vary", "Sumw2 is not created prior to this.");
        hnew->Sumw2();
    }
    if (h->GetSumOfWeights() < 1e-10){
        hdummy->Warning("Vary", "Empty histogram. Return a clone. h=%s", h->GetName());
        hdummy->SetName(newname);
        return hnew;
    }
    if (hmodel->GetSumOfWeights() < 1e-10){
        hdummy->Warning("Vary", "Empty model histogram. Return a clone. hmodel=%s", hmodel->GetName());
        hdummy->SetName(newname);
        return hnew;
    }
    UInt_t begin=1, end=nbins+1;
    if (doUOflow){  // do underflow and overflow
        begin=0; end=nbins+2;
    } else if ((h->GetBinContent(0)>1e-6) || (h->GetBinContent(nbins+1)>1e-6)){
        hdummy->Warning("Vary", "Underflow and/or overflow are not empty and not modified.");
    }
    // Normalize hmodel to equal area
    TH1F * hmodelnorm = (TH1F *) hmodel->Clone(Form("%s_norm", hmodel->GetName()));
    if (hmodel->GetSumw2N() == 0){
        hdummy->Warning("Vary", "Sumw2 is not created prior to this.");
        hmodelnorm->Sumw2();
    }
    hmodelnorm->Scale(h->GetSumOfWeights() / (hmodelnorm->GetSumOfWeights()));
    for (UInt_t i=begin; i<end; i++){
        const Double_t bincontent     = h->GetBinContent(i);
        const Double_t binerror       = h->GetBinError(i);
        const Double_t bincontentm    = hmodelnorm->GetBinContent(i);
        const Double_t binerrorm      = hmodelnorm->GetBinError(i);
        Double_t       diff           = bincontentm - bincontent;
        if (absdiff)   diff           = fabs(diff);
        const Double_t newbincontent  = TMath::Max(0., bincontent + (0.5 * diff * c));
        // the operation is h1 +/- 0.5 * (h2 - h1)
        // the bin error is simply sqrt(((1.0 +/- 0.5)**2 * e1**2) + (0.5**2 * e2**2))
        const Double_t newbinerror2   = TMath::Max(0., ((1.0+(0.5*c))*(1.0+(0.5*c)) * (binerror * binerror)) + (0.5*0.5 * (binerrorm * binerrorm)) );
        const Double_t newbinerror    = TMath::Sqrt(newbinerror2);

        hnew->SetBinContent(i, newbincontent);
        hnew->SetBinError  (i, newbinerror);
    }
    // Normalize hnew to equal area
    hnew->Scale(h->GetSumOfWeights() / (hnew->GetSumOfWeights()));
    delete hdummy;
    return hnew;
}

