#include <vector>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"

class ZbbToVBF {

  public:
    ZbbToVBF() {}
    ~ZbbToVBF() {}

    void convert(TString zbb, TString vbf,
                 const std::vector<std::vector<TString> >& zbbNames,
                 const std::vector<TString>& vbfNames,
                 TString histoName = "BDT") {
        if (!zbb.EndsWith(".root"))
            zbb += ".root";
        if (!vbf.EndsWith("/"))
            vbf += "/";
        TFile* zbbFile = (TFile*) TFile::Open(zbb, "READ");
        if (gSystem->AccessPathName(vbf))
            gSystem->mkdir(vbf);
        if (gSystem->AccessPathName("results"))
            gSystem->mkdir("results");
        assert(zbbNames.size() == vbfNames.size());
        for (unsigned int i = 0; i < zbbNames.size(); ++i) {
            TH1F* vbfHisto = 0;
            for (unsigned int j = 0; j < zbbNames.at(i).size(); ++j) {
                TString zbbName = zbbNames.at(i).at(j);
                TH1F* zbbHisto = (TH1F*) zbbFile->Get(zbbName);
                //std::cout << zbbName << std::endl;
                //assert(zbbHisto != 0);
                if (zbbHisto == 0)  continue;
                if (j == 0)
                    assert(vbfHisto == 0);
                if (vbfHisto == 0) {
                    vbfHisto = (TH1F*) zbbHisto->Clone(histoName);
                    vbfHisto->Reset();
                }
                vbfHisto->Add(zbbHisto);
            }
            if (vbfHisto == 0)  continue;
            TFile* vbfFile = (TFile*) TFile::Open(vbf+vbfNames.at(i)+".root", "RECREATE");
            vbfHisto->Write();
            vbfFile->Close();
            delete vbfHisto;
            delete vbfFile;
        }
    }

};
