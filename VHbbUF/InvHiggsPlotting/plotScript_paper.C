{

  ////////////////////////////////////////////////////////////
  // INSTRUCTIONS FOR USE                                   //
  //                                                        //
  // 1. Have ROOT files with hists to be plotted,           //
  //    e.g. ZZ.root,W+jets.root, etc each with a hist hMjj //
  //                                                        //
  // 2. Compile the StackPlot library:                      //
  //      root -q -L StackPlot.cc++                         //
  //    (Do this if you change StackPlot.cc or StackPlot.h) //
  //                                                        //
  // 3. Modify this script!                                 //
  //                                                        //
  // 4. Run this script in batch & quit mode:               //
  //      root -q -b plotScript.C                           //
  //                                                        //
  // 5. Watch as plots are made and saved as PDFs           //
  //                                                        //
  // Any questions, ask Robin: robin.aggleton@cern.ch       //
  ////////////////////////////////////////////////////////////


    // Load StackPlot library file
    gROOT->ProcessLine(".L StackPlot.cc+");

    // Convert to VBF behavior
    gROOT->ProcessLine(".L ZbbToVBF.cc+");
    gInterpreter->GenerateDictionary("vector<vector<TString> >","vector");
    ZbbToVBF zbbToVBF;
    std::vector<std::vector<TString> > zbbNames;
    std::vector<TString> zbbNames_;
    std::vector<TString> vbfNames;

        zbbNames_.push_back("VH");
    zbbNames.push_back(zbbNames_); zbbNames_.clear();
    vbfNames.push_back("VH_hbb");
        //zbbNames_.push_back("ZH");
        zbbNames_.push_back("ZH_hinv");
    zbbNames.push_back(zbbNames_); zbbNames_.clear();
    vbfNames.push_back("ZH_hinv");
        zbbNames_.push_back("Wj0b");
        zbbNames_.push_back("Wj1b");
        zbbNames_.push_back("Wj2b");
    zbbNames.push_back(zbbNames_); zbbNames_.clear();
    vbfNames.push_back("W+jets");
        zbbNames_.push_back("Zj0b");
        zbbNames_.push_back("Zj1b");
        zbbNames_.push_back("Zj2b");
    zbbNames.push_back(zbbNames_); zbbNames_.clear();
    vbfNames.push_back("Z+jets");
        zbbNames_.push_back("TT");
        zbbNames_.push_back("s_Top");
    zbbNames.push_back(zbbNames_); zbbNames_.clear();
    vbfNames.push_back("Top");
        zbbNames_.push_back("VVHF");
        zbbNames_.push_back("VVLF");
    zbbNames.push_back(zbbNames_); zbbNames_.clear();
    vbfNames.push_back("VV");
        zbbNames_.push_back("QCD");
    zbbNames.push_back(zbbNames_); zbbNames_.clear();
    vbfNames.push_back("QCD");
        zbbNames_.push_back("data_obs");
    zbbNames.push_back(zbbNames_); zbbNames_.clear();
    vbfNames.push_back("data_obs");

    TString zbb, vbf;

/*
    // MET
    zbb = "ZnunuHighPt_TT_min_METtype1corret,349_+0";
    vbf = "ZbbHinv_TT_MET";
    zbbToVBF.convert(zbb, vbf, zbbNames, vbfNames, vbf);
    StackPlot metStackPlot(vbf.Data(), "results");
    metStackPlot.setLegPos(0.65, 0.60, 0.91, 0.91);
    metStackPlot.setTextPos(0.18, 0.75, 0.63, 0.91);
    metStackPlot.addDataset("VH_hbb", "VH(b#bar{b})", kPink+1, 0);
    metStackPlot.addDataset("VV", "VV", kAzure+1, 0);
    metStackPlot.addDataset("QCD", "QCD", kMagenta-2, 0);
    metStackPlot.addDataset("Top", "t#bar{t}, tW", kBlue-4, 0);
    metStackPlot.addDataset("W+jets", "W+jets", kSpring, 0);
    metStackPlot.addDataset("Z+jets", "Z+jets", kYellow, 0);
    //metStackPlot.addDataset("ZH_hinv", "ZH(inv)", kRed, 3);  // must be after all backgrounds
    metStackPlot.addDataset("data_obs", "Data", kBlack, 1);
    metStackPlot.setLumi(18.9);
    metStackPlot.setLabel("Z(b#bar{b}) H(inv)  t#bar{t} enriched");
    metStackPlot.draw(vbf.Data(), "E_{T}^{miss} [GeV]", "Events / 15 GeV", 0, 1);
*/

    // ptjj
    zbb = "ZnunuHighPt_WjHF_min_HptReg,999_+0";
    vbf = "ZbbHinv_WjHF_ptjj";
    zbbToVBF.convert(zbb, vbf, zbbNames, vbfNames, vbf);
    StackPlot ptjjStackPlot(vbf.Data(), "results");
    ptjjStackPlot.setLegPos(0.65, 0.60, 0.91, 0.91);
    ptjjStackPlot.setTextPos(0.18, 0.75, 0.63, 0.91);
    ptjjStackPlot.addDataset("VH_hbb", "VH(b#bar{b})", kPink+1, 0);
    ptjjStackPlot.addDataset("VV", "VV", kAzure+1, 0);
    ptjjStackPlot.addDataset("QCD", "QCD", kMagenta-2, 0);
    ptjjStackPlot.addDataset("Top", "t#bar{t}, tW", kBlue-4, 0);
    ptjjStackPlot.addDataset("W+jets", "W+jets", kSpring, 0);
    ptjjStackPlot.addDataset("Z+jets", "Z+jets", kYellow, 0);
    //ptjjStackPlot.addDataset("ZH_hinv", "ZH(inv)", kRed, 3);  // must be after all backgrounds
    ptjjStackPlot.addDataset("data_obs", "Data", kBlack, 1);
    ptjjStackPlot.setLumi(18.9);
    ptjjStackPlot.setLabel("Z(b#bar{b}) H(inv)  W+b#bar{b} enriched");
    ptjjStackPlot.draw(vbf.Data(), "p_{T}^{jj} [GeV]", "Events / 20 GeV", 0, 1);


    // CSVmin
    zbb = "ZnunuHighPt_ZjHF_min_min_hJet_csv_nominal_0_,hJet_csv_nominal_1__,0999_+0";
    vbf = "ZbbHinv_ZjHF_CSVmin";
    zbbToVBF.convert(zbb, vbf, zbbNames, vbfNames, vbf);
    StackPlot csvminStackPlot(vbf.Data(), "results");
    csvminStackPlot.setLegPos(0.65, 0.60, 0.91, 0.91);
    csvminStackPlot.setTextPos(0.18, 0.75, 0.63, 0.91);
    csvminStackPlot.addDataset("VH_hbb", "VH(b#bar{b})", kPink+1, 0);
    csvminStackPlot.addDataset("VV", "VV", kAzure+1, 0);
    csvminStackPlot.addDataset("QCD", "QCD", kMagenta-2, 0);
    csvminStackPlot.addDataset("Top", "t#bar{t}, tW", kBlue-4, 0);
    csvminStackPlot.addDataset("W+jets", "W+jets", kSpring, 0);
    csvminStackPlot.addDataset("Z+jets", "Z+jets", kYellow, 0);
    //csvminStackPlot.addDataset("ZH_hinv", "ZH(inv)", kRed, 3);  // must be after all backgrounds
    csvminStackPlot.addDataset("data_obs", "Data", kBlack, 1);
    csvminStackPlot.setLumi(18.9);
    csvminStackPlot.setLabel("Z(b#bar{b}) H(inv)  Z+b#bar{b} enriched");
    csvminStackPlot.draw(vbf.Data(), "CSV_{min}", "Events / 0.07", 0, 1);


    //__________________________________________________________________________
/*
    // mjj (VH, all)
    zbb = "ZnunuHighPt_VH_HmassReg+0";
    vbf = "ZbbHinv_ZbbHighPt_mjj";
    zbbToVBF.convert(zbb, vbf, zbbNames, vbfNames, vbf);
    StackPlot mjjVH1StackPlot(vbf.Data(), "results");
    mjjVH1StackPlot.setYMin(1e-2);
    mjjVH1StackPlot.setYMax(6e5);
    mjjVH1StackPlot.setLegPos(0.67, 0.60-0.04, 0.93, 0.92);
    mjjVH1StackPlot.setTextPos(0.18, 0.75, 0.63, 0.91);
    mjjVH1StackPlot.addDataset("VH_hbb", "VH(b#bar{b})", kPink+1, 0);
    mjjVH1StackPlot.addDataset("VV", "VV", kAzure+1, 0);
    mjjVH1StackPlot.addDataset("QCD", "QCD", kMagenta-2, 0);
    mjjVH1StackPlot.addDataset("Top", "t#bar{t}, tW", kBlue-4, 0);
    mjjVH1StackPlot.addDataset("W+jets", "W+jets", kSpring, 0);
    mjjVH1StackPlot.addDataset("Z+jets", "Z+jets", kYellow, 0);
    mjjVH1StackPlot.addDataset("ZH_hinv", "#splitline{ZH m_{H}=125GeV,}{B(H #rightarrow inv)=100%}", kRed, 2);  // must be after all backgrounds
    mjjVH1StackPlot.addDataset("data_obs", "Data", kBlack, 1);
    mjjVH1StackPlot.setLumi(18.9);
    mjjVH1StackPlot.setLabel("Z(b#bar{b}) H(inv)  high E_{T}^{miss}");
    mjjVH1StackPlot.draw(vbf.Data(), "M_{jj} [GeV]", "Events / 15 GeV", 1, 1);

    // mjj (VH, after BDT)
    zbb = "ZnunuHighPt_VH_HmassReg+0+0";
    vbf = "ZbbHinv_ZbbHighPt_mjj_afterBDT";
    zbbToVBF.convert(zbb, vbf, zbbNames, vbfNames, vbf);
    StackPlot mjjVH2StackPlot(vbf.Data(), "results");
    mjjVH2StackPlot.setYMin(1e-2);
    mjjVH2StackPlot.setYMax(6e5);
    mjjVH2StackPlot.setLegPos(0.67, 0.60-0.04, 0.93, 0.92);
    mjjVH2StackPlot.setTextPos(0.18, 0.75, 0.63, 0.91);
    mjjVH2StackPlot.addDataset("VH_hbb", "VH(b#bar{b})", kPink+1, 0);
    mjjVH2StackPlot.addDataset("VV", "VV", kAzure+1, 0);
    mjjVH2StackPlot.addDataset("QCD", "QCD", kMagenta-2, 0);
    mjjVH2StackPlot.addDataset("Top", "t#bar{t}, tW", kBlue-4, 0);
    mjjVH2StackPlot.addDataset("W+jets", "W+jets", kSpring, 0);
    mjjVH2StackPlot.addDataset("Z+jets", "Z+jets", kYellow, 0);
    mjjVH2StackPlot.addDataset("ZH_hinv", "#splitline{ZH m_{H}=125GeV,}{B(H #rightarrow inv)=100%}", kRed, 2);  // must be after all backgrounds
    mjjVH2StackPlot.addDataset("data_obs", "Data", kBlack, 1);
    mjjVH2StackPlot.setLumi(18.9);
    mjjVH2StackPlot.setLabel("Z(b#bar{b}) H(inv)  high E_{T}^{miss} BDT>0.2");
    mjjVH2StackPlot.draw(vbf.Data(), "M_{jj} [GeV]", "Events / 15 GeV", 1, 1);

    // CSVmin (VH, all)
    zbb = "ZnunuHighPt_VH_min_min_hJet_csv_nominal_0_,hJet_csv_nominal_1__,0999_+0";
    vbf = "ZbbHinv_ZbbHighPt_CSVmin";
    zbbToVBF.convert(zbb, vbf, zbbNames, vbfNames, vbf);
    StackPlot csvminVH1StackPlot(vbf.Data(), "results");
    csvminVH1StackPlot.setYMin(1e-2);
    csvminVH1StackPlot.setYMax(6e5);
    csvminVH1StackPlot.setLegPos(0.67, 0.60-0.04, 0.93, 0.92);
    csvminVH1StackPlot.setTextPos(0.18, 0.75, 0.63, 0.91);
    csvminVH1StackPlot.addDataset("VH_hbb", "VH(b#bar{b})", kPink+1, 0);
    csvminVH1StackPlot.addDataset("VV", "VV", kAzure+1, 0);
    csvminVH1StackPlot.addDataset("QCD", "QCD", kMagenta-2, 0);
    csvminVH1StackPlot.addDataset("Top", "t#bar{t}, tW", kBlue-4, 0);
    csvminVH1StackPlot.addDataset("W+jets", "W+jets", kSpring, 0);
    csvminVH1StackPlot.addDataset("Z+jets", "Z+jets", kYellow, 0);
    csvminVH1StackPlot.addDataset("ZH_hinv", "#splitline{ZH m_{H}=125GeV,}{B(H #rightarrow inv)=100%}", kRed, 2);  // must be after all backgrounds
    csvminVH1StackPlot.addDataset("data_obs", "Data", kBlack, 1);
    csvminVH1StackPlot.setLumi(18.9);
    csvminVH1StackPlot.setLabel("Z(b#bar{b}) H(inv)  high E_{T}^{miss}");
    csvminVH1StackPlot.draw(vbf.Data(), "CSV_{min}", "Events / 0.07", 1, 1);

    // CSVmin (VH, after BDT)
    zbb = "ZnunuHighPt_VH_min_min_hJet_csv_nominal_0_,hJet_csv_nominal_1__,0999_+0+0";
    vbf = "ZbbHinv_ZbbHighPt_CSVmin_afterBDT";
    zbbToVBF.convert(zbb, vbf, zbbNames, vbfNames, vbf);
    StackPlot csvminVH2StackPlot(vbf.Data(), "results");
    csvminVH2StackPlot.setYMin(1e-2);
    csvminVH2StackPlot.setYMax(6e5);
    csvminVH2StackPlot.setLegPos(0.67, 0.60-0.04, 0.93, 0.92);
    csvminVH2StackPlot.setTextPos(0.18, 0.75, 0.63, 0.91);
    csvminVH2StackPlot.addDataset("VH_hbb", "VH(b#bar{b})", kPink+1, 0);
    csvminVH2StackPlot.addDataset("VV", "VV", kAzure+1, 0);
    csvminVH2StackPlot.addDataset("QCD", "QCD", kMagenta-2, 0);
    csvminVH2StackPlot.addDataset("Top", "t#bar{t}, tW", kBlue-4, 0);
    csvminVH2StackPlot.addDataset("W+jets", "W+jets", kSpring, 0);
    csvminVH2StackPlot.addDataset("Z+jets", "Z+jets", kYellow, 0);
    csvminVH2StackPlot.addDataset("ZH_hinv", "#splitline{ZH m_{H}=125GeV,}{B(H #rightarrow inv)=100%}", kRed, 2);  // must be after all backgrounds
    csvminVH2StackPlot.addDataset("data_obs", "Data", kBlack, 1);
    csvminVH2StackPlot.setLumi(18.9);
    csvminVH2StackPlot.setLabel("Z(b#bar{b}) H(inv)  high E_{T}^{miss} BDT>0.2");
    csvminVH2StackPlot.draw(vbf.Data(), "CSV_{min}", "Events / 0.07", 1, 1);
*/

    // mjj (VH, all)
    zbb = "ZnunuHighPt_VH_HmassReg-0";
    vbf = "ZbbHinv_ZbbHighPt_mjj_forReferee";
    zbbToVBF.convert(zbb, vbf, zbbNames, vbfNames, vbf);
    StackPlot mjjVH1StackPlot(vbf.Data(), "results");
    mjjVH1StackPlot.setYMin(1e-2);
    mjjVH1StackPlot.setYMax(3e6);
    mjjVH1StackPlot.setLegPos(0.67, 0.60-0.04, 0.93, 0.92);
    mjjVH1StackPlot.setTextPos(0.18, 0.75, 0.63, 0.91);
    mjjVH1StackPlot.setArrowPos(250, 2e-3, 200, 3e4);
    mjjVH1StackPlot.addDataset("VH_hbb", "VH(b#bar{b})", kPink+1, 0);
    mjjVH1StackPlot.addDataset("VV", "VV", kAzure+1, 0);
    mjjVH1StackPlot.addDataset("QCD", "QCD", kMagenta-2, 0);
    mjjVH1StackPlot.addDataset("Top", "t#bar{t}, tW", kBlue-4, 0);
    mjjVH1StackPlot.addDataset("W+jets", "W+jets", kSpring, 0);
    mjjVH1StackPlot.addDataset("Z+jets", "Z+jets", kYellow, 0);
    mjjVH1StackPlot.addDataset("ZH_hinv", "#splitline{ZH m_{H}=125GeV,}{B(H #rightarrow inv)=100%}", kRed, 2);  // must be after all backgrounds
    mjjVH1StackPlot.addDataset("data_obs", "Data", kBlack, 1);
    mjjVH1StackPlot.setLumi(18.9);
    mjjVH1StackPlot.setLabel("Z(b#bar{b}) H(inv)  high E_{T}^{miss}");
    mjjVH1StackPlot.draw(vbf.Data(), "M_{jj} [GeV]", "Events / 15 GeV", 1, 1);

    // CSVmin (VH, all)
    zbb = "ZnunuHighPt_VH_min_min_hJet_csv_nominal_0_,hJet_csv_nominal_1__,0999_-0";
    vbf = "ZbbHinv_ZbbHighPt_CSVmin_forReferee";
    zbbToVBF.convert(zbb, vbf, zbbNames, vbfNames, vbf);
    StackPlot csvminVH1StackPlot(vbf.Data(), "results");
    csvminVH1StackPlot.setYMin(1e-2);
    csvminVH1StackPlot.setYMax(3e6);
    csvminVH1StackPlot.setLegPos(0.67, 0.60-0.04, 0.93, 0.92);
    csvminVH1StackPlot.setTextPos(0.18, 0.75, 0.63, 0.91);
    csvminVH1StackPlot.setArrowPos(0.244, 2e-3, 0.4, 3e4);
    csvminVH1StackPlot.addDataset("VH_hbb", "VH(b#bar{b})", kPink+1, 0);
    csvminVH1StackPlot.addDataset("VV", "VV", kAzure+1, 0);
    csvminVH1StackPlot.addDataset("QCD", "QCD", kMagenta-2, 0);
    csvminVH1StackPlot.addDataset("Top", "t#bar{t}, tW", kBlue-4, 0);
    csvminVH1StackPlot.addDataset("W+jets", "W+jets", kSpring, 0);
    csvminVH1StackPlot.addDataset("Z+jets", "Z+jets", kYellow, 0);
    csvminVH1StackPlot.addDataset("ZH_hinv", "#splitline{ZH m_{H}=125GeV,}{B(H #rightarrow inv)=100%}", kRed, 2);  // must be after all backgrounds
    csvminVH1StackPlot.addDataset("data_obs", "Data", kBlack, 1);
    csvminVH1StackPlot.setLumi(18.9);
    csvminVH1StackPlot.setLabel("Z(b#bar{b}) H(inv)  high E_{T}^{miss}");
    csvminVH1StackPlot.draw(vbf.Data(), "CSV_{min}", "Events / 0.05", 1, 1);

    //__________________________________________________________________________
    // BDT
    zbb = "zhinv_Zbb_BDT_ZnunuHighPt_8TeV_prefit";
    vbf = "ZbbHinv_ZbbHighPt_BDT";
    zbbToVBF.convert(zbb, vbf, zbbNames, vbfNames, vbf);
    StackPlot bdtStackPlot1(vbf.Data(), "results");
    bdtStackPlot1.setYMin(1e-2);
    //bdtStackPlot1.setLegPos(0.69, 0.60-0.04, 0.93, 0.91);
    bdtStackPlot1.setLegPos(0.67, 0.60-0.04, 0.93, 0.92);
    bdtStackPlot1.setTextPos(0.18, 0.75, 0.63, 0.91);
    bdtStackPlot1.addDataset("VH_hbb", "VH(b#bar{b})", kPink+1, 0);
    bdtStackPlot1.addDataset("VV", "VV", kAzure+1, 0);
    bdtStackPlot1.addDataset("QCD", "QCD", kMagenta-2, 0);
    bdtStackPlot1.addDataset("Top", "t#bar{t}, tW", kBlue-4, 0);
    bdtStackPlot1.addDataset("W+jets", "W+jets", kSpring, 0);
    bdtStackPlot1.addDataset("Z+jets", "Z+jets", kYellow, 0);
    bdtStackPlot1.addDataset("ZH_hinv", "#splitline{ZH m_{H}=125GeV,}{B(H #rightarrow inv)=100%}", kRed, 2);  // must be after all backgrounds
    bdtStackPlot1.addDataset("data_obs", "Data", kBlack, 1);
    bdtStackPlot1.setLumi(18.9);
    bdtStackPlot1.setLabel("Z(b#bar{b}) H(inv)  high E_{T}^{miss}");
    bdtStackPlot1.draw(vbf.Data(), "BDT output", "Events / 0.1", 1, 1);

    zbb = "zhinv_Zbb_BDT_ZnunuMedPt_8TeV_prefit";
    vbf = "ZbbHinv_ZbbMedPt_BDT";
    zbbToVBF.convert(zbb, vbf, zbbNames, vbfNames, vbf);
    StackPlot bdtStackPlot2(vbf.Data(), "results");
    bdtStackPlot2.setYMin(1e-2);
    //bdtStackPlot2.setLegPos(0.69, 0.60-0.04, 0.93, 0.91);
    bdtStackPlot2.setLegPos(0.67, 0.60-0.04, 0.93, 0.92);
    bdtStackPlot2.setTextPos(0.18, 0.75, 0.63, 0.91);
    bdtStackPlot2.addDataset("VH_hbb", "VH(b#bar{b})", kPink+1, 0);
    bdtStackPlot2.addDataset("VV", "VV", kAzure+1, 0);
    bdtStackPlot2.addDataset("QCD", "QCD", kMagenta-2, 0);
    bdtStackPlot2.addDataset("Top", "t#bar{t}, tW", kBlue-4, 0);
    bdtStackPlot2.addDataset("W+jets", "W+jets", kSpring, 0);
    bdtStackPlot2.addDataset("Z+jets", "Z+jets", kYellow, 0);
    bdtStackPlot2.addDataset("ZH_hinv", "#splitline{ZH m_{H}=125GeV,}{B(H #rightarrow inv)=100%}", kRed, 2);  // must be after all backgrounds
    bdtStackPlot2.addDataset("data_obs", "Data", kBlack, 1);
    bdtStackPlot2.setLumi(18.9);
    bdtStackPlot2.setLabel("Z(b#bar{b}) H(inv)  intermediate E_{T}^{miss}");
    bdtStackPlot2.draw(vbf.Data(), "BDT output", "Events / 0.1", 1, 1);

    zbb = "zhinv_Zbb_BDT_ZnunuLowPt_8TeV_prefit";
    vbf = "ZbbHinv_ZbbLowPt_BDT";
    zbbToVBF.convert(zbb, vbf, zbbNames, vbfNames, vbf);
    StackPlot bdtStackPlot3(vbf.Data(), "results");
    bdtStackPlot3.setYMin(1e-2);
    //bdtStackPlot3.setLegPos(0.69, 0.60-0.04, 0.93, 0.91);
    bdtStackPlot3.setLegPos(0.67, 0.60-0.04, 0.93, 0.92);
    bdtStackPlot3.setTextPos(0.18, 0.75, 0.63, 0.91);
    bdtStackPlot3.addDataset("VH_hbb", "VH(b#bar{b})", kPink+1, 0);
    bdtStackPlot3.addDataset("VV", "VV", kAzure+1, 0);
    bdtStackPlot3.addDataset("QCD", "QCD", kMagenta-2, 0);
    bdtStackPlot3.addDataset("Top", "t#bar{t}, tW", kBlue-4, 0);
    bdtStackPlot3.addDataset("W+jets", "W+jets", kSpring, 0);
    bdtStackPlot3.addDataset("Z+jets", "Z+jets", kYellow, 0);
    bdtStackPlot3.addDataset("ZH_hinv", "#splitline{ZH m_{H}=125GeV,}{B(H #rightarrow inv)=100%}", kRed, 2);  // must be after all backgrounds
    bdtStackPlot3.addDataset("data_obs", "Data", kBlack, 1);
    bdtStackPlot3.setLumi(18.9);
    bdtStackPlot3.setLabel("Z(b#bar{b}) H(inv)  low E_{T}^{miss}");
    bdtStackPlot3.draw(vbf.Data(), "BDT output", "Events / 0.1", 1, 1);

}
