const std::string tagMC        = "Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC";
const std::string tagData      = "Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_DATA";
const std::string baseline     = "H.HiggsFlag==1 && (Vtype==4||Vtype==3||Vtype==2) && METtype1corr.et>80 && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5 && hJet_id[0]==1 && hJet_id[1]==1 && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && hJet_pt[0]>20 && hJet_pt[1]>20 && hJet_csv[0]>0 && hJet_csv[1]>0 && (hJet_csv_nominal[0]>0.244 || hJet_csv_nominal[1]>0.244)";
const std::string baselineZmmToZbb = "Vtype==0 && abs(vLepton_eta[0])<2.5 && abs(vLepton_eta[1])<2.5 && vLepton_pt[0]>20 && vLepton_pt[1]>20 && deltaR(vLepton_eta[0],vLepton_phi[0],vLepton_eta[1],vLepton_phi[1])>0.5";
const std::string regression   = "(Vtype==4||Vtype==3||Vtype==2) && METtype1corr.et>80 && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5 && hJet_id[0]==1 && hJet_id[1]==1 && hJet_genPt[0]>10 && hJet_genPt[1]>10 && abs(hJet_flavour[0])==5 && hJet_pt[0]>20 && hJet_pt[1]>20 && hJet_csv[0]>0 && hJet_csv[1]>0";
const std::string fjregression = "(Vtype==4||Vtype==3||Vtype==2) && METtype1corr.et>80 && FatH.FatHiggsFlag==1 && nfathFilterJets>0 && abs(fathFilterJets_eta[0])<2.5 && abs(fathFilterJets_eta[1])<2.5 && fathFilterJets_genPt[0]>10 && fathFilterJets_genPt[1]>10 && abs(fathFilterJets_flavour[0])==5 && fathFilterJets_pt[0]>15 && fathFilterJets_pt[1]>15 && fathFilterJets_csv[0]>0 && fathFilterJets_csv[1]>0";
const std::string mettrigger   = "EVENT.json && ( (190456<=EVENT.run && EVENT.run<=193752 && (triggerFlags[42]==1 || triggerFlags[49]==1 || triggerFlags[40]==1)) || (193752<=EVENT.run && EVENT.run<=208686 && (triggerFlags[42]==1 || triggerFlags[39]==1 || triggerFlags[41]==1)) )";
const std::string metfilter    = "hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && !isBadHcalEvent";
const std::string metbtagtrigger = "EVENT.json && ( (190456<=EVENT.run && EVENT.run<=193752 && (triggerFlags[49]==1)) || (193752<=EVENT.run && EVENT.run<=208686 && (triggerFlags[50]==1)) )";
const std::string singlemutrigger = "EVENT.json && ( (190456<=EVENT.run && EVENT.run<=208686 && (triggerFlags[23]==1)) )";
const std::string singleeltrigger = "EVENT.json && ( (190456<=EVENT.run && EVENT.run<=208686 && (triggerFlags[44]==1)) )";

const double lumi = 19624.0;
std::map<std::string, float> GetLumis() {
    std::map<std::string, float> values;
    values["ZbbHinv105"      ] = lumi *       0.106173 /   476181.5312;
    values["ZbbHinv115"      ] = lumi *       0.081013 /   500446.8750;
    values["ZbbHinv125"      ] = lumi *       0.062793 /   500025.3438;
    values["ZbbHinv135"      ] = lumi *       0.049276 /   500053.1562;
    values["ZbbHinv145"      ] = lumi *       0.039055 /   530427.8125;
    values["ZbbHinv150"      ] = lumi *       0.034897 /   484800.0938;
    values["ZbbHinvZmmToZbb" ] = lumi *       0.062793 /   500025.3438;
    values["ZnnH110"         ] = lumi *       0.091140 /  1000249.0000;
    values["ZnnH115"         ] = lumi *       0.075333 /  1000952.0625;
    values["ZnnH120"         ] = lumi *       0.061042 /  1001117.8125;
    values["ZnnH125"         ] = lumi *       0.047926 /  1001196.6875;
    values["ZnnH130"         ] = lumi *       0.036269 /   999389.3750;
    values["ZnnH135"         ] = lumi *       0.026333 /   998509.3750;
    values["ZnnH140"         ] = lumi *       0.018257 /   999974.0625;
    values["ZnnH145"         ] = lumi *       0.011985 /  1000448.0625;
    values["ZnnH150"         ] = lumi *       0.007247 /  1001032.5000;
    values["WlnH110"         ] = lumi *       0.259765 /   992014.7500;
    values["WlnH115"         ] = lumi *       0.212356 /  1001122.0000;
    values["WlnH120"         ] = lumi *       0.170097 /  1001201.1250;
    values["WlnH125"         ] = lumi *       0.132537 /  1000787.5625;
    values["WlnH130"         ] = lumi *       0.099348 /  1001402.9375;
    values["WlnH135"         ] = lumi *       0.071331 /  1000749.6875;
    values["WlnH140"         ] = lumi *       0.048963 /   997183.3750;
    values["WlnH145"         ] = lumi *       0.031886 /   998269.4375;
    values["WlnH150"         ] = lumi *       0.019081 /   994857.8750;
    values["WJetsPtW70"      ] = lumi *     557.570000 / 21950408.0000;
    values["WJetsPtW100"     ] = lumi *     297.570000 / 49883400.0000;
    values["WJetsPtW180"     ] = lumi *      34.292700 /  9702082.0000;
    values["WJetsHW"         ] = lumi *     297.570000 / 11796280.0000;
    values["ZJetsPtZ70"      ] = lumi *     127.610491 / 21928422.0000;
    values["ZJetsPtZ100"     ] = lumi *      83.005000 / 14171300.0000;
    values["ZJetsHT50"       ] = lumi *     495.560000 /  3978979.5000;
    values["ZJetsHT100"      ] = lumi *     208.390000 /  4383989.0000;
    values["ZJetsHT200"      ] = lumi *      51.240150 /  5054296.5000;
    values["ZJetsHT400"      ] = lumi *       6.856200 /   947465.8125;
    values["ZJetsHW"         ] = lumi *      83.005000 /  4873690.0000;
    values["DYJetsM50"       ] = lumi *    3503.710000 / 28640004.0000;
    values["DYJetsZmmToZbb"  ] = lumi *    3503.710000 / 28640004.0000;
    values["WW"              ] = lumi *      56.753200 / 10007584.0000;
    values["WZ"              ] = lumi *      33.850000 / 10008316.0000;
    values["ZZ"              ] = lumi *       8.297000 /  9809065.0000;
    values["TTFullLeptMG"    ] = lumi *      26.000000 / 11692603.0000;
    values["TTSemiLeptMG"    ] = lumi *     104.000000 / 25388066.0000;
    values["TTHadronicMG"    ] = lumi *     104.000000 / 24738606.0000;
    values["TTPowheg"        ] = lumi *     234.000000 / 26456339.0000;
    values["TTMCatNLO"       ] = lumi *     234.000000 / 27519328.0000;
    values["T_s"             ] = lumi *       3.790000 /   259923.5000;
    values["T_t"             ] = lumi *      56.400000 /  3144815.0000;
    values["T_tW"            ] = lumi *      11.100000 /   498080.5938;
    values["Tbar_s"          ] = lumi *       1.760000 /   140163.6406;
    values["Tbar_t"          ] = lumi *      30.700000 /  1933425.6250;
    values["Tbar_tW"         ] = lumi *      11.100000 /   493595.9062;
    values["QCDPt50"         ] = lumi * 8148778.000000 /  4129977.0000;
    values["QCDPt80"         ] = lumi * 1033680.000000 /  5961532.0000;
    values["QCDPt120"        ] = lumi *  156293.300000 /  5741214.5000;
    values["QCDPt170"        ] = lumi *   34138.150000 /  5391437.5000;
    values["QCDPt300"        ] = lumi *    1759.549000 /  5554703.0000;
    values["QCDPt470"        ] = lumi *     113.879100 /  2614799.0000;
    values["QCDPt600"        ] = lumi *      26.992100 /  3964562.7500;
    values["QCDPt800"        ] = lumi *       3.550036 /  3966096.2500;
    values["QCDPt1000"       ] = lumi *       0.737844 /  1962417.7500;
    values["QCDPt1400"       ] = lumi *       0.033522 /  1998678.7500;
    values["QCDPt1800"       ] = lumi *       0.001829 /   978352.1250;
    values["Data_R"          ] = lumi *       1.000000 / 25148913.0000;
    values["Data_P"          ] = lumi *       1.000000 / 46905706.0000;
    values["Data_METBTag_R"  ] = lumi *       1.000000 / 25148913.0000;
    values["Data_METBTag_P"  ] = lumi *       1.000000 / 46905706.0000;
    values["Data_SingleMu_R" ] = lumi *       1.000000 / 25148913.0000;
    values["Data_SingleMu_P" ] = lumi *       1.000000 / 46905706.0000;
    values["Data_SingleEl_R" ] = lumi *       1.000000 / 25148913.0000;
    values["Data_SingleEl_P" ] = lumi *       1.000000 / 46905706.0000;
    return values;
}

std::vector<std::string> GetLHECuts(const std::string id) {
    std::vector<std::string> values;
    if (id == "WJets") {
        values.resize(4, "");
        values[0] = "70<=lheV_pt && lheV_pt<100";
        values[1] = "100<=lheV_pt && lheV_pt<180";
        values[2] = "(180<=lheV_pt && lheV_pt<220 && lheNj>1) || (220<=lheV_pt)";
        values[3] = "180<=lheV_pt && lheV_pt<220 && lheNj==1";
        return values;
    }
    else if (id == "ZJets") {
        values.resize(5, "");
        values[0] = "50<=lheHT && lheHT<100 && lheV_pt<100";
        values[1] = "100<=lheHT && lheHT<200 && lheV_pt<100";
        values[2] = "200<=lheHT && lheHT<400 && lheV_pt<100";
        values[3] = "400<=lheHT && lheV_pt<100";
        values[4] = "100<=lheV_pt";
        return values;
    }
    return values;
}

std::map<std::string, std::string> GetLHEWeights() {
    std::map<std::string, std::string> values;
    values["WJets"           ] = "((70<=lheV_pt && lheV_pt<100) * 19624.0 * 557.570000 / 21950408.000000 * 1.000000) + ((100<=lheV_pt && lheV_pt<180) * 19624.0 * 297.570000 / 49883400.000000 * 1.000000) + (((180<=lheV_pt && lheV_pt<220 && lheNj>1) || (220<=lheV_pt)) * 19624.0 * 34.292700 / 9702082.000000 * 0.627542) + ((180<=lheV_pt && lheV_pt<220 && lheNj==1) * 19624.0 * 297.570000 / 49883400.000000 * 1.000000)";
    values["ZJets"           ] = "((50<=lheHT && lheHT<100 && lheV_pt<100) * 19624.0 * 495.560000 / 3978979.500000 * 1.000000) + ((100<=lheHT && lheHT<200 && lheV_pt<100) * 19624.0 * 208.390000 / 4383989.000000 * 1.000000) + ((200<=lheHT && lheHT<400 && lheV_pt<100) * 19624.0 * 51.240150 / 5054296.500000 * 1.000000) + ((400<=lheHT && lheV_pt<100) * 19624.0 * 6.856200 / 947465.812500 * 1.000000) + ((100<=lheV_pt) * 19624.0 * 83.005000 / 14171300.000000 * 0.775950)";
    return values;
}

