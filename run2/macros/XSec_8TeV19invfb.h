#include <map>

#include "TCut.h"

// Skim to reduce ntuple size
const TCut skimZll      = "(Vtype==0||Vtype==1) && V.pt>120 && hJet_pt[0]>20 && hJet_pt[1]>20";
const TCut skimWln      = "(Vtype==2||Vtype==3) && V.pt>130 && hJet_pt[0]>30 && hJet_pt[1]>30";
const TCut skimZnn      = "(Vtype==2||Vtype==3||Vtype==4) && METtype1corr.et>150 && hJet_pt[0]>30 && hJet_pt[1]>30";
                           //< Znn includes single lepton channels to use for control regions
const TCut skimHbb      = "H.HiggsFlag==1 && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5 && hJet_id[0]==1 && hJet_id[1]==1 && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && hJet_csv[0]>0 && hJet_csv[1]>0";

// Select events that pass triggers
//    Single muon triggers
//    "HLT_IsoMu24_v.*", #14
//    "HLT_Mu40_v.*", #21
//    "HLT_Mu40_eta2p1_v.*", #22
//    "HLT_IsoMu24_eta2p1_v.*", #23
//
//    Single electron triggers
//    "HLT_Ele27_WP80_v.*", #44
//
//    Double electron triggers
//    "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v.*", #5
//    "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v.*", #6
//
//    MET triggers
//    "HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned_v.*", #39
//    "HLT_DiCentralJet20_CaloMET65_BTagCSV07_PFMHT80_v.*", #40
//    "HLT_DiCentralPFJet30_PFMET80_BTagCSV07_v.*", #41
//    "HLT_PFMET150_v.*", #42
//    "HLT_DiCentralPFJet30_PFMHT80_v.*", #49
const TCut trigZmm      = "EVENT.json==1 && Vtype==0 && (triggerFlags[14]>0 || triggerFlags[21]>0 || triggerFlags[22]>0 || triggerFlags[23]>0) && hbhe";
const TCut trigZee      = "EVENT.json==1 && Vtype==1 && (triggerFlags[5]>0 || triggerFlags[6]>0) && hbhe";
const TCut trigWmn      = "EVENT.json==1 && Vtype==2 && (triggerFlags[14]>0 || triggerFlags[21]>0 || triggerFlags[22]>0 || triggerFlags[23]>0) && hbhe";
const TCut trigWen      = "EVENT.json==1 && Vtype==3 && (triggerFlags[44]>0) && hbhe";
const TCut trigZnn      = "EVENT.json==1 && (Vtype==2||Vtype==3||Vtype==4) && ( (190456<=EVENT.run && EVENT.run<=193752 && (triggerFlags[42]==1 || triggerFlags[49]==1 || triggerFlags[40]==1)) || (193752<=EVENT.run && EVENT.run<=208686 && (triggerFlags[42]==1 || triggerFlags[39]==1 || triggerFlags[41]==1)) ) && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && !isBadHcalEvent";
                          //< Znn also uses MET filters

// Exclude bad runs in Prompt data
const TCut exclRuns     = "(EVENT.run != 201191) && (EVENT.run < 207883 || EVENT.run > 208307)";
                          //< pass "Golden" JSON selection but exclude runs with pixel misalignment that affects b-tagging

const float lumi = 18938.;
std::map<std::string, float> GetLumis() {
    std::map<std::string, float> values;
    // Equivalent lumi = lumi * xsec / counts-with-pu
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
    values["ZllH110"         ] = lumi *       0.046026 /   999230.7500;
    values["ZllH115"         ] = lumi *       0.038043 /   996796.8750;
    values["ZllH120"         ] = lumi *       0.030826 /   990132.4375;
    values["ZllH125"         ] = lumi *       0.024202 /   969926.1250;
    values["ZllH130"         ] = lumi *       0.018316 /  1000976.6250;
    values["ZllH135"         ] = lumi *       0.013298 /   686971.1250;
    values["ZllH140"         ] = lumi *       0.009220 /   983522.9375;
    values["ZllH145"         ] = lumi *       0.006052 /  1000691.9375;
    values["ZllH150"         ] = lumi *       0.003660 /  1000272.1250;
    values["DYJetsM50"       ] = lumi *    3503.710000 / 28640004.0000;
    values["DYJetsPtZ70"     ] = lumi *      68.003000 /  1413878.7500;
    values["DYJetsPtZ100"    ] = lumi *      44.330000 /  1432132.3750;
    values["WJetsPtW100"     ] = lumi *     297.570000 / 49883400.0000;
    values["ZJetsPtZ100"     ] = lumi *      83.005000 / 14171300.0000;
    values["WW"              ] = lumi *      56.753200 / 10007584.0000;
    values["WZ"              ] = lumi *      33.850000 / 10008316.0000;
    values["ZZ"              ] = lumi *       8.297000 /  9809065.0000;
    values["TTFullLeptMG"    ] = lumi *      26.000000 / 11692603.0000;
    values["TTSemiLeptMG"    ] = lumi *     104.000000 / 25388066.0000;
    values["TTHadronicMG"    ] = lumi *     104.000000 / 24738606.0000;
    values["TTPowheg"        ] = lumi *     234.000000 / 26456339.0000;
    values["T_s"             ] = lumi *       3.790000 /   259923.5000;
    values["T_t"             ] = lumi *      56.400000 /  3144815.0000;
    values["T_tW"            ] = lumi *      11.100000 /   498080.5938;
    values["Tbar_s"          ] = lumi *       1.760000 /   140163.6406;
    values["Tbar_t"          ] = lumi *      30.700000 /  1933425.6250;
    values["Tbar_tW"         ] = lumi *      11.100000 /   493595.9062;

    values["SingleMu_ReReco" ] = 1.000000;
    values["SingleMu_Prompt" ] = 1.000000;
    values["SingleEl_ReReco" ] = 1.000000;
    values["SingleEl_Prompt" ] = 1.000000;
    values["DoubleEl_ReReco" ] = 1.000000;
    values["DoubleEl_Prompt" ] = 1.000000;
    values["MET_ReReco"      ] = 1.000000;
    values["MET_Prompt"      ] = 1.000000;

    return values;
}

