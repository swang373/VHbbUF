#include <map>

#include "TCut.h"

// Skim to reduce ntuple size
//const TCut skimZll      = "(Vtype==0||Vtype==1) && V.pt>120 && hJet_pt[0]>20 && hJet_pt[1]>20";
//const TCut skimWln      = "(Vtype==2||Vtype==3) && V.pt>130 && hJet_pt[0]>30 && hJet_pt[1]>30";
//const TCut skimZnn      = "(Vtype==2||Vtype==3||Vtype==4) && METtype1corr.et>150 && hJet_pt[0]>30 && hJet_pt[1]>30";
                           //< Znn includes single lepton channels to use for control regions
//const TCut skimHbb      = "H.HiggsFlag==1 && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5 && hJet_id[0]==1 && hJet_id[1]==1 && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && hJet_csv[0]>0 && hJet_csv[1]>0";

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

//const TCut trigZmm      = "EVENT.json==1 && Vtype==0 && (triggerFlags[14]>0 || triggerFlags[21]>0 || triggerFlags[22]>0 || triggerFlags[23]>0) && hbhe";
//const TCut trigZee      = "EVENT.json==1 && Vtype==1 && (triggerFlags[5]>0 || triggerFlags[6]>0) && hbhe";
//const TCut trigWmn      = "EVENT.json==1 && Vtype==2 && (triggerFlags[14]>0 || triggerFlags[21]>0 || triggerFlags[22]>0 || triggerFlags[23]>0) && hbhe";
//const TCut trigWen      = "EVENT.json==1 && Vtype==3 && (triggerFlags[44]>0) && hbhe";
//const TCut trigZnn      = "EVENT.json==1 && (Vtype==2||Vtype==3||Vtype==4) && ( (190456<=EVENT.run && EVENT.run<=193752 && (triggerFlags[42]==1 || triggerFlags[49]==1 || triggerFlags[40]==1)) || (193752<=EVENT.run && EVENT.run<=208686 && (triggerFlags[42]==1 || triggerFlags[39]==1 || triggerFlags[41]==1)) ) && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && !isBadHcalEvent";
                          //< Znn also uses MET filters

// Exclude bad runs in Prompt data
//const TCut exclRuns     = "(EVENT.run != 201191) && (EVENT.run < 207883 || EVENT.run > 208307)";
                          //< pass "Golden" JSON selection but exclude runs with pixel misalignment that affects b-tagging

    // Equivalent lumi = lumi * xsec / counts-with-pu
    const double lumi = 20000.0;
    std::map<std::string, float> GetLumis() {
      std::map<std::string, float> values;
      values["ZnnH125"         ] = lumi *       0.100352 /   101119.0000;
      values["WlnH125"         ] = lumi *       0.259581 /   100803.0000;
      values["monoH"           ] = lumi *       0.120000 /    10000.0000;
      values["WJetsHT100"      ] = lumi *    2234.910000 /  5262249.0000;
      values["WJetsHT200"      ] = lumi *     580.068000 /  4936055.0000;
      values["WJetsHT400"      ] = lumi *      68.400300 /  4640575.0000;
      values["WJetsIncl"       ] = lumi *   61623.000000 / 10017431.0000;
      values["ZJetsHT100"      ] = lumi *     409.860000 /  4986410.0000;
      values["ZJetsHT200"      ] = lumi *     110.880000 /  4546455.0000;
      values["ZJetsHT600"      ] = lumi *       4.524300 /  4463773.0000;
      values["TTPythia8"       ] = lumi *     809.100000 /  2991597.0000;
      values["T_s"             ] = lumi *       7.959000 /   499999.0000;
      values["T_t"             ] = lumi *     118.440000 /  3990985.0000;
      values["T_tW"            ] = lumi *      23.310000 /   986096.0000;
      values["Tbar_s"          ] = lumi *       3.696000 /   250000.0000;
      values["Tbar_t"          ] = lumi *      64.470000 /  1999793.0000;
      values["Tbar_tW"         ] = lumi *      23.310000 /   971797.0000;
      values["QCDHT100"        ] = lumi * 28730000.000000 /  4123591.0000;
      values["QCDHT250"        ] = lumi *  670500.000000 /  2668164.0000;
      values["QCDHT500"        ] = lumi *   26740.000000 /  4063331.0000;
      values["QCDHT1000"       ] = lumi *     769.700000 /  1464447.0000;
      values["WJetsHT600"      ] = lumi *      23.136300 /  4643671.0000;
      values["ZJetsHT400"      ] = lumi *      13.189000 /  4645939.0000;
      values["TTMad"           ] = lumi *     809.100000 / 25501279.0000;
      return values;
    }




