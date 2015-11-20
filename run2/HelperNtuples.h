const std::string tagMC        = "root://eoscms//eos/cms//store/user/capalmer/VHBBHeppyNtuples/V7/";
const std::string tagData      = "root://eoscms//eos/cms//store/user/capalmer/VHBBHeppyNtuples/V7/";
//const std::string baseline     = "(Vtype==4||Vtype==3||Vtype==2) && met_pt>100 && Sum$(abs(Jet_eta)<2.5 && Jet_pt>20)>1 && Sum$(Jet_btagCSV>0.7)>0";
//const std::string baseline     = "met_pt>200 &&   (Jet_btagCSV[hJidx[0]]>0.9 ||  Jet_btagCSV[hJidx[1]]>0.9) && (Vtype==4) && Jet_pt[hJidx[0]]>30  && Jet_pt[hJidx[1]] > 30  && abs(Jet_eta[hJidx[0]])<2.4 &&  abs(Jet_eta[hJidx[1]])<2.4" ;
const std::string baseline     = "(Vtype>=0 && met_pt>150)" ;



//const std::string regression   = "(Vtype==4||Vtype==3||Vtype==2) && met_pt>120 && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5 && hJet_id[0]==1 && hJet_id[1]==1 && hJet_genPt[0]>10 && hJet_genPt[1]>10 && abs(hJet_flavour[0])==5 && hJet_pt[0]>20 && hJet_pt[1]>20 && hJet_csv[0]>0 && hJet_csv[1]>0";
const std::string regression   = baseline+ "&& ( abs(Jet_mcFlavour[hJCidx[0]])==5  && abs(Jet_mcFlavour[hJCidx[1]])==5 )";
const std::string fjregression = "(Vtype==4||Vtype==3||Vtype==2) && METtype1corr.et>80 && FatH.FatHiggsFlag==1 && nfathFilterJets>0 && abs(fathFilterJets_eta[0])<2.5 && abs(fathFilterJets_eta[1])<2.5 && fathFilterJets_genPt[0]>10 && fathFilterJets_genPt[1]>10 && abs(fathFilterJets_flavour[0])==5 && fathFilterJets_pt[0]>15 && fathFilterJets_pt[1]>15 && fathFilterJets_csv[0]>0 && fathFilterJets_csv[1]>0";
const std::string mettrigger   = "1";
const std::string metfilter    = "1";
const std::string metbtagtrigger = "EVENT.json && ( (190456<=EVENT.run && EVENT.run<=193752 && (triggerFlags[49]==1)) || (193752<=EVENT.run && EVENT.run<=208686 && (triggerFlags[50]==1)) )";
const std::string singlemutrigger = "EVENT.json && ( (190456<=EVENT.run && EVENT.run<=208686 && (triggerFlags[23]==1)) )";
const std::string singleeltrigger = "EVENT.json && ( (190456<=EVENT.run && EVENT.run<=208686 && (triggerFlags[44]==1)) )";

const double lumi = 20000.0;
std::map<std::string, float> GetLumis() {
    std::map<std::string, float> values;
    values["ZnnH125"         ] = lumi *       0.100352 /   100000.0000;
    values["WlnH125"         ] = lumi *       0.259581 /   100000.0000;
    values["monoH"           ] = lumi *       0.069240  /   10000.0000;
    values["WJetsHT100"      ] = lumi *    2234.910000 /  5262249.0000;
    values["WJetsHT200"      ] = lumi *     580.068000 /  4936055.0000;
    values["WJetsHT400"      ] = lumi *      68.400300 /  4640575.0000;
    values["WJetsHT600"      ] = lumi *      23.1363  /   4643671.0000;
    values["WJetsIncl"       ] = lumi *   61623.000000 / 10017431.0000;
    values["ZJetsHT100"      ] = lumi *     409.860000 /  4986410.0000;
    values["ZJetsHT200"      ] = lumi *     110.880000 /  4546455.0000;
    values["ZJetsHT400"      ] = lumi *        13.189 /   4645939.0000;
    values["ZJetsHT600"      ] = lumi *       4.524300 /  4463773.0000;
    values["TTPythia8"       ] = lumi *     809.100000 /  2991597.0000;
    values["TTMadv2"       ] = lumi *     809.100000 /   25245138.0000;
    values["T_s"             ] = lumi *       7.959000 /   499999.0000;
    values["T_t"             ] = lumi *     118.440000 /  3990985.0000;
    values["T_tW"            ] = lumi *      23.310000 /   986096.0000;
    values["Tbar_s"          ] = lumi *       3.696000 /   250000.0000;
    values["Tbar_t"          ] = lumi *      64.470000 /  1999793.0000;
    values["Tbar_tW"         ] = lumi *      23.310000 /   971797.0000;
    values["QCDHT100"        ] = lumi * 28730000.000000 / 4123591.0000;
    values["QCDHT250"        ] = lumi *  670500.000000 /  2668164.0000;
    values["QCDHT500"        ] = lumi *   26740.000000 /  4063331.0000;
    values["QCDHT1000"       ] = lumi *     769.700000 /  1464447.0000;
    values["WJetsHT600"      ] = lumi *      23.136300 /  4643671.0000;
    values["ZJetsHT400"      ] = lumi *      13.189000 /  4645939.0000;
    values["TTMad"           ] = lumi *     809.100000 / 25501279.0000;
    return values;


}


std::vector<std::string> GetLHECuts(const std::string id) {
 
    std::vector<std::string> values;
 
         if (id == "WJets") {
           values.resize(5, "");
        values[0] = "lheHT<100";
        values[1] = "100<=lheHT && lheHT<200";
        values[2] = "200<=lheHT && lheHT<400";
        values[3] = "400<=lheHT && lheHT<600";
        values[4] = "600<=lheHT";
        return values;
    }
         else if (id == "ZJets") {
   
        values.resize(4, "");
        values[0] = "100<=lheHT && lheHT<200";
        values[1] = "200<=lheHT && lheHT<400";
        values[2] = "400<=lheHT && lheHT<600";
        values[3] = "600<=lheHT";
        return values;
    }
    
    return values;
}

std::map<std::string, std::string> GetLHEWeights() {
    std::map<std::string, std::string> values;
    
    values["WJets"           ] = "((lheHT<100) * 5000.0 * 61623.000000 / 10017431.000000 * 1.000000) + ((100<=lheHT && lheHT<200) * 5000.0 * 2234.910000 / 5262249.000000 * 0.322231) + ((200<=lheHT && lheHT<400) * 5000.0 * 580.068000 / 4936055.000000 * 0.393106) + ((400<=lheHT && lheHT<600) * 5000.0 * 68.400300 / 4640575.000000 * 0.820838) + ((600<=lheHT) * 5000.0 * 23.136300 / 4643671.000000 * 0.926380)";
    values["ZJets"           ] = "((100<=lheHT && lheHT<200) * 5000.0 * 409.860000 / 4986410.000000 * 1.000000) + ((200<=lheHT && lheHT<400) * 5000.0 * 110.880000 / 4546455.000000 * 1.000000) + ((400<=lheHT && lheHT<600) * 5000.0 * 13.189000 / 4645939.000000 * 1.000000) + ((600<=lheHT) * 5000.0 * 4.524300 / 4463773.000000 * 1.000000)";
    
    return values;
}

