# Combine Data Ntuples
hadd -f Data_MET.root Data_MET_*.root

# Combine QCD Ntuples
hadd -f QCD.root QCDHT*.root

# Combine Single Top Ntuples
hadd -f s_Top.root T_s_comb_lep.root T_t_lep.root Tbar_t_lep.root T_tW.root Tbar_tW.root

# Combine Z+Jets Ntuples
hadd -f ZJets.root ZJetsHT*.root

# Combine Diboson Ntuples
hadd -f VV.root WW.root WZ.root ZZ.root

# Combine W+Jets Ntuples (skim the inclusive sample first if using it)
hadd -f WJets.root skim_WJetsIncl.root WJetsHT*.root
