import TCutOperators as tco

"""
Analysis Settings
"""

#####################################
# Step1 Ntuple Labels and EOS Paths #
#####################################

step1_ntuples = {
# Datasets
'Data_MET_C' : '/store/group/phys_higgs/hbb/ntuples/V15/MET/VHBB_HEPPY_V15_MET__Run2015C_25ns-05Oct2015-v1',
'Data_MET_D' : '/store/group/phys_higgs/hbb/ntuples/V15/MET/VHBB_HEPPY_V15_MET__Run2015D-05Oct2015-v1',
'Data_MET_DP': '/store/group/phys_higgs/hbb/ntuples/V15/MET/VHBB_HEPPY_V15_MET__Run2015D-PromptReco-v4',
# Signals
'ZnnH125': '/store/group/phys_higgs/hbb/ntuples/V14/ZH_HToBB_ZToNuNu_M125_13TeV_amcatnloFXFX_madspin_pythia8',
'ggZH125': '/store/group/phys_higgs/hbb/ntuples/V14/ggZH_HToBB_ZToNuNu_M125_13TeV_amcatnlo_pythia8',
'WlnH125': '/store/group/phys_higgs/hbb/ntuples/V14/WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8',
# W+Jets
'WJetsHT100': '/store/group/phys_higgs/hbb/ntuples/V14/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'WJetsHT200': '/store/group/phys_higgs/hbb/ntuples/V14/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'WJetsHT400': '/store/group/phys_higgs/hbb/ntuples/V14/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'WJetsHT600': '/store/group/phys_higgs/hbb/ntuples/V14/WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'WJetsIncl' : '/store/group/phys_higgs/hbb/ntuples/V14/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
# Z+Jets
'ZJetsHT100': '/store/group/phys_higgs/hbb/ntuples/V14/ZJetsToNuNu_HT-100To200_13TeV-madgraph',
'ZJetsHT200': '/store/group/phys_higgs/hbb/ntuples/V14/ZJetsToNuNu_HT-200To400_13TeV-madgraph',
'ZJetsHT400': '/store/group/phys_higgs/hbb/ntuples/V14/ZJetsToNuNu_HT-400To600_13TeV-madgraph',
'ZJetsHT600': '/store/group/phys_higgs/hbb/ntuples/V14/ZJetsToNuNu_HT-600ToInf_13TeV-madgraph',
# TT
'TTPow': '/store/group/phys_higgs/hbb/ntuples/V14/TT_TuneCUETP8M1_13TeV-powheg-pythia8',
# Single Top s-channel
# Single Top t-channel
# Single Top Wt-channel
'T_tW'   : '/store/group/phys_higgs/hbb/ntuples/V14/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
'Tbar_tW': '/store/group/phys_higgs/hbb/ntuples/V14/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
# Single Top Leptonic
'T_s_comb_lep': '/store/group/phys_higgs/hbb/ntuples/V14/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1',
'T_t_lep'     : '/store/group/phys_higgs/hbb/ntuples/V14/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
'Tbar_t_lep'  : '/store/group/phys_higgs/hbb/ntuples/V14/ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1', 
# QCD
'QCDHT100' : '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'QCDHT200' : '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'QCDHT300' : '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'QCDHT500' : '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'QCDHT700' : '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'QCDHT1000': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'QCDHT1500': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'QCDHT2000': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
# Diboson
'WW': '/store/group/phys_higgs/hbb/ntuples/V14/WW_TuneCUETP8M1_13TeV-pythia8',
'WZ': '/store/group/phys_higgs/hbb/ntuples/V14/WZ_TuneCUETP8M1_13TeV-pythia8',
'ZZ': '/store/group/phys_higgs/hbb/ntuples/V14/ZZ_TuneCUETP8M1_13TeV-pythia8/VHBB_HEPPY_V14_ZZ_TuneCUETP8M1_13TeV-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151025_083230',}

###################
# Step2 Directory #
###################

step2_dir = '/afs/cern.ch/work/s/swang373/private/V14/'

###################
# Step2 Selection #
###################

step2_selection = '(Vtype>=0 && met_pt>150)'

##############################################
# Process Cross-Sections (in picobarns [pb]) #
##############################################

# Obtained from PDG reference pages and the following links
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
# https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV#s_13_0_TeV
# https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopRefXsec

xsec = {
# Signals
'ZnnH125': (0.8696 - 0.1057) * 0.577 * 0.2,
'ggZH125': 2 * 0.1057 * 0.577 * 0.2,
'WlnH125': 1.380 * 0.577 * 0.1080 * 3,
# W+Jets
'WJetsHT100': 1.21 * 1345,
'WJetsHT200': 1.21 * 359.7,
'WJetsHT400': 1.21 * 48.91,
'WJetsHT600': 1.21 * 18.77,
'WJetsIncl' : 61526.7,
# Z+Jets
'ZJetsHT100': 1.23 * 280.47,
'ZJetsHT200': 1.23 * 78.36,
'ZJetsHT400': 1.23 * 10.94,
'ZJetsHT600': 1.23 * 4.20,
# TT
'TTPow': 831.76,
# Single Top s-channel
# Single Top t-channel
# Single Top Wt-channel
'T_tW'   : 35.6,
'Tbar_tW': 35.6,
# Single Top Leptonic
'T_s_comb_lep': 10.32 * 0.1080 * 3,
'T_t_lep'     : 136.02 * 0.1080 * 3,
'Tbar_t_lep'  : 80.95 * 0.1080 * 3, 
# QCD
'QCDHT100' : 27850000,
'QCDHT200' : 1717000,
'QCDHT300' : 351300,
'QCDHT500' : 31630,
'QCDHT700' : 6802,
'QCDHT1000': 1206,
'QCDHT1500': 120.4,
'QCDHT2000': 25.24,
# Diboson
'WW': 118.7, # couldn't verify this one
'WZ': 47.13,
'ZZ': 16.523,
}

###############################
# Cuts and Region Definitions #
###############################

# Data Luminosity (in inverse picobarns [pb-1])
target_lumi = 1280

# General
preselection = tco.add('json', 'HCSV_reg_pt>150', 'Jet_btagCSV[hJCidx[1]]>0.2', 'abs(TVector2::Phi_mpi_pi(HCSV_reg_phi - met_phi))>1.5', 
                        'min(Jet_pt[hJCidx[0]], Jet_pt[hJCidx[1]])>30', 'HCSV_reg_mass<300', 'HCSV_mass>0')
filters = tco.add('Flag_METFilters', 'Flag_HBHENoiseFilter')
antiQCD = 

# Data Only
data_trigger = 'HLT_PFMET90_PFMHT90_IDTight_v'

# MC Only
mc_weight = tco.multiply('sign(genWeight)', target_lumi, '1./sample_lumi')
mc_trigger = 'HLT_PFMET90_PFMHT90_IDLoose_v'

