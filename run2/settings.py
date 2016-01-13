from cut import *

"""
ZnnHbb Analysis Configuration
"""

###########################
#-- Working Environment --#
###########################

WORK_DIR = '/afs/cern.ch/work/s/swang373/private/V14/testy/'

###############
#-- Samples --#
###############

# The cross-sections for the MC samples are reported in picobarns (pb).
# They were obtained from PDG reference pages and the following links:
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
# https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV#s_13_0_TeV
# https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopRefXsec

# The Step1 ntuples and their properties.
SAMPLES = { 

    # Datasets
    'Data_MET_C': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V15/MET/VHBB_HEPPY_V15_MET__Run2015C_25ns-05Oct2015-v1',
    },

    'Data_MET_D': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V15/MET/VHBB_HEPPY_V15_MET__Run2015D-05Oct2015-v1',
    },

    'Data_MET_DP': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V15a/MET/VHBB_HEPPY_V15a_MET__Run2015D-PromptReco-v4',
    },

    # Signal
    'ZnnH125': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/ZH_HToBB_ZToNuNu_M125_13TeV_amcatnloFXFX_madspin_pythia8',
        'xsec': (0.8696 - 0.1057) * 0.577 * 0.2,
    },

    'ggZH125': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/ggZH_HToBB_ZToNuNu_M125_13TeV_amcatnlo_pythia8',
        'xsec': 2 * 0.1057 * 0.577 * 0.2,
    },

    'WlnH125': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8',
        'xsec': 1.380 * 0.577 * 0.1080 * 3,
    },

    # W+Jets
    'WJetsIncl': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'xsec': 61526.7,
    },

    'WJetsHT100': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'xsec': 1.21 * 1345,
    },

    'WJetsHT200': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'xsec': 1.21 * 359.7,
    },

    'WJetsHT400': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'xsec': 1.21 * 48.91,
    },

    'WJetsHT600': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'xsec': 1.21 * 18.77,
    },

    # Z+Jets
    'ZJetsHT100': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/ZJetsToNuNu_HT-100To200_13TeV-madgraph', 
        'xsec': 1.23 * 280.47, 
    },

    'ZJetsHT200': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/ZJetsToNuNu_HT-200To400_13TeV-madgraph',
        'xsec': 1.23 * 78.36,
    },

    'ZJetsHT400': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/ZJetsToNuNu_HT-400To600_13TeV-madgraph',
        'xsec': 1.23 * 10.94,
    },

    'ZJetsHT600': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/ZJetsToNuNu_HT-600ToInf_13TeV-madgraph',
        'xsec': 1.23 * 4.20,
    },

    # TTbar
    'TTPow': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/TT_TuneCUETP8M1_13TeV-powheg-pythia8',
        'xsec': 831.76,
    },

    # Single Top
    'T_tW': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
        'xsec': 35.6,
    },

    'Tbar_tW': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1', 
        'xsec': 35.6,
    },

    'T_s_comb_lep': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1',
        'xsec': 10.32 * 0.1080 * 3,
    },

    'T_t_lep': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
        'xsec': 136.02 * 0.1080 * 3,
    },

    'Tbar_t_lep': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
        'xsec': 80.95 * 0.1080 * 3,
    },

    # QCD
    'QCDHT100': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'xsec': 27850000, 
    },
    
    'QCDHT200': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'xsec': 1717000,
    },

    'QCDHT300': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'xsec': 351300,
    },

    'QCDHT500': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'xsec': 31630,
    },

    'QCDHT700': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8', 
        'xsec': 6802,
    },

    'QCDHT1000': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'xsec': 1206,
    },

    'QCDHT1500': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'xsec': 120.4,
    },

    'QCDHT2000': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'xsec': 25.24,
    },

    # Diboson
    'WW': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/WW_TuneCUETP8M1_13TeV-pythia8',
        'xsec': 118.7, # unverified
    },

    'WZ': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/WZ_TuneCUETP8M1_13TeV-pythia8',
        'xsec': 47.13,
    },

    'ZZ': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/ZZ_TuneCUETP8M1_13TeV-pythia8/VHBB_HEPPY_V14_ZZ_TuneCUETP8M1_13TeV-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151025_083230',
        'xsec': 16.523,
    },

}

#################
#-- Processes --#
#################

# The groupings of Step2 ntuples representing a decay process.
PROCESSES = {

    'data_obs': {
        'samples': ['Data_MET_C', 'Data_MET_D', 'Data_MET_DP'],
    },

    'ZH': {
        'samples': ['ZnnH125'],
    },

    'ggZH': {
        'samples': ['ggZH125'],
    },

    'WH': {
        'samples': ['WlnH125'],
    },
    
    's_Top': {
        'samples': ['T_s_comb_lep', 'T_t_lep', 'Tbar_t_lep', 'T_tW', 'Tbar_tW'],
    },

    'TT': {
        'samples': ['TTPow'],
    },

    'Zj0b': {
        'samples': ['ZJetsHT100', 'ZJetsHT200', 'ZJetsHT400', 'ZJetsHT600'],
        'process_cut': Vudsg | Vcc,
    },

    'Zj1b': {
        'samples': ['ZJetsHT100', 'ZJetsHT200', 'ZJetsHT400', 'ZJetsHT600'],
        'process_cut': Vb,
    },

    'Zj2b': {
        'samples': ['ZJetsHT100', 'ZJetsHT200', 'ZJetsHT400', 'ZJetsHT600'],
        'process_cut': Vbb,
    },

    'QCD': {
        'samples': ['QCDHT100', 'QCDHT200', 'QCDHT300', 'QCDHT500', 'QCDHT700', 'QCDHT1000', 'QCDHT1500', 'QCDHT2000'],
    },

    'VVLF': {
        'samples': ['WW', 'WZ', 'ZZ'],
        'process_cut': Vudsg | Vcc,
    },

    'VVHF': {
        'samples': ['WW', 'WZ', 'ZZ'],
        'process_cut': Vb | Vbb,
    },

}

###############
#-- Regions --#
###############

# The signal and control region definitions.
REGIONS = {

    'Signal_Loose': {
        'cuts': [minimal, antiQCD, signal_loose],
    },

    'Signal_Tight': {
        'cuts': [minimal, antiQCD, signal_tight],
    },
    
    'TT': {
        'cuts': [minimal, antiQCD, tt],
    },

    'Z_light': {
        'cuts': [minimal, antiQCD, z_light],
    },

    'Z_bb': {
        'cuts': [minimal, antiQCD, z_bb],
    },

    'W_light': {
        'cuts': [minimal, antiQCD, w_light],
    },

    'W_bb': {
        'cuts': [minimal, antiQCD, w_bb],
    },

    'QCD': {
        'cuts': [minimal, qcd],
    },
 
}

#############
#-- Plots --#
#############

# The plots to draw and their properties.
PLOTS = {

    'H_mass': {
        'expression': 'HCSV_mass',
        'x_title': 'm_{H} [GeV]',
        'n_bins': 30,
        'x_min': 0,
        'x_max': 300,
    },

    'HCSV_pt': {
        'expression': 'HCSV_pt',
        'x_title': 'H p_{T} [GeV]',
        'n_bins': 20,
        'x_min': 150,
        'x_max': 450, 
    },

    'nAddJets': {
        'expression': 'Sum$(Jet_pt>30 && Jet_puId &&  abs(Jet_eta)<4.5)',
        'x_title': '# Add. Jets',
        'n_bins': 8,                                                          
        'x_min': 0,
        'x_max': 8,
    },

    'min_dPhi_MET_hJ': {
        'expression': 'min(abs(TVector2::Phi_mpi_pi(met_phi - Jet_phi[hJCidx[0]])), abs(TVector2::Phi_mpi_pi(met_phi - Jet_phi[hJCidx[1]])))',
        'x_title': 'Min #||{#Delta#phi(#slash{E}_{T}, H Jet)}',
        'n_bins': 32, 
        'x_min': 0,
        'x_max': 3.2,
    },

    'min_dPhi_MET_J': {
        'expression': 'MinIf$(abs(TVector2::Phi_mpi_pi(met_phi - Jet_phi)), Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)',
        'x_title': 'Min #||{#Delta#phi(#slash{E}_{T}, Jet)}',
        'n_bins': 32,
        'x_min': 0,
        'x_max': 3.2,
    },

    'ht30': {
        'expression': 'htJet30',
        'x_title': 'HT [GeV]',
        'n_bins': 50,
        'x_min': 0,
        'x_max': 1000,
    },

    'hj_maxpt': {
        'expression': 'max(Jet_pt[hJCidx[0]], Jet_pt[hJCidx[1]])',
        'x_title': 'Max H Jet p_{T} [GeV]',
        'n_bins': 15,
        'x_min': 0,
        'x_max': 300,
    },

    'hj_minpt': {
        'expression': 'min(Jet_pt[hJCidx[0]], Jet_pt[hJCidx[1]])',
        'x_title': 'Min H Jet p_{T} [GeV]',
        'n_bins': 15,
        'x_min': 0,
        'x_max': 300,
    },

    'CSV1': {
        'expression': 'Jet_btagCSV[hJCidx[0]]',
        'x_title': 'H Jet 1 CSV',
        'n_bins': 20,
        'x_min': 0,
        'x_max': 1,
    },

    'CSV2': {
        'expression': 'Jet_btagCSV[hJCidx[1]]',
        'x_title': 'H Jet 2 CSV',
        'n_bins': 20,
        'x_min': 0,
        'x_max': 1,
    },

    'MET_pt': {
        'expression': 'met_pt',
        'x_title': '#slash{E}_{T} [GeV]',
        'n_bins': 30,
        'x_min': 150,
        'x_max': 450,
    },
    
    'nPV': {
        'expression': 'nPVs',
        'x_title': '# Primary Vertices',
        'n_bins': 30,
        'x_min': 0,
        'x_max': 30,
    },

    'tkMet': {
        'expression': 'tkMet_pt',
        'x_title': 'Tracker #slash{E}_{T} [GeV]',
        'n_bins': 30,
        'x_min': 0,
        'x_max': 300,
    },

    'dPhi_MET_tkMET': {
        'expression': 'abs(TVector2::Phi_mpi_pi(met_phi - tkMet_phi))',
        'x_title': '#||{#Delta#phi(#slash{E}_{T}, Tracker #slash{E}_{T})}',
        'n_bins': 32,
        'x_min': 0,
        'x_max': 3.2,
    },

}

