from cut import *

"""
VHbb Analysis Configuration

See also cut.py
"""

#########################
#-- Working Directory --#
#########################

WORK_DIR = '/afs/cern.ch/work/s/swang373/private/V14/'

###############
#-- Samples --#
###############

"""
Properties
----------
path : str
       The path to the Step1 ntuple on EOS. The directory tree from the path
       to the sample's files must contain a single subdirectory at each level,
       as multiple subdirectories may lead to different versions of the sample.
xsec : float, optional
       The cross-section for the MC sample, reported in picobarns (pb).

BEFORE RUNNING sample.py: MAKE SURE SKIM IN cut.py IS CORRECT!

The cross-sections and k-factors were obtained from the following CERN TWiki's:
twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV#s_13_0_TeV
twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopRefXsec

The branching ratios were obtained from the PDG refernce tables.
"""

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
        'xsec': 1.380 * 0.577 * 0.1086 * 3,
    },

    # W+Jets
    #'WJetsIncl': {
    #    'path': '/store/group/phys_higgs/hbb/ntuples/V14/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    #    'xsec': 61526.7,
    #},

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
        'xsec': 1.23 * 280.35, 
    },

    'ZJetsHT200': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/ZJetsToNuNu_HT-200To400_13TeV-madgraph',
        'xsec': 1.23 * 77.67,
    },

    'ZJetsHT400': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/ZJetsToNuNu_HT-400To600_13TeV-madgraph',
        'xsec': 1.23 * 10.73,
    },

    'ZJetsHT600': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/ZJetsToNuNu_HT-600ToInf_13TeV-madgraph',
        'xsec': 1.23 * 4.116,
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
        'xsec': 3.36,
    },

    'T_t_lep': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
        'xsec': 44.33,
    },

    'Tbar_t_lep': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
        'xsec': 26.38,
    },

    # QCD
    'QCDHT100': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'xsec': 27990000, 
    },
    
    'QCDHT200': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'xsec': 1712000,
    },

    'QCDHT300': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'xsec': 347700,
    },

    'QCDHT500': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'xsec': 32100,
    },

    'QCDHT700': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8', 
        'xsec': 6831,
    },

    'QCDHT1000': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'xsec': 1207,
    },

    'QCDHT1500': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'xsec': 119.9,
    },

    'QCDHT2000': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'xsec': 25.24,
    },

    # Diboson
    'WW': {
        'path': '/store/group/phys_higgs/hbb/ntuples/V14/WW_TuneCUETP8M1_13TeV-pythia8',
        'xsec': 63.21, # unverified
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

"""
Properties
----------
samples     : list of str
              The list of samples to be combined into a decay process. 
              Sample names must be valid keys of the SAMPLES dictionary.
sample_cuts : list of str or Cut, optional
              A list of sample-specific cuts applied before combination. 
              Their ordering must match that of the samples property.
process_cut : str or Cut, optional
              A process-specific cut applied to all samples before combination. 
"""

PROCESSES = {

    'data_obs': {
        'samples': ['Data_MET_C', 'Data_MET_D', 'Data_MET_DP'],
        'types': 'data',
    },

    'ZH': {
        'samples': ['ZnnH125'],
        'types': 'MC:sig',
        'color': 632,
    },

    'ggZH': {
        'samples': ['ggZH125'],
        'types': 'MC:sig',
        'color': 632 - 7,
    },

    'WH': {
        'samples': ['WlnH125'],
        'types': 'MC:sig',
        'color': 632 + 2,
    },
    
    's_Top': {
        'samples': ['T_s_comb_lep', 'T_t_lep', 'Tbar_t_lep', 'T_tW', 'Tbar_tW'],
        'types': 'MC:bkg',
        'color': 840,
    },

    'TT': {
        'samples': ['TTPow'],
        'types': 'MC:bkg',
        'color': 600,
    },

    'Wj0b': {
        'samples': ['WJetsIncl', 'WJetsHT100', 'WJetsHT200', 'WJetsHT400', 'WJetsHT600'],
        'sample_cuts': ['lheHT<100'],
        'process_cut': Vudsg | Vcc,
        'types': 'MC:bkg',
        'color': 820 - 6,
    },

    'Wj1b': {
        'samples': ['WJetsIncl', 'WJetsHT100', 'WJetsHT200', 'WJetsHT400', 'WJetsHT600'],
        'sample_cuts': ['lheHT<100'],
        'process_cut': Vb,
        'types': 'MC:bkg',
        'color': 820 - 1,
    },

    'Wj2b': {
        'samples': ['WJetsIncl', 'WJetsHT100', 'WJetsHT200', 'WJetsHT400', 'WJetsHT600'],
        'sample_cuts': ['lheHT<100'],
        'process_cut': Vbb,
        'types': 'MC:bkg',
        'color': 820,
    },

    'Zj0b': {
        'samples': ['ZJetsHT100', 'ZJetsHT200', 'ZJetsHT400', 'ZJetsHT600'],
        'process_cut': Vudsg | Vcc,
        'types': 'MC:bkg',
        'color': 400 + 2,
    },

    'Zj1b': {
        'samples': ['ZJetsHT100', 'ZJetsHT200', 'ZJetsHT400', 'ZJetsHT600'],
        'process_cut': Vb,
        'types': 'MC:bkg',
        'color': 400 - 7,
    },

    'Zj2b': {
        'samples': ['ZJetsHT100', 'ZJetsHT200', 'ZJetsHT400', 'ZJetsHT600'],
        'process_cut': Vbb,
        'types': 'MC:bkg',
        'color': 400,
    },

    'QCD': {
        'samples': ['QCDHT100', 'QCDHT200', 'QCDHT300', 'QCDHT500', 'QCDHT700', 'QCDHT1000', 'QCDHT1500', 'QCDHT2000'],
        'types': 'MC:bkg',
        'color': 613,
    },

    'VVLF': {
        'samples': ['WW', 'WZ', 'ZZ'],
        'process_cut': Vudsg | Vcc,
        'types': 'MC:bkg',
        'color': 920,
    },

    'VVHF': {
        'samples': ['WW', 'WZ', 'ZZ'],
        'process_cut': Vb | Vbb,
        'types': 'MC:bkg',
        'color': 920 + 1,
    },

}

###############
#-- Regions --#
###############

"""
Properties
----------
cuts : list of str or Cut
       The list of cuts defining the signal or control region. Please keep 
       the individual cuts small in length, as the design choice to apply cuts
       sequentially was made because PyROOT crashes when using one large string.
"""

"""
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
        'cuts': [NoQCD],
    },
 
}
"""

#############
#-- Plots --#
#############

# The plots to draw and their properties.
PLOTS = {

    'The mass of the Higgs.': {
        'name': 'H_mass',
        'expression': 'HCSV_mass',
        'x_title': 'm_{H} [GeV]',
        'n_bins': 30,
        'x_min': 0,
        'x_max': 300,
    },

    'The Pt of the Higgs.': {
        'name': 'HCSV_pt',
        'expression': 'HCSV_pt',
        'x_title': 'H p_{T} [GeV]',
        'n_bins': 20,
        'x_min': 150,
        'x_max': 450, 
    },

    'The number of additional jets within a selection.': {
        'name': 'nAddJets',
        'expression': 'Sum$(Jet_pt>30 && Jet_puId &&  abs(Jet_eta)<4.5)',
        'x_title': '# Add. Jets',
        'n_bins': 8,                                                          
        'x_min': 0,
        'x_max': 8,
    },

    'The minimum delta phi between MET and the Higgs jets.': {
        'name': 'min_dPhi_MET_hJ',
        'expression': 'min(abs(TVector2::Phi_mpi_pi(met_phi - Jet_phi[hJCidx[0]])), abs(TVector2::Phi_mpi_pi(met_phi - Jet_phi[hJCidx[1]])))',
        'x_title': 'Min #||{#Delta#phi(#slash{E}_{T}, H Jet)}',
        'n_bins': 32, 
        'x_min': 0,
        'x_max': 3.2,
    },

    'The minimum delta phi between MET and any jet within a selection.': {
        'name': 'min_dPhi_MET_J',
        'expression': 'MinIf$(abs(TVector2::Phi_mpi_pi(met_phi - Jet_phi)), Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)',
        'x_title': 'Min #||{#Delta#phi(#slash{E}_{T}, Jet)}',
        'n_bins': 32,
        'x_min': 0,
        'x_max': 3.2,
    },

    'The HT for jets with Pt > 30.': {
        'name': 'ht30',
        'expression': 'htJet30',
        'x_title': 'HT [GeV]',
        'n_bins': 50,
        'x_min': 0,
        'x_max': 1000,
    },

    'The larger Higgs jet Pt.': {
        'name': 'hj_maxpt',
        'expression': 'max(Jet_pt[hJCidx[0]], Jet_pt[hJCidx[1]])',
        'x_title': 'Max H Jet p_{T} [GeV]',
        'n_bins': 15,
        'x_min': 0,
        'x_max': 300,
    },

    'The smaller Higgs jet Pt.': {
        'name': 'hj_minpt',
        'expression': 'min(Jet_pt[hJCidx[0]], Jet_pt[hJCidx[1]])',
        'x_title': 'Min H Jet p_{T} [GeV]',
        'n_bins': 15,
        'x_min': 0,
        'x_max': 300,
    },

    'The CSV of the first Higgs jet.': {
        'name': 'CSV1',
        'expression': 'Jet_btagCSV[hJCidx[0]]',
        'x_title': 'H Jet 1 CSV',
        'n_bins': 20,
        'x_min': 0,
        'x_max': 1,
    },

    'The CSV of the second Higgs jet.': {
        'name': 'CSV2',
        'expression': 'Jet_btagCSV[hJCidx[1]]',
        'x_title': 'H Jet 2 CSV',
        'n_bins': 20,
        'x_min': 0,
        'x_max': 1,
    },

    'The Pt of MET.': {
        'name': 'MET_pt',
        'expression': 'met_pt',
        'x_title': '#slash{E}_{T} [GeV]',
        'n_bins': 30,
        'x_min': 150,
        'x_max': 450,
    },
    
    'The number of primary vertices.': {
        'name': 'nPV',
        'expression': 'nPVs',
        'x_title': '# Primary Vertices',
        'n_bins': 30,
        'x_min': 0,
        'x_max': 30,
    },

    'The Pt of tracker MET.': {
        'name': 'tkMet',
        'expression': 'tkMet_pt',
        'x_title': 'Tracker #slash{E}_{T} [GeV]',
        'n_bins': 30,
        'x_min': 0,
        'x_max': 300,
    },

    'The delta phi between the MET and tracker MET.': {
        'name': 'dPhi_MET_tkMET',
        'expression': 'abs(TVector2::Phi_mpi_pi(met_phi - tkMet_phi))',
        'x_title': '#||{#Delta#phi(#slash{E}_{T}, Tracker #slash{E}_{T})}',
        'n_bins': 32,
        'x_min': 0,
        'x_max': 3.2,
    },

}

