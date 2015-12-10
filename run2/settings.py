import TCutOperators as tco

"""
ZnnHbb Analysis Settings

The idea is to have all configuration options
set centrally and leave the core code alone.


The cross-sections for the MC samples are reported in picobarns (pb).
Obtained from PDG reference pages and the following links
https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV#s_13_0_TeV
https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopRefXsec
"""

###########################
#-- Step2 Configuration --#
###########################

# Output directory for the Step2 ntuples.
STEP2_DIR = '/afs/cern.ch/work/s/swang373/private/V14/Step2/'

# The selection used to skim the Step1 ntuples.
STEP2_CUT = 'Vtype>=0 && met_pt>150'

# The Step1 ntuples and their properties.
SAMPLES = {
    
    # Datasets
    'Data_MET_C': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V15/MET/VHBB_HEPPY_V15_MET__Run2015C_25ns-05Oct2015-v1',
    },

    'Data_MET_D': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V15/MET/VHBB_HEPPY_V15_MET__Run2015D-05Oct2015-v1',
    },

    'Data_MET_DP': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V15a/MET/VHBB_HEPPY_V15a_MET__Run2015D-PromptReco-v4',
    },

    # Signal
    'ZnnH125': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/ZH_HToBB_ZToNuNu_M125_13TeV_amcatnloFXFX_madspin_pythia8',
        'XSEC': (0.8696 - 0.1057) * 0.577 * 0.2,
    },

    'ggZH125': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/ggZH_HToBB_ZToNuNu_M125_13TeV_amcatnlo_pythia8',
        'XSEC': 2 * 0.1057 * 0.577 * 0.2,
    },

    'WlnH125': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8',
        'XSEC': 1.380 * 0.577 * 0.1080 * 3,
    },

    # W Jets
    'WJetsIncl': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'XSEC': 61526.7,
    },

    'WJetsHT100': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'XSEC': 1.21 * 1345,
    },

    'WJetsHT200': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'XSEC': 1.21 * 359.7,
    },

    'WJetsHT400': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'XSEC': 1.21 * 48.91,
    },

    'WJetsHT600': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'XSEC': 1.21 * 18.77,
    },

    # Z Jets
    'ZJetsHT100': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/ZJetsToNuNu_HT-100To200_13TeV-madgraph', 
        'XSEC': 1.23 * 280.47, 
    },

    'ZJetsHT200': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/ZJetsToNuNu_HT-200To400_13TeV-madgraph',
        'XSEC': 1.23 * 78.36,
    },

    'ZJetsHT400': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/ZJetsToNuNu_HT-400To600_13TeV-madgraph',
        'XSEC': 1.23 * 10.94,
    },

    'ZJetsHT600': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/ZJetsToNuNu_HT-600ToInf_13TeV-madgraph',
        'XSEC': 1.23 * 4.20,
    },

    # TTbar
    'TTPow': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/TT_TuneCUETP8M1_13TeV-powheg-pythia8',
        'XSEC': 831.76,
    },

    # Single Top
    'T_tW': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
        'XSEC': 35.6,
    },

    'Tbar_tW': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1', 
        'XSEC': 35.6,
    },

    'T_s_comb_lep': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1',
        'XSEC': 10.32 * 0.1080 * 3,
    },

    'T_t_lep': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
        'XSEC': 136.02 * 0.1080 * 3,
    },

    'Tbar_t_lep': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
        'XSEC': 80.95 * 0.1080 * 3,
    },

    # QCD
    'QCDHT100': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'XSEC': 27850000, 
    },
    
    'QCDHT200': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'XSEC': 1717000,
    },

    'QCDHT300': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'XSEC': 351300,
    },

    'QCDHT500': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'XSEC': 31630,
    },

    'QCDHT700': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8', 
        'XSEC': 6802,
    },

    'QCDHT1000': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'XSEC': 1206,
    },

    'QCDHT1500': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'XSEC': 120.4,
    },

    'QCDHT2000': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'XSEC': 25.24,
    },

    # Diboson
    'WW': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/WW_TuneCUETP8M1_13TeV-pythia8',
        'XSEC': 118.7, # unverified
    },

    'WZ': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/WZ_TuneCUETP8M1_13TeV-pythia8',
        'XSEC': 47.13,
    },

    'ZZ': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V14/ZZ_TuneCUETP8M1_13TeV-pythia8/VHBB_HEPPY_V14_ZZ_TuneCUETP8M1_13TeV-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151025_083230',
        'XSEC': 16.523,
    },

}

###########################
#-- Step3 Configuration --#
###########################

# Output directory for the Step3 ntuples.
STEP3_DIR = '/afs/cern.ch/work/s/swang373/private/V14/Step3/'

# The groupings of Step2 ntuples.
GROUPS = {

    'Data_MET': {
        'SAMPLES': ['Data_MET_C', 'Data_MET_D', 'Data_MET_DP'],
    },

    'ZnnH125': {},
    
    'ggZH125': {},

    'WlnH125': {},
 
    'WJets': {
        'SAMPLES': ['WJetsIncl', 'WJetsHT100', 'WJetsHT200', 'WJetsHT400', 'WJetsHT600'],
        'CUTS': ['lheHT<100'],
    },

    'ZJets': {
        'SAMPLES': ['ZJetsHT100', 'ZJetsHT200', 'ZJetsHT400', 'ZJetsHT600'],
    },
   
    'TT': {
        'SAMPLES': ['TTPow'],
    },

    'ST': {
        'SAMPLES': ['T_s_comb_lep', 'T_t_lep', 'Tbar_t_lep', 'T_tW', 'Tbar_tW'],
    },

    'QCD': {
        'SAMPLES': ['QCDHT100', 'QCDHT200', 'QCDHT300', 'QCDHT500', 'QCDHT700', 'QCDHT1000', 'QCDHT1500', 'QCDHT2000'],
    },

    'VV': {
        'SAMPLES': ['WW', 'WZ', 'ZZ'],
    },

}

###########################
#-- Step4 Configuration --#
###########################

# Output directory for the Step3 ntuples.
STEP4_DIR = '/afs/cern.ch/work/s/swang373/private/V14/Step4/'




##############################
# Control Region Definitions #
##############################

CONTROL_REGION_DIR = '/afs/cern.ch/work/s/swang373/private/V14/CR/'

# Define Useful Cuts
Minimal = tco.add('Vtype==2||Vtype==3||Vtype==4',
                  'HCSV_pt>150',
                  'Jet_btagCSV[hJCidx[1]]>0.3',
                  'min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])>30',
                  'HCSV_mass<300')

AntiQCD = tco.add('MinIf$(abs(TVector2::Phi_mpi_pi(met_phi-Jet_phi)),Jet_pt>30&&Jet_puId&&abs(Jet_eta)<4.5)>0.7',
                  'abs(TVector2::Phi_mpi_pi(met_phi-tkMet_phi))<0.7',
                  'tkMet_pt>30')

addCenJet30m0 = '(Sum$(Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)-2)>0'

addCenJet30e0 = '(Sum$(Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)-2)==0'

addCenJet30e1 = '(Sum$(Jet_pt>30&&Jet_puId&&abs(Jet_eta)<4.5)-2)<=1'

naddGoodLeptons10e0 = '(Sum$(aLeptons_pt>10&&(aLeptons_jetBTagCSV<0.25||aLeptons_relIso03<0.4||aLeptons_looseIdSusy!=0||aLeptons_jetDR>0.3))+Sum$(vLeptons_pt>10&&(vLeptons_jetBTagCSV<0.25||vLeptons_relIso03<0.4||vLeptons_looseIdSusy!=0||vLeptons_jetDR>0.3)))==0'

naddGoodTaus20e0 = 'Sum$(TauGood_idDecayMode>=1&&TauGood_idCI3hit>=1&&TauGood_pt>20&&abs(TauGood_eta)<2.3)==0'

# Jet Flavors
LIGHT_FLAVOR = 'abs(Jet_mcFlavour[hJCidx[0]])!=5 && abs(Jet_mcFlavour[hJCidx[1]])!=5'
HEAVY_FLAVOR = 'abs(Jet_mcFlavour[hJCidx[0]])==5 || abs(Jet_mcFlavour[hJCidx[1]])==5'

# Control Regions
CONTROL_REGIONS = {

    # Signal Regions
    'Signal_Loose': [
        Minimal,
        AntiQCD,
        tco.add('Vtype==4',
                'Jet_btagCSV[hJCidx[1]]>0.605', 
                'HCSV_mass<100 || HCSV_mass>140', 
                naddGoodLeptons10e0, 
                naddGoodTaus20e0)
    ],

    'Signal_Tight': [
        Minimal,
        AntiQCD,
        tco.add('Vtype==4',
                'Jet_btagCSV[hJCidx[1]]>0.8', 
                'HCSV_mass<100||HCSV_mass>140',
                naddGoodLeptons10e0, 
                naddGoodTaus20e0, 
                addCenJet30e1)
    ],
    
    # TTbar Control Region
    'TTbar': [
        Minimal,
        AntiQCD,
        tco.add('Vtype==2 || Vtype==3',
                'vLeptons_pt>30', 
                addCenJet30m0,
                'Jet_btagCSV[hJCidx[0]]>0.97', 
                'Jet_btagCSV[hJCidx[1]]<0.97')
    ],

    # Z Jets Control Regions
    'Z_light': [
        Minimal,
        AntiQCD,
        tco.add('Vtype==4', 
                addCenJet30e0, 
                naddGoodLeptons10e0, 
                'Jet_btagCSV[hJCidx[0]]<0.97')
    ],

    'Z_bb': [
        Minimal,
        AntiQCD,
        tco.add('Vtype==4',
                'HCSV_mass<100 || HCSV_mass>140',
                addCenJet30e0, 
                naddGoodLeptons10e0, 
                'Jet_btagCSV[hJCidx[1]]>0.8')
    ],

    # W Jets Control Regions
    'W_light': [
        Minimal,
        AntiQCD,
        tco.add('Vtype==2 || Vtype==3',
                'vLeptons_pt>30',
                addCenJet30e0,
                'Jet_btagCSV[hJCidx[0]]<0.97')
    ],

    'W_bb': [
        Minimal,
        AntiQCD,
        tco.add('Vtype==2 || Vtype==3',
                'vLeptons_pt>30',
                addCenJet30e0,
                'Jet_btagCSV[hJCidx[1]]>0.8')
    ],

    # QCD Control Region
    'QCD': [
        Minimal,
        tco.add('MinIf$(abs(TVector2::Phi_mpi_pi(met_phi - Jet_phi)), Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)<0.7',
                'abs(TVector2::Phi_mpi_pi(met_phi - tkMet_phi))>0.7',
                'tkMet_pt<30')
    ]
 
}


# Plotting Categories
CATEGORIES = {

    'Data': {
        'PATH': STEP3_DIR + 'Data_MET.root',
    },

    'ZH': {
        'PATH': STEP3_DIR + 'ZnnH125.root',
    },

    'ggZH': {
        'PATH': STEP3_DIR + 'ggZH125.root',
    },

    'WH': {
        'PATH': STEP3_DIR + 'WlnH125.root',
    },

    'WjLF': {
        'PATH': STEP3_DIR + 'WJets.root',
        'CUTS': LIGHT_FLAVOR,
    },

    'WjHF': {
        'PATH': STEP3_DIR + 'WJets.root',
        'CUTS': HEAVY_FLAVOR,
    },

    'ZjLF': {
        'PATH': STEP3_DIR + 'ZJets.root',
        'CUTS': LIGHT_FLAVOR,
    },

    'ZjHF': {
        'PATH': STEP3_DIR + 'ZJets.root',
        'CUTS': HEAVY_FLAVOR,
    },

    'TT': {
        'PATH': STEP3_DIR + 'TT.root',
    },

    'ST': {
        'PATH': STEP3_DIR + 'ST.root',
    },

    'VV': {
        'PATH': STEP3_DIR + 'VV.root',
    },

    'QCD': {
        'PATH': STEP3_DIR + 'QCD.root',
    },

}


##########################
#-- Plot Configuration --#
##########################

# Output directory for the plots.
PLOT_DIR = 'plots/'

# Target luminosity of the data in inverse picobarns (pb-1).
TARGET_LUMI = 2190

# The weights to apply to the data samples.
DATA_WEIGHT = tco.mult('json', 'HLT_BIT_HLT_PFMET90_PFMHT90_IDTight_v || HLT_BIT_HLT_PFMET170_NoiseCleaned_v')

# The weights to apply to the MC samples.
MC_WEIGHT = tco.mult('sign(genWeight)', TARGET_LUMI, '1./sample_lumi', 'puWeight',
                     'HLT_BIT_HLT_PFMET90_PFMHT90_IDLoose_v || HLT_BIT_HLT_PFMET170_NoiseCleaned_v')

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
        'n_bins': 25,
        'x_min': 0,
        'x_max': 500,
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

