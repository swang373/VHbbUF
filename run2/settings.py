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

################################
# Step1 Ntuples and Properties #
################################

SAMPLES = {
    
    # Datasets
    'Data_MET_C': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V15/MET/VHBB_HEPPY_V15_MET__Run2015C_25ns-05Oct2015-v1',
    },

    'Data_MET_D': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V15/MET/VHBB_HEPPY_V15_MET__Run2015D-05Oct2015-v1',
    },

    'Data_MET_DP': {
        'EOS_DIR': '/store/group/phys_higgs/hbb/ntuples/V15/MET/VHBB_HEPPY_V15_MET__Run2015D-PromptReco-v4',
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

####################
# Step2 Parameters #
####################

STEP2_DIR = '/afs/cern.ch/work/s/swang373/private/V14/'
STEP2_CUT = 'Vtype>=0 && met_pt>150'


##############################
# Control Region Definitions #
##############################

CONTROL_REGION_DIR = '/afs/cern.ch/work/s/swang373/private/V14/CR/'

# Custom Cuts

minimal = tco.add('Vtype==2 || Vtype==3 || Vtype==4', 'HCSV_pt>150', 'met_pt>150', 'Jet_btagCSV[hJCidx[1]]>0.3',
                  'min(Jet_pt[hJCidx[0]], Jet_pt[hJCidx[1]])>30', 'HCSV_mass<300')

antiQCD = tco.add('HCSV_pt>150', 'met_pt>150', 'tkMet_pt>30', 'Jet_btagCSV[hJCidx[1]]>0.3', 'min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])>30',
                  'MinIf$(abs(TVector2::Phi_mpi_pi(met_phi - Jet_phi)), Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)>0.7',
                  'abs(TVector2::Phi_mpi_pi(met_phi - tkMet_phi))<0.7')

addCenJet30m0 = '(Sum$(Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)-2)>0'

addCenJet30e0 = '(Sum$(Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)-2)==0'

addCenJet30e1 = '(Sum$(Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)-2)<=1'

naddGoodLeptons10e0 = '(Sum$(aLeptons_pt>10 && (aLeptons_jetBTagCSV<0.25 || aLeptons_relIso03<0.4 || aLeptons_looseIdSusy!=0 || aLeptons_jetDR>0.3)) + Sum$(vLeptons_pt>10 && (vLeptons_jetBTagCSV<0.25 || vLeptons_relIso03<0.4 || vLeptons_looseIdSusy!=0 || vLeptons_jetDR>0.3)))==0'

naddGoodTaus20e0 = 'Sum$(TauGood_idDecayMode>=1 && TauGood_idCI3hit>=1 && TauGood_pt>20 && abs(TauGood_eta)<2.3)==0'

# Jet Flavors

LIGHT_FLAVOR = 'abs(Jet_mcFlavour[hJCidx[0]])!=5 && abs(Jet_mcFlavour[hJCidx[1]])!=5'
HEAVY_FLAVOR = 'abs(Jet_mcFlavour[hJCidx[0]])==5 || abs(Jet_mcFlavour[hJCidx[1]])==5'

# Control Regions

CONTROL_REGIONS = {

    # Signal Region
    'signal_loose': [
        tco.add('Vtype==4', 'min(Jet_pt[hJCidx[0]], Jet_pt[hJCidx[1]])>30', 'Jet_btagCSV[hJCidx[1]]>0.605',
                'HCSV_mass<100 || HCSV_mass>140', naddGoodLeptons10e0, naddGoodTaus20e0),
        antiQCD,
    ],

    'signal_tight': [
        tco.add('Vtype==4', 'min(Jet_pt[hJCidx[0]], Jet_pt[hJCidx[1]])>30', 'Jet_btagCSV[hJCidx[1]]>0.8',
                naddGoodLeptons10e0, naddGoodTaus20e0, addCenJet30e1),
        antiQCD,
    ],

    # TTbar Control Region
    'TTbar': [
        tco.add('Vtype==2 || Vtype==3', 'vLeptons_pt>30', addCenJet30m0, 'Jet_btagCSV[hJCidx[0]]>0.97', 'Jet_btagCSV[hJCidx[1]]<0.97'),
        antiQCD,
    ],

    # Z Jets Control Regions
    'Z_light': [
        tco.add('Vtype==4', addCenJet30e0, naddGoodLeptons10e0, 'Jet_btagCSV[hJCidx[0]]<0.97'),
        antiQCD,
    ],

    'Z_bb': [
        tco.add('Vtype==4', 'HCSV_mass<100 || HCSV_mass>140', addCenJet30e0, naddGoodLeptons10e0, 'Jet_btagCSV[hJCidx[1]]>0.8'),
        antiQCD,
    ],

    # W Jets Control Regions
    'W_light': [
        tco.add('Vtype==2 || Vtype==3', 'vLeptons_pt>30', addCenJet30e0, 'Jet_btagCSV[hJCidx[0]]<0.97'),
        antiQCD,
    ],

    'W_bb': [
        tco.add('Vtype==2 || Vtype==3', 'vLeptons_pt>30', addCenJet30e0, 'Jet_btagCSV[hJCidx[1]]>0.8'),
        antiQCD,
    ],

    # QCD Control Region
    'QCD': [
        tco.add('HCSV_pt>150', 'met_pt>150', 'tkMet_pt<30', 'Jet_btagCSV[hJCidx[1]]>0.3', 'min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])>30',
              'MinIf$(abs(TVector2::Phi_mpi_pi(met_phi - Jet_phi)), Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)<0.7', 
              'abs(TVector2::Phi_mpi_pi(met_phi - tkMet_phi))>0.7'),
    ],
}


#######################
# Plotting Categories #
#######################

CATEGORIES = {

    'Data': {
        'PATH': STEP2_DIR + 'Data_MET.root',
    },

    'ZH': {
        'PATH': STEP2_DIR + 'ZnnH125.root',
    },

    'ggZH': {
        'PATH': STEP2_DIR + 'ggZH125.root',
    },

    'WH': {
        'PATH': STEP2_DIR + 'WlnH125.root',
    },

    'WjLF': {
        'PATH': STEP2_DIR + 'WJets.root',
        'CUTS': LIGHT_FLAVOR,
    },

    'WjHF': {
        'PATH': STEP2_DIR + 'WJets.root',
        'CUTS': HEAVY_FLAVOR,
    },

    'ZjLF': {
        'PATH': STEP2_DIR + 'ZJets.root',
        'CUTS': LIGHT_FLAVOR,
    },

    'ZjHF': {
        'PATH': STEP2_DIR + 'ZJets.root',
        'CUTS': HEAVY_FLAVOR,
    },

    'TT': {
        'PATH': STEP2_DIR + 'TTPow.root',
    },

    'ST': {
        'PATH': STEP2_DIR + 's_Top.root',
    },

    'VV': {
        'PATH': STEP2_DIR + 'VV.root',
    },

    'QCD': {
        'PATH': STEP2_DIR + 'QCD.root',
    },

}


#######################
# Plotting Parameters #
#######################

TARGET_LUMI = 1280

DATA_WEIGHT = tco.mult('json', 'HLT_BIT_HLT_PFMET90_PFMHT90_IDTight_v')

MC_WEIGHT = tco.mult('sign(genWeight)', TARGET_LUMI, '1./sample_lumi', 'HLT_BIT_HLT_PFMET90_PFMHT90_IDLoose_v')






















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
'ZZ': '/store/group/phys_higgs/hbb/ntuples/V14/ZZ_TuneCUETP8M1_13TeV-pythia8/VHBB_HEPPY_V14_ZZ_TuneCUETP8M1_13TeV-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151025_083230',
}

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
minimal = tco.add('Vtype==2 || Vtype==3 || Vtype==4', 
                  'HCSV_pt>150', 
                  'met_pt>150', 
                  'Jet_btagCSV[hJCidx[1]]>0.3', 
                  'min(Jet_pt[hJCidx[0]], Jet_pt[hJCidx[1]])>30', 
                  'HCSV_mass<300')

antiQCD = tco.add('HCSV_pt>150',
                  'met_pt>150',
                  'MinIf$(abs(TVector2::Phi_mpi_pi(met_phi - Jet_phi)), Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)>0.7',
                  'abs(TVector2::Phi_mpi_pi(met_phi - tkMet_phi))<0.7',
                  'tkMet_pt>30',
                  'Jet_btagCSV[hJCidx[1]]>0.3',
                  'min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])>30')

addCenJet30m0 = '(Sum$(Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)-2)>0'
addCenJet30e0 = '(Sum$(Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)-2)==0'
addCenJet30e1 = '(Sum$(Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)-2)<=1'

naddGoodLeptons10e0 = '(Sum$(aLeptons_pt>10 && (aLeptons_jetBTagCSV<0.25 || aLeptons_relIso03<0.4 || aLeptons_looseIdSusy!=0 || aLeptons_jetDR>0.3)) + Sum$(vLeptons_pt>10 && (vLeptons_jetBTagCSV<0.25 || vLeptons_relIso03<0.4 || vLeptons_looseIdSusy!=0 || vLeptons_jetDR>0.3)))==0'

naddGoodTaus20e0 = 'Sum$(TauGood_idDecayMode>=1 && TauGood_idCI3hit>=1 && TauGood_pt>20 && abs(TauGood_eta)<2.3)==0'

heavy_flavour = 'abs(Jet_mcFlavour[hJCidx[0]])==5 || abs(Jet_mcFlavour[hJCidx[1]])==5'
light_flavour = 'abs(Jet_mcFlavour[hJCidx[0]])!=5 && abs(Jet_mcFlavour[hJCidx[1]])!=5'

# Signal Region
signal_loose = tco.add('Vtype==4',
                       'min(Jet_pt[hJCidx[0]], Jet_pt[hJCidx[1]])>30',
                       'Jet_btagCSV[hJCidx[1]]>0.605',
                       'HCSV_mass<100 || HCSV_mass>140',
                       naddGoodLeptons10e0,
                       naddGoodTaus20e0)

signal_tight = tco.add('Vtype==4',
                       'min(Jet_pt[hJCidx[0]], Jet_pt[hJCidx[1]])>30',
                       'Jet_btagCSV[hJCidx[1]]>0.8',
                       naddGoodLeptons10e0,
                       naddGoodTaus20e0,
                       addCenJet30e1)

# TTbar Control Region
TTbar = tco.add('Vtype==2 || Vtype==3',
                'vLeptons_pt>30',
                addCenJet30m0,
                'Jet_btagCSV[hJCidx[0]]>0.97',
                'Jet_btagCSV[hJCidx[1]]<0.97')

# Z+Jets Control Regions
Zlight = tco.add('Vtype==4',
                 addCenJet30e0,
                 naddGoodLeptons10e0,
                 'Jet_btagCSV[hJCidx[0]]<0.97')

Zbb = tco.add('Vtype==4',
              'HCSV_mass<100 || HCSV_mass>140',
              addCenJet30e0,
              naddGoodLeptons10e0,
              'Jet_btagCSV[hJCidx[1]]>0.8')

# W+Jets Control Regions
Wlight = tco.add('Vtype==2 || Vtype==3',
                 'vLeptons_pt>30',
                 addCenJet30e0,
                 'Jet_btagCSV[hJCidx[0]]<0.97')

Wbb = tco.add('Vtype==2 || Vtype==3',
              'vLeptons_pt>30',
              addCenJet30e0,
              'Jet_btagCSV[hJCidx[1]]>0.8')

# QCD Control Region
QCD = tco.add('HCSV_pt>150', 
              'met_pt>150', 
              'MinIf$(abs(TVector2::Phi_mpi_pi(met_phi - Jet_phi)), Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)<0.7', 
              'abs(TVector2::Phi_mpi_pi(met_phi - tkMet_phi))>0.7', 
              'Jet_btagCSV[hJCidx[1]]>0.3', 
              'min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])>30',
              'tkMet_pt<30')

# Data Event Weight
data_weight = tco.mult('json', 
                       'HLT_BIT_HLT_PFMET90_PFMHT90_IDTight_v')

# MC Event Weight
mc_weight = tco.mult('sign(genWeight)', 
                     target_lumi, 
                     '1./sample_lumi',
                     'HLT_BIT_HLT_PFMET90_PFMHT90_IDLoose_v')
"""
