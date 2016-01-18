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
SKIM    : str or Cut
          A cut used to reduce the size of the Step1 ntuples/samples.

SAMPLES : dict of 'Sample Name': {Sample Properties}

Sample Properties
-----------------
path : str
       The path to the Step1 ntuple on EOS. The directory tree from the path
       to the sample's files must contain a single subdirectory at each level,
       as multiple subdirectories may lead to different versions of the sample.
xsec : float, optional
       The MC sample cross-section reported in picobarns (pb).

The cross-sections and k-factors were obtained from the following CERN TWiki's:
twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV#s_13_0_TeV
twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopRefXsec

The branching ratios were obtained from the PDG refernce tables.
"""

SKIM = Preselection & JSON_Triggers

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
PROCESSES : dict of 'Process Name': {Process Properties}

Process Properties
------------------
samples     : list of str
              The list of samples to be combined into a decay process. 
              Sample names must be valid keys of the SAMPLES dictionary.
sample_cuts : list of str or Cut, optional
              A list of sample-specific cuts applied before combination. 
              Their ordering must match that of the samples property.
process_cut : str or Cut, optional
              A process-specific cut applied to all samples before combination.
types       : str
              A string of the following process descriptors, delimited by ':'.
              Data - Contains only data events.
              MC   - Contains only MC events.
              sig  - A signal process.
              bkg  - A background process.
color       : int
              The ROOT color wheel code with which the process will be plotted.
pid         : int, optional
              The unique process id number used by the datacard. Use a positive
              number for background and zero or a negative number for signal.
"""

PROCESSES = {

    'data_obs': {
        'samples': ['Data_MET_C', 'Data_MET_D', 'Data_MET_DP'],
        'types': 'Data',
    },

    'ZH': {
        'samples': ['ZnnH125'],
        'types': 'MC:sig',
        'color': 632,
        'pid': -2,
    },

    'ggZH': {
        'samples': ['ggZH125'],
        'types': 'MC:sig',
        'color': 632 - 7,
        'pid': -1,
    },

    'WH': {
        'samples': ['WlnH125'],
        'types': 'MC:sig',
        'color': 632 + 2,
        'pid': 0,
    },
    
    's_Top': {
        'samples': ['T_s_comb_lep', 'T_t_lep', 'Tbar_t_lep', 'T_tW', 'Tbar_tW'],
        'types': 'MC:bkg',
        'color': 840,
        'pid': 1,
    },

    'TT': {
        'samples': ['TTPow'],
        'types': 'MC:bkg',
        'color': 600,
        'pid': 2,
    },

    'Wj0b': {
        'samples': ['WJetsHT100', 'WJetsHT200', 'WJetsHT400', 'WJetsHT600'],
        'process_cut': Vudsg | Vcc,
        'types': 'MC:bkg',
        'color': 820 - 6,
        'pid': 3,
    },

    'Wj1b': {
        'samples': ['WJetsHT100', 'WJetsHT200', 'WJetsHT400', 'WJetsHT600'],
        'process_cut': Vb,
        'types': 'MC:bkg',
        'color': 820 - 1,
        'pid': 4,
    },

    'Wj2b': {
        'samples': ['WJetsHT100', 'WJetsHT200', 'WJetsHT400', 'WJetsHT600'],
        'process_cut': Vbb,
        'types': 'MC:bkg',
        'color': 820,
        'pid': 5,
    },

    'Zj0b': {
        'samples': ['ZJetsHT100', 'ZJetsHT200', 'ZJetsHT400', 'ZJetsHT600'],
        'process_cut': Vudsg | Vcc,
        'types': 'MC:bkg',
        'color': 400 + 2,
        'pid': 6,
    },

    'Zj1b': {
        'samples': ['ZJetsHT100', 'ZJetsHT200', 'ZJetsHT400', 'ZJetsHT600'],
        'process_cut': Vb,
        'types': 'MC:bkg',
        'color': 400 - 7,
        'pid': 7,
    },

    'Zj2b': {
        'samples': ['ZJetsHT100', 'ZJetsHT200', 'ZJetsHT400', 'ZJetsHT600'],
        'process_cut': Vbb,
        'types': 'MC:bkg',
        'color': 400,
        'pid': 8,
    },

    'QCD': {
        'samples': ['QCDHT100', 'QCDHT200', 'QCDHT300', 'QCDHT500', 'QCDHT700', 'QCDHT1000', 'QCDHT1500', 'QCDHT2000'],
        'types': 'MC:bkg',
        'color': 613,
        'pid': 9,
    },

    'VVLF': {
        'samples': ['WW', 'WZ', 'ZZ'],
        'process_cut': Vudsg | Vcc,
        'types': 'MC:bkg',
        'color': 920,
        'pid': 10,
    },

    'VVHF': {
        'samples': ['WW', 'WZ', 'ZZ'],
        'process_cut': Vb | Vbb,
        'types': 'MC:bkg',
        'color': 920 + 1,
        'pid': 11,
    },

}

##########################
#-- Classification BDT --#
##########################

CLASSIFICATION_SKIM = ''

CLASSIFICATION_VARS = {
    
    'M(jj)': {
        'expression': 'HCSV_reg_mass',
        'title': 'H mass',
        'unit': 'GeV',
        'type': 'F',
    },
    
    'p_T(jj)': {
        'expression': 'HCSV_reg_pt',
        'title': 'H p_{T}',
        'unit': 'GeV',
        'type': 'F',
    },

    'p_T(j1)': {
        'expression': 'max(Jet_pt_reg[hJCidx[0]], Jet_pt_reg[hJCidx[1]])',
        'title': 'H j1 p_{T}',
        'unit': 'GeV',
        'type': 'F',
    },

    'p_T(j2)': {
        'expression': 'min(Jet_pt_reg[hJCidx[0]], Jet_pt_reg[hJCidx[1]])',
        'title': 'H j2 p_{T}',
        'unit': 'GeV',
        'type': 'F',
    },
    
    'p_T(V) (Same as MET for Znn)': {
        'expression': 'met_pt',
        'title': 'Type1 Corr. PF #slash{E}_{T}',
        'unit': 'GeV',
        'type': 'F',
    },

    'H Jet Max CSV': {
        'expression': 'maxCSV := max(0, max(Jet_btagCSV[hJCidx[0]], Jet_btagCSV[hJCidx[1]]))',
        'title': 'CSV_{max}(j1,j2)',
        'unit': '',
        'type': 'F',
    },

    'H Jet Min CSV': {
        'expression': 'minCSV := max(0, min(Jet_btagCSV[hJCidx[0]], Jet_btagCSV[hJCidx[1]]))',
        'title': 'CSV_{min}(j1,j2)',
        'unit': '',
        'type': 'F',
    },

    'dPhi(H,V)': {
        'expression': 'abs(HVdPhi)',
        'title': '#||{#Delta #varphi(H,PF #slash{E}_{T}}',
        'unit': '',
        'type': 'F',
    },

    'dEta(jj)': {
        'expression': 'dEta_jj := abs(Jet_eta[hJCidx[0]] - Jet_eta[hJCidx[1]])',
        'title': '#||{#Delta #eta(j1,j2)}',
        'unit': '',
        'type': 'F',
    },

    'dR(jj)': {
        'expression': 'deltaR_jj',
        'title': '#Delta R(j1,j2)',
        'unit': '',
        'type': 'F',
    },

    'N_aj': {
        'expression': 'naJets_Znn := max(0, Sum$(Jet_pt>25 && abs(Jet_eta)<4.5 && Jet_puId==1)-2)',
        'title': '# Add. Jets {p_{T}>25}',
        'unit': '',
        'type': 'I',
    },


    #deltaPhi(pfMET,J) ADD SUPPORT
    #mindPhiMETJet_dPhi :=  MinIf$(abs(deltaPhi(met_phi,Jet_phi)), Jet_pt>25 && abs(Jet_eta)<4.5 && Jet_puId==1)
    #min #||{#Delta #varphi(pfMET,j25)}
    #''
    #F

    'Add. Jet Max CSV': {
        'expression': 'maxAddCSV := MaxIf$(max(Jet_btagCSV[aJCidx],0), Jet_pt[aJCidx]>20 && abs(Jet_eta[aJCidx])<2.5 && Jet_puId[aJCidx]==1)',
        'title': 'CSV_{max}(Add. CJ 20)',
        'unit': '',
        'type': 'F',
    },

    #mindeltaR(H,aj) ADD SUPPORT
    #mindRAddJetH := Min$(deltaR(Jet_eta[hJCidx], Jet_phi[hJCidx], Jet_eta[aJCidx], Jet_phi[aJCidx]))
    #min #DeltaR(H, add. j25)
    #''
    #F

}

###############
#-- Regions --#
###############

"""
REGION : dict of 'Region Name': list of str or Cut

Region Properties
-----------------
The cuts defining the signal or control region. PyROOT bugs out when
passed a long/complicated cut string. The fix is to split such cuts into
smaller subcuts. The cuts are applied sequentially using the order given.
"""

REGIONS = {
    
    'TT': [
        NoQCD,
        VetoLeptons >= 1,
        (addCenJet30 >= 1) & (Cut('Jet_btagCSV[hJCidx[0]]') > CSV_Medium),
    ],

    'ZLight': [
        NoQCD,
        VetoLeptons == 0,
        (addCenJet30 == 0) & (Cut('Jet_btagCSV[hJCidx[0]]') < CSV_Medium),
    ],

    'Zbb': [
        NoQCD,
        VetoLeptons == 0,
        (addCenJet30 == 0) & (Cut('Jet_btagCSV[hJCidx[0]]') > CSV_Medium) & 'HCSV_reg_mass<110 || HCSV_reg_mass>140',
    ],

    'WLight': [
        NoQCD,
        VetoLeptons == 1,
        (addCenJet30 == 0) & (Cut('Jet_btagCSV[hJCidx[0]]') < CSV_Medium),
    ],

    'Wbb': [
        NoQCD,
        VetoLeptons == 1,
        (addCenJet30 == 0) & (Cut('Jet_btagCSV[hJCidx[0]]') > CSV_Medium),
        Cut('HCSV_reg_pt>170') & 'MinIf$(abs(TVector2::Phi_mpi_pi(Jet_phi-met_phi)),Jet_pt>20 && Jet_eta<5.2)>2',
    ],


    'QCD': [
        FlagsMET & 'HCSV_pt>150',
        'MinIf$(abs(TVector2::Phi_mpi_pi(Jet_phi-met_phi)),Jet_pt>25 && Jet_eta<5.2)<0.5',
    ],
 
}

#############
#-- Plots --#
#############

"""
TARGET_LUMI : int or float
              The target luminosity of the data in inverse picobarns (pb-1).

DATA_WEIGHT : str or Cut
              The plotting weight for data events.

MC_WEIGHT   : str or Cut
              The plotting weight for MC events.

PLOTS       : dict of 'Plot description.': {Plot Properties}

Plot Properties
---------------
name       : str
             The file name for the plot.
expression : str or Cut
             A TTree.Draw expression which will be passed to TTree.Project.
x_title    : str
             The x-axis title for the plot. Supports TLatex syntax.
n_bins     : int
             The number of bins along the x-axis.
x_min      : int or float
             The lower bound of the x-axis.
x_max      : int or float
             The upper bound of the x-axis.
logy       : bool, optional
             Use a logarithmic scale for the y-axis. Default is False.
"""

TARGET_LUMI = 2200

DATA_WEIGHT = ''

# 'puWeight' is currently broken in the Heppy ntuples. Use custom reweighting.
MC_WEIGHT = Cut('sign(genWeight)') * TARGET_LUMI * '1./sample_lumi' * 'weight2(nTrueInt)'

PLOTS = {

    'nPVs': {
        'name': 'nPVs',
        'expression': 'nPVs',
        'x_title': '# of Primary Vertices',
        'n_bins': 60,
        'x_min': 0,
        'x_max': 60,
    },

    'DeltaPhiJetMet15': {
        'name': 'DeltaPhiJetMet15',
        'expression': 'MinIf$(abs(TVector2::Phi_mpi_pi(Jet_phi-met_phi)),Jet_pt>15 && Jet_eta<5.2)',
        'x_title': '#Delta#Phi(Jet,MET)',
        'n_bins': 18,
        'x_min': -0.3,
        'x_max': 3.3,
    },

    'DeltaPhiJetMet20': {
        'name': 'DeltaPhiJetMet20',
        'expression': 'MinIf$(abs(TVector2::Phi_mpi_pi(Jet_phi-met_phi)),Jet_pt>20 && Jet_eta<5.2)',
        'x_title': '#Delta#Phi(Jet,MET)',
        'n_bins': 18,
        'x_min': -0.3,
        'x_max': 3.3,
    },

    'DeltaPhiJetMet25': {
        'name': 'DeltaPhiJetMet25',
        'expression': 'MinIf$(abs(TVector2::Phi_mpi_pi(Jet_phi-met_phi)),Jet_pt>25 && Jet_eta<5.2)',
        'x_title': '#Delta#Phi(Jet,MET)',
        'n_bins': 18,
        'x_min': -0.3,
        'x_max': 3.3,
    },

    'softActivityVH_HT': {
        'name': 'softActivityVH_HT',
        'expression': 'softActivityVH_HT',
        'x_title': 'softActivity HT [GeV]',
        'n_bins': 50,
        'x_min': 0,
        'x_max': 250,
    },

    'softActivityVH_njets2': {
        'name': 'softActivityVH_njets2',
        'expression': 'softActivityVH_njets2',
        'x_title': 'softActivity, # Jets {p_{T} > 2 GeV}',
        'n_bins': 40,
        'x_min': 0,
        'x_max': 40,
    },

    'softActivityVH_njets5': {
        'name': 'softActivityVH_njets5',
        'expression': 'softActivityVH_njets5',
        'x_title': 'softActivity, # Jets {p_{T} > 5 GeV}',
        'n_bins': 20,
        'x_min': 0,
        'x_max': 20,
    },

    'softActivityVH_njets10': {
        'name': 'softActivityVH_njets10',
        'expression': 'softActivityVH_njets10',
        'x_title': 'softActivity, # Jets {p_{T} > 10 GeV}',
        'n_bins': 10,
        'x_min': 0,
        'x_max': 10,
    },

    'Vtype': {
        'name': 'Vtype',
        'expression': 'Vtype',
        'x_title': 'Vtype',
        'n_bins': 6,
        'x_min': 0,
        'x_max': 6,
    },

    'lheHT': {
        'name': 'lheHT',
        'expression': 'Alt$(lheHT,Sum$(Jet_pt))',
        'x_title': 'lheHT',
        'n_bins': 50,
        'x_min': 0,
        'x_max': 1000,
        'logy': True,
    },

    'CSVidx1': {
        'name': 'CSVidx1',
        'expression': 'min(hJCidx[0],hJCidx[1])',
        'x_title': 'Jet1 HCSV Index',
        'n_bins': 20,
        'x_min': 0,
        'x_max': 20,
    },

    'CSVidx2': {
        'name': 'CSVidx2',
        'expression': 'max(hJCidx[0],hJCidx[1])',
        'x_title': 'Jet2 HCSV Index',
        'n_bins': 20,
        'x_min': 0,
        'x_max': 20,
    },

    'MET': {
        'name': 'MET',
        'expression': 'met_pt',
        'x_title': '#slash{E}_{T} [GeV]',
        'n_bins': 50,
        'x_min': 0,
        'x_max': 500,
    },

    'minMETMHT': {
        'name': 'minMETMHT',
        'expression': 'min(met_pt,mhtJet30)',
        'x_title': 'Min(MET,MHT) [GeV]',
        'n_bins': 50,
        'x_min': 0,
        'x_max': 500,
    },

    'addCenJet30': {
        'name': 'addCenJet30',
        'expression': 'Sum$(Jet_pt>30 && abs(Jet_eta)<5.2)-2',
        'x_title': '# Add. Central Jets',
        'n_bins': 10,
        'x_min': 0,
        'x_max': 10,
    },

    'addCenJet30': {
        'name': 'addCenJet30',
        'expression': 'Sum$(Jet_pt>30 && abs(Jet_eta)<5.2)-2',
        'x_title': '# Add. Central Jets',
        'n_bins': 10,
        'x_min': 0,
        'x_max': 10,
    },

    'METpuppi': {
        'name': 'METpuppi',
        'expression': 'metPuppi_pt',
        'x_title': 'MET PUPPI [GeV]',
        'n_bins': 50,
        'x_min': 0,
        'x_max': 500,
    },

    'maxCSV': {
        'name': 'maxCSV',
        'expression': 'TMath::Max(Jet_btagCSV[hJCidx[0]],Jet_btagCSV[hJCidx[1]])',
        'x_title': 'Max CSV',
        'n_bins': 40,
        'x_min': 0,
        'x_max': 1,
    },

    'minCSV': {
        'name': 'minCSV',
        'expression': 'TMath::Min(Jet_btagCSV[hJCidx[0]],Jet_btagCSV[hJCidx[1]])',
        'x_title': 'Min CSV',
        'n_bins': 40,
        'x_min': 0,
        'x_max': 1,
    },

    'sumEt': {
        'name': 'sumEt',
        'expression': 'met_sumEt',
        'x_title': 'sumEt [GeV]',
        'n_bins': 50,
        'x_min': 0,
        'x_max': 2500,
    },

    'MET_sumEt': {
        'name': 'MET_sumEt',
        'expression': 'met_pt/(met_sumEt)',
        'x_title': 'MET/sumEt',
        'n_bins': 50,
        'x_min': 0,
        'x_max': 0.5,
    },

    'MET_sign': {
        'name': 'MET_sign',
        'expression': 'met_pt/sqrt(met_sumEt)',
        'x_title': 'MET/#sqrt{sumEt}',
        'n_bins': 50,
        'x_min': 0,
        'x_max': 20,
    },

    'DeltaPhiJet1Jet2': {
        'name': 'DeltaPhiJet1Jet2',
        'expression': 'abs(TVector2::Phi_mpi_pi(Jet_phi[0]-Jet_phi[1]))',
        'x_title': '#Delta#Phi(Jet1,Jet2)',
        'n_bins': 18,
        'x_min': -0.3,
        'x_max': 3.3,
    },

    'Hmass': {
        'name': 'Hmass',
        'expression': 'HCSV_reg_mass',
        'x_title': 'm(jj) after Regression [GeV]',
        'n_bins': 40,
        'x_min': 0,
        'x_max': 400,
    },

    'Hmass_NoReg': {
        'name': 'Hmass_NoReg',
        'expression': 'HCSV_mass',
        'x_title': 'm(jj) before Regression [GeV]',
        'n_bins': 40,
        'x_min': 0,
        'x_max': 400,
    },

    'Hpt': {
        'name': 'Hpt',
        'expression': 'HCSV_reg_pt',
        'x_title': 'p_{T}(jj) after Regression [GeV]',
        'n_bins': 20,
        'x_min': 0,
        'x_max': 400,
    },

    'Hpt_NoReg': {
        'name': 'Hpt_NoReg',
        'expression': 'HCSV_pt',
        'x_title': 'm(jj) before Regression [GeV]',
        'n_bins': 20,
        'x_min': 0,
        'x_max': 400,
    },

    'HT': {
        'name': 'HT',
        'expression': 'Sum$(Jet_pt * (Jet_pt>30))',
        'x_title': 'HT[GeV]',
        'n_bins': 50,
        'x_min': 0,
        'x_max': 1500,
    },

    'MHT': {
        'name': 'MHT',
        'expression': 'mhtJet30',
        'x_title': 'MHT [GeV]',
        'n_bins': 50,
        'x_min': 0,
        'x_max': 500,
    },

    'tkMET': {
        'name': 'tkMET',
        'expression': 'tkMet_pt',
        'x_title': 'Track MET [GeV]',
        'n_bins': 50,
        'x_min': 0,
        'x_max': 300,
    },

    'tkMetPVch': {
        'name': 'tkMetPVch',
        'expression': 'tkMetPVchs_pt',
        'x_title': 'Track MET CHS [GeV]',
        'n_bins': 50,
        'x_min': 0,
        'x_max': 300,
    },

    'METType1p2': {
        'name': 'METType1p2',
        'expression': 'metType1p2_pt',
        'x_title': 'MET1.2 [GeV]',
        'n_bins': 50,
        'x_min': 0,
        'x_max': 500,
    },

    'DeltaPhiHCSVMET': {
        'name': 'DeltaPhiHCSVMET',
        'expression': 'abs(TVector2::Phi_mpi_pi(HCSV_phi-met_phi))',
        'x_title': '#Delta#Phi(HCSV,MET)',
        'n_bins': 18,
        'x_min': -0.3,
        'x_max': 3.3,
    },

    'addJet': {
        'name': 'addJet',
        'expression': 'Sum$(Jet_pt>30 && abs(Jet_eta)<5.2)-2',
        'x_title': '# Add. Jets',
        'n_bins': 10,
        'x_min': 0,
        'x_max': 10,
    },

    'maxhjPt': {
        'name': 'maxhjPt',
        'expression': 'max(Jet_pt_reg[hJCidx[0]],Jet_pt_reg[hJCidx[1]])',
        'x_title': 'Max H Jet p_{T} after Regression [GeV]',
        'n_bins': 50,
        'x_min': 0,
        'x_max': 400,
    },

    'minhjPt': {
        'name': 'minhjPt',
        'expression': 'min(Jet_pt_reg[hJCidx[0]],Jet_pt_reg[hJCidx[1]])',
        'x_title': 'Min H Jet p_{T} after Regression [GeV]',
        'n_bins': 50,
        'x_min': 0,
        'x_max': 400,
    },

    'maxhjPtOld': {
        'name': 'maxhjPtOld',
        'expression': 'max(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])',
        'x_title': 'Max H Jet p_{T} before Regression [GeV]',
        'n_bins': 50,
        'x_min': 0,
        'x_max': 400,
    },

    'minhjPtOld': {
        'name': 'minhjPtOld',
        'expression': 'min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])',
        'x_title': 'Min H Jet p_{T} before Regression [GeV]',
        'n_bins': 50, 
        'x_min': 0,
        'x_max': 400,
    },

    'naddGoodLeptons5': {
        'name': 'naddGoodLeptons5',
        'expression': 'Sum$(aLeptons_pt>5 && (aLeptons_jetBTagCSV<0.25 || aLeptons_relIso03<0.4 || aLeptons_looseIdSusy!=0 ||  aLeptons_jetDR>0.3 ))+Sum$(vLeptons_pt>5 && (vLeptons_jetBTagCSV<0.25 || vLeptons_relIso03<0.4 || vLeptons_looseIdSusy!=0 ||  vLeptons_jetDR>0.3 ))',
        'x_title': 'naddGoodLeptons5',
        'n_bins': 5,
        'x_min': 0,
        'x_max': 5,
    },
 
    'nTightLeptons': {
        'name': 'nTightLeptons',
        'expression': '(Sum$(aLeptons_eleCutIdSpring15_25ns_v1>3 && aLeptons_pfRelIso04<0.15 && abs(aLeptons_pdgId)==11 && aLeptons_pt>20) + Sum$(vLeptons_eleCutIdSpring15_25ns_v1>3 && vLeptons_pfRelIso04<0.15 && abs(vLeptons_pdgId)==11 && vLeptons_pt>20))+(Sum$(aLeptons_tightId>0 && aLeptons_pfRelIso04<0.15 && abs(aLeptons_pdgId)==13 && aLeptons_pt>20) + Sum$(vLeptons_tightId>0 && vLeptons_pfRelIso04<0.15 && abs(vLeptons_pdgId)==13 && vLeptons_pt>20))',
        'x_title': '# Tight Leptons',
        'n_bins': 5,
        'x_min': 0,
        'x_max': 5,
    },

    'nVetoLeptons': {
        'name': 'nVetoLeptons',
        'expression': 'Sum$(aLeptons_eleCutIdSpring15_25ns_v1>0 && aLeptons_pfRelIso04<0.50 && abs(aLeptons_pdgId)==11 && aLeptons_pt>5) + Sum$(vLeptons_eleCutIdSpring15_25ns_v1>0 && vLeptons_pfRelIso04<0.50 && abs(vLeptons_pdgId)==11 && vLeptons_pt>5)+(Sum$(aLeptons_looseIdPOG>0 && aLeptons_pfRelIso04<0.50 && abs(aLeptons_pdgId)==13 && aLeptons_pt>5) + Sum$(vLeptons_looseIdPOG>0 && vLeptons_pfRelIso04<0.50 && abs(vLeptons_pdgId)==13 && vLeptons_pt>5))',
        'x_title': '# Veto Leptons',
        'n_bins': 5,
        'x_min': 0,
        'x_max': 5,
    },

    'JetPtCloseToMet03': {
        'name': 'JetPtCloseToMet03',
        'expression': 'MaxIf$(Jet_pt,abs(TVector2::Phi_mpi_pi(Jet_phi-met_phi))<0.3)',
        'x_title': 'Jet p_{T} {#Delta#Phi(Jet,MET)<0.3} [GeV]',
        'n_bins': 18,
        'x_min': 0,
        'x_max': 60,
    },

    'JetPtCloseToMet05': {
        'name': 'JetPtCloseToMet05',
        'expression': 'MaxIf$(Jet_pt,abs(TVector2::Phi_mpi_pi(Jet_phi-met_phi))<0.5)',
        'x_title': 'Jet p_{T} {#Delta#Phi(Jet,MET)<0.5} [GeV]',
        'n_bins': 18,
        'x_min': 0,
        'x_max': 60,
    },

    'JetPtCloseToMet07': {
        'name': 'JetPtCloseToMet07',
        'expression': 'MaxIf$(Jet_pt,abs(TVector2::Phi_mpi_pi(Jet_phi-met_phi))<0.7)',
        'x_title': 'Jet p_{T} {#Delta#Phi(Jet,MET)<0.7} [GeV]',
        'n_bins': 18,
        'x_min': 0,
        'x_max': 60,
    },

    'JetPtCloseToMet10': {
        'name': 'JetPtCloseToMet10',
        'expression': 'MaxIf$(Jet_pt,abs(TVector2::Phi_mpi_pi(Jet_phi-met_phi))<1.0)',
        'x_title': 'Jet p_{T} {#Delta#Phi(Jet,MET)<1.0} [GeV]',
        'n_bins': 18,
        'x_min': 0,
        'x_max': 60,
    },

    'CSV3': {
        'name': 'CSV3',
        'expression': 'Max$(Jet_btagCSV[aJCidx])',
        'x_title': 'CSV3',
        'n_bins': 40,
        'x_min': 0,
        'x_max': 1,
    },

    'HmassFSR_NoReg': {
        'name': 'HmassFSR_NoReg',
        'expression': 'HaddJetsdR08_mass',
        'x_title': 'm(jj) including FSR before Regression [GeV]',
        'n_bins': 40,
        'x_min': 0,
        'x_max': 400,
    },

    'DeltaPhiMETMHT': {
        'name': 'DeltaPhiMETMHT',
        'expression': 'abs(TVector2::Phi_mpi_pi(mhtPhiJet30-met_phi))',
        'x_title': '#Delta#Phi(MHT,MET)',
        'n_bins': 40,
        'x_min': 0,
        'x_max': 2,
    },

    'DeltaPhiMETtkMETPVchs': {
        'name': 'DeltaPhiMETtkMETPVchs',
        'expression': 'abs(TVector2::Phi_mpi_pi(tkMetPVchs_phi-met_phi))',
        'x_title': '#Delta#Phi(tkMetPVchs,MET)',
        'n_bins': 40,
        'x_min': 0,
        'x_max': 2,
    },

    'DeltaPhiMETpuppiMET': {
        'name': 'DeltaPhiMETpuppiMET',
        'expression': 'abs(TVector2::Phi_mpi_pi(metPuppi_phi-met_phi))',
        'x_title': '#Delta#Phi(PUPPI MET,MET)',
        'n_bins': 40,
        'x_min': 0,
        'x_max': 2,
    },

    'DeltaPhiMETtkMET': {
        'name': 'DeltaPhiMETtkMET',
        'expression': 'abs(TVector2::Phi_mpi_pi(tkMet_phi-met_phi))',
        'x_title': '#Delta#Phi(tkMet,MET)',
        'n_bins': 40,
        'x_min': 0,
        'x_max': 2,
    },

}

