import logging
import os

import numpy as np
import ROOT

from process import PROCESS_DIR
import settings


##########################
#-- Classification BDT --#
##########################

VARIABLES = {
    
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


DATASET_DIR = settings.WORK_DIR + 'classification/datasets/'

WEIGHT_DIR = settings.WORK_DIR + 'classification/weights/'

def make_dataset(process = '', rnd_seed = 0):

    try:
        os.makedirs(DATASET_DIR)
    except OSError:
        if not os.path.isdir(DATASET_DIR):
            raise

    infile = ROOT.TFile(PROCESS_DIR + process + '.root', 'read')
    intree = infile.Get('tree')

    # Figure out how to split the entries.
    # For uneven splits, the training set gets the remaining entries.
    n_entries = intree.GetEntriesFast()
    n_train = n_entries / 2
    n_train += n_entries % 2
 
    # Randomly shuffle the entries.
    # Seeding for now to make sure code works.
    np.random.seed(rnd_seed)
    permute_entries = np.random.permutation(n_entries)

    train_elist = ROOT.TEventList('train')
    for i in permute_entries[:n_train]:
        train_elist.Enter(i)

    test_elist  = ROOT.TEventList('test')
    for i in permute_entries[n_train:]:
        test_elist.Enter(i)

    # Start filling the new file.
    outfile = ROOT.TFile(DATASET_DIR + process + '.root', 'recreate')

    intree.SetName('train')
    intree.SetEventList(train_elist)
    train_tree = intree.CopyTree('')
    train_tree.Write()

    intree.SetName('test')
    intree.SetEventList(test_elist)
    test_tree = intree.CopyTree('')
    test_tree.Write()

    outfile.Close()
    infile.Close()


def run_classification():

    # Book Trees
    signal = ROOT.TFile.Open(DATASET_DIR + 'ZH.root', 'read')
    #sig_train_tree = signal.Get('train')
    #sig_test_tree = signal.Get('test')
    background = ROOT.TFile.Open(DATASET_DIR + 'QCD.root', 'read')
    #bkg_train_tree = background.Get('train')
    #bkg_test_tree = background.Get('test')

    # Start TMVA
    ROOT.TMVA.Tools.Instance()

    # Change the weight directory.
    ioNames = ROOT.TMVA.gConfig().GetIONames()
    ioNames.fWeightFileDir = WEIGHT_DIR
    ioNames.fWeightFileExtension = 'ZnnHbb'

    outfile = ROOT.TFile('TMVA_BDT.root', 'recreate')

    factory = ROOT.TMVA.Factory('TMVAClassification',
                                outfile,
                                #'!V:!Silent:Color:!DrawProgressBar:Transformations=I:AnalysisType=Classification'
                                '!Silent:Transformations=I:AnalysisType=Classification')

    for var in VARIABLES.itervalues():
        factory.AddVariable(var['expression'], var['title'], var['unit'], var['type']) 

    factory.AddSignalTree(signal.Get('train'), 2.0, 'train')
    factory.AddSignalTree(signal.Get('test'), 2.0, 'test')

    factory.AddBackgroundTree(background.Get('train'), 2.0, 'train')
    factory.AddBackgroundTree(background.Get('test'), 2.0, 'test')

    factory.SetSignalWeightExpression('sample_lumi')
    factory.SetBackgroundWeightExpression('sample_lumi')

    sigCut = ROOT.TCut('')
    bkgCut = ROOT.TCut('')

    factory.PrepareTrainingAndTestTree(sigCut, bkgCut, 'nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=None')

    bdt_options = '!H:V:NTrees=400:MinNodeSize=0.05:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=35:PruneMethod=NoPruning:!CreateMVAPdfs:!DoBoostMonitor'
    factory.BookMethod(ROOT.TMVA.Types.kBDT, 'BDT', bdt_options) 

    factory.TrainAllMethods()
    factory.TestAllMethods()
    factory.EvaluateAllMethods()

if __name__ == '__main__':

    run_classification()






"""
BEGIN CODE FOR RANDOM SEARCH HYPERPARAM OPTIMIZATION

import random

import numpy as np
from scipy import stats


# Dictionary of hyperparameters and their distributions.
# http://tmva.sourceforge.net/optionRef.html#MVA::BDT
# Either a list or sensible distribution should be provided.
hyperparameters = {
    'NTrees': stats.randint(200, 1001),
    'MaxDepth': stats.randint(3, 11),
    'nCuts': stats.randint(10, 31),
    'SeparationType': ['CrossEntropy', 'GiniIndex', 'GiniIndexWithLaplace', 'MisClassificationError'],
    'Shrinkage': stats.expon(loc = 0.001, scale = 0.1),
}

def rand_hyperparam():

    configuration = {}

    for option, value in hyperparameters.iteritems():
        if hasattr(value, 'rvs'):
            configuration[option] = value.rvs()
        else:
            configuration[option] = random.choice(value)

    return configuration

if __name__ == '__main__':

    # Seed the rng for repoducibility.
    np.random.seed(0)

    trials = mp.Queue()
    n_trials = 30

    for i in range(n_trials):
        trials.put(rand_hyperparam())

    trials.put(None)

    for trial in iter(trials.get, None):
        print ':'.join('{}={}'.format(*optval) for optval in trial.iteritems())
"""