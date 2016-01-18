import logging
import os

import numpy as np
import ROOT

from process import PROCESS_DIR
import settings

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
    n_test = n_entries / 2
    n_train = n_test if (n_entries % 2 == 0) else n_test + 1
    
    # Randomly shuffle the entries.
    # Seeding for now to make sure code works.
    # The entry numbers must be fed into the eventlist
    # in ASCENDING order, because of its binary search structure.
    np.random.seed(rnd_seed)
    permute_entries = np.random.permutation(n_entries)

    train_elist = ROOT.TEventList('train', '', n_train)
    for i in np.sort(permute_entries[:n_train]): 
        train_elist.Enter(i)

    test_elist  = ROOT.TEventList('test', '', n_test)
    for i in np.sort(permute_entries[n_train:]):
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
