import glob
import logging
import multiprocessing as mp
import os
import subprocess as sp
import tempfile as tf

import numpy as np
import ROOT

from process import PROCESS_DIR
from settings import WORK_DIR, PROCESSES, CLASSIFICATION_SKIM


# Output Directories
DATASET_DIR = WORK_DIR + 'classification/datasets/'
WEIGHT_DIR = WORK_DIR + 'classification/weights/'

class Classification(object):

    def __init__(self):
        self.logger = logging.getLogger('Classification')

    def _make_dataset(self, rnd_seed = 0):

        self.logger.info('Generating dataset...')

        # Output Directory
        try:
            os.makedirs(DATASET_DIR)
        except OSError:
            if not os.path.isdir(DATASET_DIR):
                raise

        # Temporary Work Directory
        tmpdir = tf.mkdtemp(dir = DATASET_DIR)
        self.tmpdir = tmpdir + '/'

        # Parallel Split
        n_tasks = 0
        tasks = mp.Queue()

        processes = []

        for process, properties in PROCESSES.iteritems():
            types = set(properties['types'].lower().split(':'))
            if 'data' in types:
                continue
            if 'sig' in types:
                tasks.put((process, 'sig'))
            elif 'bkg' in types:
                tasks.put((process, 'bkg'))
            n_tasks += 1

        _processes = [
            mp.Process(target = self._split_process, args = (tasks,))
            for i in xrange(min(n_tasks, mp.cpu_count()))
        ]

        for p in _processes:
            tasks.put(None)
            p.start()

        for p in _processes:
            p.join()

        # hadd Files
        inputfiles = glob.glob(self.tmpdir + '*.root')
        outputfile = DATASET_DIR + 'dataset.root'

        sp.check_call(['hadd', '-f', outputfile] + inputfiles)
        sp.check_call(['rm', '-r', self.tmpdir])
         
    def _split_process(self, tasks = None, rnd_seed = 0):
    
        for process, types in iter(tasks.get, None):

            infile = ROOT.TFile(PROCESS_DIR + process + '.root', 'read')
            outfile = ROOT.TFile(self.tmpdir + process + '.root', 'recreate')

            intree = infile.Get('tree')
            intree.SetName(process)

            # Apply preselection cuts and shuffle the entries which pass.
            for i, cut in enumerate(CLASSIFICATION_SKIM):
                n_entries = intree.Draw('>>{0!s}_skim_{1!s}'.format(process, i), cut)
                eventlist = ROOT.gDirectory.Get('{0!s}_skim_{1!s}'.format(process, i))
                intree.SetEventList(eventlist)

            entries = np.zeros(n_entries, dtype = np.int64)
            for i in xrange(n_entries):
                entries[i] = eventlist.GetEntry(i)

            np.random.seed(rnd_seed)
            np.random.shuffle(entries)
 
            # Split the entries into a training and test set.
            n_test = n_entries / 2
            n_train = n_test if (n_entries % 2 == 0) else n_test + 1

            train_elist = ROOT.TEventList(process + '_train', '', n_train)
            for entry in np.sort(entries[:n_train]):
                train_elist.Enter(entry)

            test_elist = ROOT.TEventList(process + '_test', '', n_test)
            for entry in np.sort(entries[n_train:]):
                test_elist.Enter(entry)

            # Create the training and test set trees.
            intree.SetEventList(train_elist)
            train_tree = intree.CopyTree('')
            train_tree.SetName('{}_{}_train'.format(process, types))
            train_tree.Write()

            intree.SetEventList(test_elist)
            test_tree = intree.CopyTree('')
            test_tree.SetName('{}_{}_test'.format(process, types))
            test_tree.Write()

            # Remove any autosaved TTrees.
            outfile.Delete(process + ';*')

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

    # Set ROOT to batch mode.
    ROOT.gROOT.SetBatch(1)

    logging.basicConfig(level = logging.INFO,
                        format = '%(name)s(%(levelname)s) - %(message)s')

    classifier = Classification()
    classifier._make_dataset()






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
