import glob
import logging
import multiprocessing as mp
import os
import subprocess as sp
import tempfile as tf

import numpy as np
import ROOT

from process import PROCESS_DIR
from settings import WORK_DIR, PROCESSES, CLASSIFICATION, VARIABLES


# Output Directories
DATASET_DIR = WORK_DIR + 'classification/datasets/'
MONITOR_DIR = WORK_DIR + 'classification/monitoring/'
WEIGHT_DIR = WORK_DIR + 'classification/weights/'

class Classification(object):

    def __init__(self, 
                 job_name = 'TMVAClassification',
                 dataset = 'test',
                 preselection = '',
                 factory = 'Silent:Transformations=I:AnalysisType=Classification',
                 model_name = 'BDT',
                 hyperparams = None,
                 **kwargs):

        self.logger = logging.getLogger('Classification')
        self.logger.info('Initialized for {}'.format(job_name))

        self.job_name = job_name
        self.dataset = dataset
        self.preselection = preselection
        self.factory = factory
        self.model_name = model_name
        if hyperparams is None:
            self.hyperparams = {
                'NTrees': 400,
                'MaxDepth': 3,
                'MinNodeSize': 0.05,
                'nCuts': 35,
                'BoostType': 'AdaBoost',
                'AdaBoostBeta': 0.5,
                'SeparationType': 'GiniIndex',
            }
        else:
            self.hyperparams = hyperparams

    def run(self):

        if not os.path.isfile(DATASET_DIR + self.dataset + '.root'):
            self._make_dataset()

        self._train_classifier()
    

    def _make_dataset(self):

        self.logger.info('Generating dataset "{}"'.format(self.dataset))

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
        results = mp.Queue()

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
            mp.Process(target = self._split_process, args = (tasks, results))
            for i in xrange(min(n_tasks, mp.cpu_count()))
        ]

        for p in _processes:
            tasks.put(None)
            p.start()

        for p in _processes:
            p.join()

        results.put(None)

        for r in iter(results.get, None):
            self.logger.info('{!s} entries in {!s} passed preselection and were split into {!s}({!s}) for training(testing).'.format(*r))

        # hadd Files
        inputfiles = glob.glob(self.tmpdir + '*.root')
        outputfile = DATASET_DIR + self.dataset + '.root'

        sp.check_call(['hadd', '-f', outputfile] + inputfiles)
        sp.check_call(['rm', '-r', self.tmpdir])
         
    def _split_process(self, tasks = None, results = None, rnd_seed = 0):
    
        for process, types in iter(tasks.get, None):

            infile = ROOT.TFile(PROCESS_DIR + process + '.root', 'read')
            outfile = ROOT.TFile(self.tmpdir + process + '.root', 'recreate')

            intree = infile.Get('tree')
            intree.SetName(process)

            # Apply preselection cuts and shuffle the entries which pass.
            for i, cut in enumerate(self.preselection):
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

            result = (n_entries, process, train_tree.GetEntriesFast(), test_tree.GetEntriesFast())

            outfile.Close()
            infile.Close()

            results.put(result)

    def _get_hyperparams(self):

        hyperparams = {}

        for name, value in self.hyperparams.iteritems():
            if hasattr(value, 'rvs'):
                hyperparams[name] = value.rvs()
            elif isinstance(value, list):
                hyperparams[name] = np.random.choice(value)
            else:
                hyperparams[name] = value

        return hyperparams

    def _train_classifier(self):

        self.logger.info('Training classifier...')

        ROOT.TMVA.Tools.Instance()

        # Change the weight directory path.
        ioNames = ROOT.TMVA.gConfig().GetIONames()
        ioNames.fWeightFileDir = WEIGHT_DIR
        #ioNames.fWeightFileExtension = 'foo' # e.g. TMVA_BDT.foo.xml

        # from datetime import datetime, include timestamp in file name.
        #timestamp = datetime.now().strftime('%d%b%Y_%H%M%S')
        outfile = ROOT.TFile(MONITOR_DIR + self.job_name + '.root', 'recreate')

        factory = ROOT.TMVA.Factory(self.job_name, outfile, self.factory)

        # Register the classification variables.
        for var in VARIABLES.itervalues():
            factory.AddVariable(var['expression'], var['title'], var['unit'], var['type'])

        # Load the dataset and add the training and validation samples.
        infile = ROOT.TFile(DATASET_DIR + self.dataset + '.root', 'read')

        for key in infile.GetListOfKeys():
            name = key.GetName()
            tree_class, tree_type = name.split('_')[-2:]
            if tree_class == 'sig':
                factory.AddSignalTree(key.ReadObj(), 1.0, tree_type)
            elif tree_class == 'bkg':
                factory.AddBackgroundTree(key.ReadObj(), 1.0, tree_type)
                    
        factory.SetSignalWeightExpression('sample_lumi')
        sigCut = ROOT.TCut('')

        factory.SetBackgroundWeightExpression('sample_lumi')
        bkgCut = ROOT.TCut('')

        factory.PrepareTrainingAndTestTree(sigCut, bkgCut, 'NormMode=None')

        hyperparams = ':'.join('{!s}={!s}'.format(*x) for x in _get_hyperparams.iteritems())
        factory.BookMethod(ROOT.TMVA.Types.kBDT, self.model_name, '!H:!V' + hyperparams)

        factory.TrainAllMethods()
        factory.TestAllMethods()
        factory.EvaluateAllMethods()

#------
# Main
#------

if __name__ == '__main__':

    # Set ROOT to batch mode.
    ROOT.gROOT.SetBatch(1)

    logging.basicConfig(level = logging.INFO,
                        format = '%(name)s(%(levelname)s) - %(message)s')

    classification = Classification(**CLASSIFICATION)
    print classification._get_hyperparams()
    print classification._get_hyperparams()
    print classification._get_hyperparams()
    classification.run()

