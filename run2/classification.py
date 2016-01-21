import logging
import multiprocessing as mp
import os

import numpy as np
import ROOT

from dataset import DATASET_DIR, Dataset
from settings import WORK_DIR, CLASSIFICATION, VARIABLES


# Output Directories
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
                 n_trials = 1,
                 **kwargs):

        self.logger = logging.getLogger('Classification')
        self.logger.info('Initialized for {}'.format(job_name))

        self.job_name = job_name

        self.dataset = dataset
        self.preselection = preselection

        self.factory = factory

        self.model_name = model_name

        if hyperparams is None:
            self.logger.info('Using default hyperparameter set.')
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

        self.n_trials = n_trials
        if n_trials > 1:
            self.logger.info('{!s} random search trials will be conducted.'.format(n_trials))
        
    def run(self):

        if not os.path.isfile(DATASET_DIR + self.dataset + '.root'):
            Dataset(self.dataset).make()

        # Parallel Train 
        tasks = mp.Queue()
        results = mp.Queue()

        if self.n_trials == 1:
            tasks.put(('', self._sample_hyperparams()))
        elif self.n_trials > 1:
            for trial_id in ('_trial{!s}'.format(i) for i in xrange(self.n_trials)):
                tasks.put((trial_id, self._sample_hyperparams()))

        _processes = [
            mp.Process(target = self._train_classifier, args = (tasks, results))
            for cpu in xrange(min(self.n_trials, mp.cpu_count()))
        ]

        for p in _processes:
            tasks.put(None)
            p.start()

        for p in _processes:
            p.join()

        results.put(None)

        for r in iter(results.get, None):
            self.logger.info('Model {!s} trained using options "{!s}".'.format(*r))

    def apply(self, model_name = ''):
    
        pass
        #reader = ROOT.TMVA.Reader()

        #variables = {}

        #for var in VARIABLES.itervalues():

    def _sample_hyperparams(self):
        sampled_params = {}
        for name, value in self.hyperparams.iteritems():
            if hasattr(value, 'rvs'):
                sampled_params[name] = value.rvs()
            elif isinstance(value, list):
                sampled_params[name] = value[np.random.randint(len(value))]
            else:
                sampled_params[name] = value
        return sampled_params

    def _train_classifier(self, tasks = None, results = None):

        for trial_id, hyperparams in iter(tasks.get, None):

            ROOT.TMVA.Tools.Instance()

            model_name = self.model_name + trial_id

	    # Change the weight directory path. TMVA makes this directory.
	    ioNames = ROOT.TMVA.gConfig().GetIONames()
	    ioNames.fWeightFileDir = WEIGHT_DIR
	    #ioNames.fWeightFileExtension = 'foo' # job_name_model_name_*.foo.xml
	    
	    # Output Directory
	    try:
                os.makedirs(MONITOR_DIR)
            except OSError:
                if not os.path.isdir(MONITOR_DIR):
                    raise

            # from datetime import datetime, include timestamp in file name.
            #timestamp = datetime.now().strftime('%d%b%Y_%H%M%S')
            outname = MONITOR_DIR + '{!s}_{!s}.root'.format(self.job_name, model_name)
            outfile = ROOT.TFile(outname, 'recreate')

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

            model_options = ':'.join('{!s}={!s}'.format(*x) for x in hyperparams.iteritems())
            factory.BookMethod(ROOT.TMVA.Types.kBDT, model_name, '!H:!V:' + model_options)

            factory.TrainAllMethods()
            factory.TestAllMethods()
            factory.EvaluateAllMethods()

            result = (model_name, model_options)
            results.put(result)

#------
# Main
#------

if __name__ == '__main__':

    # Set ROOT to batch mode.
    ROOT.gROOT.SetBatch(1)

    logging.basicConfig(level = logging.INFO,
                        format = '%(name)s(%(levelname)s) - %(message)s')

    classification = Classification(**CLASSIFICATION)
    classification.run()

