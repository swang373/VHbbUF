import logging
import os

import ROOT

from dataset import DATASET_DIR, Dataset 
from settings import WORK_DIR, VARIABLES

# Output Directories
MONITOR_DIR = WORK_DIR + 'classification/monitoring/'
WEIGHT_DIR = WORK_DIR + 'classification/weights/'

DEFAULT_HYPERPARAMS = {
    'NTrees': 400,
    'MaxDepth': 3,
    'MinNodeSize': 0.05,
    'nCuts': 35,
    'BoostType': 'AdaBoost',
    'AdaBoostBeta': 0.5,
    'SeparationType': 'GiniIndex',
}

class BDT(object):

    def __init__(self,
                 name = 'BDT',
                 job_name = 'job',
                 factory = 'Silent:!DrawProgressBar:Transformations=I:AnalysisType=Classification',
                 dataset = '',  
                 hyperparams = {},
                 **kwargs):

        self.logger = logging.getLogger('BDT')
        self.logger.info('Initialized for BDT "{!s}"'.format(name))
        
        self.name = name
        self.job_name = job_name
        self.factory = factory

        if dataset is '':
            raise ValueError, 'No dataset provided.'
        elif not os.path.isfile(DATASET_DIR + dataset + '.root'):
            raise ValueError, 'Dataset file does not exist.'
        else:
            self.dataset = dataset

        if not hyperparams:
            self.logger.info('Using default hyperparameter set.')
            self.hyperparams = DEFAULT_HYPERPARAMS
        else:
            self.hyperparams = hyperparams

    def train(self):

        # Set ROOT to batch mode.
        ROOT.gROOT.SetBatch(1)

        ROOT.TMVA.Tools.Instance() 

        # Output Directory for Monitoring File
        try:
            os.makedirs(MONITOR_DIR)
        except OSError:
            if not os.path.isdir(MONITOR_DIR):
                raise

        # Change the TMVA weight file directory.
        ioNames = ROOT.TMVA.gConfig().GetIONames()
        ioNames.fWeightFileDir = WEIGHT_DIR
        # Uncomment below to change *.weight.xml to *.foo.xml.
        #ioNames.fWeightFileExtension = 'foo'

        # from datetime import datetime, include timestamp in file name.
        #timestamp = datetime.now().strftime('%d%b%Y_%H%M%S')
        monitor_path = MONITOR_DIR + '{!s}_{!s}.root'.format(self.job_name, self.name)
        monitor_file = ROOT.TFile(monitor_path, 'recreate')

        factory = ROOT.TMVA.Factory(self.job_name, monitor_file, self.factory)

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

        bdt_options = ':'.join('{!s}={!s}'.format(*x) for x in self.hyperparams.iteritems())
        factory.BookMethod(ROOT.TMVA.Types.kBDT, self.name, '!H:!V:' + bdt_options)

        factory.TrainAllMethods()
        factory.TestAllMethods()
        factory.EvaluateAllMethods()

        self.logger.info('BDT "{!s}" trained using options "{!s}"'.format(self.name, bdt_options))

    def apply(self, model_name = ''):
    
        pass
        #reader = ROOT.TMVA.Reader()

        #variables = {}

        #for var in VARIABLES.itervalues():


#------
# Main
#------

if __name__ == '__main__':

    # Setup logging.
    logging.basicConfig(level = logging.INFO,
                        format = '%(name)s(%(levelname)s) - %(message)s')

    # Create the dataset if it does not exist.
    dataset = '20Jan2016'
    if not os.path.isfile(DATASET_DIR + dataset + '.root'):
        Dataset(dataset).make()

    # Run the BDT training and application (coming soon).
    bdt = BDT(job_name = 'test', dataset = dataset)

    bdt.train()

