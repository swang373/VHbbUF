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
                 job_id = 'TMVAClassification',
                 dataset = '',  
                 hyperparams = {},
                 **kwargs):

        self.logger = logging.getLogger('BDT')
        self.logger.info('Initialized for "{}"'.format(job_id))
        
        # Set ROOT to batch mode.
        ROOT.gROOT.SetBatch(1)

        self.job_id = job_id

        if dataset is '':
            raise ValueError, 'No dataset provided.'
        else:
            if not os.path.isfile(DATASET_DIR + dataset + '.root'):
                Dataset(dataset).make()
            self.dataset = dataset

        if not hyperparams:
            self.logger.info('Using default hyperparameter set.')
            self.hyperparams = DEFAULT_HYPERPARAMS
        else:
            self.hyperparams = hyperparams

    def train(self):

        # Output Directory for Monitoring File
        try:
            os.makedirs(MONITOR_DIR)
        except OSError:
            if not os.path.isdir(MONITOR_DIR):
                raise

        # Change the TMVA weight file directory.
        ROOT.TMVA.Tools.Instance()
        ioNames = ROOT.TMVA.gConfig().GetIONames()
        ioNames.fWeightFileDir = WEIGHT_DIR
        # Lines below change *.weight.xml to *.foo.xml.
        #ioNames.fWeightFileExtension = 'foo'

        #timestamp = datetime.now().strftime('%d%b%Y_%H%M%S')
        monitor_path = MONITOR_DIR + self.job_id + '_BDT.root'
        monitor_file = ROOT.TFile(monitor_path, 'recreate')

        factory = ROOT.TMVA.Factory(self.job_id, monitor_file, 'Silent:!DrawProgressBar:Transformations=I:AnalysisType=Classification')

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
        factory.BookMethod(ROOT.TMVA.Types.kBDT, 'BDT', '!H:!V:' + bdt_options)

        factory.TrainAllMethods()
        factory.TestAllMethods()
        factory.EvaluateAllMethods()

        monitor_file.Close()

        self.logger.info('BDT "{!s}" trained using options "{!s}"'.format(self.job_id, bdt_options))

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

    bdt = BDT(dataset = '20Jan2016')
    bdt.train()

