# sean-jiun.wang@cern.ch
# Proverbs 22:29

from functools import partial
import multiprocessing as mp
import subprocess as sp

import numpy as np
import ROOT


def copy_and_skim(subtuple = '', eos_dir = '', outdir = '', overwrite = False):
    
    if overwrite:
        cmsStage = ['cmsStage', '-f', eos_dir + subtuple, outdir]
    else:
        cmsStage = ['cmsStage', eos_dir + subtuple, outdir]

    return_code = sp.call(cmsStage)

    if return_code:
        return subtuple
    
class Step2(object):

    """
    Copy and skim the Step1 ntuples.
    """

    def __init__(self, label = ''):
    
        """
        Parameters
        ----------
        label : str
                A minimally descriptive name used to distinguish the ntuple.
        """

        # Set ROOT to batch mode.
        ROOT.gROOT.SetBatch(1)

        print '\nInitialize Step2 for {}'.format(label)
        self.label = label

    def find_subtuples(self, eos_dir = ''):

        """
        Parameters
        ----------
        eos_dir : str
                  The eos path to the ntuple directory.
        """

        print 'Collecting subtuples...'

        # Retrieve the alias for eos command. Inspired by amaltaro.
        with open('/afs/cern.ch/project/eos/installation/cms/etc/setup.sh', 'r') as f:
            for line in f:
                if 'alias eos=' in line:
                    eos = line.split('"')[1]

        # Recursively 'eos ls' until the subtuples are found.
        out = ''    

        while '.root' not in out:
            eos_dir += out.rstrip() + '/'
            eos_ls = sp.Popen([eos, 'ls', eos_dir], stdout = sp.PIPE, stderr = sp.PIPE)
            out, err = eos_ls.communicate()

        subtuples = [x for x in out.rstrip().split('\n') if '.root' in x]

        # Set subtuple attributes.
        self.eos_dir = eos_dir
        self.subtuples = subtuples

    def get_subtuples(self, step2_dir = '', selection = '', overwrite = False):
    
        """
        Parameters
        ----------
        step2_dir : str
                    The path to the Step2 ntuple directory.
        """
         
        # Make an instance-specific wrapping for parallelizing an external function.
        cas_wrapper = partial(copy_and_skim, eos_dir = self.eos_dir, 
                                             outdir = self.tmp_dir,
                                             overwrite = overwrite)

        # Copy and skim the subtuples asynchronously.
        pool = mp.Pool(processes = 4)
        results = pool.map_async(cas_wrapper, self.subtuples)
        pool.close()
        pool.join()

        self.failures = [x for x in results.get() if x is not None]
                                        

#-------------
# Main Program
#-------------

if __name__ == '__main__':

    # Create a Step2 ntuple directory if one doesn't exist.
    step2_dir = '/afs/cern.ch/work/s/swang373/private/testing/'
    if (ROOT.gSystem.AccessPathName(step2_dir)):
        ROOT.gSystem.mkdir(step2_dir)


    # Run Step2
    task = Step2('ZnnH125')
    task.find_subtuples('/store/group/phys_higgs/hbb/ntuples/V15/MET/VHBB_HEPPY_V15_MET__Run2015D-PromptReco-v4')
    task.get_subtuples(step2_dir)
    
   
    
