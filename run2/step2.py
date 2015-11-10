# sean-jiun.wang@cern.ch
# Proverbs 22:29

from functools import partial
import glob
import multiprocessing as mp
import subprocess as sp

import numpy as np
import ROOT


def copy_and_skim(infile = '', 
                  indir = '', 
                  outdir = '', 
                  label = '', 
                  overwrite = False,
                  selection = ''):
    
    # Form file paths.
    eos_path = indir + infile
    copy_path = outdir + label + '_' + infile

    if overwrite:
        cmsStage = ['cmsStage', '-f', eos_path, copy_path]
    else:
        cmsStage = ['cmsStage', eos_path, copy_path]

    return_code = sp.call(cmsStage, stdout = sp.PIPE, stderr = sp.PIPE)

    if return_code:
        # Report a copy failure.
        return infile
    elif not return_code and selection != '':
        # If a selection is given, skim the file.

        skim_path = outdir + 'skim_' + label + '_' + infile

        # Open the file and cache count histograms.
        pre_skim = ROOT.TFile(copy_path, 'read')
        tree = pre_skim.Get('tree')
        n_tree = tree.GetEntriesFast()

        Count = pre_skim.Get('Count')
        CountWeighted = pre_skim.Get('CountWeighted')
        CountPosWeight = pre_skim.Get('CountPosWeight')
        CountNegWeight = pre_skim.Get('CountNegWeight')

        # Create the skimmed file.
        post_skim = ROOT.TFile(skim_path, 'recreate')
        skim_tree = tree.CopyTree(selection)
        n_skim_tree = skim_tree.GetEntriesFast()

        skim_tree.Write()
        Count.Write()
        CountWeighted.Write()
        CountPosWeight.Write()
        CountNegWeight.Write()
      
        post_skim.Close()
        pre_skim.Close()

        print '--- {} skimmed from {!s} to {!s} entries.'.format(infile, n_tree, n_skim_tree)

        sp.check_call(['rm', copy_path])
    
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

        print 'Finding subtuples...'

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

    def get_subtuples(self, step2_dir = '', selection = '', overwrite = False, hadd = True):
    
        """
        Parameters
        ----------
        step2_dir : str
                    The path to the Step2 ntuple directory.
        """
        
        print 'Getting subtuples...'    
 
        # Make an instance-specific wrapping for parallelizing an external function.
        kwargs = {'indir': self.eos_dir,
                  'outdir': step2_dir,
                  'label': self.label,
                  'overwrite': overwrite,
                  'selection': selection}

        partial_copy_and_skim = partial(copy_and_skim, **kwargs)

        # Copy and skim the subtuples asynchronously.
        pool = mp.Pool(processes = 4)
        results = pool.map_async(partial_copy_and_skim, self.subtuples)
        pool.close()
        pool.join()

        self.failures = [x for x in results.get() if x is not None]

        for x in self.failures:
            print 'Failed to get {}'.format(x)

        if hadd:
            infiles = glob.glob(step2_dir + '*.root')
            outputfile = step2_dir + self.label + '.root'
            sp.check_call(['hadd', '-f', outputfile] + infiles, stdout = sp.PIPE, stderr = sp.PIPE)
            sp.check_call(['rm'] + infiles)
                                        

#-------------
# Main Program
#-------------

if __name__ == '__main__':

    # Create a Step2 ntuple directory if one doesn't exist.
    step2_dir = '/afs/cern.ch/work/s/swang373/private/test/'
    if (ROOT.gSystem.AccessPathName(step2_dir)):
        ROOT.gSystem.mkdir(step2_dir)

    selection = 'met_pt>150'

    # Run Step2
    task = Step2('ZnnH125')
    task.find_subtuples('/store/group/phys_higgs/hbb/ntuples/V15/MET/VHBB_HEPPY_V15_MET__Run2015D-PromptReco-v4')
    task.get_subtuples(step2_dir, selection)
    
   
    
