# sean-jiun.wang@cern.ch
# Proverbs 22:29

import glob
import multiprocessing as mp
import subprocess as sp
import tempfile as tf

import numpy as np
import ROOT


def copy_and_skim(subtuple = '', indir = '', outdir = '', overwrite = False, selection = ''):
    
    """
    Copy a Step1 subtuple from EOS and skim the events in tree.

    Parameters
    ----------
    subtuple  : str
                The name of the subtuple. This argument is positioned first 
                for the multiprocessing.Pool.apply_async() method to work.
    indir     : str
                The path to the parent directory of the input file.
    outdir    : str
                The path to the output directory.
    overwrite : bool
                Flag whether the input file is forcibly copied.
                Default is False.
    selection : str
                A TCut style string used for a preliminary cut on the ntuple events
                to reduce its file size.
    """

    # Copy the subtuple.
    if overwrite:
        cmsStage = ['cmsStage', '-f', indir + subtuple, outdir]
    else:
        cmsStage = ['cmsStage', indir + subtuple, outdir]

    return_code = sp.call(cmsStage, stdout = sp.PIPE, stderr = sp.PIPE)

    # Catch failed copies and skim the subtuple if a selection is given.
    if return_code:
        return subtuple
    else:
        if not selection:
            print '--- {} copied.'.format(subtuple)
        else:
                
            # Open the subtuple, caching the tree and count histograms.
            pre_skim = ROOT.TFile(outdir + subtuple, 'read')
            tree = pre_skim.Get('tree')
            n_tree = tree.GetEntriesFast()

            Count = pre_skim.Get('Count')
            CountWeighted = pre_skim.Get('CountWeighted')
            CountPosWeight = pre_skim.Get('CountPosWeight')
            CountNegWeight = pre_skim.Get('CountNegWeight')

            # Create an output file storing the skimmed tree and histograms.
            post_skim = ROOT.TFile(outdir + 'skim_' + subtuple, 'recreate')
            skim_tree = tree.CopyTree(selection)
            n_skim_tree = skim_tree.GetEntriesFast()

            skim_tree.Write()
            Count.Write()
            CountWeighted.Write()
            CountPosWeight.Write()
            CountNegWeight.Write()
        
            # Closing the files ensures the objects in memory are properly saved.
            post_skim.Close()
            pre_skim.Close()

            # Delete the original subtuple.
            sp.check_call(['rm', outdir + subtuple])

            print '--- {} copied and skimmed from {!s} to {!s} entries.'.format(subtuple, n_tree, n_skim_tree)

            return (n_skim_tree, n_tree)


def step2(sample = '', eos_dir = '', step2_dir = '', step2_cut = '', overwrite = False, hadd = True):

    """
    Retrieve and process a Step1 ntuple.

    Parameters
    ----------
    sample    : str
                A minimally descriptive name used to distinguish the ntuple.
    eos_dir   : str
                The Step1 ntuple's EOS directory path.
    step2_dir : str
                The path to the directory storing the Step2 ntuples.
    step2_cut : str
                A TCut style string defining a preliminary cut on the Step1 
                ntuple to reduce the file size of the Step2 ntuple.
    overwrite : bool
                Flag whether the ntuple is forcibly copied by cmsStage.
                The default is False.
    hadd      : bool
                Flag whether the ntuples should be combined into a single file.
                The default is True.
    """

    print '\n[{}]'.format(sample)

    ###################
    #-- Find Ntuple --#
    ###################

    print 'Searching EOS directory "{}"'.format(eos_dir)

    # Retrieve the alias for eos command. Inspired by amaltaro.
    with open('/afs/cern.ch/project/eos/installation/cms/etc/setup.sh', 'r') as f:
        for l in f:
            if 'alias eos=' in l:
                eos = l.split('"')[1]

    # Recursively 'eos ls' until the ntuple is found.
    out = ''    

    while '.root' not in out:

        # The path to the ntuple may not be unique, causing the function to abort.
        # Specify an unambiguous path in vhbb_config, e.g. path to most recent job.
        assert (out.count('\n') <= 1), 'Failed to find a unique EOS path.'

        # Modify the depth of the EOS path.
        eos_dir += out.rstrip() + '/'

        eos_ls = sp.Popen([eos, 'ls', eos_dir], stdout = sp.PIPE, stderr = sp.PIPE)
        out, err = eos_ls.communicate()

    # Collect the ntuple. This works even if it is split across multiple files.
    ntuples = [_ for _ in out.rstrip().split('\n') if '.root' in _]

    ##################
    #-- Get Ntuple --#
    ##################
        
    if step2_cut:
        print 'Skimming with selection "{}"'.format(step2_cut)
        
    # Collect kwargs for the external function copy_and_skim.
    cas_kwargs = {'indir'    : eos_dir,
                  'outdir'   : step2_dir,
                  'selection': step2_cut,
                  'overwrite': overwrite}

    if hadd:
        # Use a temporary directory to handle multiple files.
        tmpdir = tf.mkdtemp(prefix = sample + '_', dir = step2_dir)
        cas_kwargs['outdir'] = tmpdir + '/'

    # Copy and skim the ntuple in parallel.
    pool = mp.Pool(processes = 4)
    results = [pool.apply_async(copy_and_skim, (_,), cas_kwargs).get() for _ in ntuples]
    pool.close()
    pool.join()

    # Report any failures.
    failures = [_ for _ in results if '.root' in _]
    for fail in failures:
        print '--- Failed to get {}'.format(fail)
        results.remove(fail)

    # Report a total sum before and after skimming.
    if step2_cut:
        skim_total, total = 0, 0
        for result in results:
            skim_total += result[0]
            total += result[1]
        print 'Total Number of Entries after Skim (before Skim): {!s} ({!s})'.format(skim_total, total)
       
    # hadd the separate subtuples into a single ntuple.
    if hadd:

        inputfiles = glob.glob(tmpdir + '/*.root')
        outputfile = '{}{}.root'.format(step2_dir, sample)
            
        # Redirecting more output than sp.PIPE's buffer causes it to hang.
        # See "Subprocess Hanging: PIPE is your enemy" by Anders Pearson.
        # The solution is to redirect the output to a temporary file.
        hadd_log = tf.TemporaryFile(dir = step2_dir)

        sp.check_call(['hadd', '-f', outputfile] + inputfiles, stdout = hadd_log, stderr = hadd_log)
            
        # It is the users responsibility to delete the temporary directory.
        sp.check_call(['rm', '-r', tmpdir])
  
def write_sample_lumi(sample = '', step2_dir = '', xsec = None):
    
    print '\nWriting Sample Luminosity Branch for {}'.format(sample)
    print 'Cross Section: {} pb'.format(xsec)
    
    ntuple = ROOT.TFile('{}{}.root'.format(step2_dir, sample), 'update')
    tree = ntuple.Get('tree')
    n_entries = tree.GetEntriesFast()

    # Create the new sample luminosity branch.
    sample_lumi_address = np.array([-999], np.float32)
    sample_lumi_branch = tree.Branch('sample_lumi', sample_lumi_address, 'sample_lumi/F')

    # Calculate the sample luminosity.
    CountPosWeight = ntuple.Get('CountPosWeight')
    CountNegWeight = ntuple.Get('CountNegWeight')

    n_pos = CountPosWeight.GetBinContent(1)
    n_neg = CountNegWeight.GetBinContent(1)

    sample_lumi_address[0] = (n_pos - n_neg) / xsec
    print 'Sample Luminosity: {} pb-1'.format(sample_lumi_address[0])

    # Fill and save the new branch.
    for i in range(0, n_entries):
        tree.GetEntry(i)
        sample_lumi_branch.Fill()

    tree.Write()
    ntuple.Close()    

 
if __name__ == '__main__':

    from settings import *
 
    ##########################
    # Command Line Arguments #
    ##########################

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('samples', nargs = '*', default = [s for s in SAMPLES], 
                        help = 'The samples for which to generate Step2 ntuples.')
    args = parser.parse_args()

    #########################
    # Run the Step2 Program #
    #########################

    print "Launching Step2..."

    # Set ROOT to run in batch mode. 
    ROOT.gROOT.SetBatch(1)
 
    # Create the Step2 ntuple directory if it doesn't exist.
    if (ROOT.gSystem.AccessPathName(STEP2_DIR)):
        ROOT.gSystem.mkdir(STEP2_DIR)
    
    # Generate the Step2 ntuples.
    for sample, properties in SAMPLES.iteritems():
    
        if sample in args.samples:
            step2(sample, properties['EOS_DIR'], STEP2_DIR, STEP2_CUT)
    
            if 'XSEC' in properties:
                write_sample_lumi(sample, STEP2_DIR, properties['XSEC'])
    
    print "\nJob's done!" 
    
