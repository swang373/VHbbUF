import glob
import multiprocessing as mp
import subprocess as sp
import tempfile as tf

import numpy as np
import ROOT

from settings import *
from skim import skim


def copy_and_skim(ntuple = '', indir = '', outdir = '', overwrite = False):
    
    """
    Copy a Step1 ntuple from EOS and skim the events with the Step2 cut.

    Parameters
    ----------
    ntuple    : str
                The file name of the ntuple. This argument is positioned first 
                to conform with the usage of multiprocessing.Pool.apply_async().
    indir     : str
                The path to the parent directory of the input file.
    outdir    : str
                The path to the output directory.
    overwrite : bool
                Flag whether to forcibly copy the subtuple.
                Default is False.
    """

    # Copy the ntuple.
    if overwrite:
        cmsStage = ['cmsStage', '-f', indir + ntuple, outdir]
    else:
        cmsStage = ['cmsStage', indir + ntuple, outdir]

    return_code = sp.call(cmsStage, stdout = sp.PIPE, stderr = sp.PIPE)

    # Catch failed copies.
    if return_code:
        return ntuple
    # Skim the ntuple if a cut is specified.
    else:
        if not STEP2_CUT:
            print '--- {} copied.'.format(ntuple)
        else:
            result = skim(outdir, ntuple, STEP2_CUT)

            # Delete the original ntuple.
            sp.check_call(['rm', outdir + ntuple])

            print '--- {} copied and skimmed from {!s} to {!s} entries.'.format(ntuple, result[1], result[0])

            return result

##########################################

def step2(sample = '', overwrite = False):

    """
    Generate a Step2 ntuple.

    Parameters
    ----------
    sample    : str
                The name used to refer to the sample in settings.py.
    overwrite : bool
                Flag whether to forcibly copy the ntuple.
                The default is False.
    """

    print '\nGenerating Step2 Ntuple [{}]'.format(sample)

    # Set the sample as the Step2 ntuple's name.
    step2_file = STEP2_DIR + sample + '.root'

    # Look up the sample's EOS path and cross section.
    eos_dir = SAMPLES[sample]['EOS_DIR'] if 'EOS_DIR' in SAMPLES[sample] else ''
    xsec = SAMPLES[sample]['XSEC'] if 'XSEC' in SAMPLES[sample] else None

    if not eos_dir:
        print 'Empty EOS path!'
        return

    print 'Searching EOS directory "{}"'.format(eos_dir)

    # Retrieve the alias for the "eos" command. Inspired by amaltaro.
    with open('/afs/cern.ch/project/eos/installation/cms/etc/setup.sh', 'r') as f:
        for l in f:
            if 'alias eos=' in l:
                eos = l.split('"')[1]

    # Recursively build the full EOS path to the ntuple.
    eos_out = ''    

    while '.root' not in eos_out:

        # Abort if a unique EOS path is inaccessible. Provide a full path for
        # 'EOS_DIR' in settings.py,  e.g. full path to most recent CRAB job.
        assert (eos_out.count('\n') <= 1), 'Failed to find a unique EOS path.'

        # Append the previous "eos ls" output to the end of the EOS path.
        eos_dir += eos_out.rstrip() + '/'

        # Retrieve the "eos ls" output for the current EOS path.
        eos_ls = sp.Popen([eos, 'ls', eos_dir], stdout = sp.PIPE, stderr = sp.PIPE)
        eos_out, eos_err = eos_ls.communicate()

    # Collect the .root files, ignoring other file types.
    ntuples = [_ for _ in eos_out.rstrip().split('\n') if '.root' in _]

    # Set kwargs for the external copy_and_skim function.
    cas_kwargs = {'indir': eos_dir, 'outdir': STEP2_DIR, 'overwrite': overwrite}

    # If the ntuple is split across multiple files, use a temporary directory.
    if (len(ntuples) > 1):
        hadd = True
        tmpdir = tf.mkdtemp(prefix = sample, dir = STEP2_DIR)
        cas_kwargs['outdir'] = tmpdir + '/'
    
    # Report whether a Step2 cut was defined for skimming.
    if STEP2_CUT:
        print 'Skimming with cut "{}"'.format(STEP2_CUT)
    
    # Copy and skim the ntuple in parallel. Four processes is a safe default.
    # The (var,) syntax puts the item in a tuple for apply_async to unpack.
    results = []
    pool = mp.Pool(processes = 4)
    for ntuple in ntuples:
        pool.apply_async(copy_and_skim, (ntuple,), cas_kwargs, callback = results.append)
    pool.close()
    pool.join()

    # Report and remove files which failed to copy and skim.
    failures = [_ for _ in results if '.root' in _]
    for fail in failures:
        print '--- Failed to get {}'.format(fail)
        results.remove(fail)

    # Report a total sum before and after skimming if a Step2 cut was defined.
    if STEP2_CUT:
        skim_total, total = 0, 0
        for result in results:
            skim_total += result[0]
            total += result[1]
        print 'Total Number of Entries after Skim (before Skim): {!s} ({!s})'.format(skim_total, total)
       
    # Combine the separate files into a single ntuple.
    if hadd:

        inputfiles = glob.glob(tmpdir + '/*.root') 
        
        # Redirect output to a temporary file. Overflowing the buffer
        # of sp.PIPE causes the function call to hang. See blog post
        # "Subprocess Hanging: PIPE is your enemy" by Anders Pearson.
        hadd_log = tf.TemporaryFile(dir = STEP2_DIR)
        
        # hadd the files together.
        sp.check_call(['hadd', '-f', step2_file] + inputfiles, stdout = hadd_log, stderr = hadd_log)
            
        # It is the user's responsibility to delete the temporary directory.
        sp.check_call(['rm', '-r', tmpdir])

    # Add a sample luminosity branch if a cross section is provided.  
    if not xsec:
        return

    print 'Cross Section: {} pb'.format(xsec)

    # Open the ntuple for updating and access the tree.
    infile = ROOT.TFile(step2_file, 'update')
    tree = infile.Get('tree')

    # Create the new sample luminosity branch.
    sample_lumi_address = np.array([-999], np.float32)
    sample_lumi_branch = tree.Branch('sample_lumi', sample_lumi_address, 'sample_lumi/F')

    # Calculate the sample luminosity.
    n_pos = infile.Get('CountPosWeight').GetBinContent(1)
    n_neg = infile.Get('CountNegWeight').GetBinContent(1)
    sample_lumi_address[0] = (n_pos - n_neg) / xsec
    print 'Sample Luminosity: {} pb-1'.format(sample_lumi_address[0])

    # Fill and save the new branch.
    for i in range(0, tree.GetEntriesFast()):
        tree.GetEntry(i)
        sample_lumi_branch.Fill()

    # Write the information to file and close the ntuple.
    tree.Write()
    infile.Close()

#-------------
# Main Program
#-------------

if __name__ == '__main__':

    # Set ROOT to batch mode.
    ROOT.gROOT.SetBatch(1)

    # Parse command line arguments.
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('samples', nargs = '*', default = [s for s in SAMPLES], 
                        help = 'The samples for which to generate Step2 ntuples.')
    args = parser.parse_args()
 
    # Create the Step2 ntuple directory if it doesn't exist.
    if (ROOT.gSystem.AccessPathName(STEP2_DIR)):
        ROOT.gSystem.mkdir(STEP2_DIR)
    
    # Generate the Step2 ntuples.
    print 'Running step2.py...'

    for sample, properties in SAMPLES.iteritems():  
        if sample in args.samples:
            step2(sample)
       
    print "\nJob's done!" 
    
