from itertools import izip_longest
import subprocess as sp

import ROOT

from settings import *
from skim import skim


def step3(group = ''):

    """
    Generate a Step3 ntuple.

    Parameters
    ----------
    group   : str
              The name used to refer to the group in settings.py.
    """

    print '\nGenerating Step3 Ntuple [{}]'.format(group)

    # Set the group as the Step3 ntuple's name.
    step3_file = STEP3_DIR + group + '.root'

    # Look up the group's samples and cuts.
    samples = GROUPS[group]['SAMPLES'] if 'SAMPLES' in GROUPS[group] else []
    cuts = GROUPS[group]['CUTS'] if 'CUTS' in GROUPS[group] else []
    
    # When samples is empty, assume the group name references a Step2 ntuple.
    if not samples:

        assert group in SAMPLES, '{} is an invalid sample!'.format(group)

        # A symlink is made to avoid duplicating the Step2 ntuple.
        sp.check_call(['ln', '-s', STEP2_DIR + group + '.root', step3_file]) 

    else:
        
        for s in samples:
            assert s in SAMPLES, '{} is an invalid sample!'.format(s)

        # Only one sample provided.
        if (len(samples) == 1):

            ntuple = samples[0] + '.root'

            # Create a symlink if there are no cuts.
            if not cuts:
                sp.check_call(['ln', '-s', STEP2_DIR + ntuple, step3_file])
            # Otherwise, skim the ntuple and then move the skimmed version.
            else:
                skim(STEP2_DIR, ntuple, cuts[0])
                sp.check_call(['mv', STEP2_DIR + 'skim_' + ntuple, step3_file])
        
        # Multiple samples provided.
        elif (len(samples) > 1):

            inputfiles = []        

            # Iterate over the samples and cuts simultaneously. This is why
            # the samples and their cuts must be provided in the same order.
            for sample, cut in izip_longest(samples, cuts, fillvalue = ''):
                
                ntuple = sample + '.root'
                
                # Append the original ntuple if there are no cuts.
                if not cut:
                    inputfiles.append(STEP2_DIR + ntuple)
                # Otherwise, append the skimmed ntuple.
                else:
                    skim(STEP2_DIR, ntuple, cut)
                    inputfiles.append(STEP2_DIR + 'skim_' + ntuple)

            # hadd the files together.             
            sp.check_call(['hadd', '-f', step3_file] + inputfiles)

            # Delete the skimmed ntuples after hadding.
            for skim_file in [_ for _ in inputfiles if 'skim' in _]:
                sp.check_call(['rm', skim_file])

#-------------
# Main Program
#-------------

if __name__ == '__main__':

    # Set ROOT to batch mode.
    ROOT.gROOT.SetBatch(1)

    # Parse command line arguments.
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('groups', nargs = '*', default = [g for g in GROUPS],
                        help = 'The groups for which to generate Step3 ntuples.')
    args = parser.parse_args()

    # Create the Step3 ntuple directory if it doesn't exist.
    if (ROOT.gSystem.AccessPathName(STEP3_DIR)):
        ROOT.gSystem.mkdir(STEP3_DIR) 

    # Generate the Step3 ntuples.
    print 'Running step3.py...'

    for group in GROUPS:
        if group in args.groups:
            step3(group)

    print "\nJob's done!"

