import glob
import multiprocessing as mp
import subprocess as sp
import tempfile as tf

import ROOT

from settings import *
import TCutOperators as tco


def make_category(category = '', outdir = '', cuts = []):

    # Look up the path to the category's ntuple.
    path = CATEGORIES[category]['PATH']

    # If a cut defining a subcategory is provided, append it
    # to the list of cuts defining the signal/control region.
    if 'CUTS' in CATEGORIES[category]:
        cuts.append(CATEGORIES[category]['CUTS'])

    # Access the category's ntuple and TTree.
    infile = ROOT.TFile(path, 'read')
    tree = infile.Get('tree')
    tree.SetName(category)

    # Create an output file for the category.    
    outfile = ROOT.TFile(outdir + category + '.root', 'recreate')

    # Create the category's tree by destructively
    # iterating over the region and subcategory cuts.
    category_tree = tree.CopyTree(tco.add(*cuts.pop(0)))
    while cuts:
        category_tree = category_tree.CopyTree(tco.add(*cuts.pop(0)))
        
    print '--- Category "{}", selecting {!s} entries from "{}"'.format(category, category_tree.GetEntriesFast(), path)

    # Save the category tree.
    category_tree.Write()
    outfile.Close()
    infile.Close()
     
#############################

def make_region(region = ''):

    print '\nGenerating Signal/Control Region [{}]'.format(region)

    # Use a temporary directory to collect the different categories.
    tmpdir = tf.mkdtemp(prefix = region, dir = REGION_DIR)

    # Set kwargs for the external make_category function.
    mc_kwargs = {'outdir': tmpdir + '/', 'cuts': REGIONS[region]}

    # Add the categories in parallel. Four processes is a safe default.
    pool = mp.Pool(processes = 4)
    for category in CATEGORIES:
        pool.apply_async(make_category, (category,), mc_kwargs)
    pool.close()
    pool.join()

    # Combine the separate category files into a single ntuple.
    category_files = glob.glob(tmpdir + '/*.root')
    region_file = REGION_DIR + region + '.root'
    sp.check_call(['hadd', '-f', region_file] + category_files)

    # It is the user's responsibility to delete the temporary directory.
    sp.check_call(['rm', '-r', tmpdir])

#-------------
# Main Program
#-------------

if __name__ == '__main__':

    # Set ROOT to run in batch mode.
    ROOT.gROOT.SetBatch(1)

    # Parse command line arguments.
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('regions', nargs = '*', default = [r for r in REGIONS],
                        help = 'The signal/control regions for which to generate .root files.')
    args = parser.parse_args()

    # Create the signal/control region directory if it doesn't exist.
    if (ROOT.gSystem.AccessPathName(REGION_DIR)):
        ROOT.gSystem.mkdir(REGION_DIR)

    print 'Running make_region.py...'

    for region in REGIONS:
        if region in args.regions:
            make_region(region)
    
    print "\nJob's done!"

