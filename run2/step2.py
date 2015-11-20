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

            return (n_tree, n_skim_tree)
    
class Step2(object):

    """
    Retrieve and process a Step1 ntuple.
    """

    def __init__(self, label = ''):
    
        """
        Parameters
        ----------
        label : str
                A minimally descriptive name used to distinguish the ntuple.
        """
 
        print '\n[{}]'.format(label)
        self.label = label

    def find_subtuples(self, eos_dir = ''):

        """
        Parameters
        ----------
        eos_dir : str
                  The eos path to the ntuple directory.
        """

        print 'Searching EOS path {}'.format(eos_dir)

        # Retrieve the alias for eos command. Inspired by amaltaro.
        with open('/afs/cern.ch/project/eos/installation/cms/etc/setup.sh', 'r') as f:
            for line in f:
                if 'alias eos=' in line:
                    eos = line.split('"')[1]

        # Recursively 'eos ls' until the subtuples are found.
        out = ''    

        while '.root' not in out:

            # It may be that the path to the subtuple is not unique.
            # If so, you'll have to specify the path to use,
            # e.g. the path to the more recent crab job.
            assert (out.count('\n') <= 1), 'Branching search path.'
            
            eos_dir += out.rstrip() + '/'
            eos_ls = sp.Popen([eos, 'ls', eos_dir], stdout = sp.PIPE, stderr = sp.PIPE)
            out, err = eos_ls.communicate()

        subtuples = [_ for _ in out.rstrip().split('\n') if '.root' in _]

        # Set subtuple attributes.
        self.eos_dir = eos_dir
        self.subtuples = subtuples

    def get_subtuples(self, outdir = '', selection = '', overwrite = False, hadd = True):
    
        """
        Parameters
        ----------
        outdir    : str
                    The path to the output directory. This will usually be set as the
                    directory where the Step2 ntuples will be collected.
        selection : str
                    A TCut style string used for a preliminary cut on the ntuple events
                    to reduce its file size.
        overwrite : bool
                    Flag whether the subtuples should be forcibly copied.
                    Default is False.
        hadd      : bool
                    Flag whether the subtuples should be combined into a single file.
                    Default is True.
        """
        
        if selection:
            print 'Skimming with selection {}'.format(selection)
        
        # Collect kwargs for the external function copy_and_skim.
        kwargs = {'indir': self.eos_dir,
                  'outdir': outdir,
                  'overwrite': overwrite,
                  'selection': selection}

        if hadd:
            # Use a temporary directory instead.
            tmpdir = tf.mkdtemp(prefix = self.label + '_', dir = outdir)
            kwargs['outdir'] = tmpdir + '/'

        # Copy and skim the trees in parallel.
        pool = mp.Pool(processes = 4)
        results = [pool.apply_async(copy_and_skim, (_,), kwargs).get() for _ in self.subtuples]
        pool.close()
        pool.join()

        # Report any failures.
        failures = [_ for _ in results if '.root' in _]
        for fail in failures:
            print '--- Failed to get {}'.format(fail)
            results.remove(fail)

        # Report a total sum before and after skimming.
        if selection:
            total, skim_total = 0, 0
            for result in results:
                total += result[0]
                skim_total += result[1]
            print 'Total Number of Entries after Skim (before Skim): {!s} ({!s})'.format(skim_total, total)
       
        # hadd the separate subtuples into a single ntuple.
        if hadd:

            inputfiles = glob.glob(tmpdir + '/*.root')
            outputfile = outdir + self.label + '.root'
            
            # Redirecting the output of hadd to sp.PIPE causes it to hang.
            # "Subprocess Hanging: PIPE is your enemy" - anders pearson
            # To get around this, store the output into a temporary file.
            hadd_log = tf.TemporaryFile(dir = outdir)

            sp.check_call(['hadd', '-f', outputfile] + inputfiles, stdout = hadd_log, stderr = hadd_log)
            
            # It is the users responsibility to delete the temporary directory.
            sp.check_call(['rm', '-r', tmpdir])
                    

if __name__ == '__main__':

    #############################
    # User Specified Parameters #
    #############################

    # Step2 Ntuple Directory Path
    step2_dir = '/afs/cern.ch/work/s/swang373/private/V14/'

    # Skimming Selection
    selection = '(Vtype>0 && met_pt>150)'

    # Labels and EOS Paths for the Step1 Ntuples
    ntuples = {
    # Datasets
    'MET_RunC'       : '/store/group/phys_higgs/hbb/ntuples/V15/MET/VHBB_HEPPY_V15_MET__Run2015C_25ns-05Oct2015-v1',
    'MET_RunD'       : '/store/group/phys_higgs/hbb/ntuples/V15/MET/VHBB_HEPPY_V15_MET__Run2015D-05Oct2015-v1',
    'MET_RunD_Prompt': '/store/group/phys_higgs/hbb/ntuples/V15/MET/VHBB_HEPPY_V15_MET__Run2015D-PromptReco-v4',
    # Signals
    'ZnnH125': '/store/group/phys_higgs/hbb/ntuples/V14/ZH_HToBB_ZToNuNu_M125_13TeV_amcatnloFXFX_madspin_pythia8',
    'ggZH125': '/store/group/phys_higgs/hbb/ntuples/V14/ggZH_HToBB_ZToNuNu_M125_13TeV_amcatnlo_pythia8',
    'WlnH125': '/store/group/phys_higgs/hbb/ntuples/V14/WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8',
    # W+Jets
    'WJetsHT100': '/store/group/phys_higgs/hbb/ntuples/V14/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'WJetsHT200': '/store/group/phys_higgs/hbb/ntuples/V14/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'WJetsHT400': '/store/group/phys_higgs/hbb/ntuples/V14/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'WJetsHT600': '/store/group/phys_higgs/hbb/ntuples/V14/WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'WJetsIncl' : '/store/group/phys_higgs/hbb/ntuples/V14/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    # Z+Jets
    'ZJetsHT100': '/store/group/phys_higgs/hbb/ntuples/V14/ZJetsToNuNu_HT-100To200_13TeV-madgraph',
    'ZJetsHT200': '/store/group/phys_higgs/hbb/ntuples/V14/ZJetsToNuNu_HT-200To400_13TeV-madgraph',
    'ZJetsHT400': '/store/group/phys_higgs/hbb/ntuples/V14/ZJetsToNuNu_HT-400To600_13TeV-madgraph',
    'ZJetsHT600': '/store/group/phys_higgs/hbb/ntuples/V14/ZJetsToNuNu_HT-600ToInf_13TeV-madgraph',
    # TT
    'TTPow': '/store/group/phys_higgs/hbb/ntuples/V14/TT_TuneCUETP8M1_13TeV-powheg-pythia8',
    # Single Top s-channel
    # Single Top t-channel
    # Single Top Wt-channel
    'T_tW'   : '/store/group/phys_higgs/hbb/ntuples/V14/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
    'Tbar_tW': '/store/group/phys_higgs/hbb/ntuples/V14/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
    # Single Top Leptonic
    'T_s_comb_lep': '/store/group/phys_higgs/hbb/ntuples/V14/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1',
    'T_t_lep'     : '/store/group/phys_higgs/hbb/ntuples/V14/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
    'Tbar_t_lep'  : '/store/group/phys_higgs/hbb/ntuples/V14/ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1', 
    # QCD
    'QCDHT100' : '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'QCDHT200' : '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'QCDHT300' : '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'QCDHT500' : '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'QCDHT700' : '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'QCDHT1000': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'QCDHT1500': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'QCDHT2000': '/store/group/phys_higgs/hbb/ntuples/V14/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    # Diboson
    'WW': '/store/group/phys_higgs/hbb/ntuples/V14/WW_TuneCUETP8M1_13TeV-pythia8',
    'WZ': '/store/group/phys_higgs/hbb/ntuples/V14/WZ_TuneCUETP8M1_13TeV-pythia8',
    'ZZ': '/store/group/phys_higgs/hbb/ntuples/V14/ZZ_TuneCUETP8M1_13TeV-pythia8/VHBB_HEPPY_V14_ZZ_TuneCUETP8M1_13TeV-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151025_083230',
    }

    ##########################
    # Command Line Arguments #
    ##########################

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('labels', nargs = '*', default = [_ for _ in ntuples], 
                        help = 'The label(s) specifying which Step1 ntuples to process.')
    args = parser.parse_args()

    ######################
    # Main Step2 Program #
    ######################

    print "Launching Step2..."

    # Set ROOT to run in batch mode. 
    ROOT.gROOT.SetBatch(1)
 
    # Create the Step2 ntuple directory if one doesn't exist.
    if (ROOT.gSystem.AccessPathName(step2_dir)):
        ROOT.gSystem.mkdir(step2_dir)

    # Process the given samples.
    for label in args.labels:
        task = Step2(label)
        task.find_subtuples(ntuples[label])
        task.get_subtuples(step2_dir, selection) 

    print "\nJob's done!" 

