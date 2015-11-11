# sean-jiun.wang@cern.ch
# Proverbs 22:29

from functools import partial
import glob
import multiprocessing as mp
import subprocess as sp

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

        # Set ROOT to batch mode.
        ROOT.gROOT.SetBatch(1)

        print '\n[Run Step2 for {}]'.format(label)
        self.label = label

    def find_subtuples(self, eos_dir = ''):

        """
        Parameters
        ----------
        eos_dir : str
                  The eos path to the ntuple directory.
        """

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
         
        # Create a temporary working directory for tidiness.
        tmpdir = outdir + self.label + '/'
        if (ROOT.gSystem.AccessPathName(tmpdir)):
            ROOT.gSystem.mkdir(tmpdir)

        # Make an instance-specific wrapping for parallelizing an external function.
        inst_args = {'indir': self.eos_dir,
                     'outdir': tmpdir,
                     'overwrite': overwrite,
                     'selection': selection}

        partial_copy_and_skim = partial(copy_and_skim, **inst_args)

        # Copy and skim the subtuples asynchronously.
        pool = mp.Pool(processes = 4)
        results = pool.map_async(partial_copy_and_skim, self.subtuples)
        pool.close()
        pool.join()

        # Report any failed jobs.
        self.failures = [x for x in results.get() if x is not None]
        for x in self.failures:
            print '--- Failed to get {}'.format(x)

        # hadd the separate subtuples into a single ntuple.
        if hadd:
            inputfiles = glob.glob(tmpdir + '*.root')
            outputfile = outdir + self.label + '.root'
            sp.check_call(['hadd', '-f', outputfile] + inputfiles, stdout = sp.PIPE, stderr = sp.PIPE)
            sp.check_call(['rm', '-r', tmpdir])


if __name__ == '__main__':

    #############################
    # User Specified Parameters #
    #############################

    # Step2 Ntuple Directory Path

    step2_dir = '/afs/cern.ch/work/s/swang373/private/V14/'

    # Skimming Selection

    selection = '(Vtype>0 && met_pt>150)'

    # Step1 Ntuple Labels and EOS Paths

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
    'ZZ': '/store/group/phys_higgs/hbb/ntuples/V14/ZZ_TuneCUETP8M1_13TeV-pythia8',
    }


    ######################
    # Main Step2 Program #
    ######################

    # Create the Step2 ntuple directory if one doesn't exist.
    if (ROOT.gSystem.AccessPathName(step2_dir)):
        ROOT.gSystem.mkdir(step2_dir)

    task = Step2('MET_2015D_PromptReco')
    task.find_subtuples('/store/group/phys_higgs/hbb/ntuples/V15/MET/VHBB_HEPPY_V15_MET__Run2015D-PromptReco-v4')
    task.get_subtuples(step2_dir, selection)
   
    print "\nJob's done!" 
 
