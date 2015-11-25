# sean-jiun.wang@cern.ch
# Proverbs 22:29

import ROOT

def skim_ntuple(selection = '', indir = '', ntuple = ''):
    
    """
    This utility function will perform the common task of skimming a VHbb ntuple with
    a selection. The output ntuple will be created in the same directory as the input.
    
    Parameters
    ----------
    selection : str
                A string in the style of ROOT's TCut.
    indir     : str
                The full path to the ntuple's directory.
    ntuple    : str
                The name of the ntuple to be skimmed.
    """

    print '\nSkimming: %s' % ntuple
    print 'Selection: %s' % selection

    ROOT.gROOT.SetBatch(1)

    infile = ROOT.TFile(indir + ntuple, 'READ')
    intree = infile.Get('tree')
    n_intree = intree.GetEntriesFast()

    # Copy the tree with a selection and cache the count histograms.
    outfile = ROOT.TFile(indir + 'skim_' + ntuple, 'RECREATE')
    outtree = intree.CopyTree(selection)
    n_outtree = outtree.GetEntriesFast()
    print 'Skimmed from %s to %s entries.' % (n_intree, n_outtree)

    Count = infile.Get('Count')
    CountWeighted = infile.Get('CountWeighted')
    CountPosWeight = infile.Get('CountPosWeight')
    CountNegWeight = infile.Get('CountNegWeight')

    # Write the skimmed tree and the Step 1 count histograms.
    outtree.Write()
    Count.Write()
    CountWeighted.Write()
    CountPosWeight.Write()
    CountNegWeight.Write()

    # Properly close the files to save the objects in memory.
    outfile.Close()
    infile.Close()

if __name__ == '__main__':

    skim_ntuple('lheHT<100', '/afs/cern.ch/work/s/swang373/private/V14/', 'WJetsIncl.root')

