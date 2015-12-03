import ROOT


def skim(indir = '', ntuple = '', cut = ''):
    
    """
    Performs the common task of skimming an ntuple with a cut.
    The output is created in the same directory as the input.
    
    Parameters
    ----------
    indir  : str
             The full path to the ntuple's directory.
    ntuple : str
             The name of the ntuple to be skimmed.
    cut    : str
           A TCut style string which defines the skimming cut.
    """

    # Open the ntuple and access the tree and count histograms.
    infile = ROOT.TFile(indir + ntuple, 'read')
    tree = infile.Get('tree')
    n_tree = tree.GetEntriesFast()

    Count = infile.Get('Count')
    CountWeighted = infile.Get('CountWeighted')
    CountPosWeight = infile.Get('CountPosWeight')
    CountNegWeight = infile.Get('CountNegWeight')

    # Create an output file to store the skimmed tree and histograms. 
    outfile = ROOT.TFile(indir + 'skim_' + ntuple, 'recreate')
    skim_tree = tree.CopyTree(cut)
    n_skim_tree = skim_tree.GetEntriesFast()
   
    skim_tree.Write()
    Count.Write()
    CountWeighted.Write()
    CountPosWeight.Write()
    CountNegWeight.Write()

    # Close the files to ensure objects in memory are properly saved.
    outfile.Close()
    infile.Close()

    # Return a tuple of the number of entries after and before skimming, respectively.
    return (n_skim_tree, n_tree)

