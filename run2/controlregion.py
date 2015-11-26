import ROOT


class ControlRegion(object):

    def __init__(self, name = '', outdir = '', *cuts):

        self.name = name
        self.cuts = cuts
        self.outfile = ROOT.TFile('{}{}.root'.format(outdir, name), 'recreate')

        print '\n[{}]'.format(name)

    def add_tree(self, name = '', file = '', *add_cuts):

        print 'Adding TTree named "{}" from {}'.format(name, file)

        # Access the input file and TTree.
        infile = ROOT.TFile(file, 'read')
        tree = infile.Get('tree')

        # Change directory to the output file.
        self.outfile.cd()
    
        # Perform the control region cuts.
        if (len(self.cuts) <= 1):
            if not self.cuts:
                CR_tree = tree.CopyTree('')
            else:
                CR_tree = tree.CopyTree(self.cuts[0])
        else:
            CR_tree = tree.CopyTree(self.cuts[0])
            for cut in self.cuts[1:]:
                CR_tree = CR_tree.CopyTree(cut)
        
        # Perform any addtional cuts provided.
        if add_cuts:
            for cut in add_cuts:
                CR_tree = CR_tree.CopyTree(cut)

        print '--- Selected {!s} out of {!s} entries'.format(CR_tree.GetEntriesFast(), tree.GetEntriesFast())
        
        # Save the control region tree.
        CR_tree.SetName(name)
        CR_tree.Write()    
        infile.Close()

    def close(self):
        # Clean up extraneous TTrees and properly close.
        self.outfile.cd()
        ROOT.gDirectory.Delete('tree;*')
        self.outfile.Close()

#-------------
# Main Program
#-------------

if __name__ == '__main__':

    from settings import *

    # Set ROOT to run in batch mode.
    ROOT.gROOT.SetBatch(1)

    # Create the control region directory if it doesn't exist.
    if (ROOT.gSystem.AccessPathName(CONTROL_REGION_DIR)):
        ROOT.gSystem.mkdir(CONTROL_REGION_DIR)

    print 'Generating Control Region Ntuples...'

    for region, definition in CONTROL_REGIONS.iteritems():
 
        CR = ControlRegion(region, CONTROL_REGION_DIR, *definition)

        for category, options in CATEGORIES.iteritems():

            if 'CUTS' in options: 
                CR.add_tree(category, options['PATH'], options['CUTS'])
            else:
                CR.add_tree(category, options['PATH'])

        CR.close()
    
    print "\nJob's done!"

