import ROOT


class TChainSaw(object):

    """
    Handle a set of TChains composed of VHbb Heppy ntuples.
    Defining cuts over the TChains creates TTree attributes
    for access to events which pass a selection.
    """

    def __init__(self):
        
        self.processes = {}
        self.trees = []

    def add_process(self, name = '', ntuple = '', selection = ''):
        
        """
        Parameters
        ----------
        name      : str
                    The name for the TTree attribute. Also used by the
                    processes dictionary to key the process' TChain.
        ntuple    : str
                    The full file path of the ntuple to add to the TChain.
        selection : str
                    A TCut style string defining the cuts used to create the TTree.
                    If selection is empty, a TTree attribute is not created.
        """
            
        # Add the process.
        print 'Adding {} as process {}...'.format(ntuple, name)
        self.processes[name] = ROOT.TChain('tree')
        self.processes[name].Add(ntuple)
        
        if selection:
            # Set the TTree which passes the selection to be a named attribute. 
            setattr(self, name, self.processes[name].CopyTree(selection))
            self.trees.append(name)
            print '--- Created TTree {} with {!s} entries'.format(name, getattr(self, name).GetEntries())
        else:
            print '--- Loaded TChain'
        
    def add_subprocess(self, name = '', parent = '', selection = ''):

        """
        Parameters
        ----------
        name      : str
                    The name for the TTree attribute. Also used by the
                    processes dictionary to key its parent process' name.
        parent    : str
                    The name of the parent process whose TChain will be used.
        selection : str
                    A TCut style string defining the cuts used to create the TTree.
        """

        assert parent in self.processes, 'Invalid parent process!'
 
        # Add the subprocess.
        print 'Adding {} as subprocess of {}...'.format(name, parent)
        self.processes[name] = parent

        # Set the TTree which passes the selection to be a named attribute.
        setattr(self, name, self.processes[parent].CopyTree(selection))
        self.trees.append(name)
        print '--- Created TTree {} with {!s} entries'.format(name, getattr(self, name).GetEntries())
            

def make_plot(name = '', tchainsaw = None, expression = '', n_xbins = None, x_min = None, x_max = None):
   
    """
    Parameters
    ----------
    name       : str
                 The string used to prefix the plot's name.
    tchainsaw  : TChainSaw
                 The object managing the selected TTrees.
    expression : str
                 A TTreeFormula style string defining the expression to be passed to
                 TTree::Project(). Refer to TTree::Draw() documentation for examples.
    n_xbins    : int
                 The number of bins along the x-axis.
    x_min      : float
                 The lower bound of the x-axis.
    x_max      : float
                 The upper bound of the x-axis.
    """

    hist = {}

    #categories: ZH, ggZH, WH, WjLF, WjHF, ZjLF, ZjHF, TT, ST, VV, QCD, data

    # Book histograms and project the TTrees.
    for tree in tchainsaw.trees:
        hist[tree] = ROOT.TH1F(tree, '', n_xbins, x_min, x_max)
        getattr(tchainsaw, tree).Project(tree, expression, 'sign(genWeight)')

    # Book any additional histograms.
    hist['VH'] = ROOT.TH1F('VH', '', n_xbins, x_min, x_max)
    hist['mc_exp'] = ROOT.TH1F('mc_exp', '', n_xbins, x_min, x_max)

    # Combine the two VH plots.
    hist['VH'].Add(hist['ZH'])
    hist['VH'].Add(hist['WH'])

    # Combine all the backgrounds.
    hist['mc_exp'].Add(hist['WjLF'])
    hist['mc_exp'].Add(hist['WjHF'])
    hist['mc_exp'].Add(hist['ZjLF'])
    hist['mc_exp'].Add(hist['ZjHF'])
    hist['mc_exp'].Add(hist['TT'])
    hist['mc_exp'].Add(hist['ST'])
    hist['mc_exp'].Add(hist['VV'])
    hist['mc_exp'].Add(hist['QCD'])

    # Make the histogram stack.
    hstack = ROOT.THStack('hstack','')
    hstack.Add(hist['ST'])
    hstack.Add(hist['TT'])
    hstack.Add(hist['ZjLF'])
    hstack.Add(hist['ZjHF'])
    hstack.Add(hist['WjLF'])
    hstack.Add(hist['WjHF'])
    hstack.Add(hist['VV'])
    hstack.Add(hist['QCD'])
    hstack.Add(hist['VH'])

    # Canvas and pads.
    canvas = ROOT.TCanvas('canvas', '', 700, 700)
    upper_pad = ROOT.TPad('upper_pad', '', 0.0, 0.3, 1.0, 1.0)
    upper_pad.SetBottomMargin(0.0)
    upper_pad.Draw()
    lower_pad = ROOT.TPad('lower_pad', '', 0.0, 0.0, 1.0, 0.3)
    lower_pad.SetTopMargin(0.0)
    lower_pad.SetBottomMargin(0.35)
    lower_pad.Draw()
    upper_pad.cd()

    # Set histogram styles.
    for h in hist:

        color = {'WjLF': 814, 'WjHF': 820, 'ZjLF': 401, 'ZjHF': 5,
                 'TT': 596, 'ST': 840, 'VV': 922, 'QCD': 616}

        if (h == 'VH' or h = 'ZH', or h = 'WH'):
            hist[h].SetFillColor(2)
            hist[h].SetMarkerColor(2)
        elif (h == 'data'):
            hist[h].SetMarkerSize(0.8)
            hist[h].SetMarerStyle(20)
        elif (h in color):
            hist[h].SetLineColor(ROOT.kBlack)
            hist[h].SetFillColor(color[h])
            hist[h].SetMarkerColor(color[h])
        else:
            continue

    # Setup auxiliary histograms.
    hist['stat_unc'] = hist['mc_exp'].Clone('stat_unc')
    hist['stat_unc'].Sumw2()
    hist['stat_unc'].SetFillColor(ROOT.kGray+3)
    hist['stat_unc'].SetMarkerSize(0)
    hist['stat_unc'].SetFillStyle(3013)

    hist['ratio'] = hist['data'].Clone('ratio')
    hist['ratio'].Sumw2()
    hist['ratio'].SetMarkerSize(0.8)
    hist['ratio'].Divide(hist['data'], hist['mc_exp'], 1., 1., '')
    
    hist['ratio_stat'] = hist['mc_exp'].Clone('ratio_stat')
    hist['ratio_stat'].Sumw2()
    hist['ratio_stat'].SetStats(0)
    hist['ratio_stat'].GetYaxis().SetTitle('Data/MC')
    hist['ratio_stat'].SetMaximum(2.2)
    hist['ratio_stat'].SetMinimum(0.0)
    hist['ratio_stat'].SetMarkerSize(0)
    hist['ratio_stat'].SetFillColor(ROOT.kGray+3)
    hist['ratio_stat'].SetFillStyle(3013)
    hist['ratio_stat'].GetXaxis().SetLabelSize(0.12)
    hist['ratio_stat'].GetXaxis().SetTitleSize(0.14)
    hist['ratio_stat'].GetXaxis().SetTitleOffset(1.10)
    hist['ratio_stat'].GetYaxis().SetLabelSize(0.10)
    hist['ratio_stat'].GetYaxis().SetTitleSize(0.12)
    hist['ratio_stat'].GetYaxis().SetTitleOffset(0.6)
    hist['ratio_stat'].GetYaxis().SetNdivisions(505)

    for i in range(1, n_xbins+1):
        hist['ratio_stat'].SetBinContent(i, 1.0)
        if (hist['mc_exp'].GetBinContent(i) > 0):
            bin_error = hist['mc_exp'].GetBinError(i) / hist['mc_exp'].GetBinContent(i)
            hist['ratio_stat'].SetBinError(i, binerror)
        else:
            hist['ratio_stat'].SetBinError(i, 999.)

    ratio_unity = ROOT.TLine(x_min, 1, x_max, 1)
    ratio_unity.SetLineStyle(2)

    hist['ratio_syst'] = hist['ratio_stat'].Clone('ratio_syst')
    hist['ratio_syst'].Sumw2()
    hist['ratio_syst'].SetMarkerSize(0)
    hist['ratio_syst'].SetFillColor(ROOT.kYellow-4)
    hist['ratio_syst'].SetFillStyle(1001)

    for i in range(1, n_xbins+1):
        if (hist['mc_exp'].GetBinContent(i) > 0):
            sq_bin_error = pow(hist['mc_exp'].GetBinError(i), 2)
            sq_bin_error += pow(0.08 * hist['WjLF'].GetBinContent(i), 2)
            sq_bin_error += pow(0.20 * hist['WjHF'].GetBinContent(i), 2)
	    sq_bin_error += pow(0.08 * hist['ZjLF'].GetBinContent(i), 2)
	    sq_bin_error += pow(0.20 * hist['ZjHF'].GetBinContent(i), 2)
	    sq_bin_error += pow(0.07 * hist['TT'].GetBinContent(i), 2)
	    sq_bin_error += pow(0.25 * hist['ST'].GetBinContent(i), 2)
            sq_bin_error += pow(0.25 * hist['VV'].GetBinContent(i), 2)
            sq_bin_error += pow(0.50 * hist['QCD'].GetBinContent(i), 2)

            bin_error = sqrt(sq_bin_error)
            hist['ratio_syst'].SetBinError(i, bin_error / hist['mc_exp'].GetBinContent(i))

    



"""
if __name__ == '__main__':

    step2_dir = '/afs/cern.ch/work/s/swang373/private/V14/'

    categories = {}
    categories['ZH'] = step2_dir + 'ZnnH125.root'
    categories['WH'] = step2_dir + 'WlnH125.root'

    #ROOT.gROOT.SetBatch(1)

    t = TreeKeeper(categories)


    print dir(t)

"""
