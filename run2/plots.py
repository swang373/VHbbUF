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
        self.processes[name] = ROOT.TChain('tree')
        self.processes[name].Add(ntuple)
        print 'Added "{}" as process "{}"'.format(ntuple, name)
        
        if selection:
            # Set the TTree which passes the selection to be a named attribute. 
            setattr(self, name, self.processes[name].CopyTree(selection))
            self.trees.append(name)
            print '--- Created TTree {} with {!s} entries'.format(name, getattr(self, name).GetEntries())
        
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
        self.processes[name] = parent
        print 'Added "{}" as subprocess of "{}"'.format(name, parent)

        # Set the TTree which passes the selection to be a named attribute.
        setattr(self, name, self.processes[parent].CopyTree(selection))
        self.trees.append(name)
        print '--- Created TTree {} with {!s} entries'.format(name, getattr(self, name).GetEntries())
            

def make_plot(tchainsaw = None, expression = '', x_title = '', n_xbins = None, x_min = None, x_max = None, filename = ''):
   
    """
    Parameters
    ----------
    tchainsaw  : TChainSaw
                 The object managing the selected TTrees.
    expression : str
                 A TTreeFormula style string defining the expression to be passed to
                 TTree::Project(). Refer to TTree::Draw() documentation for examples.
    x_title    : str
                 The title of the x-axis.
    n_xbins    : int
                 The number of bins along the x-axis.
    x_min      : float
                 The lower bound of the x-axis.
    x_max      : float
                 The upper bound of the x-axis.
    filename   : str
                 The name used to save the plot.
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

        if (h == 'VH' or h == 'ZH' or h == 'WH'):
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

    # Setup legends
    legend_1 = ROOT.TLegend(0.50, 0.68, 0.72, 0.92)
    legend_1.SetFillColor(0)
    legend_1.SetLineColor(0)
    legend_1.SetShadowColor(0)
    legend_1.SetTextFont(62)
    legend_1.SetTextSize(0.03)
    legend_1.SetBorderSize(1)
    legend_1.AddEntry(hist['data'], 'Data', 'p')
    legend_1.AddEntry(hist['VH'], 'VH', 'l')
    legend_1.AddEntry(hist['TT'], 't#bar{T}', 'f')
    legend_1.AddEntry(hist['VV'], 'VV', 'f')
   
    legend_2 = ROOT.TLegend(0.72, 0.68, 0.94, 0.92)
    legend_2.SetFillColor(0)
    legend_2.SetLineColor(0)
    legend_2.SetShadowColor(0)
    legend_2.SetTextFont(62)
    legend_2.SetTextSize(0.03)
    legend_2.SetBorderSize(1)
    legend_2.AddEntry(hist['WjHF'], 'W+HF', 'f')
    legend_2.AddEntry(hist['WjLF'], 'W+LF', 'f')
    legend_2.AddEntry(hist['ZjHF'], 'Z+HF', 'f')
    legend_2.AddEntry(hist['ZjLF'], 'Z+LF', 'f')
    legend_2.AddEntry(hist['QCD'], 'QCD', 'f') 
    legend_2.AddEntry(hist['stat_unc'], 'MC Unc. (Stat)', 'f')
    
    ratio_legend_1 = ROOT.TLegend(0.72, 0.88, 0.94, 0.96)
    ratio_legend_1.SetFillColor(0)
    ratio_legend_1.SetLineColor(0)
    ratio_legend_1.SetShadowColor(0)
    ratio_legend_1.SetTextFont(62)
    ratio_legend_1.SetTextSize(0.07)
    ratio_legend_1.SetBorderSize(1)
    ratio_legend_1.AddEntry(hist['ratio_stat'], 'MC Unc. (Stat)', 'f')
    
    ratio_legend_2 = ROOT.TLegend(0.50, 0.88, 0.72, 0.96)
    ratio_legend_2.SetFillColor(0)
    ratio_legend_2.SetLineColor(0)
    ratio_legend_2.SetShadowColor(0)
    ratio_legend_2.SetTextFont(62)
    ratio_legend_2.SetTextSize(0.07)
    ratio_legend_2.SetBorderSize(1)
    ratio_legend_2.AddEntry(hist['ratio_syst'], 'MC Unc. (Syst)', 'f')

    # Draw stuff.
    hstack.Draw('hist')
    hstack.GetXaxis().SetLabelSize(0)
    hstack.GetYaxis().SetTitle('Events / {:.3f}'.format((x_max - x_min) / n_xbins))

    hist['stat_unc'].Draw('e2 same')

    hist['VH'].SetLineColor(2)
    hist['VH'].SetLineWidth(3)
    hist['VH'].SetFillColor(0)
    hist['VH'].Draw('hist same')

    hist['data'].Draw('e1 same')

    legend_1.Draw()
    legend_2.Draw()

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAlign(12)
    latex.SetTextFont(62)
    latex.SetTextSize(0.04)
    latex.DrawLatex(0.19, 0.89, 'CMS Simulation 2015')
    latex.DrawLatex(0.19, 0.84, '#sqrt{s} = 13 TeV, L = 1.28 fb^{-1}')
    latex.DrawLatex(0.19, 0.79, 'Z#nu#bar{#nu}Hb#bar{b}')

    lower_pad.cd()
    lower_pad.SetGridy(0)
    hist['ratio_stat'].Draw('e2')
    hist['ratio_syst'].Draw('e2 same')
    hist['ratio_stat'].Draw('e2 same')

    ratio_unity.Draw()
    hist['ratio'].Draw('e1 same')
    ratio_legend_1.Draw()
    ratio_legend_2.Draw()

    pave = ROOT.TPaveText(0.18, 0.86, 0.28, 0.96, 'brNDC')
    pave.SetLineColor(0)
    pave.SetFillColor(0)
    pave.SetShadowColor(0)
    pave.SetBorderSize(1)
    chi_sq = hist['data'].Chi2Test(hist['mc_exp'], 'UWCHI2/NDF')
    text = pave.AddText('#chi_{#nu}^{2} = {:.3f}'.format(chi_sq))
    text.SetTextFont(62)
    text.SetTextSize(0.07)
    pave.Draw()

    upper_pad.cd()
    upper_pad.RedrawAxis()
    upper_pad.Modified()
    upper_pad.Update()
    lower_pad.cd()
    lower_pad.RedrawAxis()
    lower_pad.Modified()
    lower_pad.Update()
    canvas.cd()


    canvas.SaveAs('plots/{}.png'.format(filename))
    canvas.SaveAs('plots/{}.pdf'.format(filename))


if __name__ == '__main__':

    from settings import *
    import TCutOperators as tco

    ROOT.gROOT.SetBatch(1)

    step2_dir = '/afs/cern.ch/work/s/swang373/private/V14/'

    #categories: ZH, ggZH, WH, WjLF, WjHF, ZjLF, ZjHF, TT, ST, VV, QCD, data

    print tco.add(antiQCD, signal_loose)

    test = TChainSaw()
    test.add_process('data', step2_dir + 'Data_MET.root', tco.add(antiQCD, signal_loose))
    test.add_process('ZH', step2_dir + 'ZnnH125.root', tco.add(antiQCD, signal_loose))
    test.add_process('ggZH', step2_dir + 'ggZH125.root', tco.add(antiQCD, signal_loose))
    test.add_process('WH', step2_dir + 'WlnH125.root', tco.add(antiQCD, signal_loose))
    test.add_process('WJets', step2_dir + 'WJets.root')
    test.add_subprocess('WjLF', 'WJets', tco.add(antiQCD, signal_loose, light_flavour))
    test.add_subprocess('WjHF', 'WJets', tco.add(antiQCD, signal_loose, heavy_flavour))
    test.add_process('ZJets', step2_dir + 'ZJets.root')
    test.add_subprocess('ZjLF', 'ZJets', tco.add(antiQCD, signal_loose, light_flavour))
    test.add_subprocess('ZjHF', 'ZJets', tco.add(antiQCD, signal_loose, heavy_flavour))
    test.add_process('TT', step2_dir + 'TTPow.root', tco.add(antiQCD, signal_loose))
    test.add_process('ST', step2_dir + 's_Top.root', tco.add(antiQCD, signal_loose))
    test.add_process('VV', step2_dir + 'VV.root', tco.add(antiQCD, signal_loose))
    test.add_process('QCD', step2_dir + 'QCD.root', tco.add(antiQCD, signal_loose))

    print test.processes
    print test.trees
    print dir(test)
