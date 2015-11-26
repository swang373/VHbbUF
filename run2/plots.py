import math

import ROOT

class ControlRegion(object):

    def __init__(self, name = '', outdir = '', *cuts):

        self.name = name
        self.cuts = cuts
        self.n_cuts = len(cuts)
        self.outfile = ROOT.TFile('{}{}.root'.format(outdir, name), 'recreate')

        print 'Control Region [{}]'.format(self.name)

    def add_tree(self, name = '', ntuple = '', *add_cuts):

        print 'Creating TTree named "{}" from {}'.format(name, ntuple)

        # Access the input file and TTree.
        infile = ROOT.TFile(ntuple, 'read')
        tree = infile.Get('tree')

        # Change directory to the output file.
        self.outfile.cd()
    
        # Perform the control region cuts.
        if (self.n_cuts <= 1):
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
        # Clean up extraneous TTrees.
        self.outfile.cd()
        ROOT.gDirectory.Delete('tree;*')
        self.outfile.Close()

def make_plot(CR_ntuple = '', expression = '', data_weight = '', mc_weight = '', x_title = '', n_xbins = None, x_min = None, x_max = None, filename = ''):
   
    """
    Parameters
    ----------
    CR_ntuple  : str
                 The path to the control region ntuple.
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

    CR_file = ROOT.TFile(CR_ntuple, 'read')

    hist = {}

    #categories: ZH, ggZH, WH, WjLF, WjHF, ZjLF, ZjHF, TT, ST, VV, QCD, data

    # Book histograms and project the TTrees.
    for name in ['Data', 'ZH', 'ggZH', 'WH', 'WjLF', 'WjHF', 'ZjLF', 'ZjHF', 'TT', 'ST', 'VV', 'QCD']:
        hist[name] = ROOT.TH1F('h_' + name, '', n_xbins, x_min, x_max)
        if name == 'Data':
            CR_file.Get(name).Project('h_' + name, expression, data_weight)
        else:
            CR_file.Get(name).Project('h_' + name, expression, mc_weight)
     
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
        elif (h == 'Data'):
            hist[h].SetMarkerSize(0.8)
            hist[h].SetMarkerStyle(20)
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
    
    hist['ratio'] = hist['Data'].Clone('ratio')
    hist['ratio'].Sumw2()
    hist['ratio'].SetMarkerSize(0.8)
    hist['ratio'].Divide(hist['Data'], hist['mc_exp'], 1., 1., '')
    
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
            hist['ratio_stat'].SetBinError(i, bin_error)
        else:
            hist['ratio_stat'].SetBinError(i, 0)

    ratio_unity = ROOT.TLine(x_min, 1, x_max, 1)
    ratio_unity.SetLineStyle(2)
    
    hist['ratio_syst'] = hist['ratio_stat'].Clone('ratio_syst')
    #hist['ratio_syst'].Sumw2()
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

            bin_error = math.sqrt(sq_bin_error)
            hist['ratio_syst'].SetBinError(i, bin_error / hist['mc_exp'].GetBinContent(i))

    # Setup legends
    legend_1 = ROOT.TLegend(0.50, 0.68, 0.72, 0.92)
    legend_1.SetFillColor(0)
    legend_1.SetLineColor(0)
    legend_1.SetShadowColor(0)
    legend_1.SetTextFont(62)
    legend_1.SetTextSize(0.03)
    legend_1.SetBorderSize(1)
    legend_1.AddEntry(hist['Data'], 'Data', 'p')
    legend_1.AddEntry(hist['VH'], 'VH', 'l')
    legend_1.AddEntry(hist['TT'], 't#bar{t}', 'f')
    legend_1.AddEntry(hist['ST'], 'Single Top', 'f')
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
    
    hist['Data'].Draw('e1 same')
    
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
    chi_sq = hist['Data'].Chi2Test(hist['mc_exp'], 'UWCHI2/NDF')
    print chi_sq
    text = pave.AddText('#chi_{{#nu}}^{{2}} = {:.3f}'.format(chi_sq))
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
    
    canvas.IsA().Destructor(canvas)
    
    CR_file.Close()
 
if __name__ == '__main__':

    from settings import *
    import tdrstyle

    ROOT.gROOT.SetBatch(1)

    tdrstyle.set_tdrStyle()

    #categories: ZH, ggZH, WH, WjLF, WjHF, ZjLF, ZjHF, TT, ST, VV, QCD, data
 
    #CR = ControlRegion('CR_Signal_Loose', step2_dir + 'CR/', antiQCD, signal_loose)
    #CR.add_tree('Data', step2_dir + 'Data_MET.root')
    #CR.add_tree('ZH', step2_dir + 'ZnnH125.root')
    #CR.add_tree('ggZH', step2_dir + 'ggZH125.root')
    #CR.add_tree('WH', step2_Dir + 'WlnH125.root')
    #CR.add_tree('WjLF', step2_Dir + 'WJets.root', light_flavour)
    #CR.add_tree('WjHF', step2_Dir + 'WJets.root', heavy_flavour)
    #CR.add_tree('ZjLF', step2_Dir + 'ZJets.root', light_flavour)
    #CR.add_tree('ZjHF', step2_Dir + 'ZJets.root', heavy_flavour)
    #CR.add_tree('TT', step2_Dir + 'TTPow.root')
    #CR.add_tree('ST', step2_Dir + 's_Top.root')
    #CR.add_tree('VV', step2_Dir + 'VV.root')
    #CR.add_tree('QCD', step2_Dir + 'QCD.root')
    #CR.close()

    make_plot(step2_dir + 'CR_Signal_Loose.root', 'HCSV_mass', data_weight, mc_weight, 'm_{jj} [GeV]', 25, 0, 250, 'test')

