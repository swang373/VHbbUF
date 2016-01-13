import math

import ROOT

from settings import *


def make_plot(CR = '', plot = '', expression = '', x_title = '', n_bins = None, x_min = None, x_max = None):
    
    infile = ROOT.TFile('{}{}.root'.format('/afs/cern.ch/work/s/swang373/private/V14/CR/', CR), 'read')

    hist = {}

    # Make histograms for each category.
    for c in PROCESSES:
        hist[c] = ROOT.TH1F('h_{}'.format(c), '', n_bins, x_min, x_max)
        if (c == 'Data_MET'):
            infile.Get(c).Project('h_{}'.format(c), expression, DATA_WEIGHT)
        else:
            infile.Get(c).Project('h_{}'.format(c), expression, MC_WEIGHT)
    
    # Make additional histogram combining all signals.
    hist['VH'] = ROOT.TH1F('VH', '', n_bins, x_min, x_max)
    hist['VH'].Add(hist['ZnnH125'])
    hist['VH'].Add(hist['WlnH125'])

    # Make additional histogram combining all backgrounds.
    hist['mc_exp'] = ROOT.TH1F('mc_exp', '', n_bins, x_min, x_max)
    hist['mc_exp'].Add(hist['ST'])
    hist['mc_exp'].Add(hist['TT'])
    hist['mc_exp'].Add(hist['ZJets'])
    hist['mc_exp'].Add(hist['WJets'])
    hist['mc_exp'].Add(hist['VV'])
    hist['mc_exp'].Add(hist['QCD'])
    
    # Make histogram stack for all signals and backgrounds.
    hstack = ROOT.THStack('hstack','')
    hstack.Add(hist['ST'])
    hstack.Add(hist['TT'])
    hstack.Add(hist['ZJets'])
    hstack.Add(hist['WJets'])
    hstack.Add(hist['VV'])
    hstack.Add(hist['QCD'])
    hstack.Add(hist['VH'])
    hstack.Add(hist['ggZH125'])
    
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

        color = {'WJets': 820, 'ZJets': 5,
                 'TT': 596, 'ST': 840, 'VV': 922, 'QCD': 616}

        if (h == 'VH' or h == 'ZnnH125' or h == 'WlnH125'):
            hist[h].SetFillColor(2)
            hist[h].SetMarkerColor(2)
        elif (h == 'ggZH125'):
            hist[h].SetFillColor(ROOT.kOrange-2)
            hist[h].SetMarkerColor(ROOT.kOrange-2)
        elif (h == 'Data_MET'):
            hist[h].SetMarkerSize(0.8)
            hist[h].SetMarkerStyle(20)
        elif (h in color):
            hist[h].SetLineColor(ROOT.kBlack)
            hist[h].SetFillColor(color[h])
            hist[h].SetMarkerColor(color[h])
        else:
            continue

    # Statistical uncertainty of all backgrounds.
    hist['stat_unc'] = hist['mc_exp'].Clone('stat_unc')
    hist['stat_unc'].Sumw2()
    hist['stat_unc'].SetFillColor(ROOT.kGray+3)
    hist['stat_unc'].SetMarkerSize(0)
    hist['stat_unc'].SetFillStyle(3013)
    
    # Data to MC ratio.
    hist['ratio'] = hist['Data_MET'].Clone('ratio')
    hist['ratio'].Sumw2()
    hist['ratio'].SetMarkerSize(0.8)
    hist['ratio'].Divide(hist['Data_MET'], hist['mc_exp'], 1., 1., '')
    
    # Statistical uncertainty of data to MC ratio.
    hist['ratio_stat'] = hist['mc_exp'].Clone('ratio_stat')
    hist['ratio_stat'].Sumw2()
    hist['ratio_stat'].SetStats(0)
    hist['ratio_stat'].GetXaxis().SetTitle(x_title)
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

    for i in range(1, n_bins + 1):
        hist['ratio_stat'].SetBinContent(i, 1.0)
        if (hist['mc_exp'].GetBinContent(i) > 0):
            bin_error = hist['mc_exp'].GetBinError(i) / hist['mc_exp'].GetBinContent(i)
            hist['ratio_stat'].SetBinError(i, bin_error)
        else:
            hist['ratio_stat'].SetBinError(i, 0)

    # Unity reference line for ratio.
    ratio_unity = ROOT.TLine(x_min, 1, x_max, 1)
    ratio_unity.SetLineStyle(2)
    
    # Systematic uncertainty of data to MC ratio.
    hist['ratio_syst'] = hist['ratio_stat'].Clone('ratio_syst')
    hist['ratio_syst'].SetMarkerSize(0)
    hist['ratio_syst'].SetFillColor(ROOT.kYellow-4)
    hist['ratio_syst'].SetFillStyle(1001)
    
    for i in range(1, n_bins + 1):
        if (hist['mc_exp'].GetBinContent(i) > 0):
	    sq_bin_error = pow(hist['mc_exp'].GetBinError(i), 2)
	    #sq_bin_error += pow(0.08 * hist['WjLF'].GetBinContent(i), 2)
	    sq_bin_error += pow(0.20 * hist['WJets'].GetBinContent(i), 2)
	    #sq_bin_error += pow(0.08 * hist['ZjLF'].GetBinContent(i), 2)
	    sq_bin_error += pow(0.20 * hist['ZJets'].GetBinContent(i), 2)
	    sq_bin_error += pow(0.07 * hist['TT'].GetBinContent(i), 2)
	    sq_bin_error += pow(0.25 * hist['ST'].GetBinContent(i), 2)
	    sq_bin_error += pow(0.25 * hist['VV'].GetBinContent(i), 2)
	    sq_bin_error += pow(0.50 * hist['QCD'].GetBinContent(i), 2)

            bin_error = math.sqrt(sq_bin_error)
            hist['ratio_syst'].SetBinError(i, bin_error / hist['mc_exp'].GetBinContent(i))

    # Setup legends.
    legend_1 = ROOT.TLegend(0.50, 0.68, 0.72, 0.92)
    legend_1.SetFillColor(0)
    legend_1.SetLineColor(0)
    legend_1.SetShadowColor(0)
    legend_1.SetTextFont(62)
    legend_1.SetTextSize(0.03)
    legend_1.SetBorderSize(1)
    legend_1.AddEntry(hist['Data_MET'], 'Data', 'p')
    legend_1.AddEntry(hist['VH'], 'VH', 'l')
    legend_1.AddEntry(hist['ggZH125'], 'ggZH', 'l')
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
    legend_2.AddEntry(hist['WJets'], 'Wj', 'f')
    #legend_2.AddEntry(hist['WjLF'], 'W+LF', 'f')
    legend_2.AddEntry(hist['ZJets'], 'Zj', 'f')
    #legend_2.AddEntry(hist['ZjLF'], 'Z+LF', 'f')
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
    
    # Scale the y-axis for viewing.
    y_max = max(hist['Data_MET'].GetMaximum(), hstack.GetMaximum())
    hstack.SetMaximum(y_max * 1.7) 

    # Draw all the things.
    hstack.Draw('hist')
    hstack.GetXaxis().SetLabelSize(0)
    hstack.GetYaxis().SetTitle('Events / {:3.3f}'.format((float(x_max) - float(x_min)) / float(n_bins)))
    
    hist['stat_unc'].Draw('e2 same')
    
    hist['VH'].SetLineColor(2)
    hist['VH'].SetLineWidth(3)
    hist['VH'].SetFillColor(0)
    hist['VH'].Draw('hist same')

    hist['ggZH125'].SetLineColor(ROOT.kOrange-2)
    hist['ggZH125'].SetLineWidth(3)
    hist['ggZH125'].SetFillColor(0)
    hist['ggZH125'].Draw('hist same')
    
    hist['Data_MET'].Draw('e1 same')
    
    legend_1.Draw()
    legend_2.Draw()
   
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAlign(12)
    latex.SetTextFont(62)
    latex.SetTextSize(0.04)
    latex.DrawLatex(0.19, 0.89, 'CMS Preliminary 2016')
    latex.DrawLatex(0.19, 0.84, '#sqrt{s} = 13 TeV, L = 2.20 fb^{-1}')
    latex.DrawLatex(0.19, 0.79, 'Z(#nu#bar{#nu})H(b#bar{b})')
    
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
    chi_sq = hist['Data_MET'].Chi2Test(hist['mc_exp'], 'UWCHI2/NDF')
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
    
    canvas.SaveAs('{}{}/{}.png'.format(PLOT_DIR, CR, plot))
    canvas.SaveAs('{}{}/{}.pdf'.format(PLOT_DIR, CR, plot))
    
    canvas.IsA().Destructor(canvas)
    
    infile.Close()
 
if __name__ == '__main__':

    import tdrstyle

    # Set ROOT to run in batch mode.
    ROOT.gROOT.SetBatch(1)

    # Change ROOT global styles to TDR style.
    tdrstyle.set_tdrStyle()

    # Create the plot directory if it doesn't exist.
    if (ROOT.gSystem.AccessPathName(PLOT_DIR)):
        ROOT.gSystem.mkdir(PLOT_DIR)
   
    
    CR = 'Signal_Loose'
 
    # Create the control region subdirectory if it doesn't exist.
    if (ROOT.gSystem.AccessPathName(PLOT_DIR + CR)):
        ROOT.gSystem.mkdir(PLOT_DIR + CR)

    for plot, options in PLOTS.iteritems():
        make_plot(CR, plot, **options)

