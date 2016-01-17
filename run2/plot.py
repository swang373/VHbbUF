import os
import sys

import ROOT

from cut import TARGET_LUMI, DATA_WEIGHT, MC_WEIGHT
from region import REGION_DIR
from settings import WORK_DIR, PROCESSES, PLOTS


# Output Directory
PLOT_DIR = WORK_DIR + 'plots/'

def set_tdrStyle():

    tdrStyle = ROOT.TStyle('tdrStyle', 'Style for P-TDR')
    
    # Canvas
    tdrStyle.SetCanvasBorderMode(0)
    tdrStyle.SetCanvasColor(ROOT.kWhite)
    tdrStyle.SetCanvasDefH(600) # Height
    tdrStyle.SetCanvasDefW(600) # Width
    tdrStyle.SetCanvasDefX(0) # Screen X Position
    tdrStyle.SetCanvasDefY(0) # Screen Y Position
    
    # Pad
    tdrStyle.SetPadBorderMode(0)
    tdrStyle.SetPadColor(ROOT.kWhite)
    tdrStyle.SetPadGridX(False)
    tdrStyle.SetPadGridY(False)
    tdrStyle.SetGridColor(0)
    tdrStyle.SetGridStyle(3)
    tdrStyle.SetGridWidth(1)
    
    # Frame
    tdrStyle.SetFrameBorderMode(0)
    tdrStyle.SetFrameBorderSize(1)
    tdrStyle.SetFrameFillColor(0)
    tdrStyle.SetFrameFillStyle(0)
    tdrStyle.SetFrameLineColor(1)
    tdrStyle.SetFrameLineStyle(1)
    tdrStyle.SetFrameLineWidth(2)
    
    # Histogram
    tdrStyle.SetHistLineColor(1)
    tdrStyle.SetHistLineStyle(0)
    tdrStyle.SetHistLineWidth(1)
    tdrStyle.SetEndErrorSize(1)
    tdrStyle.SetMarkerStyle(20)
    
    # Function/Fit
    tdrStyle.SetOptFit(1)
    tdrStyle.SetFitFormat('5.4g')
    tdrStyle.SetFuncColor(2)
    tdrStyle.SetFuncStyle(1)
    tdrStyle.SetFuncWidth(1)
    
    # Date
    tdrStyle.SetOptDate(0)
    
    # Statistics Box
    tdrStyle.SetOptFile(0)
    tdrStyle.SetOptStat(0) # Use 'mr' to display mean and RMS.
    tdrStyle.SetStatColor(ROOT.kWhite)
    tdrStyle.SetStatFont(42)
    tdrStyle.SetStatFontSize(0.025)
    tdrStyle.SetStatTextColor(1)
    tdrStyle.SetStatFormat('6.4g')
    tdrStyle.SetStatBorderSize(1)
    tdrStyle.SetStatH(0.1)
    tdrStyle.SetStatW(0.15)
    
    # Margins
    tdrStyle.SetPadTopMargin(0.05)
    tdrStyle.SetPadBottomMargin(0.13)
    tdrStyle.SetPadLeftMargin(0.15)
    tdrStyle.SetPadRightMargin(0.03)
    
    # Global Title
    tdrStyle.SetOptTitle(0)
    tdrStyle.SetTitleFont(42)
    tdrStyle.SetTitleColor(1)
    tdrStyle.SetTitleTextColor(1)
    tdrStyle.SetTitleFillColor(10)
    tdrStyle.SetTitleFontSize(0.05)
    
    # Axis Titles
    tdrStyle.SetTitleColor(1, 'XYZ')
    tdrStyle.SetTitleFont(42, 'XYZ')
    tdrStyle.SetTitleSize(0.06, 'XYZ')
    tdrStyle.SetTitleXOffset(0.9)
    tdrStyle.SetTitleYOffset(1.25)
    
    # Axis Labels
    tdrStyle.SetLabelColor(1, 'XYZ')
    tdrStyle.SetLabelFont(42, 'XYZ')
    tdrStyle.SetLabelOffset(0.007, 'XYZ')
    tdrStyle.SetLabelSize(0.05, 'XYZ')
    
    # Axis
    tdrStyle.SetAxisColor(1, 'XYZ')
    tdrStyle.SetStripDecimals(ROOT.kTRUE)
    tdrStyle.SetTickLength(0.03, 'XYZ')
    tdrStyle.SetNdivisions(510, 'XYZ')
    tdrStyle.SetPadTickX(1)  # Tick marks on opposite side of frame.
    tdrStyle.SetPadTickY(1)
    
    # Log Plots
    tdrStyle.SetOptLogx(0)
    tdrStyle.SetOptLogy(0)
    tdrStyle.SetOptLogz(0)
    
    # Postscript
    tdrStyle.SetPaperSize(20.,20.)
    
    tdrStyle.cd()

def make_plot(region = '', 
              name = '', expression = '', x_title = '', 
              n_bins = None, x_min = None, x_max = None, 
              logy = False, **kwargs):

    # Load the custom pileup reweighting function.
    ROOT.gSystem.Load('PU_C')

    # Manually set which processes are omitted. -------------------------------#
    omit = set() # E.g. set(['Wj0b', 'Wj1b', 'Wj2b'])
    #--------------------------------------------------------------------------#

    # Setup output directory.
    outdir = PLOT_DIR + region + '/'
    try:
        os.makedirs(outdir)
    except OSError:
        if not os.path.isdir(outdir):
            raise
      
    # Book histograms in dictionary and stacked histogram.
    infile = ROOT.TFile(REGION_DIR + region + '.root', 'read')

    hist = {}
    hstack = ROOT.THStack('hstack','')

    for process, properties in PROCESSES.iteritems():

        if process in omit:
            continue

        # Check process properties.
        hname = 'h_{}'.format(process)
        types = set(properties['types'].lower().split(':'))
        color = properties.get('color', 0)

        # Project histogram and set style.
        histogram = ROOT.TH1F(hname, '', n_bins, x_min, x_max)

        if 'data' in types:
            infile.Get(process).Project(hname, expression, DATA_WEIGHT)
            histogram.SetMarkerSize(0.8)
            histogram.SetMarkerStyle(20)
        else:
            infile.Get(process).Project(hname, expression, MC_WEIGHT)
            if 'sig' in types:
                histogram.SetFillColor(color)
                histogram.SetMarkerColor(color)
            elif 'bkg' in types:
                histogram.SetLineColor(ROOT.kBlack)
                histogram.SetFillColor(color)
                histogram.SetMarkerColor(color)
            hstack.Add(histogram)

        hist[process] = histogram

    # Create a reference to the sum of all MC samples.
    hist['sum_mc'] = hstack.GetStack().Last().Clone()

    # Rescale the stacked histogram for viewing.
    y_max = max(hist['data_obs'].GetMaximum(), hstack.GetMaximum())
    if logy:
        hstack.SetMaximum(y_max * 17)
    else:
        hstack.SetMaximum(y_max * 1.7)

    # Setup canvas and pads.
    canvas = ROOT.TCanvas('canvas', '', 700, 700)

    upper_pad = ROOT.TPad('upper_pad', '', 0.0, 0.3, 1.0, 1.0)
    upper_pad.SetBottomMargin(0.0)
    upper_pad.SetLogy(logy)
    upper_pad.Draw()

    lower_pad = ROOT.TPad('lower_pad', '', 0.0, 0.0, 1.0, 0.3)
    lower_pad.SetTopMargin(0.0)
    lower_pad.SetBottomMargin(0.35)
    lower_pad.Draw()

    # Setup statistical uncertainty of the MC sum.
    mc_stat_unc = ROOT.TGraphErrors(hist['sum_mc'])
    mc_stat_unc.SetFillColor(ROOT.kGray+3)
    mc_stat_unc.SetFillStyle(3013)
    hist['mc_stat_unc'] = mc_stat_unc
    
    # Setup the data to MC ratio.
    ratio = hist['data_obs'].Clone('ratio')
    ratio.Sumw2()
    ratio.SetMarkerSize(0.8)
    ratio.Divide(hist['data_obs'], hist['sum_mc'], 1., 1., '')
    hist['ratio'] = ratio

    # Setup the unity reference line for the data to MC ratio.
    ratio_unity = ROOT.TLine(x_min, 1, x_max, 1)
    ratio_unity.SetLineStyle(2)
    
    # Setup the statistical uncertainty for the data to MC ratio.
    ratio_stat_unc = hist['data_obs'].Clone('ratio_stat')

    ratio_stat_unc.Sumw2()
    ratio_stat_unc.SetStats(0)
    ratio_stat_unc.SetMaximum(2.2)
    ratio_stat_unc.SetMinimum(0.0)
    ratio_stat_unc.SetMarkerSize(0)
    ratio_stat_unc.SetFillColor(ROOT.kGray+3)
    ratio_stat_unc.SetFillStyle(3013)

    rsu_Xaxis = ratio_stat_unc.GetXaxis()
    rsu_Xaxis.SetTitle(x_title)
    rsu_Xaxis.SetLabelSize(0.12)
    rsu_Xaxis.SetTitleSize(0.14)
    rsu_Xaxis.SetTitleOffset(1.10)

    rsu_Yaxis = ratio_stat_unc.GetYaxis()
    rsu_Yaxis.SetTitle('Data/MC')
    rsu_Yaxis.SetLabelSize(0.10)
    rsu_Yaxis.SetTitleSize(0.12)
    rsu_Yaxis.SetTitleOffset(0.6)
    rsu_Yaxis.SetNdivisions(505)

    for i in xrange(1, n_bins + 1):
        ratio_stat_unc.SetBinContent(i, 1.0)
        if (hist['sum_mc'].GetBinContent(i) > 0):
            bin_error = hist['sum_mc'].GetBinError(i) / hist['sum_mc'].GetBinContent(i)
            ratio_stat_unc.SetBinError(i, bin_error)
        else:
            ratio_stat_unc.SetBinError(i, 0)

    hist['ratio_stat_unc'] = ratio_stat_unc
   
    # Removing systematic uncertainty until we know more about this.
    ''' 
    # Setup the systematic uncertainty for the data to MC ratio.
    ratio_syst_unc = hist['ratio_stat_unc'].Clone('ratio_syst_unc')
    ratio_syst_unc.SetMarkerSize(0)
    ratio_syst_unc.SetFillColor(ROOT.kYellow-4)
    ratio_syst_unc.SetFillStyle(1001)
    
    for i in xrange(1, n_bins + 1):
        if (hist['sum_mc'].GetBinContent(i) > 0): 
            bin_errors = [
                1.00 * hist['sum_mc'].GetBinError(i),
                0.08 * hist['Wj0b'].GetBinError(i),
                0.08 * hist['Wj1b'].GetBinError(i),
                0.20 * hist['Wj2b'].GetBinContent(i),
                0.08 * hist['Zj0b'].GetBinError(i),
                0.08 * hist['Zj1b'].GetBinError(i),
                0.20 * hist['Zj2b'].GetBinContent(i),
                0.07 * hist['TT'].GetBinContent(i),
                0.25 * hist['s_Top'].GetBinContent(i),
                0.25 * hist['VVLF'].GetBinContent(i),
                0.25 * hist['VVHF'].GetBinContent(i),
                0.50 * hist['QCD'].GetBinContent(i),
            ]
            total_bin_error = np.sqrt(np.sum(np.power(bin_errors, 2)))
            ratio_syst_unc.SetBinError(i, total_bin_error / hist['sum_mc'].GetBinContent(i))

    hist['ratio_syst_unc'] = ratio_syst_unc
    '''
        
    # Manually setup drawing of upper pad objects. ----------------------------#
    upper_pad.cd()

    legend_1 = ROOT.TLegend(0.61, 0.62, 0.78, 0.90)
    legend_1.SetFillColor(0)
    legend_1.SetLineColor(0)
    legend_1.SetShadowColor(0)
    legend_1.SetTextFont(62)
    legend_1.SetTextSize(0.03)
    legend_1.SetBorderSize(1)
    legend_1.AddEntry(hist['data_obs'], 'Data', 'p')
    legend_1.AddEntry(hist['ZH'], 'ZH', 'f')
    legend_1.AddEntry(hist['ggZH'], 'ggZH', 'f')
    legend_1.AddEntry(hist['WH'], 'WH', 'f')
    legend_1.AddEntry(hist['VVLF'], 'VV + udscg', 'f')
    legend_1.AddEntry(hist['VVHF'], 'VV + b, b#bar{b}', 'f')
    legend_1.AddEntry(hist['TT'], 't#bar{t}', 'f')
    legend_1.AddEntry(hist['s_Top'], 'Single Top', 'f')

    legend_2 = ROOT.TLegend(0.78, 0.62, 0.95, 0.90)
    legend_2.SetFillColor(0)
    legend_2.SetLineColor(0)
    legend_2.SetShadowColor(0)
    legend_2.SetTextFont(62)
    legend_2.SetTextSize(0.03)
    legend_2.SetBorderSize(1)
    legend_2.AddEntry(hist['Wj0b'], 'W + udscg', 'f')
    legend_2.AddEntry(hist['Wj1b'], 'W + b', 'f')
    legend_2.AddEntry(hist['Wj2b'], 'W + b#bar{b}', 'f')
    legend_2.AddEntry(hist['Zj0b'], 'Z + udscg', 'f')
    legend_2.AddEntry(hist['Zj1b'], 'Z + b', 'f')
    legend_2.AddEntry(hist['Zj2b'], 'Z + b#bar{b}', 'f')
    legend_2.AddEntry(hist['QCD'], 'QCD', 'f') 
    legend_2.AddEntry(hist['mc_stat_unc'], 'MC Unc. (Stat)', 'f')

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAlign(12)
    latex.SetTextFont(62)
    latex.SetTextSize(0.04)

    # Draw Order
    hstack.Draw('hist')
    # Stacked histogram axes exist only after drawing.
    bin_width = (float(x_max) - float(x_min)) / n_bins
    hstack.GetXaxis().SetLabelSize(0)
    hstack.GetYaxis().SetTitle('Events / {:3.3f}'.format(bin_width))
    hist['mc_stat_unc'].Draw('e2 same')
    hist['data_obs'].Draw('e1 same')
    legend_1.Draw()
    legend_2.Draw()
    latex.DrawLatex(0.19, 0.88, 'CMS Preliminary 2016')
    latex.DrawLatex(0.19, 0.83, '#sqrt{{s}} = 13 TeV, L = {:.2f} fb^{{-1}}'.format(TARGET_LUMI/1000.))
    latex.DrawLatex(0.19, 0.78, 'Z(#nu#bar{#nu})H(b#bar{b})')
   
    #hist['ggZH'].SetLineColor(ROOT.kOrange-2)
    #hist['ggZH'].SetLineWidth(3)
    #hist['ggZH'].SetFillColor(0)
    #hist['ggZH'].Draw('hist same')
    #--------------------------------------------------------------------------#
       
    # Manually setup drawing of lower pad objects. ----------------------------#
    lower_pad.cd()
    lower_pad.SetGridy(0)

    legend_3 = ROOT.TLegend(0.78, 0.88, 0.94, 0.95)
    legend_3.SetFillColor(0)
    legend_3.SetLineColor(0)
    legend_3.SetShadowColor(0)
    legend_3.SetTextFont(62)
    legend_3.SetTextSize(0.07)
    legend_3.SetBorderSize(1)
    #legend_3.SetNColumns(2)
    #legend_3.AddEntry(hist['ratio_syst_unc'], 'MC Unc. (Syst)', 'f')
    legend_3.AddEntry(hist['ratio_stat_unc'], 'MC Unc. (Stat)', 'f')
    
    pave = ROOT.TPaveText(0.18, 0.86, 0.28, 0.95, 'brNDC')
    pave.SetLineColor(0)
    pave.SetFillColor(0)
    pave.SetShadowColor(0)
    pave.SetBorderSize(1)

    chi_score = hist['data_obs'].Chi2Test(hist['sum_mc'], 'UWCHI2/NDF')
    text = pave.AddText('#chi_{{#nu}}^{{2}} = {:.3f}'.format(chi_score))
    text.SetTextFont(62)
    text.SetTextSize(0.07)

    # Draw Order
    hist['ratio_stat_unc'].Draw('e2')
    #hist['ratio_syst_unc'].Draw('e2 same')
    hist['ratio'].Draw('e1 same')
    legend_3.Draw()
    ratio_unity.Draw()
    pave.Draw()
    #--------------------------------------------------------------------------#
  
    # Refresh the pads and canvas, then save as plots.
    upper_pad.cd()
    upper_pad.RedrawAxis()
    upper_pad.Modified()
    upper_pad.Update()
    lower_pad.cd()
    lower_pad.RedrawAxis()
    lower_pad.Modified()
    lower_pad.Update()
    canvas.cd()
    
    for ftype in ['.png', '.pdf']:
        canvas.SaveAs(outdir + name + ftype)
    
    # Clean up to exit gracefully.
    canvas.IsA().Destructor(canvas)    
    infile.Close()

#------
# Main
#------
 
if __name__ == '__main__':

    # Set ROOT to run in batch mode.
    ROOT.gROOT.SetBatch(1)

    # Change ROOT global styles to TDR style.
    set_tdrStyle()

    for region in sys.argv[1:]:
 
        for plot in PLOTS.itervalues():
            make_plot(region, **plot)

