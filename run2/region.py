import glob
import logging
import multiprocessing as mp
import os
import subprocess as sp
import sys
import tempfile as tf

import ROOT

from process import Process, PROCESS_DIR
from settings import WORK_DIR, PROCESSES, REGIONS, PLOTS, TARGET_LUMI, DATA_WEIGHT, MC_WEIGHT
from tdrstyle import set_tdrStyle


# Output Directories
REGION_DIR = WORK_DIR + 'regions/'
PLOT_DIR = WORK_DIR + 'plots/'

class Region(object):

    def __init__(self, name = '', cuts = [], **kwargs):
    
        self.logger = logging.getLogger('Region')
        self.logger.info('Initialized for {}'.format(name))

        self.name = name
        self.cuts = cuts

    def make(self):

        # Output Directory
        try:
            os.makedirs(REGION_DIR)
        except OSError:
            if not os.path.isdir(REGION_DIR):
                raise

        # Prepare Processes
        self._check_processes()

        # Parallel Cut
        tasks = mp.Queue()
        results = mp.Queue()

        tmpdir = tf.mkdtemp(prefix = self.name, dir = REGION_DIR)
        self.tmpdir = tmpdir + '/'

        # _processes refers to newly spawned Python processes
        _processes = [
            mp.Process(target = self._cut_process, args = (tasks, results))
            for cpu in xrange(mp.cpu_count())
        ]

        for process in PROCESSES:
            tasks.put(process)

        for p in _processes:
            tasks.put(None)
            p.start()

        for p in _processes:
            p.join()

        results.put(None)

        for r in iter(results.get, None):
            self.logger.info('Selected {!s} out of {!s} entries in {}.'.format(*r))

        # hadd Files
        inputfiles = glob.glob(self.tmpdir + '*.root')
        self.outputfile = REGION_DIR + self.name + '.root'

        sp.check_call(['hadd', '-f', self.outputfile] + inputfiles)
        sp.check_call(['rm', '-r', self.tmpdir])

    def plot(self):

        # Output Directory     
        self.plotdir = PLOT_DIR + self.name + '/'
        try:
            os.makedirs(self.plotdir)
        except OSError:
            if not os.path.isdir(self.plotdir):
                raise

        # Change ROOT global styles to TDR style.
        set_tdrStyle()

        # Set ROOT error verbosity to warnings and higher.
        ROOT.gErrorIgnoreLevel = ROOT.kWarning

        for plot in PLOTS.itervalues():
            self.logger.info('Generating Plot: {name}'.format(**plot))
            self._make_plot(**plot)

    def _check_processes(self):
        
        for process in PROCESSES:
        
            if os.path.isfile(PROCESS_DIR + process + '.root'):
                continue
            
            self.logger.info('Getting missing process {}'.format(process))
            Process(process, **PROCESSES[process]).make()

    
    def _cut_process(self, tasks = None, results = None):

        for process in iter(tasks.get, None):

            infile = ROOT.TFile(PROCESS_DIR + process + '.root', 'read')
            outfile = ROOT.TFile(self.tmpdir + process + '.root', 'recreate')

            intree = infile.Get('tree')
            intree.SetName(process)
            
            for i, cut in enumerate(self.cuts):
               intree.Draw('>>{0!s}_elist_{1!s}'.format(process,i), cut)
               eventlist = ROOT.gDirectory.Get('{0!s}_elist_{1!s}'.format(process,i))
               intree.SetEventList(eventlist)

            outtree = intree.CopyTree('')

            result = (outtree.GetEntriesFast(), intree.GetEntriesFast(), process)

            outtree.Write()
            
            outfile.Close()
            infile.Close()

            results.put(result)

    def _make_plot(self, *args, **kwargs):
        
        # Unpack plot arguments.
        name = kwargs['name']
        expression = kwargs['expression']
        x_title = kwargs['x_title']
        n_bins = kwargs['n_bins']
        x_min = kwargs['x_min']
        x_max = kwargs['x_max']
        logy = kwargs.get('logy', False)

        # Load the custom pileup reweighting function.
        ROOT.gSystem.Load('PU_C')

        # Book histograms in dictionary and stacked histogram.
        infile = ROOT.TFile(self.outputfile, 'read')

        hist = {}
        hstack = ROOT.THStack('hstack','')

        for process, properties in PROCESSES.iteritems():

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

        # Draw upper pad objects.
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

        # Draw lower pad objects.
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

	hist['ratio_stat_unc'].Draw('e2')
	hist['ratio'].Draw('e1 same')
	legend_3.Draw()
	ratio_unity.Draw()
	pave.Draw()

	# Refresh the pads and canvas, then save the plots.
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
            canvas.SaveAs(self.plotdir + name + ftype)
    
        # Clean up to exit gracefully.
        canvas.IsA().Destructor(canvas)    
        infile.Close()

#------
# Main
#------

if __name__ == '__main__':

    # Set ROOT to run in batch mode.
    ROOT.gROOT.SetBatch(1)

    for name in sys.argv[1:]:

        logging.basicConfig(level = logging.INFO,
                            format = '%(name)s(%(levelname)s) - %(message)s')

        region = Region(name, REGIONS[name])
        region.make()
        region.plot()

