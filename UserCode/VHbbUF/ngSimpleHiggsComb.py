from ROOT import TH1F, TFile, RooWorkspace, RooArgList, RooArgSet, RooDataHist, RooHistPdf

histograms = {}
channel = "sim"
var = "var"
nbins, xmin, xmax = 5, 0., 5.
h1 = TH1F("data_obs", "", nbins, xmin, xmax); histograms[h1.GetName()] = h1
h1 = TH1F("ZH", "", nbins, xmin, xmax); histograms[h1.GetName()] = h1
h1 = TH1F("ZjLF", "", nbins, xmin, xmax); histograms[h1.GetName()] = h1
h1 = TH1F("ZjHF", "", nbins, xmin, xmax); histograms[h1.GetName()] = h1
h1 = TH1F("TT", "", nbins, xmin, xmax); histograms[h1.GetName()] = h1
h1 = TH1F("ZH_scale_jUp", "", nbins, xmin, xmax); histograms[h1.GetName()] = h1
h1 = TH1F("ZjLF_scale_jUp", "", nbins, xmin, xmax); histograms[h1.GetName()] = h1
h1 = TH1F("ZjHF_scale_jUp", "", nbins, xmin, xmax); histograms[h1.GetName()] = h1
h1 = TH1F("TT_scale_jUp", "", nbins, xmin, xmax); histograms[h1.GetName()] = h1
h1 = TH1F("ZH_scale_jDown", "", nbins, xmin, xmax); histograms[h1.GetName()] = h1
h1 = TH1F("ZjLF_scale_jDown", "", nbins, xmin, xmax); histograms[h1.GetName()] = h1
h1 = TH1F("ZjHF_scale_jDown", "", nbins, xmin, xmax); histograms[h1.GetName()] = h1
h1 = TH1F("TT_scale_jDown", "", nbins, xmin, xmax); histograms[h1.GetName()] = h1
h1 = TH1F("ZH_statZHUp", "", nbins, xmin, xmax); histograms[h1.GetName()] = h1
h1 = TH1F("ZjLF_statZjLFUp", "", nbins, xmin, xmax); histograms[h1.GetName()] = h1
h1 = TH1F("ZjHF_statZjHFUp", "", nbins, xmin, xmax); histograms[h1.GetName()] = h1
h1 = TH1F("TT_statTTUp", "", nbins, xmin, xmax); histograms[h1.GetName()] = h1
h1 = TH1F("ZH_statZHDown", "", nbins, xmin, xmax); histograms[h1.GetName()] = h1
h1 = TH1F("ZjLF_statZjLFDown", "", nbins, xmin, xmax); histograms[h1.GetName()] = h1
h1 = TH1F("ZjHF_statZjHFDown", "", nbins, xmin, xmax); histograms[h1.GetName()] = h1
h1 = TH1F("TT_statTTDown", "", nbins, xmin, xmax); histograms[h1.GetName()] = h1

histograms["data_obs"].Fill(2, 5)
histograms["ZH"].Fill(2, 1.00)
histograms["ZjLF"].Fill(2, 1.00)
histograms["ZjHF"].Fill(2, 2.00)
histograms["TT"].Fill(2, 1.00)

histograms["ZH_scale_jUp"].Fill(2, histograms["ZH"].Integral()*1.00)
histograms["ZjLF_scale_jUp"].Fill(2, histograms["ZjLF"].Integral()*1.00)
histograms["ZjHF_scale_jUp"].Fill(2, histograms["ZjHF"].Integral()*1.00)
histograms["TT_scale_jUp"].Fill(2, histograms["TT"].Integral()*1.00)
histograms["ZH_scale_jDown"].Fill(2, histograms["ZH"].Integral()*1.10)
histograms["ZjLF_scale_jDown"].Fill(2, histograms["ZjLF"].Integral()*1.00)
histograms["ZjHF_scale_jDown"].Fill(2, histograms["ZjHF"].Integral()*1.00)
histograms["TT_scale_jDown"].Fill(2, histograms["TT"].Integral()*1.00)
histograms["ZH_statZHUp"].Fill(2, histograms["ZH"].Integral()*1.50)
histograms["ZjLF_statZjLFUp"].Fill(2, histograms["ZjLF"].Integral()*1.50)
histograms["ZjHF_statZjHFUp"].Fill(2, histograms["ZjHF"].Integral()*1.50)
histograms["TT_statTTUp"].Fill(2, histograms["TT"].Integral()*1.50)
histograms["ZH_statZHDown"].Fill(2, histograms["ZH"].Integral()*0.50)
histograms["ZjLF_statZjLFDown"].Fill(2, histograms["ZjLF"].Integral()*0.50)
histograms["ZjHF_statZjHFDown"].Fill(2, histograms["ZjHF"].Integral()*0.50)
histograms["TT_statTTDown"].Fill(2, histograms["TT"].Integral()*0.50)
#histograms["ZH_statZHDown"].Fill(2, histograms["ZH"].Integral()/1.50)
#histograms["ZjLF_statZjLFDown"].Fill(2, histograms["ZjLF"].Integral()/1.50)
#histograms["ZjHF_statZjHFDown"].Fill(2, histograms["ZjHF"].Integral()/1.50)
#histograms["TT_statTTDown"].Fill(2, histograms["TT"].Integral()/1.50)

outname = "datacards/simple/simple-shape-experiment.root"
#outfile = TFile.Open(outname, "RECREATE")
ws = RooWorkspace(channel, channel)
ws.factory(var+"[0,5]")
realvar = ws.var(var)
obs = RooArgList(realvar)
for n, h in histograms.iteritems():
    datahist = RooDataHist(n, "", obs, h)
    getattr(ws,'import')(datahist)
ws.writeToFile(outname)
