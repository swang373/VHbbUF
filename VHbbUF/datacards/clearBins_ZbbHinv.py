from ROOT import TFile

channels  = ["ZnunuHighPt", "ZnunuMedPt", "ZnunuLowPt"]
firstbins = [13, 18, 19]
lastbins  = [20, 20, 20]
rootfilename = "zhinv_Zbb_J14_$CHANNEL_TH1_test.root"
newrootfilename = "zhinv_Zbb_J14_$CHANNEL_TH1.root"


for ic, channel in enumerate(channels):
    rootfile = TFile.Open(rootfilename.replace("$CHANNEL", channel))
    directory = rootfile.Get(channel)

    histograms = []
    for key in directory.GetListOfKeys():
        if key.GetClassName() == "TH1F":
            h = directory.Get(key.GetName())
            h.SetName(key.GetName())
            print h.GetName()
            #if "QCD" in h.GetName(): print "r", h.GetName()
            for b in xrange(1, h.GetNbinsX()+1):
                if b < firstbins[ic] or b > lastbins[ic]:
                    h.SetBinContent(b, 0)
                    h.SetBinError(b, 0)
            histograms.append(h)

    newrootfile = TFile.Open(newrootfilename.replace("$CHANNEL", channel), "RECREATE")
    newrootfile.mkdir(channel)
    newrootfile.cd(channel)
    for h in histograms:
        #print "w", h.GetName()
        h.Write()
    newrootfile.Close()

    rootfile.Close()
