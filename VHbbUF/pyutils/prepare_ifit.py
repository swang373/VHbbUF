from ROOT import gROOT, TFile, TH1, TH1F
import os


def draw(indir, key, var):
    histograms = {}
    
    # Ntuples to loop over
    processes = ["Zj", "Wj", "TT", "s_Top", "VV", "data_obs"]
    
    # Fit regions to loop over
    categories = {"ctrl1":1, "ctrl2":2, "ctrl3":3, "ctrl4":4, "ctrl5":5}
    
    for c, cv in categories.iteritems():
        for p in processes:
            tfile = TFile.Open(indir + "Step4_" + p + ".root")
            
            if p in ["Zj","Wj"]:
                h1 = TH1F(key % (c, p+"HF"), "; H p_{T} [GeV]", 100, 0, 500)
                tfile.tree_ZnunuHighPt_ctrl.Project(h1.GetName(), var, "weightsMC[0]*(selectFlags[%i][0] && eventFlav == 5)" % (cv), "goff")
                histograms[h1.GetName()] = h1
                
                h1 = TH1F(key % (c, p+"LF"), "; H p_{T} [GeV]", 100, 0, 500)
                tfile.tree_ZnunuHighPt_ctrl.Project(h1.GetName(), var, "weightsMC[0]*(selectFlags[%i][0] && eventFlav != 5)" % (cv), "goff")
                histograms[h1.GetName()] = h1
            
            else:
                h1 = TH1F(key % (c, p), "; H p_{T} [GeV]", 100, 0, 500)
                tfile.tree_ZnunuHighPt_ctrl.Project(h1.GetName(), var, "weightsMC[0]*(selectFlags[%i][0])" % (cv), "goff")
                histograms[h1.GetName()] = h1
    return histograms


def write(histograms, outfile):
    tfile = TFile.Open(outfile, "RECREATE")
    for h, hv in histograms.iteritems():
        hv.Write()
    tfile.Close()

##______________________________________________________________________________
if __name__ == '__main__':
    
    indir = "Step4_20130404/stitch/"
    outdir = "results/"
    outfile = "ifit_HptReg.root"
    key = "HptReg_%s_%s"
    var = "HptReg"

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Make histograms
    histograms = draw(indir, key, var)
    
    # Write histograms
    write(histograms, outdir + outfile)

