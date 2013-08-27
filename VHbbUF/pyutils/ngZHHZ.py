from ROOT import gROOT, gPad, AddressOf

gROOT.LoadMacro("HelperFunctions.h")
gROOT.LoadMacro("HelperVHbbDataFormats.h")
gROOT.LoadMacro("tdrstyle.C")
gROOT.ProcessLine("setTDRStyle()")

from ROOT import TFile, TTree, TTreeFormula, TF1, TH1F, TLatex, TLegend, TColor, TGraphErrors, TCanvas, TPad, TLatex, TLine
from math import sqrt
#gROOT.SetBatch(1)
#kGreen = 416
#kRed = 632
colZH=2
colSystZH=TColor.GetColor("#FF9A9A")
#colHZ=896
colHZ=418 # kGreen+2

fileZH = "ZH125"
fileHZ = "ZbbHinv125"

doTT=False
if doTT:
    colZH=600
    colSystZH=TColor.GetColor("#9A9AFF")
    fileZH="TT"

# ------------------------------------------------------------------------------

fZH = TFile.Open("Step4_20130404/stitch/Step4_"+fileZH+".root")
fHZ = TFile.Open("Step4_20130404/stitch/Step4_"+fileHZ+".root")

ptbins = [
    "ZnunuHighPt", 
    "ZnunuMedPt", 
    "ZnunuLowPt"
    ]
variables = [
    "HptReg;p_{T}(jj) [GeV];15,130,355;15,130,310;15,100,250",
    "METtype1corr.et;E_{T}^{miss} [GeV];12,170,350;8,130,170;6,100,130",
    "max(hJet_ptReg[0],hJet_ptReg[1]);p_{T}(j_{1}) [GeV];15,60,285;15,60,240;15,60,210",
    "min(hJet_ptReg[0],hJet_ptReg[1]);p_{T}(j_{2}) [GeV];13,30,160;13,30,160;13,30,134",
    "max(hJet_csv_nominal[0],hJet_csv_nominal[1]);CSV_{max};15,0.0,1.05;15,0.0,1.05;15,0.0,1.05",
    "min(hJet_csv_nominal[0],hJet_csv_nominal[1]);CSV_{min};15,0.0,1.05;15,0.0,1.05;15,0.0,1.05",
    ]



# Setup canvas and pads
c1 = TCanvas("c1", "c1", 700, 700)
#c1 = TCanvas("c1", "c1", 600, 600)

pad1 = TPad("pad1", "top pad"   , 0.0, 0.3, 1.0, 1.0)
pad1.SetBottomMargin(0.0)
pad1.Draw()
pad2 = TPad("pad2", "bottom pad", 0.0, 0.0, 1.0, 0.3)
pad2.SetTopMargin(0.0)
pad2.SetBottomMargin(0.35)
pad2.Draw()
pad1.cd()

def setRatioStyle(ratio):
    #ratio.Sumw2()
    ratio.SetStats(0)
    ratio.SetTitle("")
    #ratio.GetXaxis().SetTitle(xtitle)
    #ratio.GetYaxis().SetTitle("Data/MC")
    ratio.GetYaxis().SetTitle("Z(bb)H(inv)/Z(#nu#nu)H(bb)")
    #ratio.GetYaxis().SetTitle("Z(bb)H(inv)/t#bar{t}")
    ratio.SetMaximum(2.2)
    ratio.SetMinimum(0)
    ratio.SetMarkerSize(0)
    ratio.SetFillColor(923)  # kGray+3
    ratio.SetFillStyle(3013)
    ratio.GetXaxis().CenterTitle()
    ratio.GetXaxis().SetLabelSize(0.12)
    ratio.GetXaxis().SetTitleSize(0.14)
    ratio.GetXaxis().SetTitleOffset(1.10)
    ratio.GetYaxis().CenterTitle()
    ratio.GetYaxis().SetLabelSize(0.10)
    #ratio.GetYaxis().SetTitleSize(0.12)
    ratio.GetYaxis().SetTitleSize(0.10)
    ratio.GetYaxis().SetTitleOffset(0.6)
    ratio.GetYaxis().SetNdivisions(505)

#lumi = 19624.0;
#lumiZH = lumi *       0.045502 /  1001196.6875;
#lumiHZ = lumi *       0.059618 /   500025.3438;

it = 0
for ptbin in ptbins:
    tZH = fZH.Get("tree_%s_train" % ptbin)
    tHZ = fHZ.Get("tree_%s_train" % ptbin)

    for variable in variables:
        #if doTT and "CSV" not in variable:  continue
        
        tup = variable.split(";")
        ttup = tup[2+it/len(variables)].split(",")
        hZH = TH1F("hZH_%i" % it, ";%s;normalized" % tup[1], int(ttup[0]), float(ttup[1]), float(ttup[2]))
        hHZ = TH1F("hHZ_%i" % it, ";%s;normalized" % tup[1], int(ttup[0]), float(ttup[1]), float(ttup[2]))
        
        hZH.Sumw2()
        hHZ.Sumw2()
        
        hZH.SetLineColor(colZH)
        hZH.SetLineWidth(2)
        hZH.SetMarkerSize(0)
        
        hHZ.SetLineColor(colHZ)
        hHZ.SetLineWidth(2)
        hHZ.SetMarkerSize(0)
        
        exprweightZH = "2.0 * weightsMC[0] * (selectFlags[0][0]) * 19040/19624 * weightSignalEWKNew * weightSignalQCD"
        exprweightHZ = "2.0 * weightsMC[0] * (selectFlags[0][0]) * 19040/19624 * weightSignalEWK * 176432/160481 * weightSignalQCD"
        if doTT:
            exprweightZH = "2.0 * weightsMC[0] * (selectFlags[0][0]) * 19040/19624"
        
        tZH.Project(hZH.GetName(), "min("+tup[0]+",%.1f)" % float(ttup[2]), exprweightZH)
        tHZ.Project(hHZ.GetName(), "min("+tup[0]+",%.1f)" % float(ttup[2]), exprweightHZ)
        
        # Normalize
        normZH = hZH.GetSumOfWeights()
        normHZ = hHZ.GetSumOfWeights()
        hZH.Scale(1.0/normZH)
        hHZ.Scale(1.0/normHZ)
        hZH.SetMaximum(max(hZH.GetMaximum(), hHZ.GetMaximum())*1.25)
        
        # Ratio
        hratio = hHZ.Clone("hratio_%i" % it)
        hratio.Divide(hHZ, hZH, 1., 1., "")
        setRatioStyle(hratio)
        y_one_line = TLine(float(ttup[1]), 1, float(ttup[2]), 1)
        y_one_line.SetLineStyle(2)
        
        # Draw
        pad1.cd()
        hZH.Draw("e1")
        hZH.Draw("samehist")
        hHZ.Draw("same e1")
        hHZ.Draw("samehist")
        
        pad2.cd()
        hratio.Draw("e1")
        y_one_line.Draw()
        
        # Systematics
        if "CSV" in variable:
            hupZH1 = TH1F("hZH_upBC_%i" % it, ";%s;normalized" % tup[1], int(ttup[0]), float(ttup[1]), float(ttup[2]))
            hupHZ1 = TH1F("hHZ_upBC_%i" % it, ";%s;normalized" % tup[1], int(ttup[0]), float(ttup[1]), float(ttup[2]))
            hdnZH1 = TH1F("hZH_downBC_%i" % it, ";%s;normalized" % tup[1], int(ttup[0]), float(ttup[1]), float(ttup[2]))
            hdnHZ1 = TH1F("hHZ_downBC_%i" % it, ";%s;normalized" % tup[1], int(ttup[0]), float(ttup[1]), float(ttup[2]))
            
            hupZH2 = TH1F("hZH_upL_%i" % it, ";%s;normalized" % tup[1], int(ttup[0]), float(ttup[1]), float(ttup[2]))
            hupHZ2 = TH1F("hHZ_upL_%i" % it, ";%s;normalized" % tup[1], int(ttup[0]), float(ttup[1]), float(ttup[2]))
            hdnZH2 = TH1F("hZH_downL_%i" % it, ";%s;normalized" % tup[1], int(ttup[0]), float(ttup[1]), float(ttup[2]))
            hdnHZ2 = TH1F("hHZ_downL_%i" % it, ";%s;normalized" % tup[1], int(ttup[0]), float(ttup[1]), float(ttup[2]))
            
            hupZH1.Sumw2()
            hupHZ1.Sumw2()
            hdnZH1.Sumw2()
            hdnHZ1.Sumw2()
            hupZH2.Sumw2()
            hupHZ2.Sumw2()
            hdnZH2.Sumw2()
            hdnHZ2.Sumw2()
            
            tZH.Project(hupZH1.GetName(), "min("+tup[0].replace("csv_nominal","csv_upBC")+",%.1f)" % float(ttup[2]), exprweightZH.replace("selectFlags[0][0]","selectFlags[0][5]"))
            tHZ.Project(hupHZ1.GetName(), "min("+tup[0].replace("csv_nominal","csv_upBC")+",%.1f)" % float(ttup[2]), exprweightHZ.replace("selectFlags[0][0]","selectFlags[0][5]"))
            tZH.Project(hdnZH1.GetName(), "min("+tup[0].replace("csv_nominal","csv_downBC")+",%.1f)" % float(ttup[2]), exprweightZH.replace("selectFlags[0][0]","selectFlags[0][6]"))
            tHZ.Project(hdnHZ1.GetName(), "min("+tup[0].replace("csv_nominal","csv_downBC")+",%.1f)" % float(ttup[2]), exprweightHZ.replace("selectFlags[0][0]","selectFlags[0][6]"))
            
            tZH.Project(hupZH2.GetName(), "min("+tup[0].replace("csv_nominal","csv_upL")+",%.1f)" % float(ttup[2]), exprweightZH.replace("selectFlags[0][0]","selectFlags[0][7]"))
            tHZ.Project(hupHZ2.GetName(), "min("+tup[0].replace("csv_nominal","csv_upL")+",%.1f)" % float(ttup[2]), exprweightHZ.replace("selectFlags[0][0]","selectFlags[0][7]"))
            tZH.Project(hdnZH2.GetName(), "min("+tup[0].replace("csv_nominal","csv_downL")+",%.1f)" % float(ttup[2]), exprweightZH.replace("selectFlags[0][0]","selectFlags[0][8]"))
            tHZ.Project(hdnHZ2.GetName(), "min("+tup[0].replace("csv_nominal","csv_downL")+",%.1f)" % float(ttup[2]), exprweightHZ.replace("selectFlags[0][0]","selectFlags[0][8]"))
            
            # Normalize
            hupZH1.Scale(1.0/normZH)
            hupHZ1.Scale(1.0/normHZ)
            hdnZH1.Scale(1.0/normZH)
            hdnHZ1.Scale(1.0/normHZ)
            hupZH2.Scale(1.0/normZH)
            hupHZ2.Scale(1.0/normHZ)
            hdnZH2.Scale(1.0/normZH)
            hdnHZ2.Scale(1.0/normHZ)
            
            hsystZH = hZH.Clone("hsystZH_%i" % it)
            hsystHZ = hHZ.Clone("hsystHZ_%i" % it)
            
            #hsystZH.Sumw2()
            #hsystHZ.Sumw2()
            
            #hsystZH.SetFillStyle(3013)
            #hsystZH.SetFillColor(colZH)
            hsystZH.SetFillStyle(1001)
            hsystZH.SetFillColor(colSystZH)
            
            hsystHZ.SetFillStyle(3013)
            hsystHZ.SetFillColor(colHZ)
            
            for b in xrange(0,hZH.GetNbinsX()+1):
                systZH = max(abs(hupZH1.GetBinContent(b)-hZH.GetBinContent(b)),abs(hdnZH1.GetBinContent(b)-hZH.GetBinContent(b)))
                systHZ = max(abs(hupHZ1.GetBinContent(b)-hHZ.GetBinContent(b)),abs(hdnHZ1.GetBinContent(b)-hHZ.GetBinContent(b)))
                hsystZH.SetBinError(b, sqrt(hsystZH.GetBinError(b)**2 + systZH**2))
                hsystHZ.SetBinError(b, sqrt(hsystHZ.GetBinError(b)**2 + systHZ**2))
                
                systZH = max(abs(hupZH2.GetBinContent(b)-hZH.GetBinContent(b)),abs(hdnZH2.GetBinContent(b)-hZH.GetBinContent(b)))
                systHZ = max(abs(hupHZ2.GetBinContent(b)-hHZ.GetBinContent(b)),abs(hdnHZ2.GetBinContent(b)-hHZ.GetBinContent(b)))
                hsystZH.SetBinError(b, sqrt(hsystZH.GetBinError(b)**2 + systZH**2))
                hsystHZ.SetBinError(b, sqrt(hsystHZ.GetBinError(b)**2 + systHZ**2))
            
            hsystratio = hsystHZ.Clone("hsystratio_%i" % it)
            hsystratio.Divide(hsystHZ, hsystZH, 1., 1., "")
            
            # Draw
            pad1.cd()
            hsystZH.Draw("same e2")
            hsystHZ.Draw("same e2")
            hZH.Draw("samehist")
            hZH.Draw("same e1")
            hHZ.Draw("samehist")
            hHZ.Draw("same e1")
            
            pad2.cd()
            hsystratio.Draw("same e2")
            hratio.Draw("same e1")
        
        
        # Stamp
        pad1.cd()
        latex = TLatex()
        latex.SetNDC()
        latex.SetTextAlign(12)
        latex.SetTextFont(42)
        latex.SetTextSize(0.032)
        latex.DrawLatex(0.16, 0.97, "CMS Simulation  #sqrt{s} = 8 TeV")
        #latex.SetTextSize(0.045)
        latex.SetTextSize(0.038)
        if not doTT:
            latex.DrawLatex(0.19, 0.90, "#color[%i]{Z(#rightarrow E_{T}^{miss})H(#rightarrow b #bar{b})}" % colZH)
        else:
            latex.DrawLatex(0.19, 0.90, "#color[%i]{t#bar{t}(#rightarrow jets+E_{T}^{miss})}" % colZH)
        latex.DrawLatex(0.19, 0.85, "#color[%i]{Z(#rightarrow b#bar{b}) H(#rightarrow E_{T}^{miss})}" % colHZ)

        
        varname = tup[1]
        varname = varname.replace("(","")
        varname = varname.replace(")","")
        varname = varname.replace("{","")
        varname = varname.replace("}","")
        varname = varname.replace("[","")
        varname = varname.replace("]","")
        varname = varname.replace("_","")
        varname = varname.replace("^","")
        varname = varname.replace(" ","")
        varname = varname.replace("GeV","")
        if not doTT:
            c1.Print("plots/zhhz_%s_%s.png" % (ptbin,varname))
            c1.Print("plots/zhhz_%s_%s.pdf" % (ptbin,varname))
        else:
            c1.Print("plots/zhhz_TT_%s_%s.png" % (ptbin,varname))
            c1.Print("plots/zhhz_TT_%s_%s.pdf" % (ptbin,varname))
        it += 1
    
