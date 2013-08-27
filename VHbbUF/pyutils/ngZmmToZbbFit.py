from ROOT import TFile, TH1F, TH1, TLatex, TCanvas, TLine, gPad, gROOT, gStyle, kDashed, RooRealVar, RooArgList, RooArgSet, RooAbsData, RooAbsReal, RooDataHist, RooDataSet, RooWorkspace, RooPlot, RooFit
from math import sqrt
gROOT.ProcessLine(".L roofit/RooEMGaussian.cxx+")
gROOT.ProcessLine(".L roofit/RooGaussExp.cxx+")
gROOT.ProcessLine(".L roofit/RooExpGaussExp.cxx+")
from ROOT import RooEMGaussian, RooGaussExp, RooExpGaussExp
import sys


infilename = "roofit/zmmtozbb_jetres.root"
ih_ = 3
if len(sys.argv) == 2: 
    ih_ = int(sys.argv[1])
    gROOT.SetBatch(1)


# 0 = double gaussian
# 1 = exponentially modified gaussian
# 2 = gaussian + expo tail
# 3 = gaussian + expo tail x 2
# 4 = bifurcated gaussian
ifunc = 3

xlow, xup = 0.0, 3.0    # x range
fxlow, fxup = 0.4, 1.5  # fit range
#fxlow, fxup = 0.4, 1.5  # fit range
gmlow, gmup = 0.5, 1.2 # gaussian mean range
gslow, gsup = 0.0, 1.0  # gaussian sigma range
lamblow, lambup = 1.0, 10.0  # exponential lambda range
kapplow, kappup = -10.0, -0.1  # exponential kappa range
logy = False


gROOT.LoadMacro("tdrstyle.C")
gROOT.ProcessLine("setTDRStyle()")
#TH1.SetDefaultSumw2()

gStyle.SetMarkerSize(0.7)
kBlue = 600

latex = TLatex()
latex.SetNDC()
latex.SetTextAlign(12)
latex.SetTextFont(42)
latex.SetTextSize(0.035)

ws = RooWorkspace("ws")
x = RooRealVar("x", "p_{T}^{reco}/p_{T}^{parton}", xlow,xup)
obs = RooArgList(x)
c1 = TCanvas("c1", "c1",600,600) ;
if logy:  c1.SetLogy()


# Loop over input files (currently disabled)
infile = TFile.Open(infilename)
for ih in xrange(ih_,ih_+1):
    h = getattr(infile, "hh1_%i" % ih)
    data = RooDataHist("data_%i" %ih, "", obs, h)
    getattr(ws, 'import')(data)

    hl = getattr(infile, "hl_%i" % ih)
    
    frame = x.frame()


    if ifunc == 0:
        # Fit double gaussian
        ws.factory("SUM::double_gaussian(Gaussian::gaus1(x,gm1[1.0,%.1f,%.1f],gs1[0.1,%.1f,%.1f]), " 
                                    "f[0.5,0,1]*Gaussian::gaus2(x,gm2[1.0,%.1f,%.1f],gs2[0.3,%.1f,%.1f]) )"
                    % (0.9,gmup, gslow/2,gsup/2, 0.8,gmup, gslow,gsup))
        model = ws.pdf("double_gaussian")
        
        #result = model.fitTo(data, RooFit.Save())
        result = model.fitTo(data, RooFit.Range(fxlow,fxup), RooFit.Save())
        
        data.plotOn(frame,RooFit.Name("data"),RooFit.DataError(RooAbsData.SumW2))
        model.plotOn(frame,RooFit.Components("gaus1"),RooFit.LineStyle(kDashed),RooFit.LineColor(kBlue-7))
        model.plotOn(frame,RooFit.Components("gaus2"),RooFit.LineStyle(kDashed),RooFit.LineColor(kBlue-7))
        model.plotOn(frame,RooFit.Name("model"))
        frame.SetMaximum(h.GetMaximum() if logy else 6000)
        frame.SetMinimum(0.5 if logy else 0)
        frame.Draw()
        
        chi2 = frame.chiSquare("model", "data", 5)
        p1 = ws.var("gm1")
        p2 = ws.var("gs1")
        p3 = ws.var("gm2")
        p4 = ws.var("gs2")
        p5 = ws.var("f")
        result.floatParsFinal().Print("s")
        print "chi^2 = ", chi2
        
        latex.DrawLatex(0.19, 0.90, "%.0f<p_{T}^{Z}<%.0f, %.0f<p_{T}^{j}<%.0f, %.1f<|#eta^{j}|<%.1f" % tuple([hl.GetBinContent(l) for l in xrange(1,7)]))
        latex.DrawLatex(0.19, 0.85, "#color[1]{FullSim}")
        latex.DrawLatex(0.19, 0.80, "#color[4]{Fit double gaus}")
        latex.DrawLatex(0.42, 0.85, "#color[1]{#mu=%.3f, #sigma=%.3f}" % (h.GetMean(), h.GetRMS()))
        #latex.DrawLatex(0.42, 0.80, "#color[4]{#mu_{1}=%.3f #pm %.3f, #sigma_{1}=%.3f #pm %.3f}" % (p1.getVal(), p1.getError(), p2.getVal(), p2.getError() ))
        #latex.DrawLatex(0.42, 0.76, "#color[4]{#mu_{2}=%.3f #pm %.3f, #sigma_{2}=%.3f #pm %.3f}" % (p3.getVal(), p3.getError(), p4.getVal(), p4.getError() ))
        #latex.DrawLatex(0.68, 0.72, "#color[4]{f_{2/1}=%.3f #pm %.3f}" % (p5.getVal(), p5.getError() ))
        latex.DrawLatex(0.42, 0.80, "#color[4]{#mu_{1}=%.3f, #sigma_{1}=%.3f}" % (p1.getVal(), p2.getVal() ))
        latex.DrawLatex(0.42, 0.76, "#color[4]{#mu_{2}=%.3f, #sigma_{2}=%.3f, f_{2/1}=%.3f}" % (p3.getVal(), p4.getVal(), p5.getVal() ))
        latex.DrawLatex(0.19, 0.70, "#color[2]{#chi_{#nu}^{2}=%.2f}" % (chi2))
        gPad.RedrawAxis()
        gPad.Print("plots/zmmtozbb_fit%i_%s.png" % (ifunc,h.GetName() if logy else h.GetName()+"_lin"))
        #gPad.Print("plots/zmmtozbb_fit%i_%s.pdf" % (ifunc,h.GetName() if logy else h.GetName()+"_lin"))


    elif ifunc == 1:
        # Fit exponentially modified gaussian
        #ws.importClassCode(RooEMGaussian.Class(), 1)
        
        #ws.factory("Exponential::expo(x,gle[1,%.1f,%.1f])"
        #           % (lamblow,lambup) )
        #ws.factory("EXPR::expo('exp(gle*(x-gm1))', x,gm1[1.0,%.1f,%.1f],gl1[1,%.1f,%.1f])"
        #           % (gmlow,gmup, lamblow,lambup) )
        ws.factory("EMGaussian::emgaussian(x,gm1[1.0,%.1f,%.1f],gs1[0.1,%.1f,%.1f],gl1[5,%.1f,%.1f])"
                   % (gmlow,gmup, gslow,gsup, lamblow,lambup) )
        #ws.factory("SUM::emgaussian_gauss(emgaussian, f[0.5,0,1]*Gaussian::gaus2(x,gm2[1.0,%.1f,%.1f],gs2[0.5,%.1f,%.1f]) )"
        #           % (gmlow,gmup, gslow,gsup) )
                
        #model = ws.pdf("expo")
        model = ws.pdf("emgaussian")
        #model = ws.pdf("emgaussian_gauss")
        
        #result = model.fitTo(data, RooFit.Save())
        result = model.fitTo(data, RooFit.Range(fxlow,fxup), RooFit.Save())
        
        ws.factory("expr::emgm1('gm1-1/gl1',gm1,gl1)")
        ws.factory("expr::emgs1('sqrt(gs1*gs1+1/gl1/gl1)',gs1,gl1)")
        ws.factory("Gaussian::toymodel1(x,emgm1,emgs1)")
        toymodel1 = ws.pdf("toymodel1")
        #ws.factory("EXPR::toymodel2('exp(gl1*gl1 * gs1*gs1 / 2.0 + gl1*(x-gm1))', x,gm1,gs1,gl1)")
        #toymodel2 = ws.pdf("toymodel2")
        
        data.plotOn(frame,RooFit.Name("data"),RooFit.DataError(RooAbsData.SumW2))
        toymodel1.plotOn(frame,RooFit.LineStyle(kDashed),RooFit.LineColor(kBlue-7))
        #toymodel2.plotOn(frame,RooFit.LineStyle(kDashed),RooFit.LineColor(kBlue-7))
        model.plotOn(frame,RooFit.Name("model"))
        frame.SetMaximum(h.GetMaximum() if logy else 6000)
        frame.SetMinimum(0.5 if logy else 0)
        frame.Draw()
        
        chi2 = frame.chiSquare("model", "data", 3)
        p1 = ws.var("gm1")
        p2 = ws.var("gs1")
        p3 = ws.var("gl1")
        result.floatParsFinal().Print("s")
        print "chi^2 = ", chi2
        
        latex.DrawLatex(0.19, 0.90, "%.0f<p_{T}^{Z}<%.0f, %.0f<p_{T}^{j}<%.0f, %.1f<|#eta^{j}|<%.1f" % tuple([hl.GetBinContent(l) for l in xrange(1,7)]))
        latex.DrawLatex(0.19, 0.85, "#color[1]{FullSim}")
        latex.DrawLatex(0.19, 0.80, "#color[4]{Fit EM gaus}")
        latex.DrawLatex(0.42, 0.85, "#color[1]{#mu=%.3f, #sigma=%.3f}" % (h.GetMean(), h.GetRMS()))
        #latex.DrawLatex(0.42, 0.80, "#color[4]{#mu=%.3f #pm %.3f, #sigma=%.3f #pm %.3f}" % (p1.getVal(), p1.getError(), p2.getVal(), p2.getError() ))
        #latex.DrawLatex(0.42, 0.76, "#color[4]{#lambda=%.3f #pm %.3f}" % (p3.getVal(), p3.getError() ))
        latex.DrawLatex(0.42, 0.80, "#color[4]{#mu=%.3f, #sigma=%.3f, #lambda=%.3f}" % (p1.getVal(), p2.getVal(), p3.getVal() ))
        latex.DrawLatex(0.19, 0.70, "#color[2]{#chi_{#nu}^{2}=%.2f}" % (chi2))
        gPad.RedrawAxis()
        gPad.Print("plots/zmmtozbb_fit%i_%s.png" % (ifunc,h.GetName() if logy else h.GetName()+"_lin"))
        #gPad.Print("plots/zmmtozbb_fit%i_%s.pdf" % (ifunc,h.GetName() if logy else h.GetName()+"_lin"))


    elif ifunc == 2:
        # Fit gaussian + expo tail
        #ws.importClassCode(RooGaussExp.Class(), 1)
        #ws.importClassCode(RooExpGaussExp.Class(), 1)
        
        #ws.factory("GaussExp::gaussexp(x,gm1[1.0,%.1f,%.1f],gs1[0.1,%.1f,%.1f],gl1[1,%.1f,%.1f])"
        #           % (gmlow,gmup, gslow,gsup, lamblow,lambup) )
        ws.factory("GaussExp::gaussexp(x,gm1[1.0,%.1f,%.1f],gs1[0.1,%.1f,%.1f],gk1[-1,%.1f,%.1f],Flipped)"
                   % (gmlow,gmup, gslow,gsup, kapplow,kappup) )
        
        model = ws.pdf("gaussexp")
        
        #result = model.fitTo(data, RooFit.Save())
        result = model.fitTo(data, RooFit.Range(fxlow,fxup), RooFit.Save())
        
        ws.factory("Gaussian::toymodel1(x,gm1,gs1)")
        toymodel1 = ws.pdf("toymodel1")
        
        data.plotOn(frame,RooFit.Name("data"),RooFit.DataError(RooAbsData.SumW2))
        model.plotOn(frame,RooFit.Name("model"), RooFit.Invisible())
        toymodel1.plotOn(frame,RooFit.Name("toymodel1"), RooFit.Invisible())
        p1 = ws.var("gm1"); sf1 = frame.findObject("model").Eval(p1.getVal()) / frame.findObject("toymodel1").Eval(p1.getVal())
        toymodel1.plotOn(frame,RooFit.LineStyle(kDashed), RooFit.LineColor(kBlue-7),RooFit.Normalization(sf1,RooAbsReal.Relative))
        model.plotOn(frame)
        frame.SetMaximum(h.GetMaximum() if logy else 6000)
        frame.SetMinimum(0.5 if logy else 0)
        frame.Draw()
        
        chi2 = frame.chiSquare("model", "data", 3)
        p1 = ws.var("gm1")
        p2 = ws.var("gs1")
        #p3 = ws.var("gl1")
        p3 = ws.var("gk1")
        result.floatParsFinal().Print("s")
        print "chi^2 = ", chi2
        
        modelcurve = frame.findObject("model")
        lx1 = p3.getVal() * p2.getVal() + p1.getVal()
        ly1 = modelcurve.Eval(lx1)
        line1 = TLine(lx1, 0, lx1, ly1)
        line1.SetLineWidth(3); line1.SetLineStyle(kDashed); line1.SetLineColor(kBlue-7)
        line1.Draw()
        print lx1, ly1
        
        latex.DrawLatex(0.19, 0.90, "%.0f<p_{T}^{Z}<%.0f, %.0f<p_{T}^{j}<%.0f, %.1f<|#eta^{j}|<%.1f" % tuple([hl.GetBinContent(l) for l in xrange(1,7)]))
        latex.DrawLatex(0.19, 0.85, "#color[1]{FullSim}")
        latex.DrawLatex(0.19, 0.80, "#color[4]{Fit gaus+exp}")
        latex.DrawLatex(0.42, 0.85, "#color[1]{#mu=%.3f, #sigma=%.3f}" % (h.GetMean(), h.GetRMS()))
        #latex.DrawLatex(0.42, 0.80, "#color[4]{#mu=%.3f #pm %.3f, #sigma=%.3f #pm %.3f}" % (p1.getVal(), p1.getError(), p2.getVal(), p2.getError() ))
        #latex.DrawLatex(0.42, 0.76, "#color[4]{#lambda=%.3f #pm %.3f}" % (p3.getVal(), p3.getError() ))
        #latex.DrawLatex(0.42, 0.80, "#color[4]{#mu=%.3f, #sigma=%.3f, #lambda=%.3f}" % (p1.getVal(), p2.getVal(), p3.getVal() ))
        latex.DrawLatex(0.42, 0.80, "#color[4]{#mu=%.3f, #sigma=%.3f, #kappa=%.3f}" % (p1.getVal(), p2.getVal(), p3.getVal() ))
        latex.DrawLatex(0.19, 0.70, "#color[2]{#chi_{#nu}^{2}=%.2f}" % (chi2))
        gPad.RedrawAxis()
        gPad.Print("plots/zmmtozbb_fit%i_%s.png" % (ifunc,h.GetName() if logy else h.GetName()+"_lin"))
        #gPad.Print("plots/zmmtozbb_fit%i_%s.pdf" % (ifunc,h.GetName() if logy else h.GetName()+"_lin"))


    elif ifunc == 3:
        # Fit gaussian + expo tail x 2
        #ws.importClassCode(RooGaussExp.Class(), 1)
        #ws.importClassCode(RooExpGaussExp.Class(), 1)
        
        #ws.factory("GaussExp::gaussexp(x,gm1[1.0,%.1f,%.1f],gs1[0.1,%.1f,%.1f],gl1[1,%.1f,%.1f])"
        #           % (gmlow,gmup, gslow,gsup, lamblow,lambup) )
        ws.factory("ExpGaussExp::expgaussexp(x,gm1[1.0,%.1f,%.1f],gs1[0.1,%.1f,%.1f],gl1[2,%.1f,%.1f],gk1[-1,%.1f,%.1f])"
                   % (gmlow,gmup, gslow,gsup, lamblow,lambup, kapplow,kappup) )
        
        #model = ws.pdf("gaussexp")
        model = ws.pdf("expgaussexp")
        
        #result = model.fitTo(data, RooFit.Save())
        result = model.fitTo(data, RooFit.Range(fxlow,fxup), RooFit.Save())
        
        ws.factory("Gaussian::toymodel1(x,gm1,gs1)")
        toymodel1 = ws.pdf("toymodel1")
        
        data.plotOn(frame,RooFit.Name("data"),RooFit.DataError(RooAbsData.SumW2))
        model.plotOn(frame,RooFit.Name("model"), RooFit.Invisible())
        toymodel1.plotOn(frame,RooFit.Name("toymodel1"), RooFit.Invisible())
        p1 = ws.var("gm1"); sf1 = frame.findObject("model").Eval(p1.getVal()) / frame.findObject("toymodel1").Eval(p1.getVal())
        toymodel1.plotOn(frame,RooFit.LineStyle(kDashed), RooFit.LineColor(kBlue-7),RooFit.Normalization(sf1,RooAbsReal.Relative))
        model.plotOn(frame)
        frame.SetMaximum(h.GetMaximum() if logy else 6000)
        frame.SetMinimum(0.5 if logy else 0)
        frame.Draw()
        
        #chi2 = frame.chiSquare("model", "data", 3)
        chi2 = frame.chiSquare("model", "data", 4)
        p1 = ws.var("gm1")
        p2 = ws.var("gs1")
        p3 = ws.var("gk1")
        p4 = ws.var("gl1")
        result.floatParsFinal().Print("s")
        print "chi^2 = ", chi2
        
        modelcurve = frame.findObject("model")
        lx1 = p3.getVal() * p2.getVal() + p1.getVal()
        ly1 = modelcurve.Eval(lx1)
        line1 = TLine(lx1, 0, lx1, ly1)
        line1.SetLineWidth(3); line1.SetLineStyle(kDashed); line1.SetLineColor(kBlue-7)
        line1.Draw()
        print lx1, ly1
        lx2 = p4.getVal() * p2.getVal() + p1.getVal()
        ly2 = modelcurve.Eval(lx2)
        line2 = TLine(lx2, 0, lx2, ly2)
        line2.SetLineWidth(3); line2.SetLineStyle(kDashed); line2.SetLineColor(kBlue-7)
        line2.Draw()
        print lx2, ly2
        
        latex.DrawLatex(0.19, 0.90, "%.0f<p_{T}^{Z}<%.0f, %.0f<p_{T}^{j}<%.0f, %.1f<|#eta^{j}|<%.1f" % tuple([hl.GetBinContent(l) for l in xrange(1,7)]))
        latex.DrawLatex(0.19, 0.85, "#color[1]{FullSim}")
        latex.DrawLatex(0.19, 0.80, "#color[4]{Fit gaus+expx2}")
        latex.DrawLatex(0.42, 0.85, "#color[1]{#mu=%.3f, #sigma=%.3f}" % (h.GetMean(), h.GetRMS()))
        #latex.DrawLatex(0.42, 0.80, "#color[4]{#mu=%.3f #pm %.3f, #sigma=%.3f #pm %.3f}" % (p1.getVal(), p1.getError(), p2.getVal(), p2.getError() ))
        #latex.DrawLatex(0.42, 0.76, "#color[4]{#lambda=%.3f #pm %.3f}" % (p3.getVal(), p3.getError() ))
        #latex.DrawLatex(0.42, 0.80, "#color[4]{#mu=%.3f, #sigma=%.3f, #lambda=%.3f}" % (p1.getVal(), p2.getVal(), p3.getVal() ))
        latex.DrawLatex(0.42, 0.80, "#color[4]{#mu=%.3f, #sigma=%.3f, #kappa=%.3f, #lambda=%.3f}" % (p1.getVal(), p2.getVal(), p3.getVal(), p4.getVal() ))
        latex.DrawLatex(0.19, 0.70, "#color[2]{#chi_{#nu}^{2}=%.2f}" % (chi2))
        gPad.RedrawAxis()
        gPad.Print("plots/zmmtozbb_fit%i_%s.png" % (ifunc,h.GetName() if logy else h.GetName()+"_lin"))
        #gPad.Print("plots/zmmtozbb_fit%i_%s.pdf" % (ifunc,h.GetName() if logy else h.GetName()+"_lin"))


    if ifunc == 4:
        # Fit double gaussian
        ws.factory("BifurGauss::bifurgauss(x,gm1[1.0,%.1f,%.1f],gs1[0.1,%.1f,%.1f],gs2[0.2,%.1f,%.1f])" 
                    % (0.9,gmup, gslow,gsup, gslow,gsup))
        model = ws.pdf("bifurgauss")
        
        #result = model.fitTo(data, RooFit.Save())
        result = model.fitTo(data, RooFit.Range(fxlow,fxup), RooFit.Save())
        
        data.plotOn(frame,RooFit.Name("data"),RooFit.DataError(RooAbsData.SumW2))
        #model.plotOn(frame,RooFit.Components("gaus1"),RooFit.LineStyle(kDashed),RooFit.LineColor(kBlue-7))
        #model.plotOn(frame,RooFit.Components("gaus2"),RooFit.LineStyle(kDashed),RooFit.LineColor(kBlue-7))
        model.plotOn(frame,RooFit.Name("model"))
        frame.SetMaximum(h.GetMaximum() if logy else 6000)
        frame.SetMinimum(0.5 if logy else 0)
        frame.Draw()
        
        chi2 = frame.chiSquare("model", "data", 3)
        p1 = ws.var("gm1")
        p2 = ws.var("gs1")
        p3 = ws.var("gs2")
        result.floatParsFinal().Print("s")
        print "chi^2 = ", chi2
        
        latex.DrawLatex(0.19, 0.90, "%.0f<p_{T}^{Z}<%.0f, %.0f<p_{T}^{j}<%.0f, %.1f<|#eta^{j}|<%.1f" % tuple([hl.GetBinContent(l) for l in xrange(1,7)]))
        latex.DrawLatex(0.19, 0.85, "#color[1]{FullSim}")
        latex.DrawLatex(0.19, 0.80, "#color[4]{Fit bifur gaus}")
        latex.DrawLatex(0.42, 0.85, "#color[1]{#mu=%.3f, #sigma=%.3f}" % (h.GetMean(), h.GetRMS()))
        latex.DrawLatex(0.42, 0.80, "#color[4]{#mu_{1}=%.3f, #sigma_{1}=%.3f, #sigma_{2}=%.3f}" % (p1.getVal(), p2.getVal(), p3.getVal() ))
        latex.DrawLatex(0.19, 0.70, "#color[2]{#chi_{#nu}^{2}=%.2f}" % (chi2))
        gPad.RedrawAxis()
        gPad.Print("plots/zmmtozbb_fit%i_%s.png" % (ifunc,h.GetName() if logy else h.GetName()+"_lin"))
        #gPad.Print("plots/zmmtozbb_fit%i_%s.pdf" % (ifunc,h.GetName() if logy else h.GetName()+"_lin"))


# Print workspace contents
ws.Print()

