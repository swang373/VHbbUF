#!/usr/bin/env python

from ROOT import gROOT, gStyle,  TStyle


def tdrStyle():
#    gROOT.SetStyle("Plain")
    tdrStyle = TStyle("tdrStyle","Style for P-TDR")
    tdrStyle.SetPalette(1)
    
    tdrStyle.SetAxisColor(1, "XYZ")
    
    tdrStyle.SetCanvasColor(0)
    #tdrStyle.SetCanvasBorderSize(10)
    tdrStyle.SetCanvasBorderMode(0)
    tdrStyle.SetCanvasDefH(700)
    tdrStyle.SetCanvasDefW(700)
    tdrStyle.SetCanvasDefX(0)
    tdrStyle.SetCanvasDefY(0)
    
    tdrStyle.SetFitFormat("5.4g")
    tdrStyle.SetFuncColor(2)
    tdrStyle.SetFuncStyle(1)
    tdrStyle.SetFuncWidth(1)

    tdrStyle.SetFrameBorderMode(0)
    tdrStyle.SetFrameBorderSize(1)
    tdrStyle.SetFrameFillStyle(0)
    tdrStyle.SetFrameFillColor(0)
    tdrStyle.SetFrameLineColor(1)
    tdrStyle.SetFrameLineStyle(1)  # 0?
    tdrStyle.SetFrameLineWidth(1)  # 1?

    tdrStyle.SetGridColor(0)
    tdrStyle.SetGridStyle(3)
    tdrStyle.SetGridWidth(1)

    #tdrStyle.SetHistFillColor(1)
    #tdrStyle.SetHistFillStyle(0)
    tdrStyle.SetHistLineColor(1)
    tdrStyle.SetHistLineStyle(0)
    tdrStyle.SetHistLineWidth(1)
    
    tdrStyle.SetLabelColor(1, "XYZ")
    tdrStyle.SetLabelFont(42,"XYZ")
    tdrStyle.SetLabelOffset(0.007,"XYZ")  # 0.010?
    tdrStyle.SetLabelSize(0.05,"XYZ")  # 0.04?
    
    tdrStyle.SetLegendBorderSize(0)
    tdrStyle.SetLegendFillColor(0)
    tdrStyle.SetLegendFont(42)
    
    tdrStyle.SetMarkerSize(1.0)
    tdrStyle.SetMarkerStyle(20)
    
    tdrStyle.SetLineColor(1)
    tdrStyle.SetLineWidth(2)
    #tdrStyle.SetLineScalePS(2)

    tdrStyle.SetOptDate(0)
    tdrStyle.SetOptFile(0)
    tdrStyle.SetOptFit(1)
    tdrStyle.SetOptStat(0)
    tdrStyle.SetOptTitle(0)
    #tdrStyle.SetOptLogx(0)
    #tdrStyle.SetOptLogy(0)
    #tdrStyle.SetOptLogz(0)

    tdrStyle.SetPadColor(0)
    tdrStyle.SetPadBorderMode(0)
    tdrStyle.SetPadBorderSize(10)
    tdrStyle.SetPadTopMargin(0.05)  # 0.08?
    tdrStyle.SetPadBottomMargin(0.13)
    tdrStyle.SetPadLeftMargin(0.16)
    tdrStyle.SetPadRightMargin(0.03)  # 0.05?
    tdrStyle.SetPadGridX(0)
    tdrStyle.SetPadGridY(0)
    tdrStyle.SetPadTickX(1)
    tdrStyle.SetPadTickY(1)

    tdrStyle.SetStatColor(0)
    tdrStyle.SetStatFont(42)
    tdrStyle.SetStatFontSize(0.025)
    tdrStyle.SetStatTextColor(1)
    tdrStyle.SetStatFormat("6.4g")
    tdrStyle.SetStatBorderSize(1)
    tdrStyle.SetStatH(0.1)
    tdrStyle.SetStatW(0.15)
    #tdrStyle.SetStatX(0)
    #tdrStyle.SetStatY(0)

    #tdrStyle.SetTextSize(0.055)
    tdrStyle.SetTextFont(42)

    tdrStyle.SetTitleBorderSize(0)
    tdrStyle.SetTitleColor(1)
    tdrStyle.SetTitleFont(42)
    tdrStyle.SetTitleColor(1,"XYZ")
    tdrStyle.SetTitleFont(42,"XYZ")  
    tdrStyle.SetTitleSize(0.06,"XYZ")  # 0.05?
    #tdrStyle.SetTitleOffset(1.4,"XYZ")
    tdrStyle.SetTitleOffset(0.9,"X")
    tdrStyle.SetTitleOffset(1.20,"Y")
    tdrStyle.SetTitleFillColor(10)
    tdrStyle.SetTitleFontSize(0.05)
    tdrStyle.SetTitleTextColor(1)
    #tdrStyle.SetTitleH(0)
    #tdrStyle.SetTitleW(0)
    #tdrStyle.SetTitleX(0)
    #tdrStyle.SetTitleY(0.985)
    #tdrStyle.SetTitleStyle(1001)

    tdrStyle.SetPalette(1)
    #tdrStyle.SetNdivisions(510, "XYZ")  # 505?
    tdrStyle.SetNdivisions(505, "XYZ")
    tdrStyle.SetEndErrorSize(0)  # 2?
    #tdrStyle.SetErrorMarker(20)
    #tdrStyle.SetErrorX(0.)
    #tdrStyle.SetPaperSize(20.,20.)
    tdrStyle.SetStripDecimals(1)
    tdrStyle.SetTickLength(0.03, "XYZ")
    return 1

def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i+lv/3], 16) for i in range(0, lv, lv/3))

def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb
