import ROOT


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
