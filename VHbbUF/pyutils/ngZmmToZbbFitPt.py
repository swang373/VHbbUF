from ROOT import TH1, TH1F, TCanvas, gROOT, gPad

gROOT.LoadMacro("tdrstyle.C")
gROOT.ProcessLine("setTDRStyle()")

genzptbins = [
    (0,1000),
]
genptbins = [
    (20,30),
    (30,50),
    (50,70),
    (70,100),
    (100,160),
    (160,320),
    (320,640),
]
#genptbins = [
#    (30,60),
#    (60,100),
#    (100,200),
#    (200,2000),
#]
etabins = [
    (0.0,0.5),
    (0.5,1.0),
    (1.0,1.5),
    (1.5,2.5),
#NOSTAT    (2.5,5.0),
]

parvalues_ = [
0.904 ,  0.155 ,  1.109 ,  1.407 ,  0.270 , 
0.950 ,  0.115 ,  0.924 ,  1.877 ,  0.758 , 
0.982 ,  0.110 ,  0.908 ,  1.906 ,  0.579 , 
0.986 ,  0.100 ,  0.895 ,  1.924 ,  0.499 , 
0.986 ,  0.090 ,  0.875 ,  1.985 ,  0.444 , 
0.995 ,  0.080 ,  0.876 ,  2.290 ,  0.389 , 
1.004 ,  0.068 ,  0.963 ,  3.198 ,  0.468 , 
0.904 ,  0.161 ,  1.123 ,  1.396 ,  0.250 , 
0.971 ,  0.093 ,  0.922 ,  2.210 ,  0.857 , 
0.981 ,  0.110 ,  0.906 ,  1.867 ,  0.603 , 
0.984 ,  0.101 ,  0.890 ,  1.892 ,  0.488 , 
0.984 ,  0.093 ,  0.871 ,  1.922 ,  0.425 , 
0.991 ,  0.080 ,  0.887 ,  2.200 ,  0.415 , 
1.004 ,  0.062 ,  0.962 ,  3.353 ,  0.542 , 
0.893 ,  0.173 ,  1.093 ,  1.292 ,  0.332 , 
0.968 ,  0.090 ,  0.915 ,  2.392 ,  0.949 , 
0.984 ,  0.112 ,  0.915 ,  1.817 ,  0.771 , 
0.985 ,  0.107 ,  0.899 ,  1.801 ,  0.617 , 
0.985 ,  0.097 ,  0.878 ,  1.834 ,  0.521 , 
0.990 ,  0.089 ,  0.882 ,  2.045 ,  0.445 , 
1.010 ,  0.074 ,  0.960 ,  2.655 ,  0.543 , 
0.897 ,  0.185 ,  1.099 ,  1.283 ,  0.323 , 
1.001 ,  0.076 ,  0.914 ,  2.848 ,  0.970 , 
0.993 ,  0.103 ,  0.900 ,  1.955 ,  0.781 , 
0.985 ,  0.101 ,  0.879 ,  1.853 ,  0.614 , 
0.982 ,  0.087 ,  0.860 ,  1.980 ,  0.534 , 
0.987 ,  0.075 ,  0.858 ,  2.333 ,  0.459 , 
0.988 ,  0.074 ,  0.914 ,  2.973 ,  0.402 ,
]

parerrors_ = [
0.007 ,  0.005 ,  0.057 ,  0.143 ,  0.094 , 
0.005 ,  0.010 ,  0.002 ,  0.126 ,  0.054 , 
0.002 ,  0.003 ,  0.003 ,  0.033 ,  0.023 , 
0.001 ,  0.001 ,  0.002 ,  0.021 ,  0.013 , 
0.001 ,  0.001 ,  0.002 ,  0.019 ,  0.010 , 
0.001 ,  0.001 ,  0.003 ,  0.029 ,  0.010 , 
0.003 ,  0.004 ,  0.008 ,  0.156 ,  0.035 , 
0.008 ,  0.006 ,  0.064 ,  0.168 ,  0.100 , 
0.005 ,  0.007 ,  0.002 ,  0.148 ,  0.022 , 
0.002 ,  0.003 ,  0.003 ,  0.036 ,  0.023 , 
0.001 ,  0.001 ,  0.002 ,  0.021 ,  0.014 , 
0.001 ,  0.001 ,  0.002 ,  0.019 ,  0.010 , 
0.001 ,  0.001 ,  0.003 ,  0.030 ,  0.012 , 
0.003 ,  0.004 ,  0.008 ,  0.186 ,  0.034 , 
0.016 ,  0.010 ,  0.074 ,  0.158 ,  0.168 , 
0.014 ,  0.017 ,  0.002 ,  0.176 ,  0.010 , 
0.004 ,  0.005 ,  0.002 ,  0.064 ,  0.025 , 
0.002 ,  0.002 ,  0.003 ,  0.032 ,  0.019 , 
0.001 ,  0.001 ,  0.003 ,  0.025 ,  0.014 , 
0.002 ,  0.002 ,  0.004 ,  0.035 ,  0.017 , 
0.006 ,  0.009 ,  0.012 ,  0.262 ,  0.073 , 
0.017 ,  0.010 ,  0.083 ,  0.170 ,  0.187 , 
0.019 ,  0.013 ,  0.002 ,  0.461 ,  0.010 , 
0.004 ,  0.004 ,  0.002 ,  0.064 ,  0.017 , 
0.002 ,  0.002 ,  0.002 ,  0.028 ,  0.014 , 
0.001 ,  0.001 ,  0.002 ,  0.024 ,  0.010 , 
0.001 ,  0.001 ,  0.004 ,  0.043 ,  0.014 , 
0.007 ,  0.008 ,  0.025 ,  0.322 ,  0.076 ,
]

parvalues = []
npars = 5

ih = 0
for etabin in etabins:
    for genptbin in genptbins:
        for genzptbin in genzptbins:
            for i in xrange(npars):
                parvalues.append((parvalues_[ih], parerrors_[ih]))                
                ih += 1

#print parvalues

histos = []
for i in xrange(npars):
    histos_ = []
    for etabin in etabins:
        h1 = TH1F("h1_par%i_eta%i" % (i,etabins.index(etabin)), ";parton p_{T} [GeV] {%.1f#leq|#eta|<%.1f}; p%i" % (etabin[0],etabin[1],i), 200, 0, 600)
        histos_.append(h1)
    histos.append(histos_)

ip = 0
for y, ey in parvalues:
    pbin = ip % npars
    etabin = ip / (npars*len(genptbins))
    ptbin = (ip / npars) % len(genptbins)
    h = histos[pbin][etabin]
    x = 0.5*(genptbins[ptbin][0]+genptbins[ptbin][1])
    b = h.FindFixBin(x)
    h.SetBinContent(b, y)
    h.SetBinError(b, ey)
    
    if pbin == 3 and etabin == 1: print ip, x, y, ey, etabin, ptbin
    
    ip += 1


def setRatioStyle(ratio, ymin, ymax):
    #ratio.Sumw2()
    ratio.SetStats(0)
    ratio.SetTitle("")
    #ratio.GetXaxis().SetTitle(xtitle)
    #ratio.GetYaxis().SetTitle("Data/MC")
    #ratio.GetYaxis().SetTitle("Z(bb)H(inv)/Z(#nu#nu)H(bb)")
    #ratio.GetYaxis().SetTitle("Z(bb)H(inv)/t#bar{t}")
    ratio.SetMaximum(ymax)
    ratio.SetMinimum(ymin)
    #ratio.SetMarkerSize(0)
    #ratio.SetFillColor(923)  # kGray+3
    #ratio.SetFillStyle(3013)
    ratio.GetXaxis().CenterTitle()
    ratio.GetXaxis().SetLabelSize(0.10)
    ratio.GetXaxis().SetTitleSize(0.12)
    ratio.GetXaxis().SetTitleOffset(1.10)
    #ratio.GetXaxis().SetNdivisions(505)
    ratio.GetYaxis().CenterTitle()
    ratio.GetYaxis().SetLabelSize(0.10)
    #ratio.GetYaxis().SetTitleSize(0.12)
    ratio.GetYaxis().SetTitleSize(0.12)
    ratio.GetYaxis().SetTitleOffset(0.6)
    ratio.GetYaxis().SetNdivisions(505)



c1 = TCanvas("c1", "c1", 600, 210*len(etabins))
c1.Divide(1,len(etabins))


# p0
ymin,ymax = 0.8,1.1
for i in xrange(len(etabins)):
    c1.cd(i+1)
    gPad.SetBottomMargin(0.35)
    h = histos[0][i]
    setRatioStyle(h,ymin,ymax)
    h.Draw("e")
c1.Print("plots/zmmtozbb_fit0_p0.png")

# p1
ymin,ymax = 0.04,0.14
for i in xrange(len(etabins)):
    c1.cd(i+1)
    gPad.SetBottomMargin(0.35)
    h = histos[1][i]
    setRatioStyle(h,ymin,ymax)
    h.Draw("e")
c1.Print("plots/zmmtozbb_fit0_p1.png")

# p2
ymin,ymax = 0.7,1.2
for i in xrange(len(etabins)):
    c1.cd(i+1)
    gPad.SetBottomMargin(0.35)
    h = histos[2][i]
    setRatioStyle(h,ymin,ymax)
    h.Draw("e")
c1.Print("plots/zmmtozbb_fit0_p2.png")

# p3
ymin,ymax = 0.5,4
for i in xrange(len(etabins)):
    c1.cd(i+1)
    gPad.SetBottomMargin(0.35)
    h = histos[3][i]
    setRatioStyle(h,ymin,ymax)
    h.Draw("e")
c1.Print("plots/zmmtozbb_fit0_p3.png")

# p4
ymin,ymax = 0.,1.1
for i in xrange(len(etabins)):
    c1.cd(i+1)
    gPad.SetBottomMargin(0.35)
    h = histos[4][i]
    setRatioStyle(h,ymin,ymax)
    h.Draw("e")
c1.Print("plots/zmmtozbb_fit0_p4.png")
