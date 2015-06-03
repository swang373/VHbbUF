#!usr/bin/python
import ROOT
from ROOT import TFile, TH1, TCanvas, TProfile, TEfficiency, TGraphAsymmErrors, TF1, TPaveText, TPaveStats, TLegend, TLine, gROOT, gPad, gStyle
from setTDRStyle import setTDRStyle
from array import array

gROOT.Reset()
setTDRStyle()
gROOT.SetStyle('tdrStyle')
gROOT.ForceStyle()

#import optparse
#usage = "usage: %prog [options]"
#parser = optparse.OptionParser(usage)
#parser.add_option("--bit"  ,action="store",type="string",dest="bit",default='pfCor40')
#parser.add_option("--var"  ,action="store",type="string",dest="var",default='Inclusive')

#(options, args) = parser.parse_args()

#bit = options.bit
#var = options.var

COLOR = [ROOT.kRed,ROOT.kBlue]
STYLE = [21,22]

TAGS = ["baseline", "regression"  ]

canPt = TCanvas('ResVsPt','ResVsPt_',900,600)

resPt= []
gEffPt = []
fit = []
ln = []
k = 0



for tag in TAGS:
  inf      = TFile.Open('/afs/cern.ch/work/d/degrutto/public/MiniAOD/ZnnHbb_Phys14_PU20bx25/skimV11/step3/Step3_ZnnH125.root')
  tree     = inf.Get('tree')
  hRes   = TProfile('hRes_'+tag,'hRes_'+tag,8,0,400, -100, 100, "s")
  if tag == 'baseline':
#    tree.Draw('genPt>>hPtAll_'+tag, 'genPt>5')
#    tree.Draw('genPt>>hPtPass_'+tag,'pass && genMatchDR<0.2')
    tree.Draw('(Jet_pt[hJCidx] / Jet_mcPt[hJCidx]):  Jet_mcPt[hJCidx] >>hRes_'+tag, "",  "profs" )
  else:
    tree.Draw('(hJet_ptReg / Jet_mcPt[hJCidx]):  Jet_mcPt[hJCidx] >>hRes_'+tag, "",  "profs" )
  resPt.append(hRes)
  resPt[k].SetMarkerColor(COLOR[k])
  resPt[k].SetLineColor(COLOR[k])
  resPt[k].SetMarkerStyle(STYLE[k]) 

  x   = []
  y   = []
  exl = []
  exh = []
  eyl = []
  eyh = []

  for i in xrange(hRes.GetNbinsX()):
    print tag
    print hRes.GetBinCenter(i+1), hRes.GetBinContent(i+1), hRes.GetBinError(i+1)
    x.append(hRes.GetBinCenter(i+1))
    exl.append(0) 
    exh.append(0) 
    if hRes.GetBinContent(i+1)>0 :
      # resolution as sigma(hlt/gen) / <hlt/gen>
      y.append( hRes.GetBinError(i+1) / hRes.GetBinContent(i+1))
    else:
      y.append(0)
    eyl.append(0)
    eyh.append(0)

  vx = array('d',x)
  vy = array('d',y)
  vexl = array('d',exl)
  vexh = array('d',exh)
  veyl = array('d',eyl)
  veyh = array('d',eyh)

  gEffPt.append(TGraphAsymmErrors(hRes.GetNbinsX(),vx,vy,vexl,vexh,veyl,veyh))

  gEffPt[k].SetMarkerColor(COLOR[k])
  gEffPt[k].SetLineColor(COLOR[k])
  gEffPt[k].SetMarkerStyle(STYLE[k]) 

 # fit.append(TF1('fit'+str(k),'([0]-[3])/pow(1+[1]*exp(-[2]*x),[4])+[3]+[5]*(1-pow(x,-[6]))',15,200))
#  fit.append(TF1('fit'+str(k),'([0]-[3])/pow(1+[1]*exp(-[2]*x),[4])+[3]+[5]*(1-pow(x,-[6]))',0,300))
#  fit.append(TF1('fit'+str(k),'gauss',0,300))
#  fit[k].SetParameters(1,300,0.1,0,1,1,0.1)
#  fit[k].SetParNames('A','x0','#sigma')
#  fit[k].SetLineColor(COLOR[k])

  k+=1

canPt.cd()
#h = TH1F('h','h',50,0,100)
f = TF1('f','1',0,300)
f.SetLineColor(ROOT.kBlack)
f.SetLineStyle(3)
f.SetLineWidth(1)
f.GetXaxis().SetTitle('GenPt (GeV)')
f.GetYaxis().SetTitle('#sigma(pt / gen) / <pt/gen>')
#h.GetXaxis().SetRangeUser(20,70)
f.SetMinimum(0)
f.SetMaximum(2)
f.Draw()

leg = TLegend(0.7,0.2,0.9,0.45)
#leg.SetHeader(bit)
#gEffPt[0].Draw('PE1')

for k in range(len(TAGS)):
#  gEffPt[k].Fit(fit[k],'RQ')
#  x99 = fit[k].GetX(0.99)
#  x10 = fit[k].GetX(0.1)
#  x50 = fit[k].GetX(0.5)
#  x90 = fit[k].GetX(0.9)
##  print 'Turn-on: '+str(x90)
##  q = (x90-x10)/x50
##  print 'Sharpness: '+str(q)
#  gEffPt[k].Draw('samePE1)'
#  gEffPt[k].Draw('samePE')
  resPt[k].Draw('samePE')
#  #fit[k].Draw('sameL')
  gPad.Update()

#  ln.append(TLine(x90,gPad.GetFrame().GetY1(),x90,gPad.GetFrame().GetY2()))
#  ln[k].SetLineColor(COLOR[k])
#  ln[k].SetLineStyle(7)
#  ln[k].SetLineWidth(2)
#  ln[k].Draw("same")
 ## leg.AddEntry(gEffPt[k],TAGS[k]+' ('+'%.2f' % q +')','P')
  leg.AddEntry(gEffPt[k],TAGS[k],'P')

leg.SetFillColor(0)
#leg.SetBorderSize(0)
leg.SetTextFont(62)
leg.SetTextSize(0.05)
leg.Draw()

gPad.Update()

#----- keep the GUI alive ------------
if __name__ == '__main__':
  rep = ''
  while not rep in ['q','Q']:
    rep = raw_input('enter "q" to quit: ')
    if 1 < len(rep):
      rep = rep[0]
