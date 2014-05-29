import array
from math import sqrt, cos
from ROOT import gROOT, AddressOf

gROOT.ProcessLine("typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVectorM;")
gROOT.ProcessLine("typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double> > LorentzVectorE;")

gROOT.LoadMacro("HelperFunctions.h")
gROOT.LoadMacro("HelperVHbbDataFormats.h")
gROOT.LoadMacro("$HOME/style-CMSTDR.C");

from ROOT import TFile, TChain, TTree, TTreeFormula, TH1F, TLatex, TLegend, TColor, gPad, gROOT, HiggsInfo, METInfo, MHTInfo, genParticleInfo, LorentzVectorM, LorentzVectorE, setTDRStyle, deltaR, deltaPhi, triggerweight2012ABC
setTDRStyle()
gROOT.SetBatch(1)

def alphaT(pt0, phi0, pt1, phi1):
    return pt1/sqrt(2 * pt0 * pt1 * (1.0 - cos(deltaPhi(phi0, phi1))) )


dirStep3 = "Step3_20130204/"
chQCD = TChain("tree", "QCD")
chQCD.Add(dirStep3+"Step3_QCDPt80.root")
chQCD.Add(dirStep3+"Step3_QCDPt120.root")
chQCD.Add(dirStep3+"Step3_QCDPt170.root")
chQCD.Add(dirStep3+"Step3_QCDPt300.root")
chQCD.Add(dirStep3+"Step3_QCDPt470.root")
chQCD.Add(dirStep3+"Step3_QCDPt600.root")
chQCD.Add(dirStep3+"Step3_QCDPt800.root")
chQCD.Add(dirStep3+"Step3_QCDPt1000.root")
chQCD.Add(dirStep3+"Step3_QCDPt1400.root")
chQCD.Add(dirStep3+"Step3_QCDPt1800.root")

chZH = TChain("tree", "ZH125")
chZH.Add(dirStep3+"Step3_ZnnH125.root")

channel = "ZnunuHighPt"
#channel = "ZnunuLowPt"


strcutHigh = "(Vtype==4)                        && METtype1corr.et> 170 && HptReg>130 && hJet_ptReg[0]>60 && hJet_ptReg[1]>30 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.679 && min(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.244 && HmassReg<250 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && nalep_Znn==0 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && HmassReg>0"
strcutMed  = "(Vtype==4) && 130<METtype1corr.et && METtype1corr.et<=170 && HptReg>130 && hJet_ptReg[0]>60 && hJet_ptReg[1]>30 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.679 && min(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.244 && HmassReg<250 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && nalep_Znn==0 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && HmassReg>0"
strcutHigh1 = "(Vtype==4)                        && METtype1corr.et> 170 && HptReg>130 && hJet_ptReg[0]>60 && hJet_ptReg[1]>30 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.679 && min(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.244 && HmassReg<250 && mindPhiMETJet_dPhi>0.5 && abs(deltaPhi(METtype1corr.phi,METnoPUCh.phi))<0.5 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && nalep_Znn==0 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && HmassReg>0"
strcutMed1  = "(Vtype==4) && 130<METtype1corr.et && METtype1corr.et<=170 && HptReg>130 && hJet_ptReg[0]>60 && hJet_ptReg[1]>30 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.679 && min(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.244 && HmassReg<250 && mindPhiMETJet_dPhi>0.5 && abs(deltaPhi(METtype1corr.phi,METnoPUCh.phi))<0.5 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && nalep_Znn==0 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && HmassReg>0"

strcutWeight = "efflumi * PUweight * triggerweight2012ABC(METtype1corr.et)"

histograms = []

for chain in [chZH, chQCD]:
    METtype1corr = METInfo()
    chain.SetBranchAddress("METtype1corr", AddressOf(METtype1corr, "et"))
    
    METnoPUCh = METInfo()
    chain.SetBranchAddress("METnoPUCh", AddressOf(METnoPUCh, "et"))
    
    MHT = MHTInfo()
    chain.SetBranchAddress("MHT", AddressOf(MHT, "mht"))

    cut = TTreeFormula("cutHigh", strcutHigh, chain) if channel == "ZnunuHighPt" else TTreeFormula("cutMed", strcutMed, chain)
    #cut = TTreeFormula("cutHigh", strcutHigh1, chain) if channel == "ZnunuHighPt" else TTreeFormula("cutMed", strcutMed1, chain)
    cutWeight = TTreeFormula("cutWeight", strcutWeight, chain)

    nevt = chain.GetEntries()
    curTree = chain.GetTreeNumber()
    
    nbins = 1600
    h_mindphi_c50 = TH1F("%s_mindphi_c50" % (chain.GetTitle()), ";min #Delta#varphi(pfMET, cj50)", nbins, 0., 3.2)
    h_mindphi_c30 = TH1F("%s_mindphi_c30" % (chain.GetTitle()), ";min #Delta#varphi(pfMET, cj30)", nbins, 0., 3.2)
    h_mindphi_c25 = TH1F("%s_mindphi_c25" % (chain.GetTitle()), ";min #Delta#varphi(pfMET, cj25)", nbins, 0., 3.2)
    h_mindphi_c20 = TH1F("%s_mindphi_c20" % (chain.GetTitle()), ";min #Delta#varphi(pfMET, cj20)", nbins, 0., 3.2)
    h_mindphi_cf50 = TH1F("%s_mindphi_cf50" % (chain.GetTitle()), ";min #Delta#varphi(pfMET, j50)", nbins, 0., 3.2)
    h_mindphi_cf30 = TH1F("%s_mindphi_cf30" % (chain.GetTitle()), ";min #Delta#varphi(pfMET, j30)", nbins, 0., 3.2)
    h_mindphi_cf25 = TH1F("%s_mindphi_cf25" % (chain.GetTitle()), ";min #Delta#varphi(pfMET, j25)", nbins, 0., 3.2)
    h_mindphi_cf20 = TH1F("%s_mindphi_cf20" % (chain.GetTitle()), ";min #Delta#varphi(pfMET, j20)", nbins, 0., 3.2)
    h_mindphi_c20f30 = TH1F("%s_mindphi_c20f30" % (chain.GetTitle()), ";min #Delta#varphi(pfMET, cj20|fj30)", nbins, 0., 3.2)
    h_mindphi_j3 = TH1F("%s_mindphi_j3" % (chain.GetTitle()), ";min #Delta#varphi(pfMET, 3 lead j)", nbins, 0., 3.2)
    h_mindphi = TH1F("%s_mindphi" % (chain.GetTitle()), ";min #Delta#varphi(pfMET, cj30)", nbins, 0., 3.2)
    h_dphiTrk = TH1F("%s_dphiTrk_lt" % (chain.GetTitle()), ";#Delta #varphi(pfMET, trkMET)", nbins, 0., 3.2)
    h_dphiMht = TH1F("%s_dphiMht_lt" % (chain.GetTitle()), ";#Delta #varphi(pfMET, MHT)", nbins, 0., 3.2)
    h_alphaT = TH1F("%s_alphaT" % (chain.GetTitle()), ";#alpha_{T}(j1,j2)", nbins, 0., 0.8)
    h_metsig = TH1F("%s_metsig" % (chain.GetTitle()), ";pfMET/sqrt(sumET)", nbins, 0., 32.)
    h_metfancysig = TH1F("%s_metfancysig" % (chain.GetTitle()), ";fancy pfMET significance", nbins, 0., 400.)
    h_metdphisig = TH1F("%s_metdphisig" % (chain.GetTitle()), ";pfMET/(pfMET + sumET_{#Delta#varphi<0.6})", nbins, 0., 1.)
    h_najets = TH1F("%s_najets_lt" % (chain.GetTitle()), ";# add. jets", nbins, 0., 32.)
    
    for ievt in xrange(0,nevt):
        chain.LoadTree(ievt)  # used by TTreeFormula
        
        if chain.GetTreeNumber() != curTree:
            curTree = chain.GetTreeNumber()
            cut.UpdateFormulaLeaves()
            cutWeight.UpdateFormulaLeaves()
        
        cut.GetNdata()
        pass_cut = bool(cut.EvalInstance())
        evtweight = cutWeight.EvalInstance()
        
        if not pass_cut:
            continue
    
        chain.GetEntry(ievt)
    
        mindphi_c50 = 3.141593
        mindphi_c30 = 3.141593
        mindphi_c25 = 3.141593
        mindphi_c20 = 3.141593
        mindphi_cf50 = 3.141593
        mindphi_cf30 = 3.141593
        mindphi_cf25 = 3.141593
        mindphi_cf20 = 3.141593
        mindphi_c20f30 = 3.141593
        metdphisig_sumet = METtype1corr.et + 0.1
        
        for ihj in xrange(chain.nhJets):
            dphi = abs(deltaPhi(METtype1corr.phi, chain.hJet_phi[ihj]))
            if mindphi_c50 > dphi and chain.hJet_ptReg[ihj] > 50. and abs(chain.hJet_eta[ihj]) < 2.5:
                mindphi_c50 = dphi
            if mindphi_c30 > dphi and chain.hJet_ptReg[ihj] > 30. and abs(chain.hJet_eta[ihj]) < 2.5:
                mindphi_c30 = dphi
            if mindphi_c25 > dphi and chain.hJet_ptReg[ihj] > 25. and abs(chain.hJet_eta[ihj]) < 2.5:
                mindphi_c25 = dphi
            if mindphi_c20 > dphi and chain.hJet_ptReg[ihj] > 20. and abs(chain.hJet_eta[ihj]) < 2.5:
                mindphi_c20 = dphi
            if mindphi_cf50 > dphi and chain.hJet_ptReg[ihj] > 50. and abs(chain.hJet_eta[ihj]) < 5:
                mindphi_cf50 = dphi                                                            
            if mindphi_cf30 > dphi and chain.hJet_ptReg[ihj] > 30. and abs(chain.hJet_eta[ihj]) < 5:
                mindphi_cf30 = dphi                                                            
            if mindphi_cf25 > dphi and chain.hJet_ptReg[ihj] > 25. and abs(chain.hJet_eta[ihj]) < 5:
                mindphi_cf25 = dphi                                                            
            if mindphi_cf20 > dphi and chain.hJet_ptReg[ihj] > 20. and abs(chain.hJet_eta[ihj]) < 5:
                mindphi_cf20 = dphi
            if mindphi_c20f30 > dphi and chain.hJet_ptReg[ihj] > 20. and abs(chain.hJet_eta[ihj]) < 2.5:
                mindphi_c20f30 = dphi
            if dphi < 0.6:
                metdphisig_sumet += chain.hJet_ptReg[ihj]

        for iaj in xrange(chain.naJets):
            dphi = abs(deltaPhi(METtype1corr.phi, chain.aJet_phi[iaj]))
            if mindphi_c50 > dphi and chain.aJet_pt[iaj] > 50. and abs(chain.aJet_eta[iaj]) < 2.5:
                mindphi_c50 = dphi
            if mindphi_c30 > dphi and chain.aJet_pt[iaj] > 30. and abs(chain.aJet_eta[iaj]) < 2.5:
                mindphi_c30 = dphi
            if mindphi_c25 > dphi and chain.aJet_pt[iaj] > 25. and abs(chain.aJet_eta[iaj]) < 2.5:
                mindphi_c25 = dphi
            if mindphi_c20 > dphi and chain.aJet_pt[iaj] > 20. and abs(chain.aJet_eta[iaj]) < 2.5:
                mindphi_c20 = dphi
            if mindphi_cf50 > dphi and chain.aJet_pt[iaj] > 50. and abs(chain.aJet_eta[iaj]) < 5:
                mindphi_cf50 = dphi                                                            
            if mindphi_cf30 > dphi and chain.aJet_pt[iaj] > 30. and abs(chain.aJet_eta[iaj]) < 5:
                mindphi_cf30 = dphi                                                            
            if mindphi_cf25 > dphi and chain.aJet_pt[iaj] > 25. and abs(chain.aJet_eta[iaj]) < 5:
                mindphi_cf25 = dphi                                                            
            if mindphi_cf20 > dphi and chain.aJet_pt[iaj] > 20. and abs(chain.aJet_eta[iaj]) < 5:
                mindphi_cf20 = dphi
            if mindphi_c20f30 > dphi and chain.aJet_pt[iaj] > 30. and abs(chain.aJet_eta[iaj]) < 5:  # cj threshold is lower, so no need to check
                mindphi_c20f30 = dphi
            if dphi < 0.6:
                metdphisig_sumet += chain.aJet_pt[iaj]
        
        maxpt0     = chain.hJet_ptReg[0] if (chain.hJet_ptReg[0] > chain.hJet_ptReg[1]) else chain.hJet_ptReg[1]
        maxpt0_phi = chain.hJet_phi[0]   if (chain.hJet_ptReg[0] > chain.hJet_ptReg[1]) else chain.hJet_phi[1]
        maxpt1     = chain.hJet_ptReg[1] if (chain.hJet_ptReg[0] > chain.hJet_ptReg[1]) else chain.hJet_ptReg[0]
        maxpt1_phi = chain.hJet_phi[1]   if (chain.hJet_ptReg[0] > chain.hJet_ptReg[1]) else chain.hJet_phi[0]
        maxpt2     = -1.
        maxpt2_phi = 0.
        for iaj in xrange(chain.naJets):
            if maxpt0 < chain.aJet_pt[iaj]:
                maxpt2 = maxpt1
                maxpt2_phi = maxpt1_phi
                maxpt1 = maxpt0
                maxpt1_phi = maxpt0_phi
                maxpt0 = chain.aJet_pt[iaj]
                maxpt0_phi = chain.aJet_phi[iaj]
            elif maxpt1 < chain.aJet_pt[iaj]:
                maxpt2 = maxpt1
                maxpt2_phi = maxpt1_phi
                maxpt1 = chain.aJet_pt[iaj]
                maxpt1_phi = chain.aJet_phi[iaj]
            elif maxpt2 < chain.aJet_pt[iaj]:
                maxpt2 = chain.aJet_pt[iaj]
                maxpt2_phi = chain.aJet_phi[iaj]
        mindphi_j3 = min(abs(deltaPhi(METtype1corr.phi, maxpt0_phi)), abs(deltaPhi(METtype1corr.phi, maxpt1_phi)) )
        if maxpt2 > 20.:
            assert(maxpt2 < maxpt1 and maxpt1 < maxpt0)
            mindphi_j3 = min(mindphi_j3, abs(deltaPhi(METtype1corr.phi, maxpt2_phi)))
        
        h_mindphi_c50.Fill(mindphi_c50, evtweight)
        h_mindphi_c30.Fill(mindphi_c30, evtweight)
        h_mindphi_c25.Fill(mindphi_c25, evtweight)
        h_mindphi_c20.Fill(mindphi_c20, evtweight)
        h_mindphi_cf50.Fill(mindphi_cf50, evtweight)
        h_mindphi_cf30.Fill(mindphi_cf30, evtweight)
        h_mindphi_cf25.Fill(mindphi_cf25, evtweight)
        h_mindphi_cf20.Fill(mindphi_cf20, evtweight)
        h_mindphi_c20f30.Fill(mindphi_c20f30, evtweight)
        h_mindphi_j3.Fill(mindphi_j3, evtweight)
        h_mindphi.Fill(min(chain.mindPhiMETJet_dPhi,3.141593), evtweight)
        h_dphiTrk.Fill(abs(deltaPhi(METtype1corr.phi,METnoPUCh.phi)), evtweight)
        h_dphiMht.Fill(abs(deltaPhi(METtype1corr.phi,MHT.phi)), evtweight)
        h_alphaT.Fill(alphaT(maxpt0, maxpt0_phi, maxpt1, maxpt1_phi), evtweight)
        h_metsig.Fill(METtype1corr.et/sqrt(METtype1corr.sumet), evtweight)
        h_metfancysig.Fill(METtype1corr.sig, evtweight)
        h_metdphisig.Fill(METtype1corr.et/metdphisig_sumet, evtweight)
        h_najets.Fill(chain.naJets, evtweight)
        
    histograms.append(h_mindphi_c50)
    histograms.append(h_mindphi_c30)
    histograms.append(h_mindphi_c25)
    histograms.append(h_mindphi_c20)
    histograms.append(h_mindphi_cf50)
    histograms.append(h_mindphi_cf30)
    histograms.append(h_mindphi_cf25)
    histograms.append(h_mindphi_cf20)
    histograms.append(h_mindphi_c20f30)
    histograms.append(h_mindphi_j3)
    histograms.append(h_mindphi)
    histograms.append(h_dphiTrk)
    histograms.append(h_dphiMht)
    histograms.append(h_alphaT)
    histograms.append(h_metsig)
    histograms.append(h_metfancysig)
    histograms.append(h_metdphisig)
    histograms.append(h_najets)

histograms2 = []
for ih in xrange(len(histograms)/2):
    hZH = histograms[ih]
    hQCD = histograms[ih+len(histograms)/2]
    print hZH.GetName()
    print hQCD.GetName()
    sumZH = hZH.Integral(0,nbins+3)
    sumQCD = hQCD.Integral(0,nbins+3)
    
    hZH0 = hZH.Clone("hZH0")
    hQCD0 = hQCD.Clone("hQCD0")
    hZH0.Rebin(50); hZH0.SetFillColor(632); hZH0.SetLineColor(632); hZH0.Scale(1.0/sumZH)
    hQCD0.Rebin(50); hQCD0.SetLineWidth(3); hQCD0.SetLineColor(1); hQCD0.Scale(1.0/sumQCD)
    ymax = max(hZH0.GetMaximum(), hQCD0.GetMaximum()) * 1.1
    hZH0.SetMaximum(ymax)
    if "mindphi" in plotname:
        hZH0.SetMaximum(0.71)
    hZH0.Draw(); hQCD0.Draw("same")
    plotname = hZH.GetName()
    plotname = plotname.replace("ZH125_", "plots/antiQCD_%s_" % channel)
    #plotname = plotname.replace("ZH125_", "plots2/antiQCD_%s_" % channel)
    if "najets" in plotname:
        hZH0.GetXaxis().SetRangeUser(0,10)
    gPad.Modified()
    gPad.Update()
    gPad.Print(plotname+".png")
    
    hrejBvsS = TH1F("hrejBvsS_"+hZH.GetName(), "; Signal eff, effS; Background rej, 1 - effB", 150, 0.85, 1.0)
    for ibin in xrange(100):
        effS = hrejBvsS.GetBinCenter(ibin+1)
        value = hZH.GetBinCenter(nbins+1)
        for jbin in xrange(nbins+1):  # want to scan [1:nbins+1]
            if not plotname.endswith("_lt"):
                value = 1
                if hZH.Integral(nbins+1-jbin, nbins+2)/sumZH > effS: # scan [(nbins,nbins+2):(0,nbins+2)]
                    value = hZH.GetBinCenter(nbins+1-jbin)
                    break
            else:
                value = nbins+1
                if hZH.Integral(0, jbin+1)/sumZH > effS: # scan [(0,2):(0:nbins+2)]
                    value = hZH.GetBinCenter(jbin+1)
                    break
        
        effB = 1.0
        jbin = hQCD.FindFixBin(value)
        if not plotname.endswith("_lt"):
            effB = hQCD.Integral(jbin, nbins+2)/sumQCD
        else:
            effB = hQCD.Integral(0, jbin)/sumQCD
        hrejBvsS.SetBinContent(ibin+1, 1.0-effB)
    
    hrejBvsS.SetMaximum(1.1); hrejBvsS.SetMinimum(0)
    hrejBvsS.SetStats(0); hrejBvsS.SetLineWidth(2); hrejBvsS.GetXaxis().SetNdivisions(505); hrejBvsS.Draw()
    gPad.Print(plotname+"_rejBvsS"+".png")
    histograms2.append(hrejBvsS)



#histograms[1].Rebin(50)
#histograms[9].Rebin(50)
#histograms[1].Draw()
#histograms[9].SetLineColor(2)
#histograms[9].Draw("same")
#gPad.Print("temp.png")
