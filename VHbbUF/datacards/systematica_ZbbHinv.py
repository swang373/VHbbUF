import os
import re
from ROOT import TFile

# Use these commands to extract the list of nuisances
#sed -n 's/\(^\S\+\).*\(shape\|lnN\).*/("\1","\2",""),/p' zhinv_Zbb_8TeV.txt  > ! temp

nuisances = [
("FULL","shape","FULL"),  #! KEEP THIS AS THE FIRST
("CMS_scale_met","lnN","QCD"),
#("CMS_vhbb_BR_Hbb","lnN","VH2"),  # should be background
("CMS_vhbb_QCD_SF_ZnunuHighPt_8TeV","lnN","QCD"),
("CMS_vhbb_QCD_SF_ZnunuLowPt_8TeV","lnN","QCD"),
("CMS_vhbb_QCD_SF_ZnunuMedPt_8TeV","lnN","QCD"),
("CMS_vhbb_TTModel_Znn_8TeV","shape","MODEL"),
("CMS_vhbb_TT_SF_ZnunuHighPt_8TeV","lnN","SF"),
("CMS_vhbb_TT_SF_ZnunuLowPt_8TeV","lnN","SF"),
("CMS_vhbb_TT_SF_ZnunuMedPt_8TeV","lnN","SF"),
("CMS_vhbb_WJSlope_Znn_8TeV","shape","MODEL"),
("CMS_vhbb_Wj0b_SF_ZnunuHighPt_8TeV","lnN","SF"),
("CMS_vhbb_Wj0b_SF_ZnunuLowPt_8TeV","lnN","SF"),
("CMS_vhbb_Wj0b_SF_ZnunuMedPt_8TeV","lnN","SF"),
("CMS_vhbb_Wj1b_SF_ZnunuHighPt_8TeV","lnN","SF"),
("CMS_vhbb_Wj1b_SF_ZnunuLowPt_8TeV","lnN","SF"),
("CMS_vhbb_Wj1b_SF_ZnunuMedPt_8TeV","lnN","SF"),
("CMS_vhbb_Wj2b_SF_ZnunuHighPt_8TeV","lnN","SF"),
("CMS_vhbb_Wj2b_SF_ZnunuLowPt_8TeV","lnN","SF"),
("CMS_vhbb_Wj2b_SF_ZnunuMedPt_8TeV","lnN","SF"),
("CMS_vhbb_ZJSlope_Znn_8TeV","shape","MODEL"),
("CMS_vhbb_Zj0b_SF_ZnunuHighPt_8TeV","lnN","SF"),
("CMS_vhbb_Zj0b_SF_ZnunuLowPt_8TeV","lnN","SF"),
("CMS_vhbb_Zj0b_SF_ZnunuMedPt_8TeV","lnN","SF"),
("CMS_vhbb_Zj1b_SF_ZnunuHighPt_8TeV","lnN","SF"),
("CMS_vhbb_Zj1b_SF_ZnunuLowPt_8TeV","lnN","SF"),
("CMS_vhbb_Zj1b_SF_ZnunuMedPt_8TeV","lnN","SF"),
("CMS_vhbb_Zj2bModel_Znn_8TeV","shape","MODEL"),
("CMS_vhbb_Zj2b_SF_ZnunuHighPt_8TeV","lnN","SF"),
("CMS_vhbb_Zj2b_SF_ZnunuLowPt_8TeV","lnN","SF"),
("CMS_vhbb_Zj2b_SF_ZnunuMedPt_8TeV","lnN","SF"),
("CMS_vhbb_Znn_res_j","shape","JER"),
("CMS_vhbb_Znn_scale_j","shape","JES"),
("CMS_vhbb_boost_EWK_8TeV","lnN","VH"),
("CMS_vhbb_boost_QCD_8TeV","lnN","VH"),
("CMS_vhbb_eff_b","shape","BTAG"),
("CMS_vhbb_fake_b_8TeV","shape","BTAG"),
("CMS_vhbb_s_Top","lnN","ST"),
("CMS_vhbb_statQCD_ZnunuHighPt_bin1_8TeV","shape","STAT"),
("CMS_vhbb_statQCD_ZnunuHighPt_bin2_8TeV","shape","STAT"),
("CMS_vhbb_statQCD_ZnunuHighPt_bin3_8TeV","shape","STAT"),
("CMS_vhbb_statQCD_ZnunuLowPt_bin1_8TeV","shape","STAT"),
("CMS_vhbb_statQCD_ZnunuLowPt_bin2_8TeV","shape","STAT"),
("CMS_vhbb_statQCD_ZnunuLowPt_bin3_8TeV","shape","STAT"),
("CMS_vhbb_statQCD_ZnunuLowPt_bin4_8TeV","shape","STAT"),
("CMS_vhbb_statQCD_ZnunuLowPt_bin7_8TeV","shape","STAT"),
("CMS_vhbb_statQCD_ZnunuMedPt_bin14_8TeV","shape","STAT"),
("CMS_vhbb_statQCD_ZnunuMedPt_bin15_8TeV","shape","STAT"),
("CMS_vhbb_statQCD_ZnunuMedPt_bin1_8TeV","shape","STAT"),
("CMS_vhbb_statQCD_ZnunuMedPt_bin2_8TeV","shape","STAT"),
("CMS_vhbb_statQCD_ZnunuMedPt_bin3_8TeV","shape","STAT"),
("CMS_vhbb_statQCD_ZnunuMedPt_bin4_8TeV","shape","STAT"),
("CMS_vhbb_statQCD_ZnunuMedPt_bin5_8TeV","shape","STAT"),
("CMS_vhbb_statQCD_ZnunuMedPt_bin6_8TeV","shape","STAT"),
("CMS_vhbb_statQCD_ZnunuMedPt_bin7_8TeV","shape","STAT"),
("CMS_vhbb_statTT_ZnunuHighPt_bin16_8TeV","shape","STAT"),
("CMS_vhbb_statTT_ZnunuHighPt_bin17_8TeV","shape","STAT"),
("CMS_vhbb_statTT_ZnunuMedPt_bin20_8TeV","shape","STAT"),
("CMS_vhbb_statVVLF_ZnunuHighPt_bin13_8TeV","shape","STAT"),
("CMS_vhbb_statVVLF_ZnunuHighPt_bin14_8TeV","shape","STAT"),
("CMS_vhbb_statVVLF_ZnunuHighPt_bin15_8TeV","shape","STAT"),
("CMS_vhbb_statVVLF_ZnunuHighPt_bin16_8TeV","shape","STAT"),
("CMS_vhbb_statVVLF_ZnunuHighPt_bin17_8TeV","shape","STAT"),
("CMS_vhbb_statVVLF_ZnunuHighPt_bin18_8TeV","shape","STAT"),
("CMS_vhbb_statVVLF_ZnunuLowPt_bin15_8TeV","shape","STAT"),
("CMS_vhbb_statVVLF_ZnunuLowPt_bin16_8TeV","shape","STAT"),
("CMS_vhbb_statVVLF_ZnunuLowPt_bin17_8TeV","shape","STAT"),
("CMS_vhbb_statVVLF_ZnunuLowPt_bin18_8TeV","shape","STAT"),
("CMS_vhbb_statVVLF_ZnunuLowPt_bin19_8TeV","shape","STAT"),
("CMS_vhbb_statVVLF_ZnunuMedPt_bin17_8TeV","shape","STAT"),
("CMS_vhbb_statVVLF_ZnunuMedPt_bin18_8TeV","shape","STAT"),
("CMS_vhbb_statWH_SM_ZnunuLowPt_bin19_8TeV","shape","STAT"),
("CMS_vhbb_statWj0b_ZnunuHighPt_bin16_8TeV","shape","STAT"),
("CMS_vhbb_statWj0b_ZnunuHighPt_bin19_8TeV","shape","STAT"),
("CMS_vhbb_statWj0b_ZnunuLowPt_bin19_8TeV","shape","STAT"),
("CMS_vhbb_statWj0b_ZnunuMedPt_bin19_8TeV","shape","STAT"),
("CMS_vhbb_statWj1b_ZnunuHighPt_bin14_8TeV","shape","STAT"),
("CMS_vhbb_statWj1b_ZnunuHighPt_bin15_8TeV","shape","STAT"),
("CMS_vhbb_statWj1b_ZnunuHighPt_bin16_8TeV","shape","STAT"),
("CMS_vhbb_statWj1b_ZnunuHighPt_bin17_8TeV","shape","STAT"),
("CMS_vhbb_statWj1b_ZnunuHighPt_bin19_8TeV","shape","STAT"),
("CMS_vhbb_statWj1b_ZnunuLowPt_bin13_8TeV","shape","STAT"),
("CMS_vhbb_statWj1b_ZnunuLowPt_bin14_8TeV","shape","STAT"),
("CMS_vhbb_statWj1b_ZnunuLowPt_bin15_8TeV","shape","STAT"),
("CMS_vhbb_statWj1b_ZnunuLowPt_bin16_8TeV","shape","STAT"),
("CMS_vhbb_statWj1b_ZnunuLowPt_bin17_8TeV","shape","STAT"),
("CMS_vhbb_statWj1b_ZnunuLowPt_bin18_8TeV","shape","STAT"),
("CMS_vhbb_statWj1b_ZnunuMedPt_bin13_8TeV","shape","STAT"),
("CMS_vhbb_statWj1b_ZnunuMedPt_bin14_8TeV","shape","STAT"),
("CMS_vhbb_statWj1b_ZnunuMedPt_bin15_8TeV","shape","STAT"),
("CMS_vhbb_statWj1b_ZnunuMedPt_bin16_8TeV","shape","STAT"),
("CMS_vhbb_statWj1b_ZnunuMedPt_bin17_8TeV","shape","STAT"),
("CMS_vhbb_statWj1b_ZnunuMedPt_bin19_8TeV","shape","STAT"),
("CMS_vhbb_statWj2b_ZnunuHighPt_bin19_8TeV","shape","STAT"),
("CMS_vhbb_statWj2b_ZnunuHighPt_bin20_8TeV","shape","STAT"),
("CMS_vhbb_statWj2b_ZnunuLowPt_bin19_8TeV","shape","STAT"),
("CMS_vhbb_statWj2b_ZnunuLowPt_bin20_8TeV","shape","STAT"),
("CMS_vhbb_statWj2b_ZnunuMedPt_bin20_8TeV","shape","STAT"),
("CMS_vhbb_statZj0b_ZnunuHighPt_bin18_8TeV","shape","STAT"),
("CMS_vhbb_statZj0b_ZnunuLowPt_bin19_8TeV","shape","STAT"),
("CMS_vhbb_statZj0b_ZnunuLowPt_bin20_8TeV","shape","STAT"),
("CMS_vhbb_statZj0b_ZnunuMedPt_bin18_8TeV","shape","STAT"),
("CMS_vhbb_statZj0b_ZnunuMedPt_bin19_8TeV","shape","STAT"),
("CMS_vhbb_statZj1b_ZnunuHighPt_bin16_8TeV","shape","STAT"),
("CMS_vhbb_statZj1b_ZnunuHighPt_bin17_8TeV","shape","STAT"),
("CMS_vhbb_statZj1b_ZnunuHighPt_bin18_8TeV","shape","STAT"),
("CMS_vhbb_statZj1b_ZnunuLowPt_bin15_8TeV","shape","STAT"),
("CMS_vhbb_statZj1b_ZnunuLowPt_bin16_8TeV","shape","STAT"),
("CMS_vhbb_statZj1b_ZnunuLowPt_bin19_8TeV","shape","STAT"),
("CMS_vhbb_statZj1b_ZnunuMedPt_bin17_8TeV","shape","STAT"),
("CMS_vhbb_statZj1b_ZnunuMedPt_bin18_8TeV","shape","STAT"),
("CMS_vhbb_statZj1b_ZnunuMedPt_bin19_8TeV","shape","STAT"),
("CMS_vhbb_statZj2b_ZnunuHighPt_bin20_8TeV","shape","STAT"),
("CMS_vhbb_statZj2b_ZnunuMedPt_bin20_8TeV","shape","STAT"),
("CMS_vhbb_stats_Top_ZnunuHighPt_bin13_8TeV","shape","STAT"),
("CMS_vhbb_stats_Top_ZnunuHighPt_bin15_8TeV","shape","STAT"),
("CMS_vhbb_stats_Top_ZnunuHighPt_bin16_8TeV","shape","STAT"),
("CMS_vhbb_stats_Top_ZnunuLowPt_bin18_8TeV","shape","STAT"),
("CMS_vhbb_stats_Top_ZnunuMedPt_bin15_8TeV","shape","STAT"),
("CMS_vhbb_stats_Top_ZnunuMedPt_bin18_8TeV","shape","STAT"),
("CMS_vhbb_stats_Top_ZnunuMedPt_bin19_8TeV","shape","STAT"),
("CMS_vhbb_stats_Top_ZnunuMedPt_bin20_8TeV","shape","STAT"),
("CMS_vhbb_statwz3lnu_ZnunuHighPt_bin16_8TeV","shape","STAT"),
("CMS_vhbb_statwz3lnu_ZnunuHighPt_bin18_8TeV","shape","STAT"),
("CMS_vhbb_statwz3lnu_ZnunuHighPt_bin19_8TeV","shape","STAT"),
("CMS_vhbb_statwz3lnu_ZnunuLowPt_bin20_8TeV","shape","STAT"),
("CMS_vhbb_statwz3lnu_ZnunuMedPt_bin20_8TeV","shape","STAT"),
("CMS_vhbb_trigger_CSV_Znn_8TeV","shape","MET"),
("CMS_vhbb_trigger_CSV_fake_Znn_8TeV","shape","MET"),
("CMS_vhbb_trigger_MET_Znn_8TeV","shape","MET"),
("QCDscale_VH","lnN","VH2"),
("QCDscale_VV","lnN","VV"),
("lumi_8TeV","lnN","LUMI"),
("pdf_qqbar","lnN","VH2"),
]

channels = ["ZnunuHighPt", "ZnunuMedPt", "ZnunuLowPt"]
orig = "zhinv_Zbb_J14_bbb_$CHANNEL_8TeV.txt"
#orig = "zhinv_Zbb_8TeV.txt"
test = "test_8TeV.txt"
wsname = "zhinv_Zbb_8TeV.root"
realvarname = "zhinv_Zbb_BDT_$CHANNEL_8TeV"
#realvarname = "zhinv_Zbb_MJJ_$CHANNEL_8TeV"
soverb = 0.035
#soverb = -99
verbose = True
mode = "YIELD_UNCERTAINTY"

rootfile = TFile.Open(wsname)
#processes = "ZH_SM      WH_SM      zh1252lmet Wj0b       Wj1b       Wj2b       Zj0b       Zj1b       Zj2b       TT         s_Top      VVLF       zz2l2nu    wz3lnu     QCD".split()
processes = "ZH_hbb     WH_hbb     ZH_hinv    Wj0b       Wj1b       Wj2b       Zj0b       Zj1b       Zj2b       TT         s_Top      VVLF       ZZ         WZ         QCD".split()
#signal = "zh1252lmet"
signal = "ZH_hinv"


#-------------------------------------------------------------------------------

assert(nuisances[0][0] == "FULL")
nuisances_1 = []  # nuisance names
nuisances_2 = []  # nuisance types
nuisances_3 = []  # nuisance groups
for n1, n2, n3 in nuisances:
    if n3 in nuisances_3:
        nuisances_1[nuisances_3.index(n3)].append(n1)
    else:
        nuisances_1.append([n1])
        nuisances_2.append(n2)
        nuisances_3.append(n3)

if verbose:
    for i in xrange(len(nuisances_1)):
        print nuisances_3[i]
        print nuisances_1[i]

strlog1 = []
strlog2 = []


def sumStatErr(h):
    # Return sqrt of sum of stat errors squared

    sumoferrors2 = 0
    for b in xrange(h.GetNbinsX()+2):
        if h.GetBinContent(b) == 0:  continue
        sumoferrors2 += (h.GetBinError(b) * h.GetBinError(b))
    return sumoferrors2 ** 0.5


# Loop over each channel separately
for channel in channels:
    channel_8TeV = channel + "_8TeV"
    print channel_8TeV

    # Calculate sumAll, sumBkg, sumSig, then calculate S over B
    hmap = {}
    for i in xrange(len(processes)):
        h = getattr(rootfile, channel_8TeV).data("%s" % (processes[i])).createHistogram(realvarname.replace("$CHANNEL", channel))
        hmap[processes[i]] = h.Clone(processes[i])
    sumAll = hmap[signal].Clone("sumAll")
    sumBkg = hmap[signal].Clone("sumBkg")
    sumSig = hmap[signal].Clone("sumSig")
    sumAll.Reset()
    sumBkg.Reset()
    for k, v in hmap.iteritems():
        if v.Integral() == 0:  continue
        sumAll.Add(v)
        if k != signal:  sumBkg.Add(v)
    sumSoB = sumSig.Clone("sumSoB")
    sumSoB.Divide(sumBkg)

    # Find first bin in sumSoB whose bin content is < `soverb`
    firstbin = 0
    for b in xrange(sumSig.GetNbinsX()+2, 0, -1):
        bc = sumSoB.GetBinContent(b)
        #if bc != 0 and (bc < soverb or (channel == "ZnunuLowPt" and bc < (soverb * 3.5/5))):  # FIXME
        if bc != 0 and bc < soverb:
            firstbin = b
            break

    firstbin += 1  # increment it to the first bin whose bin content is >= `soverb`
    print "first bin: ", firstbin, sumSoB.GetBinContent(firstbin)
    print "integral1: ", sumBkg.Integral(), sumSig.Integral()
    print "integral2: ", sumBkg.Integral(firstbin, sumBkg.GetNbinsX()+2), sumSig.Integral(firstbin, sumSig.GetNbinsX()+2)
    print
    #print "%-6s %-48s, %.2f, %.2f" % ("STAT", "CMS_vhbb_statSIG_" + channel_8TeV + " (shape)", sumStatErr(sumSig) / sumSig.Integral(), sumStatErr(sumSig) / sumSig.Integral())
    #print "%-6s %-48s, %.2f, %.2f" % ("STAT", "CMS_vhbb_statBKG_" + channel_8TeV + " (shape)", sumStatErr(sumBkg) / sumBkg.Integral(), sumStatErr(sumBkg) / sumBkg.Integral())


    for n1, n2, n3 in zip(nuisances_1, nuisances_2, nuisances_3):
        # Copy to a temporary datacard
        os.system("cp " + orig.replace("$CHANNEL", channel) + " " + test)
        # Uses "kmax *"
        os.system("sed -i 's/^kmax.*/kmax */' " + test)

        upGrpBkg, dnGrpBkg, nmGrpBkg = 0., 0., 0.
        upGrpSig, dnGrpSig, nmGrpSig = 0., 0., 0.

        with open(test,"r") as f:
            if mode == "YIELD_UNCERTAINTY":
                for line in f:
                    if ("shape" in line or "lnN" in line) and (".root" not in line):
                        for n in n1:
                            if n in line and n != "FULL" and n3 != "STAT":
                                upBkg, dnBkg, nmBkg = 0., 0., 0.
                                upSig, dnSig, nmSig = 0., 0., 0.
                                # Regex to extract n=nprocess+2 columns
                                m = re.match(r"(\S+)\s*"*(len(processes)+2),line)
                                if m:
                                    # Loop over the columns
                                    # First two columns (nui name, type) are not numerical
                                    for i, mi in enumerate(m.groups()[2:]):
                                        hnm = hmap[processes[i]]
                                        nm = hnm.Integral(firstbin, hnm.GetNbinsX()+2)

                                        if mi.replace(".","").isdigit():
                                            if m.groups()[1] == "lnN":
                                                sf = float(mi)
                                                if sf < 1: sf = 2. - sf

                                                if processes[i] != signal:
                                                    nmBkg += nm
                                                    upBkg += (nm * sf)
                                                    dnBkg += (nm * (2. - sf))
                                                else:
                                                    nmSig += nm
                                                    upSig += (nm * sf)
                                                    dnSig += (nm * (2. - sf))

                                            elif m.groups()[1] == "shape":
                                                nnui = m.groups()[0]
                                                hup = getattr(rootfile, channel_8TeV).data("%s_%sUp" % (processes[i], nnui)).createHistogram(realvarname.replace("$CHANNEL", channel))
                                                hdn = getattr(rootfile, channel_8TeV).data("%s_%sDown" % (processes[i], nnui)).createHistogram(realvarname.replace("$CHANNEL", channel))

                                                up = hup.Integral(firstbin, hup.GetNbinsX()+2)
                                                dn = hdn.Integral(firstbin, hdn.GetNbinsX()+2)
                                                # For these nuisances, 'up' reduces yield
                                                if "res_j" in nnui or "scale_j" in nnui or "fake_b" in nnui:
                                                    up, dn = dn, up
                                                    #print nnui, processes[i], up, dn, nm, upBkg/(nmBkg if nmBkg>0 else 1e-6)

                                                if processes[i] != signal:
                                                    nmBkg += nm
                                                    upBkg += up
                                                    dnBkg += dn
                                                else:
                                                    nmSig += nm
                                                    upSig += up
                                                    dnSig += dn

                                        # Not affected by the particular systematics
                                        else:
                                            if processes[i] != signal:
                                                nmBkg += nm
                                                upBkg += nm
                                                dnBkg += nm
                                            else:
                                                nmSig += nm
                                                upSig += nm
                                                dnSig += nm

                                if nmBkg == 0.:  nmBkg = 1e-6
                                if nmSig == 0.:  nmSig = 1e-6
                                print "%-6s %-48s, %.3f, %.3f, %.3f, %.3f" % (n3, n+" ("+n2+")", upBkg / nmBkg, dnBkg / nmBkg, upSig / nmSig, dnSig / nmSig)

                                if nmGrpSig == 0.:  # initialize to nm, then add the diff
                                    nmGrpSig += nmSig
                                    upGrpSig += nmSig
                                    dnGrpSig += nmSig
                                upGrpSig += (upSig - nmSig)
                                dnGrpSig += (dnSig - nmSig)

                                if nmGrpBkg == 0.:  # initialize to nm, then add the diff
                                    nmGrpBkg += nmBkg
                                    upGrpBkg += nmBkg
                                    dnGrpBkg += nmBkg
                                upGrpBkg += (upBkg - nmBkg)
                                dnGrpBkg += (dnBkg - nmBkg)

        if nmGrpBkg == 0:  nmGrpBkg = 1e-6
        if nmGrpSig == 0:  nmGrpSig = 1e-6
        if n3 != "STAT":
            print "%-6s %-48s, %.3f, %.3f, %.3f, %.3f" % (n3, "TOTAL", upGrpBkg / nmGrpBkg, dnGrpBkg / nmGrpBkg, upGrpSig / nmGrpSig, dnGrpSig / nmGrpSig)
        else:
            print "%-6s %-48s, %.3f, %.3f, %.3f, %.3f" % ("STAT", "TOTAL", 1.0 + sumStatErr(sumBkg) / sumBkg.Integral(), 1.0 + sumStatErr(sumBkg) / sumBkg.Integral(), 1.0 + sumStatErr(sumSig) / sumSig.Integral(), 1.0 + sumStatErr(sumSig) / sumSig.Integral())

    print
