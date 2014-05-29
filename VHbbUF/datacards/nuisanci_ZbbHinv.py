import os
import re

# Use these commands to extract the list of nuisances
#sed -n 's/\(^\S\+\).*\(shape\|lnN\).*/("\1","\2",""),/p' zhinv_Zbb_8TeV.txt  > ! temp

# For stat only mu uncertainty
#combine -M MaxLikelihoodFit -m 125 --stepSize=0.05 --rMin=-5 --rMax=5 -t -1 --expectSignal=1 -S 0 vhbb_Znn_8TeV.txt

nuisances = [
("FULL","shape","FULL"),  #! KEEP THIS AS THE FIRST
("CMS_scale_met","lnN","QCD"),
("CMS_vhbb_BR_Hbb","lnN","VH2"),
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

orig = "zhinv_Zbb_8TeV.txt"
#orig = "vhbb_Znn_8TeV.txt"
#orig = "vhbb_VH_7p8TeV.txt"
test = "test_8TeV.txt"
rootfile = None
wsname = ""
processes = ""
comb1 = "combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 -t -1 --expectSignal=1 "
comb2 = "combine -M ProfileLikelihood -m 125 --signif --pvalue -t -1 --expectSignal=1 "
log1 = "syst_mu.log"
log2 = "syst_sig.log"
verbose = True
mode = "ALL_BUT_ONE"
#mode = "ONE_AT_A_TIME"
#mode = "YIELD_UNCERTAINTY"  # this should be done individually for each datacard
if mode == "YIELD_UNCERTAINTY":
    orig = "vhbb_Znn_J14_bbb_ZnunuHighPt_8TeV.txt"
    wsname = "ZnunuHighPt_8TeV"
    #orig = "vhbb_Znn_J14_bbb_ZnunuMedPt_8TeV.txt"
    #wsname = "ZnunuMedPt_8TeV"
    #orig = "vhbb_Znn_J14_bbb_ZnunuLowPt_8TeV.txt"
    #wsname = "ZnunuLowPt_8TeV"
    from ROOT import TFile
    rootfile = TFile.Open("zhinv_Zbb_8TeV.root")
    processes = "ZH_SM      WH_SM      ZH         Wj0b       Wj1b       Wj2b       Zj0b       Zj1b       Zj2b       TT         s_Top      VVLF       ZZ         WZ         QCD".split()



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

for n1, n2, n3 in zip(nuisances_1, nuisances_2, nuisances_3):
    nui = n3

    # Copy to a temporary datacard
    os.system("cp " + orig + " " + test)
    # Uses "kmax *"
    os.system("sed -i 's/^kmax.*/kmax */' " + test)

    # Write the temporary datacard
    writeme = []
    with open(test,"r") as f:
        if mode == "ALL_BUT_ONE":
            # For each line, loops over the nuisances
            # If matched, comment it out
            for line in f:
                found = False
                for n in n1:
                    if n in line:
                        found = True
                if found:
                    line = "#" + line
                writeme.append(line)

        elif mode == "ONE_AT_A_TIME":
            # For each line that contains "shape" or "lnN", loops over the nuisances
            # If not matched, comment it out
            for line in f:
                if ("shape" in line or "lnN" in line) and (".root" not in line):
                    found = False
                    for n in n1:
                        if n in line or n == "FULL":
                            found = True
                    if not found:
                        line = "#" + line
                writeme.append(line)

        elif mode == "YIELD_UNCERTAINTY":

            # For each line that contains "shape" or "lnN", loops over the nuisances
            # If matched, calculate the changes in yield
            for line in f:
                if ("shape" in line or "lnN" in line) and (".root" not in line):
                    for n in n1:
                        if n in line and n != "FULL" and n3 != "STAT":
                            up, dn = 1.00, 1.00
                            # Regex to extract n=nprocess+2 columns
                            m = re.match(r"(\S+)\s*"*(len(processes)+2),line)
                            if m:
                                # Loop over the columns, picking the largest column
                                # First two columns (nui name, type) are not numerical
                                for i,mi in enumerate(m.groups()[2:]):
                                    if mi.replace(".","").isdigit():
                                        if m.groups()[1] == "lnN":
                                            up = max(up, float(mi))
                                            dn = max(dn, float(mi))
                                        elif m.groups()[1] == "shape":
                                            up_ = getattr(rootfile,wsname).data("%s_%sUp" % (processes[i], m.groups()[0])).sumEntries()
                                            dn_ = getattr(rootfile,wsname).data("%s_%sDown" % (processes[i], m.groups()[0])).sumEntries()
                                            up_ /= getattr(rootfile,wsname).data("%s" % (processes[i])).sumEntries()
                                            dn_ /= getattr(rootfile,wsname).data("%s" % (processes[i])).sumEntries()
                                            if up_ < 1: up_ = 2. - up_
                                            if dn_ < 1: dn_ = 2. - dn_
                                            #if (n3 == "JER" or n3 == "JES" or n3 == "BTAG"):
                                            #    print n3, processes[i], up_, dn_
                                            if (n3 == "JER" or n3 == "JES" or n3 == "BTAG") and \
                                               (processes[i] == "s_Top" or processes[i] == "VVLF" or processes[i] == "ZH_SM" or processes[i] == "WH_SM"):  # no stats
                                                up = up
                                                dn = dn
                                            else:
                                                up = max(up, up_)
                                                dn = max(dn, dn_)
                            print "%-6s %-48s, %.2f, %.2f" % (n3, n+" ("+n2+")", up, dn)

            #ZnunuHighPt_8TeV->data("ZH_SM")->sumEntries()


    with open(test,"w") as f:
        # Rewrite the temporary datacard
        f.write('\n'.join(writeme))

    if mode == "YIELD_UNCERTAINTY":
        continue

    # Execute combine
    os.system(comb1 + test + " > " + log1)
    os.system(comb2 + test + " > " + log2)

    # Prepare results
    prefix = "NO " if mode == "ALL_BUT_ONE" else "ONLY "
    with open(log1) as f:
        for line in f:
            m = re.match("Best fit r: (.*)  (.*)/(.*)  \(68% CL\)", line)
            if m:
                sigma = (float(m.group(2)) * -1. + float(m.group(3))) / 2
                if verbose:
                    print nui, m.group(2), m.group(3), sigma
                if not strlog1:
                    strlog1.append((nui, sigma, sigma/sigma))  # first entry, assumed to be "FULL"
                else:
                    strlog1.append((prefix+nui, sigma, sigma/strlog1[0][1]))  # the rest, relative to the first entry

    with open(log2) as f:
        for line in f:
            m = re.match("       \(Significance = (.*)\)", line)
            if m:
                signif = float(m.group(1))
                if verbose:
                    print nui, m.group(1)
                if not strlog2:
                    strlog2.append((nui, signif, signif/signif))  # first entry, assumed to be "FULL"
                else:
                    strlog2.append((prefix+nui, signif, signif/strlog2[0][1]))  # the rest, relative to the first entry

# Print results
print "MODE:", mode
print "metric: mu uncertainty"
strlog1.sort(key=lambda tup: tup[2], reverse=True)
for q in strlog1:
    print "%-40s, %.4f, %.4f" % (q[0], q[1], q[2]*100.)

print "metric: significance"
strlog2.sort(key=lambda tup: tup[2])
for q in strlog2:
    print "%-40s, %.4f, %.4f" % (q[0], q[1], q[2]*100.)
