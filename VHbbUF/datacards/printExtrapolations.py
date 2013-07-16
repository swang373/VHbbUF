from ROOT import TFile, TH1, TH1F

processes = ["Wj0b", "Wj1b", "Wj2b", "Zj0b", "Zj1b", "Zj2b", "TT"]
systematics = ["CMS_vhbb_res_j", "CMS_vhbb_scale_j", "CMS_vhbb_eff_b", "CMS_vhbb_fake_b_8TeV"]

norms = {}

for channel in ["ZnunuHighPt_ZjLF", "ZnunuHighPt_ZjHF", "ZnunuHighPt_WjLF", "ZnunuHighPt_WjHF", "ZnunuHighPt_TT"]:

    print channel
    f1 = TFile.Open("vhbb_Znn_J11_ZnunuHighPt_TH1.root")
    f2 = TFile.Open("vhbb_Znn_SF_J11_%s_TH1.root" % channel)

    for s in systematics:
        numbers = []
        for p in processes:
            if "Zj" in p and not "Zj" in channel:  continue
            h1 = f1.Get("ZnunuHighPt/"+p)
            h1u = f1.Get("ZnunuHighPt/"+p+"_"+s+"Up")
            h1d = f1.Get("ZnunuHighPt/"+p+"_"+s+"Down")
            h2 = f2.Get(channel+"/"+p)
            h2u = f2.Get(channel+"/"+p+"_"+s+"Up")
            h2d = f2.Get(channel+"/"+p+"_"+s+"Down")
            numbers.append( (h1u.Integral() / h1.Integral(), h1d.Integral() / h1.Integral(), h2u.Integral() / h2.Integral(), h2d.Integral() / h2.Integral()) )

        what = "SF CR"
        if what == "SF SR":
            print "%28s" % (s+"Up"),
            for i, n in enumerate(numbers):
                print "%.3f   " % (n[0]),
            print
            print "%28s" % (s+"Down"),
            for i, n in enumerate(numbers):
                print "%.3f   " % (n[1]),
            print
        elif what == "SF CR":
            print "%28s" % (s+"Up"),
            for i, n in enumerate(numbers):
                print "%.3f   " % (n[2]),
            print
            print "%28s" % (s+"Down"),
            for i, n in enumerate(numbers):
                print "%.3f   " % (n[3]),
            print
        else:
            print "%28s" % (s+"Up"),
            for i, n in enumerate(numbers):
                print "%.3f   " % (n[0]/n[2]),
            print
            print "%28s" % (s+"Down"),
            for i, n in enumerate(numbers):
                print "%.3f   " % (n[1]/n[3]),
            print


        del numbers[:]
