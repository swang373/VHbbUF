import os

step = "PREPARE"
#step = "REPORT"

#masses = [125]
masses = range(105, 145+10, 10)

ref = "./125_ref"
combinecards = "run_combinecards.csh"

unblind = True


if step == "PREPARE":
    for mass in masses:
        os.system("mkdir -p %i" % mass)
        os.system("cp %s/%s %i" % (ref, combinecards, mass))
        os.system("sed -i 's/125/%i/g' %i/%s" % (mass, mass, combinecards))
    
        writeme = "cd %i\n" % mass
        writeme += "source %s\n" % combinecards
        if not unblind:
            writeme += "combine -M Asymptotic -t -1 zhinv_Zbb_8TeV.txt >! aa1.log\n"
            writeme += "combine -M Asymptotic -t -1 zhinv_Zll_7p8TeV.txt >! aa2.log\n"
            writeme += "combine -M Asymptotic -t -1 zhinv_combo_7p8TeV.txt >! aa3.log\n"
            writeme += "combine -M ProfileLikelihood -m %i --signif --pvalue -t -1 --toysFreq --expectSignal=1 zhinv_Zbb_8TeV.txt >! pp1.log\n" % mass
            writeme += "combine -M ProfileLikelihood -m %i --signif --pvalue -t -1 --toysFreq --expectSignal=1 zhinv_Zll_7p8TeV.txt >! pp2.log\n" % mass
            writeme += "combine -M ProfileLikelihood -m %i --signif --pvalue -t -1 --toysFreq --expectSignal=1 zhinv_combo_7p8TeV.txt >! pp3.log\n" % mass
        else:
            writeme += "combine -M Asymptotic zhinv_Zbb_8TeV.txt >! AA1.log\n"
            writeme += "combine -M Asymptotic zhinv_Zll_7p8TeV.txt >! AA2.log\n"
            writeme += "combine -M Asymptotic zhinv_combo_7p8TeV.txt >! AA3.log\n"
            writeme += "combine -M ProfileLikelihood -m %i --signif --pvalue --minimizerTolerance=1 zhinv_Zbb_8TeV.txt >! PP1.log\n" % mass
            writeme += "combine -M ProfileLikelihood -m %i --signif --pvalue --minimizerTolerance=1 zhinv_Zll_7p8TeV.txt >! PP2.log\n" % mass
            writeme += "combine -M ProfileLikelihood -m %i --signif --pvalue --minimizerTolerance=1 zhinv_combo_7p8TeV.txt >! PP3.log\n" % mass
        writeme += "cd -\n"
        with open("calculate_%i.csh" % mass, 'w') as f:
            f.write(writeme)
    for mass in masses:
        print "source calculate_%i.csh" % mass

elif step == "REPORT":
    import re
    for mass in masses:
        print str(mass)+":"
        if not unblind:
            alog = "aa3.log"
            plog = "pp3.log"
        else:
            alog = "AA3.log"
            plog = "PP3.log"
        with open("%i/%s" % (mass,alog)) as f:
            for line in f:
                m = re.match("Observed Limit: r < (.*)", line)
                if m:  print m.group(1)
                m = re.match("Expected  2.5%: r < (.*)", line)
                if m:  print m.group(1)
                m = re.match("Expected 16.0%: r < (.*)", line)
                if m:  print m.group(1)
                m = re.match("Expected 50.0%: r < (.*)", line)
                if m:  print m.group(1)
                m = re.match("Expected 84.0%: r < (.*)", line)
                if m:  print m.group(1)
                m = re.match("Expected 97.5%: r < (.*)", line)
                if m:  print m.group(1)
        print
        with open("%i/%s" % (mass,plog)) as f:
            for line in f:
                m = re.match("p-value of background: (.*)", line)
                if m:  print m.group(1)
                m = re.match("       \(Significance = (.*)\)", line)
                if m:  print m.group(1)
        print

        if not unblind:
            alog = "aa1.log"
            plog = "pp1.log"
        else:
            alog = "AA1.log"
            plog = "PP1.log"
        with open("%i/%s" % (mass,alog)) as f:
            for line in f:
                m = re.match("Expected 50.0%: r < (.*)", line)
                if m:  print m.group(1)
        print
        with open("%i/%s" % (mass,plog)) as f:
            for line in f:
                m = re.match("p-value of background: (.*)", line)
                if m:  print m.group(1)
        print
    
        if not unblind:
            alog = "aa2.log"
            plog = "pp2.log"
        else:
            alog = "AA2.log"
            plog = "PP2.log"
        with open("%i/%s" % (mass,alog)) as f:
            for line in f:
                m = re.match("Expected 50.0%: r < (.*)", line)
                if m:  print m.group(1)
        print
        with open("%i/%s" % (mass,plog)) as f:
            for line in f:
                m = re.match("p-value of background: (.*)", line)
                if m:  print m.group(1)
        print

