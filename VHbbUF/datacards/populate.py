import os

ref = "./125_ref"

for mass in [125]:
#for mass in range(110, 150+5, 5):
    os.system("mkdir -p %i" % mass)
    os.system("cp %s/run_*csh %i" % (ref, mass))
    os.system("sed -i 's/125/%i/g' %i/run_combinecards.csh" % (mass, mass))

    writeme = "cd %i\n" % mass
    writeme += "source run_combinecards.csh\n"
    writeme += "combine -M Asymptotic -t -1 vhbb_VH_8TeV.txt >! aa5.log\n"
    writeme += "combine -M ProfileLikelihood -m %i --signif --pvalue -t -1 --toysFreq --expectSignal=1 vhbb_VH_8TeV.txt >! pp5.log\n" % mass
    if mass <= 135:
        writeme += "combine -M Asymptotic -t -1 vhbb_VH_7p8TeV.txt >! aa8.log\n"
        writeme += "combine -M ProfileLikelihood -m %i --signif --pvalue -t -1 --toysFreq --expectSignal=1 vhbb_VH_7p8TeV.txt >! pp8.log\n" % mass
    writeme += "cd -\n"
    with open("calculate_%i.csh" % mass, 'w') as f:
        f.write(writeme)

if False:
    import re
    alog = "aa5.log"
    plog = "pp5.log"
    for mass in [110,125]:
    #for mass in range(110, 150+5, 5):
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
    
