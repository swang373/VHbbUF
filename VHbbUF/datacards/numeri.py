import subprocess
import os
import re
import glob

header1 = "HelperBDTShape.h"
macro1 = "BDTShapeJ12.C"
macro2 = "BDTShapeWorkspaceJ12.C"
#combinecards = "run_combinecards.csh"
#combineasymp = "run_combineasymptoticJ12.csh"
#combineprofi = "run_combineprofilelikelihoodJ12.csh"
#combinemaxli = "run_combinemaxlikelihoodJ12.csh"
combinecards = "run_combinecards_bbb.csh"
combineasymp = "run_combineasymptoticJ12_bbb.csh"
combineprofi = "run_combineprofilelikelihoodJ12_bbb.csh"
combinemaxli = "run_combinemaxlikelihoodJ12_bbb.csh"
observeasymp = "run_observeasymptoticJ12_bbb.csh"
observeprofi = "run_observeprofilelikelihoodJ12_bbb.csh"
observemaxli = "run_observemaxlikelihoodJ12_bbb.csh"
injectsignal = "run_injectsignalJ12_bbb.csh"

outdir = "res_20140407/"
afsdir = "/afs/cern.ch/user/j/jiafulow/public/vhbb_LHCP_20140407/Znn/"

#------------------------------------------------------------------------------

if not outdir.endswith("/"):
    outdir += "/"
header1 = header1.replace(".h", "")
macro1 = macro1.replace(".C", "")
macro2 = macro2.replace(".C", "")
newheader1 = header1+"_new"
newmacro1 = macro1+"_new"
newmacro2 = macro2+"_new"
arg = ""
argws = ""


def make(var, analysis, data):

    global arg, argws
    outsubdir = outdir+var+"_"+data+"/"+str(analysis)+"/"
    outsubdir = outsubdir.replace("_LHCP", "")
    if not os.path.exists(outsubdir):
        os.makedirs(outsubdir)

    subprocess.call(["cp", header1+".h", newheader1+".h"])
    subprocess.call(["cp", macro1+".C", newmacro1+".C"])
    subprocess.call(["cp", macro2+".C", newmacro2+".C"])
    subprocess.call(["sed", "-i", "s@%s@%s@" % (macro1+"(", newmacro1+"("), newmacro1+".C"])
    subprocess.call(["sed", "-i", "s@%s@%s@" % (macro2+"(", newmacro2+"("), newmacro2+".C"])
    subprocess.call(["sed", "-i", "s@%s@%s@" % (header1, newheader1), newmacro1+".C"])
    subprocess.call(["sed", "-i", "s@%s@%s@" % (header1, newheader1), newmacro2+".C"])

    if var == "MJJ":
        old = "//#define MJJANALYSIS"
        new = "#define MJJANALYSIS"
        subprocess.call(["sed", "-i", "s@%s@%s@" % (old, new), newmacro1+".C"])
        subprocess.call(["sed", "-i", "s@%s@%s@" % (old, new), newmacro2+".C"])
        old = "Step4_20130404/reload_20130401/"
        new = "Step4_20130404/stitch/"
        subprocess.call(["sed", "-i", "s@%s@%s@" % (old, new), newmacro1+".C"])

    elif var == "mBDT":
        arg = "500,100"

    elif var == "nBDT":
        arg = "500,-100"

        old = "//#define REBINFIRSTTHREE"
        new = "#define REBINFIRSTTHREE"
        subprocess.call(["sed", "-i", "s@%s@%s@" % (old, new), newheader1+".h"])

    elif var == "BDT":
        arg = ""


    if analysis == "VV":
        old = "//#define VVANALYSIS"
        new = "#define VVANALYSIS"
        subprocess.call(["sed", "-i", "s@%s@%s@" % (old, new), newmacro1+".C"])

        old = "Step4_20130404/reload_20130401/"
        new = "Step4_20130404/stitch/"
        subprocess.call(["sed", "-i", "s@%s@%s@" % (old, new), newmacro1+".C"])

    elif analysis in range(110, 150+5, 5):
        massH = analysis
        #massH = min(135, analysis)
        old = "const int massH             = 125;"
        new = "const int massH             = %i;" % massH
        subprocess.call(["sed", "-i", "s@%s@%s@" % (old, new), newmacro1+".C"])

        old = "Form(\"ZH%i\", massH)"
        new = "Form(\"ZH%i\", "+ ("%i)" % analysis)
        subprocess.call(["sed", "-i", "s@%s@%s@" % (old, new), newmacro1+".C"])

        old = "Form(\"WH%i\", massH)"
        new = "Form(\"WH%i\", "+ ("%i)" % analysis)
        subprocess.call(["sed", "-i", "s@%s@%s@" % (old, new), newmacro1+".C"])


    if data == "HCP":
        old = "//#define HCPANALYSIS"
        new = "#define HCPANALYSIS"
        subprocess.call(["sed", "-i", "s@%s@%s@" % (old, new), newmacro1+".C"])

    # Run
    p1 = subprocess.Popen(["root", "-b", "-l", "-q", "%s+(%s)" % (newmacro1+".C", arg)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = p1.communicate()
    #print output[1]

    # Signal injection
    p2 = subprocess.Popen(["root", "-b", "-l", "-q", "%s+(%s)" % (newmacro2+".C", "\"injectsignal\"")], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = p2.communicate()
    #print output[1]
    subprocess.call("cp vhbb_Znn_SI_J12_bbb_Znunu*Pt_8TeV.txt vhbb_Znn_SI_J12_8TeV.root " + outsubdir, shell=True)  # use shell for wildcard expansion

    # Default
    p3 = subprocess.Popen(["root", "-b", "-l", "-q", "%s+(%s)" % (newmacro2+".C", "")], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = p3.communicate()
    #print output[1]
    subprocess.call("cp vhbb_Znn_J12_bbb_Znunu*Pt_8TeV.txt vhbb_Znn_J12_8TeV.root vhbb_Znn_J12_Znunu*Pt_TH1.* diffNuisances.py " + outsubdir, shell=True)  # use shell for wildcard expansion
    subprocess.call(("cp %s " + outsubdir) % (" ".join([combinecards, combineasymp, combineprofi, combinemaxli, observeasymp, observeprofi, observemaxli, injectsignal])), shell=True)

    # Keep .pdf files
    pdfs = glob.glob("plotsJ12/*pdf")
    pdfs = sorted(pdfs, key=os.path.getmtime)
    for pdf in pdfs[-3:]:
        subprocess.call(["cp", pdf, outsubdir])


def numeri(var, analysis, data):

    outsubdir = outdir+var+"_"+data+"/"+str(analysis)+"/"
    outsubdir = outsubdir.replace("_LHCP", "")

    cwd = os.getcwd()
    os.chdir(outsubdir)
    p1 = subprocess.Popen(["/bin/tcsh", combinecards], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = p1.communicate()
    #print output[0]
    p2 = subprocess.Popen(["/bin/tcsh", combineasymp], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = p2.communicate()
    #print output[0]
    p3 = subprocess.Popen(["/bin/tcsh", combineprofi], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = p3.communicate()
    #print output[0]
    p4 = subprocess.Popen(["/bin/tcsh", combinemaxli], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = p4.communicate()
    #print output[0]
    p5 = subprocess.Popen(["/bin/tcsh", observeasymp], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = p5.communicate()
    #print output[0]
    p6 = subprocess.Popen(["/bin/tcsh", observeprofi], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = p6.communicate()
    #print output[0]
    p7 = subprocess.Popen(["/bin/tcsh", observemaxli], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = p7.communicate()
    #print output[0]
    p8 = subprocess.Popen(["/bin/tcsh", injectsignal], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = p8.communicate()
    #print output[0]

    # Make copies
    subprocess.call("cp vhbb_Znn_J12_bbb_combo_8TeV.txt vhbb_Znn_J12_combo_8TeV.txt", shell=True)
    subprocess.call("cp vhbb_Znn_SI_J12_bbb_combo_8TeV.txt vhbb_Znn_SI_J12_combo_8TeV.txt", shell=True)
    subprocess.call("cp A5.log asymp.log", shell=True)
    subprocess.call("cat P5.log p5.log > profi.log", shell=True)
    subprocess.call("cp M5.log maxli.log", shell=True)
    subprocess.call("cp N1.html nuisances_ZnunuHighPt.html", shell=True)
    subprocess.call("cp N6.html nuisances_combo2.html", shell=True)
    subprocess.call("cp N5.html nuisances_combo.html", shell=True)
    subprocess.call("cp SI1.log inject.log", shell=True)

    # Clean up
    subprocess.call("rm roostats-*", shell=True)

    print outsubdir
    with open("asymp.log") as f:
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
    with open("profi.log") as f:
        for line in f:
            m = re.match("p-value of background: (.*)", line)
            if m:  print m.group(1)
            m = re.match("       \(Significance = (.*)\)", line)
            if m:  print m.group(1)
    print
    with open("maxli.log") as f:
        for line in f:
            m = re.match("Best fit r: (.*)  (.*)/(.*)  \(68% CL\)", line)
            if m:
                print m.group(1)
                print m.group(2)
                print m.group(3)
    print

    os.chdir(cwd)


def transfer(outdir, afsdir):

    writeme_ = []
    for dirname, dirnames, filenames in os.walk(outdir):
    # make all subdirectories first.
        for subdirname in dirnames:
            path1 = os.path.join(dirname, subdirname)
            path2 = os.path.join(dirname.replace(outdir, afsdir), subdirname)
            if path1.count('/') != 2:  continue
            writeme_.append("mkdir -p " + path2)
            writeme_.append("cd " + path1)
            writeme_.append("cp vhbb_Znn_J12_bbb_Znunu*Pt_8TeV.txt vhbb_Znn_J12_bbb_combo_8TeV.txt vhbb_Znn_J12_8TeV.root vhbb_Znn_J12_Znunu*Pt_TH1.* vhbb_Znn_J12_combo_8TeV.txt asymp.log profi.log maxli.log nuisances*.html " + path2)
            writeme_.append("cp vhbb_Znn_SI_J12_* " + path2)
            writeme_.append("cd -")
    writeme = "\n".join(writeme_)
    with open("afstransfer_"+outdir.strip("/")+".csh", "w") as f:
        f.write(writeme)
    print "To transfer, do: \nsource afstransfer_"+outdir.strip("/")+".csh"


if __name__ == "__main__":

    #var = "BDT"
    var = "mBDT"
    #var = "nBDT"
    #var = "MJJ"
    data = "LHCP"
    #data = "HCP"
    for analysis in [125]:
    #for analysis in ["VV"] + range(110, 150+5, 5):
        make(var, analysis, data)
        numeri(var, analysis, data)

    transfer(outdir, afsdir)

