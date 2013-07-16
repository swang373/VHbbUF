import os

# This specifies the set of NTrees parameters
#r = range(400,700,100)
r = [400,500,600,800]

# This specifies the number of optimization points for a given NTrees.
# It should agree with the number in TrainBDT.C
n = 10

# Name of the weight folders
w = "weights_V6FullB"
w = w.replace("weights","skim")

# Name of the log folder
logs = "logs_V6FullB_20130329"

if not os.path.exists(logs):
    os.makedirs(logs)

from ROOT import gROOT
flag = gROOT.ProcessLine('.L BDTShapeJ12_reload.C++O')

for i in r:
    for j in xrange(n):
        os.system(r'root -b -l -q BDTShapeJ12_reload.C+\(\"skim/%sNTrees%i\",%i\) >& %s/a1_%i_%i' % (w,i,j,logs,i,j))
        os.system(r'root -b -l -q BDTShapeWorkspaceJ12.C+\(\) >& %s/a2_%i_%i' % (logs,i,j))
        os.system(r'combine -M Asymptotic -t -1 vhbb_Znn_J12_ZnunuHighPt_8TeV.txt >& %s/a3_%i_%i' % (logs,i,j))
        os.system(r'combine -M Asymptotic -t -1 vhbb_Znn_J12_ZnunuMedPt_8TeV.txt >& %s/a4_%i_%i' % (logs,i,j))
        os.system(r'combine -M Asymptotic -t -1 vhbb_Znn_J12_ZnunuLowPt_8TeV.txt  >& %s/a5_%i_%i' % (logs,i,j))
        #os.system(r'combine -M Asymptotic -t -1 vhbb_Znn_J12_ZnunuLowCSV_8TeV.txt >& %s/a6_%i_%i' % (logs,i,j))
        
        os.system(r'combine -m 125 -M ProfileLikelihood --signif --pvalue -t -1 --toysFreq --expectSignal=1 vhbb_Znn_J12_ZnunuHighPt_8TeV.txt >& %s/b3_%i_%i' % (logs,i,j))
        os.system(r'combine -m 125 -M ProfileLikelihood --signif --pvalue -t -1 --toysFreq --expectSignal=1 vhbb_Znn_J12_ZnunuMedPt_8TeV.txt  >& %s/b4_%i_%i' % (logs,i,j))
        os.system(r'combine -m 125 -M ProfileLikelihood --signif --pvalue -t -1 --toysFreq --expectSignal=1 vhbb_Znn_J12_ZnunuLowPt_8TeV.txt  >& %s/b5_%i_%i' % (logs,i,j))
        #os.system(r'combine -m 125 -M ProfileLikelihood --signif --pvalue -t -1 --toysFreq --expectSignal=1 vhbb_Znn_J12_ZnunuLowCSV_8TeV.txt >& %s/b6_%i_%i' % (logs,i,j))
