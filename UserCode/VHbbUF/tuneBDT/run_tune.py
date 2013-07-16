import os

# This specifies the set of NTrees parameters
#r = range(400,700,100)
r = [400, 500, 600, 800]

# This specifies the number of optimization points for a given NTrees.
# It should agree with the number in TrainBDT.C
n = 10

# Name of the weight folders
w = "weights_V6FullB"

#------------------------------------------------------------------------------

ww = w.replace("weights","")
ww = ww.lower()

print "# TRAIN"
for i in r:
    print r'root -l -b -q TrainBDT.C+\(125,\"\",\"\",\"\",\"%sNTrees%i\",1,%i\) >&! train%s_%i.log &' % (w,i,i,ww,i)

print "# RELOAD"
ss = []
for i in r:
    ss.append(r'root -b -l -q TrimTree.C+\(\"Step4_ZH125\",\"BDT:%sNTrees%i:%i\"\)' % (w,i,n))
    ss.append(r'root -b -l -q TrimTree.C+\(\"Step4_WH125\",\"BDT:%sNTrees%i:%i\"\)' % (w,i,n))
    ss.append(r'root -b -l -q TrimTree.C+\(\"Step4_Zj\",\"BDT:%sNTrees%i:%i\"\)' % (w,i,n))
    ss.append(r'root -b -l -q TrimTree.C+\(\"Step4_ZjHW\",\"BDT:%sNTrees%i:%i\"\)' % (w,i,n))
    ss.append(r'root -b -l -q TrimTree.C+\(\"Step4_Wj\",\"BDT:%sNTrees%i:%i\"\)' % (w,i,n))
    ss.append(r'root -b -l -q TrimTree.C+\(\"Step4_WjHW\",\"BDT:%sNTrees%i:%i\"\)' % (w,i,n))
    ss.append(r'root -b -l -q TrimTree.C+\(\"Step4_data_obs\",\"BDT:%sNTrees%i:%i\"\)' % (w,i,n))
    ss.append(r'root -b -l -q TrimTree.C+\(\"Step4_TT\",\"BDT:%sNTrees%i:%i\"\)' % (w,i,n))
    ss.append(r'root -b -l -q TrimTree.C+\(\"Step4_TTPowheg\",\"BDT:%sNTrees%i:%i\"\)' % (w,i,n))
    ss.append(r'root -b -l -q TrimTree.C+\(\"Step4_s_Top\",\"BDT:%sNTrees%i:%i\"\)' % (w,i,n))
    ss.append(r'root -b -l -q TrimTree.C+\(\"Step4_ZZ\",\"BDT:%sNTrees%i:%i\"\)' % (w,i,n))
    ss.append(r'root -b -l -q TrimTree.C+\(\"Step4_VV\",\"BDT:%sNTrees%i:%i\"\)' % (w,i,n))
    ss.append(r'root -b -l -q TrimTree.C+\(\"Step4_QCD\",\"BDT:%sNTrees%i:%i\"\)' % (w,i,n))
    s = '\n'.join(ss)
    with open('run_reload%s_%i.sh' % (ww,i), 'w') as f:
        f.write(s)
    del ss[:]
    print r'source run_reload%s_%i.sh >&! reload%s_%i.log &' % (ww,i,ww,i)
