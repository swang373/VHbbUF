#! /usr/bin/env python

import sys
import os
import commands
import string

numFilesPerJob = 100

def fib(n):
	os.system('echo "cd /afs/cern.ch/user/m/madfish/scratch0/batch/CMSSW_4_1_5" > Run%d.sh' % n)
	os.system('echo "source STARTUP" >> Run%d.sh' % n)
	os.system('echo "cd -" >> Run%d.sh' % n)
	os.system('echo "cmsRun /afs/cern.ch/user/m/madfish/scratch0/batch/CMSSW_4_1_5/src/MyCode/SigEventFilter/test/Job%d.py" >> Run%d.sh' % (n,n))
	os.system('echo "scp CSCTF%d.root madfish@pcufl2:/exports/2a/madfish/CSCTF-BADFILES" >> Run%d.sh' % (n,n))
        os.system("sed 's/NNN/%d/g' SigFilter.py > Job%d.py" % (n,n))
	os.system('echo "bsub -q cmscaf1nd -J JOB_RUN%d < Run%d.sh" >> submit' % (n,n))


num_lines = sum(1 for line in open('files'))
count = 0
os.system('rm submit')
os.system('echo " " > submit')
fib(count)

for x in range(1,num_lines+1):
    os.system("head -n 1 files >> Job%d.py" % count)
    os.system("sed -i '1,+0d' files")
    if x%numFilesPerJob == 0:
       os.system("echo '))' >> Job%d.py" % count)
       count += 1
       fib(count)
            
       
os.system("echo '))' >> Job%d.py" % count)
