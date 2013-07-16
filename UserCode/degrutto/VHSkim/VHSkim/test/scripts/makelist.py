#! /usr/bin/env python


import sys
import os
import commands
import string


dir = raw_input("What Directory (eg oct9Ntuple)?\n")
os.system("source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.sh")
os.system("source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.sh")
os.system("voms-proxy-init -voms cms")
os.system("echo '\n\nGetting Filenames in Folder: %s... (may take a couple of minutes)'" % dir)
os.system("lcg-ls 'srm://cmsdcache.pi.infn.it:8443/srm/managerv2?SFN=/pnfs/pi.infn.it/data/cms/store/user/arizzi/%s' > %s.txt" % (dir,dir))
os.system("sed -i 's@/pnfs/pi.infn.it/data/cms/store/user/arizzi/%s/@@g' %s.txt" % (dir,dir))
os.system("sed -i 's@Test@\"Test@g' %s.txt" % dir) 
os.system("sed -i 's@DiJet@\"DiJet@g' %s.txt" % dir) 
os.system("sed -i 's@BestCSV@\"BestCSV@g' %s.txt" % dir) 
os.system("sed -i 's@root@root\",@g' %s.txt" % dir)
os.system("echo '#! /usr/bin/env python' > GET%s.py" % dir)
os.system("echo 'import sys' >> GET%s.py" % dir)
os.system("echo 'import os' >> GET%s.py" % dir)
os.system("echo 'import commands' >> GET%s.py" % dir)
os.system("echo 'import string' >> GET%s.py" % dir)
os.system("echo '\n\n' >> GET%s.py" % dir)
os.system("echo 'def fib(n,dir):' >> GET%s.py" % dir)
os.system("echo '    os.system(\"echo Copying File: %%s\" %% n)' >> GET%s.py" % dir)
os.system("echo '    os.system(\"lcg-cp -V cms --srm-timeout=6000 -n 1 --connect-timeout=6000 -T srmv2 -U srmv2 \'srm://cmsdcache.pi.infn.it:8443/srm/managerv2?SFN=/pnfs/pi.infn.it/data/cms/store/user/arizzi/%%s//%%s\' %%s\" %% (dir,n,n))' >> GET%s.py" % dir)
os.system("echo '    os.system(\"echo Finished File: %%s\" %% n)\n\n' >> GET%s.py" % dir)
os.system("echo 'os.system(\"voms-proxy-init -voms cms\")' >> GET%s.py" % dir)
os.system("echo 'dir = \"%s\" \n\n' >> GET%s.py" % (dir,dir))
os.system("echo 'L = [' >> GET%s.py" % dir)
os.system("more %s.txt >> GET%s.py" % (dir,dir))
os.system("echo ']\n\n' >> GET%s.py" % dir)
os.system("echo 'for item in L:' >> GET%s.py" % dir)
os.system("echo '   fib(item,dir)' >> GET%s.py" % dir)
os.system("chmod +x GET%s.py" % dir)
os.system("echo 'GET%s.py written'" % dir)
os.system("rm %s.txt" % dir)
