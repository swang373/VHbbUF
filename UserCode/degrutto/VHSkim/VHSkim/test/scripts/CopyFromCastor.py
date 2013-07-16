#! /usr/bin/env python

import sys
import os
import commands
import string


castorFolder = raw_input('Enter Full Source Castor Path (ie /castor/cern.ch/user/m/madfish/Data2011/): ')
destFolder = raw_input('Enter Full Destination Path (will be created if not present): ')
print('To Begin Copying, Execute The File %s/temp' % destFolder) 

os.system('mkdir -p %s' % destFolder)
os.system('nsls %s | grep Edm > %s/temp ' % (castorFolder,destFolder))
os.system('chmod +x %s/temp ' % destFolder)
os.system('sed -i \'s@Edm@rfcp %s/Edm@g\' %s/temp' % (castorFolder,destFolder))
os.system('sed -i \'s@root@root .@g\' %s/temp' % destFolder)
os.system('more  %s/temp | wc -l' % destFolder)

os.system('echo "hadd sample.root *.root" >> %s/temp' % destFolder)
os.system('echo "rm Edm*.root" >> %s/temp' % destFolder)
