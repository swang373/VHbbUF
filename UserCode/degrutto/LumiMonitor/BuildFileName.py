import os, sys
usage = "usage: %s fileName line"   %         os.path.basename(sys.argv[0])

if len(sys.argv) < 2:
      print usage
      sys.exit(2)
else:
      argv = sys.argv
      infile = argv[1]
      line = argv[2]
#      print argv[1]
      f = open(infile, 'read')
               


##f = open('runlog.log', 'r')

#for in in f.readlines():
#   run= i.rstrip()
#   fileloc = 'https://cmsweb.cern.ch/dqm/offline/data/browse/ROOT/OfflineData/Run2010/Mu/000'+ str(run[0:3]) + 'xx/DQM_V0001_R000' + str(run) + '__Mu__Run2010A-PromptReco-v1__RECO.root'
#   print fileloc

filelocv1 = 'x'
filelocv2 = 'x'
filelocv4 = 'x'
filenamev1 = 'x'
filenamev2 = 'x'
filenamev4 = 'x'

run= f.readlines()[int(line)].rstrip()

if( int(run) < 135800):
       filelocv1 = 'https://cmsweb.cern.ch/dqm/offline/data/browse/ROOT/OfflineData/Commissioning10/MinimumBias/000'+ str(run[0:4]) + 'xx/DQM_V0001_R000' + str(run) + '__MinimumBias__Commissioning10-DQM-May27thReReco_v1__DQM.root'
       filenamev1 = '/tmp/degrutto/DQM_V0001_R000' + str(run) + '__MinimumBias__Commissioning10-DQM-May27thReReco_v1__DQM.root'

if (135800< int(run)<136000):
   filelocv1 = 'https://cmsweb.cern.ch/dqm/offline/data/browse/ROOT/OfflineData/Run2010/Mu/000'+ str(run[0:3]) + 'xxx/DQM_V0001_R000' + str(run) + '__Mu__Run2010A-PromptReco-v1__RECO.root'
   filenamev1 = '/tmp/degrutto/DQM_V0001_R000' + str(run) + '__Mu__Run2010A-PromptReco-v1__RECO.root'    

if (int(run)>136000):
      filelocv1 = 'https://cmsweb.cern.ch/dqm/offline/data/browse/ROOT/OfflineData/Run2010/Mu/000'+ str(run[0:4]) + 'xx/DQM_V0001_R000' + str(run) + '__Mu__Run2010A-PromptReco-v1__RECO.root'
      filenamev1 = '/tmp/degrutto/DQM_V0001_R000' + str(run) + '__Mu__Run2010A-PromptReco-v1__RECO.root'


      filelocv2 = 'https://cmsweb.cern.ch/dqm/offline/data/browse/ROOT/OfflineData/Run2010/Mu/000'+ str(run[0:4]) + 'xx/DQM_V0001_R000' + str(run) + '__Mu__Run2010A-PromptReco-v2__RECO.root'
      filenamev2 = '/tmp/degrutto/DQM_V0001_R000' + str(run) + '__Mu__Run2010A-PromptReco-v2__RECO.root'
#print filelocv1 , filenamev1, filelocv2,  filenamev2, run

      filelocv4 = 'https://cmsweb.cern.ch/dqm/offline/data/browse/ROOT/OfflineData/Run2010/Mu/000'+ str(run[0:4]) + 'xx/DQM_V0001_R000' + str(run) + '__Mu__Run2010A-PromptReco-v4__RECO.root'
      filenamev4 = '/tmp/degrutto/DQM_V0001_R000' + str(run) + '__Mu__Run2010A-PromptReco-v4__RECO.root'

print filelocv1 , filenamev1, filelocv2,  filenamev2, filelocv4,  filenamev4, run



 

#import os
#l= os.popen("rfdir /castor/cern.ch/cms/store/express/Run2010A/ExpressPhysics/FEVT/v2/000/136/088").readlines()
#for i in range(len(l)):
   #print l 
#   lstrip =l[i].split(' ')
#   file= lstrip[len(lstrip) - 1].rstrip()
   #### selectong the hour
#   if lstrip[len(lstrip) - 2][:3] == '18:' :  
#     print'"rfio:/castor/cern.ch/cms/store/express/Run2010A/ExpressPhysics/FEVT/v2/000/136/088/' + file + '",'

   
