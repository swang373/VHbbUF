import os
l= os.popen("rfdir /castor/cern.ch/cms/store/express/Run2010A/ExpressPhysics/FEVT/v2/000/136/088").readlines()
for i in range(len(l)):
   #print l 
   lstrip =l[i].split(' ')
   file= lstrip[len(lstrip) - 1].rstrip()
   #### selectong the hour
   if lstrip[len(lstrip) - 2][:3] == '18:' :  
     print'"rfio:/castor/cern.ch/cms/store/express/Run2010A/ExpressPhysics/FEVT/v2/000/136/088/' + file + '",'

   
