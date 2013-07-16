 
import FWCore.ParameterSet.Config as cms

process = cms.Process("FWLitePlots")

###fileName   = cms.string('file:/data/1a/degrutto/data/data11/METv2/sample.root'),         ## mandatory
process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(),
    PUmcfileName = cms.string(''),
    PUdatafileName = cms.string('PUdist_160404-167151_7TeV_PromptReco_Collisions11.root'),
    maxEvents   = cms.int32(0),                             ## optional
    outputEvery = cms.uint32(0),                            ## optional
    )


 
##channel =  "METBTagMay10"
#channel =  "Zinv100"
##channel =  "WJetsMad"
##channel =  "ZinvHbb115"
#T-t  T-tWDR  TTbar  Tbar-s 
channel =  "T-t"
import os
##data
##dirnameOld = "/exports/1a/degrutto/data/data11/42X/"
## mc
dirnameOld = "/exports/1a/degrutto/data/Summer11/"






#for i in range(len(channels)):
 

dirname =  dirnameOld + channel
process.fwliteInput.PUmcfileName = dirname + '/PUMC.root' 

dirlist = os.listdir(dirname)
basenamelist = os.listdir(dirname + "/")
for basename in basenamelist:
###for basename in basenamelist:
    if (basename[0:3]== "Edm" and basename[-4:]== "root") :
#        process.fwliteInput.fileName = "file:" + dirname + "/" + basename
         process.fwliteInput.fileNames.append("file:" + dirname + "/" + basename)
         #print "file to process is %s" %process.fwliteInput.fileNames
print "Number of files to process is %s" % (len(process.fwliteInput.fileNames))    
    



#


fname = 'ZinvH42X' + channel + '.root'

process.fwliteOutput = cms.PSet(
    
    fileName  = cms.string(fname),## mandatory
    )

process.Analyzer = cms.PSet(
   lumi       =     cms.double(10.),
       sampleName =     cms.string(channel),
       isData     =     cms.bool(False),
       isZinv     =     cms.bool(True),
       isZee      =     cms.bool(False),
       isZmm      =     cms.bool(False),
       isWln      =     cms.bool(False),
       doJID      =     cms.bool(True),
       isMC      =     cms.bool(True)
   
    )

    
  
    
