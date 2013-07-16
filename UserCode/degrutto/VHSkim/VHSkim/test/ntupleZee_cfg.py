import FWCore.ParameterSet.Config as cms

process = cms.Process("FWLitePlots")

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(),
    maxEvents   = cms.int32(-1),                             ## optional
    outputEvery = cms.uint32(0),                            ## optional
    )

#ZllHbb115

channel =  "TTbar"

import os
#dirnameOld = "/data/4/madfish/VHSkim/"
dirnameOld = "/data/1a/degrutto/data/Spring11/ZllHbb/"
#for i in range(len(channels)):
 

dirname =  dirnameOld + channel
dirlist = os.listdir(dirname)
basenamelist = os.listdir(dirname + "/")
for basename in basenamelist:
    process.fwliteInput.fileNames.append("file:" + dirname + "/" + basename)
print "Number of files to process is %s" %(len(process.fwliteInput.fileNames)) 

fname = 'ZeeH' + channel + '.root'

process.fwliteOutput = cms.PSet(
    
    fileName  = cms.string(fname),## mandatory
    )

process.Analyzer = cms.PSet(
    lumi       =     cms.double(10.),
    sampleName =     cms.string(channel),
    isData     =     cms.bool(True),
    isZinv     =     cms.bool(False),
    isZee      =     cms.bool(True),
    isZmm      =     cms.bool(False),
    isWln      =     cms.bool(False),
    )

