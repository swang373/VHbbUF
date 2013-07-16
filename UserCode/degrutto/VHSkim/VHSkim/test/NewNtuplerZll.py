#! /usr/bin/env python


import sys
import os
import commands
import string

#make sure no ending / 
dirname = "/exports/1a/degrutto/data/Summer11"
channel = "XXX"
pufile = dirname + "/" + channel + "/PUMC.root"

numFile = 0
numFiles = 0

def process(n):
        global numFile
        numFile += 1
        print("Processing File (%i/%i): %s" % (numFile,numFiles,n)) 
        os.system('echo \"import FWCore.ParameterSet.Config as cms\" > ntupleZll-%s.py' % channel)
        os.system('echo \'process = cms.Process(\"FWLitePlots\")\' >> ntupleZll-%s.py' % channel)
        os.system('echo \"\" >> ntupleZll-%s.py' % channel)
        os.system('echo \"process.fwliteInput = cms.PSet(\" >> ntupleZll-%s.py' % channel)
        os.system('echo \"    fileName  = cms.string(\'%s/%s/%s\'),\" >> ntupleZll-%s.py' % (dirname, channel, n,channel))
        os.system('echo \'    PUmcfileName = cms.string(\"%s\"),\' >> ntupleZll-%s.py' % (pufile,channel))
        os.system("echo \'    PUdatafileName = cms.string(\"PUdist_160404-167151_7TeV_PromptReco_Collisions11.root\"),\' >> ntupleZll-%s.py" % channel)
        os.system('echo \"    maxEvents   = cms.int32(0),                             ## optional\" >> ntupleZll-%s.py' % channel)
        os.system('echo \"    outputEvery = cms.uint32(0),                            ## optional\" >> ntupleZll-%s.py' % channel)
        os.system('echo \"    )\" >> ntupleZll-%s.py' % channel)
        os.system('echo \"\" >> ntupleZll-%s.py' % channel)
        os.system('echo \"import os\" >> ntupleZll-%s.py' % channel)
        os.system('echo \"\" >> ntupleZll-%s.py' % channel)
        os.system('echo \"\" >> ntupleZll-%s.py' % channel)
        os.system("echo \"fname = \'%s-%s\' \" >> ntupleZll-%s.py" % (channel,n,channel))
        os.system('echo \"\" >> ntupleZll-%s.py' % channel)
        os.system('echo \"process.fwliteOutput = cms.PSet(\" >> ntupleZll-%s.py' % channel)
        os.system('echo \"    \" >> ntupleZll-%s.py' % channel)
        os.system('echo "    fileName  = cms.string(fname),## mandatory" >> ntupleZll-%s.py' % channel)
        os.system('echo "    )" >> ntupleZll-%s.py' % channel)
        os.system('echo "" >> ntupleZll-%s.py' % channel)
        os.system('echo "process.Analyzer = cms.PSet(" >> ntupleZll-%s.py' % channel)
        os.system('echo "    isData    =     cms.bool(False)," >> ntupleZll-%s.py' % channel)
        os.system('echo "    isZinv    =     cms.bool(False)," >> ntupleZll-%s.py' % channel)
        os.system('echo "    isZll     =     cms.bool(True)," >> ntupleZll-%s.py' % channel)
        os.system('echo "    isWln     =     cms.bool(False)," >> ntupleZll-%s.py' % channel)
        os.system('echo "    doJID     =     cms.bool(False)," >> ntupleZll-%s.py' % channel)
        os.system('echo "    isMC      =     cms.bool(True)" >> ntupleZll-%s.py' % channel)
        os.system('echo "    )" >> ntupleZll-%s.py' % channel)
        os.system('echo "" >> ntupleZll-%s.py' % channel)
        os.system('Ntupler ntupleZll-%s.py' % channel)
        os.system('rm ntupleZll-%s.py' % channel)



dirlist = os.listdir(dirname+"/"+channel)
numFiles = len(dirlist)


for filename in dirlist:
    if (filename[0:3]== "Edm" and filename[-4:]== "root") :
      process(filename)

os.system("hadd 42-%s.root  %s-EdmNtuples*.root" % (channel,channel))
