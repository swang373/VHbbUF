import FWCore.ParameterSet.Config as cms

process = cms.Process("FWLitePlots")

process.fwliteInput = cms.PSet(
    fileName  = cms.string("DDD/CCC/FFF"),
    PUmcfileName = cms.string('PPP'),
    PUdatafileName = cms.string('PUdist_160404-167151_7TeV_PromptReco_Collisions11.root'),
    maxEvents   = cms.int32(0),                             ## optional
    outputEvery = cms.uint32(0),                            ## optional
    )

import os


fname = 'm-CCC-FFF' 

process.fwliteOutput = cms.PSet(
    
    fileName  = cms.string(fname),## mandatory
    )

process.Analyzer = cms.PSet(
    isData    =     cms.bool(False),
    isZinv    =     cms.bool(False),
    isZll     =     cms.bool(True),
    isWln     =     cms.bool(False),
    doJID     =     cms.bool(False),
    isMC      =     cms.bool(True)
    )

