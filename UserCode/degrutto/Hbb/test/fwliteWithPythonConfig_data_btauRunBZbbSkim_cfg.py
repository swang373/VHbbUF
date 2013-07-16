import FWCore.ParameterSet.Config as cms

process = cms.Process("FWLitePlots")


##process.fwliteParameters = cms.PSet(
process.fwliteInput = cms.PSet(
#  fileNames   = cms.vstring('rfio:/castor/cern.ch/user/d/degrutto/test/jetRunB/zjjEdmNtuples_90_1_OCA.root'),         ## mandatory

    fileNames   = cms.vstring(),
    maxEvents   = cms.int32(-1),                             ## optional
    outputEvery = cms.uint32(1000),                            ## optional
)

import os
#dirname = "/data/1a/degrutto/data/MuRunBZbbSkim/"
dirname = "/data/1a/degrutto/data/CorrJet/btauRunB/"
dirlist = os.listdir(dirname)
basenamelist = os.listdir(dirname + "/")
for basename in basenamelist:
             process.fwliteInput.fileNames.append("file:" + dirname + "/" + basename)
print "Number of files to process is %s" % (len(process.fwliteInput.fileNames))



#



process.fwliteOutput = cms.PSet(
    fileName  = cms.string('analyzeFWLiteHistograms_data_btauRunBZbbSkim.root'),## mandatory
    )

process.Analyzer = cms.PSet(
    ## using now tight cut on jbProb_cut
    ### no cut  on CSVMVA 
    jbProb_cut   = cms.double(1.95),
    jCSVMVA_cut  = cms.double(0.9),
    j1pt_cut = cms.double(30.),
    j2pt_cut = cms.double(20.),
    jjDeltaRmin_cut = cms.double(2.5),
    jjDeltaRmax_cut = cms.double(100.5),
    pfMet_cut = cms.double(30.),
    jjDistSVtx_cut =     cms.double(0.),
    )
