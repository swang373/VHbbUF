import FWCore.ParameterSet.Config as cms

process = cms.Process("FWLitePlots")


##process.fwliteParameters = cms.PSet(
process.fwliteInput = cms.PSet(
  fileNames   = cms.vstring('file:2l2bMetEdmNtuples.root'),         ## mandatory

#    fileNames   = cms.vstring(),
    maxEvents   = cms.int32(-1),                             ## optional
    outputEvery = cms.uint32(1000),                            ## optional
)

import os
#dirname = "/data/1a/degrutto/data/MuRunBZbbSkim/"
#dirname = "//btauRunB/"
#dirlist = os.listdir(dirname)
#basenamelist = os.listdir(dirname + "/")
#for basename in basenamelist:
#             process.fwliteInput.fileNames.append("file:" + dirname + "/" + basename)
#print "Number of files to process is %s" % (len(process.fwliteInput.fileNames))



#



process.fwliteOutput = cms.PSet(
    fileName  = cms.string('analyzeFWLiteHistograms_data_jet20_btauRunBZbbSkim.root'),## mandatory
    )

process.Analyzer = cms.PSet(
      ### no cut  on CSVMVA and jbProb
    jbProb_cut   = cms.double(1.95),
 ##   jCSVMVA_cut  = cms.double(0.9),
##    jbProb_cut   = cms.double(0.99),
    jTCHP_cut   = cms.double(3.14),
    jCSVMVA_cut  = cms.double(-10000),
    j1pt_cut = cms.double(30.),
    j2pt_cut = cms.double(30.),
    jjDeltaRmin_cut = cms.double(2.5),
    jjDeltaPhimin_cut = cms.double(2.5),
    jjDeltaRmax_cut = cms.double(100.5),
    pfMet_cut = cms.double(9999.),
    jjDistSVtx_cut =     cms.double(-99999),
    )
