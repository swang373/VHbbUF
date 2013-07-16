
import FWCore.ParameterSet.Config as cms

process = cms.Process("FWLitePlots")

#fileNames   = cms.vstring('file:2l2bMetEdmNtuples.root'),         ## mandatory
process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(),
    maxEvents   = cms.int32(-1),                             ## optional
    outputEvery = cms.uint32(0),                            ## optional
    )


channel =  "TTbar"
import os
dirnameOld = "/data/1a/madfish/VHSkim/"







#for i in range(len(channels)):
 

dirname =  dirnameOld + channel
dirlist = os.listdir(dirname)
basenamelist = os.listdir(dirname + "/")
for basename in basenamelist:
    process.fwliteInput.fileNames.append("file:" + dirname + "/" + basename)
print "Number of files to process is %s" %(len(process.fwliteInput.fileNames)) 
    
    



#


fname = 'WH' + channel + '.root'

process.fwliteOutput = cms.PSet(
    
    fileName  = cms.string(fname),## mandatory
    )

process.Analyzer = cms.PSet(
    #### bbbar piece
#    jbProb_cut   = cms.double(1.9),
    jbProb_cut   = cms.double(0.0),
    jTCHP_cut   = cms.double(0.0),
    jCSVL_cut  = cms.double(0.6),
    jCSVH_cut  = cms.double(0.6),
    jCSVMVA_cut  = cms.double(0.0),
    j1pt_cut = cms.double(20.),
    j2pt_cut = cms.double(20.),
    jjDeltaRmin_cut = cms.double(0.0),
    jjDeltaPhimin_cut = cms.double(-99999),
    zjjmetDeltaPhimin_cut = cms.double(-99999),
    jjDeltaRmax_cut = cms.double(1.5),
    pfMet_cut = cms.double(0),
    pfMetSig_cut = cms.double(0),
    jjDistSVtx_cut =     cms.double(-99999),
    zjjpt_cut =     cms.double(150),

    
    isZinv =     cms.bool(False),
    isWln  =     cms.bool(True),
    isZll =     cms.bool(False),
    thirdJetPt_cut =     cms.double(20),
    
    nZjj_cut =     cms.uint32(9999),
    
    nofLeptons_cut =     cms.uint32(999),
    

    ###anticuts for suppressing ttbar
    ### Zmm piece
    
    
    zmmmassMin_cut = cms.double(0),
    zmmmassMax_cut = cms.double(1000),
    zmmptMin_cut = cms.double(150),
    zjjzmmDeltaPhimin_cut = cms.double(-9999),


    ### Zee piece
    zeemassMin_cut = cms.double(0),
    zeemassMax_cut = cms.double(1000),
    zeeptMin_cut = cms.double(150),
    zjjzeeDeltaPhimin_cut = cms.double(-99999),


    ### wmn piece
    wmnmtMin_cut = cms.double(0),
    wmnptMin_cut = cms.double(150),
    zjjwmnDeltaPhimin_cut = cms.double(3.0),

    ### wen piece
    wenmtMin_cut = cms.double(0),
    wenptMin_cut = cms.double(150),
    zjjwenDeltaPhimin_cut = cms.double(3.0)


    )

    
  
    
