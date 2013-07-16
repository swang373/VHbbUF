import FWCore.ParameterSet.Config as cms

process = cms.Process("MuEleMetBbbarSkimv2")

# Setup PAT
from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.pfTools import *
from PhysicsTools.PatAlgos.tools.trigTools import *

process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load('Configuration.StandardSequences.Services_cff')

# Specify the Global Tag
process.GlobalTag.globaltag = 'START311_V2::All'

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.EventContent.EventContent_cff')


import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
#process.hltFilter.HLTPaths = ["HLT_MET100_v1"]
#process.hltFilter.HLTPaths = ["HLT_PFMHT150_v2"]
#process.hltFilter.HLTPaths = ["HLT_2CentralJet20_BtagIP_SingleTrackTC_MET65_v1"]
#process.hltFilter.HLTPaths = ["HLT_MET45_2CentralJet20_2BTagIP_1LooseTC_1Tight_v1"]
#process.hltFilter.HLTPaths = ["HLT_MET45_2CentralJet20_1BTagIP_SingleTrackTC_v1"]
#process.hltFilter.HLTPaths = ["HLT_MET45_2CentralJet20_2BTagIP_SingleTrackTC_v1"]
#process.hltFilter.HLTPaths = ["HLT_MET45_2CentralJet20_2BTagIP_1BTag6IP_SingleTrackTC_v1"]
#process.hltFilter.HLTPaths = ["HLT_MET45_2CentralJet20_1BTagIP3D_SingleTrackTC_v1"]
#process.hltFilter.HLTPaths = ["HLT_MET45_2CentralJet20_2BTagIP3D_SingleTrackTC_v1"]#
#process.hltFilter.HLTPaths = ["HLT_MET45_2CentralJet20_2BTagIP3D_1BTag6IP_SingleTrackTC_v1"]


#HLT_2CentralJet20_BtagIP_MET65_v1"]
#HLT_2CentralJet20_MET80_v1"]
#HLT_CentralJet55CentralJet35_MET80_v1"]
#process.hltFilter.HLTPaths = ["HLT_PFMHT150_v2"]
#process.hltFilter.HLTPaths =["HLT_CentralJet55CentralJet35_MET80_v1"]
#process.hltFilter.HLTPaths =["HLT_2CentralJet20_BtagIP_MET65_v1"]

#process.hltFilter.HLTPaths=["HLT_MET45_2CentralJet20_2BTagIP_SingleTrackTC_v1"]
# Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load('L1Trigger.Skimmer.l1Filter_cfi')
#process.l1Filter.algorithms = cms.vstring('L1_ETM30')
process.l1Filter.algorithms = cms.vstring('L1_ETM50')



# Source file : To be run on a Full RECO sample
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'file://data/4/madfish/step2_RAW2DIGI_L1Reco_RECOApr14.root'
#file:/uscms_data/d2/dlopes/CMSSW_4_2_0_pre8/src/GeneratorInterface/MCatNLOInterface/test/step2_RAW2DIGI_L1Reco_RECOApr14.root',
#file:/uscms_data/d2/dlopes/CMSSW_4_2_0_pre8/src/GeneratorInterface/MCatNLOInterface/test/step2_RAW2DIGI_L1Reco_RECOApr14_2.root',
#file:/uscms_data/d2/dlopes/CMSSW_4_2_0_pre8/src/GeneratorInterface/MCatNLOInterface/test/step2_RAW2DIGI_L1Reco_RECOApr14_3.root','
#"dcache:/pnfs/cms/WAX/11/store/user/dlopes/testZnunuHGen413_2/step2_RAW2DIGI_L1Reco_RECO.root",
#"file:../../../step2_RAW2DIGI_L1Reco_RECO_v2.root",
#"file:../../../step2_RAW2DIGI_L1Reco_RECO_v3.root",
#"file:../../../step2_RAW2DIGI_L1Reco_RECO_v4.root",
#"file:../../../step2_RAW2DIGI_L1Reco_RECO_v5.root"


    

    )
)


isData = False 

removeMCMatching(process, ['All'])

removeAllPATObjectsBut(process, ['METs', 'Jets'])



  
## now selecting a good MET!!!
# Switch to using PFMET 
switchToPFMET(
    process, 
    cms.InputTag('pfMet'), 
    ""
)



#process.selectedPatMET = cms.EDFilter("PATMETSelector",
#     src = cms.InputTag("patMETs"),
#     cut = cms.string("et > 0 && mEtSig>0 ")
# )



#process.METFilter = cms.EDFilter("CandViewCountFilter",
#                                     src = cms.InputTag("selectedPatMET"),
#                                     minNumber = cms.uint32(1)
#                                 )


switchOnTrigger( process )

process.metTriggerMatchHLTMET45 = cms.EDProducer(
      "PATTriggerMatcherDRLessByR"                    # match by DeltaR only, best match by DeltaR
      , src     = cms.InputTag( 'patMETs' )
      , matched = cms.InputTag( 'patTrigger' )          # default producer label as defined in PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
      #, matchedCuts = cms.string( 'type("TriggerMET")')
      , matchedCuts = cms.string( 'path( "HLT_MET45_2CentralJet20_2BTagIP_SingleTrackTC_v1", 0) && type("TriggerMET")')
      , maxDPtRel = cms.double( 3.0 )
      , maxDeltaR = cms.double( 0.5 )
      , resolveAmbiguities    = cms.bool( True )        # only one match per trigger object
      , resolveByMatchQuality = cms.bool( True )        # take best match found per reco object: by DeltaR here (s. above)
      )

process.selectedMETTriggerMatch = cms.EDProducer( "PATTriggerMatchMETEmbedder",
                                                   src     = cms.InputTag( "patMETs" ),
                                                   matches = cms.VInputTag( "metTriggerMatchHLTMET45" )

)

# Switch to using PFJets
switchJetCollection(process,cms.InputTag('ak5PFJets'),
    doJTA        = True,
    doBTagging   = True,
    jetCorrLabel = ('AK5PF', cms.vstring(['L2Relative', 'L3Absolute', 'L2L3Residual'])),
    doType1MET   = True,
 #   genJetCollection=cms.InputTag("ak5GenJets"),
    doJetID      = True
)

if isData!=True: 
    switchJetCollection(process,cms.InputTag('ak5PFJets'),
                            doJTA        = True,
                            doBTagging   = True,
                            jetCorrLabel = ('AK5PF', cms.vstring(['L2Relative', 'L3Absolute' ])),
                            doType1MET   = True,
                            genJetCollection=cms.InputTag("ak5GenJets"),
                            doJetID      = True
                        )
    


process.patJets.addTagInfos  = False# ### needed for secondary vertex variables....

### ntuple
process.load("METbb.VHSkim.VHEdmNtuples_cff")







# Clean the Jets from the selected leptons, applying here only pt cut pf 20 
process.cleanPatJets = cms.EDProducer("PATJetCleaner",
                                      src = cms.InputTag("selectedPatJets"),
                                     # preselection = cms.string('pt > 20.0 && abs(eta) < 2.5 && ( (neutralEmEnergy/energy < 0.99) &&  (neutralHadronEnergy/energy < 0.99) && numberOfDaughters>1) '),
                                         preselection = cms.string('pt > 10.0 & (neutralEmEnergy/energy < 0.99) &&  (neutralHadronEnergy/energy < 0.99) && numberOfDaughters>1'),
                                      checkOverlaps = cms.PSet(
####
    
    ),
                                      ###applying 2 loose cutson the b-tagged jets
                                      
                                      finalCut = cms.string('')
                                      ##finalCut = cms.string('bDiscriminator(\"trackCountingHighEffBJetTags\")>1.7 ')
                                      )

### select pat jets, applying eta loose ID and  b-tag
process.myPatJets = cms.EDFilter("PATJetSelector",
                               src= cms.InputTag("cleanPatJets"),
                               cut = cms.string("abs(eta)< 2.5 & pt > 20.0" ) #&& bDiscriminator(\"trackCountingHighEffBJetTags\")>1.7 " ) ### try without loose b-tag now..... 
                               )
                                      

process.jetFilter = cms.EDFilter("CandViewCountFilter",
                                     src = cms.InputTag("myPatJets"),
                                     minNumber = cms.uint32(2),

                                 )





process.selJetTriggerMatch = cms.EDProducer(
      "PATTriggerMatcherDRLessByR"                 # match by DeltaR only, best match by DeltaR
      , src     = cms.InputTag( "myPatJets" )
      , matched = cms.InputTag( "patTrigger" )          # default producer label as defined in PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
      , matchedCuts = cms.string( 'path("HLT_MET45_2CentralJet20_2BTagIP_SingleTrackTC_v1", 0) && type("TriggerJet")')
      , maxDPtRel = cms.double( 3.0 )
      , maxDeltaR = cms.double(0.4)
      , resolveAmbiguities    = cms.bool( True)        # only one match per trigger object
      , resolveByMatchQuality = cms.bool( True )        # take best match found per reco object: by DeltaR here (s. above)
      )

process.selectedJetTriggerMatch = cms.EDProducer( "PATTriggerMatchJetEmbedder",
                                                   src     = cms.InputTag( "myPatJets" ),
                                                   matches = cms.VInputTag( "selJetTriggerMatch" )
                                               )



process.hjj = cms.EDProducer("CandViewShallowCloneCombiner",
    checkCharge = cms.bool(False),
    cut = cms.string('mass>0'),
    decay = cms.string("selectedJetTriggerMatch selectedJetTriggerMatch")
)   



process.hjjmet = cms.EDProducer("CandViewCombiner",
                            checkCharge = cms.bool(False),
                            cut = cms.string(''),
                            decay = cms.string("hjj selectedMETTriggerMatch")
                            )


process.HLTInfoMetJetEdmNtuple= cms.EDProducer(
        "HLTInfoDumper",
            hjj = cms.InputTag("hjj"),
            met = cms.InputTag("selectedMETTriggerMatch"),
            hltPath = cms.string("HLT_MET45_2CentralJet20_2BTagIP_SingleTrackTC_v1")
            )





process.ZinvHbbPath = cms.Path(
#       process.hltFilter +
#        process.l1Filter +
#     process.patTrigger +
#           process.patTriggerSequence *
        
#     process.makePatJets +
#     process.makePatMETs +

    process.patDefaultSequence +



        
#     process.selectedPatMET +
     process.metTriggerMatchHLTMET45+
        process.selectedMETTriggerMatch+
#     process.METFilter +

     process.cleanPatJets +
     process.myPatJets +

     process.jetFilter +
     process.selJetTriggerMatch +
     process.selectedJetTriggerMatch +

        
     process.hjj *
     process.HjjEdmNtuple *
        process.HLTInfoMetJetEdmNtuple*
     process.MetEdmNtuple *
     process.hjjmet *
     process.HjjMetEdmNtuple 


     
 )



#switchOnTriggerMatching( process,  [ 'selJetTriggerMatch' ] )
 # Switch to selected PAT objects in the trigger matching

#removeCleaningFromTriggerMatching( process )

#switchOnTriggerMatchEmbedding( process, [ 'selectedJetTriggerMatch' ] )


process.MuEleMetBbbarSkimEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring(
    'keep *_selectedPatElectrons_*_*',
    'keep *_selectedPatMuons_*_*',
    'keep *_selectedPatMET_*_*',
    'keep *_cleanPatJets_*_*',
    'keep *_zee_*_*',
    'keep *_zmm_*_*',
    'keep *_hjj_*_*',
    'keep *_wmn_*_*',
    'keep *_wen_*_*',
    

    )
)

# to add the AOD output uncomment the following
#process.MuEleMetBbbarSkimEventContent.outputCommands.extend(process.AODEventContent.outputCommands)


process.out = cms.OutputModule("PoolOutputModule",
    process.MuEleMetBbbarSkimEventContent,
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('ZinvHbbPath')
        ),
                                                       
       dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('muelemetbbarSkim'),
        dataTier = cms.untracked.string('USER')
   ),
   fileName = cms.untracked.string('MuEleMetBbbarSkim.root')
)







process.edmNtuplesOut.fileName = cms.untracked.string('14AprnoHLTFilter_MET45_2CentralJet20_2BTagIP_SingleTrackTC_v1_matching.root')
process.edmNtuplesOut.SelectEvents = cms.untracked.PSet(
            SelectEvents = cms.vstring( 'ZinvHbbPath',
                                  
                                       )
            )                            

process.endPath = cms.EndPath(process.out +
                              process.edmNtuplesOut)


