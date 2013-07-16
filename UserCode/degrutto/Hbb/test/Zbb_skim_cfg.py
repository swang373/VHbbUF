
import FWCore.ParameterSet.Config as cms

process = cms.Process("ZbbSkim")


# Setup PAT


from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.pfTools import *


# Load some generic cffs  
process.load('Configuration.StandardSequences.Services_cff')
#process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
#process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# Specify the Global Tag
process.GlobalTag.globaltag = 'START38_V13::All'

# Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Source file : To be run on a Full RECO sample
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
"file:/data/1a/degrutto/data/Fall10/DYToBB_M_50/DYToBB_M_50_TuneZ2_1.root"
# "rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_8_6/RelValTTbar/GEN-SIM-RECO/START38_V13-v1/0068/98EA8C65-25E8-DF11-B0A0-0018F3D095F8.root"
 
    )
)






isData = True
isHLT6e31 = True
isHLT2e32 = False



if isData!=True: 

    process.load("PhysicsTools.HepMCCandAlgos.flavorHistoryProducer_cfi")
    process.load("PhysicsTools.HepMCCandAlgos.flavorHistoryFilter_cfi")

removeMCMatching(process, ['All'])




import FWCore.ParameterSet.Config as cms

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
if isHLT6e31:
     process.hltFilter.HLTPaths = ["HLT_BTagMu_Jet20U"]
     #     process.hltFilter.HLTPaths = ["HLT_Mu14"]
     process.hltFilter.throw=cms.bool(True)



if isHLT2e32:
    process.hltFilter.HLTPaths = ["HLT_BTagMu_DiJet30U_Mu5_v3","HLT_BTagMu_DiJet30U_v3", "HLT_BTagMu_DiJet20U_v3" ]     

if  isData!=True:
    process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults", "", "REDIGI38XPU")




# Muon Selection
process.selectedPatMuons.cut = (
        "pt > 15 && isGlobalMuon && isTrackerMuon && globalTrack().normalizedChi2 < 10 &&" +
            "innerTrack().hitPattern().numberOfValidTrackerHits > 10 && "                      +
            "innerTrack().hitPattern().numberOfValidPixelHits > 0 && "                         +
            "globalTrack().hitPattern().numberOfValidMuonHits > 0 && "                         +
            "dB < 0.2 && "                                                                     +
            "trackIso + caloIso < 0.15 * pt && "                                               +
            "numberOfMatches > 1 && abs(eta) < 2.4"
        )

process.muonFilter = cms.EDFilter("CandViewCountFilter",
                                     src = cms.InputTag("selectedPatMuons"),
                                     minNumber = cms.uint32(1)
                                 )


# Switch to using PFMET 
switchToPFMET(
    process, 
    cms.InputTag('pfMet'), 
    ""
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
                            jetCorrLabel = ('AK5PF', cms.vstring(['L2Relative', 'L3Absolute', 'L2L3Residual'])),
                            doType1MET   = True,
                            genJetCollection=cms.InputTag("ak5GenJets"),
                            doJetID      = True
                        )
    

process.patJets.addTagInfos  = True # ### needed for secondary vertex variables....


# Clean the Jets from the seleted leptons, and apply loose(???) btag cuts
process.cleanPatJets = cms.EDProducer("PATJetCleaner",
                                      src = cms.InputTag("patJets"),
                                      preselection = cms.string('pt > 20.0 && abs(eta) < 2.1'),
                                      checkOverlaps = cms.PSet(
    muons = cms.PSet(
        src       = cms.InputTag("selectedPatMuons"),
    algorithm = cms.string("byDeltaR"),
    preselection        = cms.string(""),
   deltaR              = cms.double(0.5),
   checkRecoComponents = cms.bool(False),
   pairCut             = cms.string(""),
    requireNoOverlaps   = cms.bool(True),
    )
    ),
###applying 2 loose cutson the b-tagged jets

#finalCut = cms.string('')
#                                      finalCut = cms.string('bDiscriminator(\"trackCountingHighEffBJetTags\")>1.7 && bDiscriminator(\"jetProbabilityBJetTags\")>0.988')
                                      finalCut = cms.string('bDiscriminator(\"trackCountingHighEffBJetTags\")>1.7')
                                
                                      )


                                      

process.jetFilter = cms.EDFilter("CandViewCountFilter",
                                     src = cms.InputTag("cleanPatJets"),
                                     minNumber = cms.uint32(2),
#                                     maxNumber = cms.uint32(4)
                                 )

# Z Candidates 
process.zjj = cms.EDProducer("CandViewShallowCloneCombiner",
    checkCharge = cms.bool(False),
    cut = cms.string(''),
    decay = cms.string("cleanPatJets cleanPatJets")
)   


### ntuple
process.load("ZjjEdmNtuples_cff")



##  NSVtexJets 

#EventNjetAndMetInfoNtupleDumper

process.SVtxEdmNtuple= cms.EDProducer(
    "SVtxInfoDumper",
    zjj = cms.InputTag("zjj"),
    )


## # path for dumping svtx info in the ntuple
##   generalEventInfoPath = cms.Path(
##  	  	   process.SVtxEdmNtuple 
##  	  	     )



# Define the relevant paths and schedule them
process.analysisPath = cms.Path(
#    process.cFlavorHistoryProducer +
#    process.bFlavorHistoryProducer +
    #process.flavorHistoryFilter +
 #   process.hltFilter +
    process.makePatMuons +

    process.makePatJets +
    process.makePatMETs +
    
    
    process.selectedPatMuons +
#    process.muonFilter +
    process.cleanPatJets +
    process.jetFilter +
    process.zjj +
    process.SVtxEdmNtuple+ 
    process.ZjjEdmNtuple+
    process.MetEdmNtuple 
  
    )



if isData!=True:
    process.analysisPath = cms.Path(
   # process.hltFilter +
    process.cFlavorHistoryProducer +
    process.bFlavorHistoryProducer +
    process.flavorHistoryFilter +
    
    process.makePatMuons +

    process.makePatJets +
    process.makePatMETs +
    
    
    process.selectedPatMuons +
 #   process.muonFilter +
    process.cleanPatJets +
    process.jetFilter +
    process.zjj +
    process.SVtxEdmNtuple +
    process.ZjjEdmNtuple+
    process.MetEdmNtuple 

)

# Setup for a basic filtering






# Output Module : Hopefully we keep all we need
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('zjj.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring("analysisPath")
        ),
    outputCommands =  cms.untracked.vstring(
        'drop *_*_*_*',
        'keep *_patJets_*_*',
        'keep *_zjj_*_*',
        'keep *_flavorHistoryFilter_*_*',
    )
#    verbose = cms.untracked.bool(True)
)

#process.out.SelectEvents = cms.untracked.PSet(
#        SelectEvents = cms.vstring("filterPath",
#                                   )
#        )



process.edmNtuplesOut.fileName = cms.untracked.string('zjjEdmNtuples.root')
process.edmNtuplesOut.SelectEvents = cms.untracked.PSet(
            SelectEvents = cms.vstring("analysisPath",
                                       )
            )                            

process.endPath = cms.EndPath(process.out +
                              process.edmNtuplesOut)




