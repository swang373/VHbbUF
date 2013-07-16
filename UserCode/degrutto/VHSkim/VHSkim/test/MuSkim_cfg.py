
import FWCore.ParameterSet.Config as cms

process = cms.Process("MuSkim")


process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 100




import HLTrigger.HLTfilters.hltHighLevel_cfi
process.muHLTFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.muHLTFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.muHLTFilter.HLTPaths = ["HLT_Mu15v1"]

# Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Source file : To be run on a Full RECO sample
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'rfio:/castor/cern.ch/cms/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/291/A4F73A5C-CFE4-DF11-B0F8-001D09F2512C.root',

    )
)




  
process.goodGlobalMuons = cms.EDFilter("MuonViewRefSelector",
                                 src = cms.InputTag("muons"),
                                 cut = cms.string('isGlobalMuon = 1 & isTrackerMuon = 1 &  pt > 15 & abs(eta)<2.4 & dB<1.0 & globalTrack().hitPattern().numberOfValidTrackerHits>7  & globalTrack().hitPattern().numberOfValidMuonHits>0'),
                                 filter = cms.bool(True)
                               )

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.MuFilterPath = cms.Path(
    #process.muHLTFilter *
        process.goodGlobalMuons

        )



process.load("Configuration.EventContent.EventContent_cff")



process.MuSkimEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring()
)

# to add the AOD output uncomment the following
process.MuSkimEventContent.outputCommands.extend(process.AODEventContent.outputCommands)


process.MuSkimOutputModule = cms.OutputModule("PoolOutputModule",
    process.MuSkimEventContent,
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
           'MuFilterPath',
           )
    ),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('muSkim'),
        dataTier = cms.untracked.string('USER')
   ),
   fileName = cms.untracked.string('MuSkim.root')
)




process.endPath = cms.EndPath(
    process.MuSkimOutputModule

    )

