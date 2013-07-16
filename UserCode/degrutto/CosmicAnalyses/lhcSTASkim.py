
import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process("lhcSTASkim")



process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
    )

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
    )

process.load("ElectroWeakAnalysis.ZReco.SP_cff")

## process.source = cms.Source(
##     "PoolSource",
##     fileNames = cms.untracked.vstring(
##     "rfio:/castor/cern.ch/user/d/degrutto/cosmics/DimuonSkimCosmicsSP_100K.root"
##      )
##     )

## RawEventContent = cms.PSet(
##    outputCommands = cms.untracked.vstring(
##         'drop *_*_*_PU'
##             )
## )


process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring(
        'keep *_lhcStandAloneMuonsBarrelOnly_*_*',
        'keep *_lhcStandAloneMuonsEndCapsOnly_*_*' 
        ),
    fileName= cms.untracked.string('SkimCosmiclhcSta.root')
)

## process.out = cms.OutputModule("PoolOutputModule",
##     outputCommands = cms.untracked.vstring(
##         'drop *_*_*_PU'
##                ),
##     fileName= cms.untracked.string('SkimCosmiclhcSta.root')
## )
## process.path = cms.Path(
##        process.out
##       process.out
##        )

#process.outpath = cms.EndPath(process.out)
