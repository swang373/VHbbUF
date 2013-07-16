
import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process("CosmicLHCAnalyzer")

process.TFileService=cms.Service(
    "TFileService",
    fileName=cms.string("CosmicsLHC_lhcSTABarrel66692_PtEtaCuts_v3.root")

    )

#process.load("ElectroWeakAnalysis.ZReco.SP_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
    )

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    #"rfio:/castor/cern.ch/user/d/degrutto/cosmics/DimuonSkimCosmicsSP_100K.root"
    "rfio:/castor/cern.ch/user/d/degrutto/cosmics/DimuonSkimCosmics66692SP.root"
     )
    )

process.cosmicAnalyzer = cms.EDAnalyzer(
    "CosmicLHCAnalyzer",
#    src_muons = cms.InputTag("lhcStandAloneMuonsBarrelOnly:UpdatedAtVtx"),
    src_muons = cms.InputTag("lhcStandAloneMuonsBarrelOnly"),
    ptCut = cms.double(20.0),
    etaCut = cms.double(2.0),
    massCut = cms.double(20.0),
    dzCut = cms.double(1000.0),
    dxyCut = cms.double(1000.0)
    )

process.path = cms.Path(
    process.cosmicAnalyzer 
    )

## process.out = cms.OutputModule("PoolOutputModule",
##     outputCommands = cms.untracked.vstring(
##         'keep *_lhcStandAloneMuonsBarrelOnly_*_*',
##         'keep *_lhcStandAloneMuonsEndCapsOnly_*_*' 
##         ),
##     fileName= cms.untracked.string('lhcSta.root')
## )

## process.outpath = cms.EndPath(process.out)
