import FWCore.ParameterSet.Config as cms


process = cms.Process("GEN2")

# Number of events to be generated
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


process.filt = cms.EDFilter("SigEventFilter")
process.f = cms.Path(process.filt)
# To write out events (not need: FastSimulation _is_ fast!)
process.o1 = cms.OutputModule(
    "PoolOutputModule",
    #outputCommands = cms.untracked.vstring('keep *'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('f')
    ),

    fileName = cms.untracked.string('CSCTFNNN.root')
)
process.outpath = cms.EndPath(process.o1)


# Make the job crash in case of missing product
process.options = cms.untracked.PSet( Rethrow = cms.untracked.vstring('ProductNotFound') )



# Famos with tracks

process.schedule = cms.Schedule(process.f,process.outpath)

process.source = cms.Source("PoolSource",
                fileNames = cms.untracked.vstring(



