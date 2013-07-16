
import FWCore.ParameterSet.Config as cms

process = cms.Process("MuEleMetBbbarSkimv2")

# Setup PAT
from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.pfTools import *


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
process.hltFilter.HLTPaths = ["HLT_Mu11", "HLT_Ele17_SW_CaloEleId_L1R"]

# Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Source file : To be run on a Full RECO sample
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
  #  "rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_9_7/RelValTTbar/GEN-SIM-RECO/START39_V8-v1/0048/DC84E63C-B90D-E011-BF62-002618943836.root"
#        'rfio:/castor/cern.ch/cms/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/291/A4F73A5C-CFE4-DF11-B0F8-001D09F2512C.root',
#"file:/data/1a/degrutto/data/Fall10/DYToBB_M_50/DYToBB_M_50_TuneZ2_1.root"
"rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_1_4//RelValTTbar/GEN-SIM-RECO/START311_V2-v1/0013/2A878D65-7D60-E011-9308-00261894393D.root"
#    "file:/data/1a/degrutto/data/Fall10/VH/ZH_ZToLL_HToBB_M-125.root"
   # "file:/data/1a/degrutto/data/Fall10/VH/ZH_ZToNuNu_HToBB_M-125.root"
#     "file:/data/1a/degrutto/data/Fall10/VH/WH_WToLNu_HToBB_M-125.root"


    

    

    )
)


isData = True

if isData:
    removeMCMatching(process, ['All'])

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




# Modules for the Cut-based Electron ID in the VBTF prescription
import ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi as vbtfid

process.eidVBTFRel95 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '95relIso' )
process.eidVBTFRel80 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '80relIso' )
process.eidVBTFCom95 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '95cIso'   )
process.eidVBTFCom80 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '80cIso'   )

process.eidSequence = cms.Sequence(
    process.eidVBTFRel95 +
    process.eidVBTFRel80 +
    process.eidVBTFCom95 +
    process.eidVBTFCom80 
)


# Electron Selection
process.patElectrons.electronIDSources = cms.PSet(
    eidVBTFRel95 = cms.InputTag("eidVBTFRel95"),
    eidVBTFRel80 = cms.InputTag("eidVBTFRel80"),
    eidVBTFCom95 = cms.InputTag("eidVBTFCom95"),
    eidVBTFCom80 = cms.InputTag("eidVBTFCom80")
)

process.selectedPatElectrons.cut = ( 
    "pt > 15.0 && abs(eta) < 2.5 &&"                               +
    "(isEE || isEB) && !isEBEEGap &&"                              +
    "electronID('eidVBTFCom95') == 7" 
)






  
## now selecting a good MET!!!
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
                            jetCorrLabel = ('AK5PF', cms.vstring(['L2Relative', 'L3Absolute' ])),
                            doType1MET   = True,
                            genJetCollection=cms.InputTag("ak5GenJets"),
                            doJetID      = True
                        )
    


####process.patJets.addTagInfos  = True # ### needed for secondary vertex variables....

### ntuple
process.load("HSkim.VHSkim.DiLeptEdmNtuples_cff")




# Clean the Jets from the selected leptons, applying here only pt cut pf 20 
process.cleanPatJets = cms.EDProducer("PATJetCleaner",
                                      src = cms.InputTag("patJets"),
                                     # preselection = cms.string('pt > 20.0 && abs(eta) < 2.5 && ( (neutralEmEnergy/energy < 0.99) &&  (neutralHadronEnergy/energy < 0.99) && numberOfDaughters>1) '),
                                         preselection = cms.string('pt > 20.0 & (neutralEmEnergy/energy < 0.99) &&  (neutralHadronEnergy/energy < 0.99) && numberOfDaughters>1'),
                                      checkOverlaps = cms.PSet(
    muons = cms.PSet(
    src       = cms.InputTag("selectedPatMuons"),
    algorithm = cms.string("byDeltaR"),
    preselection        = cms.string(""),
    deltaR              = cms.double(0.5),
    checkRecoComponents = cms.bool(False),
    pairCut             = cms.string(""),
    requireNoOverlaps   = cms.bool(True),
    ),
    electrons = cms.PSet(
    src       = cms.InputTag("selectedPatElectrons"),
    algorithm = cms.string("byDeltaR"),
    preselection        = cms.string(""),
    deltaR              = cms.double(0.5),
    checkRecoComponents = cms.bool(False),
    pairCut             = cms.string(""),
    requireNoOverlaps   = cms.bool(True),
    ),
    ),
                                      ###applying 2 loose cutson the b-tagged jets
                                      
                                      finalCut = cms.string('')
                                      ##finalCut = cms.string('bDiscriminator(\"trackCountingHighEffBJetTags\")>1.7 ')
                                      )
                                      ### select pat jets, applying eta loose ID and  b-tag
process.myPatJets = cms.EDFilter("PATJetSelector",
                               src= cms.InputTag("cleanPatJets"),
                               cut = cms.string("abs(eta)< 10. & pt > 20.0" ) #&& bDiscriminator(\"trackCountingHighEffBJetTags\")>1.7 " ) ### try without loose b-tag now..... 
                               )
                                      



# Z Candidates and Higgs Candidates
process.zee = cms.EDProducer("CandViewShallowCloneCombiner",
    checkCharge = cms.bool(True),
    cut = cms.string('mass > 0'),
    decay = cms.string("selectedPatElectrons@+ selectedPatElectrons@-")
)

process.zmm = cms.EDProducer("CandViewShallowCloneCombiner",
    checkCharge = cms.bool(True),
    cut = cms.string('mass > 0'),
    decay = cms.string("selectedPatMuons@+ selectedPatMuons@-")
)


process.zem = cms.EDProducer("CandViewShallowCloneCombiner",
    checkCharge = cms.bool(True),
    cut = cms.string('mass > 0'),
    decay = cms.string("selectedPatMuons selectedPatElectrons")
)

# Setup for a basic filtering
process.zll = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag("zee", "zmm", "zem")
) 

process.zllFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("zll"),
    minNumber = cms.uint32(1),
)




process.JetLepEdmNtuple= cms.EDProducer(
    "JetLepInfoDumper",
    j = cms.InputTag("myPatJets"), ### looking after the selections 
    m = cms.InputTag("selectedPatMuons"), 
    e = cms.InputTag("selectedPatElectrons"), 
    )


process.ZllHbbPath = cms.Path(
 #    process.hltFilter +
     process.eidSequence + 
     process.patDefaultSequence+
     
     process.zmm +
     process.zee +
     process.zem +
     process.zll +
     process.zllFilter + 
     process.cleanPatJets +
     process.myPatJets +  
     process.ZmmEdmNtuple +
     process.ZeeEdmNtuple+
     process.ZemEdmNtuple+
     process.MetEdmNtuple +
     process.jEdmNtuple + 
     process.JetLepEdmNtuple 
 )




process.MuEleMetBbbarSkimEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring(
    'keep *_selectedPatElectrons_*_*',
    'keep *_selectedPatMuons_*_*',
    'keep *_selectedPatMET_*_*',
    'keep *_cleanPatJets_*_*',
    'keep *_zee_*_*',
    'keep *_zmm_*_*',
    'keep *_zem_*_*',

    

    )
)

# to add the AOD output uncomment the following
process.MuEleMetBbbarSkimEventContent.outputCommands.extend(process.AODEventContent.outputCommands)


process.out = cms.OutputModule("PoolOutputModule",
    process.MuEleMetBbbarSkimEventContent,
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('ZllHbbPath' )
        ),
                                                       
       dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('muelemetbbarSkim'),
        dataTier = cms.untracked.string('USER')
   ),
   fileName = cms.untracked.string('MuEleMetBbbarSkim.root')
)







process.edmNtuplesOut.fileName = cms.untracked.string('dileptonsNtuples.root')
process.edmNtuplesOut.SelectEvents = cms.untracked.PSet(
            SelectEvents = cms.vstring('ZllHbbPath',
                                       )
            )                            

process.endPath = cms.EndPath(process.out +
                              process.edmNtuplesOut)


