

import FWCore.ParameterSet.Config as cms

process = cms.Process("MuEleMetBbbarSkimv2")

# Setup PAT
from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.pfTools import *


process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.load('Configuration.StandardSequences.Services_cff')

# Specify the Global Tag
process.GlobalTag.globaltag = 'START39_V9::All'

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
#        'rfio:/castor/cern.ch/cms/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/291/A4F73A5C-CFE4-DF11-B0F8-001D09F2512C.root',
#"file:/data/1a/degrutto/data/Fall10/DYToBB_M_50/DYToBB_M_50_TuneZ2_1.root"
    "file:/data/1a/degrutto/data/Fall10/VH/ZH_ZToLL_HToBB_M-125.root"
#"rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_9_7/RelValTTbar/GEN-SIM-RECO/START39_V8-v1/0048/DC84E63C-B90D-E011-BF62-002618943836.root"
   # "file:/data/1a/degrutto/data/Fall10/VH/ZH_ZToNuNu_HToBB_M-125.root"
#     "file:/data/1a/degrutto/data/Fall10/VH/WH_WToLNu_HToBB_M-125.root"

)
 )   


process.printList = cms.EDAnalyzer("ParticleListDrawer",
                                 src = cms.InputTag("genParticles"),
                                 maxEventToPrint = cms.untracked.int32(1000)
                                 )

process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
                                            src = cms.InputTag("genParticles"),
                                            printP4 = cms.untracked.bool(False),
                                            printPtEtaPhi = cms.untracked.bool(False),
                                            printVertex = cms.untracked.bool(True),
                                            printStatus = cms.untracked.bool(True),
                                            printIndex = cms.untracked.bool(True),
                                            status = cms.untracked.vint32(1, 2, 3)
                                        )



process.genLep = cms.EDFilter("GenParticleSelector",
  filter = cms.bool(False),
  src = cms.InputTag("genParticles"),
  cut = cms.string("abs(pdgId)= 13  &  status=1")
                              )



process.selGenLep = cms.EDFilter("GenParticleSelector",
                                  src = cms.InputTag("genLep"),
#                                  etaMin= cms.double(-2.4),
#                                  etaMax= cms.double(2.4),
                                  cut = cms.string("abs(eta)<2.4 & pt>10. ")

)

process.genBfromH = cms.EDFilter("GenParticleSelector",
  filter = cms.bool(False),
  src = cms.InputTag("genB"),
  cut = cms.string("mother().mother().pdgId=25")
        )




process.genZstar = cms.EDFilter("PdgIdAndStatusCandViewSelector",
 
  src = cms.InputTag("genParticles"),
  pdgId = cms.vint32(23),
  status = cms.vint32(3)
                                
)







process.zllCands = cms.EDProducer("CandViewCombiner",
    checkCharge = cms.bool(False),
    cut = cms.string('mass>86 && mass<96 && (daughter(0).pdgId = - daughter(1).pdgId) '),
    decay = cms.string("selGenLep selGenLep")
)


# dimuon filter
process.dimuonsFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("zllCands"),
    minNumber = cms.uint32(1)
)


# Switch to using PFMET
switchToPFMET(
        process,
            cms.InputTag('pfMet'),
            ""
        )

removeMCMatching(process, ['All'])

switchJetCollection(process,cms.InputTag('ak5GenJets'),
                                                doJTA        = True,
                                                doBTagging   = True,
                                                jetCorrLabel = ('AK5PF', cms.vstring(['L2Relative', 'L3Absolute' ])),
                                                doType1MET   = True,
                                                genJetCollection=cms.InputTag("ak5GenJets"),
                                                doJetID      = True
                                            )




     

# Clean the Jets from the selected leptons, applying here only pt cut pf 20 
process.cleanPatJets = cms.EDProducer("PATJetCleaner",
                                      src = cms.InputTag("patJets"),
                                      # preselection = cms.string('pt > 20.0 && abs(eta) < 2.5 && ( (neutralEmEnergy/energy < 0.99) &&  (neutralHadronEnergy/energy < 0.99) && numberOfDaughters>1) '),
                                      preselection = cms.string('pt > 0.0'),
                                      checkOverlaps = cms.PSet(
     muons = cms.PSet(
     src       = cms.InputTag("selGenLep"),
     algorithm = cms.string("byDeltaR"),
     preselection        = cms.string(""),
     deltaR              = cms.double(0.5),
     checkRecoComponents = cms.bool(True),
     pairCut             = cms.string(""),
     requireNoOverlaps   = cms.bool(True),
    )
    ),
                                      finalCut = cms.string('')
                                      
                                      )



  


### select pat jets, applying eta loose ID and  b-tag
process.selectedPatJets = cms.EDFilter("PATJetSelector",
                               src= cms.InputTag("cleanPatJets"),
                               cut = cms.string("pt>0" )
                                       #&& bDiscriminator(\"trackCountingHighEffBJetTags\")>1.7 " ) ### try without loose b-tag now..... 
                               )

                                                    

process.jetFilter = cms.EDFilter("CandViewCountFilter",
                                     src = cms.InputTag("selectedPatJets"),
                                     minNumber = cms.uint32(2),
                                   
                                 )


process.myPartons = cms.EDProducer("PartonSelector",
                                           withLeptons = cms.bool(True),
                                   src =  cms.InputTag("genParticles"),
                                       )
 
process.flavourByRef = cms.EDProducer("JetPartonMatcher",
                                              jets = cms.InputTag("selectedPatJets"),
                                              coneSizeToAssociate = cms.double(0.3),
                                              partons = cms.InputTag("myPartons")
                                          )
 
process.flavourByVal = cms.EDProducer("JetFlavourIdentifier",
                                         srcByReference = cms.InputTag("flavourByRef"),
                                          physicsDefinition = cms.bool(False)
                                      )
 

process.genB = cms.EDFilter("GenParticleSelector",
  filter = cms.bool(False),
  src = cms.InputTag("genParticles"),
  cut = cms.string("abs(pdgId)= 5 &  status=2 & mother().numberOfMothers()>0 & pt>0 & abs(eta)<2.5")
                            )

process.genBfromH = cms.EDFilter("GenParticleSelector",
  filter = cms.bool(False),
  src = cms.InputTag("genB"),
  cut = cms.string("mother().mother().pdgId=25")
        )                                         


process.bbFilter = cms.EDFilter("CandViewCountFilter",
                                     src = cms.InputTag("genBfromH"),
                                     minNumber = cms.uint32(2),
                                    
                                 )

process.patJetPartonMatch = cms.EDProducer("MCMatcherByPt",        # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                   src         = cms.InputTag("selectedPatJets"),       # RECO objects to match
                                   matched     = cms.InputTag("genBfromH"),      # mc-truth particle collection
                                   mcPdgId     = cms.vint32(5,-5),     # one or more PDG ID (quarks except top; gluons)
                                   mcStatus    = cms.vint32(2),                     # PYTHIA status code (3 = hard scattering)
                                   checkCharge = cms.bool(False),                   # False = any value of the charge of MC and RECO is ok
                                   maxDeltaR   = cms.double(0.1),                   # Minimum deltaR for the match
                                   maxDPtRel   = cms.double(0.5),                   # Minimum deltaPt/Pt for the match
                                   resolveAmbiguities    = cms.bool(True),          # Forbid two RECO objects to match to the same GEN object
                                   resolveByMatchQuality = cms.bool(True),         # False = just match input in order; True = pick lowest deltaR/deltaPt? pair first
                                   )

process.BNtuple = cms.EDProducer("BInfoDumper",
                                map = cms.InputTag("patJetPartonMatch"),
                                jet  = cms.InputTag("selectedPatJets"),
                                z = cms.InputTag("zllCands"),
                                allJet =   cms.InputTag("patJets"),
                                met = cms.InputTag("patMETs")
                                 )

 



process.hbbCands = cms.EDProducer("CandViewCombiner",
    checkCharge = cms.bool(True),
    cut = cms.string('mass>0'),
    decay = cms.string("patJetPartonMatch patJetPartonMatch")
)


process.bFilter = cms.EDFilter("CandViewCountFilter",
                                     src = cms.InputTag("hbbCands"),
                                     minNumber = cms.uint32(1)

                                 )


process.hjjCands = cms.EDProducer("CandViewShallowCloneCombiner",
    checkCharge = cms.bool(False),
    cut = cms.string('mass>0'),
    decay = cms.string("selectedPatJets selectedPatJets")
)   


process.hjjzllCands = cms.EDProducer("CandViewCombiner",
    checkCharge = cms.bool(False),
    cut = cms.string(''),
    decay = cms.string("hjjCands zllCands")
)


process.hbbzllCands = cms.EDProducer("CandViewCombiner",
    checkCharge = cms.bool(False),
    cut = cms.string(''),
    decay = cms.string("hbbCands zllCands")
)   


process.hjjzllNtuple = cms.EDProducer(
    "CandViewNtpProducer",
    src = cms.InputTag("hjjzllCands"),
    variables = cms.VPSet(
    cms.PSet(
    tag = cms.untracked.string("jjllmass"),
    quantity = cms.untracked.string("mass")
    ),
    cms.PSet(
    tag = cms.untracked.string("jjllpt"),
    quantity = cms.untracked.string("pt")
    ),
    )
    )

process.hbbzllNtuple = cms.EDProducer(
    "CandViewNtpProducer",
    src = cms.InputTag("hbbzllCands"),
    variables = cms.VPSet(
    cms.PSet(
    tag = cms.untracked.string("bbllmass"),
    quantity = cms.untracked.string("mass")
    ),
    cms.PSet(
    tag = cms.untracked.string("bbllpt"),
    quantity = cms.untracked.string("pt")
    ),
    )
    )


process.hjjNtuple = cms.EDProducer(
    "CandViewNtpProducer",
    src = cms.InputTag("hjjCands"),
    variables = cms.VPSet(
    cms.PSet(
    tag = cms.untracked.string("jjmass"),
    quantity = cms.untracked.string("mass")
    ),
    cms.PSet(
    tag = cms.untracked.string("jjpt"),
    quantity = cms.untracked.string("pt")
    ),
    cms.PSet(
    tag = cms.untracked.string("jjDau1pt"),
    quantity = cms.untracked.string("daughter(0).pt")

    ),
    cms.PSet(
    tag = cms.untracked.string("jjDau2pt"),
    quantity = cms.untracked.string("daughter(1).pt")
    ),

    
    )
    )



process.hbbNtuple = cms.EDProducer(
    "CandViewNtpProducer",
    src = cms.InputTag("hbbCands"),
    variables = cms.VPSet(
    cms.PSet(
    tag = cms.untracked.string("bbmass"),
    quantity = cms.untracked.string("mass")
    ),
    cms.PSet(
    tag = cms.untracked.string("bbpt"),
    quantity = cms.untracked.string("pt")
    ),
    cms.PSet(
    tag = cms.untracked.string("bbDau1pt"),
    quantity = cms.untracked.string("daughter(0).pt")
    ),
    cms.PSet(
    tag = cms.untracked.string("bbDau2pt"),
    quantity = cms.untracked.string("daughter(1).pt")
    ),
    
    )
    )


process.zllNtuple = cms.EDProducer(
    "CandViewNtpProducer",
    src = cms.InputTag("zllCands"),
    variables = cms.VPSet(
    cms.PSet(
    tag = cms.untracked.string("llmass"),
    quantity = cms.untracked.string("mass")
    ),
    cms.PSet(
    tag = cms.untracked.string("llpt"),
    quantity = cms.untracked.string("pt")
    ),
        cms.PSet(
    tag = cms.untracked.string("lleta"),
    quantity = cms.untracked.string("eta")
    ),
        cms.PSet(
    tag = cms.untracked.string("llphi"),
    quantity = cms.untracked.string("phi")
    ),
    )
    )
        
process.genZstarNtuple  = cms.EDProducer(
    "CandViewNtpProducer",
    src = cms.InputTag("genZstar"),
    variables = cms.VPSet(
    cms.PSet(
    tag = cms.untracked.string("genZstarMass"),
    quantity = cms.untracked.string("?mass>200?mass:-1")
    ),
     cms.PSet(
    tag = cms.untracked.string("genZstarPt"),
    quantity = cms.untracked.string("?mass>200?pt:-1")
    )
    )
)



process.ZllHbbPath = cms.Path(
    # process.printList *    
     process.genLep +
     process.selGenLep +
     process.zllCands +
     process.dimuonsFilter +
      process.makePatMETs +
     process.makePatJets +
     process.cleanPatJets +
     process.selectedPatJets +
     process.jetFilter +
     process.hjjCands *
     process.hjjzllCands *
     process.genB *
     process.genBfromH *
     process.bbFilter *
#     process.selGenB *
     process.myPartons *
  #    process.flavourByRef *
  #   process.flavourByVal*
     process.patJetPartonMatch *
     process.BNtuple *
     process.genZstar *
     process.genZstarNtuple *
     
     #process.hbbCands *
    # process.bFilter *
    # process.hbbzllCands *
     process.hjjzllNtuple  *
     process.hjjNtuple *
    # process.hbbNtuple *
    # process.hbbzllNtuple *     
     process.zllNtuple *
    #  process.printTree *
    process.printList
     
 )



process.out = cms.OutputModule(
        "PoolOutputModule",
            fileName = cms.untracked.string('Ntuple_test_v3.root'),
        outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_hjjzllNtuple_*_*",
    "keep *_hbbzllNtuple_*_*",
    "keep *_hjjNtuple_*_*",
    "keep *_BNtuple_*_*",
    "keep *_zllNtuple_*_*",
    "keep *_genZstarNtuple_*_*"
    ),
            SelectEvents = cms.untracked.PSet(
          SelectEvents = cms.vstring(
            "ZllHbbPath",
            )
          )
        )

process.endPath = cms.EndPath(
    
     process.out
    )

