
import FWCore.ParameterSet.Config as cms
import copy

### Zjj  variables


ZjjEdmNtuple =  cms.EDProducer(
    "CandViewNtpProducer",
    src=cms.InputTag("zjj"),
    lazyParser=cms.untracked.bool(True),
    prefix=cms.untracked.string("zjj"),
    eventInfo=cms.untracked.bool(True),
    variables = cms.VPSet(
    cms.PSet(
    tag = cms.untracked.string("Mass"),
    quantity = cms.untracked.string("mass")
    ),
    cms.PSet(
    tag = cms.untracked.string("Pt"),
    quantity = cms.untracked.string("pt")
    ),
    cms.PSet(
    tag = cms.untracked.string("Eta"),
    quantity = cms.untracked.string("eta")
    ),
    cms.PSet(
    tag = cms.untracked.string("Phi"),
    quantity = cms.untracked.string("phi")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1Pt"),
    quantity = cms.untracked.string("daughter(0).pt ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2Pt"),
    quantity = cms.untracked.string("daughter(1).pt")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1Eta"),
    quantity = cms.untracked.string("daughter(0).eta ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2Eta"),
    quantity = cms.untracked.string("daughter(1).eta ")
    ),cms.PSet(
    tag = cms.untracked.string("Jet1Phi"),
    quantity = cms.untracked.string("daughter(0).phi ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2Phi"),
    quantity = cms.untracked.string("daughter(1).phi ")
    ),

    ### b tagging ###
    cms.PSet(
    tag = cms.untracked.string("Jet1CSV"),
    quantity = cms.untracked.string("daughter(0).bDiscriminator(\"combinedSecondaryVertexBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2CSV"),
    quantity = cms.untracked.string("daughter(1).bDiscriminator(\"combinedSecondaryVertexBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1CSVMVA"),
    quantity = cms.untracked.string("daughter(0).bDiscriminator(\"combinedSecondaryVertexMVABJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2CSVMVA"),
    quantity = cms.untracked.string("daughter(1).bDiscriminator(\"combinedSecondaryVertexMVABJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1JProb"),
    quantity = cms.untracked.string("daughter(0).bDiscriminator(\"jetProbabilityBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2JProb"),
    quantity = cms.untracked.string("daughter(1).bDiscriminator(\"jetProbabilityBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1JbProb"),
    quantity = cms.untracked.string("daughter(0).bDiscriminator(\"jetBProbabilityBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2JbProb"),
    quantity = cms.untracked.string("daughter(1).bDiscriminator(\"jetBProbabilityBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1SSVHE"),
    quantity = cms.untracked.string("daughter(0).bDiscriminator(\"simpleSecondaryVertexHighEffBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2SSVHE"),
    quantity = cms.untracked.string("daughter(1).bDiscriminator(\"simpleSecondaryVertexHighEffBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1SSVHP"),
    quantity = cms.untracked.string("daughter(0).bDiscriminator(\"simpleSecondaryVertexHighPurBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2SSVHP"),
    quantity = cms.untracked.string("daughter(1).bDiscriminator(\"simpleSecondaryVertexHighPurBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1ElPt"),
    quantity = cms.untracked.string("daughter(0).bDiscriminator(\"softElectronByPtBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2ElPt"),
    quantity = cms.untracked.string("daughter(1).bDiscriminator(\"softElectronByPtBJetTags\") ")
    ),
        cms.PSet(
    tag = cms.untracked.string("Jet1ElIp"),
    quantity = cms.untracked.string("daughter(0).bDiscriminator(\"softElectronByIP3dBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2ElIp"),
    quantity = cms.untracked.string("daughter(1).bDiscriminator(\"softElectronByIP3dBJetTags\") ")
    ),
        cms.PSet(
    tag = cms.untracked.string("Jet1Mu"),
    quantity = cms.untracked.string("daughter(0).bDiscriminator(\"softMuonBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2Mu"),
    quantity = cms.untracked.string("daughter(1).bDiscriminator(\"softMuonBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1MuPt"),
    quantity = cms.untracked.string("daughter(0).bDiscriminator(\"softMuonByPtBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2MuPt"),
    quantity = cms.untracked.string("daughter(1).bDiscriminator(\"softMuonByPtBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1MuIp"),
    quantity = cms.untracked.string("daughter(0).bDiscriminator(\"softMuonByIP3dBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2MuIp"),
    quantity = cms.untracked.string("daughter(1).bDiscriminator(\"softMuonByIP3dBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1TKHE"),
    quantity = cms.untracked.string("daughter(0).bDiscriminator(\"trackCountingHighEffBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2TKHE"),
    quantity = cms.untracked.string("daughter(1).bDiscriminator(\"trackCountingHighEffBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1TKHP"),
    quantity = cms.untracked.string("daughter(0).bDiscriminator(\"trackCountingHighPurBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2TKHP"),
    quantity = cms.untracked.string("daughter(1).bDiscriminator(\"trackCountingHighPurBJetTags\") ")
    ),


    )
    )


MetEdmNtuple =  cms.EDProducer(
        "CandViewNtpProducer",
        src=cms.InputTag("pfMet"),
        lazyParser=cms.untracked.bool(True),
        prefix=cms.untracked.string("pfMet"),
        eventInfo=cms.untracked.bool(False),
        variables = cms.VPSet(
    cms.PSet(
    tag = cms.untracked.string("Pt"),
    quantity = cms.untracked.string("pt")
    ),
    )
        )


edmNtuplesOut = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('zbbEdmNtuples.root'),
    outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_ZjjEdmNtuple_*_*",
    "keep *_MetEdmNtuple_*_*",
    "keep *_SVtxEdmNtuple_*_*",
    
    )
    )
