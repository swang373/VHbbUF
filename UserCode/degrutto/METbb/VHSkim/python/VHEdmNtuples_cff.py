import FWCore.ParameterSet.Config as cms
import copy






###  hjj standard and bDiscriminator variables

hjj = (
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
    tag = cms.untracked.string("Jet1Mass"),
    quantity = cms.untracked.string("daughter(0).mass")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2Mass"),
    quantity = cms.untracked.string("daughter(1).mass")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1Eta"),
    quantity = cms.untracked.string("daughter(0).eta ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2Eta"),
    quantity = cms.untracked.string("daughter(1).eta ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1Phi"),
    quantity = cms.untracked.string("daughter(0).phi ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2Phi"),
    quantity = cms.untracked.string("daughter(1).phi ")
    ),
#    cms.PSet(
#    tag = cms.untracked.string("Jet1HLTBit"),
#    quantity = cms.untracked.string("daughter(0).masterClone.triggerObjectMatchesByPath('HLT_MET45_2CentralJet20_2BTagIP_SingleTrackTC_v1').size()")
#    ),
#    cms.PSet(
#    tag = cms.untracked.string("Jet2HLTBit"),
#    quantity = cms.untracked.string("daughter(1).masterClone.triggerObjectMatchesByPath('HLT_MET45_2CentralJet20_2BTagIP_SingleTrackTC_v1').size()")
#    ),
  
    ### b tagging ###
    cms.PSet(
    tag = cms.untracked.string("Jet1CSV"),
    quantity = cms.untracked.string("daughter(0).masterClone.bDiscriminator(\"combinedSecondaryVertexBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2CSV"),
    quantity = cms.untracked.string("daughter(1).masterClone.bDiscriminator(\"combinedSecondaryVertexBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1CSVMVA"),
    quantity = cms.untracked.string("daughter(0).masterClone.bDiscriminator(\"combinedSecondaryVertexMVABJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2CSVMVA"),
    quantity = cms.untracked.string("daughter(1).masterClone.bDiscriminator(\"combinedSecondaryVertexMVABJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1JProb"),
    quantity = cms.untracked.string("daughter(0).masterClone.bDiscriminator(\"jetProbabilityBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2JProb"),
    quantity = cms.untracked.string("daughter(1).masterClone.bDiscriminator(\"jetProbabilityBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1JbProb"),
    quantity = cms.untracked.string("daughter(0).masterClone.bDiscriminator(\"jetBProbabilityBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2JbProb"),
    quantity = cms.untracked.string("daughter(1).masterClone.bDiscriminator(\"jetBProbabilityBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1SSVHE"),
    quantity = cms.untracked.string("daughter(0).masterClone.bDiscriminator(\"simpleSecondaryVertexHighEffBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2SSVHE"),
    quantity = cms.untracked.string("daughter(1).masterClone.bDiscriminator(\"simpleSecondaryVertexHighEffBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1SSVHP"),
    quantity = cms.untracked.string("daughter(0).masterClone.bDiscriminator(\"simpleSecondaryVertexHighPurBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2SSVHP"),
    quantity = cms.untracked.string("daughter(1).masterClone.bDiscriminator(\"simpleSecondaryVertexHighPurBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1ElPt"),
    quantity = cms.untracked.string("daughter(0).masterClone.bDiscriminator(\"softElectronByPtBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2ElPt"),
    quantity = cms.untracked.string("daughter(1).masterClone.bDiscriminator(\"softElectronByPtBJetTags\") ")
    ),
        cms.PSet(
    tag = cms.untracked.string("Jet1ElIp"),
    quantity = cms.untracked.string("daughter(0).masterClone.bDiscriminator(\"softElectronByIP3dBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2ElIp"),
    quantity = cms.untracked.string("daughter(1).masterClone.bDiscriminator(\"softElectronByIP3dBJetTags\") ")
    ),
        cms.PSet(
    tag = cms.untracked.string("Jet1Mu"),
    quantity = cms.untracked.string("daughter(0).masterClone.bDiscriminator(\"softMuonBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2Mu"),
    quantity = cms.untracked.string("daughter(1).masterClone.bDiscriminator(\"softMuonBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1MuPt"),
    quantity = cms.untracked.string("daughter(0).masterClone.bDiscriminator(\"softMuonByPtBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2MuPt"),
    quantity = cms.untracked.string("daughter(1).masterClone.bDiscriminator(\"softMuonByPtBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1MuIp"),
    quantity = cms.untracked.string("daughter(0).masterClone.bDiscriminator(\"softMuonByIP3dBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2MuIp"),
    quantity = cms.untracked.string("daughter(1).masterClone.bDiscriminator(\"softMuonByIP3dBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1TKHE"),
    quantity = cms.untracked.string("daughter(0).masterClone.bDiscriminator(\"trackCountingHighEffBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2TKHE"),
    quantity = cms.untracked.string("daughter(1).masterClone.bDiscriminator(\"trackCountingHighEffBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1TKHP"),
    quantity = cms.untracked.string("daughter(0).masterClone.bDiscriminator(\"trackCountingHighPurBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2TKHP"),
    quantity = cms.untracked.string("daughter(1).masterClone.bDiscriminator(\"trackCountingHighPurBJetTags\") ")
    ),


    )
    
    

template =  cms.EDProducer(
    "CandViewNtpProducer",
    src=cms.InputTag("test"),
    lazyParser=cms.untracked.bool(True),
    prefix=cms.untracked.string("template"),
    eventInfo=cms.untracked.bool(True),
    variables = cms.VPSet(
    ),
   )
  



HjjEdmNtuple = copy.deepcopy(template)
HjjEdmNtuple.variables += hjj
HjjEdmNtuple.src = cms.InputTag("hjj")
HjjEdmNtuple.prefix = cms.untracked.string("hjj")


### met info
MetEdmNtuple =  cms.EDProducer(
        "CandViewNtpProducer",
        src=cms.InputTag("selectedMETTriggerMatch"),
        lazyParser=cms.untracked.bool(True),
        prefix=cms.untracked.string("pfMet"),
        eventInfo=cms.untracked.bool(False),
        variables = cms.VPSet(
    cms.PSet(
    tag = cms.untracked.string("Et"),
    quantity = cms.untracked.string("et")
    ),
    cms.PSet(
    tag = cms.untracked.string("Sig"),
    quantity = cms.untracked.string("mEtSig"),
    ),
     cms.PSet(
    tag = cms.untracked.string("Phi"),
    quantity = cms.untracked.string("phi"),
    ),
    
        )
        )


    
HjjMetEdmNtuple = cms.EDProducer(
    "CandViewNtpProducer",
    src=cms.InputTag("hjjmet"),
    lazyParser=cms.untracked.bool(True),
    prefix=cms.untracked.string("hjjmet"),
    eventInfo=cms.untracked.bool(False),
    variables = cms.VPSet(
    cms.PSet(
    tag = cms.untracked.string("MT"),
    quantity = cms.untracked.string("sqrt((daughter(0).pt+daughter(1).pt)*(daughter(0).pt+daughter(1).pt)-pt*pt)")
    ),
    cms.PSet(
    tag = cms.untracked.string("Et"),
    quantity = cms.untracked.string("et"),
    ),
    cms.PSet(
    tag = cms.untracked.string("Eta"),
    quantity = cms.untracked.string("eta")
    ),
    cms.PSet(
    tag = cms.untracked.string("Y"),
    quantity = cms.untracked.string("y"),
    ),
    cms.PSet(
    tag = cms.untracked.string("Phi"),
    quantity = cms.untracked.string("phi"),
    )
    )
    )




       
edmNtuplesOut = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('2l2bMetEdmNtuples.root'),
    outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_HjjEdmNtuple_*_*",
    "keep *_HjjMetEdmNtuple_*_*",
    "keep *_MetEdmNtuple_*_*",
    "keep *_HLTInfoMetJetEdmNtuple_*_*",
    
    )
    )
