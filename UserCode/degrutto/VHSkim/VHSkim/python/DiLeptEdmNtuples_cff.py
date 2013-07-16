import FWCore.ParameterSet.Config as cms
import copy





### Zll variables

zll = (
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
    tag = cms.untracked.string("Y"),
    quantity = cms.untracked.string("y")
    ),
    cms.PSet(
    tag = cms.untracked.string("Phi"),
    quantity = cms.untracked.string("phi")
    ),
    ## pt
    cms.PSet(
    tag = cms.untracked.string("LeptDau1Pt"),
    quantity = cms.untracked.string("daughter(0).pt ")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2Pt"),
    quantity = cms.untracked.string("daughter(1).pt")
    ),
    ## eta
    cms.PSet(
    tag = cms.untracked.string("LeptDau1Eta"),
    quantity = cms.untracked.string("daughter(0).eta ")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2Eta"),
    quantity = cms.untracked.string("daughter(1).eta ")
    ),
    ## phi
    cms.PSet(
    tag = cms.untracked.string("LeptDau1Phi"),
    quantity = cms.untracked.string("daughter(0).phi ")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2Phi"),
    quantity = cms.untracked.string("daughter(1).phi ")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1GlobalMuonBit"),
    quantity = cms.untracked.string("daughter(0).isGlobalMuon")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2GlobalMuonBit"),
    quantity = cms.untracked.string("daughter(1).isGlobalMuon")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1StandAloneBit"),
    quantity = cms.untracked.string("daughter(0).isStandAloneMuon")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2StandAloneBit"),
    quantity = cms.untracked.string("daughter(1).isStandAloneMuon")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1TrackerMuonBit"),
    quantity = cms.untracked.string("daughter(0).isTrackerMuon")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2TrackerMuonBit"),
    quantity = cms.untracked.string("daughter(1).isTrackerMuon")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1NofMuonHits"),
    quantity = cms.untracked.string("?daughter(0).isGlobalMuon?daughter(0).masterClone.globalTrack.hitPattern.numberOfValidMuonHits: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2NofMuonHits"),
    quantity = cms.untracked.string("?daughter(1).isGlobalMuon?daughter(1).masterClone.globalTrack.hitPattern.numberOfValidMuonHits: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1NofStripHits"),
    quantity = cms.untracked.string("?daughter(0).isGlobalMuon?daughter(0).masterClone.globalTrack.hitPattern.numberOfValidStripHits: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2NofStripHits"),
    quantity = cms.untracked.string("?daughter(1).isGlobalMuon?daughter(1).masterClone.globalTrack.hitPattern.numberOfValidStripHits: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1NofPixelHits"),
    quantity = cms.untracked.string("?daughter(0).isGlobalMuon?daughter(0).masterClone.globalTrack.hitPattern.numberOfValidPixelHits: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2NofPixelHits"),
    quantity = cms.untracked.string("?daughter(1).isGlobalMuon?daughter(1).masterClone.globalTrack.hitPattern.numberOfValidPixelHits: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1NormChi2"),
    quantity = cms.untracked.string("?daughter(0).isGlobalMuon?daughter(0).masterClone.globalTrack.normalizedChi2: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2NormChi2"),
    quantity = cms.untracked.string("?daughter(1).isGlobalMuon?daughter(1).masterClone.globalTrack.normalizedChi2: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1NofChambers"),
    quantity = cms.untracked.string("?daughter(0).isGlobalMuon?daughter(0).masterClone.numberOfChambers: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2NofChambers"),
    quantity = cms.untracked.string("?daughter(1).isGlobalMuon?daughter(1).masterClone.numberOfChambers: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1dB"),
    quantity = cms.untracked.string("daughter(0).masterClone.dB")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2dB"),
    quantity = cms.untracked.string("daughter(1).masterClone.dB")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1TrkIso"),
    quantity = cms.untracked.string("daughter(0).masterClone.trackIso")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2TrkIso"),
    quantity = cms.untracked.string("daughter(1).masterClone.trackIso")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1EcalIso"),
    quantity = cms.untracked.string("daughter(0).masterClone.ecalIso")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2EcalIso"),
    quantity = cms.untracked.string("daughter(1).masterClone.ecalIso")
    ),
      cms.PSet(
    tag = cms.untracked.string("LeptDau1HcalIso"),
    quantity = cms.untracked.string("daughter(0).masterClone.hcalIso")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2HcalIso"),
    quantity = cms.untracked.string("daughter(1).masterClone.hcalIso")
    ),
     cms.PSet(
    tag = cms.untracked.string("LeptDau1CombRelIso"),
    quantity = cms.untracked.string("(daughter(0).masterClone.hcalIso + daughter(0).masterClone.ecalIso + daughter(0).masterClone.trackIso )/ daughter(0).pt")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2CombrelIso"),
    quantity = cms.untracked.string("(daughter(1).masterClone.hcalIso + daughter(1).masterClone.ecalIso + daughter(1).masterClone.trackIso )/ daughter(1).pt")
    ),

    )

zee =(

    cms.PSet(
    tag = cms.untracked.string("EleDau1VBTF80CombID"),
    quantity = cms.untracked.string("daughter(0).masterClone.electronID(\"eidVBTFCom80\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("EleDau2VBTF80CombID"),
    quantity = cms.untracked.string("daughter(1).masterClone.electronID(\"eidVBTFCom80\")")
    )
    )

zem =(
    cms.PSet(
    tag = cms.untracked.string("MuonCharge"),
    quantity = cms.untracked.string("daughter(0).charge()")
    ),
    cms.PSet(
    tag = cms.untracked.string("EleCharge"),
    quantity = cms.untracked.string("daughter(1).charge()")
    ),
    cms.PSet(
    tag = cms.untracked.string("EleDau2VBTF80CombID"),
    quantity = cms.untracked.string("daughter(1).masterClone.electronID(\"eidVBTFCom80\")")
    )
    )

j = (
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

    ### b tagging ###
    cms.PSet(
    tag = cms.untracked.string("JetCSV"),
    quantity = cms.untracked.string("bDiscriminator(\"combinedSecondaryVertexBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetTCHE"),
    quantity = cms.untracked.string("bDiscriminator(\"trackCountingHighEffBJetTags\") ")
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
  


## zll.variables += baseKinematics


ZmmEdmNtuple = copy.deepcopy(template)
ZmmEdmNtuple.variables += zll
ZmmEdmNtuple.src = cms.InputTag("zmm")
ZmmEdmNtuple.prefix = cms.untracked.string("zmm")


ZeeEdmNtuple = copy.deepcopy(template)
ZeeEdmNtuple.variables += zll
ZeeEdmNtuple.variables += zee
ZeeEdmNtuple.src = cms.InputTag("zee")
ZeeEdmNtuple.prefix = cms.untracked.string("zee")



ZemEdmNtuple = copy.deepcopy(template)
ZemEdmNtuple.variables += zll
ZemEdmNtuple.variables += zem
ZemEdmNtuple.src = cms.InputTag("zem")
ZemEdmNtuple.prefix = cms.untracked.string("zem")


jEdmNtuple = copy.deepcopy(template)
jEdmNtuple.variables += j
jEdmNtuple.src = cms.InputTag("selectedPatJets")
jEdmNtuple.prefix = cms.untracked.string("j")


### met info
MetEdmNtuple =  cms.EDProducer(
        "CandViewNtpProducer",
        src=cms.InputTag("patMETs"),
        lazyParser=cms.untracked.bool(True),
        prefix=cms.untracked.string("pfMet"),
        eventInfo=cms.untracked.bool(False),
        variables = cms.VPSet(
    cms.PSet(
    tag = cms.untracked.string("Et"),
    quantity = cms.untracked.string("et")
    ),
    cms.PSet(
    tag = cms.untracked.string("MEtSig"),
    quantity = cms.untracked.string("mEtSig"),
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
    "keep *_jEdmNtuple_*_*",
    "keep *_ZmmEdmNtuple_*_*",
    "keep *_ZeeEdmNtuple_*_*",
    "keep *_ZemEdmNtuple_*_*",
    "keep *_MetEdmNtuple_*_*",
    "keep *_JetLepEdmNtuple_*_*",

    
    )
    )
