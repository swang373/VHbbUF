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
    tag = cms.untracked.string("LeptDau1aodCombRelIso"),
    quantity = cms.untracked.string("(daughter(0).masterClone.hcalIso + daughter(0).masterClone.ecalIso + daughter(0).masterClone.trackIso )/ daughter(0).pt")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2aodCombrelIso"),
    quantity = cms.untracked.string("(daughter(1).masterClone.hcalIso + daughter(1).masterClone.ecalIso + daughter(1).masterClone.trackIso )/ daughter(1).pt")
    ),
     cms.PSet(
    tag = cms.untracked.string("LeptDau1pfCombRelIso"),
    quantity = cms.untracked.string("(daughter(0).masterClone.chargedHadronIso + daughter(0).masterClone.neutralHadronIso + daughter(0).masterClone.photonIso )/ daughter(0).pt")
    ),
     cms.PSet(
    tag = cms.untracked.string("LeptDau2pfCombRelIso"),
    quantity = cms.untracked.string("(daughter(1).masterClone.chargedHadronIso + daughter(1).masterClone.neutralHadronIso + daughter(1).masterClone.photonIso )/ daughter(1).pt")
    ),

    )
    

zmm = (  
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
    tag = cms.untracked.string("Lept1MuEnergyEm"),
    quantity = cms.untracked.string("daughter(0).masterClone.calEnergy().em")
    ),
    cms.PSet(
    tag = cms.untracked.string("Lept2MuEnergyEm"),
    quantity = cms.untracked.string("daughter(1).masterClone.calEnergy().em")
    ),
    cms.PSet(
    tag = cms.untracked.string("Lept1MuEnergyHad"),
    quantity = cms.untracked.string("daughter(0).masterClone.calEnergy().had")
    ),
    cms.PSet(
    tag = cms.untracked.string("Lept2MuEnergyHad"),
    quantity = cms.untracked.string("daughter(1).masterClone.calEnergy().had")
    ),

    )

zee =(
##### we can more stuff for electrons
  #  cms.PSet(
  #  tag = cms.untracked.string("EleDau1VBTF80CombID"),
  #  quantity = cms.untracked.string("daughter(0).masterClone.electronID(\"eidVBTFCom80\")")
  #  ),
  #  cms.PSet(
  #  tag = cms.untracked.string("EleDau2VBTF80CombID"),
  #  quantity = cms.untracked.string("daughter(1).masterClone.electronID(\"eidVBTFCom80\")")
  #  )
    )

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
    ### MC parton flavour ###
    cms.PSet(
    tag = cms.untracked.string("Jet1PartonFlavour"),
    quantity = cms.untracked.string("daughter(0).partonFlavour")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2PartonFlavour"),
    quantity = cms.untracked.string("daughter(1).partonFlavour")
    ),
    ### MC genJet ###
    ## cms.PSet(
##     tag = cms.untracked.string("Jet1genJetPt"),
##     quantity = cms.untracked.string("daughter(0).genJet.pt")
##     ),
##     cms.PSet(
##     tag = cms.untracked.string("Jet1genJetEta"),
##     quantity = cms.untracked.string("daughter(0).genJet.eta")
##     ),
##     cms.PSet(
##     tag = cms.untracked.string("Jet1genJetPhi"),
##     quantity = cms.untracked.string("daughter(0).genJet.phi")
##     ),
##     cms.PSet(
##     tag = cms.untracked.string("Jet2genJetPt"),
##     quantity = cms.untracked.string("daughter(1).genJet.pt")
##     ),
##     cms.PSet(
##     tag = cms.untracked.string("Jet2genJetEta"),
##     quantity = cms.untracked.string("daughter(1).genJet.eta")
##     ),
##     cms.PSet(
##     tag = cms.untracked.string("Jet2genJetPhi"),
##     quantity = cms.untracked.string("daughter(1).genJet.phi")
##     ),



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
ZmmEdmNtuple.variables += zmm
ZmmEdmNtuple.src = cms.InputTag("zmm")
ZmmEdmNtuple.prefix = cms.untracked.string("zmm")

ZeeEdmNtuple = copy.deepcopy(template)
ZeeEdmNtuple.variables += zll
ZeeEdmNtuple.variables += zee
ZeeEdmNtuple.src = cms.InputTag("zee")
ZeeEdmNtuple.prefix = cms.untracked.string("zee")


HjjEdmNtuple = copy.deepcopy(template)
HjjEdmNtuple.variables += hjj
HjjEdmNtuple.src = cms.InputTag("hjj")
HjjEdmNtuple.prefix = cms.untracked.string("hjj")


### met info
MetEdmNtuple =  cms.EDProducer(
        "CandViewNtpProducer",
        src=cms.InputTag("patMETs"),
        lazyParser=cms.untracked.bool(True),
        prefix=cms.untracked.string("pfMet"),
        eventInfo=cms.untracked.bool(True),
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
    tag = cms.untracked.string("MEtSignificance"),
    quantity = cms.untracked.string("significance"),
    ),
     cms.PSet(
    tag = cms.untracked.string("Phi"),
    quantity = cms.untracked.string("phi"),

    )
        )
        )


### met info
JetEdmNtuple =  cms.EDProducer(
        "CandViewNtpProducer",
        src=cms.InputTag("cleanPatJets"),
        lazyParser=cms.untracked.bool(True),
        prefix=cms.untracked.string("jets"),
        eventInfo=cms.untracked.bool(False),
        variables = cms.VPSet(
    cms.PSet(
    tag = cms.untracked.string("Pt"),
    quantity = cms.untracked.string("pt")
    ),
    cms.PSet(
    tag = cms.untracked.string("Eta"),
    quantity = cms.untracked.string("eta"),
    ),
    cms.PSet(
    tag = cms.untracked.string("Phi"),
    quantity = cms.untracked.string("phi"),
    ),
    cms.PSet(
    tag = cms.untracked.string("Mass"),
    quantity = cms.untracked.string("mass")
    ),
    cms.PSet(
    tag = cms.untracked.string("flavour"),
    quantity = cms.untracked.string("partonFlavour")
    ),
#### Mc genJet
  ##   cms.PSet(
##     tag = cms.untracked.string("genJetPt"),
##     quantity = cms.untracked.string("genJet().pt")
##     ),
##     cms.PSet(
##     tag = cms.untracked.string("genJetEta"),
##     quantity = cms.untracked.string("genJet().eta")
##     ),
##     cms.PSet(
##     tag = cms.untracked.string("genJetPhi"),
##     quantity = cms.untracked.string("genJet().phi")
##     ),

    
    ### b tagging ###
    cms.PSet(
    tag = cms.untracked.string("CSV"),
    quantity = cms.untracked.string("bDiscriminator(\"combinedSecondaryVertexBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("CSVMVA"),
    quantity = cms.untracked.string("bDiscriminator(\"combinedSecondaryVertexMVABJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("JProb"),
    quantity = cms.untracked.string("bDiscriminator(\"jetProbabilityBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("JbProb"),
    quantity = cms.untracked.string("bDiscriminator(\"jetBProbabilityBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("SSVHE"),
    quantity = cms.untracked.string("bDiscriminator(\"simpleSecondaryVertexHighEffBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("SSVHP"),
    quantity = cms.untracked.string("bDiscriminator(\"simpleSecondaryVertexHighPurBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("ElPt"),
    quantity = cms.untracked.string("bDiscriminator(\"softElectronByPtBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("ElIp"),
    quantity = cms.untracked.string("bDiscriminator(\"softElectronByIP3dBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("Mu"),
    quantity = cms.untracked.string("bDiscriminator(\"softMuonBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("MuPt"),
    quantity = cms.untracked.string("bDiscriminator(\"softMuonByPtBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("MuIp"),
    quantity = cms.untracked.string("bDiscriminator(\"softMuonByIP3dBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("TKHE"),
    quantity = cms.untracked.string("bDiscriminator(\"trackCountingHighEffBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("TKHP"),
    quantity = cms.untracked.string("bDiscriminator(\"trackCountingHighPurBJetTags\") ")
    ),
      ### JET ID
  ##       chf = patJet->chargedHadronEnergyFraction();
##         nhf = ( patJet->neutralHadronEnergy() + patJet->HFHadronEnergy() ) / patJet->energy();
##         cef = patJet->chargedEmEnergyFraction();
##         nef = patJet->neutralEmEnergyFraction();
##         nch = patJet->chargedMultiplicity();
##         nconstituents = patJet->numberOfDaughters();

    cms.PSet(
    tag = cms.untracked.string("chf"),
    quantity = cms.untracked.string("chargedHadronEnergyFraction")
    ),
    cms.PSet(
    tag = cms.untracked.string("nhf"),
    quantity = cms.untracked.string("(neutralHadronEnergy() + HFHadronEnergy() ) / energy()")
    ),
    cms.PSet(
    tag = cms.untracked.string("cef"),
    quantity = cms.untracked.string("chargedEmEnergyFraction")
    ),
    cms.PSet(
    tag = cms.untracked.string("nef"),
    quantity = cms.untracked.string("neutralEmEnergyFraction")
    ),
    cms.PSet(
    tag = cms.untracked.string("nch"),
    quantity = cms.untracked.string("chargedMultiplicity")
    ),
    cms.PSet(
    tag = cms.untracked.string("nconstituents"),
    quantity = cms.untracked.string("numberOfDaughters")
    ),
    
    )
        )



### met info
MuEdmNtuple =  cms.EDProducer(
        "CandViewNtpProducer",
        src=cms.InputTag("selectedPatMuons"),
        lazyParser=cms.untracked.bool(True),
        prefix=cms.untracked.string("mu"),
        eventInfo=cms.untracked.bool(False),
        variables = cms.VPSet(
    cms.PSet(
    tag = cms.untracked.string("Pt"),
    quantity = cms.untracked.string("pt")
    ),
    cms.PSet(
    tag = cms.untracked.string("Eta"),
    quantity = cms.untracked.string("eta"),
    ),
    cms.PSet(
    tag = cms.untracked.string("Phi"),
    quantity = cms.untracked.string("phi"),
    ),
    cms.PSet(
    tag = cms.untracked.string("aodcombIso"),
    quantity = cms.untracked.string("hcalIso + ecalIso + trackIso")
    ),
    cms.PSet(
    tag = cms.untracked.string("particleIso"),
    quantity = cms.untracked.string("particleIso")
    ),
     cms.PSet(
    tag = cms.untracked.string("chargeHadronIso"),
    quantity = cms.untracked.string("chargedHadronIso")
    ),
     cms.PSet(
    tag = cms.untracked.string("neutralHadronIso"),
    quantity = cms.untracked.string("neutralHadronIso")
    ),
    cms.PSet(
    tag = cms.untracked.string("photonIso"),
    quantity = cms.untracked.string("photonIso")
    ),
)
)

ElecEdmNtuple = copy.deepcopy(MuEdmNtuple)
ElecEdmNtuple.src= cms.InputTag("selectedPatElectrons")
ElecEdmNtuple.prefix=cms.untracked.string("elec")



HjjZllEdmNtuple = cms.EDProducer(
    "CandViewNtpProducer",
    src=cms.InputTag("hjjzll"),
    lazyParser=cms.untracked.bool(True),
    prefix=cms.untracked.string("hjjzll"),
    eventInfo=cms.untracked.bool(False),
    variables = cms.VPSet(
    cms.PSet(
    tag = cms.untracked.string("Mass"),
    quantity = cms.untracked.string("mass")
    ),
    cms.PSet(
    tag = cms.untracked.string("Pt"),
    quantity = cms.untracked.string("pt"),
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




WmnEdmNtuple =  cms.EDProducer(
    "CandViewNtpProducer",
    src=cms.InputTag("wmn"),
    lazyParser=cms.untracked.bool(True),
    prefix=cms.untracked.string("wmn"),
    eventInfo=cms.untracked.bool(False),
    variables = cms.VPSet(
     cms.PSet(
    tag = cms.untracked.string("MT"),
    quantity = cms.untracked.string("sqrt((daughter(0).pt+daughter(1).pt)*(daughter(0).pt+daughter(1).pt)-pt*pt)")
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
    ## eta
    cms.PSet(
    tag = cms.untracked.string("LeptDau1Eta"),
    quantity = cms.untracked.string("daughter(0).eta ")
    ),
    ## phi
    cms.PSet(
    tag = cms.untracked.string("LeptDau1Phi"),
    quantity = cms.untracked.string("daughter(0).phi ")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1GlobalMuonBit"),
    quantity = cms.untracked.string("daughter(0).isGlobalMuon")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1StandAloneBit"),
    quantity = cms.untracked.string("daughter(0).isStandAloneMuon")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1TrackerMuonBit"),
    quantity = cms.untracked.string("daughter(0).isTrackerMuon")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1NofMuonHits"),
    quantity = cms.untracked.string("?daughter(0).isGlobalMuon?daughter(0).masterClone.globalTrack.hitPattern.numberOfValidMuonHits: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1NofStripHits"),
    quantity = cms.untracked.string("?daughter(0).isGlobalMuon?daughter(0).masterClone.globalTrack.hitPattern.numberOfValidStripHits: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1NofPixelHits"),
    quantity = cms.untracked.string("?daughter(0).isGlobalMuon?daughter(0).masterClone.globalTrack.hitPattern.numberOfValidPixelHits: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1NormChi2"),
    quantity = cms.untracked.string("?daughter(0).isGlobalMuon?daughter(0).masterClone.globalTrack.normalizedChi2: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1NofChambers"),
    quantity = cms.untracked.string("?daughter(0).isGlobalMuon?daughter(0).masterClone.numberOfChambers: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1dB"),
    quantity = cms.untracked.string("daughter(0).masterClone.dB")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1TrkIso"),
    quantity = cms.untracked.string("daughter(0).masterClone.trackIso")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1EcalIso"),
    quantity = cms.untracked.string("daughter(0).masterClone.ecalIso")
    ),
      cms.PSet(
    tag = cms.untracked.string("LeptDau1HcalIso"),
    quantity = cms.untracked.string("daughter(0).masterClone.hcalIso")
    ),
     cms.PSet(
    tag = cms.untracked.string("LeptDau1CombRelIso"),
    quantity = cms.untracked.string("(daughter(0).masterClone.hcalIso + daughter(0).masterClone.ecalIso + daughter(0).masterClone.trackIso )/ daughter(0).pt")
    ),
     cms.PSet(
    tag = cms.untracked.string("MEtEt"),
    quantity = cms.untracked.string("daughter(1).masterClone.et")
    ),
    cms.PSet(
    tag = cms.untracked.string("MEtSig"),
    quantity = cms.untracked.string("daughter(1).masterClone.mEtSig"),
    ),
     cms.PSet(
    tag = cms.untracked.string("MEtPhi"),
    quantity = cms.untracked.string("daughter(1).masterClone.phi"),
    
    )
    )
    )
        
wen =(
###we can add some for electrons
   # cms.PSet(
   # tag = cms.untracked.string("EleDau1VBTF80CombID"),
   # quantity = cms.untracked.string("daughter(0).masterClone.electronID(\"eidVBTFCom80\")")
   # ),

    )

WenEdmNtuple = copy.deepcopy(WmnEdmNtuple)
WenEdmNtuple.variables += wen
WenEdmNtuple.src= cms.InputTag("wen")
WenEdmNtuple.prefix=cms.untracked.string("wen")


#### now adding also vertext constrain info ####
#hjjVtxed = cms.EDProducer(
#    "KalmanVertexFitCompositeCandProducer",
#    src = cms.InputTag("hjj")
#)


#hjjzllVtxed = cms.EDProducer(
#    "KalmanVertexFitCompositeCandProducer",
#    src = cms.InputTag("hjjzll")
#)

#hjjVtxedEdmNtuple = cms.EDProducer(
#    "CandViewNtpProducer",
#    src = cms.InputTag("hjjVtxed"),
#    lazyParser=cms.untracked.bool(True),
#    prefix=cms.untracked.string("hjjVtxed"),
#    variables = cms.VPSet(
#      cms.PSet(
#        tag = cms.untracked.string("mass"),
#        quantity = cms.untracked.string("mass")
#      ),
#       cms.PSet(
#        tag = cms.untracked.string("vertexNdof"),
#        quantity = cms.untracked.string("vertexNdof")
#      ),
#       cms.PSet(
#        tag = cms.untracked.string("vertexNormalizedChi2"),
#        quantity = cms.untracked.string("vertexNormalizedChi2")
#      ),
#    )
#)

#hjjzllVtxedEdmNtuple = copy.deepcopy(hjjVtxedEdmNtuple)
#hjjzllVtxedEdmNtuple.src= cms.InputTag("hjjzllVtxed")
#hjjzllVtxedEdmNtuple.prefix=cms.untracked.string("hjjzllVtxed")



       
edmNtuplesOut= cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('2l2bMetEdmNtuples.root'),
    outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_HjjEdmNtuple_*_*",
    "keep *_ZmmEdmNtuple_*_*",
    "keep *_ZeeEdmNtuple_*_*",
    "keep *_WmnEdmNtuple_*_*",
    "keep *_WenEdmNtuple_*_*",
    "keep *_HjjZllEdmNtuple_*_*",
    "keep *_HjjMetEdmNtuple_*_*",
    "keep *_MetEdmNtuple_*_*",
    "keep *_JetLepEdmNtuple_*_*",
    "keep *_JetEdmNtuple_*_*",
    "keep *_MuEdmNtuple_*_*",
    "keep *_ElecEdmNtuple_*_*",
 #   "keep *_hjjVtxedEdmNtuple_*_*",
 #   "keep *_hjjzllVtxedEdmNtuple_*_*",
    "keep *_SVtxEdmNtuple_*_*",
    "keep *_BMCTruthEdmNtuple_*_*",
    "keep *_HLTInfoMetJetEdmNtuple_*_*"
    )
    )
