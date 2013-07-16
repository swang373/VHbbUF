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
isData = False 
process.GlobalTag.globaltag = 'GR_R_41_V0::All'
if isData==False:
    process.GlobalTag.globaltag = 'START41_V0::All'
    
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.EventContent.EventContent_cff')


import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
#process.hltFilter.HLTPaths = ["HLT_Mu11", "HLT_Ele17_SW_CaloEleId_L1R"]
#HLT_Ele25_CaloIdVT_TrkIdT_CentralJet30_BTagIP_v*
#HLT_DoubleMu7_v*
#HLT_IsoMu17_v*
#HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v3
# Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Source file : To be run on a Full RECO sample
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
"file:/data/1a/degrutto/data/Spring11/ZH_ZToLL_HToBB_M-115_AODSIM_PU_S1_START311_V1G1-v1.root"
#file:/data/1a/degrutto/data/Spring11/ZBB2JetsToLNu_TuneZ2_7TeV-alpgen-tauola_AODSIM_PU_S1_START311_V1G1-v1.root"
#"file:/data/1a/degrutto/data/Spring11/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_AODSIM_PU_S1_START311_V1G1-v1.root"
###ZH_ZToLL_HToBB_M-115_AODSIM_PU_S1_START311_V1G1-v1.root



    )
)




if isData:
    removeMCMatching(process, ['All'])

##removeAllPATObjectsBut(process, ['METs', 'Jets', 'Electrons', 'Muons'])


### sequence to use PFnoPU
#Create good primary vertices to be used for PF association
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

process.goodOfflinePrimaryVertices = cms.EDFilter(
        "PrimaryVertexObjectFilter",
            filterParams = pvSelector.clone( minNdof = cms.double(7.0), maxZ = cms.double(24.0) ),
            src=cms.InputTag('offlinePrimaryVertices')
            )

#Create the "top-down projection" for the PF2PAT sequence
from PhysicsTools.PatAlgos.tools.pfTools import *
postfix = "PFlow"
if isData:
#    usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=not options.useData, postfix=postfix)
    usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=True, postfix=postfix)
else:
    usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=True, postfix=postfix)
process.pfPileUpPFlow.Enable = True
process.pfPileUpPFlow.Vertices = 'goodOfflinePrimaryVertices'
process.pfJetsPFlow.doAreaFastjet = True

process.pfJetsPFlow.doRhoFastjet = False

process.pfJetsPFlow.Rho_EtaMax = cms.double(5.0)
process.pfJetsPFlow.Ghost_EtaMax = cms.double(4.5)


process.patJetCorrFactorsPFlow.rho = cms.InputTag("kt6PFJetsPFlow", "rho")


###NOTA BENE! In 4.2.x you also have to enable this switch, but not in 4.1.x:

process.pfPileUpPFlow.checkClosestZVertex = cms.bool(False)

        
#Compute the mean pt per unit area ("rho") using KT6 Jets with the Voronoi areas method (which is equivalent to the passive area for KT ONLY).

from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
process.kt6PFJetsPFlow = kt4PFJets.clone(
        rParam = cms.double(0.6),
            src = cms.InputTag('pfNoElectron'+postfix),
            doAreaFastjet = cms.bool(True),
            doRhoFastjet = cms.bool(True),
            voronoiRfact = cms.double(0.9)
           ,Rho_EtaMax = cms.double(5.0)
            ,Ghost_EtaMax = cms.double(4.5)
            )

###needed for efficeinct ghost subtruction
process.kt6PFJetsPFlow.jetPtMin      = 0.5

getattr(process,"patPF2PATSequence"+postfix).replace( getattr(process,"pfNoElectron"+postfix), getattr(process,"pfNoElectron"+postfix)*process.kt6PFJetsPFlow )


process.patJetCorrFactors.rho = cms.InputTag('kt6PFJetsPFlow','rho')





# Muon Selection
process.selectedPatMuons.cut = (
        "pt > 20 && isGlobalMuon && isTrackerMuon && globalTrack().normalizedChi2 < 10 &&" +
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
    "pt > 20.0 && abs(eta) < 2.5 &&"                               +
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

process.selectedPatMET = cms.EDFilter("PATMETSelector",
     src = cms.InputTag("patMETs"),
     cut = cms.string("et > 20 && mEtSig>1 ")
 )



process.METFilter = cms.EDFilter("CandViewCountFilter",
                                     src = cms.InputTag("selectedPatMET"),
                                     minNumber = cms.uint32(1)

                                 )



# Add PF jets



### 41X JEC
inputJetCorrLabel = ('AK5PFchs', ['L1FastJet','L2Relative', 'L3Absolute', 'L2L3Residual'])
if isData==False:
      inputJetCorrLabel = ('AK5PFchs', ['L1FastJet','L2Relative', 'L3Absolute' ])

      
# Switch to using PFJets
switchJetCollection(process,cms.InputTag('ak5PFJets'),
    doJTA        = True,
    doBTagging   = True,
    jetCorrLabel = inputJetCorrLabel,
    doType1MET   = True,
    genJetCollection=cms.InputTag(""),
    doJetID      = True
)

if isData!=True: 
    switchJetCollection(process,cms.InputTag('ak5PFJets'),
                            doJTA        = True,
                            doBTagging   = True,
                            jetCorrLabel = inputJetCorrLabel,
                            doType1MET   = True,
                            genJetCollection=cms.InputTag("ak5GenJets"),
                            doJetID      = True
                        )
    




process.patJets.addTagInfos  = True # ### needed for secondary vertex variables....
process.patJets.tagInfoSources  = cms.VInputTag(
        cms.InputTag("secondaryVertexTagInfosAOD"),
        )

    
### ntuple
process.load("HSkim.VHSkim.VHEdmNtuples_cff")



##  NSVtexJets 



process.SVtxEdmNtuple= cms.EDProducer(
    "SVtxInfoDumper",
    hjj = cms.InputTag("hjj"),
    )

if (isData!=0):
    process.BMCTruthEdmNtuple= cms.EDProducer(
    "BMCTruthDumper",
      ptPairCut = cms.untracked.double(5.)
    )


# Clean the Jets from the selected leptons, applying here only pt cut pf 20 
process.cleanPatJets = cms.EDProducer("PATJetCleaner",
                                      src = cms.InputTag("patJets"),
                                     # preselection = cms.string('pt > 20.0  && ( (neutralEmEnergy/energy < 0.99) &&  (neutralHadronEnergy/energy < 0.99) && numberOfDaughters>1) '),
                                         preselection = cms.string('abs(eta)< 5.0 && pt > 20.0 & (neutralEmEnergy/energy < 0.99) &&  (neutralHadronEnergy/energy < 0.99) && numberOfDaughters>1'),
                                      checkOverlaps = cms.PSet(
  ),
                                      finalCut = cms.string('')
                                      ##finalCut = cms.string('bDiscriminator(\"trackCountingHighEffBJetTags\")>1.7 ')
                                      )
                                      ### select pat jets, applying eta loose ID and  b-tag
process.myPatJets = cms.EDFilter("PATJetSelector",
                                 src= cms.InputTag("cleanPatJets"),
                                 cut = cms.string("abs(eta)< 2.5 & pt > 20.0 ") ##  && bDiscriminator(\"trackCountingHighEffBJetTags\")>1.7 " ) ### 
                                 )
                                      

process.jetFilter = cms.EDFilter("CandViewCountFilter",
                                     src = cms.InputTag("myPatJets"),
                                     minNumber = cms.uint32(2)

                                 )


# Z Candidates and Higgs Candidates
process.zee = cms.EDProducer("CandViewShallowCloneCombiner",
    checkCharge = cms.bool(True),
    cut = cms.string('mass > 20'),
    decay = cms.string("selectedPatElectrons@+ selectedPatElectrons@-")
)

process.zmm = cms.EDProducer("CandViewShallowCloneCombiner",
    checkCharge = cms.bool(True),
    cut = cms.string('mass > 20'),
    decay = cms.string("selectedPatMuons@+ selectedPatMuons@-")
)

# Setup for a basic filtering
process.zll = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag("zee", "zmm")
) 

process.zllFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("zll"),
    minNumber = cms.uint32(1),
)




process.hjj = cms.EDProducer("CandViewCombiner",
    checkCharge = cms.bool(False),
    cut = cms.string('mass>0'),
    decay = cms.string("myPatJets myPatJets")
)   

process.wmn = cms.EDProducer("CandViewShallowCloneCombiner",
                            checkCharge = cms.bool(False),
                            cut = cms.string('sqrt((daughter(0).pt+daughter(1).pt)*(daughter(0).pt+daughter(1).pt)-pt*pt)>10'),
                            decay = cms.string("selectedPatMuons selectedPatMET")
                            )

process.wen = cms.EDProducer("CandViewShallowCloneCombiner",
                            checkCharge = cms.bool(False),
                            cut = cms.string('sqrt((daughter(0).pt+daughter(1).pt)*(daughter(0).pt+daughter(1).pt)-pt*pt)>10'),
                            decay = cms.string("selectedPatElectrons selectedPatMET")
                            )

# W filter
process.wln = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag("wen", "wmn")
) 

process.wlnFilter = cms.EDFilter("CandViewCountFilter",
                           src = cms.InputTag("wln"),
                           minNumber = cms.uint32(1)
                           )



process.hjjzll = cms.EDProducer("CandViewCombiner",
    checkCharge = cms.bool(False),
    cut = cms.string(''),
    decay = cms.string("hjj zll")
)   

process.hjjmet = cms.EDProducer("CandViewCombiner",
                            checkCharge = cms.bool(False),
                            cut = cms.string(''),
                            decay = cms.string("hjj selectedPatMET")
                            )


process.SVtxEdmNtuple= cms.EDProducer(
    "SVtxInfoDumper",
    hjj = cms.InputTag("hjj"),
    )

process.JetLepEdmNtuple= cms.EDProducer(
    "JetLepInfoDumper",
    j = cms.InputTag("cleanPatJets"), ### looking at jets after the cleaning...but before the selections (so no eta....), but with corrected energy!!!!!
    m = cms.InputTag("selectedPatMuons"), 
    e = cms.InputTag("selectedPatElectrons"), 
    )


#process.patseq = cms.Sequence(
#        process.goodOfflinePrimaryVertices*
#            getattr(process,"patPF2PATSequence"+postfix)
#            )

process.ZllHbbPath = cms.Path(
  #   process.hltFilter +
     process.eidSequence +
     process.goodOfflinePrimaryVertices +
     getattr(process,"patPF2PATSequence"+postfix)*
     process.patDefaultSequence*
     process.zmm +
     process.zee +
     process.zll +
     process.zllFilter + 
     process.cleanPatJets +
     process.myPatJets +  
 #    process.jetFilter +
     process.hjj +
     process.hjjzll +
     process.SVtxEdmNtuple +
     process.HjjEdmNtuple +
     process.ZmmEdmNtuple +
     process.ZeeEdmNtuple+
     process.HjjZllEdmNtuple +
     process.MetEdmNtuple +
     process.JetLepEdmNtuple 
 )
if (isData!=0):
    process.ZllHbbPath.__iadd__(process.BMCTruthEdmNtuple)



process.ZinvHbbPath = cms.Path(
#    process.hltFilter +
     process.eidSequence +
     process.goodOfflinePrimaryVertices +
     getattr(process,"patPF2PATSequence"+postfix)*
     process.patDefaultSequence*
     process.selectedPatMET +
     process.METFilter +
     process.cleanPatJets +
     process.myPatJets +
#     process.jetFilter +
     process.hjj *
     process.HjjEdmNtuple *
    process.SVtxEdmNtuple *
     process.MetEdmNtuple *
     process.JetLepEdmNtuple *
     process.hjjmet *
     process.HjjMetEdmNtuple  
 )


process.WlnHbbPath = cms.Path(
    # process.hltFilter +
     process.eidSequence +
     process.goodOfflinePrimaryVertices *
     getattr(process,"patPF2PATSequence"+postfix)*
     process.patDefaultSequence*
     process.selectedPatMET +
     process.METFilter +
     process.wmn +
     process.wen +
     process.wln +
     process.wlnFilter + 
     process.cleanPatJets +
     process.myPatJets +
  #   process.jetFilter +
     process.hjj *
     process.SVtxEdmNtuple *
     process.HjjEdmNtuple *
     process.WmnEdmNtuple *
     process.WenEdmNtuple*
     process.MetEdmNtuple *
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
    'keep *_hjj_*_*',
    'keep *_wmn_*_*',
    'keep *_wen_*_*',
    

    )
)

# to add the AOD output uncomment the following
process.MuEleMetBbbarSkimEventContent.outputCommands.extend(process.AODEventContent.outputCommands)


process.out = cms.OutputModule("PoolOutputModule",
    process.MuEleMetBbbarSkimEventContent,
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('ZllHbbPath', 'ZinvHbbPath',
                                   'WlnHbbPath')
        ),
                                                       
       dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('muelemetbbarSkim'),
        dataTier = cms.untracked.string('USER')
   ),
   dropMetaData = cms.untracked.string("ALL"),                            
   fileName = cms.untracked.string('MuEleMetBbbarSkim.root')
)







process.edmNtuplesOut.fileName = cms.untracked.string('EdmNtuples_Zll.root')
process.edmNtuplesOut.SelectEvents = cms.untracked.PSet(
            SelectEvents = cms.vstring('ZllHbbPath', 'ZinvHbbPath',
                                   'WlnHbbPath',
                                       )
            )                            

process.edmNtuplesOut.dropMetaData = cms.untracked.string("ALL")

process.endPath = cms.EndPath(process.out +
                              process.edmNtuplesOut)


