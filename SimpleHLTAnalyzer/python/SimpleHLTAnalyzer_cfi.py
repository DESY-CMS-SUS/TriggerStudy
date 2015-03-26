import FWCore.ParameterSet.Config as cms

hltana = cms.EDAnalyzer('SimpleHLTAnalyzer',
    hltProcess = cms.string('TEST'),
    trigFilter = cms.InputTag('hltPFMET70Filter', '', 'TEST'),
    trigJetFilter = cms.InputTag('hltPreIsoMu24eta2p1TriCentralPFJet30', '', 'TEST'),
    trigMuonFilterAux = cms.InputTag('hltL2fL1sMu5L1f0L2Filtered0Q','','TEST'),
    trigElectronFilterAux = cms.InputTag('hltEle15VVVLGsfTrackIsoFilter','TEST'),
    trigSummary = cms.InputTag('hltTriggerSummaryAOD', '', 'TEST'),
    trigResults = cms.InputTag('TriggerResults::TEST'),
    genMet = cms.InputTag('genMetTrue::HLT'),
    genJets = cms.InputTag('ak4GenJets::HLT'),
    genPart = cms.InputTag('genParticles::HLT'),
)
