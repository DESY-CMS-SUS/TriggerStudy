import FWCore.ParameterSet.Config as cms

hltana = cms.EDAnalyzer('SimpleHLTAnalyzer',
    hltProcess = cms.string('TEST'),
    trigFilter = cms.InputTag('hltPFMET120Filter', '', 'TEST'),
    trigJetFilter = cms.InputTag('hltDiCentralPFJet70', '', 'TEST'),
    trigMuonFilterAux = cms.InputTag('hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09','','TEST'),
    trigSummary = cms.InputTag('hltTriggerSummaryAOD', '', 'TEST'),
    trigResults = cms.InputTag('TriggerResults::TEST'),
    genMet = cms.InputTag('genMetTrue::HLT'),
    genJets = cms.InputTag('ak4GenJets::HLT'),
)
