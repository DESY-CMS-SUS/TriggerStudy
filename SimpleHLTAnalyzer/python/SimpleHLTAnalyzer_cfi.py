import FWCore.ParameterSet.Config as cms

hltana = cms.EDAnalyzer('SimpleHLTAnalyzer',
    hltProcess = cms.string('TEST'),
    trigFilter = cms.InputTag('hltPFMET70Filter', '', 'TEST'),
    trigJetFilter = cms.InputTag('hltPreIsoMu24eta2p1TriCentralPFJet30', '', 'TEST'),
    trigMuonFilterAux = cms.InputTag('hltL3fL1sMu20Eta2p1L1f0L2f10QL3f24QL3pfecalIsoRhoFilteredEB0p13EE0p10','','TEST'),
    trigElectronFilterAux = cms.InputTag('L1SingleEG40ORL1SingleIsoEG30er','TEST'),
    trigSummary = cms.InputTag('hltTriggerSummaryAOD', '', 'TEST'),
    trigResults = cms.InputTag('TriggerResults::TEST'),
    genMet = cms.InputTag('genMetTrue::HLT'),
    genJets = cms.InputTag('ak4GenJets::HLT'),
    genPart = cms.InputTag('genParticles::HLT'),
)
