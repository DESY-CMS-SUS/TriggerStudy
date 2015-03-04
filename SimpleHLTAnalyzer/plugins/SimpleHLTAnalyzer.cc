// -*- C++ -*-
//
// Package:    HLTrigger/SimpleHLTAnalyzer
// Class:      SimpleHLTAnalyzer
// 
/**\class SimpleHLTAnalyzer SimpleHLTAnalyzer.cc HLTrigger/SimpleHLTAnalyzer/plugins/SimpleHLTAnalyzer.cc

 Description: save trigger decisions into a TTree, as well as physics objects for computing efficiencies

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Valentina Dutta
//         Created:  Tue, 03 Mar 2015 10:42:08 GMT
//
//


#include <memory>

#include "TTree.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"


class SimpleHLTAnalyzer : public edm::EDAnalyzer {
  public:
    explicit SimpleHLTAnalyzer(const edm::ParameterSet&);
    ~SimpleHLTAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    virtual void endRun(edm::Run const&, edm::EventSetup const&) override;

    // ----------member data ---------------------------

    HLTConfigProvider hltConfig_;
    std::string hltProcess_;

    edm::InputTag trigFilter_;
    edm::InputTag trigJetFilter_;
    edm::InputTag trigMuonFilterAux_;

    edm::EDGetTokenT<edm::TriggerResults>    trigresultsToken_;
    edm::EDGetTokenT<trigger::TriggerEvent>  trigsummaryToken_;
    edm::EDGetTokenT<reco::GenMETCollection> genmetToken_;
    edm::EDGetTokenT<reco::GenJetCollection> genjetsToken_;

    unsigned int run_;
    unsigned int lumi_;
    unsigned int evt_;

    bool firstEvent_;

    const int maxResults_;

    unsigned int *passtrig_;
    float genmetpt_, genmetphi_;
    int ngenjets_;
    float *genjetpt_, *genjeteta_, *genjetphi_, *genjetmass_;
    float hltmetpt_, hltmetphi_;
    int nhltjets_;
    float *hltjetpt_, *hltjeteta_, *hltjetphi_, *hltjetmass_;
    int nhltmuons_;
    float *hltmuonpt_, *hltmuoneta_, *hltmuonphi_, *hltmuonmass_;

    TTree* outTree_;

};

SimpleHLTAnalyzer::SimpleHLTAnalyzer(const edm::ParameterSet& iConfig) :
  hltProcess_       (iConfig.getParameter<std::string>("hltProcess")),
  trigFilter_       (iConfig.getParameter<edm::InputTag>("trigFilter")),
  trigJetFilter_    (iConfig.getParameter<edm::InputTag>("trigJetFilter")),
  trigMuonFilterAux_(iConfig.getParameter<edm::InputTag>("trigMuonFilterAux")),
  trigresultsToken_ (consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("trigResults"))),
  trigsummaryToken_ (consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("trigSummary"))),
  genmetToken_      (consumes<reco::GenMETCollection>(iConfig.getParameter<edm::InputTag>("genMet"))),
  genjetsToken_     (consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets"))),
  run_  (0),
  lumi_ (0),
  evt_  (0),
  firstEvent_ (true),
  maxResults_ (10000)
{

  passtrig_ = new unsigned int[maxResults_];
  genmetpt_ = 0.0;
  genmetphi_ = 0.0;
  ngenjets_ = 0;
  genjetpt_ = new float[maxResults_];
  genjeteta_ = new float[maxResults_];
  genjetphi_ = new float[maxResults_];
  genjetmass_ = new float[maxResults_];
  hltmetpt_ = 0.0;
  hltmetphi_ = 0.0;
  nhltjets_ = 0;
  hltjetpt_ = new float[maxResults_];
  hltjeteta_ = new float[maxResults_];
  hltjetphi_ = new float[maxResults_];
  hltjetmass_ = new float[maxResults_];
  nhltmuons_ = 0;
  hltmuonpt_ = new float[maxResults_];
  hltmuoneta_ = new float[maxResults_];
  hltmuonphi_ = new float[maxResults_];
  hltmuonmass_ = new float[maxResults_];

  edm::Service<TFileService> fs;
  outTree_ = fs->make<TTree>("HLTAnalysis","");

  outTree_->Branch("run",  &run_,  "run/i");
  outTree_->Branch("lumi", &lumi_, "lumi/i");
  outTree_->Branch("evt",  &evt_,  "evt/i");
  outTree_->Branch("genmet",     &genmetpt_,  "genmet/F");
  outTree_->Branch("genmetphi",  &genmetphi_, "genmetphi/F");
  outTree_->Branch("ngenjets",   &ngenjets_,  "ngenjets/I");
  outTree_->Branch("genjetpt",   genjetpt_,   "genjetpt[ngenjets]/F");
  outTree_->Branch("genjeteta",  genjeteta_,  "genjeteta[ngenjets]/F");
  outTree_->Branch("genjetphi",  genjetphi_,  "genjetphi[ngenjets]/F");
  outTree_->Branch("genjetmass", genjetmass_, "genjetmass[ngenjets]/F");
  outTree_->Branch("hltmet",     &hltmetpt_,  "hltmet/F");
  outTree_->Branch("hltmetphi",  &hltmetphi_, "hltmetphi/F");
  outTree_->Branch("nhltjets",   &nhltjets_,  "nhltjets/I");
  outTree_->Branch("hltjetpt",   hltjetpt_,   "hltjetpt[nhltjets]/F");
  outTree_->Branch("hltjeteta",  hltjeteta_,  "hltjeteta[nhltjets]/F");
  outTree_->Branch("hltjetphi",  hltjetphi_,  "hltjetphi[nhltjets]/F");
  outTree_->Branch("hltjetmass", hltjetmass_, "hltjetmass[nhltjets]/F");
  outTree_->Branch("nhltmuons",   &nhltmuons_,  "nhltmuons/I");
  outTree_->Branch("hltmuonpt",   hltmuonpt_,   "hltmuonpt[nhltmuons]/F");
  outTree_->Branch("hltmuoneta",  hltmuoneta_,  "hltmuoneta[nhltmuons]/F");
  outTree_->Branch("hltmuonphi",  hltmuonphi_,  "hltmuonphi[nhltmuons]/F");
  outTree_->Branch("hltmuonmass", hltmuonmass_, "hltmuonmass[nhltmuons]/F");

}


SimpleHLTAnalyzer::~SimpleHLTAnalyzer()
{}


// ------------ method called for each event  ------------
void
SimpleHLTAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // event header
  run_  = iEvent.id().run();
  lumi_ = iEvent.luminosityBlock();
  evt_  = iEvent.id().event();

  // trigger results
  edm::Handle<edm::TriggerResults> trigresults;
  iEvent.getByToken(trigresultsToken_, trigresults);

  edm::TriggerNames const& triggerNames = iEvent.triggerNames(*trigresults);

  for(unsigned int itrig = 0; itrig < trigresults->size(); ++itrig) {
    TString trigname = triggerNames.triggerName(itrig);
    if(firstEvent_) {
      outTree_->Branch(trigname,&passtrig_[itrig],trigname+"/i");
    }
    bool accept = trigresults->accept(itrig);
    if(accept) passtrig_[itrig] = 1;
    else       passtrig_[itrig] = 0;
  }

  if(firstEvent_) firstEvent_ = false;

  edm::Handle<trigger::TriggerEvent> trigsummary;
  iEvent.getByToken(trigsummaryToken_, trigsummary);

  // filter results ... access trigger objects
  trigger::size_type trigFilterIndex = trigsummary->filterIndex( trigFilter_ );
  trigger::size_type trigJetFilterIndex = trigsummary->filterIndex( trigJetFilter_ );
  trigger::size_type trigMuonFilterAuxIndex = trigsummary->filterIndex( trigMuonFilterAux_ );

  trigger::TriggerObjectCollection triggerObjects = trigsummary->getObjects();

  if(trigFilterIndex < trigsummary->sizeFilters()) {
    const trigger::Keys& trigKeys = trigsummary->filterKeys(trigFilterIndex); 
    for(unsigned int ik = 0; ik < trigKeys.size(); ++ik){ 
      const trigger::TriggerObject& obj = triggerObjects[trigKeys[ik]];
      if(ik == 0) {
        hltmetpt_ = obj.pt();
        hltmetphi_ = obj.phi();
      }
    }
  }

  if(trigJetFilterIndex < trigsummary->sizeFilters()) {
    const trigger::Keys& trigKeys = trigsummary->filterKeys(trigJetFilterIndex); 
    nhltjets_ = trigKeys.size();
    for(unsigned int ik = 0; ik < trigKeys.size(); ++ik){ 
      const trigger::TriggerObject& obj = triggerObjects[trigKeys[ik]];
      hltjetpt_[ik] = obj.pt();
      hltjeteta_[ik] = obj.eta();
      hltjetphi_[ik] = obj.phi();
      hltjetmass_[ik] = obj.mass();
    }
  }

  if(trigMuonFilterAuxIndex < trigsummary->sizeFilters()) {
    const trigger::Keys& trigKeys = trigsummary->filterKeys(trigMuonFilterAuxIndex);
    nhltmuons_ = trigKeys.size();
    for(unsigned int ik = 0; ik < trigKeys.size(); ++ik){ 
      const trigger::TriggerObject& obj = triggerObjects[trigKeys[ik]];
      hltmuonpt_[ik] = obj.pt();
      hltmuoneta_[ik] = obj.eta();
      hltmuonphi_[ik] = obj.phi();
      hltmuonmass_[ik] = obj.mass();
    }
  }

  // gen quantities
  edm::Handle<reco::GenMETCollection> genmets;
  iEvent.getByToken(genmetToken_, genmets);

  const reco::GenMET &genmet = genmets->front();

  genmetpt_  = genmet.pt();
  genmetphi_ = genmet.phi();

  edm::Handle<reco::GenJetCollection> genjets;
  iEvent.getByToken(genjetsToken_, genjets);

  ngenjets_ = genjets->size();
  unsigned int igj = 0;
  for(const reco::GenJet &gj : *genjets) {
    genjetpt_[igj]   = gj.pt();
    genjeteta_[igj]  = gj.eta();
    genjetphi_[igj]  = gj.phi();
    genjetmass_[igj] = gj.mass();
    igj++;
  }

  outTree_->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
SimpleHLTAnalyzer::beginJob()
{}

// ------------ method called once each job just after ending the event loop  ------------
void 
SimpleHLTAnalyzer::endJob() 
{}

// ------------ method called when starting to processes a run  ------------
void 
SimpleHLTAnalyzer::beginRun(edm::Run const &run, edm::EventSetup const &es)
{

  bool changed;

  if (!hltConfig_.init(run, es, hltProcess_, changed)) {
    edm::LogError("SimpleHLTAnalyzer") << "Initialization of HLTConfigProvider failed!!";
    return;
  }

}

// ------------ method called when ending the processing of a run  ------------
void 
SimpleHLTAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SimpleHLTAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimpleHLTAnalyzer);
