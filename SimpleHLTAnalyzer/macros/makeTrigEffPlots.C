#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TMath.h"
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include <vector>
#include "Style.hh"
#endif

float deltaR(float eta1, float phi1, float eta2, float phi2)
{
  float deltaPhi = TMath::Abs(phi1-phi2);
  float deltaEta = eta1-eta2;
  if(deltaPhi > TMath::Pi())
    deltaPhi = TMath::TwoPi() - deltaPhi;
  return TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}

TGraphAsymmErrors* makeEffGraph(TH1F* pass, TH1F* total, bool debug=false) {

  // make sure <pass> and <total> have the same binning!

  int npoints = total->GetNbinsX();

  float x[npoints], y[npoints], errx[npoints], erryl[npoints], erryh[npoints];

  float npass = 0.0;
  float ntotal = 0.0;

  for(int ibin = 1; ibin < npoints+1; ibin++) {
    x[ibin-1] = total->GetBinCenter(ibin);
    npass = pass->GetBinContent(ibin);
    ntotal = total->GetBinContent(ibin);
    y[ibin-1] = ntotal < 1.0 ? 0.0 : npass/ntotal;
    errx[ibin-1] = 0.0;
    if(y[ibin-1]==0.0) {
      erryl[ibin-1] = 0.0; erryh[ibin-1] = 0.0;
    } else {
      if(debug) printf("npass = %3.1f, ntotal = %3.1f, eff = %4.2f", npass, ntotal, y[ibin-1]);
      erryl[ibin-1] = y[ibin-1] - TEfficiency::ClopperPearson((unsigned int)ntotal, (unsigned int)npass, 0.683, false);
      erryh[ibin-1] = TEfficiency::ClopperPearson((unsigned int)ntotal, (unsigned int)npass, 0.683, true) - y[ibin-1];
    }
  }

  TGraphAsymmErrors *gr = new TGraphAsymmErrors(npoints, x, y, errx, errx, erryl, erryh);

  return gr;

}

void runTrigEffPlots(std::vector<TString> infileNames)
{

  TChain *chain = new TChain("hltana/HLTAnalysis");

  // add input file names
  for(const auto fname : infileNames) {
    chain->Add(fname);
  }

  int nevents = chain->GetEntries();

  printf("Will process %d events\n",nevents);

  const int maxsize = 1000;

  // trigger decision
  unsigned int passtrig = 0;  chain->SetBranchAddress("HLT_DiCentralPFJet70_PFMET120_NoiseCleaned_v1", &passtrig);
  // auxiliary trigger decision
  unsigned int passaux = 0;   chain->SetBranchAddress("HLT_IsoMu24_eta2p1_v1", &passaux);
  // GenMET
  float genmet = 0.0;         chain->SetBranchAddress("genmet", &genmet);
  float genmetphi = 0.0;      chain->SetBranchAddress("genmetphi", &genmetphi);
  // GenJets
  int ngenjets = 0;           chain->SetBranchAddress("ngenjets", &ngenjets);
  float genjetpt[maxsize];    chain->SetBranchAddress("genjetpt", &genjetpt);
  float genjeteta[maxsize];   chain->SetBranchAddress("genjeteta", &genjeteta);
  float genjetphi[maxsize];   chain->SetBranchAddress("genjetphi", &genjetphi);
  // Trigger objects (jets)
  int nhltjets = 0;           chain->SetBranchAddress("nhltjets", &nhltjets);
  float hltjetpt[maxsize];    chain->SetBranchAddress("hltjetpt", &hltjetpt);
  float hltjeteta[maxsize];   chain->SetBranchAddress("hltjeteta", &hltjeteta);
  float hltjetphi[maxsize];   chain->SetBranchAddress("hltjetphi", &hltjetphi);
  // Trigger objects (muons)
  int nhltmuons = 0;          chain->SetBranchAddress("nhltmuons", &nhltmuons);
  float hltmuonpt[maxsize];   chain->SetBranchAddress("hltmuonpt", &hltmuonpt);
  float hltmuoneta[maxsize];  chain->SetBranchAddress("hltmuoneta", &hltmuoneta);
  float hltmuonphi[maxsize];  chain->SetBranchAddress("hltmuonphi", &hltmuonphi);

  // use a variable binning to get sufficient statistics everywhere
  int nbins = 136;
  float xbins[nbins+1];
  for(int ibin = 0; ibin < nbins+1; ibin++) {
    if(ibin < 61) {
      xbins[ibin] = 5.0*ibin;
    }
    if(ibin >= 61 && ibin < 91) {
      xbins[ibin] = xbins[60] + 10.0*(ibin-60);
    }
    if(ibin >= 91 && ibin < 121) {
      xbins[ibin] = xbins[90] + 20.0*(ibin-90);
    }
    if(ibin >= 121) {
      xbins[ibin] = xbins[120] + 50.0*(ibin-120);
    }
  }

  // histograms to store denominators and numerators
  // vs MET
  TH1F* hGenMet = new TH1F("hGenMet","",nbins,xbins);                                 // denominator
  TH1F* hGenMet_PassTrig = new TH1F("hGenMet_PassTrig","",nbins,xbins);               // numerator for full trigger
  // vs sub-leading jet pT
  TH1F* hGenJet2Pt = new TH1F("hGenJet2Pt","",nbins,xbins);                                 // denominator
  TH1F* hGenJet2Pt_PassDijetFilter = new TH1F("hGenJet2Pt_PassDijetFilter","",nbins,xbins); // numerator for jet part of trigger

  for(int ievent = 0; ievent < nevents; ievent++) {
    chain->GetEntry(ievent);

    int ncenjets30 = 0, ncenjets70 = 0;
    float jet1pt = 0.0, jet2pt = 0.0;
    for(int ijet = 0; ijet < ngenjets; ijet++) {
      if(fabs(genjeteta[ijet]) < 2.4 && genjetpt[ijet] > 30.0) {
        bool overlap = false;
        // check whether there is overlap with the muon trigger object which fires the auxiliary trigger ... if so, don't count the jet
        for(int imu = 0; imu < nhltmuons; imu++) {
          if(deltaR(hltmuoneta[imu], hltmuonphi[imu], genjeteta[ijet], genjetphi[ijet]) < 0.4) {
            overlap = true;
          }
        }
        if(!overlap) {
          ncenjets30++;
        }
        if(genjetpt[ijet] > jet1pt) {
          jet1pt = genjetpt[ijet];
        }
        if(genjetpt[ijet] > jet2pt && genjetpt[ijet] < jet1pt) {
          jet2pt = genjetpt[ijet];
        }
        if(genjetpt[ijet] > 70.0) {
          if(!overlap) ncenjets70++;
        }
      }
    }

    // check whether the event passed the DiJet filter
    int npassjets = 0;
    for(int ijet = 0; ijet < nhltjets; ijet++) {
      if(hltjetpt[ijet] > 70.0 && fabs(hltjeteta[ijet]) < 2.6) npassjets++;
    }
    bool passdijetfilter = npassjets > 1;

    if(passaux && ncenjets30 > 4) { // require auxiliary trigger + offline selection
      // vs met
      hGenMet->Fill(genmet);
      if(passtrig) hGenMet_PassTrig->Fill(genmet);
      // vs sub-leading jet pT
      hGenJet2Pt->Fill(jet2pt);
      if(passdijetfilter) hGenJet2Pt_PassDijetFilter->Fill(jet2pt);
    }

  }

  // make efficiency graphs
  TGraphAsymmErrors *effvsGenMet = makeEffGraph(hGenMet_PassTrig, hGenMet);
  TGraphAsymmErrors *effDijetvsGenJet2Pt = makeEffGraph(hGenJet2Pt_PassDijetFilter, hGenJet2Pt);

  // make plots
  SetStyle();
  TCanvas *c = MakeCanvas("effcanvas","",600,600);
  InitGraph(effvsGenMet,"","Gen #slash{E}_{T} [GeV]","Efficiency",0.0,800.0,0.0,1.1,kRed);
  effvsGenMet->SetMarkerStyle(kFullCircle);
  effvsGenMet->SetMarkerSize(0.9);
  effvsGenMet->Draw("AP");
  c->SaveAs("eff_trig_vs_genmet.png");
  InitGraph(effDijetvsGenJet2Pt,"","Sub-leading genjet p_{T} [GeV]","Efficiency",0.0,800.0,0.0,1.1,kRed);
  effDijetvsGenJet2Pt->SetMarkerStyle(kFullCircle);
  effDijetvsGenJet2Pt->SetMarkerSize(0.9);
  effDijetvsGenJet2Pt->Draw("AP");
  c->SaveAs("eff_dijetfilter_vs_genjet2pt.png");

}

void makeTrigEffPlots()
{

  std::vector<TString> filelist;

  filelist.push_back("/afs/cern.ch/work/v/vdutta/public/TriggerTutorial_040315/hltanaoutput_ttbar.root");

  runTrigEffPlots(filelist);

}
