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
  unsigned int passtrig = 0;  chain->SetBranchAddress("HLT_PFHT350_PFMET120_NoiseCleaned_v1", &passtrig);
  unsigned int passtrig0 = 0;  chain->SetBranchAddress("HLT_PFHT900_v1", &passtrig0);
  unsigned int passtrig1 = 0;  chain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT600_v1", &passtrig1);
  unsigned int passtrig2 = 0;  chain->SetBranchAddress("HLT_IsoMu24_eta2p1_TriCentralPFJet30_v1", &passtrig2);
  unsigned int passtrig3 = 0;  chain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT400_PFMET70_v1", &passtrig3);
  unsigned int passtrig4 = 0;  chain->SetBranchAddress("HLT_PFMET120_NoiseCleaned_Mu5_v1", &passtrig4);
  unsigned int passtrig5 = 0;  chain->SetBranchAddress("HLT_IsoMu20_eta2p1_TriCentralPFJet50_40_30_v1", &passtrig5);

  // auxiliary trigger decision
  unsigned int passaux = 0;   chain->SetBranchAddress("HLT_Mu45_eta2p1_v1", &passaux);
  //unsigned int passaux = 0;   chain->SetBranchAddress("HLT_IsoMu27_v1", &passaux);
  //unsigned int passaux = 0;   chain->SetBranchAddress("HLT_PFHT900_v1", &passaux);
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

  int nmuons = 0;          chain->SetBranchAddress("nmuons", &nmuons);
  float muonpt[maxsize];   chain->SetBranchAddress("muonpt", &muonpt);
  float muoneta[maxsize];  chain->SetBranchAddress("muoneta", &muoneta);
  float muonphi[maxsize];  chain->SetBranchAddress("muonphi", &muonphi);


  // use a variable binning to get sufficient statistics everywhere
  int nbins = 136;
  float xbins[nbins+1];
  for(int ibin = 0; ibin < nbins+1; ibin++) {
    if(ibin < 61) {
      xbins[ibin] = 40.0*ibin;
    }
    if(ibin >= 61 && ibin < 91) {
      xbins[ibin] = xbins[60] + 60.0*(ibin-60);
    }
    if(ibin >= 91 && ibin < 121) {
      xbins[ibin] = xbins[90] + 100.0*(ibin-90);
    }
    if(ibin >= 121) {
      xbins[ibin] = xbins[120] + 200.0*(ibin-120);
    }
  }

  // histograms to store denominators and numerators
  // vs MET
  TH1F* hMet = new TH1F("hMet","",nbins,xbins);                                 // denominator
  TH1F* hMet_PassTrig = new TH1F("hMet_PassTrig","",nbins,xbins);               // numerator for full trigger
  TH1F* hST = new TH1F("hST","",nbins,xbins);                                 // denominator
  TH1F* hST_PassTrig = new TH1F("hST_PassTrig","",nbins,xbins);               // numerator for full trigger
  TH1F* hHT = new TH1F("hHT","",nbins,xbins);                                 // denominator
  TH1F* hHT_PassTrig = new TH1F("hHT_PassTrig","",nbins,xbins);               // numerator for full trigger

TH1F* hMet_0 = new TH1F("hMet_0","",nbins,xbins);                                 // denominator
  TH1F* hMet_PassTrig_0 = new TH1F("hMet_PassTrig_0","",nbins,xbins);               // numerator for full trigger
  TH1F* hST_0 = new TH1F("hST_0","",nbins,xbins);                                 // denominator
  TH1F* hST_PassTrig_0 = new TH1F("hST_PassTrig_0","",nbins,xbins);               // numerator for full trigger
  TH1F* hHT_0 = new TH1F("hHT_0","",nbins,xbins);                                 // denominator
  TH1F* hHT_PassTrig_0 = new TH1F("hHT_PassTrig_0","",nbins,xbins);    

  TH1F* hMet_1 = new TH1F("hMet_1","",nbins,xbins);                                 // denominator
  TH1F* hMet_PassTrig_1 = new TH1F("hMet_PassTrig_1","",nbins,xbins);               // numerator for full trigger
  TH1F* hST_1 = new TH1F("hST_1","",nbins,xbins);                                 // denominator
  TH1F* hST_PassTrig_1 = new TH1F("hST_PassTrig_1","",nbins,xbins);               // numerator for full trigger
  TH1F* hHT_1 = new TH1F("hHT_1","",nbins,xbins);                                 // denominator
  TH1F* hHT_PassTrig_1 = new TH1F("hHT_PassTrig_1","",nbins,xbins);   


  TH1F* hMet_2 = new TH1F("hMet_2","",nbins,xbins);                                 // denominator
  TH1F* hMet_PassTrig_2 = new TH1F("hMet_PassTrig_2","",nbins,xbins);               // numerator for full trigger
  TH1F* hST_2 = new TH1F("hST_2","",nbins,xbins);                                 // denominator
  TH1F* hST_PassTrig_2 = new TH1F("hST_PassTrig_2","",nbins,xbins);               // numerator for full trigger
  TH1F* hHT_2 = new TH1F("hHT_2","",nbins,xbins);                                 // denominator
  TH1F* hHT_PassTrig_2 = new TH1F("hHT_PassTrig_2","",nbins,xbins);    

  TH1F* hMet_3 = new TH1F("hMet_3","",nbins,xbins);                                 // denominator
  TH1F* hMet_PassTrig_3 = new TH1F("hMet_PassTrig_3","",nbins,xbins);               // numerator for full trigger
  TH1F* hST_3 = new TH1F("hST_3","",nbins,xbins);                                 // denominator
  TH1F* hST_PassTrig_3 = new TH1F("hST_PassTrig_3","",nbins,xbins);               // numerator for full trigger
  TH1F* hHT_3 = new TH1F("hHT_3","",nbins,xbins);                                 // denominator
  TH1F* hHT_PassTrig_3 = new TH1F("hHT_PassTrig_3","",nbins,xbins);    
 
  TH1F* hMet_4 = new TH1F("hMet_4","",nbins,xbins);                                 // denominator
  TH1F* hMet_PassTrig_4 = new TH1F("hMet_PassTrig_4","",nbins,xbins);               // numerator for full trigger
  TH1F* hST_4 = new TH1F("hST_4","",nbins,xbins);                                 // denominator
  TH1F* hST_PassTrig_4 = new TH1F("hST_PassTrig_4","",nbins,xbins);               // numerator for full trigger
  TH1F* hHT_4 = new TH1F("hHT_4","",nbins,xbins);                                 // denominator
  TH1F* hHT_PassTrig_4 = new TH1F("hHT_PassTrig_4","",nbins,xbins);    

  TH1F* hMet_5 = new TH1F("hMet_5","",nbins,xbins);                                 // denominator
  TH1F* hMet_PassTrig_5 = new TH1F("hMet_PassTrig_5","",nbins,xbins);               // numerator for full trigger
  TH1F* hST_5 = new TH1F("hST_5","",nbins,xbins);                                 // denominator
  TH1F* hST_PassTrig_5 = new TH1F("hST_PassTrig_5","",nbins,xbins);               // numerator for full trigger
  TH1F* hHT_5 = new TH1F("hHT_5","",nbins,xbins);                                 // denominator
  TH1F* hHT_PassTrig_5 = new TH1F("hHT_PassTrig_5","",nbins,xbins);    
  //TH1F* hMupT = new TH1F("hMupT","",100,5,500);                                 // denominator
  TH1F* hMupT = new TH1F("hMupT","",nbins,xbins);                                 // denominator
  TH1F* hMupT_PassTrig = new TH1F("hMupT_PassTrig","",nbins,xbins);                                 // denominator
  TH1F* hGenST = new TH1F("hGenST","",nbins,xbins);                                 // denominator
  TH1F* hGenST_PassTrig = new TH1F("hGenST_PassTrig","",nbins,xbins);               // numerator for full trigger
  TH1F* hGenHT = new TH1F("hGenHT","",nbins,xbins);                                 // denominator
  TH1F* hGenHT_PassTrig = new TH1F("hGenHT_PassTrig","",nbins,xbins);               // numerator for full trigger
  // vs sub-leading jet pT
  TH1F* hGenJet2Pt = new TH1F("hGenJet2Pt","",nbins,xbins);                                 // denominator
  TH1F* hGenJet2Pt_PassDijetFilter = new TH1F("hGenJet2Pt_PassthreeFilter","",nbins,xbins); // numerator for jet part of trigger

  for(int ievent = 0; ievent < nevents; ievent++) {
    chain->GetEntry(ievent);

    int ncenjets30 = 0, ncenjets70 = 0;
    float jet1pt = 0.0, jet2pt = 0.0, HT =0.0; 
    float Mu1pt = 0.0 , Mu1eta =0.0;
    for(int ijet = 0; ijet < ngenjets; ijet++) {
      if(fabs(genjeteta[ijet]) < 2.4 && genjetpt[ijet] > 30.0) {
        bool overlap = false;
       HT = HT+genjetpt[ijet];
        // check whether there is overlap with the muon trigger object which fires the auxiliary trigger ... if so, don't count the jet
        for(int imu = 0; imu < nhltmuons; imu++) {
          if(deltaR(hltmuoneta[imu], hltmuonphi[imu], genjeteta[ijet], genjetphi[ijet]) < 0.4) {
            overlap = true;
          }
        }
       // if(!overlap) {
          ncenjets30++;
       // }
        if(genjetpt[ijet] > jet1pt) {
          jet1pt = genjetpt[ijet];
        }
        if(genjetpt[ijet] > jet2pt && genjetpt[ijet] < jet1pt) {
          jet2pt = genjetpt[ijet];
        }
        if(genjetpt[ijet] > 30.0) {
        //  if(!overlap) 
            ncenjets70++;
        }
      }
    }
      int Nmu =0;
      int Nmu_V =0;
      for(int imu = 0; imu < nmuons; imu++) {
       if(fabs(muoneta[imu])< 2.1){
       if(muonpt[imu] > Mu1pt) {
          Mu1pt = muonpt[imu];
          }
       if(muonpt[imu] > 25) 
          Nmu ++;
         else if (muonpt[imu] > 10) Nmu_V ++;
     }

       } 
     //std::cout<<HT<<std::endl;
     // check whether the event passed the DiJet filter
    int npassjets = 0;
    for(int ijet = 0; ijet < nhltjets; ijet++) {
      if(hltjetpt[ijet] > 50.0 && fabs(hltjeteta[ijet]) < 2.6) { npassjets++; }
    }
    int npassmu = 0;
    for(int imu = 0; imu < nhltmuons; imu++) {
      if(hltmuonpt[imu] > 50.0 && fabs(hltmuoneta[imu]) < 2.1) { npassmu++; }
    }
    bool passdijetfilter = npassjets > 1;
    bool passmufilter = npassmu > 1;
    float ST = genmet+Mu1pt ;
    // if (Mu1pt > 45){
    //if( Nmu == 1 && ncenjets30 >= 3 && Nmu_V ==0  ) {  //offline selection
    if(  ncenjets30 >= 3   ) {  //offline selection
    if(passaux) { // offline selection
      hMet->Fill(genmet);
      hST->Fill(ST);     
      hHT->Fill(HT);
      if(  passtrig) {
       hMet_PassTrig->Fill(genmet);
      hST_PassTrig->Fill(ST);
     hHT_PassTrig->Fill(HT);
    }
 if( passtrig0) {
       hMet_PassTrig_0->Fill(genmet);
      hST_PassTrig_0->Fill(ST);
       hHT_PassTrig_0->Fill(HT);
    }

   if(  passtrig1) {
       hMet_PassTrig_1->Fill(genmet);
      hST_PassTrig_1->Fill(ST);
       hHT_PassTrig_1->Fill(HT);
    }

if(  passtrig1 || passtrig3) {
       hMet_PassTrig_2->Fill(genmet);
      hST_PassTrig_2->Fill(ST);
       hHT_PassTrig_2->Fill(HT);
    }

   if(  passtrig3) {
       hMet_PassTrig_3->Fill(genmet);
      hST_PassTrig_3->Fill(ST);
       hHT_PassTrig_3->Fill(HT);
    }
 
 if(  passtrig4) {
       hMet_PassTrig_4->Fill(genmet);
      hST_PassTrig_4->Fill(ST);
       hHT_PassTrig_4->Fill(HT);
    }

   if(  passtrig1|| passtrig3 || passtrig4 ) {
       hMet_PassTrig_5->Fill(genmet);
      hST_PassTrig_5->Fill(ST);
       hHT_PassTrig_5->Fill(HT);
    }
    }
   }
  }

  // make efficiency graphs
  TGraphAsymmErrors *effvsMet = makeEffGraph(hMet_PassTrig, hMet);
  TGraphAsymmErrors *effvsST = makeEffGraph(hST_PassTrig, hST);
  TGraphAsymmErrors *effvsHT = makeEffGraph(hHT_PassTrig, hHT);
  TGraphAsymmErrors *effvsMupT = makeEffGraph(hMupT_PassTrig, hMupT);
  TGraphAsymmErrors *effDijetvsGenJet2Pt = makeEffGraph(hGenJet2Pt_PassDijetFilter, hGenJet2Pt);
  TGraphAsymmErrors *effvsMet_0 =makeEffGraph(hMet_PassTrig_0, hMet);
  TGraphAsymmErrors *effvsST_0 =makeEffGraph(hST_PassTrig_0, hST);
  TGraphAsymmErrors *effvsHT_0 =makeEffGraph(hHT_PassTrig_0, hHT);
  TGraphAsymmErrors *effvsMet_1 =makeEffGraph(hMet_PassTrig_1, hMet);
  TGraphAsymmErrors *effvsST_1 =makeEffGraph(hST_PassTrig_1, hST);
  TGraphAsymmErrors *effvsHT_1 =makeEffGraph(hHT_PassTrig_1, hHT);
   TGraphAsymmErrors *effvsMet_2 =makeEffGraph(hMet_PassTrig_2, hMet);
  TGraphAsymmErrors *effvsST_2 =makeEffGraph(hST_PassTrig_2, hST);
  TGraphAsymmErrors *effvsHT_2 =makeEffGraph(hHT_PassTrig_2, hHT);
 TGraphAsymmErrors *effvsMet_3 =makeEffGraph(hMet_PassTrig_3, hMet);
  TGraphAsymmErrors *effvsST_3 =makeEffGraph(hST_PassTrig_3, hST);
  TGraphAsymmErrors *effvsHT_3 =makeEffGraph(hHT_PassTrig_3, hHT);
TGraphAsymmErrors *effvsMet_4 =makeEffGraph(hMet_PassTrig_4, hMet);
  TGraphAsymmErrors *effvsST_4 =makeEffGraph(hST_PassTrig_4, hST);
  TGraphAsymmErrors *effvsHT_4 =makeEffGraph(hHT_PassTrig_4, hHT);
  TGraphAsymmErrors *effvsMet_5 =makeEffGraph(hMet_PassTrig_5, hMet);
  TGraphAsymmErrors *effvsST_5 =makeEffGraph(hST_PassTrig_5, hST);
  TGraphAsymmErrors *effvsHT_5 =makeEffGraph(hHT_PassTrig_5, hHT);


  // make plots
  SetStyle();
  TCanvas *c = MakeCanvas("effcanvas","HLT_PFHT350_PFMET120_NoiseCleaned_v1",600,600);
  InitGraph(effvsMet,"HLT_PFHT350_PFMET120_NoiseCleaned_v1","#slash{E}_{T} [GeV]","Efficiency",0.0,800.0,0.0,1.1,kRed);
  effvsMet->SetMarkerStyle(kFullCircle);
  effvsMet->SetMarkerSize(0.9); 
  effvsMet->SetMarkerStyle(kFullCircle);
  effvsMet->SetMarkerSize(0.9);
  effvsMet->Draw("AP");
 c->SaveAs("./fig/eff_trig_vs_MET.pdf");

 InitGraph(effvsST,"HLT_PFHT350_PFMET120_NoiseCleaned_v1","S_{T} [GeV]","Efficiency",0.0,800.0,0.0,1.1,kRed);
  effvsST->SetMarkerStyle(kFullCircle);
  effvsST->SetMarkerSize(0.9);
  effvsST->Draw("AP");
 c->SaveAs("./fig/eff_trig_vs_ST.pdf");

  InitGraph(effvsHT,"HLT_PFHT350_PFMET120_NoiseCleaned_v1","HT [GeV]","Efficiency",300,3000.0,0.0,1.1,kRed);
  effvsHT->SetMarkerStyle(kFullCircle);
  effvsHT->SetMarkerSize(0.9);
  effvsHT->Draw("AP");
 c->SaveAs("./fig/eff_trig_vs_HT.pdf");

 InitGraph(effvsMet_0,"HLT_PFHT900_v1","#slash{E}_{T} [GeV]","Efficiency",0.0,800.0,0.0,1.1,kRed);
  effvsMet_0->SetMarkerStyle(kFullCircle);
  effvsMet_0->SetMarkerSize(0.9); 
  effvsMet_0->SetMarkerStyle(kFullCircle);
  effvsMet_0->SetMarkerSize(0.9);
  effvsMet_0->Draw("AP");
 c->SaveAs("./fig/eff_trig_vs_MET_0.pdf");

 InitGraph(effvsST_0,"HLT_PFHT900_v1","S_{T} [GeV]","Efficiency",0.0,800.0,0.0,1.1,kRed);
  effvsST_0->SetMarkerStyle(kFullCircle);
  effvsST_0->SetMarkerSize(0.9);
  effvsST_0->Draw("AP");
 c->SaveAs("./fig/eff_trig_vs_ST_0.pdf");

  InitGraph(effvsHT_0,"HLT_PFHT900_v1","HT [GeV]","Efficiency",900,3000.0,0.0,1.1,kRed);
  effvsHT_0->SetMarkerStyle(kFullCircle);
  effvsHT_0->SetMarkerSize(0.9);
  effvsHT_0->Draw("AP");
 c->SaveAs("./fig/eff_trig_vs_HT_0.pdf");


InitGraph(effvsMet_1,"HLT_Mu15_IsoVVVL_PFHT600_v1","#slash{E}_{T} [GeV]","Efficiency",0.0,800.0,0.0,1.1,kRed);
  effvsMet_1->SetMarkerStyle(kFullCircle);
  effvsMet_1->SetMarkerSize(0.9); 
  effvsMet_1->SetMarkerStyle(kFullCircle);
  effvsMet_1->SetMarkerSize(0.9);
  effvsMet_1->Draw("AP");
 c->SaveAs("./fig/eff_trig_vs_MET_1.pdf");

 InitGraph(effvsST_1,"HLT_Mu15_IsoVVVL_PFHT600_v1","S_{T} [GeV]","Efficiency",0.0,800.0,0.0,1.1,kRed);
  effvsST_1->SetMarkerStyle(kFullCircle);
  effvsST_1->SetMarkerSize(0.9);
  effvsST_1->Draw("AP");
 c->SaveAs("./fig/eff_trig_vs_ST_1.pdf");

  InitGraph(effvsHT_1,"HLT_Mu15_IsoVVVL_PFHT600_v1","HT [GeV]","Efficiency",600,3000.0,0.0,1.1,kRed);
  effvsHT_1->SetMarkerStyle(kFullCircle);
  effvsHT_1->SetMarkerSize(0.9);
  effvsHT_1->Draw("AP");
 c->SaveAs("./fig/eff_trig_vs_HT_1.pdf");

InitGraph(effvsMet_2,"HLT_Mu15_IsoVVVL_PFHT600 || HLT_Mu15_IsoVVVL_PFHT400_PFMET70","#slash{E}_{T} [GeV]","Efficiency",0.0,800.0,0.0,1.1,kRed);
  effvsMet_2->SetMarkerStyle(kFullCircle);
  effvsMet_2->SetMarkerSize(0.9); 
  effvsMet_2->SetMarkerStyle(kFullCircle);
  effvsMet_2->SetMarkerSize(0.9);
  effvsMet_2->Draw("AP");
 c->SaveAs("./fig/eff_trig_vs_MET_2.pdf");

 InitGraph(effvsST_2,"HLT_Mu15_IsoVVVL_PFHT600 || HLT_Mu15_IsoVVVL_PFHT400_PFMET70","S_{T} [GeV]","Efficiency",0.0,800.0,0.0,1.1,kRed);
  effvsST_2->SetMarkerStyle(kFullCircle);
  effvsST_2->SetMarkerSize(0.9);
  effvsST_2->Draw("AP");
 c->SaveAs("./fig/eff_trig_vs_ST_2.pdf");

  InitGraph(effvsHT_2,"HLT_Mu15_IsoVVVL_PFHT600 || HLT_Mu15_IsoVVVL_PFHT400_PFMET70","HT [GeV]","Efficiency",400,3000.0,0.0,1.1,kRed);
  effvsHT_2->SetMarkerStyle(kFullCircle);
  effvsHT_2->SetMarkerSize(0.9);
  effvsHT_2->Draw("AP");
 c->SaveAs("./fig/eff_trig_vs_HT_2.pdf");


InitGraph(effvsMet_3,"HLT_Mu15_IsoVVVL_PFHT400_PFMET70_v1","#slash{E}_{T} [GeV]","Efficiency",0.0,800.0,0.0,1.1,kRed);
  effvsMet_3->SetMarkerStyle(kFullCircle);
  effvsMet_3->SetMarkerSize(0.9); 
  effvsMet_3->SetMarkerStyle(kFullCircle);
  effvsMet_3->SetMarkerSize(0.9);
  effvsMet_3->Draw("AP");
 c->SaveAs("./fig/eff_trig_vs_MET_3.pdf");

 InitGraph(effvsST_3,"HLT_Mu15_IsoVVVL_PFHT400_PFMET70_v1","S_{T} [GeV]","Efficiency",0.0,800.0,0.0,1.1,kRed);
  effvsST_3->SetMarkerStyle(kFullCircle);
  effvsST_3->SetMarkerSize(0.9);
  effvsST_3->Draw("AP");
 c->SaveAs("./fig/eff_trig_vs_ST_3.pdf");

  InitGraph(effvsHT_3,"HLT_Mu15_IsoVVVL_PFHT400_PFMET70_v1","HT [GeV]","Efficiency",600,3000.0,0.0,1.1,kRed);
  effvsHT_3->SetMarkerStyle(kFullCircle);
  effvsHT_3->SetMarkerSize(0.9);
  effvsHT_3->Draw("AP");
 c->SaveAs("./fig/eff_trig_vs_HT_3.pdf");

InitGraph(effvsMet_4,"HLT_PFMET120_NoiseCleaned_Mu5_v1","#slash{E}_{T} [GeV]","Efficiency",0.0,800.0,0.0,1.1,kRed);
  effvsMet_4->SetMarkerStyle(kFullCircle);
  effvsMet_4->SetMarkerSize(0.9); 
  effvsMet_4->SetMarkerStyle(kFullCircle);
  effvsMet_4->SetMarkerSize(0.9);
  effvsMet_4->Draw("AP");
 c->SaveAs("./fig/eff_trig_vs_MET_4.pdf");

 InitGraph(effvsST_4,"HLT_PFMET120_NoiseCleaned_Mu5_v1","S_{T} [GeV]","Efficiency",0.0,800.0,0.0,1.1,kRed);
  effvsST_4->SetMarkerStyle(kFullCircle);
  effvsST_4->SetMarkerSize(0.9);
  effvsST_4->Draw("AP");
 c->SaveAs("./fig/eff_trig_vs_ST_4.pdf");

  InitGraph(effvsHT_4,"HLT_PFMET120_NoiseCleaned_Mu5_v1","HT [GeV]","Efficiency",600,3000.0,0.0,1.1,kRed);
  effvsHT_4->SetMarkerStyle(kFullCircle);
  effvsHT_4->SetMarkerSize(0.9);
  effvsHT_4->Draw("AP");
 c->SaveAs("./fig/eff_trig_vs_HT_4.pdf");

InitGraph(effvsMet_5,"HLT_Mu15_IsoVVVL_PFHT600 || HLT_Mu15_IsoVVVL_PFHT400_PFMET70 || HLT_PFMET120_NoiseCleaned_Mu5","#slash{E}_{T} [GeV]","Efficiency",0.0,800.0,0.0,1.1,kRed);
  effvsMet_5->SetMarkerStyle(kFullCircle);
  effvsMet_5->SetMarkerSize(0.9); 
  effvsMet_5->SetMarkerStyle(kFullCircle);
  effvsMet_5->SetMarkerSize(0.9);
  effvsMet_5->Draw("AP");
 c->SaveAs("./fig/eff_trig_vs_MET_5.pdf");

 InitGraph(effvsST_5,"HLT_Mu15_IsoVVVL_PFHT600 || HLT_Mu15_IsoVVVL_PFHT400_PFMET70|| HLT_PFMET120_NoiseCleaned_Mu5 ","S_{T} [GeV]","Efficiency",0.0,800.0,0.0,1.1,kRed);
  effvsST_5->SetMarkerStyle(kFullCircle);
  effvsST_5->SetMarkerSize(0.9);
  effvsST_5->Draw("AP");
 c->SaveAs("./fig/eff_trig_vs_ST_5.pdf");

  InitGraph(effvsHT_5,"HLT_Mu15_IsoVVVL_PFHT600 || HLT_Mu15_IsoVVVL_PFHT400_PFMET70|| HLT_PFMET120_NoiseCleaned_Mu5","HT [GeV]","Efficiency",400,3000.0,0.0,1.1,kRed);
  effvsHT_5->SetMarkerStyle(kFullCircle);
  effvsHT_5->SetMarkerSize(0.9);
  effvsHT_5->Draw("AP");
 c->SaveAs("./fig/eff_trig_vs_HT_5.pdf");



  TCanvas *c1 = MakeCanvas("MupT","",600,600);
 // hMupT->Draw();  
  hHT->Draw("");
  c1->SaveAs("HT.pdf");


}

void makeTrigEffPlots_mu()
{

  std::vector<TString> filelist;

  filelist.push_back("merge.root");

  runTrigEffPlots(filelist);

}
