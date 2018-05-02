#include "Settings.h"
#include "Tools.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TAxis.h"
#include "TVector3.h"
#include <iostream>
#include <cmath>
#include <string>
#include "../../DataProcessing/include/eventData.h"
#include "../../DataProcessing/include/particleData.h"
#include "../../DataProcessing/include/alephTrkEfficiency.h"

//#include "include/smartJetName.h"

//inline Float_t jtp(Float_t pt, Float_t eta){ return pt*TMath::CosH(eta); }

void Analyzer(){
  TH1::SetDefaultSumw2();

  alephTrkEfficiency a;

  //set up plots
  Settings s = Settings();

  TFile * output = TFile::Open("Output.root","recreate");
  TH2F * signal2PC[s.nMultBins]; 
  TH2F * bkgrnd2PC[s.nMultBins];
  TH2F * ratio2PC[s.nMultBins]; 
  TH1F * longRangeYield[s.nMultBins]; 
  Float_t nSignalEvts[s.nMultBins] = {0};
  Float_t nBkgrndEvts[s.nMultBins] = {0};

  TH1D * nEvtSigHist = new TH1D("nEvtSigHisto","nEvtSigHisto",10,0,10);
  TH1D * nEvtBkgHist = new TH1D("nEvtBkgHisto","nEvtBkgHisto",10,0,10);
  
  for(int i = 0; i<s.nMultBins; i++){
    signal2PC[i] = new TH2F(Form("signal2PC_%d_%d",s.multBinsLow[i],s.multBinsHigh[i]),";#Delta#eta;#Delta#Phi",s.dEtaBins,-s.etaPlotRange,s.etaPlotRange,s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
    bkgrnd2PC[i] = new TH2F(Form("bkgrnd2PC_%d_%d",s.multBinsLow[i],s.multBinsHigh[i]),";#Delta#eta;#Delta#Phi",s.dEtaBins,-s.etaPlotRange,s.etaPlotRange,s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
    longRangeYield[i] = new TH1F(Form("longRangeYield_%d_%d",s.multBinsLow[i],s.multBinsHigh[i]),";#Delta#phi;Y(#Delta#Phi)",s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0); 
  }
  TH1F * multiplicity = new TH1F("multiplicity",";nTrk;nEvents",200,0,200);

  //files and variables for input
  TFile * f = TFile::Open(s.inputFile.c_str(),"read");
  TTree * t = (TTree*)f->Get("t"); 
  TFile * fMix = TFile::Open(s.inputFileMix.c_str(),"read");
  TTree * tMix = (TTree*)fMix->Get("t");

  eventData ed, edm;
  particleData pd, pdm;
 
  std::vector<std::string> inList;
  ed.SetStatusAndAddressRead(t,inList);
  pd.SetStatusAndAddressRead(t,inList);
  edm.SetStatusAndAddressRead(tMix,inList); 
  pdm.SetStatusAndAddressRead(tMix,inList); 

  t->GetEntry(0);
  tMix->GetEntry(0);


  for(int i = 0; i< (s.doAllData?t->GetEntries():s.nEvts); i++){
    if(i%1000==0) std::cout << i << "/" << (s.doAllData?t->GetEntries():s.nEvts)<< std::endl;  
    t->GetEntry(i);

    //if(!ed.passesAll) continue;
    if(!ed.passesLEP1TwoPC) continue;
    tMix->GetEntry(i);
 
    int nTrk = ed.nChargedHadronsHP;
    Float_t nTrig = 0;
    for(int t = 0; t<pd.nParticle; t++){
      if(!pd.highPurity[t]) continue;
      if(pd.pwflag[t]>2) continue;
      if(TMath::Abs(pd.eta[t])>s.etaCut) continue;
      if(!s.doThrust && (pd.pt[t]<s.trigPt[0] || pd.pt[t]>=s.trigPt[1])) continue;//nTrig calculation
      if(s.doThrust && (pd.pt_wrtThr[t]<s.trigPt[0] || pd.pt_wrtThr[t]>=s.trigPt[1])) continue;//nTrig calculation
      Float_t corr = 1.0/a.efficiency(pd.theta[t],pd.phi[t],pd.pt[t]);
      nTrig += corr;
    }
    multiplicity->Fill(nTrk);
 
    for(int k = 0; k<s.nMultBins; k++){
      if(s.isInMultBin(nTrk,k))  nSignalEvts[k]++;
      if(s.isInMultBin(nTrk,k))  nBkgrndEvts[k] ++;
    }
 
      //fill signal histogram
    for(int j1 = 0; j1<pd.nParticle; j1++){
      if(!pd.highPurity[j1]) continue;
      if(pd.pwflag[j1]>2) continue;
      if(TMath::Abs(pd.eta[j1])>s.etaCut) continue;
      if(!s.doThrust && (pd.pt[j1]<s.trigPt[0] || pd.pt[j1]>s.trigPt[1])) continue;
      if(s.doThrust && (pd.pt_wrtThr[j1]<s.trigPt[0] || pd.pt_wrtThr[j1]>s.trigPt[1])) continue;
      Float_t corr1 = 1.0/a.efficiency(pd.theta[j1],pd.phi[j1],pd.pt[j1]);

      //signal histogram
      for(int j2 = 0; j2<j1; j2++){
        if(!pd.highPurity[j2]) continue;
        if(pd.pwflag[j2]>2) continue;
        if(TMath::Abs(pd.eta[j2])>s.etaCut) continue;
        if(!s.doThrust && (pd.pt[j2]<s.assocPt[0] || pd.pt[j2]>s.assocPt[1])) continue;
        if(s.doThrust && (pd.pt_wrtThr[j2]<s.assocPt[0] || pd.pt_wrtThr[j2]>s.assocPt[1])) continue;
        Float_t corr2 = 1.0/a.efficiency(pd.theta[j2],pd.phi[j2],pd.pt[j2]);
       
        //correct for both particles and also divide by hte bin widths
        for(int k = 0; k<s.nMultBins; k++){
          if(s.isInMultBin(nTrk,k)){
            if(!s.doThrust) fillAllQuadrants( signal2PC[k], dEta(pd.eta[j1],pd.eta[j2]), dPhi(pd.phi[j1],pd.phi[j2]), corr1*corr2/(2*s.etaPlotRange/(Float_t)s.dEtaBins)/(2*TMath::Pi()/(Float_t)s.dPhiBins)/nTrig);
            if(s.doThrust) fillAllQuadrants( signal2PC[k], dEta(pd.eta_wrtThr[j1],pd.eta_wrtThr[j2]), dPhi(pd.phi_wrtThr[j1],pd.phi_wrtThr[j2]), corr1*corr2/(2*s.etaPlotRange/(Float_t)s.dEtaBins)/(2*TMath::Pi()/(Float_t)s.dPhiBins)/nTrig);
          }
        }//end mult bin loop
      }//end 2nd particle loop
      
      //background mixed event histogram 
      for(int j2 = 0; j2<pdm.nParticle; j2++){
        if(!pdm.highPurity[j2]) continue;
        if(pdm.pwflag[j2]>2) continue;
        if(TMath::Abs(pdm.eta[j2])>s.etaCut) continue;
        if(!s.doThrust && (pdm.pt[j2]<s.assocPt[0] || pdm.pt[j2]>s.assocPt[1])) continue;
        if(s.doThrust && (pdm.pt_wrtThr[j2]<s.assocPt[0] || pdm.pt_wrtThr[j2]>s.assocPt[1])) continue;
        Float_t corr2 = pdm.particleWeight/a.efficiency(pdm.theta[j2],pdm.phi[j2],pdm.pt[j2]);
        //correct for both particles and also divide by hte bin widths
        for(int k = 0; k<s.nMultBins; k++){
          if(s.isInMultBin(nTrk,k)){
            if(!s.doThrust) fillAllQuadrants( bkgrnd2PC[k], dEta(pd.eta[j1],pdm.eta[j2]), dPhi(pd.phi[j1],pdm.phi[j2]), corr1*corr2/(2*s.etaPlotRange/(Float_t)s.dEtaBins)/(2*TMath::Pi()/(Float_t)s.dPhiBins)/nTrig);
            if(s.doThrust) fillAllQuadrants( bkgrnd2PC[k], dEta(pd.eta_wrtThr[j1],pdm.eta_wrtThr[j2]), dPhi(pd.phi_wrtThr[j1],pdm.phi_wrtThr[j2]), corr1*corr2/(2*s.etaPlotRange/(Float_t)s.dEtaBins)/(2*TMath::Pi()/(Float_t)s.dPhiBins)/nTrig);
          }
        }//end mult bin loop
      }//end mixed particle loop
    }//end trigger particle loop
  }//end main evt loop
  
  output->cd();
  for(int k = 0; k<s.nMultBins; k++){
    std::cout << nSignalEvts[k] <<" " << std::endl;
    if(!s.doParallel) signal2PC[k]->Scale(1.0/(Float_t)nSignalEvts[k]);
    if(!s.doParallel) bkgrnd2PC[k]->Scale(1.0/(Float_t)nBkgrndEvts[k]);

    nEvtSigHist->Fill(k,nSignalEvts[k]);
    nEvtBkgHist->Fill(k,nBkgrndEvts[k]);

    ratio2PC[k] = (TH2F*)signal2PC[k]->Clone(Form("ratio2PC_%d_%d",s.multBinsLow[k],s.multBinsHigh[k]));
    ratio2PC[k]->Divide(bkgrnd2PC[k]);
    ratio2PC[k]->Scale(bkgrnd2PC[k]->GetBinContent(bkgrnd2PC[k]->FindBin(0,0)));
    getLongRangeYield(s,ratio2PC[k],longRangeYield[k]);
  }

  f->Close();
  fMix->Close();
  output->Write();
  exit(1); 
}
