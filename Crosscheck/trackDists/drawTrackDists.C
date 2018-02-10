#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

#include <iostream>

void drawTrackDists(){

  TFile * f = TFile::Open("/data/cmcginn/StudyMultSamples/ALEPH/LEP1/20180119/LEP1Data1995_recons_aftercut-MERGED.root","read");
  TTree * t = (TTree*)f->Get("t");
  int nParticle;
  float d0[10000];
  float z0[10000];
  int ntpc[10000];
  t->SetBranchAddress("nParticle",&nParticle);
  t->SetBranchAddress("d0",d0);
  t->SetBranchAddress("z0",z0);
  t->SetBranchAddress("ntpc",ntpc);
 
  TH1D * ntpc_h = new TH1D("ntpc_h",";n TPC Clusters;N_{trk}",25,0,25);
  TH1D * d0_h = new TH1D("d0_h",";d0;N_{trk}",50,-2,2);
  TH1D * z0_h = new TH1D("z0_h",";z0;N_{trk}",50,-5,5);

  for(int i = 0; i<t->GetEntries(); i++){
    t->GetEntry(i);
    for(int t = 0; t<nParticle; t++){
      if(ntpc[t]<0) continue;
      ntpc_h->Fill(ntpc[t]);
      d0_h->Fill(d0[t]);
      z0_h->Fill(z0[t]);
    }
  }

  TCanvas * c = new TCanvas("c","c",800,600);
  ntpc_h->Draw();
  c->SaveAs("img/ntpc.pdf");
  c->SaveAs("img/ntpc.png");
  c->SaveAs("img/ntpc.C");
  d0_h->Draw();
  c->SaveAs("img/d0.pdf");
  c->SaveAs("img/d0.png");
  c->SaveAs("img/d0.C");
  z0_h->Draw();
  c->SaveAs("img/z0.pdf");
  c->SaveAs("img/z0.png");
  c->SaveAs("img/z0.C");
}
