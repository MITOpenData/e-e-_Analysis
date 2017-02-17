#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>

using namespace std;


#define PI 3.1415926
enum SIMPLEPID {PHOTON, ELECTRON, PION, MUON, KAON, PROTON};

void analysis(int maxevt=0,int mult=0,int nbin=20){

  TFile *f = new TFile("output-2.root");
  TTree *t1 = (TTree*)f->Get("t");
  Int_t nParticle;
  Float_t pt[10000];
  Float_t eta[10000];
  Float_t pid[10000];
  Float_t phi[10000];
  Float_t mass[10000];
  
  
  TH1F*hpt[3]=new 
  
  hpt[0]
  
  
  t1->SetBranchAddress("nParticle",&nParticle);
  t1->SetBranchAddress("pt",pt);
  t1->SetBranchAddress("eta",eta);
  t1->SetBranchAddress("pid",pid);
  t1->SetBranchAddress("phi",phi);
  t1->SetBranchAddress("mass",mass);

  // two histograms
  TH1F *hpt = new TH1F ( "hpt", "hpt ",100,0.,20.);

  // all entries and fill the histograms
  Int_t nevent = (Int_t)t1->GetEntries();
  
  int nevent_process = nevent;
  if( maxevt>0 && maxevt<nevent ) nevent_process = maxevt;  
  
  for (Int_t i=0;i<nevent_process;i++) {
    if (i%1000==0) cout <<i<<"/"<<nevent_process<<endl;
    t1->GetEntry(i);
    
    int nparticles = nParticle;
    if (nparticles<mult) continue;

    for ( int j=0;j<nparticles;j++ ) {
    
      int pid1 = pid[j];
      if (pid1!=PION&&pid1!=PROTON&&pid1!=KAON) continue;
      float eta1 = eta[j];
      float phi1 = phi[j];
      float pt1 = pt[j];
      float mass1 = mass[j];
      if(pt1<0.||pt1>4.) continue;
      if(pid1==PION)hpt[0]->Fill(pt1);      
      }//end of second loop 
  }// end of loop over events  
      hpt->Draw();

}