#include "TFile.h"
#include "TTree.h"
#include "TAttLine.h"
#include "THelix.h"
#include "TView.h"
#include <string>

void formatHelix(THelix * h){
  h->SetLineColor(kRed);
  h->SetLineWidth(2);
}

void EventDisplay(std::string inputFile = "/data/abaty/ALEPHTrees/mergedLEP1_20171022.root" ,int eventIndx = 0){

  TFile * f = TFile::Open(inputFile.c_str(),"read");
  TTree * t = (TTree*)f->Get("t");

  int nParticle;
  float px[500];
  float py[500];
  float pz[500];
  float charge[500];
  int pwflag[500], pwflagMix[500];

  t->SetBranchAddress("nParticle",&nParticle);
  t->SetBranchAddress("px",&px);              
  t->SetBranchAddress("py",&py);              
  t->SetBranchAddress("pz",&pz);              
  t->SetBranchAddress("pwflag",&pwflag);     
  t->SetBranchAddress("charge",&charge); 
  t->GetEntry(eventIndx);

  TCanvas* helix_example_c1 = new TCanvas("helix_example_c1");
  TView *view = TView::CreateView(1);
  view->SetRange(-1,-1,-1,1,1,1);

  THelix * helix[1000];
  int nHelix = 0;
  for(int i = 0; i<nParticle; i++){
    if(pwflag[i]!=0) continue;
    std::cout << i<<" "  << px[i]<<" "  << py[i]<<" " << pz[i]<<" " << charge[i] << std::endl; 
    helix[nHelix] = new THelix(0,0,0,px[i],py[i],pz[i],charge[i]*1.5);
    formatHelix(helix[nHelix]);
    helix[nHelix]->Draw("same");
    nHelix++;
  }

}
