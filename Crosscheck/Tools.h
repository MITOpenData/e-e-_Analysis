#ifndef TOOLS
#define TOOLS

#include "TMath.h"
#include "TH2F.h"

inline float dPhi(float phi1, float phi2){
  return TMath::ACos(TMath::Cos(phi1-phi2));
}

inline float dEta(float eta1, float eta2){
  return TMath::Abs(eta1-eta2);
}

void symmetrizeDetaDphi(TH2F * h, int etaBins, int phiBins){
  for(int i = etaBins/2; i<etaBins; i++){
    for(int j = phiBins/4 ; j<3*phiBins/4; j++){
      int bin = phiBins/2-j;
      if(bin<1) bin = phiBins-(j-phiBins/2);
      
      h->SetBinContent(i+1,bin,h->GetBinContent(i+1,j+1));
      h->SetBinError(i+1,bin,h->GetBinError(i+1,j+1));
      
      h->SetBinContent(etaBins-i,j+1,h->GetBinContent(i+1,j+1));
      h->SetBinError(etaBins-i,j+1,h->GetBinError(i+1,j+1));
      
      h->SetBinContent(etaBins-i,bin,h->GetBinContent(i+1,j+1));
      h->SetBinError(etaBins-i,bin,h->GetBinError(i+1,j+1));
    }
  } 
}

#endif
