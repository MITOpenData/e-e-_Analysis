#ifndef TOOLS
#define TOOLS

#include "TMath.h"
#include "TH2F.h"
#include "Settings.h"

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

void getLongRangeYield(Settings s, TH2F * h, TH1F * g){
  float binArea = (4*s.etaCut/(float)s.dEtaBins)*(2*TMath::Pi()/(float)s.dPhiBins); 
 
  for(int i = 0; i<s.dPhiBins; i++){
    float sum = 0; 
    float err2 = 0;
    for(int j = 0; j<s.dEtaBins; j++){
      float center = h->GetXaxis()->GetBinCenter(j+1); 
      if(TMath::Abs(center)>=s.dEtaRangeToIntegrate[0] && TMath::Abs(center)<s.dEtaRangeToIntegrate[1]){
        sum += h->GetBinContent(j+1,i+1);
        err2+= h->GetBinError(j+1,i+1)*h->GetBinError(j+1,i+1);
      }
    }
    g->SetBinContent(i+1,sum);
    g->SetBinError(i+1,TMath::Power(err2,0.5));
  }
                      //scale by bin area to get yield
  g->Scale(2*binArea);//scale by 2 because we only integrated the positive side, adn there is also yield in the negative eta side
}

#endif
