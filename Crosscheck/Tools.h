#ifndef TOOLS
#define TOOLS

#include "TMath.h"
#include "TH2F.h"
#include "TVector3.h"
#include "Settings.h"
#include "TH2D.h"

inline float dPhi(float phi1, float phi2){
  return TMath::Abs(TMath::ACos(TMath::Cos(phi1-phi2)));
}

inline float dEta(float eta1, float eta2){
  return TMath::Abs(eta1-eta2);
}

void fillAllQuadrants(TH2F * h, float deta, float dphi, float weight){
  h->Fill(deta,dphi,weight); 
  h->Fill(-deta,dphi,weight);
  if(-dphi<-TMath::Pi()/2.0) dphi -=TMath::Pi()*2; 
  h->Fill(deta,-dphi,weight); 
  h->Fill(-deta,-dphi,weight); 
}


void getLongRangeYield(Settings s, TH2F * h, TH1F * g){
  float binArea = (4*s.etaPlotRange/(float)s.dEtaBins)*(2*TMath::Pi()/(float)s.dPhiBins); 
 
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

float getEff(Settings s, float pt, float eta){
  if(s.experiment == 0){
    return 1;
  }
  if(s.experiment == 3){
    return 1;
  }
  return 1;
}


#endif
