#ifndef TOOLS
#define TOOLS

#include "TMath.h"
#include "TH2F.h"
#include "TVector3.h"
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

inline double ptFromThrust(TVector3 thrust, TVector3 p){
  return p.Perp(thrust); 
}
inline double thetaFromThrust(TVector3 thrust, TVector3 p){
  return p.Angle(thrust);
}

//this is actually rapidity w/ pion mass assumption
inline double etaFromThrust(TVector3 thrust, TVector3 p){
  float pl = p*thrust.Unit();//logitudinal momentum component
  float E = TMath::Power(p.Mag2()+0.13957*0.13957,0.5);//energy w/ mass assumption
  return 0.5*TMath::Log((E+pl)/(E-pl));//rapidity
}
inline double phiFromThrust(TVector3 thrust, TVector3 p){
  TVector3 pt = p-((p*thrust.Unit())*(thrust.Unit()));//pt vector
  TVector3 z = TVector3(0,0,1);
  TVector3 phiOrigin = thrust.Unit().Cross((thrust.Unit().Cross(z)));//vector that will be phi=0 (in plane of thrust and beam line
  double phi = pt.Angle(phiOrigin);//get phi from 0 to pi

  //determine sign of phi based on cross product of pt and origin
  if( (phiOrigin.Cross(pt.Unit()))*thrust >= 0) return phi;
  else return -phi;
}

#endif
