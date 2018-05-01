/**************************************************************************************/
// Common Utilities
/**************************************************************************************/
#include <TCanvas.h>
#include <TVector3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>
#include <TMath.h>

#include "Selection.h"

#define PI 3.14159265358979


// Calculate Delta Phi
double dphi(double phi1,double phi2)
{
    double a=phi1-phi2;
    while (a<-PI) a+=2*PI;
    while (a>PI) a-=2*PI;
    
    if (a<-PI/2) a=2*PI+a;
    return a;
}

/**************************************************************************************/
// Austin Thrust Angle Rotators
/**************************************************************************************/

//this is actually rapidity w/ pion mass assumption
inline double ptFromThrust(TVector3 thrust, TVector3 p){
  return p.Perp(thrust); 
}

inline double thetaFromThrust(TVector3 thrust, TVector3 p){
  return p.Angle(thrust);
}

inline double etaFromThrust(TVector3 thrust, TVector3 p){
  float pl = p*thrust.Unit();//logitudinal momentum component
  float E = TMath::Power(p.Mag2()+0.13957*0.13957,0.5);//energy w/ mass assumption
  return 0.5*TMath::Log((E+pl)/(E-pl));//rapidity
}

inline double phiFromThrust(TVector3 thrust, TVector3 p){
  TVector3 pt = p-((p*thrust.Unit())*thrust.Unit());//pt vector
  TVector3 z = TVector3(0,0,1);
  TVector3 phiOrigin = thrust.Unit().Cross((thrust.Unit().Cross(z)));//vector that will be phi=0 (in plane of thrust and beam line
  double phi = pt.Angle(phiOrigin);//get phi from 0 to pi
  
  //determine sign of phi based on cross product of pt and origin
  if( (phiOrigin.Cross(pt.Unit()))*thrust >= 0) return phi;
  else return -phi;
}

/**************************************************************************************/
// BELLE Particle Definition 
/**************************************************************************************/
enum SIMPLEPID {BELLE_PHOTON, BELLE_ELECTRON, BELLE_PION, BELLE_MUON, BELLE_KAON, BELLE_PROTON};

/**************************************************************************************/
// ALEPH Particle Flow Classification
/**************************************************************************************/
enum SIMPLEPWFLAG {ALEPH_CHARGED_TRACK, ALEPH_CHARGED_LEPTONS1, ALEPH_CHARGED_LEPTONS2, ALEPH_V0, ALEPH_PHOTON, ALEPH_NEUTRAL_HADRON};


/**************************************************************************************/
// Calculate 2D correlation function ratios
/**************************************************************************************/
int calculateRatio(TH2F *h_2D, TH2F *h_2Dmix, TH2F *h_ratio)
{
    
    //double normalization = h_2D->GetXaxis()->GetBinWidth(1)*h_2D->GetYaxis()->GetBinWidth(1);
    double ratio = 0;
    double errrel_ratio = 0;
    double errrel_num = 0;
    double errrel_den = 0;

    double b00_x=h_2Dmix->GetXaxis()->FindBin(0.);
    double b00_y=h_2Dmix->GetYaxis()->FindBin(0.);
    double B00=h_2Dmix->GetBinContent(b00_x,b00_y);
    //double errrel_B00=h_2Dmix->GetBinError(b00_x,b00_y)/B00;
    
    std::cout<<"value of B(0,0)="<<B00<<std::endl;
    std::cout<<"x axis "<<h_2Dmix->GetXaxis()->GetBinCenter(b00_x)<<std::endl;
    std::cout<<"y axis "<<h_2Dmix->GetYaxis()->GetBinCenter(b00_y)<<std::endl;
    
    for (Int_t x=0;x<=h_2D->GetNbinsX();x++){        
        for (Int_t y=0;y<=h_2D->GetNbinsY();y++){
            if(h_2Dmix->GetBinContent(x,y)>0){
                
                ratio=B00*(h_2D->GetBinContent(x,y)/h_2Dmix->GetBinContent(x,y));
                errrel_num=h_2D->GetBinError(x,y)/h_2D->GetBinContent(x,y);
                errrel_den=h_2Dmix->GetBinError(x,y)/h_2Dmix->GetBinContent(x,y);
                // Yen-Jie: Take out the error of B00 for the moment since this is a global uncertainty
                errrel_ratio=TMath::Sqrt(errrel_num*errrel_num+errrel_den*errrel_den);   // +errrel_B00*errrel_B00
                //h_ratio->SetBinContent(x,y,ratio/normalization);
                h_ratio->SetBinContent(x,y,ratio);
                h_ratio->SetBinError(x,y,errrel_ratio*ratio);
                //h_ratio->SetBinError(x,y,errrel_ratio*ratio/normalization);
            } else {
	      h_ratio->SetBinContent(x,y,0);
	      h_ratio->SetBinError(x,y,0);
	    }
        }
    }
    
    return 0;
}

void getLongRangeYield(Selection s, TH2F * h, TH1F * g){
    
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
    g->Scale(2*s.getDifferential());//scale by 2 because we only integrated the positive side, adn there is also yield in the negative eta side
}
