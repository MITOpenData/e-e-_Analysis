#ifndef BOOSTTOOLS
#define BOOSTTOOLS
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "particleData.h"
#include "thrustTools.h"
#include "boostedEvtData.h"
#include <iostream>

TVector3 findBack2BackBoost(TLorentzVector a, TLorentzVector b){
  //find plane between 2 vectors
  TVector3 a3 = a.Vect();
  TVector3 b3 = b.Vect();
  float angle = a3.Angle(b3);
  float phi   = 2*TMath::Pi()-angle;

  //calculate angle from 1st vector that the boost vector is in the plane of a and b
  float theta = TMath::ATan(TMath::Sin(phi)/(a3.Mag()/b3.Mag()+TMath::Cos(phi)));
  TVector3 axis  = a3.Cross(b3).Unit();
  a3.Rotate(-theta,axis);
  TVector3 boost = a3.Unit();

  //boost magnitude in order to get the momenta to balanced out each other
  float boostMag = (a3.Mag()*TMath::Cos(theta)+b3.Mag()*TMath::Cos(phi-theta))/(a.E()+b.E());
  boost.SetMag(-boostMag); 

  return boost;
}


void setBoostedVariables(bool doBoost, particleData *p, boostedEvtData *b, TLorentzVector mainAxis = TLorentzVector(0,0,0,0), TVector3 boost = TVector3(0,0,0)){
  if(doBoost==0){
    b->WTAAxis_Theta = -999;
    b->WTAAxis_ThetaPerp = -999;
    b->WTAAxis_Phi = -999;
    b->boostx = 0; 
    b->boosty = 0; 
    b->boostz = 0; 
    b->boost = 0; 
 
    for(int i = 0; i< p->nParticle; i++){
      b->pt[i] = -999;
      b->pmag[i] = -999;
      b->eta[i] = -999;
      b->phi[i] = -999;
      b->theta[i] = -999;

      b->pt_Perp[i] = -999;
      b->pmag_Perp[i] = -999;
      b->eta_Perp[i] = -999;
      b->phi_Perp[i] = -999;
      b->theta_Perp[i] = -999;
    } 
  }
  else{
    mainAxis.Boost(boost);

    b->WTAAxis_Theta = mainAxis.Vect().Theta();
    b->WTAAxis_ThetaPerp = getPerpVector(mainAxis.Vect()).Theta();
    b->WTAAxis_Phi = mainAxis.Vect().Phi();

    b->boostx = boost.X(); 
    b->boosty = boost.Y(); 
    b->boostz = boost.Z(); 
    b->boost = boost.Mag(); 

    //    std::cout << "SETTING BOOST, " << p->nParticle << std::endl;

    for(int i = 0; i< p->nParticle; i++){
      TLorentzVector part;
      part.SetXYZM(p->px[i], p->py[i], p->pz[i], p->mass[i]);
      part.Boost(boost);
      
      //      std::cout << " " << p->pt[i] << ", " << p->phi[i] << ", " << p->eta[i] << std::endl;

      b->pt[i] = ptFromThrust(mainAxis.Vect().Unit(), part.Vect());
      b->pmag[i] = part.Vect().Mag();
      b->eta[i] = etaFromThrust(mainAxis.Vect().Unit(),part.Vect());
      b->theta[i] = thetaFromThrust(mainAxis.Vect().Unit(), part.Vect());
      b->phi[i] = phiFromThrust(mainAxis.Vect().Unit(), part.Vect());

      //      std::cout << "  " << b->pt[i] << ", " << b->phi[i] << ", " << b->eta[i] << ", " << b->theta[i] << std::endl;

      b->pt_Perp[i] = ptFromThrust(mainAxis.Vect(), part.Vect(), true);
      b->pmag_Perp[i] = part.Vect().Mag();
      b->eta_Perp[i] = etaFromThrust(mainAxis.Vect(), part.Vect(), true);
      b->theta_Perp[i] = thetaFromThrust(mainAxis.Vect(), part.Vect(), true);
      b->phi_Perp[i] = phiFromThrust(mainAxis.Vect(), part.Vect(), true);
    }
  }
}


//used for testing
void boostTools(){
  TLorentzVector a = TLorentzVector(4,0,-1,7);
  TLorentzVector b = TLorentzVector(-1,0.2,1,3);
  TVector3 boost = findBack2BackBoost(a,b);
  //  std::cout << "boost magnitude: " << boost.Mag() << std::endl;
  //  std::cout << "boost phi: " << boost.Phi()*360.0/2.0/TMath::Pi() << std::endl;
  a.Boost(boost);
  b.Boost(boost);
  //  std::cout << "post-boost dPhi: " << a.Vect().Angle(b.Vect()) << std::endl;
  //  std::cout << "post-boost difference in momentum (want 0): " << a.Vect().Mag() - b.Vect().Mag() << std::endl;
}

#endif

