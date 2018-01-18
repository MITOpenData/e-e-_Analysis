#ifndef BOOSTTOOLS
#define BOOSTTOOLS
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include <iostream>

TVector3 findBack2BackBoost(TLorentzVector a, TLorentzVector b){
  //find plane between 2 vectors
  TVector3 a3 = a.Vect();
  TVector3 b3 = b.Vect();
  float angle = a3.Angle(b3);
  float phi   = 2*TMath::Pi()-angle;

  float theta = TMath::ATan(TMath::Sin(phi)/(a3.Mag()/b3.Mag()+TMath::Cos(phi)));
  TVector3 axis  = a3.Cross(b3).Unit();
  a3.Rotate(-theta,axis);
  TVector3 boost = a3.Unit();
  //boost magnitude
  float boostMag = (a3.Mag()*TMath::Cos(theta)+b3.Mag()*TMath::Cos(phi-theta))/(a.E()+b.E());
  boost.SetMag(-boostMag); 

  return boost;
}

void boostTools(){
  TLorentzVector a = TLorentzVector(4,0,-1,7);
  TLorentzVector b = TLorentzVector(-1,0.2,1,3);
  TVector3 boost = findBack2BackBoost(a,b);
  std::cout << "boost magnitude: " << boost.Mag() << std::endl;
  std::cout << "boost phi: " << boost.Phi()*360.0/2.0/TMath::Pi() << std::endl;
  a.Boost(boost);
  b.Boost(boost);
  std::cout << "post-boost dPhi: " << a.Vect().Angle(b.Vect()) << std::endl;
  std::cout << "post-boost difference in momentum (want 0): " << a.Vect().Mag() - b.Vect().Mag() << std::endl;
}

#endif

