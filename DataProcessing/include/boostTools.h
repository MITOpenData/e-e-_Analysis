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
  float phi   = TMath::Pi()*2-angle;

  std::cout << a3.Mag() << " " << b3.Mag() << std::endl;
  float theta = TMath::ATan(TMath::Sin(phi)/(a3.Mag()/b3.Mag()+TMath::Cos(phi)));
  std::cout << a3.Mag() << " " << b3.Mag() << std::endl;
  std::cout << a3.Phi() << " "<< phi << " "  << theta << std::endl;
  TVector3 axis  = a3.Cross(b3).Unit();
  a3.Rotate(-theta,axis);
  std::cout << a3.Phi() << std::endl;
  TVector3 boost = a3.Unit();
  //boost magnitude
  float boostMag = (a3.Mag()*TMath::Cos(theta)+b3.Mag()*TMath::Cos(phi-theta))/(a.E()+b.E());
  boost.SetMag(-boostMag); 

  return boost;
}

void testCaseB2BBoost(){
  TLorentzVector a = TLorentzVector(2,-1.3,1,5);
  TLorentzVector b = TLorentzVector(-2,1,0,4);
  TVector3 boost = findBack2BackBoost(a,b);
  std::cout << boost.Mag() << std::endl;
  std::cout << boost.Phi() << std::endl;
  a.Boost(boost);
  b.Boost(boost);
  std::cout << boost.Phi()*360./2./TMath::Pi() << " " << a.Vect().Angle(b.Vect()) << std::endl;;
}

#endif

