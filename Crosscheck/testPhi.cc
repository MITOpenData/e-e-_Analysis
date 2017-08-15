#include "Tools.h"
#include "Settings.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TAxis.h"
#include "TVector3.h"
#include <iostream>
#include <cmath>
#include <vector>

//code for checking kinematics w/ respect to thrust axis
//throws particles on a uniform sphere and then calculates the kinematics on a new random thrust axis

void testPhi(){

  TRandom2 randGen = TRandom2();
  std::vector< TVector3 > v;

  TFile * f = TFile::Open("testPhi.root","recreate");

  TH1D * hphi = new TH1D("phi","phi",100,0,2*TMath::Pi());
  TH1D * htheta = new TH1D("theta","theta",100,0,TMath::Pi());
  TH1D * hphi2 = new TH1D("phi2","phi2",100,-TMath::Pi(),TMath::Pi());
  TH1D * htheta2 = new TH1D("theta2","theta2",100,0,TMath::Pi());

  for(int i = 0; i<10000; i++){
    if(i%1000==0) std::cout << i << std::endl;
    TVector3 temp = TVector3(0,0,0);
    float phi = randGen.Rndm()*TMath::Pi()*2;
    std::cout << phi << std::endl;
    hphi->Fill(phi);    


    float theta = TMath::ACos(randGen.Rndm()*2-1);
    std::cout << theta << std::endl;
    htheta->Fill(theta);
    temp.SetMagThetaPhi(1,theta,phi);
    v.push_back(temp);
  }

  float tPhi = randGen.Rndm()*TMath::Pi()*2;
  float tTheta = TMath::ACos(randGen.Rndm()*2-1);
  TVector3 thrust = TVector3(0,0,0);
  thrust.SetMagThetaPhi(1,tTheta,tPhi);

  std::cout << tTheta << " " << tPhi << std::endl;
  for(int i = 0; i<10000; i++){
    TVector3 p = TVector3(v.at(i).Px(),v.at(i).Py(),v.at(i).Pz());
    hphi2->Fill(phiFromThrust(thrust,p));
    htheta2->Fill(thetaFromThrust(thrust,p));
  }
  f->Write();

}
