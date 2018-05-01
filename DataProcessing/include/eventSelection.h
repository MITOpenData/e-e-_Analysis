#ifndef EVENTSELECTION_H
#define EVENTSELECTION_H

#include <vector>

#include "TMath.h"
#include "TLorentzVector.h"

#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include "include/trackSelection.h"
#include "include/particleData.h"
#include "include/eventData.h"

class eventSelection{
 public:
  std::vector<fastjet::PseudoJet> particles;

  int nJISR = 2;
  fastjet::JetDefinition jDefISR;

  int nJWW = 4;
  fastjet::JetDefinition jDefWW;

  Double_t pxSumInit_;
  Double_t pySumInit_;
  Double_t pzSumInit_;
  Double_t eSumInit_;

  std::vector<TLorentzVector> initFourJet;

  Double_t masslessRescale_;

  Double_t factorA_;
  Double_t factorB_;
  Double_t factorC_;
  Double_t factorD_;

  Double_t pxSumFinal_;
  Double_t pySumFinal_;
  Double_t pzSumFinal_;
  Double_t eSumFinal_;

  std::vector<TLorentzVector> finalFourJet;

  Bool_t passesTotalChgEnergyMin = false;
  Bool_t passesNTrkMin = false;
  Bool_t passesSTheta = false;
  Bool_t passesMissP = false;

  Bool_t passesISR = false;
  Bool_t passesWW = false;
  Bool_t passesD2 = false;
  Bool_t passesCW = false;

  Float_t d2 = 9999999.;
  Float_t cW = 9999999.;

  eventSelection(){
    jDefISR = fastjet::JetDefinition(fastjet::ee_kt_algorithm);
    jDefWW = fastjet::JetDefinition(fastjet::ee_kt_algorithm);
  }
  eventSelection(particleData* inPart, eventData* inData){
    jDefISR = fastjet::JetDefinition(fastjet::ee_kt_algorithm);
    jDefWW = fastjet::JetDefinition(fastjet::ee_kt_algorithm);
    setEventSelection(inPart, inData);
  }
  Bool_t checkDoRescale(std::vector<TLorentzVector>* inJets, std::vector<double>* factors, Float_t cmEnergy);
  void setEventSelection(particleData* inPart, eventData* inData);
  Bool_t getPassesTotalChgEnergyMin(){return passesTotalChgEnergyMin;}
  Bool_t getPassesNTrkMin(){return passesNTrkMin;}
  Bool_t getPassesSTheta(){return passesSTheta;}
  Bool_t getPassesMissP(){return passesMissP;}

  Bool_t getPassesISR(){return passesISR;}
  Bool_t getPassesWW(){return passesWW;}

  private:
    TrackSelection trkSel = TrackSelection();

};


Bool_t eventSelection::checkDoRescale(std::vector<TLorentzVector>* inJets, std::vector<double>* factors, Float_t cmEnergy)
{
  Float_t px0 = inJets->at(0).Px();
  Float_t py0 = inJets->at(0).Py();
  Float_t pz0 = inJets->at(0).Pz();
  Float_t p0 = TMath::Sqrt(px0*px0 + py0*py0 + pz0*pz0);

  Float_t px1 = inJets->at(1).Px();
  Float_t py1 = inJets->at(1).Py();
  Float_t pz1 = inJets->at(1).Pz();
  Float_t p1 = TMath::Sqrt(px1*px1 + py1*py1 + pz1*pz1);
  
  Float_t px2 = inJets->at(2).Px();
  Float_t py2 = inJets->at(2).Py();
  Float_t pz2 = inJets->at(2).Pz();
  Float_t p2 = TMath::Sqrt(px2*px2 + py2*py2 + pz2*pz2);

  Float_t px3 = inJets->at(3).Px();
  Float_t py3 = inJets->at(3).Py();
  Float_t pz3 = inJets->at(3).Pz();
  Float_t p3 = TMath::Sqrt(px3*px3 + py3*py3 + pz3*pz3);

  Float_t p10xy = py1 - py0*px1/px0;
  Float_t p20xy = py2 - py0*px2/px0;
  Float_t p30xy = py3 - py0*px3/px0;

  Float_t ct = pz2 - pz1*p20xy/p10xy - pz0*px2/px0 + pz0*px1*p20xy/(px0*p10xy);
  Float_t dt = pz3 - pz1*p30xy/p10xy - pz0*px3/px0 + pz0*px1*p30xy/(px0*p10xy);

  Float_t d = 1.;
  Float_t c = -d*dt/ct;
  Float_t b = -(c*p20xy + d*p30xy)/p10xy;
  Float_t a = -1.*(b*px1 + c*px2 + d*px3)/px0;

  d = 1.;
  c = 1.;
  b = 1.;
  a = 1.;

  p0 = a*TMath::Sqrt(px0*px0 + py0*py0 + pz0*pz0);
  p1 = b*TMath::Sqrt(px1*px1 + py1*py1 + pz1*pz1);
  p2 = c*TMath::Sqrt(px2*px2 + py2*py2 + pz2*pz2);
  p3 = d*TMath::Sqrt(px3*px3 + py3*py3 + pz3*pz3);

  Float_t cmFactor = cmEnergy/(p0 + p1 + p2 + p3);
  p0 *= cmFactor;
  p1 *= cmFactor;
  p2 *= cmFactor;
  p3 *= cmFactor;
  
  a *= cmFactor;
  b *= cmFactor;
  c *= cmFactor;
  d *= cmFactor;

  inJets->at(0).SetPxPyPzE(a*px0, a*py0, a*pz0, p0);
  inJets->at(1).SetPxPyPzE(b*px1, b*py1, b*pz1, p1);
  inJets->at(2).SetPxPyPzE(c*px2, c*py2, c*pz2, p2);
  inJets->at(3).SetPxPyPzE(d*px3, d*py3, d*pz3, p3);

  factors->push_back(a);
  factors->push_back(b);
  factors->push_back(c);
  factors->push_back(d);
  //  std::cout << " a,b,c,d: " << a << ", " << b << ", " << c << ", " << d << std::endl;

  return false;
}

void eventSelection::setEventSelection(particleData* inPart, eventData* inData)
{
  particles.clear();

  passesTotalChgEnergyMin = false;
  passesNTrkMin = false;
  passesSTheta = false;
  passesMissP = false;

  passesISR = false;
  passesWW = false;
  passesD2 = false;
  passesCW = false;

  initFourJet.clear();
  finalFourJet.clear();

  initFourJet.push_back({0,0,0,0});
  initFourJet.push_back({0,0,0,0});
  initFourJet.push_back({0,0,0,0});
  initFourJet.push_back({0,0,0,0});

  finalFourJet.push_back({0,0,0,0});
  finalFourJet.push_back({0,0,0,0});
  finalFourJet.push_back({0,0,0,0});
  finalFourJet.push_back({0,0,0,0});

  if(inPart->nParticle < 4) return;

  Float_t TotalChgEnergy = 0;
  Int_t NTrk = 0;
  Float_t STheta = inData->STheta;
  Float_t MissP = inData->missP;

  for(Int_t pI = 0; pI < inPart->nParticle; ++pI){
    Double_t e = TMath::Sqrt(inPart->pmag[pI]*inPart->pmag[pI] + inPart->mass[pI]*inPart->mass[pI]);
    particles.push_back(fastjet::PseudoJet(inPart->px[pI], inPart->py[pI], inPart->pz[pI], e));

    if(!trkSel.highPurity(inPart, pI)) continue;
    TotalChgEnergy += TMath::Sqrt(inPart->pmag[pI]*inPart->pmag[pI] + inPart->mass[pI]*inPart->mass[pI]);
    NTrk += 1;
  }

  passesTotalChgEnergyMin = TotalChgEnergy >= 15;
  passesNTrkMin = NTrk >= 5;
  passesSTheta = TMath::Abs(TMath::Cos(STheta)) <= .82;
  passesMissP = MissP < 20;

  fastjet::ClusterSequence csISR(particles, jDefISR);
  std::vector<fastjet::PseudoJet> twoJetsFJ = csISR.exclusive_jets(nJISR);
  std::vector<TLorentzVector> twoJets;

  //two jet implementation follows https://arxiv.org/pdf/hep-ex/9810047.pdf, 1,2,3 eq.
  for(unsigned int i = 0; i < twoJetsFJ.size(); ++i){
    TLorentzVector temp(twoJetsFJ.at(i).px(), twoJetsFJ.at(i).py(), twoJetsFJ.at(i).pz(), twoJetsFJ.at(i).E());
    twoJets.push_back(temp);
  }
  Double_t absBeta = TMath::Abs(TMath::Sin(twoJets.at(0).Theta() + twoJets.at(1).Theta()));
  absBeta /= (TMath::Sin(twoJets.at(0).Theta()) + TMath::Abs(twoJets.at(1).Theta()));
  Double_t x = 2*absBeta/(1+absBeta);
  Double_t sPrime = inPart->Energy*inPart->Energy*(1-x);
  Double_t mass = (twoJets.at(0) + twoJets.at(1)).M();

  //cuts from http://cds.cern.ch/record/690637/files/ep-2003-084.pdf, page2, middle paragraph (before d^2 equation
  passesISR = mass/inPart->Energy > .7 || sPrime/(inPart->Energy*inPart->Energy) > .81;

  if(!passesISR) return;

  fastjet::ClusterSequence csWW(particles, jDefWW);
  std::vector<fastjet::PseudoJet> fourJetsFJ = csWW.exclusive_jets(nJWW);
  std::vector<TLorentzVector> fourJets;

  for(unsigned int i = 0; i < fourJetsFJ.size(); ++i){
    TLorentzVector temp(fourJetsFJ.at(i).px(), fourJetsFJ.at(i).py(), fourJetsFJ.at(i).pz(), fourJetsFJ.at(i).E());
    fourJets.push_back(temp);
    initFourJet.at(i) = temp;
  }

  for(unsigned int i = 0; i < fourJets.size(); ++i){
    for(unsigned int j = i+1; j < fourJets.size(); ++j){
      if(fourJets.at(i).E() > fourJets.at(j).E()){
	TLorentzVector temp = fourJets.at(j);
	fourJets.at(j) = fourJets.at(i);
	fourJets.at(i) = temp;

	initFourJet.at(j) = initFourJet.at(i);
	initFourJet.at(i) = temp;
      }
    }
  }

  pxSumInit_ = 0;
  pySumInit_ = 0;
  pzSumInit_ = 0;
  eSumInit_ = 0;

  pxSumFinal_ = 0;
  pySumFinal_ = 0;
  pzSumFinal_ = 0;
  eSumFinal_ = 0;
  
  for(unsigned int i = 0; i < fourJets.size(); ++i){
    pxSumInit_ += fourJets.at(i).Px();
    pySumInit_ += fourJets.at(i).Py();
    pzSumInit_ += fourJets.at(i).Pz();
    eSumInit_ += fourJets.at(i).E();
  }
  
  bool doRescale = (TMath::Abs(pxSumInit_) > 1. || TMath::Abs(pySumInit_) > 1. || TMath::Abs(pzSumInit_) > 1. || TMath::Abs(eSumInit_ - inPart->Energy) > 1.);
  std::vector<double> factors;

  factorA_ = 1.;
  factorB_ = 1.;
  factorC_ = 1.;
  factorD_ = 1.;
  
  while(doRescale){
    doRescale = checkDoRescale(&fourJets, &factors, inPart->Energy);
    factorA_ = factors.at(0);
    factorB_ = factors.at(1);
    factorC_ = factors.at(2);
    factorD_ = factors.at(3);
  }


  /*
  Float_t pxSumMid_ = 0;
  Float_t pySumMid_ = 0;
  Float_t pzSumMid_ = 0;
  Float_t eSumMid_ = 0;
  for(unsigned int i = 0; i < fourJets.size(); ++i){
    pxSumMid_ += fourJets.at(i).Px();
    pySumMid_ += fourJets.at(i).Py();
    pzSumMid_ += fourJets.at(i).Pz();
    eSumMid_ += fourJets.at(i).E();
  }

  masslessRescale_ = inPart->Energy/eSumMid_;
  //  std::cout << pxSumMid_ << ", " << pySumMid_ << ", " << pzSumMid_ << std::endl;

  for(unsigned int i = 0; i < fourJets.size(); ++i){fourJets.at(i) *= masslessRescale_;}
  */

  for(unsigned int i = 0; i < fourJets.size(); ++i){
    pxSumFinal_ += fourJets.at(i).Px();
    pySumFinal_ += fourJets.at(i).Py();
    pzSumFinal_ += fourJets.at(i).Pz();
    eSumFinal_ += fourJets.at(i).E();
    finalFourJet.at(i) = fourJets.at(i);
  }

  d2 = 999999999.;
  Double_t smallestAngle = 999999.;
  
  for(unsigned int i = 0; i < fourJets.size()-1; ++i){
    for(unsigned int j = i+1; j < fourJets.size(); ++j){
      Double_t tempAngle = fourJets.at(i).Angle(fourJets.at(j).Vect());

      if(tempAngle < smallestAngle) smallestAngle = tempAngle;
    }
  }
  
  cW = TMath::Cos(smallestAngle);

  std::vector<double> twoJetMasses;
  for(unsigned int i = 0; i < fourJets.size()-1; ++i){
    for(unsigned int j = i+1; j < fourJets.size(); ++j){
      TLorentzVector temp = fourJets.at(i) + fourJets.at(j);
      twoJetMasses.push_back(temp.M());
    }
  }

  for(unsigned int i = 0; i < twoJetMasses.size()-1; ++i){
    for(unsigned int j = i+1; j < twoJetMasses.size(); ++j){
      Double_t d2Temp = (twoJetMasses.at(i) - 80.4)*(twoJetMasses.at(i) - 80.4);
      d2Temp += (twoJetMasses.at(j) - 80.4)*(twoJetMasses.at(j) - 80.4);
      d2Temp /= (80.4*80.4);
      
      if(d2Temp < d2) d2 = d2Temp;
    }
  }

  passesD2 = d2 >= .1;
  passesCW = cW >= .9;
  passesWW = passesD2 || passesCW;

  return;
}

#endif
