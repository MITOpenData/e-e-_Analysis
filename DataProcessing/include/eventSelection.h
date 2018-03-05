#ifndef EVENTSELECTION_H
#define EVENTSELECTION_H

#include <vector>

#include "TMath.h"
#include "TLorentzVector.h"

#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include "include/particleData.h"

class eventSelection{
 public:
  std::vector<fastjet::PseudoJet> particles;
  int nJWW = 4;
  fastjet::JetDefinition jDefWW;

  Double_t pxSumInit_;
  Double_t pySumInit_;
  Double_t pzSumInit_;
  Double_t eSumInit_;

  std::vector<TLorentzVector> initFourJet;

  Double_t masslessRescale_;

  Double_t pxSumFinal_;
  Double_t pySumFinal_;
  Double_t pzSumFinal_;
  Double_t eSumFinal_;

  std::vector<TLorentzVector> finalFourJet;

  Bool_t passesWW = false;
  Bool_t passesD2 = false;
  Bool_t passesCW = false;

  Float_t d2 = 9999999.;
  Float_t cW = 9999999.;

  eventSelection(){jDefWW = fastjet::JetDefinition(fastjet::ee_kt_algorithm);}
  eventSelection(particleData* inPart){jDefWW = fastjet::JetDefinition(fastjet::ee_kt_algorithm); setEventSelection(inPart);}
  Bool_t checkDoRescale(std::vector<TLorentzVector>* inJets, Float_t eScale);
  void setEventSelection(particleData* inPart);
  Bool_t getPassesWW(){return passesWW;}
  Bool_t getPassesD2(){return passesD2;}
  Bool_t getPassesCW(){return passesCW;}
  Float_t getD2(){return d2;}
  Float_t getCW(){return cW;}

  Double_t getInitFourJetSumPx(){return pxSumInit_;}
  Double_t getInitFourJetSumPy(){return pySumInit_;}
  Double_t getInitFourJetSumPz(){return pzSumInit_;}
  Double_t getInitFourJetSumE(){return eSumInit_;}

  std::vector<TLorentzVector> getInitFourJet(){return initFourJet;}

  Double_t getFinalFourJetSumPx(){return pxSumFinal_;}
  Double_t getFinalFourJetSumPy(){return pySumFinal_;}
  Double_t getFinalFourJetSumPz(){return pzSumFinal_;}
  Double_t getFinalFourJetSumE(){return eSumFinal_;}

  std::vector<TLorentzVector> getFinalFourJet(){return finalFourJet;}
};


/*
Bool_t eventSelection::checkDoRescale(std::vector<TLorentzVector>* inJets, Float_t eScale)
{
  Float_t currPx = 0;
  Float_t currPy = 0;
  Float_t currPz = 0;
  Float_t currE = 0;

  for(unsigned int i = 0 ; i < inJets.size(); ++i){
    currPx += inJets.at(i).Px();
    currPy += inJets.at(i).Py();
    currPz += inJets.at(i).Pz();

    inJets.at(i).SetE(TMath::Sqrt(inJets.at(i).Px()*inJets.at(i).Px() + inJets.at(i).Py()*inJets.at(i).Py() + inJets.at(i).Pz()*inJets.at(i).Pz()));
    
    currE += inJets.at(i).E();
  }

  if(TMath::Abs(currPx) > 1.){
    
    return true;
  }
  else if(TMath::Abs(currPy) > 1.){

    return true;
  }
  else if(TMath::Abs(currPz) > 1.){

    return true;
  }
  else if(TMath::Abs(currE) > 1.){

    return true;
  }

  return false;
}
*/

Bool_t eventSelection::checkDoRescale(std::vector<TLorentzVector>* inJets, Float_t eScale)
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

  Float_t d = eScale;
  Float_t divD = -p0*(p20xy*px1*dt/(p10xy*ct) - px1*p30xy/p10xy - px2*dt/ct + px3)/px0;
  divD += p1*(p20xy*dt/(p10xy*ct) - p30xy/p10xy) - p2*dt/ct + p3;

  d /= divD;

  Float_t c = -d*dt/ct;
  Float_t b = -(c*p20xy + d*p30xy)/p10xy;
  Float_t a = -1.*(b*px1 + c*px2 + d*px3)/px0;

  inJets->at(0).SetPxPyPzE(a*px0, a*py0, a*pz0, a*p0);
  inJets->at(1).SetPxPyPzE(b*px1, b*py1, b*pz1, b*p1);
  inJets->at(2).SetPxPyPzE(c*px2, c*py2, c*pz2, c*p2);
  inJets->at(3).SetPxPyPzE(d*px3, d*py3, d*pz3, d*p3);

  return false;
}

void eventSelection::setEventSelection(particleData* inPart)
{
  particles.clear();
  passesWW = false;
  passesD2 = false;
  passesCW = false;

  if(inPart->nParticle < 4) return;

  for(Int_t pI = 0; pI < inPart->nParticle; ++pI){
    Double_t e = TMath::Sqrt(inPart->pmag[pI]*inPart->pmag[pI] + inPart->mass[pI]*inPart->mass[pI]);
    particles.push_back(fastjet::PseudoJet(inPart->px[pI], inPart->py[pI], inPart->pz[pI], e));
  }

  fastjet::ClusterSequence cs(particles, jDefWW);
  std::vector<fastjet::PseudoJet> fourJetsFJ = cs.exclusive_jets(nJWW);
  std::vector<TLorentzVector> fourJets;

  initFourJet.clear();
  finalFourJet.clear();

  for(unsigned int i = 0; i < fourJetsFJ.size(); ++i){
    TLorentzVector temp(fourJetsFJ.at(i).px(), fourJetsFJ.at(i).py(), fourJetsFJ.at(i).pz(), fourJetsFJ.at(i).E());
    fourJets.push_back(temp);
    initFourJet.push_back(temp);
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
  
  bool doRescale = (TMath::Abs(pxSumInit_) > 1. || TMath::Abs(pySumInit_) > 1. || TMath::Abs(pzSumInit_) > 1.);

  while(doRescale){doRescale = checkDoRescale(&fourJets, inPart->Energy);}

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
    finalFourJet.push_back(fourJets.at(i));
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
