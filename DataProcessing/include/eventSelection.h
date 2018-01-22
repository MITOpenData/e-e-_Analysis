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

  Double_t masslessRescale_;

  Double_t pxSumFinal_;
  Double_t pySumFinal_;
  Double_t pzSumFinal_;
  Double_t eSumFinal_;

  Bool_t passesWW = false;

  eventSelection(){jDefWW = fastjet::JetDefinition(fastjet::ee_kt_algorithm);}
  eventSelection(particleData inPart){jDefWW = fastjet::JetDefinition(fastjet::ee_kt_algorithm); setEventSelection(inPart);}
  void setEventSelection(particleData inPart);
  Bool_t getPassesWW(){return passesWW;}
};

void eventSelection::setEventSelection(particleData inPart)
{
  particles.clear();
  passesWW = false;

  if(inPart.nParticle < 4) return;

  for(Int_t pI = 0; pI < inPart.nParticle; ++pI){
    Double_t e = TMath::Sqrt(inPart.pmag[pI]*inPart.pmag[pI] + inPart.mass[pI]*inPart.mass[pI]);
    particles.push_back(fastjet::PseudoJet(inPart.px[pI], inPart.py[pI], inPart.pz[pI], e));
  }

  fastjet::ClusterSequence cs(particles, jDefWW);
  std::vector<fastjet::PseudoJet> fourJetsFJ = cs.exclusive_jets(nJWW);
  std::vector<TLorentzVector> fourJets;

  for(unsigned int i = 0; i < fourJetsFJ.size(); ++i){
    TLorentzVector temp(fourJetsFJ.at(i).px(), fourJetsFJ.at(i).py(), fourJetsFJ.at(i).pz(), fourJetsFJ.at(i).E());
    fourJets.push_back(temp);
  }

  for(unsigned int i = 0; i < fourJets.size(); ++i){
    pxSumInit_ += fourJets.at(i).Px();
    pySumInit_ += fourJets.at(i).Py();
    pzSumInit_ += fourJets.at(i).Pz();
    eSumInit_ += fourJets.at(i).E();
  }

  masslessRescale_ = inPart.Energy/eSumInit_;

  for(unsigned int i = 0; i < fourJets.size(); ++i){fourJets.at(i) *= masslessRescale_;}

  for(unsigned int i = 0; i < fourJets.size(); ++i){
    pxSumFinal_ += fourJets.at(i).Px();
    pySumFinal_ += fourJets.at(i).Py();
    pzSumFinal_ += fourJets.at(i).Pz();
    eSumFinal_ += fourJets.at(i).E();
  }

  Double_t d2 = 999999999.;
  Double_t smallestAngle = 999999.;
  
  for(unsigned int i = 0; i < fourJets.size()-1; ++i){
    for(unsigned int j = i+1; j < fourJets.size(); ++j){
      Double_t tempAngle = fourJets.at(i).Angle(fourJets.at(j).Vect());

      if(tempAngle < smallestAngle) smallestAngle = tempAngle;
    }
  }
  
  Double_t cW = TMath::Cos(smallestAngle);

  if(cW < .9){
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
	d2Temp /= 80.4;
      
	if(d2Temp < d2) d2 = d2Temp;
      }
    }
  }

  passesWW = d2 >= .1 || cW >= .9;

  return;
}

#endif
