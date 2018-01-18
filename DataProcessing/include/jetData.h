#ifndef JETDATA_H
#define JETDATA_H

#include "TLorentzVector.h"

class jetData{
 public:
  static const Int_t nMaxJet = 500;
  static const int pwFlag = 6;

  int nref;
  Float_t jtpt[nMaxJet];
  Float_t jteta[nMaxJet];
  Float_t jtphi[nMaxJet];
  Int_t jtN[nMaxJet];
  Int_t jtNPW[nMaxJet][pwFlag];
  Float_t jtptFracPW[nMaxJet][pwFlag];

  TLorentzVector fourJet[nMaxJet];
};

#endif
