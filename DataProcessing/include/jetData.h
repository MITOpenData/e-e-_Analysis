#ifndef JETDATA_H
#define JETDATA_H

class jetData{
 public:
  static const Int_t nMaxJet = 500;

  int nref;
  Float_t jtpt[nMaxJet];
  Float_t jteta[nMaxJet];
  Float_t jtphi[nMaxJet];
};

#endif
