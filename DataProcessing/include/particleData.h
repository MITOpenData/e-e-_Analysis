#ifndef PARTICLEDATA_H
#define PARTICLEDATA_H

class particleData{
 public:
  static const Int_t nMaxPart = 100000;

  int nParticle;
  int EventNo;
  int RunNo;
  float Energy;
  Float_t px[nMaxPart];
  Float_t py[nMaxPart];
  Float_t pz[nMaxPart];
  Float_t pt[nMaxPart];
  Float_t pmag[nMaxPart];
  Float_t eta[nMaxPart];
  Float_t theta[nMaxPart];
  Float_t phi[nMaxPart];
  Float_t mass[nMaxPart];
  Float_t charge[nMaxPart];
  Int_t pwflag[nMaxPart];
  Int_t pid[nMaxPart];
};

#endif
