#ifndef PARTICLEDATA_H
#define PARTICLEDATA_H

class particleData{
 public:
  static const int nMaxPart = 10000;

  int nParticle;
  int EventNo;
  int RunNo;
  int year;
  int process;
  float Energy;
  int bFlag;
  float bx;
  float by;
  float ebx;
  float eby;
  float px[nMaxPart];
  float py[nMaxPart];
  float pz[nMaxPart];
  float pt[nMaxPart];
  float pmag[nMaxPart];
  float eta[nMaxPart];
  float theta[nMaxPart];
  float phi[nMaxPart];
  float mass[nMaxPart];
  float charge[nMaxPart];
  // Starting from 0, pwflag (via Marcello) - CHARGED_TRACK, CHARGED_LEPTONS1, CHARGED_LEPTONS2, V0, PHOTON, NEUTRAL_HADRON
  int pwflag[nMaxPart];
  int pid[nMaxPart];
  float d0[nMaxPart];
  float z0[nMaxPart];
  int ntpc[nMaxPart];
  
  //thrust axis variables
  float pt_wrtThr[nMaxPart];
  float eta_wrtThr[nMaxPart];
  float theta_wrtThr[nMaxPart];
  float phi_wrtThr[nMaxPart];
  float pt_wrtChThr[nMaxPart];
  float eta_wrtChThr[nMaxPart];
  float theta_wrtChThr[nMaxPart];
  float phi_wrtChThr[nMaxPart];
};

#endif
