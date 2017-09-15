#ifndef EVENTDATA_H
#define EVENTDATA_H

class eventData{
 public:
  Float_t missP;
  Float_t missPt;
  Float_t missTheta;
  Float_t missPhi;
  Float_t missChargedP;
  Float_t missChargedPt;
  Float_t missChargedTheta;
  Float_t missChargedPhi;
  Int_t nChargedHadrons;
  Int_t nChargedHadrons_GT0p4;
  Int_t nChargedHadrons_GT0p4Thrust;

  //thrust axis variables
  Float_t TTheta;
  Float_t TPhi;
  Float_t TTheta_charged;
  Float_t TPhi_charged;
};

#endif
