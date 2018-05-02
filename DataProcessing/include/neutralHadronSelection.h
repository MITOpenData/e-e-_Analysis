#ifndef NeutralHadronSELECTION_H
#define NeutralHadronSELECTION_H

#include <vector>

#include "TMath.h"
#include "TLorentzVector.h"
#include "DataProcessing/include/particleData.h"

/*example do getting highPurity requirement for ith particle from pData:
NeutralHadronSelection trkSel = NeutralHadronSelection();
particleData pData;
int i;
bool isGoodTrk = NeutralHadronSelection.highPurity(&pData,i);
*/

class NeutralHadronSelection{
 public:
   NeutralHadronSelection();
   ~NeutralHadronSelection();

   void setAbsCosThCut(float cut);
   void setECut(float cut);
   
   inline bool passesAbsCosThCut(float theta);
   inline bool passesE(float p);
   inline bool passesPWFlag(Short_t pwflag);

   bool highPurity(particleData * p, int indx);
   void fillHighPurity(particleData * p); 
 private:
   float ECut = 0.4;
   float absCosThCut = 0.98;                       //maximum abs(cos(th)) of charged NeutralHadrons
};

bool NeutralHadronSelection::highPurity(particleData * p, int indx){
  if(!passesPWFlag(p->pwflag[indx])) return false;
  if(!passesAbsCosThCut(p->theta[indx]))   return false;
  if(!passesE(sqrt(p->pmag[indx]*p->pmag[indx]+p->mass[indx]*p->mass[indx])))         return false;

  return true;
}

void NeutralHadronSelection::fillHighPurity(particleData * p){
  for(int i = 0; i<p->nParticle; i++){
    if(highPurity(p,i)) p->highPurity[i] = true;
    else p->highPurity[i] = false;
  }
}

inline bool NeutralHadronSelection::passesAbsCosThCut(float theta){ return TMath::Abs(cos(theta))<=absCosThCut; }
inline bool NeutralHadronSelection::passesE(float e){ return e>=ECut;}
inline bool NeutralHadronSelection::passesPWFlag(Short_t pwflag){ return (pwflag==4||pwflag==5);}

void NeutralHadronSelection::setAbsCosThCut(float cut){ absCosThCut = cut;}
void NeutralHadronSelection::setECut(float cut){ ECut = cut;}

NeutralHadronSelection::NeutralHadronSelection(){}

NeutralHadronSelection::~NeutralHadronSelection(){}

#endif
