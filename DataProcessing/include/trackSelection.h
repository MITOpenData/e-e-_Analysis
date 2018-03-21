#ifndef TRACKSELECTION_H
#define TRACKSELECTION_H

#include <vector>

#include "TMath.h"
#include "TLorentzVector.h"
#include "include/particleData.h"

/*example do getting highPurity requirement for ith particle from pData:
TrackSelection trkSel = TrackSelection();
particleData pData;
int i;
bool isGoodTrk = TrackSelection.highPurity(&pData,i);
*/

class TrackSelection{
 public:
   TrackSelection();
   ~TrackSelection();

   void setnTPCCut(Short_t cut);
   void setThetaCutLow(float cut);
   void setThetaCutHigh(float cut);
   void setPCut(float cut);
   void setPtCut(float cut);
   void setD0Cut(float cut);
   void setZ0Cut(float cut);

   inline bool passesNTPC(Short_t ntpc);
   inline bool passesTheta(float theta);
   inline bool passesP(float p);
   inline bool passesPt(float pt);
   inline bool passesD0(float d0);
   inline bool passesZ0(float z0);
   inline bool passesPWFlag(Short_t pwflag);

   bool highPurity(particleData * p, int indx);
 
 private:
   Short_t nTPCcut = 4;
   float thetaCutLow = 20.*TMath::Pi()/180.;
   float thetaCutHigh = 160.*TMath::Pi()/180.;
   float pCut = 0.2;//currently not used in highPurity
   float ptCut = 0.2;
   float d0Cut = 3;
   float z0Cut = 5;
};

bool TrackSelection::highPurity(particleData * p, int indx){
  if(!passesPWFlag(p->pwflag[indx])) return false;
  if(!passesTheta(p->theta[indx]))   return false;
  if(!passesP(p->pmag[indx]))        return false;
  if(!passesD0(p->d0[indx]))         return false;
  if(!passesZ0(p->z0[indx]))         return false;
  if(!passesNTPC(p->ntpc[indx]))     return false;

  return true;

}

inline bool TrackSelection::passesNTPC(Short_t ntpc){ return ntpc>=nTPCcut; }
inline bool TrackSelection::passesTheta(float theta){ return (theta>=thetaCutLow) && (theta<=thetaCutHigh); }
inline bool TrackSelection::passesP(float p){ return p>=pCut;}
inline bool TrackSelection::passesPt(float pt){ return pt>=ptCut;}
inline bool TrackSelection::passesD0(float d0){ return TMath::Abs(d0)<=d0Cut;}
inline bool TrackSelection::passesZ0(float z0){ return TMath::Abs(z0)<=z0Cut;}
inline bool TrackSelection::passesPWFlag(Short_t pwflag){ return pwflag==0;}

void TrackSelection::setnTPCCut(Short_t cut){ nTPCcut = cut;}
void TrackSelection::setThetaCutLow(float cut){ thetaCutLow = cut;}
void TrackSelection::setThetaCutHigh(float cut){ thetaCutHigh = cut;}
void TrackSelection::setPCut(float cut){ pCut = cut;}
void TrackSelection::setPtCut(float cut){ ptCut = cut;}
void TrackSelection::setD0Cut(float cut){ d0Cut = cut;}
void TrackSelection::setZ0Cut(float cut){ z0Cut = cut;}

TrackSelection::TrackSelection(){}

TrackSelection::~TrackSelection(){}

#endif
