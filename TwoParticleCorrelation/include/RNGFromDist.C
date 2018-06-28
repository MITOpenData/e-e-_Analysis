#include "TRandom3.h"
#include "TF1.h"
#include "TMath.h"

class RNGFromDist
{
  public:
    RNGFromDist(TF1* f);
    ~RNGFromDist();
    float getRand();

  private:
    TRandom3 * rng3;
    TF1 * f;
    float normalization;

};

float RNGFromDist::getRand(){
  while(true){
    float a = rng3->Rndm()*2*TMath::Pi()-TMath::Pi()/2;
    float b = rng3->Rndm()*normalization;
  
    if(b<f->Eval(a)) return a;
  }
}

RNGFromDist::RNGFromDist(TF1* f1){
  rng3 = new TRandom3();
  f = f1;
  normalization = f->GetMaximum(-TMath::Pi()/2,3*TMath::Pi()/2.0);
}

RNGFromDist::~RNGFromDist(){
  delete rng3;
}

