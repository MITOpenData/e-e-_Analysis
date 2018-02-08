#ifndef SPHERETOOLS
#define SPHERETOOLS
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include <iostream>

class Sphericity{

  public:
    Sphericity(int n, float *px, float *py, float *pz);
    ~Sphericity();

    inline TVector3 sphericityAxis();
    inline TVector3 linSphericityAxis();
    inline float sphericity(); 
    inline float aplanarity(); 
    inline float linSphericity(); 
    inline float linAplanarity(); 
    inline float linC(); 
    inline float linD(); 

  private:
    float l1, l2, l3;   //eigenvalues
    TVector3 v1, v2, v3;//eigenvectors

    //linearized quantities
    float linl1, linl2, linl3;  //eigenvalues
    TVector3 linv1, linv2, linv3;//eigenvectors

    void calculateSphericity(int n, float *px, float *py, float *pz, float r);
    inline float p2(float px,float py,float pz);
};

inline float Sphericity::p2(float px,float py,float pz){
  return px*px+py*py+pz*pz;
}

inline TVector3 Sphericity::sphericityAxis(){
  return v1;
}

inline TVector3 Sphericity::linSphericityAxis(){
  return linv1;
}

inline float Sphericity::sphericity(){
  return 1.5*(l2+l3);
}

inline float Sphericity::linSphericity(){
  return 1.5*(linl2+linl3);
}

inline float Sphericity::aplanarity(){
  return 1.5*l3;
}

inline float Sphericity::linAplanarity(){
  return 1.5*linl3;
}

inline float Sphericity::linC(){
  return 3*(linl1*linl2+linl2*linl3+linl1*linl3);
}

inline float Sphericity::linD(){
  return 27*linl1*linl2*linl3;
}

void Sphericity::calculateSphericity(int n, float *px, float *py, float *pz, float r){
  float norm = 0;

  TMatrixD m = TMatrixD(3,3);
 
  //calculate matrix elements
  for(int i = 0; i<n; i++){
    m(0,0) += px[i]*px[i]*TMath::Power(p2(px[i],py[i],pz[i]),(r-2)/2.0);   
    m(1,1) += py[i]*py[i]*TMath::Power(p2(px[i],py[i],pz[i]),(r-2)/2.0);   
    m(2,2) += pz[i]*pz[i]*TMath::Power(p2(px[i],py[i],pz[i]),(r-2)/2.0);   
    m(1,0) += px[i]*py[i]*TMath::Power(p2(px[i],py[i],pz[i]),(r-2)/2.0);   
    m(2,0) += py[i]*pz[i]*TMath::Power(p2(px[i],py[i],pz[i]),(r-2)/2.0);   
    m(1,2) += px[i]*pz[i]*TMath::Power(p2(px[i],py[i],pz[i]),(r-2)/2.0);
    norm += TMath::Power(p2(px[i],py[i],pz[i]),r/2.0);   
  } 

  //normalize
  for(int i = 0; i<3; i++){
    for(int j = 0; j<3; j++){
      m(i,j) = m(i,j)/norm;
    }
  }

  //symmetrize other side of the matrix
  m(0,1) = m(1,0);
  m(0,2) = m(2,0);
  m(2,1) = m(1,2);
  
   
}

Sphericity::Sphericity(int n, float *px, float *py, float *pz){
  calculateSphericity(n, px, py, pz, 2);
  calculateSphericity(n, px, py, pz, 1);
}

Sphericity::~Sphericity(){}

#endif
