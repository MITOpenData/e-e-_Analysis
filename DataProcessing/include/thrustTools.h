#ifndef THRUSTTOOLS
#define THRUSTTOOLS

#include "TVector3.h"
#include "TMath.h"
#include <vector>

double dphi(double phi1,double phi2)
{
    double a=phi1-phi2;
    while (a<-TMath::Pi()) a+=2*TMath::Pi();
    while (a>TMath::Pi()) a-=2*TMath::Pi();
    
    if (a<-TMath::Pi()/2) a=2*TMath::Pi()+a;
    return a;
}

inline double ptFromThrust(TVector3 thrust, TVector3 p){
  return p.Perp(thrust); 
}

inline double thetaFromThrust(TVector3 thrust, TVector3 p){
  return p.Angle(thrust);
}

//this is actually rapidity w/ pion mass assumption
inline double etaFromThrust(TVector3 thrust, TVector3 p){
  float pl = p*thrust.Unit();//logitudinal momentum component
  float E = TMath::Power(p.Mag2()+0.13957*0.13957,0.5);//energy w/ mass assumption
  return 0.5*TMath::Log((E+pl)/(E-pl));//rapidity
}

inline double phiFromThrust(TVector3 thrust, TVector3 p){
  TVector3 pt = p-((p*thrust.Unit())*(thrust.Unit()));//pt vector
  TVector3 z = TVector3(0,0,1);
  TVector3 phiOrigin = thrust.Unit().Cross((thrust.Unit().Cross(z)));//vector that will be phi=0 (in plane of thrust and beam line
  double phi = pt.Angle(phiOrigin);//get phi from 0 to pi

  //determine sign of phi based on cross product of pt and origin
  if( (phiOrigin.Cross(pt.Unit()))*thrust >= 0) return phi;
  else return -phi;
}

//based on code from herwig: http://herwig.hepforge.org/svn/tags/herwig-2-0-beta/Analysis/EventShapes.cc
//ported by A. Baty
//n is number of particles, px,py,pz are arrays of momentum components
TVector3 getThrust(int n, float *px, float *py, float *pz){
  TVector3 thrust = TVector3(0,0,0);

  if(n==1){//thrust is just the particle
    return TVector3(px[0],py[0],pz[0]);   
  }

  if(n==2){//special case for 2 particles
    if(TMath::Power(px[0],2)+TMath::Power(py[0],2)+TMath::Power(pz[0],2) >= TMath::Power(px[1],2)+TMath::Power(py[1],2)+TMath::Power(pz[1],2)){
      return TVector3(px[0],py[0],pz[0]);
    }
    else{
      return TVector3(px[1],py[1],pz[1]);
    }
  }

  if(n==3){//combine lowest 2 magnitude momentum, then use same algo as n=2
    if(TMath::Power(px[0],2)+TMath::Power(py[0],2)+TMath::Power(pz[0],2) >= TMath::Power(px[1],2)+TMath::Power(py[1],2)+TMath::Power(pz[1],2)){
      if(TMath::Power(px[0],2)+TMath::Power(py[0],2)+TMath::Power(pz[0],2) >= TMath::Power(px[2],2)+TMath::Power(py[2],2)+TMath::Power(pz[2],2)){
        return TVector3(px[0],py[0],pz[0]);//0 is largest momentum
      }
      else{
        return TVector3(px[2],py[2],pz[2]);//2 is largest momentum
      } 
    }
    else{
      if(TMath::Power(px[1],2)+TMath::Power(py[1],2)+TMath::Power(pz[1],2) >= TMath::Power(px[2],2)+TMath::Power(py[2],2)+TMath::Power(pz[2],2)){
        return TVector3(px[1],py[1],pz[1]);//1 is largest momentum
      }
      else{
        return TVector3(px[2],py[2],pz[2]);//2 is largest momentum
      } 
    }
  }


  if(n>3){
    //make vector of TVector3's of each particle
    std::vector< TVector3 > pVec;
    for(int i = 0; i<n; i++){
      TVector3 v = TVector3(px[i],py[i],pz[i]);
      pVec.push_back(v); 
    }
    //std::cout << pVecs.at(0).x();
  
    TVector3 cross; 
    float t = 0;
    for(int i = 1; i<n; i++){//loop through all possible cross products of 2 unique vectors
      for(int j = 0; j<i; j++){
        cross = pVec.at(i).Cross(pVec.at(j));
        TVector3 ptot = TVector3(0,0,0);
 
        for(int k = 0; k<n; k++){ //loop through all 3rd particles not used for the cross product
          if(k!=i && k!=j){
            if(pVec.at(k)*cross > 0){//if dot product is >0
              ptot += pVec.at(k); 
            }
            else{
              ptot -= pVec.at(k);
            }
          }
        }

        std::vector< TVector3 > cpm;//add or subtract in last 2 vectors used for cross product
        cpm.push_back(ptot - pVec.at(j) - pVec.at(i));
        cpm.push_back(ptot - pVec.at(j) + pVec.at(i));
        cpm.push_back(ptot + pVec.at(j) - pVec.at(i));
        cpm.push_back(ptot + pVec.at(j) + pVec.at(i));
        for(vector<TVector3>::iterator it = cpm.begin(); it != cpm.end(); it++){
          float tval = (*it).Mag2();
          if(tval > t){
            t = tval;
            thrust = *it;
          }
        }
      }
    } 
  }
  return thrust;
}




#endif
