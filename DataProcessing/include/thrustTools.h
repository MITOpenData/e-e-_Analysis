#ifndef THRUSTTOOLS
#define THRUSTTOOLS

#include "TVector3.h"
#include "TMath.h"
#include <vector>
#include "particleData.h"
#include "eventData.h"

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
inline double rapFromThrust(TVector3 thrust, TVector3 p, float mass){
  float pl = p*(thrust.Unit());//logitudinal momentum component
  float E = TMath::Power(p.Mag2()+mass*mass,0.5);//energy
  return 0.5*TMath::Log((E+pl)/(E-pl));//rapidity
}

inline double etaFromThrust(TVector3 thrust, TVector3 p){
  return -TMath::Log( TMath::Tan( thetaFromThrust(thrust,p)/2.0));
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

void setThrustVariables(particleData *p, eventData *e, TVector3 thrust, TVector3 chThrust){
  int nTrk = 0;
  for(int i = 0; i< p->nParticle; i++){
    TVector3 part = TVector3(p->px[i], p->py[i], p->pz[i]);
    p->pt_wrtThr[i] = ptFromThrust(thrust, part);
    p->eta_wrtThr[i] = etaFromThrust(thrust, part);
    p->theta_wrtThr[i] = thetaFromThrust(thrust, part);
    p->phi_wrtThr[i] = phiFromThrust(thrust, part);
    
    p->pt_wrtChThr[i] = ptFromThrust(chThrust, part);
    p->eta_wrtChThr[i] = etaFromThrust(chThrust, part);
    p->theta_wrtChThr[i] = thetaFromThrust(chThrust, part);
    p->phi_wrtChThr[i] = phiFromThrust(chThrust, part);
    
    if(p->pwflag[i]==0 && p->pt_wrtThr[i]>0.4) nTrk++;
  }
  e->nChargedHadrons_GT0p4Thrust = nTrk;
}


/* Almost a direct copy and paste from Belle standard Thrust axis algorithm */
template <class Iterator, class Function>
TVector3 thrustBelleSTD(Iterator begin, Iterator end, Function func)
{
  // Temporary variables
  Iterator p, q;
  TVector3 rvec, Axis;

  double sump = 0;
  for (p = begin; p != end; p++)
    sump += ((TVector3)func(*p)).Mag();

  // Thrust and thrust vectors

  double Thru = 0;
  for (p = begin; p != end; p++) {
    TVector3 rvec(func(*p));
    if (rvec.z() <= 0.0) rvec = -rvec;

    double s = rvec.Mag();
    if (s != 0.0) rvec *= (1/s);

    for (Iterator loopcount = begin; loopcount != end; loopcount++) {
      TVector3 rprev(rvec);
      rvec = TVector3(); // clear

      for (q = begin; q != end; q++) {
        const TVector3 qvec(func(*q));
        rvec += (qvec.Dot(rprev) >= 0) ? qvec : - qvec;
      }

      for (q = begin; q != end; q++) {
        const TVector3 qvec(func(*q));
        if (qvec.Dot(rvec) * qvec.Dot(rprev) < 0) break;
      }

      if (q == end) break;
    }

    double ttmp = 0.0;
    for (q = begin; q != end; q++) {
      const TVector3 qvec = func(*q);
      ttmp += std::fabs(qvec.Dot(rvec));
    }
    ttmp /= (sump * rvec.Mag());
    rvec *= 1/rvec.Mag();
    if (ttmp > Thru) {
      Thru = ttmp;
      Axis = rvec;
    }
  }
  Axis *= Thru;
  return Axis;
}

// ----------------------------------------------------------------------
// SelfFunc - retrieve the pointer to a function which returns itself
// example:
//   list<Vector3> vl;
//   Vector3 t = thrust(vl.begin(), vp.end(), SelfFunc(Vector3()));
//
//   list<Vector3 *> vp;
//   Vector3 t = thrust(vp.begin(), vp.end(), SelfFunc(Vector3()));
// ----------------------------------------------------------------------

template <class T>
class ptr_to_self_func {
 protected:
	typedef T (*pfun)(T &);
 public:
 ptr_to_self_func() : ptr(NULL) {};
	T operator()(T& t) const { return t; };
	T operator()(T *t) const { return *t; };
	const T operator()(const T& t) const { return t; };
	const T operator()(const T *t) const { return *t; };
 protected:
	const pfun ptr;
};

template <class T>
ptr_to_self_func<T> SelfFunc(const T &) {
	return ptr_to_self_func<T>();
}

/* wrapper of the Belle thrust axis algorithm */
TVector3 getThrustBelle(int n, float *px, float *py, float *pz){
  std::vector<TVector3> momenta;
  for (int i = 0; i < n; ++i) {
    momenta.push_back(TVector3(px[i], py[i], pz[i]));
  }
  return thrustBelleSTD(momenta.begin(), momenta.end(), SelfFunc(TVector3()));
}

TVector3 getChargedThrustBelle(int n, float *px, float *py, float *pz, int *pwflag){
  std::vector<TVector3> momenta;
  for (int i = 0; i < n; ++i) {
    if(pwflag[i]!=0) continue;
    momenta.push_back(TVector3(px[i], py[i], pz[i]));
  }
  return thrustBelleSTD(momenta.begin(), momenta.end(), SelfFunc(TVector3()));
}

//based on code from herwig: http://herwig.hepforge.org/svn/tags/herwig-2-0-beta/Analysis/EventShapes.cc
//ported by A. Baty
//n is number of particles, px,py,pz are arrays of momentum components
TVector3 getThrustHerwig(int n, float *px, float *py, float *pz){
  TVector3 thrust = TVector3(0,0,0);
  float pSum = 0;
  for(int t = 0; t<n; t++) pSum += TVector3(px[t],py[t],pz[t]).Mag();
 
  if(n<=0) return thrust;

  if(n==1){//thrust is just the particle
    TVector3 thrust = TVector3(px[0],py[0],pz[0]);   
    thrust.SetMag(thrust.Mag()/pSum);
    return thrust;
  }

  if(n==2){//special case for 2 particles
    if(TMath::Power(px[0],2)+TMath::Power(py[0],2)+TMath::Power(pz[0],2) >= TMath::Power(px[1],2)+TMath::Power(py[1],2)+TMath::Power(pz[1],2)){
      TVector3 thrust = TVector3(px[0],py[0],pz[0]);   
      thrust.SetMag(thrust.Mag()/pSum);
      return thrust;
    }
    else{
      TVector3 thrust = TVector3(px[1],py[1],pz[1]);   
      thrust.SetMag(thrust.Mag()/pSum);
      return thrust;
    }
  }

  if(n==3){//combine lowest 2 magnitude momentum, then use same algo as n=2
    if(TMath::Power(px[0],2)+TMath::Power(py[0],2)+TMath::Power(pz[0],2) >= TMath::Power(px[1],2)+TMath::Power(py[1],2)+TMath::Power(pz[1],2)){
      if(TMath::Power(px[0],2)+TMath::Power(py[0],2)+TMath::Power(pz[0],2) >= TMath::Power(px[2],2)+TMath::Power(py[2],2)+TMath::Power(pz[2],2)){
        TVector3 thrust = TVector3(px[0],py[0],pz[0]);   
        thrust.SetMag(thrust.Mag()/pSum);
        return thrust;
      }
      else{
        TVector3 thrust = TVector3(px[2],py[2],pz[2]);   
        thrust.SetMag(thrust.Mag()/pSum);
        return thrust;
      } 
    }
    else{
      if(TMath::Power(px[1],2)+TMath::Power(py[1],2)+TMath::Power(pz[1],2) >= TMath::Power(px[2],2)+TMath::Power(py[2],2)+TMath::Power(pz[2],2)){
        TVector3 thrust = TVector3(px[1],py[1],pz[1]);   
        thrust.SetMag(thrust.Mag()/pSum);
        return thrust;
      }
      else{
        TVector3 thrust = TVector3(px[2],py[2],pz[2]);   
        thrust.SetMag(thrust.Mag()/pSum);
        return thrust;
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
        for(std::vector<TVector3>::iterator it = cpm.begin(); it != cpm.end(); it++){
          float tval = (*it).Mag2();
          if(tval > t){
            t = tval;
            thrust = *it;
          }
        }
      }
    } 
  }
  thrust.SetMag(thrust.Mag()/pSum);
  return thrust;
}

//almost a straight copy of above, but filter for pwflag==0 (tracks)
TVector3 getChargedThrustHerwig(int n, float *px, float *py, float *pz, int *pwflag){
  float nTrk = 0;
  float pSum = 0;
  for(int t = 0; t<n; t++){
    if(pwflag[t]!=0) continue;
    pSum += TVector3(px[t],py[t],pz[t]).Mag();
    nTrk++;
  }
  
  TVector3 thrust = TVector3(0,0,0);
  if(nTrk<=0) return thrust;

  if(nTrk==1){//thrust is just the particle
    for(int t = 0; t<n; t++){
      if(pwflag[t]==0){
        TVector3 thrust = TVector3(px[t],py[t],pz[t]);   
        thrust.SetMag(thrust.Mag()/pSum);
        return thrust;
      }
    }
  }
 
  if(nTrk==2){//special case for 2 particles
    int n1 = -1, n2 = -1;
    for(int t = 0; t<n; t++){
      if(pwflag[t]==0 && n1==-1){ n1 = t; continue;}
      if(pwflag[t]==0 && n1!=-1) n2 = t;  
    }

    if(TMath::Power(px[n1],2)+TMath::Power(py[n1],2)+TMath::Power(pz[n1],2) >= TMath::Power(px[n2],2)+TMath::Power(py[n2],2)+TMath::Power(pz[n2],2)){
      TVector3 thrust = TVector3(px[n1],py[n1],pz[n1]);   
      thrust.SetMag(thrust.Mag()/pSum);
      return thrust;
    }
    else{
      TVector3 thrust = TVector3(px[n2],py[n2],pz[n2]);   
      thrust.SetMag(thrust.Mag()/pSum);
      return thrust;
    }
  }
 
  if(nTrk==3){//combine lowest 2 magnitude momentum, then use same algo as n=2
    int n1 = -1, n2 = -1, n3 = -1;
    for(int t = 0; t<n; t++){
      if(pwflag[t]==0 && n1==-1){ n1 = t; continue;}
      if(pwflag[t]==0 && n1!=-1 && n2==-1){ n2 = t;  continue;}
      if(pwflag[t]==0 && n1!=-1 && n2!=-1) n3 = t;  
    }
    if(TMath::Power(px[n1],2)+TMath::Power(py[n1],2)+TMath::Power(pz[n1],2) >= TMath::Power(px[n2],2)+TMath::Power(py[n2],2)+TMath::Power(pz[n2],2)){
      if(TMath::Power(px[n1],2)+TMath::Power(py[n1],2)+TMath::Power(pz[n1],2) >= TMath::Power(px[n3],2)+TMath::Power(py[n3],2)+TMath::Power(pz[n3],2)){
        TVector3 thrust = TVector3(px[n1],py[n1],pz[n1]);//n1 is largest p   
        thrust.SetMag(thrust.Mag()/pSum);
        return thrust;
      }
      else{
        TVector3 thrust = TVector3(px[n3],py[n3],pz[n3]);//n3 is largest p   
        thrust.SetMag(thrust.Mag()/pSum);
        return thrust;
      } 
    }
    else{
      if(TMath::Power(px[n2],2)+TMath::Power(py[n2],2)+TMath::Power(pz[n2],2) >= TMath::Power(px[n3],2)+TMath::Power(py[n3],2)+TMath::Power(pz[n3],2)){
        TVector3 thrust = TVector3(px[n2],py[n2],pz[n2]);//n3 is largest p   
        thrust.SetMag(thrust.Mag()/pSum);
        return thrust;
      }
      else{
        TVector3 thrust = TVector3(px[n3],py[n3],pz[n3]);//n3 is largest p   
        thrust.SetMag(thrust.Mag()/pSum);
        return thrust;
      } 
    }
  }
 
  if(nTrk>3){
    //make vector of TVector3's of each particle
    std::vector< TVector3 > pVec;
    for(int i = 0; i<n; i++){
      if(pwflag[i]!=0) continue;
      TVector3 v = TVector3(px[i],py[i],pz[i]);
      pVec.push_back(v); 
    }
    //std::cout << pVecs.at(0).x();
  
    TVector3 cross; 
    float t = 0;
    for(unsigned int i = 1; i<pVec.size(); i++){//loop through all possible cross products of 2 unique vectors
      for(unsigned int j = 0; j<i; j++){
        cross = pVec.at(i).Cross(pVec.at(j));
        TVector3 ptot = TVector3(0,0,0);
 
        for(unsigned int k = 0; k<pVec.size(); k++){ //loop through all 3rd particles not used for the cross product
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
        for(std::vector<TVector3>::iterator it = cpm.begin(); it != cpm.end(); it++){
          float tval = (*it).Mag2();
          if(tval > t){
            t = tval;
            thrust = *it;
          }
        }
      }
    } 
  }
  thrust.SetMag(thrust.Mag()/pSum);
  return thrust;
}


struct THRUST
{
  enum algorithm {HERWIG, BELLE, OPTIMAL};
};

// Interface to the different algorithm to compute thrust axis
TVector3 getThrust(int n, float *px, float *py, float *pz, THRUST::algorithm algo=THRUST::HERWIG){
  TVector3 thrustAxis(0, 0, 0);
  switch (algo) {
  case THRUST::HERWIG: thrustAxis = getThrustHerwig(n, px, py, pz); break;
  case THRUST::BELLE: thrustAxis = getThrustBelle(n, px, py, pz); break;
  case THRUST::OPTIMAL: {
    if (n < 4) {
      thrustAxis = getThrustBelle(n, px, py, pz);
    } else {
      thrustAxis = getThrustHerwig(n, px, py, pz);
    }
    break;
  }
  default:
    break;
  }
  return thrustAxis;
}

TVector3 getChargedThrust(int n, float *px, float *py, float *pz, int *pwflag, THRUST::algorithm algo=THRUST::HERWIG){
  float nTrk = 0;
  for(int t = 0; t<n; t++){
    if(pwflag[t]!=0) continue;
    nTrk++;
  }

  TVector3 thrustAxis(0, 0, 0);
  switch (algo) {
  case THRUST::HERWIG: thrustAxis = getChargedThrustHerwig(n, px, py, pz, pwflag); break;
  case THRUST::BELLE: thrustAxis = getChargedThrustBelle(n, px, py, pz, pwflag); break;
  case THRUST::OPTIMAL: {
    if (nTrk < 4) {
      thrustAxis = getChargedThrustBelle(n, px, py, pz,pwflag);
    } else {
      thrustAxis = getChargedThrustHerwig(n, px, py, pz, pwflag);
    }
    break;
  }
  default:
    break;
  }
  return thrustAxis;
}

#endif
