/**************************************************************************************/
// Common Utilities
/**************************************************************************************/
#include <TCanvas.h>
#include <TVector3.h>

#define PI 3.14159265358979


// Calculate Delta Phi
double dphi(double phi1,double phi2)
{
    double a=phi1-phi2;
    while (a<-PI) a+=2*PI;
    while (a>PI) a-=2*PI;
    
    if (a<-PI/2) a=2*PI+a;
    return a;
}

TCanvas *CFViewer(char *canvasName, char *title, double x, double y)
{
   TCanvas *c = new TCanvas(canvasName,title,x,y);
   // The MAGIC angle
   c->SetTheta(60.839);
   c->SetPhi(38.0172);
   return c;
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
  double phi = (thrust.Cross(p)).Angle(thrust.Orthogonal());
  if( ((thrust.Orthogonal().Unit()).Cross(((thrust.Cross(p).Unit()))))*thrust >= 0) return phi;
  else return -phi;
}


// BELLE Particle Definition 
enum SIMPLEPID {BELLE_PHOTON, BELLE_ELECTRON, BELLE_PION, BELLE_MUON, BELLE_KAON, BELLE_PROTON};

// ALEPH Particle Flow Classification
enum SIMPLEPWFLAG {ALEPH_CHARGED_TRACK, ALEPH_CHARGED_LEPTONS1, ALEPH_CHARGED_LEPTONS2, ALEPH_V0, ALEPH_PHOTON, ALEPH_NEUTRAL_HADRON};


