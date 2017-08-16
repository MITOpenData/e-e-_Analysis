/**************************************************************************************/
// Common Utilities
/**************************************************************************************/
#include <TCanvas.h>
#include <TVector3.h>
#include <TH1F.h>
#include <TH2F.h>


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

/**************************************************************************************/
// Handsome ROOT Style
/**************************************************************************************/

// Handsome Canvas for 2D histogram
TCanvas *CFViewer(const char *canvasName, const char *title, double x, double y)
{
   TCanvas *c = new TCanvas(canvasName,title,x,y);
   // The MAGIC angle
   c->SetTheta(60.839);
   c->SetPhi(38.0172);
   return c;
}

// Handsome Correlation Function Histogram
TH2F *CFTH2F(const char *name, const char *title, Int_t nBinX, Float_t xL, Float_t xH, Int_t nBinY, Float_t yL, Float_t yH, Float_t offset=1.5, Float_t offsetZ=2) {
   TH2F *h = new TH2F(name,title,nBinX,xL,xH,nBinY,yL,yH);
   h->Sumw2();
   h->SetTitleOffset(offset,"X");
   h->SetTitleOffset(offset,"Y");
   h->SetTitleOffset(offsetZ,"Z");
   h->GetXaxis()->CenterTitle();
   h->GetYaxis()->CenterTitle();
   h->GetZaxis()->CenterTitle();
   h->SetNdivisions(505,"X");
   h->SetNdivisions(505,"Y");
   h->SetNdivisions(505,"Z");
   return h;   
}
inline double ptFromThrust(TVector3 thrust, TVector3 p){
  return p.Perp(thrust); 
}
inline double thetaFromThrust(TVector3 thrust, TVector3 p){
  return p.Angle(thrust);
}


/**************************************************************************************/
// Austin Thrust Angle Rotators
/**************************************************************************************/

//this is actually rapidity w/ pion mass assumption
inline double etaFromThrust(TVector3 thrust, TVector3 p){
  float pl = p*thrust.Unit();//logitudinal momentum component
  float E = TMath::Power(p.Mag2()+0.13957*0.13957,0.5);//energy w/ mass assumption
  return 0.5*TMath::Log((E+pl)/(E-pl));//rapidity
}
inline double phiFromThrust(TVector3 thrust, TVector3 p){
  TVector3 pt = p-((p*thrust.Unit())*thrust.Unit());//pt vector
  TVector3 z = TVector3(0,0,1);
  TVector3 phiOrigin = thrust.Unit().Cross((thrust.Unit().Cross(z)));//vector that will be phi=0 (in plane of thrust and beam line
  double phi = pt.Angle(phiOrigin);//get phi from 0 to pi
  
  //determine sign of phi based on cross product of pt and origin
  if( (phiOrigin.Cross(pt.Unit()))*thrust >= 0) return phi;
  else return -phi;
}


// BELLE Particle Definition 
enum SIMPLEPID {BELLE_PHOTON, BELLE_ELECTRON, BELLE_PION, BELLE_MUON, BELLE_KAON, BELLE_PROTON};

// ALEPH Particle Flow Classification
enum SIMPLEPWFLAG {ALEPH_CHARGED_TRACK, ALEPH_CHARGED_LEPTONS1, ALEPH_CHARGED_LEPTONS2, ALEPH_V0, ALEPH_PHOTON, ALEPH_NEUTRAL_HADRON};


