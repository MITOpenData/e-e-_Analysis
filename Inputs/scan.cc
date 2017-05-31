// The numbers of particles in each event are stored in a histogram.
// All particle masses in the analyzed events are calculated and stored in a histogram.
// All pion masses in the analyzed events are calculated and stored in a histogram.
// The neutral pion (pizero) and the charged pion are differentiated and respective
// mass distributions are produced.
//+
// File : pion.cc
// Description : Analyze a data file
// 
// Author : Ryosuke Itoh, IPNS, KEK
// Date : 1 - Mar - 2004
//-
// Modified by T.Nozaki  (2005/09/21)
// Modified by S.Nishida (2007/02/03)
// Modified by T.Nozaki (2009/07/01) Comments in English
#include <cstdio>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TClonesArray.h"
#include "TBrowser.h"
#include "TLorentzVector.h"
#include <iostream>

using namespace std;


#define PI 3.1415926
class particleData
{
    public:
    int nParticle;
    Float_t pt[100000];
    Float_t pid[100000];
    Float_t eta[100000];
    Float_t theta[100000];
    Float_t phi[100000];
    Float_t etathrust[100000];
    Float_t thetathrust[100000];
    Float_t phithrust[100000];

    Float_t mass[10000];
};
double dphi(double phi1,double phi2)
{
    double a=phi1-phi2;
    while (a<-PI) a+=2*PI;
    while (a>PI) a-=2*PI;
    
    if (a<-PI/2) a=2*PI+a;
    return a;
}
void analysis (  ) {

void calculateThrust(float&,float&, int,float[],float[],float[],float[]);

  TFile *f = new TFile("../Inputs/output_2_withtheta.root");
  TTree *tentry = (TTree*)f->Get("t");
  
  Int_t nParticle_input;
  Float_t pt_input[100000];
  Float_t eta_input[100000];
  Float_t theta_input[100000];
  Float_t pid_input[100000];
  Float_t phi_input[100000];
  Float_t mass_input[100000];
  
  tentry->SetBranchAddress("nParticle",&nParticle_input);
  tentry->SetBranchAddress("pt",pt_input);
  tentry->SetBranchAddress("eta",eta_input);
  tentry->SetBranchAddress("theta",theta_input);
  tentry->SetBranchAddress("pid",pid_input);
  tentry->SetBranchAddress("phi",phi_input);
  tentry->SetBranchAddress("mass",mass_input);
   
  TFile *hf = new TFile("../Inputs/output_2_withtheta_skimmed.root", "RECREATE" );
  TTree *tout = new TTree("t_thrust","");
  
  particleData pData;
  tout->Branch("nParticle", &pData.nParticle,"nParticle/I");
  tout->Branch("pt", pData.pt,"pt[nParticle]/F");
  tout->Branch("pid", pData.pid,"pid[nParticle]/F");
  tout->Branch("eta", pData.eta,"eta[nParticle]/F");
  tout->Branch("theta", pData.theta,"theta[nParticle]/F");
  tout->Branch("phi", pData.phi,"phi[nParticle]/F");
  tout->Branch("etathrust", pData.etathrust,"etathrust[nParticle]/F");
  tout->Branch("thetathrust", pData.thetathrust,"thetathrust[nParticle]/F");
  tout->Branch("phithrust", pData.phithrust,"phithrust[nParticle]/F");
  tout->Branch("mass", pData.mass,"mass[nParticle]/F");
  
  int nevent = (int)tentry->GetEntries();
  for ( int i=0;i<10000;i++ ) {
    if ( i%100 == 0 ) printf ( "***** Event %d\n", i );
    int nb = tentry->GetEntry ( i );
    pData.nParticle=nParticle_input;
    
    float theta_trust,phi_trust;
    //void calculateThrust(float& thetaThrust , float& phiThrust, int n, float mypt[], float mytheta[], float myphi[], float eta[]){
    calculateThrust(theta_trust,phi_trust,nParticle_input,pt_input, theta_input, phi_input, eta_input);
    
    pData.nParticle=nParticle_input;

    
    for ( int j=0;j<nParticle_input;j++ ) {
    
      pData.pt[j]=pt_input[j];
      pData.theta[j]=theta_input[j];
      pData.eta[j]=eta_input[j];
      pData.phi[j]=phi_input[j];
      pData.mass[j]=mass_input[j];
      
      pData.etathrust[j]=0.;
      pData.thetathrust[j]=theta_trust;
      pData.phithrust[j]=phi_trust;
      
    }
    tout->Fill();
    // Clear at the end of particle loop.
  }
  
  hf->Write();
}   


void calculateThrust(float& thetaThrust , float& phiThrust, int _mynparticles, float _mypt[], float _mytheta[], float _myphi[], float _eta[]){

  int nstepstheta=100;
  double mintheta=0.;
  double maxtheta=3.14/2;  
  
  int nstepsphi=100;
  double minphi=0.;
  double maxphi=3.14;
  
  double sizetheta=(maxtheta-mintheta)/(double)nstepstheta;
  double sizephi=(maxphi-minphi)/(double)nstepsphi;

  double mytheta,myphi;
  double x_versor,y_versor,z_versor;
  TVector3 unity_versor;
  TVector3 particle_vector;
  
  double maxscalar,theta_max,phi_max;
  maxscalar=-1;
  
  for ( int itheta=0; itheta<nstepstheta; itheta++ ) {
    mytheta=sizetheta*(double)itheta;
    for ( int iphi=0; iphi<nstepsphi; iphi++ ) {
      myphi=sizephi*(double)iphi;

      x_versor=sin(mytheta)*cos(myphi);
      y_versor=sin(mytheta)*sin(myphi);
      z_versor=cos(mytheta);
        
      unity_versor.SetXYZ(x_versor,y_versor,z_versor);
      double scalar_prod;
      scalar_prod=0.;
       
      for ( int j=0;j<_mynparticles;j++ ) {
        float pt1 = _mypt[j];
        float theta1 = _mytheta[j];
        float phi1 = _myphi[j];
        float px1=pt1*sin(theta1)*cos(phi1);
        float py1=pt1*sin(theta1)*sin(phi1);
        float pz1=pt1*sin(phi1);
        particle_vector.SetXYZ(px1,py1,pz1);
        scalar_prod=scalar_prod+fabs(particle_vector.Dot(unity_versor));         
      }//end of loop over particles

      if (scalar_prod>maxscalar) { maxscalar=scalar_prod; theta_max=mytheta; phi_max=myphi;} 
    }//loop over theta
  }//loop over phi
  
  thetaThrust=theta_max;
  phiThrust=phi_max;
  cout<<"thetaThrust="<<thetaThrust<<endl;
  cout<<"phiThrust="<<phiThrust<<endl;
}

