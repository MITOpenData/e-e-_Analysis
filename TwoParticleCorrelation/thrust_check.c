#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <iostream>
#include <stdlib.h>
#include <TMath.h>
#include <TLorentzVector.h>

using namespace std;


#define PI 3.1415926
enum SIMPLEPWFLAG {CHARGED_TRACK, CHARGED_LEPTONS1, CHARGED_LEPTONS2, V0, PHOTON, NEUTRAL_HADRON};


double dphi(double phi1,double phi2)
{
    double a=phi1-phi2;
    while (a<-PI) a+=2*PI;
    while (a>PI) a-=2*PI;
    
    if (a<-PI/2) a=2*PI+a;
    return a;
}


void analysis(){
  TString filename;
  //if(isBelle) filename="/data/flowex/Datasamples/Belle/output_2_withtheta.root";
  //else filename="/data/flowex/Datasamples/LEP2_MAIN/ROOTfiles/cleaned_ALEPH_Data-all.aleph.root";
  
  filename = "/mnt/c/Users/Bibek Kumar Pandit/Desktop/Root_Directory/StudyMult/DataFiles/ROOTfiles/cleaned_ALEPH_Data2-all.aleph.root";
  
  TFile *f = new TFile(filename.Data());
  TTree *t1 = (TTree*)f->Get("t");
  Int_t nParticle;
  Float_t pt[50000];
  Float_t eta[50000];
  Float_t pid[50000];
  Float_t phi[50000];
  Float_t mass[50000];
  Float_t pwflag[50000];
  Float_t px[100000];
  Float_t py[100000];
  Float_t pz[100000];
  Float_t TTheta;
  Float_t TPhi;
  
  t1->SetBranchAddress("nParticle",&nParticle);
  t1->SetBranchAddress("pt",pt);
  t1->SetBranchAddress("eta",eta);
  t1->SetBranchAddress("pid",pid);
  t1->SetBranchAddress("phi",phi);
  t1->SetBranchAddress("mass",mass);
  t1->SetBranchAddress("pwflag",pwflag);
  t1->SetBranchAddress("px",px);
  t1->SetBranchAddress("py",py);
  t1->SetBranchAddress("pz",pz);
  t1->SetBranchAddress("TTheta", &TTheta);
  t1->SetBranchAddress("TPhi", &TPhi);

  TFile *f_mix = new TFile(filename.Data());
  TTree *t1_mix = (TTree*)f_mix->Get("t");
  Int_t nParticle_mix;
  Float_t pt_mix[50000];
  Float_t eta_mix[50000];
  Float_t pid_mix[50000];
  Float_t phi_mix[50000];
  Float_t mass_mix[50000];
  Float_t pwflag_mix[50000];
  Float_t px_mix[100000];
  Float_t py_mix[100000];
  Float_t pz_mix[100000];
  Float_t TTheta_mix[100000];
  Float_t TPhi_mix[100000];
  
  t1_mix->SetBranchAddress("nParticle",&nParticle_mix);
  t1_mix->SetBranchAddress("pt",pt_mix);
  t1_mix->SetBranchAddress("eta",eta_mix);
  t1_mix->SetBranchAddress("pid",pid_mix);
  t1_mix->SetBranchAddress("phi",phi_mix);
  t1_mix->SetBranchAddress("mass",mass_mix);
  t1_mix->SetBranchAddress("pwflag",pwflag_mix);
  t1_mix->SetBranchAddress("px",px_mix);
  t1_mix->SetBranchAddress("py",py_mix);
  t1_mix->SetBranchAddress("pz",pz_mix);
  t1_mix->SetBranchAddress("TTheta", TTheta_mix);
  t1_mix->SetBranchAddress("TPhi", TPhi_mix);
  
  Int_t nevent = (Int_t)t1->GetEntries();
  
  for (Int_t i=0;i<5;i++) {

  //if (i%10000==0) cout <<i<<"/"<<nevent_process<<endl;
  t1->GetEntry(i);
  cout<<"-------------------------------------------------"<<endl;
  float theta1 = TTheta;
  float phi1 = TPhi;
  
  TVector3 v(1,1,1);
  
  TVector3 v1;
  v.SetPhi(phi1);
  v.SetTheta(theta1);
  v1.SetX(-v.Y());
  v1.SetY(v.X());
  v.Rotate(-theta1, v1);
 
  cout<<v.X()<<endl;
  cout<<v.Y()<<endl;
  cout<<v.Z()<<endl;
 
  cout<<endl;
  }
}