#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <iostream>
#include <stdlib.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>
using namespace std;


#define PI 3.1415926
enum SIMPLEPWFLAG {CHARGED_TRACK, CHARGED_LEPTONS1, CHARGED_LEPTONS2, V0, PHOTON, NEUTRAL_HADRON};



void analysis(int maxevt=10){

  TString filename="/mnt/c/Users/Bibek Kumar Pandit/Desktop/Root_Directory/StudyMult/LEP2/ROOTfiles/cleaned_ALEPH_Data-all.aleph.root";
  
  TFile *f = new TFile(filename.Data());
  TTree *t1 = (TTree*)f->Get("t");
  Int_t nParticle;
  Float_t pt[50000];
  Float_t px[50000];
  Float_t py[50000];
  Float_t pz[50000];
  Float_t eta[50000];
  Float_t theta[50000];
  Float_t pid[50000];
  Float_t phi[50000];
  Float_t mass[50000];
  Float_t pwflag[50000];
  
  t1->SetBranchAddress("nParticle",&nParticle);
  t1->SetBranchAddress("pt",pt);
  t1->SetBranchAddress("px",px);
  t1->SetBranchAddress("py",py);
  t1->SetBranchAddress("pz",pz);
  t1->SetBranchAddress("eta",eta);
  t1->SetBranchAddress("theta",theta);
  t1->SetBranchAddress("pid",pid);
  t1->SetBranchAddress("phi",phi);
  t1->SetBranchAddress("mass",mass);
  t1->SetBranchAddress("pwflag",pwflag);

  /*
  x=sin(theta)*cos(phi)
  y=sin(theta)*sin(phi)
  z=cos(theta)
  */
  
  int nstepstheta=100;
  double mintheta=0.;
  double maxtheta=3.14;  
  
  int nstepsphi=100;
  double minphi=-3.14;
  double maxphi=3.14;
  
  double sizetheta=(maxtheta-mintheta)/(double)nstepstheta;
  double sizephi=(maxphi-minphi)/(double)nstepsphi;

  double mytheta,myphi;
  double x_versor,y_versor,z_versor;
  TVector3 unity_versor;
  TVector3 particle_vector;
  
  double maxscalar,theta_max,phi_max;
  
    // all entries and fill the histograms
  Int_t nevent = (Int_t)t1->GetEntries();
  int nevent_process = nevent;
  
  if (maxevt > 0 && maxevt < nevent) nevent_process = maxevt;
  //const int nevent_process = 10;
  
  TH2F*hprojection[nevent_process];
  TH2F*hparticles[nevent_process];
  
  TCanvas *c1 = new TCanvas("particles_on_thrust", "particles_on_thrust", 600,  600);
  /*
  px[0] = 1.5;
  py[0] = 2.5;
  pz[0] = 0.1;
  
  px[1] = 1.6;
  py[1] = 2.6;
  pz[1] = 0.0;
  
  px[2] = -1.6;
  py[2] = -2.5;
  pz[2] = -0.1;
  
  px[3] = -1.5;
  py[3] = -2.6;
  pz[3] = 0.0;
  */
  for (Int_t i=0;i<nevent_process;i++) {  
  if (i%10==0) cout <<i<<"/"<<nevent_process<<endl;
  t1->GetEntry(i);
  
  int nparticles = nParticle;
  
  hprojection[i]=new TH2F(Form("hprojection%d",i),Form("hprojection%d",i),nstepstheta,mintheta,maxtheta,nstepsphi,minphi,maxphi);
  hprojection[i]->GetXaxis()->SetTitle("#theta");
  hprojection[i]->GetYaxis()->SetTitle("#phi");
  
  hparticles[i] = new TH2F (Form("hprojection%d",i),Form("hprojection%d",i), nstepstheta,mintheta,maxtheta,nstepsphi,minphi,maxphi);

  
  
  maxscalar=-1;

  for ( int itheta=0; itheta<nstepstheta; itheta++ ) {
    mytheta=sizetheta*(double)itheta;
    myphi = -3.14;
    for ( int iphi=0; iphi<nstepsphi; iphi++ ) {
      //myphi+=sizephi*(double)iphi;

      x_versor=sin(mytheta)*cos(myphi);
      y_versor=sin(mytheta)*sin(myphi);
      z_versor=cos(mytheta);
      
      unity_versor.SetXYZ(x_versor,y_versor,z_versor);
      double scalar_prod;
      
      scalar_prod=0.;
     
      for ( int j=0;j<nparticles;j++ ) {
       
       float px1=px[j];
       float py1=py[j];
       float pz1=pz[j];
       /*
       float px1=pt1*sin(theta1)*cos(phi1);
       float py1=pt1*sin(theta1)*sin(phi1);
       float pz1=pt1*cos(theta1);
       */
       particle_vector.SetXYZ(px1,py1,pz1);
       scalar_prod = scalar_prod + fabs(particle_vector.Dot(unity_versor));         
      }//end of loop over particles

     hprojection[i]->SetBinContent(itheta,iphi,scalar_prod);  
     
     myphi+=sizephi;
     if (scalar_prod>maxscalar) { maxscalar=scalar_prod; theta_max=mytheta; phi_max=myphi;} 
    }//loop over theta
 }//loop over phi
 float x_momentum = 0.;
 float y_momentum = 0.;
 float z_momentum = 0.;
 for (int j = 0; j< nparticles; j++)
 {
  x_momentum += px[j];
  y_momentum += py[j];
  z_momentum += pz[j];
  hparticles[i]->SetMarkerStyle(20);
  /*
  TVector3 v;
  v.SetXYZ(px[j], py[j], pz[j]);
  float th = v.Theta();
  float ph = v.Phi();
  cout << "thta:" << th << "phi:" << ph <<endl;
  hparticles[i]->Fill(th, ph);
  */
  hparticles[i]->Fill(theta[j], phi[j]);         
 }
 x_momentum = x_momentum/(float)nparticles;
 y_momentum = y_momentum/(float)nparticles;
 z_momentum = z_momentum/(float)nparticles;
 cout<<"X-momentum ="<< x_momentum << endl;
 cout<<"Y-momentum ="<< y_momentum << endl;
 cout<<"Z-momentum ="<< z_momentum << endl;
 
 c1->cd();
 char *h = new char[100];
 sprintf(h, "Particles_on_thrust_%d.pdf", i);
 hprojection[i]->Draw("colz");
 hparticles[i]->Draw("same");
 gStyle->SetOptStat(0);
 c1->SaveAs(h);
}// end of loop over events


}
