#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <stdlib.h>
#include <TMath.h>
#include <TVector3.h>


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



void thrust(int isBelle, float x[], float y[]){

  TString filename;
  if(isBelle) filename="/data/flowex/Datasamples/Belle/output_2_withtheta.root";
  else filename="/data/flowex/Datasamples/LEP2_MAIN/ROOTfiles/cleaned_ALEPH_Data-all.aleph.root";

  
  //filename ="/mnt/c/Users/Bibek Kumar Pandit/Desktop/Root_Directory/StudyMult/DataFiles/ROOTfiles/cleaned_ALEPH_Data-all";
  
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
  
  int nevent_process = nevent; //59874
  
  for (Int_t i=0;i<nevent_process;i++) 
  {  
    if (i%10==0) cout <<i<<"/"<<nevent_process<<endl;
    t1->GetEntry(i);
    
    int nparticles = nParticle;
    
    maxscalar=-1;

    for ( int itheta=0; itheta<nstepstheta; itheta++ ) 
    {
      mytheta=sizetheta*(double)itheta;
      myphi = -3.14;
      for ( int iphi=0; iphi<nstepsphi; iphi++ ) 
      {

        x_versor=sin(mytheta)*cos(myphi);
        y_versor=sin(mytheta)*sin(myphi);
        z_versor=cos(mytheta);
        
        unity_versor.SetXYZ(x_versor,y_versor,z_versor);
        double scalar_prod;
        
        scalar_prod=0.;
       
        for ( int j=0;j<nparticles;j++ ) 
        {         
          float px1=px[j];
          float py1=py[j];
          float pz1=pz[j];
          particle_vector.SetXYZ(px1,py1,pz1);
          scalar_prod = scalar_prod + fabs(particle_vector.Dot(unity_versor));         
        }//end of loop over particles

        myphi+=sizephi;
        if (scalar_prod>maxscalar) { maxscalar=scalar_prod; theta_max=mytheta; phi_max=myphi;} 
      }//loop over phi
    }//loop over theta
    //int binmax = hprojection[i]->GetMaximumBin();
    x[i] = theta_max;
    y[i] = phi_max;
  }// end of loop over events
}




void scan_thrust(){
  TString filename;
  //filename="/mnt/c/Users/Bibek Kumar Pandit/Desktop/Root_Directory/StudyMult/DataFiles/ROOTfiles/cleaned_ALEPH_Data-all.aleph.root";
  filename = "/data/flowex/Datasamples/LEP2_MAIN/ROOTfiles/cleaned_ALEPH_Data-all.aleph.root";


  
  TFile *f = new TFile(filename.Data());
  TTree *t1 = (TTree*)f->Get("t");
  
  
  
  int nParticle;
  int EventNo;
  int RunNo;
  float Energy;
  float TTheta; //to be added
  float TPhi; //to be added

  Float_t px[100000];
  Float_t py[100000];
  Float_t pz[100000];
  Float_t pt[100000];
  Float_t eta[100000];
  Float_t theta[100000];
  Float_t phi[100000];
  Float_t mass[100000];
  Float_t charge[10000];
  Float_t pwflag[10000];
  Float_t pid[10000];
  
  t1->SetBranchAddress("nParticle",&nParticle);
  t1->SetBranchAddress("EventNo",&EventNo);
  t1->SetBranchAddress("RunNo",&RunNo);
  t1->SetBranchAddress("Energy",&Energy);
  t1->SetBranchAddress("px",px);
  t1->SetBranchAddress("py",py);
  t1->SetBranchAddress("pz",pz);
  t1->SetBranchAddress("pt",pt);
  t1->SetBranchAddress("mass",mass);
  t1->SetBranchAddress("eta",eta);
  t1->SetBranchAddress("theta",theta);
  t1->SetBranchAddress("phi",phi);
  t1->SetBranchAddress("charge",charge);
  t1->SetBranchAddress("pwflag", pwflag);
  t1->SetBranchAddress("pid",pid);
  

  
  


  //TFile *g = new TFile("/mnt/c/Users/Bibek Kumar Pandit/Desktop/Root_Directory/StudyMult/DataFiles/ROOTfiles/cleaned_ALEPH_Data2-all.aleph.root", "RECREATE");

  TFile *g = new TFile("/data/flowex/Datasamples/LEP2_MAIN/ROOTfiles/cleaned_ALEPH_Data2-all.aleph.root","RECREATE");


 

  TTree *newtree = t1->CloneTree(0); // Do no copy the data yet
  // DO NOT delete the old tree and close the old file
  //add the branch to the new tree and try to fill it
  newtree->Branch("TTheta", &TTheta, "TTheta/F");
  newtree->Branch("TPhi", &TPhi, "TPhi/F");
  
  
  //int nentries=newtree->GetEntries();
  Int_t nevent = (Int_t)t1->GetEntries();

  int nevent_process = nevent;

  //cout<<nentries<<endl;
  //int nentries = 100;
  float x[nevent_process];
  memset(x, 0, nevent*sizeof(float));
  
  float y[nevent_process];
  memset(y, 0, nevent*sizeof(float)); 
 
  thrust(0,x,y);

  
  for( int i=0; i < nevent_process; i++){
     t1->GetEntry(i);
     TTheta = x[i];
     TPhi = y[i];
     newtree->Fill();
  }
  g->Write();
}
