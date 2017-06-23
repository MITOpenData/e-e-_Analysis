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

void thrust(int isBelle, float x[]){


  if(isBelle) filename="/data/flowex/Datasamples/Belle/output_2_withtheta.root";
  else filename="/data/flowex/Datasamples/LEP2_MAIN/ROOTfiles/cleaned_ALEPH_Data-all.aleph.root";
  
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
  
    // all entries and fill the histograms
  Int_t nevent = (Int_t)t1->GetEntries();
  
  //int nevent_process = nevent;
  //cout<<"n"<<nevent<<endl;
  int nevent_process = nevent; //59874
  //TH2F*hprojection[nevent_process];
  //float x[nevent_process];
  
    for (Int_t i=0;i<nevent_process;i++) {  
    if (i%10==0) cout <<i<<"/"<<nevent_process<<endl;
    t1->GetEntry(i);
    
    //hprojection[i]=new TH2F(Form("hprojection%d",i),Form("hprojection%d",i),nstepstheta,mintheta,maxtheta,nstepsphi,minphi,maxphi);

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
       
       for ( int j=0;j<nParticle;j++ ) {
         //float pt1 = pt[j];
         //float theta1 = theta[j];
         //float phi1 = phi[j];
         float px1=px[j];
         float py1=py[j];
         float pz1=pz[j];
         particle_vector.SetXYZ(px1,py1,pz1);
         scalar_prod=scalar_prod+fabs(particle_vector.Dot(unity_versor));         
       }//end of loop over particles

       //hprojection[i]->SetBinContent(itheta,iphi,scalar_prod);       
       if (scalar_prod>maxscalar) { maxscalar=scalar_prod; theta_max=mytheta; phi_max=myphi;} 
      }//loop over theta
   }//loop over phi
   //int binmax = hprojection[i]->GetMaximumBin();
   x[i] = theta_max;
   
}// end of loop over events
}



void scan_thrust(){
   TString filename;
  //filename="/mnt/c/Users/Bibek Kumar Pandit/Desktop/Root_Directory/StudyMult/LEP2/ROOTfiles/cleaned_ALEPH_Data-all.aleph.root";
  filename = "/data/flowex/Datasamples/LEP2_MAIN/ROOTfiles/cleaned_ALEPH_Data-all.aleph.root";
  
  TFile *f = new TFile(filename.Data());
  TTree *t1 = (TTree*)f->Get("t");
  
  
  
  int nParticle;
  int EventNo;
  int RunNo;
  float Energy;
  float tAngle;

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
  

  
  

  //TFile *g = new TFile("/mnt/c/Users/Bibek Kumar Pandit/Desktop/Root_Directory/StudyMult/LEP2/ROOTfiles/cleaned_ALEPH_Data2-all.aleph.root", "RECREATE");
  TFile *g = new TFile("/data/flowex/Datasamples/LEP2_MAIN/ROOTfiles/cleaned_ALEPH_Data2-all.aleph.roo","RECREATE");
 

  TTree *newtree = t1->CloneTree(0); // Do no copy the data yet
  // DO NOT delete the old tree and close the old file
  //add the branch to the new tree and try to fill it
  newtree->Branch("tAngle", &tAngle, "tAngle/F");
  
  
  //int nentries=newtree->GetEntries();
  Int_t nevent = (Int_t)t1->GetEntries();
    int nevent_process = nevent;
  //cout<<nentries<<endl;
  //int nentries = 100;
  float x[nevent_process];
  memset(x, 0, nevent*sizeof(float) );
  cout<<"elo"<<endl; 
    
  thrust(0,x);
  
  for( int i=0; i < nevent_process; i++){
     t1->GetEntry(i);
     tAngle = x[i];
     newtree->Fill();
  }
  g->Write();
}
