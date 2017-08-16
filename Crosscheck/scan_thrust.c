#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <stdlib.h>
#include <TMath.h>
#include <TVector3.h>
#include <vector>
#include <iostream>

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


//based on code from herwig: http://herwig.hepforge.org/svn/tags/herwig-2-0-beta/Analysis/EventShapes.cc
//ported by A. Baty
void thrust(int isBelle, float x[], float y[], float xch[], float ych[]){

  TString filename;
  if(isBelle) filename="/data/flowex/Datasamples/Belle/output_2_withtheta.root";
  else filename="/data/flowex/Datasamples/LEP2_MAIN/ROOTfiles/cleaned_ALEPH_DATA-all.aleph.root";
  
  //filename ="/mnt/c/Users/Bibek Kumar Pandit/Desktop/Root_Directory/StudyMult/DataFiles/ROOTfiles/cleaned_ALEPH_Data-all";
  
  TFile *f = new TFile(filename.Data());
  TTree *t1 = (TTree*)f->Get("t");
  int n;
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
  
  t1->SetBranchAddress("nParticle",&n);
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

  for (int i=0;i<1337;i++) 
  {  
    if (i%100==0) cout <<i<<"/"<<t1->GetEntries()<<endl;
    t1->GetEntry(i);

    
    TVector3 thrust = TVector3(0,0,0);
    if(n==1){//thrust is just the particle
      thrust = TVector3(px[0],py[0],pz[0]);   
    }
  
    if(n==2){//special case for 2 particles
      if(TMath::Power(px[0],2)+TMath::Power(py[0],2)+TMath::Power(pz[0],2) >= TMath::Power(px[1],2)+TMath::Power(py[1],2)+TMath::Power(pz[1],2)){
        thrust = TVector3(px[0],py[0],pz[0]);
      }
      else{
        thrust = TVector3(px[1],py[1],pz[1]);
      }
    }
  
    if(n==3){//combine lowest 2 magnitude momentum, then use same algo as n=2
      if(TMath::Power(px[0],2)+TMath::Power(py[0],2)+TMath::Power(pz[0],2) >= TMath::Power(px[1],2)+TMath::Power(py[1],2)+TMath::Power(pz[1],2)){
        if(TMath::Power(px[0],2)+TMath::Power(py[0],2)+TMath::Power(pz[0],2) >= TMath::Power(px[2],2)+TMath::Power(py[2],2)+TMath::Power(pz[2],2)){
          thrust = TVector3(px[0],py[0],pz[0]);//0 is largest momentum
        }
        else{
          thrust = TVector3(px[2],py[2],pz[2]);//2 is largest momentum
        } 
      }
      else{
        if(TMath::Power(px[1],2)+TMath::Power(py[1],2)+TMath::Power(pz[1],2) >= TMath::Power(px[2],2)+TMath::Power(py[2],2)+TMath::Power(pz[2],2)){
          thrust = TVector3(px[1],py[1],pz[1]);//1 is largest momentum
        }
        else{
          thrust = TVector3(px[2],py[2],pz[2]);//2 is largest momentum
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

    x[i] = thrust.Theta();
    y[i] = thrust.Phi();
  }

  for (int i=0;i<1337;i++) 
  {  
    if (i%100==0) cout <<i<<"/"<<t1->GetEntries()<<endl;
    t1->GetEntry(i);
    float nTrk = 0;
    for(int t = 0; t<n; t++){
      if(pwflag[t]!=0) continue;
      nTrk++;
    }
    
    TVector3 thrust = TVector3(0,0,0);
    if(nTrk==1){//thrust is just the particle
      for(int t = 0; t<n; t++){
        if(pwflag[t]==0) thrust = TVector3(px[t],py[t],pz[t]);   
      }
    }
  
    if(nTrk==2){//special case for 2 particles
      int n1 = -1, n2 = -1;
      for(int t = 0; t<n; t++){
        if(pwflag[t]==0 && n1==-1){ n1 = t; continue;}
        if(pwflag[t]==0 && n1!=-1) n2 = t;  
      }

      if(TMath::Power(px[n1],2)+TMath::Power(py[n1],2)+TMath::Power(pz[n1],2) >= TMath::Power(px[n2],2)+TMath::Power(py[n2],2)+TMath::Power(pz[n2],2)){
        thrust = TVector3(px[n1],py[n1],pz[n1]);
      }
      else{
        thrust = TVector3(px[n2],py[n2],pz[n2]);
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
          thrust = TVector3(px[n1],py[n1],pz[n1]);//n1 is largest momentum
        }
        else{
          thrust = TVector3(px[n3],py[n3],pz[n3]);//2 is largest momentum
        } 
      }
      else{
        if(TMath::Power(px[n2],2)+TMath::Power(py[n2],2)+TMath::Power(pz[n2],2) >= TMath::Power(px[n3],2)+TMath::Power(py[n3],2)+TMath::Power(pz[n3],2)){
          thrust = TVector3(px[n2],py[n2],pz[n2]);//1 is largest momentum
        }
        else{
          thrust = TVector3(px[n3],py[n3],pz[n3]);//2 is largest momentum
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
      for(int i = 1; i<pVec.size(); i++){//loop through all possible cross products of 2 unique vectors
        for(int j = 0; j<i; j++){
          cross = pVec.at(i).Cross(pVec.at(j));
          TVector3 ptot = TVector3(0,0,0);
 
          for(int k = 0; k<pVec.size(); k++){ //loop through all 3rd particles not used for the cross product
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

    xch[i] = thrust.Theta();
    ych[i] = thrust.Phi();
  }
}

//legacy thrust code written using a scan over theta/phi
void thrustLegacy(int isBelle, float x[], float y[], float xch[], float ych[]){

  TString filename;
  if(isBelle) filename="/data/flowex/Datasamples/Belle/output_2_withtheta.root";
  else filename="/data/flowex/Datasamples/LEP2_MAIN/ROOTfiles/cleaned_ALEPH_DATA-all.aleph.root";

  
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
  Int_t nevent = 1337;
  
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
  filename = "/data/flowex/Datasamples/LEP2_MAIN/ROOTfiles/cleaned_ALEPH_DATA-all.aleph.root";


  
  TFile *f = new TFile(filename.Data());
  TTree *t1 = (TTree*)f->Get("t");
  
  
  
  int nParticle;
  int EventNo;
  int RunNo;
  float Energy;
  float TTheta; //to be added
  float TPhi; //to be added
  float TchTheta;
  float TchPhi;


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

  TFile *g = new TFile("/data/abaty/EpEmStudies/ALEPH_Thrust/thrust_Aug16.root","RECREATE");
 

  TTree *newtree = t1->CloneTree(0); // Do no copy the data yet
  // DO NOT delete the old tree and close the old file
  //add the branch to the new tree and try to fill it
  newtree->Branch("TTheta", &TTheta, "TTheta/F");
  newtree->Branch("TPhi", &TPhi, "TPhi/F");
  newtree->Branch("TchTheta", &TchTheta, "TchTheta/F");
  newtree->Branch("TchPhi", &TchPhi, "TchPhi/F");
  
  
  //int nentries=newtree->GetEntries();
  //Int_t nevent = (Int_t)t1->GetEntries();
  Int_t nevent = 1337;
 
  int nevent_process = nevent;

  //cout<<nentries<<endl;
  //int nentries = 100;
  float x[nevent_process];
  memset(x, 0, nevent*sizeof(float));
  
  float y[nevent_process];
  memset(y, 0, nevent*sizeof(float)); 
  
  float xch[nevent_process];
  memset(xch, 0, nevent*sizeof(float));
  
  float ych[nevent_process];
  memset(ych, 0, nevent*sizeof(float)); 
 
  thrust(0,x,y,xch,ych);

  
  for( int i=0; i < nevent_process; i++){
     t1->GetEntry(i);
     TTheta = x[i];
     TPhi = y[i];
     TchTheta = xch[i];
     TchPhi = ych[i];
     newtree->Fill();
  }
  g->Write();
}
