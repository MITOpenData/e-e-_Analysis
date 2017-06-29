#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <iostream>
#include <stdlib.h>
#include <TMath.h>
#include <TVector3.h>

using namespace std;


#define PI 3.1415926
enum SIMPLEPID {PHOTON, ELECTRON, PION, MUON, KAON, PROTON};


double dphi(double phi1,double phi2)
{
    double a=phi1-phi2;
    while (a<-PI) a+=2*PI;
    while (a>PI) a-=2*PI;
    
    if (a<-PI/2) a=2*PI+a;
    return a;
}

void thrust(int isBelle = 0){
  
  TString filename;
  if(isBelle) filename="/data/flowex/Datasamples/Belle/output_2_withtheta.root"; 			
  else filename="/data/flowex/Datasamples/LEP2_MAIN/ROOTfiles/cleaned_ALEPH_Data2-all.aleph.root";
  
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
  t1->SetBranchAddress("px",pt);
  t1->SetBranchAddress("py",pt);
  t1->SetBranchAddress("pz",pt);
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
  const int nevent_process = 100;
  TH2F *hprojection[nevent_process];
  
    for (Int_t i=0;i<nevent_process;i++) {  
    if (i%10==0) cout <<i<<"/"<<nevent_process<<endl;
    t1->GetEntry(i);
    
    hprojection[i]=new TH2F(Form("hprojection%d",i),Form("hprojection%d",i),nstepstheta,mintheta,maxtheta,nstepsphi,minphi,maxphi);

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
         float pt1 = pt[j];
         float theta1 = theta[j];
         float phi1 = phi[j];
         float px1=px[j];
         float py1=py[j];
         float pz1=pz[j];
         /*
         float px1=pt1*sin(theta1)*cos(phi1);
         float py1=pt1*sin(theta1)*sin(phi1);
         float pz1=pt1*cos(theta1);
         */
         
         particle_vector.SetXYZ(px1,py1,pz1);
         scalar_prod=scalar_prod+fabs(particle_vector.Dot(unity_versor));         
       }//end of loop over particles

       hprojection[i]->SetBinContent(itheta,iphi,scalar_prod);       
       if (scalar_prod>maxscalar) { maxscalar=scalar_prod; theta_max=mytheta; phi_max=myphi;} 
      }//loop over theta
   }//loop over phi
   
   cout<<"max theta="<<theta_max<<endl;
   cout<<"max phi="<<phi_max<<endl;
   cout<<"max scaler="<<maxscalar<<endl;
   cout<<"max histo="<<hprojection[i]->GetMaximum()<<endl;
   hprojection[i]->Draw("colz");
}// end of loop over events



if (isBelle) TFile*fout=new TFile("data/flowex/Datasamples/Belle/lookatthisfile.root","recreate"); 
else TFile*fout=new TFile("/data/flowex/Datasamples/LEP2_MAIN/ROOTfiles/lookatthisfile.root","recreate"); 
fout->cd();
for (int m=0;m<nevent_process;m++)hprojection[m]->Write();
fout->Close();

}
