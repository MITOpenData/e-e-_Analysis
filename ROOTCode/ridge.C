#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>

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

void analysis(int maxevt=0,int mult=50,int nbin=15){

  TFile *f = new TFile("../Inputs/output-2.root");
  TTree *t1 = (TTree*)f->Get("t");
  Int_t nParticle;
  Float_t pt[50000];
  Float_t eta[50000];
  Float_t pid[50000];
  Float_t phi[50000];
  Float_t mass[50000];
  
  t1->SetBranchAddress("nParticle",&nParticle);
  t1->SetBranchAddress("pt",pt);
  t1->SetBranchAddress("eta",eta);
  t1->SetBranchAddress("pid",pid);
  t1->SetBranchAddress("phi",phi);
  t1->SetBranchAddress("mass",mass);

  TFile *f_mix = new TFile("../Inputs/output-2.root");
  TTree *t1_mix = (TTree*)f_mix->Get("t");
  Int_t nParticle_mix;
  Float_t pt_mix[50000];
  Float_t eta_mix[50000];
  Float_t pid_mix[50000];
  Float_t phi_mix[50000];
  Float_t mass_mix[50000];
  
  t1_mix->SetBranchAddress("nParticle",&nParticle_mix);
  t1_mix->SetBranchAddress("pt",pt_mix);
  t1_mix->SetBranchAddress("eta",eta_mix);
  t1_mix->SetBranchAddress("pid",pid_mix);
  t1_mix->SetBranchAddress("phi",phi_mix);
  t1_mix->SetBranchAddress("mass",mass_mix);

  // two histograms
  TH2F *h_2D = new TH2F ( "h_2D", "eta-phi of all particles ",nbin, -3.5, 3.5,nbin, -3.1416/2., 3.1416*1.5);
  TH2F *h_2Dmix = new TH2F ( "h_2Dmix", "eta-phi of all particles ",nbin, -3.5, 3.5,nbin, -3.1416/2., 3.1416*1.5);
  TH2F *h_ratio = new TH2F ( "h_ratio", "eta-phi of all particles ", nbin, -3.5, 3.5,nbin, -3.1416/2.,3.1416*1.5);


  // all entries and fill the histograms
  Int_t nevent = (Int_t)t1->GetEntries();
  
  int nevent_process = nevent;
  if( maxevt>0 && maxevt<nevent ) nevent_process = maxevt;  
  
  for (Int_t i=0;i<nevent_process;i++) {
    if (i%1000==0) cout <<i<<"/"<<nevent_process<<endl;
    t1->GetEntry(i);
    
    int nparticles = nParticle;
    if (nparticles<mult) continue;
  
    int nparticles2=-1000;
    int selected=i+1;  //questo diventa il numero di evento estratto nel file +1
    int flag=0;    //definisco una flag a zero
    if (nparticles>500) continue;
    while ((abs(nparticles2-nparticles)>5&&nparticles<100)||i==selected){
       //cout <<nparticles<<" "<<selected<<endl;
       selected++;
       if (selected>nevent_process&&flag==1) break;
       if (selected>nevent_process) flag=1;
       selected = selected % nevent_process;
       t1_mix->GetEntry ( selected );
       nparticles2= nParticle_mix;
    }

    for ( int j=0;j<nparticles;j++ ) {
    
      int pid1 = pid[j];
      if (pid1!=PION&&pid1!=PROTON&&pid1!=KAON) continue;
      float eta1 = eta[j];
      float phi1 = phi[j];
      float pt1 = pt[j];
      float mass1 = mass[j];
      if(pt1<0||pt1>4) continue;
      
      //cout<<"particle1, pid="<<pid1<<", eta="<<eta1<<", phi="<<phi1<<", eta="<<eta1<<", pt="<<pt1<<", mass="<<mass1<<endl;;
      
      for ( int k=j+1;k<nparticles;k++ ) {
        int pid2 = pid[k];
        if (pid2!=PION&&pid2!=PROTON&&pid2!=KAON) continue;
        float eta2 = eta[k];
        float phi2 = phi[k];
        float pt2 = pt[k];
        float mass2 = mass[k];
        if(pt2<0||pt2>4) continue;
        
        h_2D->Fill(eta1-eta2,dphi(phi1,phi2));    
        h_2D->Fill(eta1-eta2,dphi(phi2,phi1));    
        h_2D->Fill(eta2-eta1,dphi(phi1,phi2));    
        h_2D->Fill(eta2-eta1,dphi(phi2,phi1));    
        ///cout<<"Fill"<<endl;
        //cout<<"particle2, pid="<<pid2<<", eta="<<eta2<<", phi="<<phi2<<", eta="<<eta2<<", pt="<<pt2<<", mass="<<mass2<<endl;;
      }//end of second loop 


      for ( int k=0;k<nparticles2;k++ ) {
        int pidmix = pid_mix[k];
        if (pidmix!=PION&&pidmix!=PROTON&&pidmix!=KAON) continue;
        float etamix = eta_mix[k];
        float phimix = phi_mix[k];
        float ptmix = pt_mix[k];
        float massmix = mass_mix[k];
        if(ptmix<0||ptmix>4) continue;
        
        h_2Dmix->Fill(eta1-etamix,dphi(phi1,phimix));    
        h_2Dmix->Fill(eta1-etamix,dphi(phimix,phi1));    
        h_2Dmix->Fill(etamix-eta1,dphi(phi1,phimix));    
        h_2Dmix->Fill(etamix-eta1,dphi(phimix,phi1));    
        //cout<<"particle2, pid="<<pid2<<", eta="<<eta2<<", phi="<<phi2<<", eta="<<eta2<<", pt="<<pt2<<", mass="<<mass2<<endl;;
      }//end of second loop 

    }// end of first loop
  }// end of loop over events

  for (int x=0;x<=h_2D->GetNbinsX();x++){
  for (int y=0;y<=h_2D->GetNbinsY();y++){
     if (h_2Dmix->GetBinContent(x,y)>0) {h_ratio->SetBinContent(x,y,h_2D->GetBinContent(x,y)/h_2Dmix->GetBinContent(x,y));}
  }
  }

  TFile*fout=new TFile("myoutput.root","recreate");
  fout->cd();
  h_2D->Write();
  h_2Dmix->Write();
  h_ratio->Write();
  fout->Close();
  delete fout;
  
  }