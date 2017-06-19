//original code is high_low-eta_phi.c

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <iostream>
#include <stdlib.h>
#include <TMath.h>

using namespace std;


#define PI 3.1415926
enum SIMPLEPWFLAG {CHARGED_TRACK, CHARGED_LEPTONS1, CHARGED_LEPTONS2, V0, PHOTON, NEUTRAL_HADRON};

void analysis(int nplots = 10, int maxevt = 0){

  TString filename;
  filename="/mnt/c/Users/Bibek Kumar Pandit/Desktop/Root_Directory/StudyMult/LEP2/ROOTfiles/cleaned_ALEPH_Data2-all.aleph.root";			
  
  TFile *f = new TFile(filename.Data());
  TTree *t1 = (TTree*)f->Get("t");
  Int_t nParticle;
  Float_t pt[50000];
  Float_t eta[50000];
  Float_t theta[50000];
  Float_t pid[50000];
  Float_t phi[50000];
  Float_t mass[50000];
  Float_t pwflag[50000];
  
  t1->SetBranchAddress("nParticle",&nParticle);
  t1->SetBranchAddress("pt",pt);
  t1->SetBranchAddress("eta",eta);
  t1->SetBranchAddress("theta",theta);
  t1->SetBranchAddress("pid",pid);
  t1->SetBranchAddress("phi",phi);
  t1->SetBranchAddress("mass",mass);
  t1->SetBranchAddress("pwflag",pwflag);
  
  TCanvas *c1 = new TCanvas("eta_phi", "eta_phi", 600,  600);
  
  
  
  Int_t nevent = (Int_t)t1->GetEntries();
  int nevent_process = nevent;
  
  if (maxevt > 0 && maxevt < nevent) nevent_process = maxevt;
  int count_h = 0;
  int count_l = 0;
  
    for (Int_t i = 0; i<nevent; i++)
    {
      if (count_h >nplots && count_l>nplots) break;
      
      if (i%10000 == 0) cout<<i<<"/"<<nevent_process<<endl;
      t1->GetEntry(i);
      
      int nparticles = nParticle;
      char *h_2D_name = new char[100];
      if(count_h<nplots)
      {   
        sprintf(h_2D_name, "h_2D_name_%d", i);
        TH2F *h_2D = new TH2F (h_2D_name, "eta_phi of Charged Tracks", 50, -2.5, 2.5, 50, -PI, PI);
        h_2D->GetXaxis()->SetTitle("#eta");
        h_2D->GetYaxis()->SetTitle("#phi");
        
        
        for (int j = 0; j< nparticles; j++)
          {
              if (pwflag[j]==CHARGED_TRACK)
              {
                
                h_2D->SetMarkerStyle(20);
                h_2D->Fill(eta[j], phi[j]);
              }
              
          }
        
        count_h++;
        
        c1->cd();
        char *h = new char[100];
        sprintf(h, "Eta-Phi-Charged_Track_%d.pdf", i);
        h_2D->Draw();
        c1->SaveAs(h);
      }
      
    }
}



void analysis2(int maxevt = 10)
{
  TString filename;
  filename="/mnt/c/Users/Bibek Kumar Pandit/Desktop/Root_Directory/StudyMult/LEP2/ROOTfiles/cleaned_ALEPH_Data2-all.aleph.root";			
  
  TFile *f = new TFile(filename.Data());
  TTree *t1 = (TTree*)f->Get("t");
  Int_t nParticle;
  Float_t eta[50000];
  Float_t theta[50000];
  Float_t phi[50000];
  Float_t pwflag[50000];
  
  t1->SetBranchAddress("nParticle",&nParticle);
  t1->SetBranchAddress("eta",eta);
  t1->SetBranchAddress("theta",theta);
  t1->SetBranchAddress("phi",phi);;
  t1->SetBranchAddress("pwflag",pwflag);
  
  TCanvas *c1 = new TCanvas("eta_phi_all", "eta_phi_all", 600,  600);
  TCanvas *c2 = new TCanvas("eta_phi_charged", "eta_phi_charged", 600,  600);
  TCanvas *c3 = new TCanvas("theta_phi_all", "theta_phi_all", 600,  600);
  TCanvas *c4 = new TCanvas("theta_phi_charged", "theta_phi_charged", 600,  600);
  
  Int_t nevent = (Int_t)t1->GetEntries();
  int nevent_process = nevent;
  
  if (maxevt > 0 && maxevt < nevent) nevent_process = maxevt;
   
  for (Int_t i = 0; i<nevent_process; i++)
  {
   
    
    if (i%10000 == 0) cout<<i<<"/"<<nevent_process<<endl;
    t1->GetEntry(i);
    
    int nparticles = nParticle;
    
    char *h1_name = new char[100];
    char *h2_name = new char[100];
    char *h3_name = new char[100];
    char *h4_name = new char[100];
    
    
    sprintf(h1_name, "h1_name_%d", i);
    TH2F *h1 = new TH2F (h1_name, "eta_phi of All Particles", 50, -2.5, 2.5, 50, -PI, PI);
    h1->GetXaxis()->SetTitle("#eta");
    h1->GetYaxis()->SetTitle("#phi");
    
    sprintf(h2_name, "h2_name_%d", i);
    TH2F *h2 = new TH2F (h2_name, "eta_phi of Charged Tracks", 50, -2.5, 2.5, 50, -PI, PI);
    h2->GetXaxis()->SetTitle("#eta");
    h2->GetYaxis()->SetTitle("#phi");
    
    sprintf(h3_name, "h3_name_%d", i);
    TH2F *h3 = new TH2F (h3_name, "theta_phi of All Particles", 50, 0, PI, 50, -PI, PI);
    h3->GetXaxis()->SetTitle("#theta");
    h3->GetYaxis()->SetTitle("#phi");
    
    sprintf(h4_name, "h4_name_%d", i);
    TH2F *h4 = new TH2F (h4_name, "theta_phi of Charged Tracks", 50, 0, PI, 50, -PI, PI);
    h4->GetXaxis()->SetTitle("#theta");
    h4->GetYaxis()->SetTitle("#phi");
    
    for (int j = 0; j< nparticles; j++)
      {
        h1->SetMarkerStyle(20);
        h1->Fill(eta[j], phi[j]);
        
        h3->SetMarkerStyle(20);
        h3->Fill(theta[j], phi[j]);
        
        if (pwflag[j]==CHARGED_TRACK)
        {
          
          h2->SetMarkerStyle(20);
          h2->Fill(eta[j], phi[j]);
          
          h4->SetMarkerStyle(20);
          h4->Fill(theta[j], phi[j]);
        }
          
      }
    
    
    
    c1->cd();
    char *ha = new char[100];
    sprintf(ha, "Eta-Phi-All_%d.pdf", i);
    h1->Draw();
    c1->SaveAs(ha);
    
    c2->cd();
    char *hb = new char[100];
    sprintf(hb, "Eta-Phi-Charged_%d.pdf", i);
    h2->Draw();
    c2->SaveAs(hb);
    
    c3->cd();
    char *hc = new char[100];
    sprintf(hc, "Theta-Phi-All_%d.pdf", i);
    h3->Draw();
    c3->SaveAs(hc);
    
    c4->cd();
    char *hd = new char[100];
    sprintf(hd, "Theta-Phi-Charged_%d.pdf", i);
    h4->Draw();
    c4->SaveAs(hd);
    
    
    
  }
  
}