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


double dphi(double phi1,double phi2)
{
    double a=phi1-phi2;
    while (a<-PI) a+=2*PI;
    while (a>PI) a-=2*PI;
    
    if (a<-PI/2) a=2*PI+a;
    return a;
}

<<<<<<< HEAD
void analysis(int isBelle=1, int maxevt=0,int mult = 0, int mult_upper_bound =0,int nbin=40,bool verbose=0){
    
=======
void analysis(int isBelle=1, int maxevt=0,int mult=100, int mult_upper_bound = 1000, int nbin=40,bool verbose=0){
  
>>>>>>> origin/master
  
  TString filename;
  if(isBelle) filename="/mnt/c/Users/Bibek Kumar Pandit/Desktop/Root_Directory/StudyMult/LEP2/ROOTfiles/cleaned_ALEPH_DATA-all.aleph.root";
  //if(isBelle) filename="/Users/anthony/Documents/StudyMult/LEP2/ROOTfiles/cleaned_ALEPH_DATA-all.aleph.root";
  //else filename="../LEP2dataMarcello/myALEPH.root";
  //else filename="../LEP2dataMarcello/ROOTfiles/final_ALEPH.root";
  
  TFile *f = new TFile(filename.Data());
  TTree *t1 = (TTree*)f->Get("t");
  Int_t nParticle;
  Float_t pt[50000];
  Float_t eta[50000];
  Float_t pid[50000];
  Float_t phi[50000];
  Float_t mass[50000];
  Float_t pwflag[50000];
  
  t1->SetBranchAddress("nParticle",&nParticle);
  t1->SetBranchAddress("pt",pt);
  t1->SetBranchAddress("eta",eta);
  t1->SetBranchAddress("pid",pid);
  t1->SetBranchAddress("phi",phi);
  t1->SetBranchAddress("mass",mass);
  t1->SetBranchAddress("pwflag",pwflag);

  TFile *f_mix = new TFile(filename.Data());
  TTree *t1_mix = (TTree*)f_mix->Get("t");
  Int_t nParticle_mix;
  Float_t pt_mix[50000];
  Float_t eta_mix[50000];
  Float_t pid_mix[50000];
  Float_t phi_mix[50000];
  Float_t mass_mix[50000];
  Float_t pwflag_mix[50000];
  
  t1_mix->SetBranchAddress("nParticle",&nParticle_mix);
  t1_mix->SetBranchAddress("pt",pt_mix);
  t1_mix->SetBranchAddress("eta",eta_mix);
  t1_mix->SetBranchAddress("pid",pid_mix);
  t1_mix->SetBranchAddress("phi",phi_mix);
  t1_mix->SetBranchAddress("mass",mass_mix);
  t1_mix->SetBranchAddress("pwflag",pwflag_mix);


  // two histograms
  double detaRange = 2.4;
  double normalization = detaRange*2/nbin*2*3.14159/nbin;
  TH2F *h_2D = new TH2F ( "h_2D", "eta-phi of all particles ",nbin, -detaRange, detaRange,nbin, -3.1416/2., 3.1416*1.5);
  TH2F *h_2Dmix = new TH2F ( "h_2Dmix", "eta-phi of all particles ",nbin, -detaRange, detaRange,nbin, -3.1416/2., 3.1416*1.5);
  TH2F *h_ratio = new TH2F ( "h_ratio", "eta-phi of all particles ", nbin, -detaRange, detaRange,nbin, -3.1416/2.,3.1416*1.5);
  
  h_2D->Sumw2();
  h_2Dmix->Sumw2();
  h_ratio->Sumw2();
  

  // all entries and fill the histograms
  Int_t nevent = (Int_t)t1->GetEntries();
  
  int nevent_process = nevent;
  if( maxevt>0 && maxevt<nevent ) nevent_process = maxevt;  
  
    double averageN=0;
    double nEventProcessed=0;
    double nEventInMultBin=0;
    for (Int_t i=0;i<nevent_process;i++) {
  
    if (i%10000==0) cout <<i<<"/"<<nevent_process<<endl;
    t1->GetEntry(i);
    
    int nparticles = nParticle;
    if (verbose) cout<<"nparticles="<<nparticles<<endl;
  
    int nparticles2=-1000;
    int selected=i+1;  //questo diventa il numero di evento estratto nel file +1
    int flag=0;    //definisco una flag a zero
    
    // Yen-Jie: cut on maximum number of particles to avoid infinite loop, for the moment it is 100
    if (nparticles>mult_upper_bound) continue;
    

    //double N=0;
    double ptMin=0.4;
    double ptMax=4;
    
    
    //int count = 0;
    
    // calculate the number of tracks in the passing selection
    double N=0;
    for (int p = 0; p < 10; p++)
    {
      
      for ( int j=0;j<nparticles;j++ ) {
        float pt1 = pt[j];
        int pid1 = pid[j];
        int pwflag1 = pwflag[j];
        if (pwflag1 != CHARGED_TRACK) continue;
        if(pt1<ptMin||pt1>ptMax) continue;
        N++;
      }
      averageN+=N;
      nEventProcessed++;

      if (nparticles<mult) continue;
      //if (N<mult) continue;
      
      nEventInMultBin++;

      // find a mixed event
      //	cout <<N<<endl;
      while ((fabs(nparticles2-nparticles)>5&&nparticles<mult_upper_bound)||i==selected){
         //cout <<nparticles<<" "<<selected<<endl;
         selected++;
         if (selected>nevent_process&&flag==1) break;
         if (selected>nevent_process) flag=1;
         selected = selected % nevent_process;
         t1_mix->GetEntry ( selected );
         nparticles2= nParticle_mix;
      }
      
      double N2=0;
      // calculate the number of tracks in the mixed event passing selection
      for ( int j=0;j<nparticles;j++ ) {
        float pt1 = pt_mix[j];
        int pid1 = pid_mix[j];
        int pwflag1 = pwflag_mix[j];
        if (pwflag1 != CHARGED_TRACK) continue;
        if(pt1<ptMin||pt1>ptMax) continue;
        N2++;
      }
      
      
      for ( int j=0;j<nparticles;j++ ) {
        int pid1 = pid[j];
        float eta1 = eta[j];
        float phi1 = phi[j];
        float pt1 = pt[j];
        float mass1 = mass[j];
        float pwflag1 = pwflag[j];
        if (pwflag1 != CHARGED_TRACK) continue;
        if(pt1<ptMin||pt1>ptMax) continue;
        
        // Signal loop, calculate S correlation function
        if (p == 0)
        {
          for ( int k=j+1;k<nparticles;k++ ) {
            int pid2 = pid[k];
            float eta2 = eta[k];
            float phi2 = phi[k];
            float pt2 = pt[k];
            float mass2 = mass[k];
            float pwflag2 = pwflag[k];
            if (pwflag2 != CHARGED_TRACK) continue;
            if(pt2<ptMin||pt2>ptMax) continue;
            
            h_2D->Fill(eta1-eta2,dphi(phi1,phi2),1./N);    
            h_2D->Fill(eta1-eta2,dphi(phi2,phi1),1./N);    
            h_2D->Fill(eta2-eta1,dphi(phi1,phi2),1./N);    
            h_2D->Fill(eta2-eta1,dphi(phi2,phi1),1./N);    
          }//end of second loop 
        }

        // Background loop, calculate B correlation function from mixed event
        for ( int k=0;k<nparticles2;k++ ) {
          int pidmix = pid_mix[k];
          float etamix = eta_mix[k];
          float phimix = phi_mix[k];
          float ptmix = pt_mix[k];
          float massmix = mass_mix[k];
          float pwflagmix = pwflag_mix[k];
          if (pwflagmix != CHARGED_TRACK) continue;
          if(ptmix<ptMin||ptmix>ptMax) continue;
            
          h_2Dmix->Fill(eta1-etamix,dphi(phi1,phimix));
          h_2Dmix->Fill(eta1-etamix,dphi(phimix,phi1));
          h_2Dmix->Fill(etamix-eta1,dphi(phi1,phimix));
          h_2Dmix->Fill(etamix-eta1,dphi(phimix,phi1));
        }//end of second loop 

    }// end of first loop
    }
    // NOW DIVIDE
    h_2Dmix->Scale(1./N);
    
  }// end of loop over events
  
  
  h_2Dmix->Scale(1./nEventInMultBin);
  h_2D->Scale(1./nEventInMultBin);
  
  averageN=averageN/nEventProcessed;
  cout <<"Average N = "<<averageN<<endl;
  double ratio;
  double errrel_ratio;
  double errrel_num;
  double errrel_den;
  
    // calculate the  correlation function
  
  double b00_x=h_2Dmix->GetXaxis()->FindBin(0.);
  double b00_y=h_2Dmix->GetYaxis()->FindBin(0.);
  double B00=h_2Dmix->GetBinContent(b00_x,b00_y);
  double errrel_B00=h_2Dmix->GetBinError(b00_x,b00_y)/B00;
  
  cout<<"value of B(0,0)="<<B00<<endl;


  cout<<"x axis "<<h_2Dmix->GetXaxis()->GetBinCenter(b00_x);
  cout<<"y axis "<<h_2Dmix->GetYaxis()->GetBinCenter(b00_y);
  
  for (int x=0;x<=h_2D->GetNbinsX();x++){
     for (int y=0;y<=h_2D->GetNbinsY();y++){
        if(h_2Dmix->GetBinContent(x,y)>0){
          ratio=B00*(h_2D->GetBinContent(x,y)/h_2Dmix->GetBinContent(x,y));
          errrel_num=h_2D->GetBinError(x,y)/h_2D->GetBinContent(x,y);
          errrel_den=h_2Dmix->GetBinError(x,y)/h_2Dmix->GetBinContent(x,y);
          errrel_ratio=TMath::Sqrt(errrel_num*errrel_num+errrel_den*errrel_den+errrel_B00*errrel_B00);
          h_ratio->SetBinContent(x,y,ratio);
          h_ratio->SetBinError(x,y,ratio*errrel_ratio);
        }
     }
  }

  int etaranges[4]={0,1,2,3};
  int minbin,maxbin;
  
  TH1D*h_deltaphi[3];

  for (int i=0;i<3;i++){
    minbin =  h_ratio->GetXaxis()->FindBin(etaranges[i]);
    maxbin =  h_ratio->GetXaxis()->FindBin(etaranges[i+1]);
    //h_deltaphi[i]  = (TH1D*) h_ratio->ProjectionY(Form("h_deltaphi_etamin%d_max%d",etaranges[i],etaranges[i+1]),minbin,maxbin);
    h_deltaphi[i]  = (TH1D*) h_ratio->ProjectionY(Form("h_deltaphi%d",i),minbin,maxbin);
    h_deltaphi[i]->Sumw2();
    //h_deltaphi[i]->SetName(Form("h_deltaphi_etamin%d_max%d",etaranges[i],etaranges[i+1]));
    h_deltaphi[i]->GetZaxis()->CenterTitle();
    h_deltaphi[i]->GetXaxis()->SetTitle("#Delta#phi");
    h_deltaphi[i]->Scale(1./(maxbin-minbin+1));
  }  


  cout<<"error"<<h_deltaphi[0]->GetBinError(5)<<endl;

  h_ratio->GetXaxis()->SetTitle("#Delta#eta");
  h_ratio->GetYaxis()->SetTitle("#Delta#phi");
  

  TFile*fout=new TFile(Form("ROOTfiles/myoutput_isBelle%d_minMult%d.root",isBelle,mult),"recreate");
  fout->cd();
  h_2D->Write();
  h_2Dmix->Write();
  h_ratio->Write();
  TH1D *h_proj = h_ratio->ProjectionY("",2,4);
  h_proj->SetTitle(Form("Multiplicity Between %d and %d",mult,mult_upper_bound));
  h_proj->Write();
  for (int i=0;i<3;i++) h_deltaphi[i]->Write();
  fout->Close();
  delete fout;
    
    
  TCanvas *c1 = new TCanvas("c1","",600,600);
    c1->SetTheta(60.839);
    c1->SetPhi(38.0172);
  h_ratio->SetTitle(Form("Ratio for Multiplicity between %d and %d",mult,mult_upper_bound));
  h_ratio->Draw("LEGO2");
  c1->SaveAs(Form("ratio%d_%d.pdf",mult,mult_upper_bound));
  
  
  }

void goofy()
{
  int low[2] = {0, 40};
  int high[2] = {20, 100};
  for (int i = 0; i<3; i++){
    analysis(1, 0, low[i], high[i], 40, 0);
  }
<<<<<<< HEAD



void bar()
{
    int mult[1][2] = {{0,20}};//,{20,40},{40,1000},{50,1000},{60,1000},{70,1000},{80,1000},{90,1000}};
    cout<<mult[0]<<"/n";
    analysis(1, 0,0,20,40,0);
    /**
    for (int i =0; i<1;i++)
    {
        analysis(1, 0,mult[i][0],mult[i][1],40,0);
    }**/
    
=======
  //return 0;
>>>>>>> origin/master
}
