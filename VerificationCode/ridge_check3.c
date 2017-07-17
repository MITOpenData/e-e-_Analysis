#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>
#include <stdlib.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TFormula.h>
#include "fourier.h"
using namespace std;


#define PI 3.1415926
//enum SIMPLEPID {PHOTON, ELECTRON, PION, MUON, KAON, PROTON};
enum SIMPLEPWLAG {CHARGED_TRACK, CHARGED_LEPTONS1, CHARGED_LEPTONS2, V0, PHOTON, NEUTRAL_HADRON};

double dphi(double phi1,double phi2)
{
    double a=phi1-phi2;
    while (a<-PI) a+=2*PI;
    while (a>PI) a-=2*PI;
    
    if (a<-PI/2) a=2*PI+a;
    return a;
}

void analysis(int isBelle = 0,int isThrust = 0, int maxevt = 0,int mult_low = 0,int mult_high = 10,int nbin = 50,bool verbose = 0,int num_runs= 1){
  
  TString filename;
  //if(isBelle) filename="/data/flowex/Datasamples/Belle/output_2_withtheta.root";/Users/anthony/Desktop/ROOTUsersGuideLetter.pdf
  //else filename="/data/flowex/Datasamples/LEP2_MAIN/ROOTfiles/cleaned_ALEPH_Data2-all.aleph.root";
  
  //filename = "/Users/anthony/documents/StudyMult/DataFiles/ROOTfiles/cleaned_ALEPH_Data2-all.aleph.root";
  filename = "/mnt/c/Users/Bibek Kumar Pandit/Desktop/Root_Directory/StudyMult/DataFiles/ROOTfiles/cleaned_ALEPH_Data2-all.aleph.root";
  
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
    
  Float_t px[100000];
  Float_t py[100000];
  Float_t pz[100000];
  Float_t TTheta;
  Float_t TPhi;
    
  t1->SetBranchAddress("nParticle",&nParticle);
  t1->SetBranchAddress("pt",pt);
  t1->SetBranchAddress("eta",eta);
  t1->SetBranchAddress("theta",theta);
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
  Float_t theta_mix[50000];
  Float_t pid_mix[50000];
  Float_t phi_mix[50000];
  Float_t mass_mix[50000];
  Float_t pwflag_mix[50000];
  
  Float_t px_mix[100000];
  Float_t py_mix[100000];
  Float_t pz_mix[100000];
  Float_t TTheta_mix;
  Float_t TPhi_mix;
    
  t1_mix->SetBranchAddress("nParticle",&nParticle_mix);
  t1_mix->SetBranchAddress("pt",pt_mix);
  t1_mix->SetBranchAddress("eta",eta_mix);
  t1_mix->SetBranchAddress("theta",theta_mix);
  t1_mix->SetBranchAddress("pid",pid_mix);
  t1_mix->SetBranchAddress("phi",phi_mix);
  t1_mix->SetBranchAddress("mass",mass_mix);
  t1_mix->SetBranchAddress("pwflag",pwflag_mix);

  t1_mix->SetBranchAddress("px",px_mix);
  t1_mix->SetBranchAddress("py",py_mix);
  t1_mix->SetBranchAddress("pz",pz_mix);
  t1_mix->SetBranchAddress("TTheta",&TTheta_mix);
  t1_mix->SetBranchAddress("TPhi",&TPhi);
  // two histograms
  
  double detaRange = 4;//1.8
  double dthetaRange = PI;
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
    
    
    // S calculation using multiplicity cut
    for (Int_t i=0;i<nevent_process;i++) {
      
    if (i%10000==0) cout <<i<<"/"<<nevent_process<<endl;
    t1->GetEntry(i);
    
    int nparticles = nParticle;
    float ttheta = TTheta;
    float tphi = TPhi;
    if (verbose) cout<<"nparticles="<<nparticles<<endl;
  
    int nparticles2=-1000;
    float ttheta_mix;
    float tphi_mix;
    int selected=i+1;  //questo diventa il numero di evento estratto nel file +1
    int flag=0;    //definisco una flag a zero
    
    // Yen-Jie: cut on maximum number of particles to avoid infinite loop, for the moment it is 100
    if (nparticles>mult_high) continue;
    

    double N=0;
    double ptMin=0.1; //0.4
    double ptMax=4;
        for (int p = 0;p<num_runs;p++)
        {
    // calculate the number of tracks in the passing selection
    for ( int j=0;j<nparticles;j++ ) {
      float pt1 = pt[j];
      int pid1 = pid[j];
      int pwflag1 = pwflag[j];
      //if (pid1!=PION&&pid1!=PROTON&&pid1!=KAON&&!(!isBelle&&pwflag1==0)) continue;
      if (pwflag1!=CHARGED_TRACK)continue;
      if(pt1<ptMin||pt1>ptMax) continue;
      N++;
    }
    
    averageN+=N;
    nEventProcessed++;

    //if (nparticles<mult_low) continue;
    if (N<mult_low) continue;
    
    nEventInMultBin++;

    // find a mixed event
    while ((fabs(nparticles2-nparticles)>5&&nparticles<mult_high)||i==selected){
       selected++;
       if (selected>nevent_process&&flag==1) break;
       if (selected>nevent_process) flag=1;
       selected = selected % nevent_process;
       t1_mix->GetEntry ( selected );
       nparticles2= nParticle_mix;
       ttheta_mix = TTheta_mix;
       tphi_mix =TPhi_mix;
    }
    
    double N2=0;
    // calculate the number of tracks in the mixed event passing selection
    for ( int j=0;j<nparticles;j++ ) {
      float pt1 = pt_mix[j];
      int pid1 = pid_mix[j];
      int pwflag1 = pwflag_mix[j];
      if (pwflag1!=CHARGED_TRACK) continue;
      //if (pid1!=PION&&pid1!=PROTON&&pid1!=KAON&&!(!isBelle&&pwflag1==0)) continue;
      if(pt1<ptMin||pt1>ptMax) continue;
      N2++;
    }
    
    TVector3 thrust(1, 1, 1);
    thrust.SetTheta(ttheta);    
    thrust.SetPhi(tphi);

    
    TVector3 rotVec;
    rotVec.SetX(-thrust.Y());
    rotVec.SetY(thrust.X());
    
    TVector3 thrust_mix(1, 1, 1);
    thrust_mix.SetTheta(ttheta_mix);
    thrust_mix.SetPhi(tphi_mix);
    
    TVector3 rotVec_mix;
    rotVec_mix.SetX(-thrust_mix.Y());
    rotVec_mix.SetY(thrust_mix.X());
    
    for ( int j=0;j<nparticles;j++ ) {
      int pid1 = pid[j];
      float eta1 = eta[j];
      float phi1 = phi[j];
      float pt1 = pt[j];
      float mass1 = mass[j];
      float pwflag1 = pwflag[j];
        float theta1 = theta[j];
      if (isThrust == 1)
      {
          TVector3 v1;
          float px1 = px[j];
          float py1 = py[j];
          float pz1 = pz[j];
          v1.SetXYZ(px1,py1,pz1);
          v1.Rotate(-ttheta, rotVec);
          phi1 = v1.Phi();
          eta1 = v1.PseudoRapidity();
          pt1 = v1.Pt();
          theta1 = v1.Theta();
          //
      }
      //cout<<phi1<<endl;
      if(pwflag1!=CHARGED_TRACK)continue;
      //if (pid1!=PION&&pid1!=PROTON&&pid1!=KAON&&!(!isBelle&&pwflag1==0)) continue;
      if(pt1<ptMin||pt1>ptMax) continue;
      if(p==0)
        {
      // Signal loop, calculate S correlation function
      for ( int k=j+1;k<nparticles;k++ ) {
        int pid2 = pid[k];
        float eta2 = eta[k];
        float phi2 = phi[k];
        float pt2 = pt[k];
        float mass2 = mass[k];
        float pwflag2 = pwflag[k];
          float theta2 = theta[k];
        if (isThrust == 1)
        {
            TVector3 v2;
            float px2 = px[k];
            float py2 = py[k];
            float pz2 = pz[k];
            v2.SetXYZ(px2,py2,pz2);
            v2.Rotate(-ttheta, rotVec);
            phi2 = v2.Phi();
            eta2 = v2.PseudoRapidity();
            pt2 = v2.Pt();
            theta2 = v2.Theta();
            //
        }
        if(pwflag2!=CHARGED_TRACK)continue;
        if(pt2<ptMin||pt2>ptMax) continue;
        if (N!=0)
          {
              //cout<<"hi"<<endl;
        h_2D->Fill(eta1-eta2,dphi(phi1,phi2),1./N);
        h_2D->Fill(eta1-eta2,dphi(phi2,phi1),1./N);
        h_2D->Fill(eta2-eta1,dphi(phi1,phi2),1./N);
        h_2D->Fill(eta2-eta1,dphi(phi2,phi1),1./N);
          }
      }
        }//end of second loop
    }
    }// end of first loop
  }// end of loop over events
  
    
    
    // BACKGROUND LOOP using no multiplicity cut
    for (Int_t i=0;i<nevent_process;i++) {
        
        if (i%10000==0) cout <<i<<"/"<<nevent_process<<endl;
        t1->GetEntry(i);
        
        int nparticles = nParticle;
        float ttheta = TTheta;
        float tphi = TPhi;
        if (verbose) cout<<"nparticles="<<nparticles<<endl;
        
        int nparticles2=-1000;
        float ttheta_mix;
        float tphi_mix;
        int selected=i+1;  //questo diventa il numero di evento estratto nel file +1
        int flag=0;    //definisco una flag a zero
        
        // Yen-Jie: cut on maximum number of particles to avoid infinite loop, for the moment it is 100
        if (nparticles>1000) continue;
        
        
        double N=0;
        double ptMin=0.1; //0.4
        double ptMax=4;
        for (int p = 0;p<num_runs;p++)
        {
            // calculate the number of tracks in the passing selection
            for ( int j=0;j<nparticles;j++ ) {
                float pt1 = pt[j];
                int pid1 = pid[j];
                int pwflag1 = pwflag[j];
                //if (pid1!=PION&&pid1!=PROTON&&pid1!=KAON&&!(!isBelle&&pwflag1==0)) continue;
                if (pwflag1!=CHARGED_TRACK)continue;
                if(pt1<ptMin||pt1>ptMax) continue;
                N++;
            }
            
            averageN+=N;
            nEventProcessed++;
            
            //if (nparticles<mult_low) continue;
            if (N<0) continue;
            
            nEventInMultBin++;
            
            // find a mixed event
            //	cout <<N<<endl;
            while ((fabs(nparticles2-nparticles)>5&&nparticles<1000)||i==selected){
                //cout <<nparticles<<" "<<selected<<endl;
                selected++;
                if (selected>nevent_process&&flag==1) break;
                if (selected>nevent_process) flag=1;
                selected = selected % nevent_process;
                t1_mix->GetEntry ( selected );
                nparticles2= nParticle_mix;
                ttheta_mix = TTheta_mix;
                tphi_mix =TPhi_mix;
            }
            
            double N2=0;
            // calculate the number of tracks in the mixed event passing selection
            for ( int j=0;j<nparticles;j++ ) {
                float pt1 = pt_mix[j];
                int pid1 = pid_mix[j];
                int pwflag1 = pwflag_mix[j];
                if (pwflag1!=CHARGED_TRACK) continue;
                //if (pid1!=PION&&pid1!=PROTON&&pid1!=KAON&&!(!isBelle&&pwflag1==0)) continue;
                if(pt1<ptMin||pt1>ptMax) continue;
                N2++;
            }
            
            TVector3 thrust(1, 1, 1);
            thrust.SetTheta(ttheta);
            thrust.SetPhi(tphi);
            
            
            TVector3 rotVec;
            rotVec.SetX(-thrust.Y());
            rotVec.SetY(thrust.X());
            
            TVector3 thrust_mix(1, 1, 1);
            thrust_mix.SetTheta(ttheta_mix);
            thrust_mix.SetPhi(tphi_mix);
            
            TVector3 rotVec_mix;
            rotVec_mix.SetX(-thrust_mix.Y());
            rotVec_mix.SetY(thrust_mix.X());
            
            for ( int j=0;j<nparticles;j++ ) {
                int pid1 = pid[j];
                float eta1 = eta[j];
                float phi1 = phi[j];
                float pt1 = pt[j];
                float mass1 = mass[j];
                float pwflag1 = pwflag[j];
                float theta1 = theta[j];
                if (isThrust == 1)
                {
                    TVector3 v1;
                    float px1 = px[j];
                    float py1 = py[j];
                    float pz1 = pz[j];
                    v1.SetXYZ(px1,py1,pz1);
                    v1.Rotate(-ttheta, rotVec);
                    phi1 = v1.Phi();
                    eta1 = v1.PseudoRapidity();
                    pt1 = v1.Pt();
                    theta1 = v1.Theta();
                    //
                }
                //cout<<phi1<<endl;
                if(pwflag1!=CHARGED_TRACK)continue;
                //if (pid1!=PION&&pid1!=PROTON&&pid1!=KAON&&!(!isBelle&&pwflag1==0)) continue;
                if(pt1<ptMin||pt1>ptMax) continue;
                
                // Background loop, calculate B correlation function from mixed event
                for ( int k=0;k<nparticles2;k++ ) {
                    int pidmix = pid_mix[k];
                    float etamix = eta_mix[k];
                    float phimix = phi_mix[k];
                    float ptmix = pt_mix[k];
                    float massmix = mass_mix[k];
                    float pwflagmix = pwflag_mix[k];
                    float thetamix = theta_mix[k];
                    if (isThrust == 1)
                    {
                        TVector3 vmix;
                        float pxmix = px_mix[j];
                        float pymix = py_mix[j];
                        float pzmix = pz_mix[j];
                        vmix.SetXYZ(pxmix,pymix,pzmix);
                        vmix.Rotate(-ttheta_mix, rotVec_mix);
                        phimix = vmix.Phi();
                        etamix = vmix.PseudoRapidity();
                        ptmix = vmix.Pt();
                        thetamix = vmix.Theta();
                        
                    }
                    if(pwflagmix!= CHARGED_TRACK)continue;
                    //if (pidmix!=PION&&pidmix!=PROTON&&pidmix!=KAON&&!(!isBelle&&pwflagmix==0)) continue;
                    if(ptmix<ptMin||ptmix>ptMax) continue;
                    if(N!=0)
                    {
                        h_2Dmix->Fill(eta1-etamix,dphi(phi1,phimix),1./N);    
                        h_2Dmix->Fill(eta1-etamix,dphi(phimix,phi1),1./N);    
                        h_2Dmix->Fill(etamix-eta1,dphi(phi1,phimix),1./N);    
                        h_2Dmix->Fill(etamix-eta1,dphi(phimix,phi1),1./N);
                        
                        //h_2Dmix_theta->Fill(eta1-etamix,dphi(phi1,phimix),1./N);
                        //h_2Dmix_theta->Fill(eta1-etamix,dphi(phimix,phi1),1./N);
                        //h_2Dmix_theta->Fill(etamix-eta1,dphi(phi1,phimix),1./N);
                        //h_2Dmix_theta->Fill(etamix-eta1,dphi(phimix,phi1),1./N);
                    }
                }//end of second loop
                
            }
        }// end of first loop
    }
    
  
  cout<<"nEventInMultBin"<<nEventInMultBin<<endl;
  
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


  cout<<"x axis "<<h_2Dmix->GetXaxis()->GetBinCenter(b00_x)<<endl;
  cout<<"y axis "<<h_2Dmix->GetYaxis()->GetBinCenter(b00_y)<<endl;

  TH2F *h_ratio_theta = new TH2F ( "h_ratio_theta", "theta-phi of all particles ", nbin, -PI*1.1, PI/6.,nbin, -3.1416/2.,3.1416*1.5);
    
  h_ratio_theta->Sumw2();
  /*
  for (int x=0;x<=h_2D->GetNbinsX();x++){
     for (int y=0;y<=h_2D->GetNbinsY();y++){
        if(h_2Dmix->GetBinContent(x,y)>0){
            
          ratio=B00*(h_2D->GetBinContent(x,y)/h_2Dmix->GetBinContent(x,y));
          errrel_num=h_2D->GetBinError(x,y)/h_2D->GetBinContent(x,y);
          errrel_den=h_2Dmix->GetBinError(x,y)/h_2Dmix->GetBinContent(x,y);
          errrel_ratio=TMath::Sqrt(errrel_num*errrel_num+errrel_den*errrel_den+errrel_B00*errrel_B00);
            //cout<<ratio*errrel_ratio<<endl;
          h_ratio->SetBinContent(x,y,ratio);
          h_ratio->SetBinError(x,y,ratio*errrel_ratio);
          
        }
     }
  
  }
  */
  double ratio_phi;
  double ratio_eta;
   
  //float deta_count = 2*detaRange/(float)nbin;
  //float dphi_count = (3.1416*1.5 - -3.1416/2.)/(float)nbin;
  
    for (int x=0;x<=h_2D->GetNbinsX();x++){

     for (int y=0;y<=h_2D->GetNbinsY();y++){

       
       if(h_2Dmix->GetBinContent(x,y)>0){
            
          ratio=B00*(h_2D->GetBinContent(x,y)/h_2Dmix->GetBinContent(x,y));
          errrel_num=h_2D->GetBinError(x,y)/h_2D->GetBinContent(x,y);
          errrel_den=h_2Dmix->GetBinError(x,y)/h_2Dmix->GetBinContent(x,y);
          errrel_ratio=TMath::Sqrt(errrel_num*errrel_num+errrel_den*errrel_den+errrel_B00*errrel_B00);
          
          ratio_eta = h_2D->GetXaxis()->GetBinCenter(x);
          ratio_phi = h_2D->GetYaxis()->GetBinCenter(y);
            //cout<<ratio*errrel_ratio<<endl;
          h_ratio->Fill(ratio_eta,ratio_phi,ratio);
          //h_ratio->SetBinError(x,y,ratio*errrel_ratio);

        }
     }

    }
        
  
  
  for(int binIterX = 0; binIterX < h_ratio->GetNbinsX(); ++binIterX){
      for(int binIterY = 0; binIterY < h_ratio->GetNbinsX(); ++binIterY){
          double eta = h_ratio->GetXaxis()->GetBinCenter(binIterX+1);
          double phi = h_ratio->GetYaxis()->GetBinCenter(binIterY+1);
          double val = h_ratio->GetBinContent(binIterX,binIterY);
          double err = h_ratio->GetBinError(binIterX,binIterY);
          double theta = -2*TMath::ATan(exp(-eta));
          h_ratio_theta->Fill(theta,phi,val);
      }
      
    }
    cout<<"bin"<<h_ratio_theta->GetBinContent(0,0)<<endl;
  double etaranges[8]={1.5,4,1.5,3, 2,3, 2,4};
  int minbin,maxbin;
    int thetaranges[4] = {-3,-2,-1,0};
  
  TH1D*h_deltaphi[7];
  TH1D*h_deltaphi_theta[3];

  //h_deltaphi  = (TH1D*) h_ratio->ProjectionY("h_deltaphi",0,-1);
  //h_deltaphi->Sumw2();
  
  for (int i=0;i<7;i =i+2){
    
    minbin =  h_ratio->GetXaxis()->FindBin(etaranges[i]);
    maxbin =  h_ratio->GetXaxis()->FindBin(etaranges[i+1]);
    cout<<"1"<<" "<<minbin<<" "<<maxbin<<endl;
    //h_deltaphi[i]  = (TH1D*) h_ratio->ProjectionY(Form("h_deltaphi_etamin%d_max%d",etaranges[i],etaranges[i+1]),minbin,maxbin);
    h_deltaphi[i]  = (TH1D*) h_ratio->ProjectionY(Form("h_deltaphi%d",i),minbin,maxbin);
    h_deltaphi[i]->Sumw2();
    /*if (i == 0) 
     {
       histogram  = (TH1D*) h_ratio->ProjectionY(Form("h_deltaphi_thetamin%d_max%d",thetaranges[i],thetaranges[i+1]),minbin,maxbin);
       histogram->Sumw2();
     }*/
    h_deltaphi[i]->SetName(Form("h_deltaphi_%d",i));
    //h_deltaphi[i]->GetZaxis()->CenterTitle();
    
    h_deltaphi[i]->GetXaxis()->SetTitle("#Delta#phi");
    h_deltaphi[i]->SetTitle(Form("#Delta#phi, #Delta#eta (%f, %f), Multipliplicity (%d, %d)",etaranges[i],etaranges[i+1], mult_low, mult_high));
    h_deltaphi[i]->GetYaxis()->SetTitle("Y(#Delta#phi)");
    cout<<1./(maxbin-minbin+1)<<endl;
    h_deltaphi[i]->Scale(1./(maxbin-minbin+1));
    /*
    TCanvas *c = new TCanvas("c","",600,600);
    c->SetTheta(60.839);
    c->SetPhi(38.0172);
    
    h_deltaphi[i]->Draw("hist");
    c->SaveAs(Form("likeATLAS/projection%d_%d_%d_t%d.pdf",mult_low,mult_high,i,isThrust));
    */
  }  
  

  for (int i=0;i<3;i++){
     minbin =  h_ratio_theta->GetXaxis()->FindBin(thetaranges[i]);
     maxbin =  h_ratio_theta->GetXaxis()->FindBin(thetaranges[i+1]);
     cout<<minbin<<" "<<maxbin<<endl;
     h_deltaphi_theta[i]  = (TH1D*) h_ratio_theta->ProjectionY(Form("h_deltaphi_thetamin%d_max%d",thetaranges[i],thetaranges[i+1]),minbin,maxbin);
     
     //h_deltaphi_theta[i]  = (TH1D*) h_ratio_theta->ProjectionY(Form("h_deltaphi_theta%d",i),minbin,maxbin);
     h_deltaphi_theta[i]->Sumw2();
     //h_deltaphi[i]->SetName(Form("h_deltaphi_etamin%d_max%d",etaranges[i],etaranges[i+1]));
     
     h_deltaphi_theta[i]->GetXaxis()->SetTitle("#Delta#phi For #theta");
      h_deltaphi_theta[i]->SetTitle(Form("deltaphi theta (%d, %d)",thetaranges[i],thetaranges[i+1]));
      
     h_deltaphi_theta[i]->Scale(1./(maxbin-minbin+1));
    }
    
  //cout<<"error"<<h_deltaphi[0]->GetBinError(5)<<endl;
  
  h_ratio->GetXaxis()->SetTitle("#Delta#eta");
  h_ratio->GetYaxis()->SetTitle("#Delta#phi");
  
  //TFile *fout;
  //if (isBelle) fout=new TFile(Form("/data/flowex/Datasamples/Belle/myoutput_minMult%d.root",mult_low),"recreate");
  //else fout=new TFile(Form("/data/flowex/Outputs/myoutput_minMult%d_modified.root",mult_low),"recreate");
  
  //fout = new TFile(Form("/Users/anthony/desktop/myoutput_minMult%d_modified1.root",mult_low),"recreate");
  
  /*
  fout = new TFile(Form("/mnt/c/Users/Bibek Kumar Pandit/Desktop/Root_Directory/StudyMult/VerificationCode/likeATLAS/myoutput_minMult%d_modified_t%d.root",mult_low, isThrust),"recreate");
  
  fout->cd();
  h_2D->Write();
  h_2Dmix->Write();
  h_ratio->Write();
  h_ratio_theta->Write();
  //h_deltaphi->Write();
  for (int i=0;i<3;i++)
    {//h_deltaphi[i]->Write();
        h_deltaphi_theta[i]->Write();
    }
    
  fout->Close();
  delete fout;

    /
    TCanvas *c1 = new TCanvas("c1","",600,600);
    c1->SetTheta(60.839);
    c1->SetPhi(38.0172);
    h_2D->SetTitle(Form("S correlation between %d and %d",mult_low,mult_high));
    h_2D->Draw("LEGO2");
    c1->SaveAs(Form("ridge_plots/S%d_%d_t%d.pdf",mult_low,mult_high,isThrust));
    
    TCanvas *c2 = new TCanvas("c2","",600,600);
    c2->SetTheta(60.839);
    c2->SetPhi(38.0172);
    h_2Dmix->SetTitle(Form("B correlation between %d and %d",mult_low,mult_high));
    h_2Dmix->Draw("LEGO2");
    c2->SaveAs(Form("ridge_plots/B%d_%d_t%d.pdf",mult_low,mult_high,isThrust));
    
    TCanvas *c3 = new TCanvas("c3","",600,600);
    c3->SetTheta(60.839);
    c3->SetPhi(38.0172);
    h_ratio->SetTitle(Form("Ratio of correlations between %d and %d",mult_low,mult_high));
    h_ratio->Draw("LEGO2");
    c3->SaveAs(Form("ridge_plots/ratio%d_%d_t%d.pdf",mult_low,mult_high,isThrust));
    

    TCanvas *c4 = new TCanvas("c4","",600,600);
    c4->SetTheta(60.839);
    c4->SetPhi(38.0172);
    h_ratio_theta->GetXaxis()->SetTitle("#Delta#theta");
    h_ratio_theta->GetYaxis()->SetTitle("#Delta#phi");
    h_ratio_theta->SetTitle(Form("Ratio of correlations (in #Delta#theta) between %d and %d",mult_low,mult_high));
    h_ratio_theta->Draw("LEGO2");
    c4->SaveAs(Form("ridge_plots/ratio_theta%d_%d_t%d.pdf",mult_low,mult_high,isThrust));
    
    histogram->FillRandom("gaus", 10000);
    
    TCanvas *c = new TCanvas("c","",600,600);
    c->cd();
    cout<<"hi"<<" "<<"a"<<endl;
    histogram->Draw();
  
    c->SaveAs("fit1.pdf");
    */
    TFile *background = new TFile(Form("correlation_%d_%d.root", mult_low, mult_high), "recreate");
    h_deltaphi[0]->Write();
    background->Close();
    
  }

TH1D *background;// = analysis(0, 0,  0, 30, 40, 50, 0, 1); 
// = analysis(0, 0,  0, 0, 10, 50, 0, 1);

Double_t ftotal(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = background->GetXaxis()->FindBin(xx);
   Int_t bin2 = background->GetXaxis()->FindBin(0.);
   Double_t br = par[1] + par[2]*background->GetBinContent(bin);
   //Double_t arg = (2*xx);
   Double_t sr = background->GetBinContent(bin2)+par[0]*TMath::Cos(2*xx);
   return br + sr;
}

void subtract()
{
  analysis();
  TFile *file = new TFile("correlation_0_10.root");
  background = (TH1D*)file->Get("h_deltaphi_0");
  //TH1D* hist2 = new TH1D("blablabla2", "",50,-3.1416/2., 3.1416*1.5);
  
  analysis(0, 0,  0, 30, 40, 50, 0, 1);
  
  TFile *file2 = new TFile("correlation_30_40.root");
  TH1D* result = (TH1D*)file2->Get("h_deltaphi_0");

  /*
  TCanvas *c1 = new TCanvas("c1","",600,600);
  c1->cd();
  cout<<"hi"<<" "<<"3"<<endl;
  hist1->Draw();
  
  c1->SaveAs("fit.pdf");
  cout<<"hi"<<" "<<"4"<<endl;
  */
  


   //fit function ftotal to signal + background

//   histgen();

   //TFile *f = new TFile("background.root");
   //background = (TH1F*)f->Get("background"); //pointer used in ftotal
   //TH1F *result = (TH1F*)f->Get("result");

  TF1 *ftot = new TF1("ftot",ftotal,0,10,4);
  Double_t norm = result->GetMaximum();
  ftot->SetParameters(0.5*norm,0.5*norm, 0.1);
  //ftot->SetParLimits(2,0.001,0.5);
    
  result->Fit("ftot", "b");
  
  TCanvas *c1 = new TCanvas("c1","",600,600);
  c1->cd();
  cout<<"hi"<<" "<<"3"<<endl;
  result->Draw();
  
  c1->SaveAs("fit.pdf");
  cout<<"hi"<<" "<<"4"<<endl;
  
  file->Close();
  file2->Close();
}