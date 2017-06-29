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


void analysis(int isBelle, int maxevt,int mult, int mult_upper_bound, int nbin,bool verbose,int isThrust){
  TString filename;
  //if(isBelle) filename="/data/flowex/Datasamples/Belle/output_2_withtheta.root";
  //else filename="/data/flowex/Datasamples/LEP2_MAIN/ROOTfiles/cleaned_ALEPH_Data-all.aleph.root";
  
  filename = "/mnt/c/Users/Bibek Kumar Pandit/Desktop/Root_Directory/StudyMult/LEP2/ROOTfiles/cleaned_ALEPH_Data2-all.aleph.root";
  
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
  Float_t TTheta[100000];
  Float_t TPhi[100000];
  
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
  t1->SetBranchAddress("TTheta", TTheta);
  t1->SetBranchAddress("TPhi",TPhi);

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
  
  // two histograms
  double detaRange = 3.0;
  double normalization = detaRange*2/nbin*2*3.14159/nbin;
  TH2F *h_2D = new TH2F ( "h_2D", "eta-phi of all particles ",nbin*1.5, -detaRange, detaRange,nbin, -3.1416/4., 3.1416*0.75);
  TH2F *h_2Dmix = new TH2F ( "h_2Dmix", "eta-phi of all particles ",nbin*1.5, -detaRange, detaRange,nbin, -3.1416/4., 3.1416*0.75);
  TH2F *h_ratio = new TH2F ( "h_ratio", "eta-phi of all particles ", nbin*1.5, -detaRange, detaRange,nbin, -3.1416/4.,3.1416*0.75);
  
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
    double num_runs = 10;
    
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
        
        
    
    // calculate the number of tracks in the passing selection
    double N=0;

    for (int p = 0; p < num_runs; p++)
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

      //if (nparticles<mult) continue;
      if (N<mult) continue;
      
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
      
        
      for ( int j=0;j<nparticles;j++ )
      {
        float mass1 = mass[j];
        int pid1 = pid[j];
        float pwflag1 = pwflag[j];
        
        float eta1;
        float phi1;
        float pt1;
        // If isThrust == 0 then do not rotate by thrust angle
        // If isThrust == 1 then rotate by thrust angle
        if (isThrust == 0)
        {
            eta1 = eta[j];
            phi1 = phi[j];
            pt1 = pt[j];
            //if(pt1<ptMin||pt1>ptMax) continue;
        }
          
        else //(isThrust == 1)
        {
            TLorentzVector v1;
            float px1 = px[j];
            float py1 = py[j];
            float pz1 = pz[j];
            float thrust_angle1 = TTheta[i];//EDIT
            v1.SetXYZM(px1,py1,pz1,mass1);
            v1.RotateZ(thrust_angle1);
            phi1 = v1.Phi();
            eta1 = v1.PseudoRapidity();
            
            pt1 = v1.Pt();
            //if(pt1<ptMin||pt1>ptMax) continue;
        }

        if (pwflag1 != CHARGED_TRACK) continue;
        if(pt1<ptMin||pt1>ptMax) continue;
        
        // Signal loop, calculate S correlation function
        if (p == 0)
        {
          for ( int k=j+1;k<nparticles;k++ )
          {
            int pid2 = pid[k];
            float mass2 = mass[k];
            float pwflag2 = pwflag[k];
            
            float eta2;
            float phi2;
            float pt2;
            
            if( isThrust == 0)
            {
                eta2 = eta[k];
                phi2 = phi[k];
                pt2 = pt[k];
            }
              
            if( isThrust == 1)
            {
                TLorentzVector v2;
                float px2 = px[k];
                float py2 = py[k];
                float pz2 = pz[k];
                float thrust_angle2 = TTheta[i];
                v2.SetXYZM(px2,py2,pz2,mass2);
                v2.RotateZ(thrust_angle2);
                eta2 = v2.PseudoRapidity();
                phi2 = v2.Phi();
                pt2 = v2.Pt();
            }
            
            
            if (pwflag2 != CHARGED_TRACK) continue;
            if(pt2<ptMin||pt2>ptMax) continue;
            
            h_2D->Fill(eta1-eta2,dphi(phi1,phi2),1./N);    
            h_2D->Fill(eta1-eta2,dphi(phi2,phi1),1./N);    
            h_2D->Fill(eta2-eta1,dphi(phi1,phi2),1./N);    
            h_2D->Fill(eta2-eta1,dphi(phi2,phi1),1./N);
              
           /* TEST IDEA we want to be able to see two peaks one in the front and one in the back. try adding an if statement if eta is greater than 180 then take absolute value
          h_2D->Fill(fabs(eta1)-fabs(eta2),dphi(phi1,phi2),1./N);
          h_2D->Fill(fabs(eta1)-fabs(eta2),dphi(phi2,phi1),1./N);
          h_2D->Fill(fabs(eta2)-fabs(eta1),dphi(phi1,phi2),1./N);
          h_2D->Fill(fabs(eta2)-fabs(eta1),dphi(phi2,phi1),1./N);
            */
          }//end of second loop 
        }

        // Background loop, calculate B correlation function from mixed event
        for ( int k=0;k<nparticles2;k++ )
        {
          int pidmix = pid_mix[k];
          float massmix = mass_mix[k];
          float pwflagmix = pwflag_mix[k];
          
          float etamix;
          float phimix;
          float ptmix;
          
          if (isThrust == 0)
          {
              etamix = eta_mix[k];
              phimix = phi_mix[k];
              ptmix = pt_mix[k];
          }
          
          if (isThrust == 1)
          {
              TLorentzVector vmix;
              float pxmix = px_mix[k];
              float pymix = py_mix[k];
              float pzmix = pz_mix[k];
              float thrust_angle_mix = TTheta_mix[selected];
              vmix.SetXYZM(pxmix,pymix,pzmix,massmix);
              vmix.RotateZ(thrust_angle_mix);
              etamix = vmix.PseudoRapidity();
              phimix = vmix.Phi();
              ptmix = vmix.Pt();
          }
          
          
          if (pwflagmix != CHARGED_TRACK) continue;
          if(ptmix<ptMin||ptmix>ptMax) continue;

          h_2Dmix->Fill(eta1-etamix,dphi(phi1,phimix),1./N);
          h_2Dmix->Fill(eta1-etamix,dphi(phimix,phi1),1./N);
          h_2Dmix->Fill(etamix-eta1,dphi(phi1,phimix),1./N);
          h_2Dmix->Fill(etamix-eta1,dphi(phimix,phi1),1./N);
        }//end of second loop 

    }// end of first loop
    }
    // NOW DIVIDE
    //h_2Dmix->Scale(1./N);
  }//end of second loop
    
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


  cout<<"x axis"<<h_2Dmix->GetXaxis()->GetBinCenter(b00_x);
  cout<<"y axis"<<h_2Dmix->GetYaxis()->GetBinCenter(b00_y);
  
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
    
  h_ratio->Scale(num_runs);
    
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
  
  TFile *fout;
  //if (isBelle) fout= new TFile(Form("/data/flowex/Outputs/myoutput_minMult%d.root",mult),"recreate");
  //else fout=new TFile(Form("/data/flowex/Outputs/myoutput_minMult%d.root",mult),"recreate");
  fout = new TFile(Form("/mnt/c/Users/Bibek Kumar Pandit/Desktop/Root_Directory/StudyMult/DataFiles/ROOTfiles_from_ROOTCode/myoutput_minMult%d.root", mult), "recreate");
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
  h_2D->SetTitle(Form("S correlation between %d and %d",mult,mult_upper_bound));
  h_2D->Draw("LEGO2");
  c1->SaveAs(Form("S%d_%d.pdf",mult,mult_upper_bound));
  
  TCanvas *c2 = new TCanvas("c2","",600,600);
    c2->SetTheta(60.839);
    c2->SetPhi(38.0172);
  h_2Dmix->SetTitle(Form("B correlation between %d and %d",mult,mult_upper_bound));
  h_2Dmix->Draw("LEGO2");
  c2->SaveAs(Form("B%d_%d.pdf",mult,mult_upper_bound));
    
  TCanvas *c3 = new TCanvas("c3","",600,600);
    c3->SetTheta(60.839);
    c3->SetPhi(38.0172);
  h_ratio->SetTitle(Form("Ratio of correlations between %d and %d",mult,mult_upper_bound));
  h_ratio->Draw("LEGO2");
  c3->SaveAs(Form("ratio%d_%d.pdf",mult,mult_upper_bound));

  }

void goofy()
{
  int low[1] = {0};
  int high[1] = {1000};
  for (int i = 0; i<1; i++){
    analysis(0, 0, low[i], high[i], 40, 0,0);
  
//void analysis(int isBelle, int maxevt,int mult, int mult_upper_bound, int nbin,bool verbose,int isThrust)
  }
}

