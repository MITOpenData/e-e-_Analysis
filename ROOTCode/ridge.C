#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
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

void analysis(int isBelle=1, int maxevt=100000,int mult=70,int nbin=20){

  TString filename;
  if(isBelle) filename="../Inputs/output-2.root";
  else filename="../LEP2dataMarcello/myALEPH.root";
  
  TFile *f = new TFile(filename.Data());
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

  TFile *f_mix = new TFile(filename.Data());
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
  h_2D->Sumw2();
  h_2Dmix->Sumw2();
  h_ratio->Sumw2();
  

  // all entries and fill the histograms
  Int_t nevent = (Int_t)t1->GetEntries();
  
  int nevent_process = nevent;
  if( maxevt>0 && maxevt<nevent ) nevent_process = maxevt;  
  
  for (Int_t i=0;i<nevent_process;i++) {
    if (i%100==0) cout <<i<<"/"<<nevent_process<<endl;
    t1->GetEntry(i);
    
    int nparticles = nParticle;
//    cout<<"nparticles="<<nparticles<<endl;
    if (nparticles<mult) continue;
  
    int nparticles2=-1000;
    int selected=i+1;  //questo diventa il numero di evento estratto nel file +1
    int flag=0;    //definisco una flag a zero
    if (nparticles>100) continue;
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
     // if (pid1!=PION&&pid1!=PROTON&&pid1!=KAON) continue;
      float eta1 = eta[j];
      float phi1 = phi[j];
      float pt1 = pt[j];
      float mass1 = mass[j];
      if(pt1<0||pt1>4) continue;
      
      //cout<<"particle1, pid="<<pid1<<", eta="<<eta1<<", phi="<<phi1<<", eta="<<eta1<<", pt="<<pt1<<", mass="<<mass1<<endl;;
      
      for ( int k=j+1;k<nparticles;k++ ) {
        int pid2 = pid[k];
        //if (pid2!=PION&&pid2!=PROTON&&pid2!=KAON) continue;
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
        //if (pidmix!=PION&&pidmix!=PROTON&&pidmix!=KAON) continue;
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

    int etaranges[4]={0,1,2,3};
    int minbin,maxbin;
    
	TH1D*h_deltaphi[3];

	for (int i=0;i<3;i++){
	  
      minbin =  h_ratio->GetXaxis()->FindBin(etaranges[i]);
      maxbin =  h_ratio->GetXaxis()->FindBin(etaranges[i+1]);
      
        h_deltaphi[i]  = (TH1D*) h_ratio->ProjectionY(Form("h_deltaphi_etamin%d_max%d",etaranges[i],etaranges[i+1]),minbin,maxbin);
        h_deltaphi[i]->Sumw2();
		h_deltaphi[i]->SetName(Form("h_deltaphi_etamin%d_max%d",etaranges[i],etaranges[i+1]));
	    h_deltaphi[i]->GetZaxis()->CenterTitle();
	    h_deltaphi[i]->GetXaxis()->SetTitle("#Delta#phi");
	}  


  cout<<"error"<<h_deltaphi[0]->GetBinError(5)<<endl;

  h_ratio->GetXaxis()->SetTitle("#Delta#eta");
  h_ratio->GetYaxis()->SetTitle("#Delta#phi");
  
  TCanvas * c1 = new  TCanvas("c1","c1",1000,500); 
  h_ratio->Draw("lego2");

  TCanvas * c2 = new TCanvas("c2");
  c2->Divide(3,1);
  c2->cd(1);
  h_deltaphi[0]->Draw("ep");
  c2->cd(2);
  h_deltaphi[1]->Draw("ep");
  c2->cd(3);
  h_deltaphi[2]->Draw("ep");


  TFile*fout=new TFile(Form("myoutput_isBelle%d_minMult%d.root",isBelle,mult),"recreate");
  fout->cd();
  h_2D->Write();
  h_2Dmix->Write();
  h_ratio->Write();
  for (int i=0;i<3;i++) h_deltaphi[i]->Write();
  fout->Close();
  delete fout;

  }
