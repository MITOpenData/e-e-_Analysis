// The numbers of particles in each event are stored in a histogram.
// All particle masses in the analyzed events are calculated and stored in a histogram.
// All pion masses in the analyzed events are calculated and stored in a histogram.
// The neutral pion (pizero) and the charged pion are differentiated and respective
// mass distributions are produced.
//+
// File : pion.cc
// Description : Analyze a data file
// 
// Author : Ryosuke Itoh, IPNS, KEK
// Date : 1 - Mar - 2004
//-
// Modified by T.Nozaki  (2005/09/21)
// Modified by S.Nishida (2007/02/03)
// Modified by T.Nozaki (2009/07/01) Comments in English

#include <cstdio>
#include <cmath>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TClonesArray.h"
#include "TBrowser.h"
#include "TLorentzVector.h"
#include "BParticle.cc"
#include "BEvent.h"
#include <iostream>
//enum SIMPLEPID {PHOTON, ELECTRON, PION, MUON, KAON, PROTON};

using namespace std;

#define PI 3.1415926

double dphi(double phi1,double phi2)
{
    double a=phi1-phi2;
    while (a<-PI) a+=2*PI;
    while (a>PI) a-=2*PI;
    
    if (a<-PI/2) a=2*PI+a;
    return a;
}

void analysis ( char* file, int mult=30,int nbin=15, int maxevt = 0 ) {
   /* These are initialization codes.
     You can ignore them as a black box. */

  // Open a data file
  TFile *f = new TFile( file );
  TFile *f2 = new TFile( file );

  // Obtain a pointer to a series of "event" data in the file
  TTree *t = (TTree*)f->Get("T");
  TTree *t2 = (TTree*)f2->Get("T");

  // Create a pointer to "BEvent" object where data are loaded from the file
  BEvent* event = new BEvent();
  BEvent* event2 = new BEvent();

  // Obtain a branch for "BEvent" in the tree
  TBranch *branch = t->GetBranch("BEvent");
  TBranch *branch2 = t2->GetBranch("BEvent");

  // register created "BEvent" object to the branch to load event data
  branch->SetAddress(&event);
  branch2->SetAddress(&event2);

   /*  Initialization codes finish here.*/

  // ---------------------------------------------------------------------- //

   /* We define the  histograms.*/

  // First, we define the name of the file where histograms will be saved. 
  // We define "histfile" as the histogram file name, as shown below. 
  TFile *hf = new TFile("histfile", "RECREATE" );

 // We define each histogram.  
  // The first one is the distribution of the number of particles in each event.
  // The longitudinal axis is divided into 50 bins between 0 and 50. 
  // The first argument in the blacket is the name of the histogram and 
  // the second one is its explanation.

  TH1F *h_nprt = new TH1F ( "h_nprt", "No. of particles", 300, 0, 300 );
  
  // The distribution of all particle masses.
  TH1F *h_mass = new TH1F ( "h_mass", "mass of all particles ",
			    1000, 0.0, 1.0);
  // The distribution of all pion masses.
  TH1F *h_pimass = new TH1F ( "h_pimass", "mass of pion ",
			      1000, 0.1, 0.2);
  // The distribution of all charged pion  masses.
  TH1F *h_pipmass = new TH1F ( "h_pipmass", "mass of charged pion ",
			       1000, 0.1, 0.2 );
  // The distribution of all neutral pion  masses.
  TH1F *h_pi0mass = new TH1F ( "h_pi0mass", "mass of neutral pion ",
			       1000, 0.1, 0.2 );
//  TH1F *h_mass = new TH1F ( "h_mass", "mass of all particles ",
//			    1000, 0.0, 1.0);
  TH1F *h_eta = new TH1F ( "h_eta", "eta of all particles ",
			    1000, -4.0, 4.0);

  TH2F *h_2D = new TH2F ( "h_2D", "eta-phi of all particles ",
			    nbin, -3.5, 3.5,nbin, -3.1416/2., 3.1416*1.5);
  TH2F *h_mix = new TH2F ( "h_mix", "eta-phi of all particles ",
			    nbin, -3.5, 3.5,nbin, -3.1416/2.,3.1416*1.5);
  TH2F *h_ratio;
  // TH2F *h_ratio = new TH2F ( "h_ratio", "eta-phi of all particles ",
			    // nbin, -4.5, 4.5,nbin, -3.1416/2.,3.1416*1.5);

  // ---------------------------------------------------------------------- //

   /* The analysis codes start here. */

  // Let us examine how many events exists in the data file.
  int nevent = (int)t->GetEntries();

  // Let us determine how many events should be analyzed.
  // If maxevent is not specified all the events in the data are analyzed.
  int nevent_process = nevent;
  if( maxevt>0 && maxevt<nevent ) nevent_process = maxevt;

  // Analysis is carried out event by event.
  // loop of events.
  int multiplicitycut = 50;
  for ( int i=0;i<nevent_process;i++ ) {
    //  Let us write an event number every 10000 events.
    // It helps for debugging when the program is dead.
    if ( i%10000 == 0 ) printf ( "***** Event %d\n", i );

    // An event is read.
    int nb = t->GetEntry ( i );
//    int nb2 = t2->GetEntry ( i-1 );

    // If the comment delimeter in the following sentence is removed, 
    // the detailed information of the first event is shown.

    /* -------- The main body of the analysis codes--------------- */

    // Obtain the number of particles in each event.
    int nparticles = event->NParticles();
//    cout <<nparticles<<endl;
    if (nparticles<mult) continue;

    int nparticles2=-1000;
    int selected=i+1;
    int flag=0;
    if (nparticles<multiplicitycut) continue;
    while (abs(nparticles2-nparticles)>5&&nparticles<100||i==selected)
    {
     //cout <<nparticles<<" "<<selected<<endl;
     selected++;
     if (selected>nevent_process&&flag==1) break;
     if (selected>nevent_process) flag=1;
     selected = selected % nevent_process;
     t2->GetEntry ( selected );
     nparticles2= event2->NParticles();
    }
    
     // It is stored in the histogram named h_nprt.
    h_nprt->Fill ( nparticles );

    // Particles are treated as an array.
    TClonesArray& plist = *(event->GetParticleList());
    TClonesArray& plist2 = *(event2->GetParticleList());
    TLorentzVector v;
    TLorentzVector v2;
    // loop of all particles.
    for ( int j=0;j<nparticles;j++ ) {
      // The j-th particle is named  p1.
      BParticle* p1 = (BParticle*)plist[j];
      int pid1 = p1->pid();
      if (pid1!=PION&&pid1!=PROTON&&pid1!=KAON) continue;
      // Obtain the momentum and energy of a particle.
      float px1 = p1->px();  // x component of momentum.
      float py1 = p1->py();  // y component of momentum.
      float pz1 = p1->pz();  // z component of momentum.
      float  e1 = p1->e();   // Energy.
      float mass1 = sqrt ( e1*e1-px1*px1-py1*py1-pz1*pz1 ); // Calculate the mass.
      v.SetXYZM(px1,py1,pz1,mass1);
      float eta1 = v.PseudoRapidity();
      float phi1 = v.Phi();
      if (v.Pt()<0||v.Pt()>4) continue;
      
      for ( int k=j+1;k<nparticles;k++ ) {
        BParticle* p2 = (BParticle*)plist[k];
        int pid2 = p2->pid();
        if (pid2!=PION&&pid2!=PROTON&&pid2!=KAON) continue;
        // Obtain the momentum and energy of a particle.
        float px2 = p2->px();  // x component of momentum.
        float py2 = p2->py();  // y component of momentum.
        float pz2 = p2->pz();  // z component of momentum.
        float  e2 = p2->e();   // Energy.
        float mass2 = sqrt ( e2*e2-px2*px2-py2*py2-pz2*pz2 ); // Calculate the mass.
        v2.SetXYZM(px2,py2,pz2,mass2);
        float eta2 = v2.PseudoRapidity();
        float phi2 = v2.Phi();
        if (v2.Pt()<0||v2.Pt()>4) continue;
        h_2D->Fill(eta1-eta2,dphi(phi1,phi2));    
        h_2D->Fill(eta1-eta2,dphi(phi2,phi1));    
        h_2D->Fill(eta2-eta1,dphi(phi1,phi2));    
        h_2D->Fill(eta2-eta1,dphi(phi2,phi1));    
      } 

      for ( int k=0;k<nparticles2;k++ ) {
        BParticle* p2 = (BParticle*)plist2[k];
        int pid2 = p2->pid();
        if (pid2!=PION&&pid2!=PROTON&&pid2!=KAON) continue;
        // Obtain the momentum and energy of a particle.
        float px2 = p2->px();  // x component of momentum.
        float py2 = p2->py();  // y component of momentum.
        float pz2 = p2->pz();  // z component of momentum.
        float  e2 = p2->e();   // Energy.
        float mass2 = sqrt ( e2*e2-px2*px2-py2*py2-pz2*pz2 ); // Calculate the mass.
        v2.SetXYZM(px2,py2,pz2,mass2);
        if (v.Pt()<0||v.Pt()>4) continue;
        float eta2 = v2.PseudoRapidity();
        float phi2 = v2.Phi();
        h_mix->Fill(eta1-eta2,dphi(phi1,phi2));    
        h_mix->Fill(eta1-eta2,dphi(phi2,phi1));    
        h_mix->Fill(eta2-eta1,dphi(phi1,phi2));    
        h_mix->Fill(eta2-eta1,dphi(phi2,phi1));    
      } 

    }
    // Clear at the end of particle loop.
    event->Clear();
    event2->Clear();
  }
  h_2D->Sumw2();
  h_ratio = (TH2F*) h_2D->Clone("h_ratio");
  h_ratio->Sumw2();
  h_mix->Sumw2();

  h_ratio->Divide(h_mix);
  // for (int x=0;x<=h_2D->GetNbinsX();x++){
  // for (int y=0;y<=h_2D->GetNbinsY();y++){
     // if (h_mix->GetBinContent(x,y)>0) {h_ratio->SetBinContent(x,y,h_2D->GetBinContent(x,y)/h_mix->GetBinContent(x,y));}
  // }
  // }
  
  
  TFile*foutput=new TFile("foutput.root'","RECREATE");
  foutput->cd();
  h_2D->Write();
  h_ratio->Write();
  h_mix->Write();
  foutput->Close();
      
  // Let us show how many events are processed.
  printf( "***** %d events processed. Exit.\n", nevent_process );
  // Let us store  all the histogram data in the file named hist
  hf->Write();
  int minbin =  h_ratio->GetXaxis()->FindBin(2);
  int maxbin =  h_ratio->GetXaxis()->FindBin(4);

  TH1D * h_deltaphi = (TH1D*) h_ratio->ProjectionY("h_deltaphi",minbin, maxbin);
  
  h_ratio->GetXaxis()->SetTitle("#Delta#eta");
  h_ratio->GetYaxis()->SetTitle("#Delta#phi");
  h_deltaphi->GetZaxis()->CenterTitle();
  
  h_deltaphi->GetXaxis()->SetTitle("#Delta#phi");
  h_deltaphi->GetYaxis()->CenterTitle();
  h_deltaphi->Draw();
  
  TCanvas * c1 = new  TCanvas("c1"); 
  h_ratio->SetTitle(Form("eta-phi of all particles(mult_cut: %d)", multiplicitycut));
  h_ratio->Draw("lego2");
  //c1->SaveAs(Form("eta-phi_%d.pdf",multiplicitycut));

  TCanvas * c2 = new TCanvas("c2");
  h_deltaphi->SetTitle(Form("projection of #Delta#eta on #Delta#phi(mult_cut: %d)", multiplicitycut));
  h_deltaphi->Draw("pe");
  //c2->SaveAs(Form("delta-phi_%d.pdf",multiplicitycut));
  //h_ratio->Draw("lego2");
  // Display histograms.
//  showhist();
}   
// The function showhist is used when the histograms are displayed.
void showhist() {
  TFile *f = new TFile ( "histfile" );
  TBrowser *t = new TBrowser;
}
