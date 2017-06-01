//
//  aleph_analysis.c
//  
//
//  Created by Anthony Badea and Bibek K Pandit on 5/31/17.
//
//

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
#define plank (4.135667516 * pow(10,-18))
enum SIMPLEPWFLAG {CHARGED_TRACK, CHARGED_LEPTONS1, CHARGED_LEPTONS2, V0, PHOTON, NEUTRAL_HADRON};
// LEP2mcggbby200e


void analysis(int isBelle=1, int maxevt=0,int mult=50,int nbin=40,bool verbose=0){
    
    TString filename;
    //if(isBelle) filename="/Users/anthony/Documents/StudyMult/LEP2/ROOTfiles/cleaned_ALEPH_DATA-all.aleph.root";
    if(isBelle) filename="ROOTfiles/cleaned_ALEPH_DATA-all.aleph.root";
    
    TFile *f = new TFile(filename.Data());
    TTree *t1 = (TTree*)f->Get("t");
    Int_t nParticle;
    Float_t pt[50000];
    Float_t eta[50000];
    Float_t pid[50000];
    Float_t phi[50000];
    Float_t mass[50000];
    Float_t pwflag[50000];
    Float_t Energy[50000];
    
    t1->SetBranchAddress("nParticle",&nParticle);
    t1->SetBranchAddress("pt",pt);
    t1->SetBranchAddress("eta",eta);
    t1->SetBranchAddress("pid",pid);
    t1->SetBranchAddress("phi",phi);
    t1->SetBranchAddress("mass",mass);
    t1->SetBranchAddress("pwflag",pwflag);
    t1->SetBranchAddress("Energy",Energy);
    
    TCanvas *c1 = new TCanvas("eta_photon", "eta_photon", 600,  600);
    TH1F *h1 = new TH1F("eta_photon", "eta_photon", 100, -PI,PI);
    TCanvas *c2 = new TCanvas("eta_charged", "eta_charged", 600,  600);
    TH1F *h2 = new TH1F("eta_charged", "eta_charged", 100, -PI,PI);
    TCanvas *c3 = new TCanvas("eta_neutral", "eta_neutral", 600,  600);
    TH1F *h3 = new TH1F("eta_neutral", "eta_neutral", 100, -PI,PI);
    
    TCanvas *c4 = new TCanvas("phi_photon", "phi_photon", 600,  600);
    TH1F *h4 = new TH1F("phi_photon", "phi_photon", 100, -PI,PI);
    TCanvas *c5 = new TCanvas("phi_charged", "phi_charged", 600,  600);
    TH1F *h5 = new TH1F("phi_charged", "phi_charged", 100, -PI,PI);
    TCanvas *c6 = new TCanvas("phi_neutral", "phi_neutral", 600,  600);
    TH1F *h6 = new TH1F("phi_neutral", "phi_neutral", 100, -PI,PI);
    
    TCanvas *c7 = new TCanvas("2d_photon", "2d_photon", 600,  600);
    TH2F *h7 = new TH2F("2d_photon", "2d_photon", 100, -PI,PI,100, -PI, PI);
    TCanvas *c8 = new TCanvas("2d_charged", "2d_charged", 600,  600);
    TH2F *h8 = new TH2F("2d_charged", "2d_charged", 100, -PI,PI,100, -PI, PI);
    TCanvas *c9 = new TCanvas("2d_neutral", "2d_neutral", 600,  600);
    TH2F *h9 = new TH2F("2d_neutral", "2d_neutral", 100, -PI,PI,100, -PI,PI);
    
    char *h_photon_energy = new char[100];
    char *h_charged_energy = new char[100];
    char *h_neutral_energy = new char[100];
    char *h_all_energy = new char[100];
    
    TCanvas *c10 = new TCanvas("mult", "mult", 600,  600);
    TH1F *h10 = new TH1F("mult", "mult", 100, 0, 100);
    
    
    
    // all entries and fill the histograms
    Int_t nevent = (Int_t)t1->GetEntries();

    if (maxevt > 0 && maxevt < nevent) nevent_process = maxevt;
    int nevent_process = nevent;
    
    for (int e = 130; e<=210; i = i + 10)
    {
      sprintf(h_photon_energy, "h_photon_%d", i);
      TH1F *h_photon = new TH1F (h_photon, h_photon, 100, 0, 100);
      
    }
    
    for (Int_t i = 0; i<nevent_process; i++)
    {
        t1->GetEntry(i);
        int nparticles = nParticle;
        for (Int_t j =0;j<nparticles;j++)
        {
            if(pwflag[j]==PHOTON)
            {
                h1->Fill(eta[j]);
                h4->Fill(phi[j]);
                h7->Fill(eta[j],phi[j]);
            }
            
            if(pwflag[j] == CHARGED_TRACK || pwflag[j] == CHARGED_LEPTONS1 || pwflag[j] == CHARGED_LEPTONS2)
            {
                h2->Fill(eta[j]);
                h5->Fill(phi[j]);
                h8->Fill(eta[j],phi[j]);
            }
            
            if(pwflag[j] == NEUTRAL_HADRON) //pwflag[j] == V0 ||
            {
                h3->Fill(eta[j]);
                h6->Fill(phi[j]);
                h9->Fill(eta[j],phi[j]);
            }
        }
        
    }
    /*
    c1->cd();
    h1->Draw();
    c2->cd();
    h2->Draw();
    c3->cd();
    h3->Draw();
    c4->cd();
    h4->Draw();
    c5->cd();
    h5->Draw();
    c6->cd();
    h6->Draw();
    c7->cd();
    h7->Draw();
    c8->cd();
    h8->Draw();
    c9->cd();
    h9->Draw();
    */
    c10->cd();
    t1->Draw("Sum$(pwflag == CHARGED_TRACK):Energy>>h10","","prof");
    //h10->D
}

