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
    
    TCanvas *c1 = new TCanvas("mult", "mult", 600,  600);
    TH1F *h1 = new TH1F("mult", "mult", 100, 0 , 100);
    
    // all entries and fill the histograms
    Int_t nevent = (Int_t)t1->GetEntries();

    if (maxevt > 0 && maxevt < nevent) nevent_process = maxevt;
    int nevent_process = nevent;
    
    for (Int_t i = 0; i<nevent_process; i++)
    {
      t1->GetEntry(i);
      int nparticles = nParticle;
      
      if (i % 100 == 0){
      cout<<i<<endl;
      }
      
      int N =0;
      
      for (Int_t j = 0; j<nparticles; j++)
      {
        //pwflag1 = pwflag[j];
        //cout<<pwflag1<<endl;
        if (pwflag[j] == CHARGED_TRACK)
        {
          //cout<<"hi"<<endl;
          N++;
        }
        
      }
      //if (N == nparticles)
      //{
        //cout<<"hi"<<endl;
      //}
      //cout<<nparticles - N<<endl;
      for (int k = 0; k<=N; k++)
      {
        h1->Fill(k);
      }
    }
    c1->cd();
    c1->SetLogy();
    h1->Draw();
    
    //c1->SaveAs(Form("Images/eta_photon.pdf"));
}

