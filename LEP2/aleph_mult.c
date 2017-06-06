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
    
    TCanvas *c1 = new TCanvas("h", "h", 600,  600);
    TCanvas *c2 = new TCanvas("i", "i", 600,  600);
    TCanvas *c3 = new TCanvas("j", "j", 600,  600);
    TH1F *h1 = new TH1F("photon","h",100,0,100);
    TH1F *h2 = new TH1F("charged","h",100,0,100);
    TH1F *h3 = new TH1F("neutral","h",100,0,100);
    
    // all entries and fill the histograms
    Int_t nevent = (Int_t)t1->GetEntries();
    
    if (maxevt > 0 && maxevt < nevent) nevent_process = maxevt;
    int nevent_process = nevent;
    
    for (Int_t i = 0; i<nevent_process; i++)
    {
        t1->GetEntry(i);
        int nparticles = nParticle;
        int P =0;
        int C = 0;
        int N = 0;
        for (Int_t j =0;j<nparticles;j++)
        {
            if(pwflag[j]==PHOTON)
            {
                P+=1;
            }
            
            if(pwflag[j] == CHARGED_TRACK || pwflag[j] == CHARGED_LEPTONS1 || pwflag[j] == CHARGED_LEPTONS2)
            {
                C+=1;
            }
            
            if(pwflag[j] == NEUTRAL_HADRON) //pwflag[j] == V0 ||
            {
                N+=1;
            }
            
        }
        for (int j = 0; j<=P;j++){h1->Fill(j);}
        for (int k = 0; k<=C;k++){h2->Fill(k);}
        for (int l = 0; l<=N;l++){h3->Fill(l);}
    }
    c1->cd();
    c1->SetLogy();
    h1->Draw();
    c2->cd();
    c2->SetLogy();
    h2->Draw();
    c3->cd();
    c3->SetLogy();
    h3->Draw();
    

}

