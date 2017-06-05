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
   
}

