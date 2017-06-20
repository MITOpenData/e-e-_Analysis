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

void analysis(){

  TString filename1;
  filename1="/data/flowex/Datasamples/Belle/myoutput_isBelle1_minMult40.root";
  
  TString filename2;
  filename2="/data/flowex/Datasamples/Belle/myoutput_isBelle1_minMult0.root";
  
  TFile *f1 = new TFile(filename1.Data());
  //TTree *t1 = (TTree*)f1->Get("t");
  
  TFile *f2 = new TFile(filename2.Data());
  //TTree *t2 = (TTree*)f2->Get("t");
  
  TH2F * h1 = (TH2F*)f1->Get("h_ratio");
  TH2F * h2 = (TH2F*)f2->Get("h_ratio");
  
  
  TCanvas *c = new TCanvas("mult", "mult", 600,  600);
  
  
  // Define the ratio plot
  TH2F *h3 = (TH2F*)h1->Clone("h3");
  //h3->SetLineColor(kBlack);
  //h3->SetMinimum(0.8);  // Define Y ..
  //h3->SetMaximum(1.35); // .. range
  h3->Sumw2();
  h3->SetStats(0);      // No statistics on lower plot
  h3->Divide(h2);
  //h3->SetMarkerStyle(21);
  
  c->cd();
  h3->Draw("ep");       // Draw the ratio plot


}
