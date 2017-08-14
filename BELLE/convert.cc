// The numbers of particles in each event are stored in a histogram.
// All particle masses in the analyzed events are calculated and stored in a histogram.
// All pion masses in the analyzed events are calculated and stored in a histogram.
// The neutral pion (pizero) and the charged pion are differentiated and respective
// mass distributions are produced.
//+
// File : convert.cc
// Description : convert BELLE open data format to TPCNtuple
// 
// Author : Ryosuke Itoh, IPNS, KEK
// Date : 1 - Mar - 2004
//-
// Modified by T.Nozaki  (2005/09/21)
// Modified by S.Nishida (2007/02/03)
// Modified by T.Nozaki (2009/07/01) Comments in English
#include <cstdio>
#include <cmath>
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
class particleData
{
    public:
    int nParticle;
    Float_t pt[100000];
    Float_t pid[100000];
    Float_t eta[100000];
    Float_t theta[100000];
    Float_t phi[100000];
    Float_t mass[10000];
};

void analysis ( char* file="../Inputs/hadron-2.root", int maxevt = 0 ) {
   /* These are initialization codes.
     You can ignore them as a black box. */
  // Open a data file

  TFile *f = new TFile( file );
  // Obtain a pointer to a series of "event" data in the file
  TTree *t = (TTree*)f->Get("T");
  // Create a pointer to "BEvent" object where data are loaded from the file
  BEvent* event = new BEvent();
  // Obtain a branch for "BEvent" in the tree
  TBranch *branch = t->GetBranch("BEvent");
  // register created "BEvent" object to the branch to load event data
  branch->SetAddress(&event);
   
  // output file
  TFile *hf = new TFile("../Inputs/output_2_withtheta.root", "RECREATE" );
  TTree *tout = new TTree("t","");
  
  particleData pData;
  tout->Branch("nParticle", &pData.nParticle,"nParticle/I");
  tout->Branch("pt", pData.pt,"pt[nParticle]/F");
  tout->Branch("eta", pData.eta,"eta[nParticle]/F");
  tout->Branch("theta", pData.theta,"theta[nParticle]/F");
  tout->Branch("pid", pData.pid,"pid[nParticle]/F");
  tout->Branch("phi", pData.phi,"phi[nParticle]/F");
  tout->Branch("mass", pData.mass,"mass[nParticle]/F");
  
  int nevent = (int)t->GetEntries();
  int nevent_process = nevent;
  if( maxevt>0 && maxevt<nevent ) nevent_process = maxevt;

  // Analysis is carried out event by event.
  // loop of events.
  for ( int i=0;i<nevent_process;i++ ) {
    //  Let us write an event number every 10000 events.
    // It helps for debugging when the program is dead.
    if ( i%10000 == 0 ) printf ( "***** Event %d\n", i );
    // An event is read.
    int nb = t->GetEntry ( i );
    // If the comment delimeter in the following sentence is removed, 
    // the detailed information of the first event is shown.
    /* -------- The main body of the analysis codes--------------- */
    // Obtain the number of particles in each event.
    int nparticles = event->NParticles();
  
     // It is stored in the histogram named h_nprt.
     // Particles are treated as an array.
    TClonesArray& plist = *(event->GetParticleList());
    TLorentzVector v;
    // loop of all particles.
    
    pData.nParticle=nparticles;
    
    for ( int j=0;j<nparticles;j++ ) {
      // The j-th particle is named  p1.
      BParticle* p1 = (BParticle*)plist[j];
      int pid1 = p1->pid();
      pData.pid[j]=pid1;
      // Obtain the momentum and energy of a particle.
      float px1 = p1->px();  // x component of momentum.
      float py1 = p1->py();  // y component of momentum.
      float pz1 = p1->pz();  // z component of momentum.
      float  e1 = p1->e();   // Energy.
      float mass1 = sqrt ( e1*e1-px1*px1-py1*py1-pz1*pz1 ); // Calculate the mass.
      v.SetXYZM(px1,py1,pz1,mass1);
      float eta1 = v.PseudoRapidity();
      float phi1 = v.Phi();
      pData.pt[j]=v.Pt();
      pData.theta[j]=v.Theta();
      pData.eta[j]=v.PseudoRapidity();
      pData.phi[j]=v.Phi();    
      pData.mass[j]=mass1;  
    }
    tout->Fill();
    // Clear at the end of particle loop.
    event->Clear();
  }
  
  // Let us show how many events are processed.
  printf( "***** %d events processed. Exit.\n", nevent_process );
  // Let us store  all the histogram data in the file named hist
  hf->Write();
  // Display histograms.
//  showhist();
}   
// The function showhist is used when the histograms are displayed.
void showhist() {
  TFile *f = new TFile ( "histfile" );
  TBrowser *t = new TBrowser;
}

// for ROOT interactive job
void convert ( char* file="../Inputs/hadron-2.root", int maxevt = 0 ) {
   analysis(file,maxevt);
}
