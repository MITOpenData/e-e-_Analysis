#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include <iostream>
#include "TLorentzVector.h"

#define MAX_NPARTICLE 1000

class EvtInfoBranches {
public:
  double beam_px;
  double beam_pz;
  double beam_e;

  void Register(TTree *root) {
    root->SetBranchAddress("EvtInfo.beam_px"     , &beam_px     );
    root->SetBranchAddress("EvtInfo.beam_pz"     , &beam_pz     );
    root->SetBranchAddress("EvtInfo.beam_e"      , &beam_e      );
  };
};

class particleData
{
public:
  int nParticle;
  Float_t pt[MAX_NPARTICLE];
  Float_t pid[MAX_NPARTICLE];
  Float_t eta[MAX_NPARTICLE];
  Float_t theta[MAX_NPARTICLE];
  Float_t phi[MAX_NPARTICLE];
  Float_t mass[MAX_NPARTICLE];
};

void boost ( TString file="charmmc-e07-00.root", int maxevt = 0 ) {
  /* These are initialization codes.
     You can ignore them as a black box. */

  // Open a data file
  TFile *f = new TFile( file );
  // Obtain a pointer to a series of "event" data in the file
  TTree *t = (TTree*)f->Get("t");

  particleData inP;
  t->SetBranchAddress("nParticle", &inP.nParticle);
  t->SetBranchAddress("pt", inP.pt);
  t->SetBranchAddress("eta", inP.eta);
  t->SetBranchAddress("theta", inP.theta);
  t->SetBranchAddress("pid", inP.pid);
  t->SetBranchAddress("phi", inP.phi);
  // t->SetBranchStatus("mass", inP.mass);

  TString foutname = file(0, file.Length() - 5) + "_labframe.root";
  // Open an output root file
  TFile *fout = new TFile(foutname, "RECREATE");
  // output tree
  TTree *tout = new TTree("t", "momentums in lab frame ");

  particleData pData;
  tout->Branch("nParticle", &pData.nParticle,"nParticle/I");
  tout->Branch("pt", pData.pt,"pt[nParticle]/F");
  tout->Branch("eta", pData.eta,"eta[nParticle]/F");
  tout->Branch("theta", pData.theta,"theta[nParticle]/F");
  tout->Branch("pid", pData.pid,"pid[nParticle]/F");
  tout->Branch("phi", pData.phi,"phi[nParticle]/F");
  tout->Branch("mass", pData.mass,"mass[nParticle]/F");

  memset(&pData, 0, sizeof(pData));

  EvtInfoBranches evt;
  evt.Register(t);

  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("EvtInfo.beam*", 1);
  t->SetBranchStatus("nParticle", 1);
  t->SetBranchStatus("pt", 1);
  t->SetBranchStatus("eta", 1);
  t->SetBranchStatus("theta", 1);
  t->SetBranchStatus("pid", 1);
  t->SetBranchStatus("phi", 1);

  // end of initialization section

  // Loop over all events
  for (auto ievt = 0; ievt < t->GetEntries(); ++ievt) {
    //  Let us write an event number every 10000 events.
    // It helps for debugging when the program is dead.
    if ( ievt % 10000 == 0 ) {
      std::cout << "***** Event " << ievt << std::endl;
    }

    t->GetEntry(ievt);

    // Get the beam 4-vector.
    TLorentzVector beam(evt.beam_px, 0, evt.beam_pz, evt.beam_e);

    pData.nParticle = inP.nParticle;
    for (auto i = 0; i < inP.nParticle; ++i) {
      // Get the frame-independent variables
      pData.pid[i] = inP.pid[i];

      TLorentzVector p;
      p.SetPtEtaPhiM(inP.pt[i], inP.eta[i], inP.phi[i], inP.mass[i]);
      // Boost from the CM frame to the lab frame
      p.Boost(-beam.BoostVector());

      // Get the kinematic in the lab frame
      pData.pt[i] = p.Pt();
      pData.eta[i] = p.Eta();
      pData.theta[i] = p.Theta();
      pData.phi[i] = p.Phi();
      pData.mass[i] = p.M();
    }

    tout->Fill();
  }
  fout->Write();

  std::cout << "output file " << foutname << " saved." << std::endl;
}
