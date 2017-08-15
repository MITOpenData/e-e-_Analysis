#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include <iostream>

#define MAX_NGEN 2048    // GEN : generator info

class GenInfoBranches {
public:
  int idhep[MAX_NGEN];
  int daFirst[MAX_NGEN];
  int daLast[MAX_NGEN];

  void Register(TTree *root){
    root->SetBranchAddress("GenInfo.idhep"   , &idhep[0]   );
    root->SetBranchAddress("GenInfo.daFirst" , &daFirst[0] );
    root->SetBranchAddress("GenInfo.daLast"  , &daLast[0]  );

    root->SetBranchStatus("GenInfo.idhep", 1);
    root->SetBranchStatus("GenInfo.daFirst", 1);
    root->SetBranchStatus("GenInfo.daLast", 1);
  };
};

  void flavor ( TString file="../qqmc-e07-00.root", int maxevt = 0 ) {
  /* These are initialization codes.
     You can ignore them as a black box. */

  // Open a data file
  TFile *f = new TFile( file );
  // Obtain a pointer to a series of "event" data in the file
  TTree *t = (TTree*)f->Get("T");

  TString foutname = file(0, file.Length() - 5) + "_flavor.root";
  // Open an output root file
  TFile *fout = new TFile(foutname, "RECREATE");
  // output tree
  TTree *tout = new TTree("t", "extended data");

  // flavor of the qqbar event
  Int_t flavor = 0;

  tout->Branch("flavor", &flavor, "flavor/I");

  root->SetBranchStatus("*", 0);
  // Get the generator information
  GenInfoBranches gen;
  gen.Register(t);

  // Loop over all events
  for (auto evt = 0; evt < t->GetEntries(); ++evt) {
    //  Let us write an event number every 10000 events.
    // It helps for debugging when the program is dead.
    if ( evt % 10000 == 0 ) std::cout << "***** Event " << evt << std::endl;

    t->GetEntry(evt);
    // Check whether this looks like a qqbar event
    if (gen.idhep[0] != 10022) {
      std::cout << "The first generated particle is not a virtual photon (" << gen.idhep[0] << "). Is this a qqbar MC sample?" << std::endl;
    }

    // Check the flavor of the daughters of the virtual photon
    flavor = 0;
    for (auto i = gen.daFirst[0]; i <= gen.daLast[0]; i++) {
      if (gen.idhep[i] == 4) {
        flavor = 1;
        break;
      }
    }

    tout->Fill();
  }
  fout->Write();

  std::cout << "output file " << foutname << " saved." << std::endl;
}
