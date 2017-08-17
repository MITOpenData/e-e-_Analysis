//
//  convertCMSTrackTree.cc
//  
//
//  Created by Anthony Badea on 8/9/17.
//
//

#include <stdio.h>
#include <TFile.h>
#include <TTree.h>
#include "TLorentzVector.h"
#include <TString.h>
void copytree()
{
    //gSystem->Load("$ROOTSYS/test/libEvent");
    
  const Double_t piPlusMass = 0.13957018;

    TFile *oldfile = new TFile("MinBias_TuneCUETP8M1_5p02TeV-pythia8-HINppWinter16DR-NoPU_75X_mcRun2_asymptotic_ppAt5TeV_forest_v2_track.root");
    TTree *t1 = (TTree*)oldfile->Get("ppTrack/trackTree");
    // Now Open the branches and then put them into a new file
    const Int_t nMaxTrk = 50000;
    Int_t nTrk;
    Int_t nEv;
    Int_t nRun;
    Float_t trkPt[nMaxTrk];
    Float_t trkPhi[nMaxTrk];
    Float_t trkEta[nMaxTrk];
    Int_t trkCharge[nMaxTrk];
    
    t1->SetBranchAddress("nTrk",&nTrk);
    t1->SetBranchAddress("nEv",&nEv);
    t1->SetBranchAddress("nRun",&nRun);
    t1->SetBranchAddress("trkPt",trkPt);
    t1->SetBranchAddress("trkPhi",trkPhi);
    t1->SetBranchAddress("trkEta",trkEta);
    t1->SetBranchAddress("trkCharge",trkCharge);
    
    Int_t nevent = (Int_t)t1->GetEntries();
    TFile *newfile = new TFile("FORMATTED_MinBias_TuneCUETP8M1_5p02TeV-pythia8-HINppWinter16DR-NoPU_75X_mcRun2_asymptotic_ppAt5TeV_forest_v2_track.root","recreate");
    TTree *tout = new TTree("t","");
    //cout<<"step1"<<endl;
    const Int_t nMaxPart = 100000;
    Int_t nParticles, EventNo, RunNo;
    Double_t Energy;
    Float_t px[nMaxPart],py[nMaxPart],pz[nMaxPart],pt[nMaxPart],mass[nMaxPart],theta[nMaxPart],eta[nMaxPart],phi[nMaxPart],charge[nMaxPart],pwflag[nMaxPart],pid[nMaxPart];

    tout->Branch("nParticle", &nParticles,"nParticle/I");
    tout->Branch("EventNo", &EventNo,"EventNo/I");
    tout->Branch("RunNo", &RunNo,"RunNo/I");
    tout->Branch("Energy", &Energy,"Energy/F");
    tout->Branch("px", px,"px[nParticle]/F");
    tout->Branch("py", py,"py[nParticle]/F");
    tout->Branch("pz", pz,"pz[nParticle]/F");
    tout->Branch("pt", pt,"pt[nParticle]/F");
    tout->Branch("mass", mass,"mass[nParticle]/F");
    tout->Branch("theta",theta,"theta[nParticle]/F");
    tout->Branch("eta", eta,"eta[nParticle]/F");
    tout->Branch("phi", phi,"phi[nParticle]/F");
    tout->Branch("charge", charge,"charge[nParticle]/F");
    tout->Branch("pwflag", pwflag,"pwflag[nParticle]/F");
    tout->Branch("pid", pid,"pid[nParticle]/F");
    cout<<nevent<<endl;
    //cout<<nTrk[53360]<<endl;
    for (Int_t i=0;i<nevent;i++)
    {
        t1->GetEntry(i);
        nParticles = nTrk;
        EventNo = nEv;
        RunNo = nRun;
        Energy = 5.02; // TeV

	for(int trkI = 0; trkI < nTrk; ++trkI)
    {
	  pt[trkI] = trkPt[trkI];
	  phi[trkI] = trkPhi[trkI];
	  eta[trkI] = trkEta[trkI];
	  theta[trkI] = 2*TMath::ATan(TMath::Exp(-trkEta[trkI]));
	  charge[trkI] = trkCharge[trkI];
	  pwflag[trkI] = 0;
	  // pData.pid =
	  mass[trkI] = piPlusMass;
	  TLorentzVector v;
	  v.SetPtEtaPhiM(pt[trkI], eta[trkI], phi[trkI], piPlusMass);
	  px[trkI] = v.Px();
	  py[trkI] = v.Py();
	  pz[trkI] = v.Pz();
    }
        tout->Fill();
    }
    tout->Write();
    newfile->Write();
    
}
