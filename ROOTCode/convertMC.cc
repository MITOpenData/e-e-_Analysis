#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TString.h>
#include <stdio.h> 
#include <TLegendEntry.h>
#include <TGraphAsymmErrors.h>

class particleData
{
    public:
    int nParticle;
    Float_t pt[100000];
    Float_t pid[100000];
    Float_t eta[100000];
    Float_t phi[100000];
    Float_t mass[10000];
};


void convertMC(){

  const TString forestinput = "HiForestAOD.root";
  TFile *lFile = TFile::Open(forestinput);
  TTree *gentree = (TTree*)lFile->Get("HiGenParticleAna/hi");
  
  Int_t _mult;
  
  vector<float> *_pt;
  vector<float> *_eta;
  vector<float> *_pdg;
  vector<float> *_phi;
  
  gentree->SetBranchAddress("mult",&_mult);
  gentree->SetBranchAddress("pt",&_pt);
  gentree->SetBranchAddress("eta",&_eta);
  gentree->SetBranchAddress("phi",&_phi);
  gentree->SetBranchAddress("pdg",&_pdg);
  
  
  TFile *hf = new TFile("outputMC.root", "RECREATE" );
  TTree *tout = new TTree("t","");
  
  particleData pData;
  tout->Branch("nParticle", &pData.nParticle,"nParticle/I");
  tout->Branch("pt", pData.pt,"pt[nParticle]/F");
  tout->Branch("eta", pData.eta,"eta[nParticle]/F");
  tout->Branch("pid", pData.pid,"pid[nParticle]/F");
  tout->Branch("phi", pData.phi,"phi[nParticle]/F");
  tout->Branch("mass", pData.mass,"mass[nParticle]/F");


  int entries = gentree->GetEntries();
  
  for (Int_t i=0;i<entries;i++) {
    gentree->GetEntry(i);
    pData.nParticle=_mult;
    
      for ( int j=0;j<_mult;j++ ) {
     
      pData.pid[j]=_pdg.at(j);
      pData.pt[j]=_pt.at(j);
      pData.eta[j]=_eta.at(j);
      pData.phi[j]=_phi.at(j);
      pData.mass[j]=0;
      cout<<"pt="<<_pt.at(j)<<endl;
     }
     tout->Fill();
  }
  hf->Write();
}