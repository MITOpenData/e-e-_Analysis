//this takes in the virtual n-tuple created by rand_angle.cxx 
//and adjusts the px, py and py values so that the mean of these
//parameters accross an event is 0

#include <random>
#include <iostream>
#include <math.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#define PI 3.1415926

using namespace std;

class particleData
{
	public:
		int nParticle;
		int EventNo;
		int RunNo;
		float Energy;

		Float_t px[100000];
		Float_t py[100000];
		Float_t pz[100000];
		Float_t pt[100000];
    Float_t pmag[100000];
		Float_t eta[100000];
    Float_t theta[100000];
		Float_t phi[100000];
		Float_t mass[100000];
		Float_t charge[100000];
		Float_t pwflag[10000];
		Float_t pid[10000];
};



void scan()
{
  TString filename = "new_data8.root";
  
  TFile *f = new TFile(filename.Data());
  TTree *t1 = (TTree*)f->Get("t");
  
  int nParticle1;
  Float_t px1[10000];
  Float_t py1[10000];
  Float_t pz1[10000];
  
  //Note, we use px1, py1 and pz1 and not px, py and pz.
  //We have reserved px, py and pz for the new n-tuple we seek to create
  
  t1->SetBranchAddress("nParticle",&nParticle1);  
  t1->SetBranchAddress("px", px1);
  t1->SetBranchAddress("py", py1);
  t1->SetBranchAddress("pz", pz1);
  
  t1->GetEntries(0);
  
  /*
  double mean_px1 = 0;
  double mean_py1 = 0;
  double mean_pz1 = 0;
  
  for (int i = 0; i <1000; i++)
  {
   mean_px1 += (double) px1[i];
   //cout<<px1[i]<<endl;
   mean_py1 += py1[i];
   mean_pz1 += pz1[i];
  }
  mean_px1 = mean_px1/(float)1000;
  mean_py1 = mean_py1/(float)1000;
  mean_pz1 = mean_pz1/(float)1000;
  
  */
 
  t1->Draw("px>>h1");
  TH1F *h1 = (TH1F*)gDirectory->Get("h1");

  t1->Draw("py>>h2");
  TH1F *h2 = (TH1F*)gDirectory->Get("h2");

  t1->Draw("pz>>h3");
  TH1F *h3 = (TH1F*)gDirectory->Get("h3"); 


  float mean_px1 = h1->GetMean();
  float mean_py1 = h2->GetMean();  
  float mean_pz1 = h3->GetMean();
  
  TFile *hf = new TFile("good_data8.root", "RECREATE" );
	TTree *tout = new TTree("t","");
	
  
	
  particleData pData;
	tout->Branch("nParticle", &pData.nParticle,"nParticle/I");
	//tout->Branch("EventNo", &pData.EventNo,"EventNo/I");
	//tout->Branch("RunNo", &pData.RunNo,"RunNo/I");
	//tout->Branch("Energy", &pData.Energy,"Energy/F");
	tout->Branch("px", pData.px,"px[nParticle]/F");
	tout->Branch("py", pData.py,"py[nParticle]/F");
	tout->Branch("pz", pData.pz,"pz[nParticle]/F");
	tout->Branch("pt", pData.pt,"pt[nParticle]/F");
  tout->Branch("pmag", pData.pmag,"pmag[nParticle]/F");
	//tout->Branch("mass", pData.mass,"mass[nParticle]/F");
	tout->Branch("eta", pData.eta,"eta[nParticle]/F");
  tout->Branch("theta", pData.theta,"theta[nParticle]/F");
	tout->Branch("phi", pData.phi,"phi[nParticle]/F");
	//tout->Branch("charge", pData.charge,"charge[nParticle]/F");
	//tout->Branch("pwflag", pData.pwflag,"pwflag[nParticle]/F");
	//tout->Branch("pid", pData.pid,"pid[nParticle]/F");
  
  pData.nParticle = nParticle1;
  
  for (int j=0; j<nParticle1; j++)
  {
    pData.px[j] = px1[j] - mean_px1;
    pData.py[j] = py1[j] - mean_py1;
    pData.pz[j] = pz1[j] - mean_pz1;
  }

  TVector3 v;
  
  for (int counterParticles = 0; counterParticles < pData.nParticle; counterParticles++)
  {
    
    v.SetXYZ(pData.px[counterParticles], pData.py[counterParticles], pData.pz[counterParticles]);
    
    pData.pt[counterParticles]= v.Pt();
    pData.pmag[counterParticles]=v.Mag();
    pData.eta[counterParticles]=v.PseudoRapidity();
    pData.theta[counterParticles]=v.Theta();
    pData.phi[counterParticles]=v.Phi();
    
	}
  tout->Fill();
  hf->Write();
  cout<<"complete"<<endl;
}
