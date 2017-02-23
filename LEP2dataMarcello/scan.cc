#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include "TLorentzVector.h"


using namespace std;

class particleData
{
	public:
		int nParticle;

		Float_t px[100000];
		Float_t py[100000];
		Float_t pz[100000];
		Float_t pt[100000];
		Float_t eta[100000];
		Float_t phi[100000];
		Float_t mass[100000];
		Float_t charge[10000];
		Float_t pwflag[10000];
		Float_t pid[10000];
};

void scan (TString infile="cleaned_ALEPH_Data1998_189GeV_V0.txt"){

	FILE *fp=fopen(infile.Data(),"r");

	float _px,_py,_pz,_m,_charge,_pwflag;

	TFile *hf = new TFile("myALEPH.root", "RECREATE" );
	TTree *tout = new TTree("t","");
	TLorentzVector v;

	particleData pData;
	tout->Branch("nParticle", &pData.nParticle,"nParticle/I");
	tout->Branch("px", pData.px,"px[nParticle]/F");
	tout->Branch("py", pData.py,"py[nParticle]/F");
	tout->Branch("pz", pData.pz,"pz[nParticle]/F");
	tout->Branch("pt", pData.pt,"pt[nParticle]/F");
	tout->Branch("mass", pData.mass,"mass[nParticle]/F");
	tout->Branch("eta", pData.eta,"eta[nParticle]/F");
	tout->Branch("phi", pData.phi,"phi[nParticle]/F");
	tout->Branch("charge", pData.charge,"charge[nParticle]/F");
	tout->Branch("pwflag", pData.pwflag,"pwflag[nParticle]/F");
	tout->Branch("pid", pData.pid,"pid[nParticle]/F");

	int counterEntries=0;
	int counterParticles=0;

	while(fscanf(fp,"%f %f %f %f %f %f",&_px,&_py,&_pz,&_m,&_charge,&_pwflag)!=EOF) {
		if (_px==-999.&&_py==-999.&&_pz==-999.) { 
		   //cout<<"counterEntries="<<counterEntries<<endl;
			pData.nParticle=counterParticles; 
            if(counterEntries>0) tout->Fill(); 
			counterParticles=0;   
			//cout<<"------------------------------------------------------------------------"<<endl; 
			continue;
		}  
		pData.px[counterParticles]=_px;
		pData.py[counterParticles]=_py;
		pData.pz[counterParticles]=_pz;
		pData.mass[counterParticles]=_m;
		v.SetXYZM(_px,_py,_pz,_m);
        pData.pt[counterParticles]=v.Pt();
        pData.eta[counterParticles]=v.PseudoRapidity();
        pData.phi[counterParticles]=v.Phi();    
		pData.charge[counterParticles]=_charge;
		pData.pwflag[counterParticles]=_pwflag;
		pData.pid[counterParticles]=-1;
		counterParticles++;	
		counterEntries++;	
	}
	hf->Write();
}


void check(){

  TFile *f = new TFile("../LEP2dataMarcello/myALEPH.root");
  TTree *t1 = (TTree*)f->Get("t");
  Int_t nParticle;
  Float_t px[50000];
  Float_t py[50000];
  Float_t pz[50000];
  Float_t pt[50000];
  Float_t eta[50000];
  Float_t pid[50000];
  Float_t phi[50000];
  Float_t mass[50000];
  
  t1->SetBranchAddress("nParticle",&nParticle);
  t1->SetBranchAddress("px",px);
  t1->SetBranchAddress("py",py);
  t1->SetBranchAddress("pz",pz);
  t1->SetBranchAddress("pt",pt);
  t1->SetBranchAddress("eta",eta);
  t1->SetBranchAddress("pid",pid);
  t1->SetBranchAddress("phi",phi);
  t1->SetBranchAddress("mass",mass);

  Int_t nevent = (Int_t)t1->GetEntries();
    
  for (Int_t i=0;i<nevent;i++) {
    t1->GetEntry(i);
    
    int nparticles = nParticle;
    cout<<"------------------------------------------------------------------------"<<endl;
    cout<<"nparticles="<<nparticles<<endl;

    for ( int j=0;j<nparticles;j++ ) {
    
      int pid1 = pid[j];
     // if (pid1!=PION&&pid1!=PROTON&&pid1!=KAON) continue;
      float eta1 = eta[j];
      float phi1 = phi[j];
      float px1 = px[j];
      float py1 = py[j];
      float pz1 = pz[j];
      float pt1 = pt[j];
      float mass1 = mass[j];
      cout<<"px="<<px1<<"py="<<py1<<"pz="<<pz1<<"mass="<<mass1<<endl;
   }
 }
}