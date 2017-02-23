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
		Float_t m[100000];
		Float_t charge[10000];
		Float_t pwflag[10000];
		Int_t pid[10000];
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
	tout->Branch("m", pData.m,"m[nParticle]/F");
	tout->Branch("eta", pData.eta,"eta[nParticle]/F");
	tout->Branch("phi", pData.phi,"phi[nParticle]/F");
	tout->Branch("charge", pData.charge,"charge[nParticle]/F");
	tout->Branch("pwflag", pData.pwflag,"pwflag[nParticle]/F");
	tout->Branch("pid", pData.pid,"pid[nParticle]/I");

	int counterEntries=0;
	int counterParticles=0;

	while(fscanf(fp,"%f %f %f %f %f %f",&_px,&_py,&_pz,&_m,&_charge,&_pwflag)!=EOF) {
		if (_px==-999.&&_py==-999.&&_pz==-999.) { 
		   cout<<"counterEntries="<<counterEntries<<endl;
		    if(counterEntries>0) tout->Fill(); 
			pData.nParticle=counterParticles; 
			counterParticles=0;   
			cout<<"------------------------------------------------------------------------"<<endl; 
			continue;
		}  
		pData.px[counterParticles]=_px;
		pData.py[counterParticles]=_py;
		pData.pz[counterParticles]=_pz;
		pData.m[counterParticles]=_m;
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
