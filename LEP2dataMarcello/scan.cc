#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>

using namespace std;

class particleData
{
    public:
    int nParticle;
    Float_t px[100000];
    Float_t py[100000];
    Float_t pz[100000];
    Float_t m[100000];
    Float_t charge[10000];
    Float_t pwflag[10000];
};

void scan (TString infile="final.txt",TString output="PredictionsCUJET3_pt_0_10.root"){

	FILE *fp=fopen(infile.Data(),"r");

	float _px,_py,_pz,_m,_charge,_pwflag;
	
    TFile *hf = new TFile("myALEPH.root", "RECREATE" );
    TTree *tout = new TTree("t","");
  
    particleData pData;
    tout->Branch("nParticle", &pData.nParticle,"nParticle/I");
    tout->Branch("px", pData.px,"px[nParticle]/F");
    tout->Branch("py", pData.py,"py[nParticle]/F");
    tout->Branch("pz", pData.pz,"pz[nParticle]/F");
    tout->Branch("m", pData.m,"m[nParticle]/F");
    tout->Branch("charge", pData.charge,"charge[nParticle]/F");
    tout->Branch("pwflag", pData.pwflag,"pwflag[nParticle]/F");
    
    
    int counterEntries=0;
    int counterParticles=0;

	while(fscanf(fp,"%f %f %f %f %f %f",&_px,&_py,&_pz,&_m,&_charge,&_pwflag)!=EOF) {
	if (_px==-999.&&_py==-999.&&_pz==-999.) { 
	  pData.nParticle=counterParticles; counterParticles=0;  
	  if(counterEntries>0) tout->Fill(); 
	  cout<<"------------------------------------------------------------------------"<<endl; 
	  continue;
	  }  
      pData.px[counterParticles]=_px;
      pData.py[counterParticles]=_py;
      pData.pz[counterParticles]=_pz;
      pData.m[counterParticles]=_m;
      pData.charge[counterParticles]=_charge;
      pData.pwflag[counterParticles]=_pwflag;
      cout<<"px="<<_px<<endl;
      counterParticles++;	
      counterEntries++;	
}
hf->Write();
}