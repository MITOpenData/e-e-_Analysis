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
		int EventNo;
		int RunNo;
		float Energy;

		Float_t px[100000];
		Float_t py[100000];
		Float_t pz[100000];
		Float_t pt[100000];
    Float_t pmag[100000];//Added later on
		Float_t eta[100000];
    Float_t theta[100000];
		Float_t phi[100000];
		Float_t mass[100000];
		Float_t charge[10000];
		Float_t pwflag[10000];
		Float_t pid[10000];
};

void scan (TString infile="LEP2MCGGCCY1997E183_recons_aftercut-001.aleph"){


  // "/data/flowex/Datasamples/LEP2_MAIN/ROOTfiles/cleaned_ALEPH_Data-all.aleph"
	FILE *fp=fopen(Form("%s",infile.Data()),"r");

	float _px,_py,_pz,_m,_charge,_pwflag;

	TFile *hf = new TFile(Form("%s.root",infile.Data()), "RECREATE" );
	TTree *tout = new TTree("t","");
	TLorentzVector v;
    cout<<"step1"<<endl;
	particleData pData;
	tout->Branch("nParticle", &pData.nParticle,"nParticle/I");
	tout->Branch("EventNo", &pData.EventNo,"EventNo/I");
	tout->Branch("RunNo", &pData.RunNo,"RunNo/I");
	tout->Branch("Energy", &pData.Energy,"Energy/F");
	tout->Branch("px", pData.px,"px[nParticle]/F");
	tout->Branch("py", pData.py,"py[nParticle]/F");
	tout->Branch("pz", pData.pz,"pz[nParticle]/F");
	tout->Branch("pt", pData.pt,"pt[nParticle]/F");
  tout->Branch("pmag", pData.pmag,"pmag[nParticle]/F");//Added later on
	tout->Branch("mass", pData.mass,"mass[nParticle]/F");
	tout->Branch("eta", pData.eta,"eta[nParticle]/F");
  tout->Branch("theta", pData.theta,"theta[nParticle]/F");
	tout->Branch("phi", pData.phi,"phi[nParticle]/F");
	tout->Branch("charge", pData.charge,"charge[nParticle]/F");
	tout->Branch("pwflag", pData.pwflag,"pwflag[nParticle]/F");
	tout->Branch("pid", pData.pid,"pid[nParticle]/F");

	int counterEntries=0;
	int counterParticles=0;

  cout<<"step2"<<endl;
	while(fscanf(fp,"%f %f %f %f %f %f",&_px,&_py,&_pz,&_m,&_charge,&_pwflag)!=EOF) {
		if (_px==-999.&&_py==-999.&&_pz==-999.) { 
		  cout<<"counterEntries="<<counterEntries<<endl;
			pData.nParticle=counterParticles; 
      if(counterEntries>0) tout->Fill(); 
		  pData.RunNo=_m; 
		  pData.EventNo=_charge; 
			pData.Energy=_pwflag; 
			counterParticles=0;   
			cout<<"------------------------------------------------------------------------"<<endl; 
			continue;
		}  
		cout<<_px<<"    "<<_py<<"    "<<_pz<<"    "<<endl;
		pData.px[counterParticles]=_px;
		pData.py[counterParticles]=_py;
		pData.pz[counterParticles]=_pz;
		pData.mass[counterParticles]=_m;
		v.SetXYZM(_px,_py,_pz,_m);
        pData.pt[counterParticles]=v.Pt();
        pData.pmag[counterParticles]=v.Rho(); //Added later on
        pData.eta[counterParticles]=v.PseudoRapidity();
        pData.theta[counterParticles]=v.Theta();
        pData.phi[counterParticles]=v.Phi();
        
		pData.charge[counterParticles]=_charge;
		pData.pwflag[counterParticles]=_pwflag;
		pData.pid[counterParticles]=0;
		counterParticles++;	
		counterEntries++;	
	}
	hf->Write();
}


void check(){

  // THIS FILE LOCATION NEEDS TO BE UPDATED TO SVMITHI02 ADDRESS
  // THIS WAS NOT DONE BECAUSE THE FILE WAS NOT IN THE ORIGINAL GIT FILES
  TFile *f = new TFile("/data/flowex/Datasamples/LEP2_MAIN/ROOTfiles/cleaned_ALEPH_Data1999_200GeV_V0.root");
  TTree *t1 = (TTree*)f->Get("t");
  Int_t nParticle,EventNo,RunNo;
  Float_t Energy;
  Float_t px[100000];
  Float_t py[100000];
  Float_t pz[100000];
  Float_t pt[100000];
  Float_t pmag[100000];//Added later on
  Float_t eta[100000];
  Float_t pid[100000];
  Float_t phi[100000];
  Float_t mass[100000];
  
  t1->SetBranchAddress("nParticle",&nParticle);
  t1->SetBranchAddress("EventNo",&EventNo);
  t1->SetBranchAddress("RunNo",&RunNo);
  t1->SetBranchAddress("Energy",&Energy);
  t1->SetBranchAddress("px",px);
  t1->SetBranchAddress("py",py);
  t1->SetBranchAddress("pz",pz);
  t1->SetBranchAddress("pt",pt);
  t1->SetBranchAddress("pmag",pmag);//Added later on
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
      float pmag1 = pmag[j];
      float mass1 = mass[j];
      cout<<"EventNo="<<EventNo<<" px="<<px1<<"py="<<py1<<"pz="<<pz1<<"pt="<<pt1<<"pmag="<<pmag1<<"mass="<<mass1<<endl;
   }
 }
}



int main(int argc, char *argv[])
{
  if((argc != 2))
  {
    std::cout << "Wrong number of inputs" << std::endl;
    return 1;
  }
  
  if(argc == 2)
    std::cout << "Running for the sake of god" << std::endl;
    scan(argv[1]);
  return 0;
}
