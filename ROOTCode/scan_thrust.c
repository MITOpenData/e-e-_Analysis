#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include "TLorentzVector.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <stdlib.h>
#include <TMath.h>
#include <TVector3.h>


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
		Float_t eta[100000];
		Float_t phi[100000];
		Float_t mass[100000];
		Float_t charge[10000];
		Float_t pwflag[10000];
		Float_t pid[10000];
};

void scan (TString infile="cleaned_LEP2MCGGUD-MCtrue-all.aleph"){

	FILE *fp=fopen(Form("%s",infile.Data()),"r");

	float _px,_py,_pz,_m,_charge,_pwflag;

	TFile *hf = new TFile(Form("ROOTfiles/%s.root",infile.Data()), "RECREATE" );
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
	tout->Branch("mass", pData.mass,"mass[nParticle]/F");
	tout->Branch("eta", pData.eta,"eta[nParticle]/F");
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
        pData.eta[counterParticles]=v.PseudoRapidity();
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

  TFile *f = new TFile("ROOTfiles/cleaned_ALEPH_Data1999_200GeV_V0.root");
  TTree *t1 = (TTree*)f->Get("t");
  Int_t nParticle,EventNo,RunNo;
  Float_t Energy;
  Float_t px[100000];
  Float_t py[100000];
  Float_t pz[100000];
  Float_t pt[100000];
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
      cout<<"EventNo="<<EventNo<<" px="<<px1<<"py="<<py1<<"pz="<<pz1<<"pt="<<pt1<<"mass="<<mass1<<endl;
   }
 }
}





#define PI 3.1415926
enum SIMPLEPID {PHOTON, ELECTRON, PION, MUON, KAON, PROTON};


double dphi(double phi1,double phi2)
{
    double a=phi1-phi2;
    while (a<-PI) a+=2*PI;
    while (a>PI) a-=2*PI;
    
    if (a<-PI/2) a=2*PI+a;
    return a;
}

float thrust(float x[]){

  TString filename="/mnt/c/Users/Bibek Kumar Pandit/Desktop/Root_Directory/StudyMult/LEP2/ROOTfiles/cleaned_ALEPH_Data-all.aleph.root";
  
  TFile *f = new TFile(filename.Data());
  TTree *t1 = (TTree*)f->Get("t");
  Int_t nParticle;
  Float_t pt[50000];
  Float_t eta[50000];
  Float_t theta[50000];
  Float_t pid[50000];
  Float_t phi[50000];
  Float_t mass[50000];
  Float_t pwflag[50000];
  
  t1->SetBranchAddress("nParticle",&nParticle);
  t1->SetBranchAddress("pt",pt);
  t1->SetBranchAddress("eta",eta);
  t1->SetBranchAddress("theta",theta);
  t1->SetBranchAddress("pid",pid);
  t1->SetBranchAddress("phi",phi);
  t1->SetBranchAddress("mass",mass);
  t1->SetBranchAddress("pwflag",pwflag);

  /*
  x=sin(theta)*cos(phi)
  y=sin(theta)*sin(phi)
  z=cos(theta)
  */
  
  int nstepstheta=100;
  double mintheta=0.;
  double maxtheta=3.14/2;  
  
  int nstepsphi=100;
  double minphi=0.;
  double maxphi=3.14;
  
  double sizetheta=(maxtheta-mintheta)/(double)nstepstheta;
  double sizephi=(maxphi-minphi)/(double)nstepsphi;

  double mytheta,myphi;
  double x_versor,y_versor,z_versor;
  TVector3 unity_versor;
  TVector3 particle_vector;
  
  double maxscalar,theta_max,phi_max;
  
    // all entries and fill the histograms
  Int_t nevent = (Int_t)t1->GetEntries();
  
  //int nevent_process = nevent;
  //cout<<"n"<<nevent<<endl;
  const int nevent_process = 59874; //59874
  TH2F*hprojection[nevent_process];
  //float x[nevent_process];
  
    for (Int_t i=0;i<nevent_process;i++) {  
    if (i%10==0) cout <<i<<"/"<<nevent_process<<endl;
    t1->GetEntry(i);
    
    hprojection[i]=new TH2F(Form("hprojection%d",i),Form("hprojection%d",i),nstepstheta,mintheta,maxtheta,nstepsphi,minphi,maxphi);

    maxscalar=-1;

    for ( int itheta=0; itheta<nstepstheta; itheta++ ) {
      mytheta=sizetheta*(double)itheta;
      for ( int iphi=0; iphi<nstepsphi; iphi++ ) {
        myphi=sizephi*(double)iphi;

        x_versor=sin(mytheta)*cos(myphi);
        y_versor=sin(mytheta)*sin(myphi);
        z_versor=cos(mytheta);
        
        unity_versor.SetXYZ(x_versor,y_versor,z_versor);
        double scalar_prod;
        
       scalar_prod=0.;
       
       for ( int j=0;j<nParticle;j++ ) {
         float pt1 = pt[j];
         float theta1 = theta[j];
         float phi1 = phi[j];
         float px1=pt1*sin(theta1)*cos(phi1);
         float py1=pt1*sin(theta1)*sin(phi1);
         float pz1=pt1*sin(phi1);
         particle_vector.SetXYZ(px1,py1,pz1);
         scalar_prod=scalar_prod+fabs(particle_vector.Dot(unity_versor));         
       }//end of loop over particles

       hprojection[i]->SetBinContent(itheta,iphi,scalar_prod);       
       if (scalar_prod>maxscalar) { maxscalar=scalar_prod; theta_max=mytheta; phi_max=myphi;} 
      }//loop over theta
   }//loop over phi
   int binmax = hprojection[i]->GetMaximumBin();
   x[i] = hprojection[i]->GetYaxis()->GetBinCenter(binmax);
   
}// end of loop over events
 return x[5];
}



void scan_thrust(){
  
 
  //TFile *f = new TFile(filename.Data());
  //TTree *t1 = (TTree*)f->Get("t");
  
 	//TFile *hf = TFile(Form("ROOTfiles/%s.root",infile.Data()), "RECREATE" );
	//TTree *tout = new TTree("t","");
	//TLorentzVector v;
  //  cout<<"step1"<<endl;
  
  /*
  TFile *f = new TFile("/mnt/c/Users/Bibek Kumar Pandit/Desktop/Root_Directory/StudyMult/LEP2/ROOTfiles/cleaned_ALEPH_Data-all.aleph.root","update");
  TTree *tout = (TTree*)f->Get("ntuple");
	
  float tAngle;
	TBranch *etAngle = tout->Branch("tAngle", &tAngle,"tAngle/F");
  Int_t nevent = (Int_t)tout->GetEntries();
  
  
  float x[nevent];
  memset(x, 0, nevent*sizeof(float) );
  cout<<"elo"<<endl;
  
  thrust(x);
  
  //tAngle_array = thrust();
  
  for (Int_t i =0; i < nevent; i++)
  {
    if (i%10 == 0) 
    {   
      cout<<i<<endl;
    }
    tAngle = x[i];
    etAngle->Fill();
  }
  tout->Print();
  tout->Write();
  delete f;
  
  
 */
 



  /*
  TFile *f = new TFile("/mnt/c/Users/Bibek Kumar Pandit/Desktop/Root_Directory/StudyMult/LEP2/ROOTfiles/cleaned_ALEPH_Data-all.aleph.root","RECREATE");
   TTree *T = (TTree*)f->Get("t");
  // float px,py;
   float tAngle;
   TBranch *bpt = T->Branch("tAngle",&tAngle,"tAngle/F");
  // T->SetBranchAddress("px",&px);
   //T->SetBranchAddress("py",&py);
   //Long64_t nentries = T->GetEntries
   //Int_t nevent = (Int_t)T->GetEntries();
  
  
   float x[100];
   //memset(x, 0, nevent*sizeof(float) );
   cout<<"elo"<<endl; 
    
    thrust(x);
   // Int_t nevent = (Int_t)tout->GetEntries();
   for (i=0;i<100;i++) {
      if (i%10 == 0) 
    {   
      cout<<i<<endl;
    }
    tAngle = x[i];
    bpt->Fill();
   }
   T->Print();
   T->Write();
   delete f;
   
   */
  
  //Float_t px, py, energy, transmom;
  TString filename;
  filename="/mnt/c/Users/Bibek Kumar Pandit/Desktop/Root_Directory/StudyMult/LEP2/ROOTfiles/cleaned_ALEPH_Data-all.aleph.root";
  
  TFile *f = new TFile(filename.Data());
  //TFile *f = TFile("/mnt/c/Users/Bibek Kumar Pandit/Desktop/Root_Directory/StudyMult/LEP2/ROOTfiles/cleaned_ALEPH_Data-all.aleph.root");
  TTree *t1 = (TTree*)f->Get("t");
  
  
  
  int nParticle;
  int EventNo;
  int RunNo;
  float Energy;
  float tAngle;

  Float_t px[100000];
  Float_t py[100000];
  Float_t pz[100000];
  Float_t pt[100000];
  Float_t eta[100000];
  Float_t theta[100000];
  Float_t phi[100000];
  Float_t mass[100000];
  Float_t charge[10000];
  Float_t pwflag[10000];
  Float_t pid[10000];
  
  t1->SetBranchAddress("nParticle",&nParticle);
  t1->SetBranchAddress("EventNo",&EventNo);
  t1->SetBranchAddress("RunNo",&RunNo);
  t1->SetBranchAddress("Energy",&Energy);
  t1->SetBranchAddress("px",px);
  t1->SetBranchAddress("py",py);
  t1->SetBranchAddress("pz",pz);
  t1->SetBranchAddress("pt",pt);
  t1->SetBranchAddress("mass",mass);
  t1->SetBranchAddress("eta",eta);
  t1->SetBranchAddress("theta",theta);
  t1->SetBranchAddress("phi",phi);
  t1->SetBranchAddress("charge",charge);
  t1->SetBranchAddress("pwflag", pwflag);
  t1->SetBranchAddress("pid",pid);
  

  
  

  TFile *g = new TFile("/mnt/c/Users/Bibek Kumar Pandit/Desktop/Root_Directory/StudyMult/LEP2/ROOTfiles/cleaned_ALEPH_Data2-all.aleph.root", "RECREATE");
  TTree *newtree = t1->CloneTree(0); // Do no copy the data yet
  // DO NOT delete the old tree and close the old file
  //add the branch to the new tree and try to fill it
  newtree->Branch("tAngle", &tAngle, "tAngle/F");
  
  int nentries=newtree->GetEntries();
  //cout<<nentries<<endl;
  //int nentries = 100;
  float x[59874];
   //memset(x, 0, nevent*sizeof(float) );
  cout<<"elo"<<endl; 
    
  thrust(x);
  
  for( int i=0; i < 59874; i++){
     t1->GetEntry(i);
     tAngle = x[i];
     newtree->Fill();
  }
  g->Write();
}
