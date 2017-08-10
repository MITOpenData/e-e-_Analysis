//This macro here creates a virtual event. It starts by picking a random thrust
//direction and sampling from the normal distribution surrounding the thrust angle.
//These samples represent the trajectories the particles. Then it assigns these 
//particles a momentum from the real data. And  finally create the n-tuple of the event.

#include <random>
#include <iostream>
#include <math.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#define PI 3.1415926

using namespace std;

double dphi(double alpha, double lower_bound)
{
  if (alpha < lower_bound)
  {
    while (alpha <lower_bound) alpha += PI - lower_bound;
  }
  if (alpha > PI)
  {
    while (alpha > PI) alpha -= PI - lower_bound;
  }
  
  return alpha;
}



//can we put this in and get rid of the second line?
random_device r;  //Will be used to obtain a seed for the random number engine
mt19937 gen(r());//Standard mersenne_twister_engine seeded with rd()

double binary_generator()
{
  seed_seq s{r(), r(), r(), r(), r(), r(), r(), r(), r(), r()}; 
  mt19937 e(s);
  bernoulli_distribution d;

  binary = -1 + 2*d(e);
  return (double)binary;
}

class particleData
{
	public:
		int nParticle;
		//int EventNo;
		//int RunNo;
		//float Energy;

		Float_t px[100000];
		Float_t py[100000];
		Float_t pz[100000];
		Float_t pt[100000];
    Float_t pmag[100000];
		Float_t eta[100000];
    Float_t theta[100000];
		Float_t phi[100000];
		//Float_t mass[100000];
		//Float_t charge[100000];
		//Float_t pwflag[10000];
		//Float_t pid[10000];
};

double random_angle(double lower_bound, double upper_bound, int nparticles, double sigma, float angle_array[], int istheta)
{
  

  uniform_real_distribution<> dis(lower_bound, upper_bound);//theta and phi are chosen from uniform distribution
  double angle = dis(gen);
  
  random_device generator;
  // distribution(mean,std dev)
  normal_distribution<double> distribution(angle, sigma);//normal distribution centered on theta/phi
  
  for (int i = 0; i < nparticles; i++)
  {
    if (istheta == 1) angle_array[i] = (float)dphi((1-binary_generator())*PI + (-1. + 2.*binary_generator())*distribution(generator), lower_bound);
    else angle_array[i] = (float)dphi(distribution(generator) - (1. - binary_generator())*PI, lower_bound);
  }
  return angle;

}
    

void goofy(particleData& pData)
{
  
  
  int nparticles = pData.nParticle;
  double sigma_phi = 0.3;
  double sigma_theta = 0.3;
  
  double thrust_phi = random_angle(-PI, PI, nparticles, sigma_phi, pData.phi, 0);
  double thrust_theta = random_angle(0., PI, nparticles, sigma_theta, pData.theta, 1);
  
  
  TString filename;
  //filename="/data/flowex/Datasamples/Belle/myoutput_isBelle1_minMult40.root";
  filename = "cleaned_ALEPH_Data-all.aleph.root";
  
  
  
  TFile *f = new TFile(filename.Data());
  TTree *t1 = (TTree*)f->Get("t");
  
  
  t1->Draw("pmag>>h");
  TH1F *h = (TH1F*)gDirectory->Get("h");
  
  for (int i = 0; i < nparticles; i++)
  {
    pData.pmag[i] = h->GetRandom();
    cout<<pData.pmag[i]<<endl;
  }
  

}
//test() function was used to plot the created particles on theta-phi space to see how well distributed they were
/*
void test()
{
  particleData pData;
  goofy(pData);
  nparticles = 1000;
  TH2F *h_theta_phi = new  TH2F ("theta-phi of generated particles", "#theta-#phi of generated particles", 100, 0, PI, 100, -PI,PI);
  
  for (int i=0; i<nparticles; i++)
  {
    h_theta_phi->SetMarkerStyle(20);
    h_theta_phi->Fill(pData.theta[i], pData.phi[i]);
  }
  
  TCanvas *c1 = new TCanvas("theta-phi of generated particles", "theta_phi of generated particles", 600,  600);
  
  c1->cd();
  h_theta_phi->Draw();
  
}
*/

void scan()
{
  TFile *hf = new TFile("new_data8.root", "RECREATE" );
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
  
  pData.nParticle = 40;//this can be adjusted as per wish, it can also be made random
  
  goofy(pData);
    
  
  for (int counterParticles = 0; counterParticles < pData.nParticle; counterParticles++)
  {
    
    TVector3 v(1, 1, 1); 
    //we use v(1,1,1) and not just v because v by default is (0,0,0) and there are problems stretching it
    v.SetMag(pData.pmag[counterParticles]);
    v.SetPhi(pData.phi[counterParticles]);
    v.SetTheta(pData.theta[counterParticles]);
    
    pData.px[counterParticles]= v.X();
    pData.py[counterParticles]=v.Y();
    pData.pz[counterParticles]=v.Z();
        
    pData.pt[counterParticles]=v.Pt();
    pData.eta[counterParticles]=v.PseudoRapidity();
    
	}
  //
  tout->Fill();
  hf->Write();
  cout<<"complete"<<endl;
  
}
