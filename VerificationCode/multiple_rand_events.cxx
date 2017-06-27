//This is basically doing what rand_angle does, but it can create multiple events at once.

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




double dphi(double alpha, double lower_bound)
{
  if (alpha < lower_bound)
  {
    while (alpha <lower_bound) alpha += PI - lower_bound; //it's upper_bound - lower_bound
  }
  if (alpha > PI) //PI is the upper bound for both theta and phi
  {
    while (alpha > PI) alpha -= PI - lower_bound;
  }
  
  return alpha;
}

/*
Standard Procedure of Generating Random Numbers:
  random_device r;
  seed_seq s{r(), r(), ...} or we don't need this line at all
  mt19937 e(s) or mt19937 e(r()) if the above line not present
  TYPE_OF_DISTRIBUTION <> dist
  random_number = dist(e)
*/





double binary_generator() //This funciton can be used to generate binary values, 0 and 1, at random
{
  random_device r;  //Will be used to obtain a seed for the random number engine
  seed_seq s{r(), r(), r(), r(),r(),r(),r(),r(),r(),r()}; 
  mt19937 e(s);//Standard mersenne_twister_engine seeded with rd()
  bernoulli_distribution d;

  return (double)d(e);
}






double random_angle(double lower_bound, double upper_bound, int nparticles, double sigma, float angle_array[])
{
  random_device r;  
  mt19937 gen(r());
  uniform_real_distribution<> dis(lower_bound, upper_bound);//theta and phi are chosen from uniform distribution
  double angle = dis(gen);
  
  
  random_device r2;
  mt19937 generator(r2());
  // distribution(mean,std dev)
  normal_distribution<double> distribution(angle, sigma);//normal distribution centered on theta/phi
  
  
  for (int i = 0; i < nparticles; i = i+2) //we do i=i+2 because index i=i+1 is not random, but complimentary. we set that part in next function
  {
    angle_array[i] = (float)dphi(distribution(generator), lower_bound);
  }
  return angle;
}




void event_generator(particleData& pData) //This funcition calls on other function and generates a virtual event
{  
  int nparticles = pData.nParticle;
  double sigma_phi = 0.1;
  double sigma_theta = 0.1;
  
  double thrust_phi = random_angle(-PI, PI, nparticles, sigma_phi, pData.phi); //thrust_phi is the phi of thrust vector
  double thrust_theta = random_angle(0., PI, nparticles, sigma_theta, pData.theta); // thrust_theta is the theta of thrust vector
  
  for (int i = 1; i < nparticles; i = i+2) //Note here that i starts at 1. Fixing the complimentary part.
  { 
    //rotating vector by 180 degrees in 3D corresponds to changing: phi --> phi - pi, theta --> pi - theta, and we do that below
    pData.phi[i] = (float)dphi((double)pData.phi[i-1] - PI, -PI); 
    pData.theta[i] = (float)dphi(PI - (double)pData.theta[i-1], 0.);
  }
  
  TString filename;
  //filename="/data/flowex/Datasamples/Belle/myoutput_isBelle1_minMult40.root";
  filename = "cleaned_ALEPH_Data-all.aleph.root";
  
  
  
  TFile *f = new TFile(filename.Data());
  TTree *t1 = (TTree*)f->Get("t");
  
  
  t1->Draw("pmag>>h");
  TH1F *h = (TH1F*)gDirectory->Get("h");
  
  for (int i = 0; i < nparticles; i = i+2)
  {
    float momentum = h->GetRandom();
    pData.pmag[i] = momentum;
    if (i+1<nparticles) pData.pmag[i+1] = momentum;
  }
  

}

void scan()
{
  TFile *hf = new TFile("new_mult_data.root", "RECREATE" );
	TTree *tout = new TTree("t","");
	
  
	
  particleData pData;
	tout->Branch("nParticle", &pData.nParticle,"nParticle/I");
	tout->Branch("EventNo", &pData.EventNo,"EventNo/I");
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
  
  const int no_of_events = 5; //As many events you want in your data
  int particle_array[no_of_events] = {30, 40, 26, 23, 33}; //Array of no. of particles in each event
  
  for (int i =0; i<no_of_events; i++)
  {
    pData.nParticle = particle_array[i];
    pData.EventNo = i;
    event_generator(pData);
      
    
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
  }
  hf->Write();
  
  cout<<"Complete"<<endl; //Marks everything went well
}
