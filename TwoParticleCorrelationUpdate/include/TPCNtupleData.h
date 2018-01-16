#include <TTree.h>
#include <TVector3.h>
#include "utilities.h"

// DataFormat
class TPCNtupleData{
    public:
    
    Int_t nParticle;
    Float_t pt[5000];
    Float_t eta[5000];
    Float_t theta[5000];
    Int_t pid[5000];
    Float_t phi[5000];
    Float_t mass[5000];
    Int_t pwflag[5000];
    
    Float_t px[5000];
    Float_t py[5000];
    Float_t pz[5000];
    Float_t TTheta;
    Float_t TPhi;
    
    Float_t N;
    Float_t N_TP;
    
    // Jet Tree
    Int_t nref;
    Float_t jtpt;
    Float_t jteta[5000];
    Float_t jtphi[5000];
    
    
    bool isBelle;
    int doThrust;
    Float_t memory;
    TVector3 thrust;
    TVector3 p;
    
    TPCNtupleData(bool ana=0, int thrustAna=0){
       isBelle = ana;
       doThrust = thrustAna;
       thrust.SetXYZ(1,0,0);
       p.SetXYZ(1,0,0);
    }
    
    // Decide if paricle j is a charged hadron    
    bool isChargedHadron(int j){
       if (isBelle) {
          // for BELLE analysis
          if (pid[j]!=BELLE_PION&&pid[j]!=BELLE_PROTON&&pid[j]!=BELLE_KAON) return 0; 
       } else {
          // for ALEPH analysis
         if (pwflag[j]!=ALEPH_CHARGED_TRACK) return 0;
       }
       
       return 1;
    }

    // Return Transverse Momentum
    Float_t getPt(int j){
     selfCheck();
     if (doThrust>0){
      p.SetXYZ(px[j],py[j],pz[j]);
      return ptFromThrust(thrust,p);
     }
     return pt[j];
    }
    
    // Return Pseudorapidity
    Float_t getEta(int j){
     selfCheck();
     if (doThrust>0){
      p.SetXYZ(px[j],py[j],pz[j]);
      return etaFromThrust(thrust,p);
     }
     return eta[j];
    }

    // Return phi angle
    Float_t getPhi(int j){
     selfCheck();
     if (doThrust>0){
      p.SetXYZ(px[j],py[j],pz[j]);
      return phiFromThrust(thrust,p);
     }
     return phi[j];
    }

    // Return theta angle
    Float_t getTheta(int j){
     selfCheck();
     if (doThrust>0){
      p.SetXYZ(px[j],py[j],pz[j]);
      return thetaFromThrust(thrust,p);
     }
     return theta[j];
    }

    // 
    void update(){
       memory = pt[0]*1000+eta[0];
       if (doThrust==1) {
         thrust.SetTheta(TTheta);
         thrust.SetPhi(TPhi);
       } else if (doThrust==2) {
         thrust.SetTheta(2*atan(exp(-jteta[0])));
	 thrust.SetPhi(jtphi[0]);
       }
    }
    
    // check if the class user remember to update the Thrust axis in the analysis code when getting a new entry.
    void selfCheck(){
        if (memory!=pt[0]*1000+eta[0]) std::cout <<"Bug in the code!"<<std::endl;
    }
    
    void setTPCTreeStatus(TTree *t1){
      t1->SetBranchStatus("*", 0);
      t1->SetBranchStatus("nParticle", 1);
      t1->SetBranchStatus("pt", 1);
      t1->SetBranchStatus("eta", 1);
      t1->SetBranchStatus("theta", 1);
      t1->SetBranchStatus("pid", 1);
      t1->SetBranchStatus("phi", 1);
      t1->SetBranchStatus("mass", 1);

      if (!isBelle) {
        t1->SetBranchStatus("pwflag", 1);

        t1->SetBranchStatus("px", 1);
        t1->SetBranchStatus("py", 1);
        t1->SetBranchStatus("pz", 1);
        t1->SetBranchStatus("TTheta", 1);
        t1->SetBranchStatus("TPhi", 1);
      }
    }
};


// Set the branch addresses
void setupTPCTree(TTree *t1, TTree *t2, TPCNtupleData &data)
{
    t1->SetBranchAddress("nParticle",&data.nParticle);
    t1->SetBranchAddress("pt",data.pt);
    t1->SetBranchAddress("eta",data.eta);
    t1->SetBranchAddress("theta",data.theta);
    t1->SetBranchAddress("pid",data.pid);
    t1->SetBranchAddress("phi",data.phi);
    t1->SetBranchAddress("mass",data.mass);

    t2->SetBranchAddress("jteta",data.jteta);
    t2->SetBranchAddress("jtphi",data.jtphi);
    t2->SetBranchAddress("nref",&data.nref);
    

    if (!data.isBelle) {
        t1->SetBranchAddress("pwflag",data.pwflag);

        t1->SetBranchAddress("px",data.px);
        t1->SetBranchAddress("py",data.py);
        t1->SetBranchAddress("pz",data.pz);
        t1->SetBranchAddress("TTheta", &data.TTheta);
        t1->SetBranchAddress("TPhi", &data.TPhi);
    }
}

