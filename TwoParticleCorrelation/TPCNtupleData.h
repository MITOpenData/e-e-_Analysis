#include <TTree.h>
#include <TVector3.h>
#include "utilities.h"

// DataFormat
class TPCNtupleData{
    public:
    
    Int_t nParticle;
    Float_t pt[50000];
    Float_t eta[50000];
    Float_t theta[50000];
    Float_t pid[50000];
    Float_t phi[50000];
    Float_t mass[50000];
    Float_t pwflag[50000];
    
    Float_t px[100000];
    Float_t py[100000];
    Float_t pz[100000];
    Float_t TTheta;
    Float_t TPhi;
    bool isBelle;
    bool doThrust;
    Float_t memory;
    TVector3 thrust;
    TVector3 p;
    
    TPCNtupleData(bool ana=0, bool thrustAna=0){
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
     if (doThrust){
      p.SetXYZ(px[j],py[j],pz[j]);
      return ptFromThrust(thrust,p);
     }
     return pt[j];
    }
    
    // Return Pseudorapidity
    Float_t getEta(int j){
     selfCheck();
     if (doThrust){
      p.SetXYZ(px[j],py[j],pz[j]);
      return etaFromThrust(thrust,p);
     }
     return eta[j];
    }

    // Return phi angle
    Float_t getPhi(int j){
     selfCheck();
     if (doThrust){
      p.SetXYZ(px[j],py[j],pz[j]);
      return phiFromThrust(thrust,p);
     }
     return phi[j];
    }

    // Return theta angle
    Float_t getTheta(int j){
     selfCheck();
     return theta[j];
    }

    // 
    void update(){
       memory = pt[0]*1000+eta[0];
       thrust.SetTheta(TTheta);
       thrust.SetPhi(TPhi);
    }
    
    // check if the class user remember to update the Thrust axis in the analysis code when getting a new entry.
    void selfCheck(){
         if (memory!=pt[0]*1000+eta[0]) cout <<"Bug in the code!"<<endl;
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
void setupTPCTree(TTree *t1, TPCNtupleData &data)
{
    t1->SetBranchAddress("nParticle",&data.nParticle);
    t1->SetBranchAddress("pt",data.pt);
    t1->SetBranchAddress("eta",data.eta);
    t1->SetBranchAddress("theta",data.theta);
    t1->SetBranchAddress("pid",data.pid);
    t1->SetBranchAddress("phi",data.phi);
    t1->SetBranchAddress("mass",data.mass);
    t1->SetBranchAddress("pwflag",data.pwflag);
    
    t1->SetBranchAddress("px",data.px);
    t1->SetBranchAddress("py",data.py);
    t1->SetBranchAddress("pz",data.pz);
    t1->SetBranchAddress("TTheta", &data.TTheta);
    t1->SetBranchAddress("TPhi", &data.TPhi);
}
