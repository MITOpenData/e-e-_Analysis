#include <TTree.h>
#include <TVector3.h>
#include "utilities.h"

// DataFormat
class TPCNtupleData{
    public:
    static const int nMaxPart = 10000;
    
    Int_t nParticle;
    bool passesWW;
    Float_t missP;
    
    Float_t pt[nMaxPart];
    Float_t eta[nMaxPart];
    Float_t theta[nMaxPart];
    Int_t pid[nMaxPart];
    Float_t phi[nMaxPart];
    Float_t mass[nMaxPart];
    Int_t pwflag[nMaxPart];
    Int_t nTPC[nMaxPart];
    
    Float_t px[nMaxPart];
    Float_t py[nMaxPart];
    Float_t pz[nMaxPart];
    
    Float_t N;
    Float_t N_TP;
    
    // Jet Tree
    Int_t nref;
    Float_t jtpt[nMaxPart];
    Float_t jteta[nMaxPart];
    Float_t jtphi[nMaxPart];
    
    //thrust axis variables
    Float_t pt_wrtThr[nMaxPart];
    Float_t eta_wrtThr[nMaxPart];
    Float_t theta_wrtThr[nMaxPart];
    Float_t phi_wrtThr[nMaxPart];
    Float_t TTheta;
    Float_t TPhi;
    
    // WTA axis variables
    Float_t WTAAxis_Theta;
    Float_t WTAAxis_Phi;
    Float_t boostx;
    Float_t boosty;
    Float_t boostz;
    Float_t boost;
    
    bool isBelle;
    bool isWTA;
    int doThrust;
    Float_t memory;
    TVector3 thrust;
    TVector3 p;
    
    TPCNtupleData(bool ana=0, int thrustAna=0, bool tree=0)
    {
       isBelle = ana;
       doThrust = thrustAna;
       isWTA = tree;
       thrust.SetXYZ(1,0,0);
       p.SetXYZ(1,0,0);
    }
    
    // Decide if paricle j is a charged hadron    
    bool isChargedHadron(int j)
    {
       if (isBelle)
       {
          // for BELLE analysis
          if (pid[j]!=BELLE_PION&&pid[j]!=BELLE_PROTON&&pid[j]!=BELLE_KAON) return 0; 
       }
        
       else
       {
          // for ALEPH analysis
          if (pwflag[j]!=ALEPH_CHARGED_TRACK) return 0;
       }
       
       return 1;
    }

    // Return Transverse Momentum
    Float_t getPt(int j)
    {
     selfCheck();
     if (doThrust>0)
     {
         return pt_wrtThr[j];
         //p.SetXYZ(px[j],py[j],pz[j]);
         //return ptFromThrust(thrust,p);
     }
     return pt[j];
    }
    
    // Return Pseudorapidity
    Float_t getEta(int j)
    {
     selfCheck();
     if (doThrust>0)
     {
         return eta_wrtThr[j];
         //p.SetXYZ(px[j],py[j],pz[j]);
         //return etaFromThrust(thrust,p);
     }
     return eta[j];
    }

    // Return phi angle
    Float_t getPhi(int j)
    {
     selfCheck();
     if (doThrust>0)
     {
         return phi_wrtThr[j];
         //p.SetXYZ(px[j],py[j],pz[j]);
         //return phiFromThrust(thrust,p);
     }
     return phi[j];
    }

    // Return theta angle
    Float_t getTheta(int j)
    {
     selfCheck();
     if (doThrust>0)
     {
         return theta_wrtThr[j];
         //p.SetXYZ(px[j],py[j],pz[j]);
         //return thetaFromThrust(thrust,p);
     }
     return theta[j];
    }

    // 
    void update()
    {
       memory = pt[0]*1000+eta[0];
       if (doThrust==1)
       {
         thrust.SetTheta(TTheta);
         thrust.SetPhi(TPhi);
       }
       else if (doThrust==2)
       {
         thrust.SetTheta(2*atan(exp(-jteta[0])));
         thrust.SetPhi(jtphi[0]);
       }
    }
    
    // check if the class user remember to update the Thrust axis in the analysis code when getting a new entry.
    void selfCheck()
    {
        if (memory!=pt[0]*1000+eta[0]) std::cout <<"Bug in the code!"<<std::endl;
    }
    
    void setTPCTreeStatus(TTree *t1)
    {
      t1->SetBranchStatus("*", 0);
      t1->SetBranchStatus("nParticle", 1);
      t1->SetBranchStatus("pt", 1);
      t1->SetBranchStatus("eta", 1);
      t1->SetBranchStatus("theta", 1);
      t1->SetBranchStatus("pid", 1);
      t1->SetBranchStatus("phi", 1);
      t1->SetBranchStatus("mass", 1);
      t1->SetBranchStatus("pt_wrtThr",1);
      t1->SetBranchStatus("eta_wrtThr",1);
      t1->SetBranchStatus("theta_wrtThr",1);
      t1->SetBranchStatus("phi_wrtThr",1);
      t1->SetBranchStatus("ntpc",1);
      t1->SetBranchStatus("passesWW",1);
      t1->SetBranchStatus("missP",1);
        
      if (!isBelle)
      {
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
    t1->SetBranchAddress("passesWW",&data.passesWW);
    t1->SetBranchAddress("missP",&data.missP);
    t1->SetBranchAddress("pt",data.pt);
    t1->SetBranchAddress("eta",data.eta);
    t1->SetBranchAddress("theta",data.theta);
    t1->SetBranchAddress("pid",data.pid);
    t1->SetBranchAddress("phi",data.phi);
    t1->SetBranchAddress("mass",data.mass);
    t1->SetBranchAddress("pt_wrtThr",data.pt_wrtThr);
    t1->SetBranchAddress("eta_wrtThr",data.eta_wrtThr);
    t1->SetBranchAddress("theta_wrtThr",data.theta_wrtThr);
    t1->SetBranchAddress("phi_wrtThr",data.phi_wrtThr);
    t1->SetBranchAddress("ntpc",data.nTPC);
    
    t2->SetBranchAddress("jtpt",data.jtpt);
    t2->SetBranchAddress("jteta",data.jteta);
    t2->SetBranchAddress("jtphi",data.jtphi);
    t2->SetBranchAddress("nref",&data.nref);

    
    if (!data.isBelle)
    {
        t1->SetBranchAddress("pwflag",data.pwflag);
        t1->SetBranchAddress("px",data.px);
        t1->SetBranchAddress("py",data.py);
        t1->SetBranchAddress("pz",data.pz);
        t1->SetBranchAddress("TTheta", &data.TTheta);
        t1->SetBranchAddress("TPhi", &data.TPhi);
    }
    
    if(data.isWTA)
    {
        t1->SetBranchAddress("WTAAxis_Theta",&data.WTAAxis_Theta);
        t1->SetBranchAddress("WTAAxis_Phi",&data.WTAAxis_Phi);
    }
}

