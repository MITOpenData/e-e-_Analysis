#include <TTree.h>
#include <TVector3.h>
#include "utilities.h"
#include "Selection.h"
// DataFormat
class TPCNtupleData{
    public:
    static const int nMaxPart = 10000;
    
    Selection s;

    Int_t nParticle;
    bool passesWW;
    Float_t missP;
    Float_t pthatWeight;

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
    int sideTree;
    int doThrust;
    Float_t memory;
    TVector3 thrust;
    TVector3 p;
    
    TPCNtupleData(bool ana=0, int thrustAna=0, int tree=0)
    {
       isBelle = ana;
       doThrust = thrustAna;
       sideTree = tree;
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
    
    // main tree, side tree
    void setTPCTreeStatus(TTree *t1, TTree *t2)
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
        
      if(s.doPP)
      {
        t1->SetBranchStatus("pthatWeight",1);
      }
      if (!isBelle)
      {
        t1->SetBranchStatus("pwflag", 1);
        t1->SetBranchStatus("px", 1);
        t1->SetBranchStatus("py", 1);
        t1->SetBranchStatus("pz", 1);
        t1->SetBranchStatus("TTheta", 1);
        t1->SetBranchStatus("TPhi", 1);
      }
    
      if (sideTree == 1)
      {
          t2->SetBranchStatus("WTAAxis_Theta",1);
          t2->SetBranchStatus("WTAAxis_Phi",1);
          t2->SetBranchStatus("pt", 1);
          t2->SetBranchStatus("eta", 1);
          t2->SetBranchStatus("theta", 1);
          t2->SetBranchStatus("phi", 1);
      }
    }
};


// Set the branch addresses
// main tree, side tree, jet tree, data
void setupTPCTree(TTree *t1, TTree *t2, TTree *t3, TPCNtupleData &data)
{
    Selection s;

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
    
    t3->SetBranchAddress("jtpt",data.jtpt);
    t3->SetBranchAddress("jteta",data.jteta);
    t3->SetBranchAddress("jtphi",data.jtphi);
    t3->SetBranchAddress("nref",&data.nref);

    if(s.doPP)
    {
        t1->SetBranchAddress("pthatWeight",&data.pthatWeight);
    }
    
    if (!data.isBelle)
    {
        t1->SetBranchAddress("pwflag",data.pwflag);
        t1->SetBranchAddress("px",data.px);
        t1->SetBranchAddress("py",data.py);
        t1->SetBranchAddress("pz",data.pz);
        t1->SetBranchAddress("TTheta", &data.TTheta);
        t1->SetBranchAddress("TPhi", &data.TPhi);
    }
    
    if(data.sideTree == 1) // WTA Axis
    {
        t2->SetBranchAddress("WTAAxis_Theta",&data.WTAAxis_Theta);
        t2->SetBranchAddress("WTAAxis_Phi",&data.WTAAxis_Phi);
        
        t2->SetBranchAddress("pt_wrtThr",data.pt);
        t2->SetBranchAddress("eta_wrtThr",data.eta);
        t2->SetBranchAddress("theta_wrtThr",data.theta);
        t2->SetBranchAddress("phi_wrtThr",data.phi);
    }
}

