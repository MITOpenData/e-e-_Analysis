#include <TTree.h>
#include <TVector3.h>
#include "utilities.h"
#include "Selection.h"
#include "../../DataProcessing/include/particleData.h"
#include "../../DataProcessing/include/boostedEvtData.h"
#include "../../DataProcessing/include/eventData.h"
#include "../../DataProcessing/include/jetData.h"

// DataFormat
class TPCNtupleData{
    public:
    static const int nMaxPart = 10000;
    
    Selection s;
    particleData particle;
    jetData jet;
    boostedEvtData boosted;
    eventData event;
    
    bool isBelle;
    bool doThrust;
    bool doWTA;
     
       
    TPCNtupleData(bool ana=0, bool thrust=0, bool wta=0) {
       isBelle = ana;
       doWTA=wta;
       doThrust=thrust;
    }
    
    void setupTPCTree(TTree *tData, TTree *tBoost, TTree *tJet) {
       std::vector<std::string> listData;
       std::vector<std::string> list;
       listData.push_back("nParticle");
       listData.push_back("Energy");
       listData.push_back("pt");
       listData.push_back("pmag");
       listData.push_back("mass");
       listData.push_back("pwflag");
       listData.push_back("eta");
       listData.push_back("theta");
       listData.push_back("phi");
       listData.push_back("d0");
       listData.push_back("z0");
       listData.push_back("ntpc");
       
       if (doThrust) {
          listData.push_back("pt_wrtThr");
          listData.push_back("eta_wrtThr");
          listData.push_back("phi_wrtThr");
          listData.push_back("theta_wrtThr");
       }
       
       std::vector<std::string> listJet;
       listJet.push_back("nref");
       listJet.push_back("jtpt");
       listJet.push_back("jteta");
       
       particle.SetStatusAndAddressRead(tData,listData);
       event.SetStatusAndAddressRead(tData,list);
       if (doWTA) boosted.SetStatusAndAddressRead(tBoost,list);
       jet.SetStatusAndAddressRead(tJet,listJet);
    }	
    
    // Decide if paricle j is a charged hadron    
    bool isChargedHadron(int j)
    {
       if (isBelle)
       {
          // for BELLE analysis
          if (particle.pid[j]!=BELLE_PION&&particle.pid[j]!=BELLE_PROTON&&particle.pid[j]!=BELLE_KAON) return 0; 
       }
        
       else
       {
          // for ALEPH analysis
          if (particle.pwflag[j]!=ALEPH_CHARGED_TRACK) return 0;
       }
       
       return 1;
    }


    Float_t getPt(int j)
    {
     if (doThrust)
     {
         return particle.pt_wrtThr[j];
     } else if (doWTA){
         return boosted.pt[j];
     }
     return particle.pt[j];
    }
        
    // Return Pseudorapidity
    Float_t getEta(int j)
    {
     if (doThrust)
     {
         return particle.eta_wrtThr[j];
     } else if (doWTA)
     {
         return boosted.eta[j];
     }
     return particle.eta[j];
    }
    
    // Return phi angle
    Float_t getPhi(int j)
    {
     if (doThrust)
     {
         return particle.phi_wrtThr[j];
     } else if (doWTA)
     {
         return boosted.phi[j];
     }
     return particle.phi[j];
    }

    // Return theta angle
    Float_t getTheta(int j)
    {
     if (doThrust)
     {
         return particle.theta_wrtThr[j];
     } else if (doWTA) 
     {
        return boosted.theta[j];
     }
     return particle.theta[j];
    }

    // return pmag
    Float_t getPmag(int j)
    {
     if (doWTA) return boosted.pmag[j];
     return particle.pmag[j];
    }
    
    // return TTheta
    Float_t getTTheta()
    {
     if (doWTA) return boosted.WTAAxis_Theta;
     return event.TTheta;
    }

    // return TPhi
    Float_t getTPhi()
    {
     if (doWTA>0) return boosted.WTAAxis_Phi;
     return event.TPhi;
    }
};

