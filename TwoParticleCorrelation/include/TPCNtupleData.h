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
    particleData genParticle;
    jetData jet;
    boostedEvtData boosted;
    eventData event;
    
    bool isBelle;
    bool doThrust;
    bool doWTA;
    bool doPerp;
    bool doTGen;
    bool isMix;
     
       
    TPCNtupleData(bool ana=0, bool thrust=0, bool wta=0, bool perp=0, bool mix=0, bool tgen=0) {
       isBelle = ana;
       doWTA=wta;
       doThrust=thrust;
       doPerp=perp;
       doTGen=tgen;
       isMix=0;
    }
    
    void setupTPCTree(TTree *tData, TTree *tBoost, TTree *tJet, TTree *tGen) {
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
       //listData.push_back("passesLEP1TwoPC");
       //listData.push_back("nChargedHadronsHP");
       //listData.push_back("d0");
       //listData.push_back("z0");
       //listData.push_back("ntpc");
       listData.push_back("highPurity");
       listData.push_back("particleWeight");
       
       if (doThrust) {
          if (doPerp) {
             listData.push_back("pt_wrtThrPerp");
             listData.push_back("eta_wrtThrPerp");
             listData.push_back("phi_wrtThrPerp");
             listData.push_back("theta_wrtThrPerp");
	  } else {
             listData.push_back("pt_wrtThr");
             listData.push_back("eta_wrtThr");
             listData.push_back("phi_wrtThr");
             listData.push_back("theta_wrtThr");
	  }   
       }
       
       std::vector<std::string> listJet;
       listJet.push_back("nref");
       listJet.push_back("jtpt");
       listJet.push_back("jteta");
       listJet.push_back("jtphi");
       
       particle.SetStatusAndAddressRead(tData,listData);
       event.SetStatusAndAddressRead(tData,list);
       boosted.SetStatusAndAddressRead(tBoost,list);
       jet.SetStatusAndAddressRead(tJet,listJet);
       if (doTGen) genParticle.SetStatusAndAddressRead(tGen,listData);
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
         if (doPerp) {
	    if (doTGen) return genParticle.pt_wrtThrPerp[j];
	    else return particle.pt_wrtThrPerp[j];
	 } else {
	    if (doTGen) return genParticle.pt_wrtThr[j];
	    else return particle.pt_wrtThr[j];
	 }
     } else if (doWTA){
         if (doPerp) {
	    return boosted.pt_Perp[j];
	 } else {
            return boosted.pt[j];
	 }
     }
     if (doTGen) return genParticle.pt[j];
     else return particle.pt[j];
    }
        
    // Return Pseudorapidity
    Float_t getEta(int j)
    {
     if (doThrust)
     {
         if (doPerp) {
            if (doTGen) return genParticle.eta_wrtThrPerp[j];
            else return particle.eta_wrtThrPerp[j];
         } else {
            if (doTGen) return genParticle.eta_wrtThr[j];
            else return particle.eta_wrtThr[j];
         }
     } else if (doWTA)
     {
         if (doPerp) {
            return boosted.eta_Perp[j];
	 } else {   
            return boosted.eta[j];
	 }   
     }
     if (doTGen) return genParticle.eta[j];
     else return particle.eta[j];
    }
    
    // Return phi angle
    Float_t getPhi(int j)
    {
     if (doThrust)
     {
         if (doPerp) {
            if (doTGen) return genParticle.phi_wrtThrPerp[j];
            else return particle.phi_wrtThrPerp[j];
         } else {
            if (doTGen) return genParticle.phi_wrtThr[j];
            else return particle.phi_wrtThr[j];
         }
     } else if (doWTA)
     {
         if (doPerp) {
            return boosted.phi_Perp[j];
	 } else { 
	    return boosted.phi[j];
	 }     
     }
     if (doTGen) return genParticle.phi[j];
     else return particle.phi[j];
    }

    // Return theta angle
    Float_t getTheta(int j)
    {
     if (doThrust)
     {
         if (doPerp) {
            if (doTGen) return genParticle.theta_wrtThrPerp[j];
            else return particle.theta_wrtThrPerp[j];
	 } else {
	    if (doTGen) return genParticle.theta_wrtThr[j];
            else return particle.theta_wrtThr[j];
	 }    
     } else if (doWTA) 
     {
         if (doPerp) {
            return boosted.theta_Perp[j];
	 } else {   
            return boosted.theta[j];
         }
     }
     if (doTGen) return genParticle.theta[j];
     else return particle.theta[j];
    }

    // return pmag
    Float_t getPmag(int j)
    {
     if (doWTA) {
        if (doPerp) {
	   return boosted.pmag_Perp[j];
        } else {
	   return boosted.pmag[j];
	}
     }	 
     if (doTGen) return genParticle.pmag[j];
     else return particle.pmag[j];
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

