//
//  ridge_eventSelection.h
//  
//
//  Created by Anthony Badea on 1/22/18.
//

#ifndef SELECTION
#define SELECTION

//c and c++ dependencies
#include <sstream>
#include <string>
#include <iostream>
#include <stdlib.h>

//root dependencies
#include <TMath.h>
#define PI 3.14159265358979

inline float jtp(float pt, float eta){ return pt*TMath::CosH(eta);}

class Selection
{
    public:
    
        // event cuts
        Float_t missPCut = 20;
        bool domissPCut = false;
        bool doWW = false;
        static const Int_t nMultBins = 3;
        Int_t multBinsLow[nMultBins]  = {0 , 20, 30};
        Int_t multBinsHigh[nMultBins] = {20, 30, 999};
    
        // particle cuts
        Float_t ptMin = 0.4;
        Float_t ptMax = 100.0;
        Float_t etaCut = 1.8;
        Int_t nTPCMin = 3; // completely inclusive, maybe tighten this
        Int_t nTPCMax = 23 ; 
    
        // jet cuts
        Float_t AjCut = 0.1;
        Float_t thirdJetCut = 0.03;
        bool doAjCut = false;
        bool do3jetEvtCut = false;
        bool doMixedJetCut = false;
    
        // plotting
        Float_t etaPlotRange = 1.6;//this gets multiplied by 2
        Float_t dEtaBins = 20;//keep even
        Float_t dPhiBins = 20;//keep factor of 4
        Float_t differential = (etaPlotRange/(float)dEtaBins)*(2*TMath::Pi()/(float)dPhiBins);
        Float_t fillAj = 0.0; // used for plotting h_Aj
    
        // other
        bool doThrust = false;
        bool doBelle = false;
        bool doTheta = false;
        //Int_t maxevt = ;
        Int_t nbin = 20;
        Int_t bkgrd_runs = 1; // decide if we should actually be using this
    
                                                                
        // BELLE Particle Definition
        enum SIMPLEPID {BELLE_PHOTON, BELLE_ELECTRON, BELLE_PION, BELLE_MUON, BELLE_KAON, BELLE_PROTON};
        // ALEPH Particle Flow Classification
        enum SIMPLEPWFLAG {ALEPH_CHARGED_TRACK, ALEPH_CHARGED_LEPTONS1, ALEPH_CHARGED_LEPTONS2, ALEPH_V0, ALEPH_PHOTON, ALEPH_NEUTRAL_HADRON};
    
        Selection();
        int ridge_trackSelection(Float_t pt, Float_t eta, Int_t nTPC, Int_t pwflag);
        int ridge_eventSelection(bool passesWW, Int_t nParticle, Float_t missP,Float_t pt[], Float_t eta[], Int_t nTPC[], Int_t pwflag[], Int_t nref, Float_t jtpt[], Float_t jteta[]);
        bool mixedEvent(Int_t nParticle, Int_t nParticle_mix, Float_t jteta, Float_t jteta_mix);
        int histNum(Int_t N);
    
    private:
};

Selection::Selection()
{
    std::cout << "Getting settings.." << std::endl;
    return;
}

// return 1 if the track passes the selection
// otherwise return 0
int Selection::ridge_trackSelection
(
 // particle variables
 Float_t pt,
 Float_t eta,
 Int_t nTPC,
 Int_t pwflag
 )
{
    if (pwflag != ALEPH_CHARGED_TRACK) return 0;
    if (pt < ptMin || pt > ptMax) return 0;
    if (nTPC <= nTPCMin || nTPC >= nTPCMax) return 0;
    if (TMath::Abs(eta) >= etaCut) return 0;
    
    return 1;
}

// return -1 if does not pass event selection
// return multiplicity
int Selection::ridge_eventSelection
(
    // QCD paper
    bool passesWW,
 
    // particle variables
    Int_t nParticle,
    Float_t missP,
    Float_t pt[],
    Float_t eta[],
    Int_t nTPC[],
    Int_t pwflag[],
 
    // jet variables
    Int_t nref,
    Float_t jtpt[],
    Float_t jteta[]

)
{
    
    ///////// QCD Paper Selection /////////
    if (doWW && !passesWW) return -1;
    
    ///////// Missing Momentum /////////
    if ( domissPCut && missP > missPCut) return -1;
    
    ///////// 3-Jet /////////
    Float_t j12 = jtp(jtpt[0],jteta[0])+jtp(jtpt[1],jteta[1]);
    fillAj = TMath::Abs(jtp(jtpt[0],jteta[0])-jtp(jtpt[1],jteta[1])) / j12;
    if (nref>=2)
    {
        
        // require 2 jets to be pretty equally balanced in momentum
        if(doAjCut && fillAj > AjCut) return -1;
        // require 3rd jet has low momentum relative to first two jets (take average momentum of jet 1 and 2)
        if(do3jetEvtCut && nref > 2 && (2*jtp(jtpt[2],jteta[2]) / j12) > thirdJetCut ) return -1;
    }
    
    ///////// CHARGED PARTICLES /////////
    Int_t N = 0;
    for (Int_t j=0;j<nParticle;j++)
    {
        N += ridge_trackSelection(pt[j],eta[j],nTPC[j],pwflag[j]);
    }
    if (N <= 1) return -1;
    // this is to ensure that we are able to find a mixed event
    if( N > 1000) return -1;
    
    return N;
}

// return true if the event passes the criteria for a mixed event
// Currently matching the total multiplicity in the event and leading jet eta
bool Selection::mixedEvent(Int_t nParticle, Int_t nParticle_mix, Float_t jteta, Float_t jteta_mix)
{
    // Original matched event criteria
    // (fabs(mix.nParticle-data.nParticle)>4&&data.nParticle<1000&&fabs(mix.jteta[0]-data.jteta[0])>0.2)
    if (TMath::Abs(nParticle - nParticle_mix) <= 4) return false;
    if (doMixedJetCut && TMath::Abs(jteta - jteta_mix) <= 0.2) return false;
    
    return true;
}
// return histogram number of the histogram corresponding to the multiplicity rather than the multiplicity itself
int Selection::histNum(Int_t N)
{
    // inclusive on low end to help statistics in higher multiplicity bins
    if (N >= multBinsLow[0] && N < multBinsHigh[0]) return 0;
    if (N >= multBinsLow[1] && N < multBinsHigh[1]) return 1;
    if (N >= multBinsLow[2] && N < multBinsHigh[2]) return 2;
    
    return -1;
}
#endif /* Selection.h */
