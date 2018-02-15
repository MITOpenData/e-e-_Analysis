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
    
        Int_t experiment = 0; // 0 ALEPH , 1 DELPHI, 2 BELLE, 3 CMS pp
    
        // Initial Setup
        int side_tree = 0;  // 0 t, 1 BoostedWTAR8Evt,
        int jttree = 0; // 0 ak4ESchemeJetTree, 1 ak4WTAmodpSchemeJetTree, 2 ak8ESchemeJetTree, 3 ak8WTAmodpSchemeJetTree
        bool doThrust = true;  // used for determining thrust/beam angles to use in filling histograms
        bool donTrkThrust = false; // used for calculation of nTrk (true = Beam)
        bool doBelle = false;
        bool doTheta = false;
        bool doPP = false;
        bool doParallel = true;
        // for testing purposes
        bool doOneEvent = false;     int numEvents = 50000;
    
        // event cuts
        Float_t missPCut = 20;
        bool domissPCut = false;
        bool doWW = false;
        static const Int_t nMultBins = 3;
        Int_t multBinsLow[nMultBins]  = {0 , 20, 30};
        Int_t multBinsHigh[nMultBins] = {20, 30, 999};
        bool doMixedMultCut = false;
    
        // particle cuts
        Float_t ptMin = 0.4;    Float_t ptMin_wrtThr = 0.4;
        Float_t ptMax = 100.0;  Float_t ptMax_wrtThr = 100.0;
        Float_t etaCut = 1.8;   Float_t etaCut_wrtThr = 5.0;
        Float_t thetaMin = 20;
        Float_t thetaMax = 160;
        Float_t pMin = 0.2;
        Float_t dxyCut = 3;
        Float_t dzCut = 5;
        Int_t nTPCMin = 4; // completely inclusive, maybe tighten this
        Int_t nTPCMax = 100;
        bool doNTPC = 0;
    
        // jet cuts
        Float_t AjCut = 0.1;
        Float_t thirdJetCut = 0.03;
        bool doAjCut = false;
        bool do3jetEvtCut = false;
        bool doMixedJetCut = false;
    
        // plotting
        Float_t etaPlotRange = 3.2;     Float_t etaPlotRange_wrtThr = 6.0; //this gets multiplied by 2
        Float_t dEtaBins = 20;//keep even
        Float_t dPhiBins = 20;//keep factor of 4
        Float_t differential = (2*etaPlotRange/(float)dEtaBins)*(2*TMath::Pi()/(float)dPhiBins);
        Float_t fillAj = 0.0; // used for plotting h_Aj
        Float_t dEtaRangeToIntegrate[2] = {2.0,3.6}; //try to make this correspond with bin edges based on above parameters
    
        //Int_t maxevt = ;
        Int_t nbin = 20;
        Int_t bkgrd_runs = 1; // decide if we should actually be using this
    
        //kinematics (if trig != assoc cuts, make sure doExcludeNTrigLT2 is set to false)
        float trigPt[2] = {0.4,100};
        float assocPt[2] = {0.4,100};
        float nTrkPt[2] = {0.4,100};
    
        // BELLE Particle Definition
        enum SIMPLEPID {BELLE_PHOTON, BELLE_ELECTRON, BELLE_PION, BELLE_MUON, BELLE_KAON, BELLE_PROTON};
        // ALEPH Particle Flow Classification
        enum SIMPLEPWFLAG {ALEPH_CHARGED_TRACK, ALEPH_CHARGED_LEPTONS1, ALEPH_CHARGED_LEPTONS2, ALEPH_V0, ALEPH_PHOTON, ALEPH_NEUTRAL_HADRON};
    
        Selection();
        Float_t getEtaPlotRange();
        Float_t getDifferential();
        int ridge_trackSelection(Float_t pt, Float_t eta, Int_t nTPC, Int_t pwflag, bool isThrust);
        int ridge_eventSelection(bool passesWW, Int_t nParticle, Float_t missP,Float_t pt[], Float_t eta[], Int_t nTPC[], Int_t pwflag[], Int_t nref, Float_t jtpt[], Float_t jteta[], bool isThrust);
        bool isMixedEvent(Int_t nParticle, Int_t nParticle_mix, Float_t jteta, Float_t jteta_mix);
        int histNum(Int_t N);
    
    private:
};

Selection::Selection()
{
    std::cout << "Getting settings.." << std::endl;
    if(doPP)
    {
        multBinsLow[2] = 110;
        multBinsHigh[1] = 110;
        ptMax = 3.0;
        side_tree = 0;  // 0 t, 1 BoostedWTAR8Evt,
        jttree = 0; // 0 ak4ESchemeJetTree, 1 ak4WTAmodpSchemeJetTree, 2 ak8ESchemeJetTree, 3 ak8WTAmodpSchemeJetTree
        doThrust = false;  // used for determining thrust/beam angles to use in filling histograms
        donTrkThrust = false; // used for calculation of nTrk (true = Beam)
    }
    return;
}

Float_t Selection::getEtaPlotRange()
{
    if(doThrust) return etaPlotRange_wrtThr;
    return etaPlotRange;
}

Float_t Selection::getDifferential()
{
    if(doThrust) return (2*etaPlotRange_wrtThr/(float)dEtaBins)*(2*TMath::Pi()/(float)dPhiBins);
    return (2*etaPlotRange/(float)dEtaBins)*(2*TMath::Pi()/(float)dPhiBins);
}
// return 1 if the track passes the selection
// otherwise return 0
int Selection::ridge_trackSelection
(
 // particle variables
 Float_t pt,
 Float_t eta,
 Int_t nTPC,
 Int_t pwflag,
 
 // thrust vs beam
 bool isThrust
 )
{
    
    if (pwflag != ALEPH_CHARGED_TRACK) return 0;
    
    Float_t ptMin_temp = 0;
    Float_t ptMax_temp = 0;
    Float_t etaCut_temp = 0;
    
    if(!isThrust)
    {
        ptMin_temp = ptMin;
        ptMax_temp = ptMax;
        etaCut_temp = etaCut;
    }
    else if (isThrust)
    {
        ptMin_temp = ptMin_wrtThr;
        ptMax_temp = ptMax_wrtThr;
        etaCut_temp = etaCut_wrtThr;
    }
    if(doPP) ptMax_temp = 3.0;
    if (pt < ptMin_temp || pt > ptMax_temp) return 0;
    if (doNTPC && (nTPC <= nTPCMin || nTPC >= nTPCMax) ) return 0;
    if (TMath::Abs(eta) >= etaCut_temp) return 0;
    
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
    Float_t jteta[],
 
    // thrust vs beam
    bool isThrust

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
        N += ridge_trackSelection(pt[j],eta[j],nTPC[j],pwflag[j],isThrust);
    }
    if (N <= 1) return -1;
    // this is to ensure that we are able to find a mixed event
    if( N > 1000) return -1;
    
    return N;
}

// return true if the event passes the criteria for a mixed event
// Currently matching the total multiplicity in the event and leading jet eta
bool Selection::isMixedEvent(Int_t nParticle, Int_t nParticle_mix, Float_t jteta, Float_t jteta_mix)
{
    // Original matched event criteria
    // (fabs(mix.nParticle-data.nParticle)>4&&data.nParticle<1000&&fabs(mix.jteta[0]-data.jteta[0])>0.2)
    
    if (doMixedMultCut && TMath::Abs(nParticle - nParticle_mix) > 4) return false;
    if (doMixedJetCut && TMath::Abs(jteta - jteta_mix) > 0.2) return false;
    
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
