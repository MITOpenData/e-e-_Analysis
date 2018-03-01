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
#include <math.h>

//root dependencies
#include <TMath.h>
#define PI 3.14159265358979

inline float jtp(float pt, float eta){ return pt*TMath::CosH(eta);}

class Selection
{
    public:
    
        // Initial Setup
        bool doParallel = true;
        bool doOneEvent = false;     int numEvents = 50000;
        bool doBelle = false;
        Int_t experiment = 0; // 0 ALEPH , 1 DELPHI, 2 BELLE, 3 CMS pp
        Int_t jttree = 0; // 0 ak4ESchemeJetTree, 1 ak4WTAmodpSchemeJetTree, 2 ak8ESchemeJetTree, 3 ak8WTAmodpSchemeJetTree
        Int_t nbin = 20;
        Int_t bkgrd_runs = 1;
        enum SIMPLEPID {BELLE_PHOTON, BELLE_ELECTRON, BELLE_PION, BELLE_MUON, BELLE_KAON, BELLE_PROTON};    // BELLE Particle Definition
        enum SIMPLEPWFLAG {ALEPH_CHARGED_TRACK, ALEPH_CHARGED_LEPTONS1, ALEPH_CHARGED_LEPTONS2, ALEPH_V0, ALEPH_PHOTON, ALEPH_NEUTRAL_HADRON};  // ALEPH Particle Flow Classification
        bool doGen = true; // turn of event and track selections
        
        /* Detector Specific Cuts */
        
        // From 1990 "Properties of Hadronic Events in e+e- Annihilation at sqrt(s) = 91 GeV" ALEPH Collaboration paper
            // Track Cuts
        bool doNTPC = 0;    Int_t nTPCMin = 4;  Int_t nTPCMax = 999;
        bool doTheta = true;   Float_t thetaMin = 20*TMath::Pi()/180.;  Float_t thetaMax = 160*TMath::Pi()/180.;  // measured in radians
        bool doP = true;   Float_t pMin = 0.2; // measured in GeV
        bool dod0 = true;     Float_t d0Cut = 3; // measured in cm
        bool doz0 = true;   Float_t z0Cut = 5;  // measured in cm
            // Event Cuts
        bool doE = true; Float_t TotalChrgEnergyMin = 15; // measured in GeV
        Float_t nTrkMin = 5;
        bool doSTheta = true; Float_t SThetaMin = 35;  Float_t SThetaMax = 145;  // measured in degrees
        // end of 1990 cuts
        bool domissPCut = true;  Float_t missPCut = 20;  // measured in GeV
        /* Frame Dependent Cuts */

        bool doPt = true;
        
        bool doWW = true;  // impose WW cut
        bool doMixedMultCut = false; // cut on |nTrk - nTrkMix|
        // jet cuts
        bool doAjCut = false;   Float_t AjCut = 0.1;
        bool do3jetEvtCut = false;  Float_t thirdJetCut = 0.03;
        bool doMixedJetCut = false;

        bool doMixedThrustAxisCut = false; Float_t maxTTheta_TPhi_Cut = 0.3;
        Float_t fillAj = 0.0; // used for plotting h_Aj
        Int_t nptBins = 0;
        Int_t netaBins = 0;
        Float_t etaPlotRange;

        // Beam Axis
        static const Int_t nptBins_wrtBeam = 2;
        Float_t ptBinsLow_wrtBeam[nptBins_wrtBeam]  = {1.0,0.4};  // measured in GeV
        Float_t ptBinsHigh_wrtBeam[nptBins_wrtBeam] = {3.0,100};
        static const Int_t netaBins_wrtBeam = 2;
        Float_t etaBinsLow_wrtBeam[netaBins_wrtBeam]  = {1.6, 1.8};
        Float_t missPCut_wrtBeam = 20;
        Float_t etaPlotRange_wrtBeam = 3.2;

        Float_t ptMin = 0.4; Float_t ptMax = 100.0;
        Float_t etaCut = 1.8;
    
        // Thrust Axis
        bool doThrust = false;  bool donTrkThrust = false; // false = use beam axis for nTrk calculation
        static const Int_t nptBins_wrtThr = 1;
        Float_t ptBinsLow_wrtThr[nptBins_wrtThr]  = {1.0,0.4};  // measured in GeV {1.0,0.4}
        Float_t ptBinsHigh_wrtThr[nptBins_wrtThr] = {3.0,100}; // {3.0,100.0}
        static const Int_t netaBins_wrtThr = 1;
        Float_t etaBinsLow_wrtThr[netaBins_wrtThr]  = {1.8}; //{4.5,5.0}
        Float_t missPCut_wrtThr = 20;
        Float_t etaPlotRange_wrtThr = 6.0;

        Float_t ptMin_wrtThr = 0.4; Float_t ptMax_wrtThr = 100.0; // measured in GeV
        Float_t etaCut_wrtThr = 5.0;

        // WTA Axis
        bool doWTA = false;
        static const Int_t nptBins_wrtWTA = 2;
        Float_t ptBinsLow_wrtWTA[nptBins_wrtWTA]  = {1.0,0.4};  // measured in GeV
        Float_t ptBinsHigh_wrtWTA[nptBins_wrtWTA] = {3.0,100};
        static const Int_t netaBins_wrtWTA = 3;
        Float_t etaBinsLow_wrtWTA[netaBins_wrtWTA]  = {10.0,20.0,40.0};
        Float_t missPCut_wrtWTA = 20;
        Float_t etaPlotRange_wrtWTA = 6.0;

        /* Independent Cuts */
        static const Int_t nEnergyBins = 2;
        Float_t energyBinsLow[nEnergyBins] = {0,100};
        Float_t energyBinsHigh[nEnergyBins] = {100,999};

        static const Int_t nMultBins = 3;
        Int_t multBinsLow[nMultBins]  = {0 , 20, 30};
        Int_t multBinsHigh[nMultBins] = {20, 30, 999};

        /* Plotting */
        Float_t dEtaBins = 20; //keep even
        Float_t dPhiBins = 20; //keep factor of 4
        Float_t dEtaRangeToIntegrate[2] = {2.0,3.6}; //used for Austin's implementation of getLongRangeYield(.) in utilities.h
    
        //kinematics (if trig != assoc cuts, make sure doExcludeNTrigLT2 is set to false)
        //float trigPt[2] = {0.4,100};
        //float assocPt[2] = {0.4,100};
        //float nTrkPt[2] = {0.4,100};
    
        /* pp Cross-Check */
        bool doPP = false;
        static const Int_t nptBins_pp = 2;
        Float_t ptBinsLow_pp[nptBins_wrtWTA]  = {0.4, 1.0};  // measured in GeV
        Float_t ptBinsHigh_pp[nptBins_wrtWTA] = {100, 3.0};
        static const Int_t netaBins_pp = 2;
        Float_t etaBinsLow_pp[netaBins_wrtWTA]  = {4.5,5.0};
        Float_t missPCut_pp = 20;
        Float_t etaPlotRange_pp = 6.0;
    
        Selection();
        Float_t getEtaPlotRange();
        Float_t getDifferential();
        int ridge_trackSelection(Int_t nTPC, Float_t theta, Float_t p, Float_t d0, Float_t z0, Int_t pwflag);
        int ridge_eventSelection(bool passesWW, Float_t missP, Int_t nParticle, Int_t nref, Float_t jtpt[], Float_t jteta[], Float_t STheta, Float_t mass[], Int_t nTPC[], Float_t theta[], Float_t pmag[], Float_t d0[], Float_t z0[], Int_t pwflag[]);
        bool isMixedEvent(Int_t nParticle, Int_t nParticle_mix, Float_t jteta, Float_t jteta_mix, Float_t TTheta, Float_t TTheta_mix, Float_t TPhi, Float_t TPhi_mix);
        int histNtrk(Int_t N);
        int histPt(Float_t pt);
        int histEta(Float_t eta);
        int histEnergy(Float_t Energy);
    
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
        doWTA = false;  // 0 t, 1 BoostedWTAR8Evt,
        jttree = 0; // 0 ak4ESchemeJetTree, 1 ak4WTAmodpSchemeJetTree, 2 ak8ESchemeJetTree, 3 ak8WTAmodpSchemeJetTree
        doThrust = false;  // used for determining thrust/beam angles to use in filling histograms
        donTrkThrust = false; // used for calculation of nTrk (true = Beam)
    }

    if(doThrust) { nptBins = nptBins_wrtThr; netaBins = netaBins_wrtThr; etaPlotRange = etaPlotRange_wrtThr;}
    else if(doWTA) { nptBins = nptBins_wrtWTA; netaBins = netaBins_wrtWTA; etaPlotRange = etaPlotRange_wrtWTA;}
    else { nptBins = nptBins_wrtBeam; netaBins = netaBins_wrtBeam; etaPlotRange = etaPlotRange_wrtBeam;}

    if(doGen)
    {
        // nTrk is just the number of charged particles
        doNTPC = false;
        doTheta = false;
        doP = false;
        dod0 = false;
        doz0 = false;
        doWW = false;
        nTrkMin = 1;
        doE = false;
        doSTheta = false;
        domissPCut = false;
        doAjCut = false;
        do3jetEvtCut = false;
    }
    return;
}

Float_t Selection::getEtaPlotRange()
{
    if(doThrust) return etaPlotRange_wrtThr;
    if(doWTA) return etaPlotRange_wrtWTA;
    return etaPlotRange_wrtBeam;
}

Float_t Selection::getDifferential()
{
    if(doThrust) return (2*etaPlotRange_wrtThr/(float)dEtaBins)*(2*TMath::Pi()/(float)dPhiBins);
    if(doWTA) return (2*etaPlotRange_wrtWTA/(float)dEtaBins)*(2*TMath::Pi()/(float)dPhiBins);
    return (2*etaPlotRange_wrtBeam/(float)dEtaBins)*(2*TMath::Pi()/(float)dPhiBins);
}
// return 1 if the track passes the selection
// otherwise return 0
int Selection::ridge_trackSelection(Int_t nTPC, Float_t theta, Float_t p, Float_t d0, Float_t z0, Int_t pwflag)
{
    // only charged tracks
    if (pwflag != ALEPH_CHARGED_TRACK) return 0;
    
    // From 1990 "Properties of Hadronic Events in e+e- Annihilation at sqrt(s) = 91 GeV" ALEPH Collaboration paper
    if (doNTPC && (nTPC <= nTPCMin || nTPC >= nTPCMax) ) return 0;
    if (doTheta && (theta <= thetaMin || theta >= thetaMax)) return 0;
    if (doP && ( TMath::Abs(p) < pMin)) return 0;
    if (dod0 && (TMath::Abs(d0) > d0Cut)) return 0;
    if (doz0 && (TMath::Abs(z0) > z0Cut)) return 0;

    return 1;
}

// return -1 if does not pass event selection
// return multiplicity



/////// MUST UPDATE THE EVENT SELECTION BASED ON AXIS/THE 1990 PAPER/////////
int Selection::ridge_eventSelection(bool passesWW, Float_t missP, Int_t nParticle, Int_t nref, Float_t jtpt[], Float_t jteta[], Float_t STheta, Float_t mass[],/* track selection */ Int_t nTPC[], Float_t theta[], Float_t pmag[], Float_t d0[], Float_t z0[], Int_t pwflag[])
{
    
    ///////// QCD Paper Selection /////////
    if (doWW && !passesWW) return -1;

    // From 1990 "Properties of Hadronic Events in e+e- Annihilation at sqrt(s) = 91 GeV" ALEPH Collaboration paper
    Int_t N = 0;
    Float_t E = 0;
    int i = 0;
    for (Int_t j=0;j<nParticle;j++)
    {
        if (ridge_trackSelection(nTPC[j],theta[j],pmag[j], d0[j], z0[j], pwflag[j]))
        {
            N += 1;
            E += sqrt(pow(pmag[j],2) + pow(mass[j],2));
        }
    }
    if (N < nTrkMin) return -1;
    // this is to ensure that we are able to find a mixed event
    if( N >= 1000) return -1;

    if (doE && E < TotalChrgEnergyMin) return -1;
    if(doSTheta && (STheta < SThetaMin && STheta > SThetaMax)) return -1;
    //// End of 1990 cuts ////

    ///////// Missing Momentum /////////
    if (domissPCut && missP > missPCut) return -1;
    
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
    
    return N;
}

// return true if the event passes the criteria for a mixed event
// Currently matching the total multiplicity in the event and leading jet eta
bool Selection::isMixedEvent(Int_t nParticle, Int_t nParticle_mix, Float_t jteta, Float_t jteta_mix, Float_t TTheta, Float_t TTheta_mix, Float_t TPhi, Float_t TPhi_mix)
{
    // Original matched event criteria
    // (fabs(mix.nParticle-data.nParticle)>4&&data.nParticle<1000&&fabs(mix.jteta[0]-data.jteta[0])>0.2)
    
    if (doMixedMultCut && TMath::Abs(nParticle - nParticle_mix) > 4) return false;
    if (doMixedJetCut && TMath::Abs(jteta - jteta_mix) > 0.2) return false;
    if (doMixedThrustAxisCut && ( sqrt( pow(TMath::Abs(TTheta - TTheta_mix),2) +  pow(TMath::Abs(TPhi - TPhi_mix),2) ) >= maxTTheta_TPhi_Cut)) return false;
    
    return true;
}
// return histogram number of the histogram corresponding to the multiplicity range
int Selection::histNtrk(Int_t N)
{
    // inclusive on low end to help statistics in higher multiplicity bins
    for (Int_t i = 0; i < nMultBins; ++i)
    {
        if (N >= multBinsLow[i] && N < multBinsHigh[i]) return i;
    }
    return -1;
}

int Selection::histEnergy(Float_t Energy)
{
    for (Int_t i = 0; i < nEnergyBins; ++i)
    {
        if (Energy >= energyBinsLow[i] && Energy < energyBinsHigh[i]) return i;
    }
    return -1;
}

// return histogram number of the histogram corresponding to the pt range
int Selection::histPt(Float_t pt)
{
    for (Int_t i = 0; i < nptBins; ++i)
    {
        if (doThrust && (pt >= ptBinsLow_wrtThr[i] && pt < ptBinsHigh_wrtThr[i])) return i;
        else if (doWTA && (pt >= ptBinsLow_wrtThr[i] && pt < ptBinsHigh_wrtThr[i])) return i;
        else if (pt >= ptBinsLow_wrtBeam[i] && pt < ptBinsHigh_wrtBeam[i]) return i;
    }
    return -1;
}
// return histogram number of the histogram corresponding to the pt range
int Selection::histEta(Float_t eta)
{
    int numEta = 0;
    // we iterate backwards because we want to fill all of the bins with an eta fill less than or equal to the min 
    for (Int_t i = netaBins-1; i >=0; --i)
    {
        if (doThrust && (eta < etaBinsLow_wrtThr[i])) return i;
        else if (doWTA && (eta < etaBinsLow_wrtThr[i])) return i;
        else if (eta < etaBinsLow_wrtBeam[i]) return i;
    }
    return -1;
}
#endif /* Selection.h */
