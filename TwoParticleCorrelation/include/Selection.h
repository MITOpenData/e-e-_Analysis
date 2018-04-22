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
#include <TH1F.h>
#define PI 3.14159265358979

//DataProcessing header
#include "DataProcessing/include/trackSelection.h"
#include "DataProcessing/include/neutralHadronSelection.h"
#include "DataProcessing/include/particleData.h"
#include "DataProcessing/include/eventData.h"
#include "DataProcessing/include/jetData.h"

using namespace std;

inline float jtp(float pt, float eta){ return pt*TMath::CosH(eta);}

class Selection
{
    public:
    
        // Initial Setup
	TrackSelection trackSelector;
	NeutralHadronSelection neutralHadronSelector;
        bool doParallel = true;
        bool doOneEvent = false;     int numEvents = 50000;
        bool doBelle = false;
        Int_t experiment = 0; // 0 ALEPH , 1 DELPHI, 2 BELLE, 3 CMS pp
        Int_t jttree = 0; // 0 ak4ESchemeJetTree, 1 ak4WTAmodpSchemeJetTree, 2 ak8ESchemeJetTree, 3 ak8WTAmodpSchemeJetTree, 4 ktN2WTAmodpSchemeJetTree
        Int_t nbin = 20;
        Int_t bkgrd_runs = 1;
        enum SIMPLEPID {BELLE_PHOTON, BELLE_ELECTRON, BELLE_PION, BELLE_MUON, BELLE_KAON, BELLE_PROTON};    // BELLE Particle Definition
        enum SIMPLEPWFLAG {ALEPH_CHARGED_TRACK, ALEPH_CHARGED_LEPTONS1, ALEPH_CHARGED_LEPTONS2, ALEPH_V0, ALEPH_PHOTON, ALEPH_NEUTRAL_HADRON};  // ALEPH Particle Flow Classification
        bool doGen = false; // turn of event and track selections
        bool getThetaAngle = false;
        
        /* Detector Specific Cuts */
        
        // From 1990 "Properties of Hadronic Events in e+e- Annihilation at sqrt(s) = 91 GeV" ALEPH Collaboration paper
        // Event Cuts
        bool doE = true; Float_t TotalChrgEnergyMin = 15; // measured in GeV
        Float_t nTrkMin = 5;
        bool doSTheta = false; Float_t SThetaMax = 0.4;  
        // end of 1990 cuts
        bool domissPCut = false;  Float_t missPCut = 20;  // measured in GeV
        
	/* Frame Dependent Cuts */
        bool doWW = false;  // impose WW cut
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
        static const Int_t nptBins_wrtBeam = 3;
        Float_t ptBinsLow_wrtBeam[nptBins_wrtBeam]  = {0.4,0.4,1.0};  // measured in GeV
        Float_t ptBinsHigh_wrtBeam[nptBins_wrtBeam] = {100.0,4.0,3.0};
        static const Int_t netaBins_wrtBeam = 2;
        Float_t etaBinsLow_wrtBeam[netaBins_wrtBeam]  = {1.6, 1.8};
        Float_t missPCut_wrtBeam = 20;
        Float_t etaPlotRange_wrtBeam = 6.0;

        // Thrust Axis
        bool doThrust = true;  bool donTrkThrust = false; // false = use beam axis for nTrk calculation
        static const Int_t nptBins_wrtThr = 3;
        Float_t ptBinsLow_wrtThr[nptBins_wrtThr]  = {0.4,0.4,1.0};  // measured in GeV {1.0,0.4}
        Float_t ptBinsHigh_wrtThr[nptBins_wrtThr] = {100.0,4.0,3.0}; // {3.0,100.0}
        static const Int_t netaBins_wrtThr = 1;
        Float_t etaBinsLow_wrtThr[netaBins_wrtThr]  = {4.5}; //{4.5,5.0}
        Float_t missPCut_wrtThr = 20;
        Float_t etaPlotRange_wrtThr = 6.0;

        // WTA Axis
        bool doWTA = false;
        static const Int_t nptBins_wrtWTA = 3;
        Float_t ptBinsLow_wrtWTA[nptBins_wrtWTA]  = {0.4,0.4,1.0};  // measured in GeV
        Float_t ptBinsHigh_wrtWTA[nptBins_wrtWTA] = {100,4.0,3.0};
        static const Int_t netaBins_wrtWTA = 1;
        Float_t etaBinsLow_wrtWTA[netaBins_wrtWTA]  = {5.0};
        Float_t missPCut_wrtWTA = 20;
        Float_t etaPlotRange_wrtWTA = 6.0;

        // Use the perpendicular direction of WTA or Thrust axis?
	bool doPerp = true;

        /* Independent Cuts */
        static const Int_t nEnergyBins = 1;
        Float_t energyBinsLow[nEnergyBins] = {0}; // {0,100}
        Float_t energyBinsHigh[nEnergyBins] = {100}; // {100,999}

        static const Int_t nMultBins = 4;
        Int_t multBinsLow[nMultBins]  = {0 , 20, 30, 35};
        Int_t multBinsHigh[nMultBins] = {20, 30, 999, 999};

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
        Float_t ptBinsLow_pp[nptBins_pp]  = {0.4, 1.0};  // measured in GeV
        Float_t ptBinsHigh_pp[nptBins_pp] = {100, 3.0};
        static const Int_t netaBins_pp = 2;
        Float_t etaBinsLow_pp[netaBins_pp]  = {4.5,5.0};
        Float_t missPCut_pp = 20;
        Float_t etaPlotRange_pp = 6.0;
        

        /*In this section we define the variables used to perform an analysis inside and outside what we call "barrel" that is a cilindrical region around the thrust axis.
        The variable listanalysisregion is a selector that allows to decide whether to perform an thrust analysis in the entire space, in a barrel region around of thrust axis of eta limits = etabarrelcut*/
        
        enum listanalysisregion{ kAnaTypeAllDetector, kAnaTypeInBarrelThrust, kAnaTypeOutBarrelThrust};
        int anatyperegion=kAnaTypeInBarrelThrust;
	    double etabarrelcut=-1;

        /*In this section we define the variables needed to perform a selection with respect to the ratio of the energy of particles (charged and neutral) inside the barrel divided by the total energy
        The variable is indeed EnergyInBarrel/Total energy. The barrel around the thrust axis is defined by 	double etabarrelcutforEselection. The cut on the relative energy is defined by maxrelenergyinsidebarrel.
        Three options are available:
        1) kNoEBarrelThrustSel = no selection
        2) kEInBarrelThrustSel selecting events with a small fraction of energy in barrel
        3) kEOutBarrelThrustSel selecting events with a large fraction of energy in barrel 
        With typemultiplicity, one can choose to use as multiplicity variable to study the multiplicity dependent of the right the total event multiplicity of charged track or the total multiplicity of charged track in the 
        barrel defined by etabarrelcutforEselection as well.*/
          
        enum listtypeEnergyBarrelThrustSel{ kNoEBarrelThrustSel, kEInBarrelThrustSel, kEOutBarrelThrustSel};        
        int typeEnergyBarrelSel=kNoEBarrelThrustSel;
	    double maxrelenergyinsidebarrel=0.;        
	    double etabarrelcutforEselection=2.0;
        enum listtypemultiplicity{ kMultAllEvent, kMultInBarrelThrust};
        int typemultiplicity=kMultAllEvent;

        Selection();
        Float_t getEtaPlotRange();
        Float_t getDifferential();
        int ridge_eventSelection(eventData *event, jetData *jet, particleData *particle);
        bool isMixedEvent(Int_t nParticle, Int_t nParticle_mix, Float_t jteta, Float_t jteta_mix, Float_t TTheta, Float_t TTheta_mix, Float_t TPhi, Float_t TPhi_mix);
        void histEnergy(std::vector<Int_t> &hists, Float_t Energy);
        void histNtrk(std::vector<Int_t> &hists, Int_t N);
        void histPt(std::vector<Int_t> &hists, Float_t pt);
        void histEta(std::vector<Int_t> &hists, Float_t eta);
	void writeMetaData(TH1F *h);
        
    private:
};

Selection::Selection()
{
    std::cout << "Getting settings.." << std::endl;
    if (getThetaAngle) std::cout << "THIS IS A THETA ANALYSIS!!!"<<std::endl;
    if(doPP)
    {
        multBinsLow[2] = 110;
        multBinsHigh[1] = 110;
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
        std::cout<<"Turning off event and track selection criteria..."<<std::endl;
        // nTrk is just the number of charged particles
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

/*
Yen-Jie: merged with DataProcessing/include/trackSelection.h

int Selection::ridge_trackSelection(Int_t j)
{
    // only charged tracks
    if (pwflag != ALEPH_CHARGED_TRACK) return 0;
    
    // From 1990 "Properties of Hadronic Events in e+e- Annihilation at sqrt(s) = 91 GeV" ALEPH Collaboration paper
    if (doNTPC && (nTPC <= nTPCMin || nTPC >= nTPCMax) ) return 0;
    if (doTheta && (theta <= thetaMin || theta >= thetaMax)) return 0;
    if (doP && ( TMath::Abs(p) < pMin)) return 0;
    if (doPt && ( TMath::Abs(pt) < ptMin)) return 0;
    if (dod0 && (TMath::Abs(d0) > d0Cut)) return 0;
    if (doz0 && (TMath::Abs(z0) > z0Cut)) return 0;

    return 1;
}
*/

// return -1 if does not pass event selection
// return multiplicity



/////// MUST UPDATE THE EVENT SELECTION BASED ON AXIS/THE 1990 PAPER/////////
int Selection::ridge_eventSelection(eventData *event, jetData *jet, particleData *particle)
{
    if (doGen) return particle->nParticle;
    ///////// QCD Paper Selection /////////
    if (particle->nParticle<nTrkMin || particle->nParticle<13) return -1;

    // Yen-Jie: probably need to update and apply neutral particle selection

    if (doWW && !event->passesWW) return -1;
    if (doSTheta && TMath::Abs(cos(event->STheta))>=SThetaMax) return -1;

    // From 1990 "Properties of Hadronic Events in e+e- Annihilation at sqrt(s) = 91 GeV" ALEPH Collaboration paper
    Int_t returnNch=-1;
    Int_t Nch = 0;
    Int_t Neu = 0;
    Float_t E = 0;

    for (Int_t j=0;j<particle->nParticle;j++) {
        if (trackSelector.highPurity(particle,j)){
            Nch++;
            E += sqrt(particle->pmag[j]*particle->pmag[j] + particle->mass[j]*particle->mass[j]);
        }
	if (neutralHadronSelector.highPurity(particle,j)) {
	   Neu++;
	}
    }

    if (Nch < nTrkMin) return -1;
    // this is to ensure that we are able to find a mixed event
    if( Nch >= 1000) return -1;
    if ((Neu+Nch)<13) return -1;

    if (doE && E < TotalChrgEnergyMin) return -1;


    //// End of 1990 cuts ////

    ///////// Missing Momentum /////////
    if (domissPCut && event->missP > missPCut) return -1;
    
    ///////// 3-Jet /////////
    if (doAjCut&& jet->nref>=2) {
       Float_t p1 = jtp(jet->jtpt[0],jet->jteta[0]);
       Float_t p2 = jtp(jet->jtpt[1],jet->jteta[1]);
       Float_t j12 = p1+p2;
       fillAj = TMath::Abs(p1-p2) / j12;
       // require 2 jets to be pretty equally balanced in momentum
       if(doAjCut && fillAj > AjCut) return -1;
       // require 3rd jet has low momentum relative to first two jets (take average momentum of jet 1 and 2)
       if(do3jetEvtCut && jet->nref > 2 && (2*jtp(jet->jtpt[2],jet->jteta[2]) / j12) > thirdJetCut ) return -1;
    }
    
    ///////////// barrel energy cut //////////////
    
    if(typeEnergyBarrelSel==1 || typeEnergyBarrelSel==2){
      double totalenergy=0.;
      double totalenergyinbarrel=0.;
    
      for (Int_t j=0;j<particle->nParticle;j++) {
        if (trackSelector.highPurity(particle,j) || neutralHadronSelector.highPurity(particle,j)){
          totalenergy += sqrt(particle->pmag[j]*particle->pmag[j] + particle->mass[j]*particle->mass[j]);
          //std::cout<<"eta"<<particle->eta_wrtThr[j]<<std::endl;
          if(fabs(particle->eta_wrtThr[j])<etabarrelcutforEselection){
            totalenergyinbarrel += sqrt(particle->pmag[j]*particle->pmag[j] + particle->mass[j]*particle->mass[j]);
          }
        }
      }
      //std::cout<<"ratio"<<totalenergyinbarrel/totalenergy<<std::endl;
      if (typeEnergyBarrelSel==1 && totalenergyinbarrel/totalenergy>maxrelenergyinsidebarrel) return -1;
      if (typeEnergyBarrelSel==2 && totalenergyinbarrel/totalenergy<maxrelenergyinsidebarrel) return -1;
    }
    
    returnNch=Nch;
    
    if (typemultiplicity==kMultInBarrelThrust) { 
    /* here we define a new version of the multiplicity that considers only tracks inside the barrel region defined by etabarrelcutforEselection*/
    Int_t NchBarrel = 0;
      for (Int_t j=0;j<particle->nParticle;j++) {
        if (trackSelector.highPurity(particle,j)){
          if(fabs(particle->eta_wrtThr[j])<etabarrelcutforEselection){
            NchBarrel++;
          }
        }  
     }
     returnNch=NchBarrel;
   }
   return returnNch;
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

// return int of histogram for event energy
void Selection::histEnergy(vector<Int_t> &hists, Float_t Energy)
{
    for (Int_t i = 0; i < nEnergyBins; ++i)
    {
        if (Energy >= energyBinsLow[i] && Energy < energyBinsHigh[i]) hists.push_back(i);
    }
}

// return histogram number of the histogram corresponding to the multiplicity range
void Selection::histNtrk(std::vector<Int_t> &hists, Int_t N)
{
    // inclusive on low end to help statistics in higher multiplicity bins
    for (Int_t i = 0; i < nMultBins; ++i)
    {
        if (N >= multBinsLow[i] && N < multBinsHigh[i]) hists.push_back(i);
    }
}

// return vector containing the histogram numbers to fill for pt
void Selection::histPt(std::vector<Int_t> &hists, Float_t pt)
{
    for (Int_t i = 0; i < nptBins; ++i) {
        if (doThrust && (pt >= ptBinsLow_wrtThr[i] && pt < ptBinsHigh_wrtThr[i])) hists.push_back(i);
        else if (doWTA && (pt >= ptBinsLow_wrtThr[i] && pt < ptBinsHigh_wrtThr[i])) hists.push_back(i);
        else if (pt >= ptBinsLow_wrtBeam[i] && pt < ptBinsHigh_wrtBeam[i]) hists.push_back(i);
    }
}

// return vector containing the histogram numbers to fill for eta
void Selection::histEta(std::vector<Int_t> &hists, Float_t eta)
{
    // we iterate backwards because we want to fill all of the bins with an eta fill less than or equal to the min 
    for (Int_t i = 0; i < netaBins; ++i)
    {
        if (doThrust && (TMath::Abs(eta) < etaBinsLow_wrtThr[i])) hists.push_back(i);
        else if (doWTA && (TMath::Abs(eta) < etaBinsLow_wrtThr[i])) hists.push_back(i);
        else if (TMath::Abs(eta) < etaBinsLow_wrtBeam[i]) hists.push_back(i);
    }
}

void Selection::writeMetaData(TH1F *h)
{
   h->SetBinContent(1,experiment);
   h->SetBinContent(2,doThrust);
   h->SetBinContent(3,doWTA);
   h->SetBinContent(4,doPerp);
   h->SetBinContent(5,doGen);
   h->SetBinContent(6,getThetaAngle);
   h->SetBinContent(7,nEnergyBins);
   h->SetBinContent(8,nMultBins);
   h->SetBinContent(9,nptBins);
   h->SetBinContent(10,netaBins);
      
}

/*
void Selection::readMetaData(TH1F *h)
{
   experiment = h->GetBinContent(1);
   doThrust   = h->SetBinContent(2);
   doWTA      = h->GetBinContent(3);
   doPerp     = h->GetBinContent(4);
   doGen      = h->GetBinContent(5);
   getThetaAngle = h->GetBinContent(6);
   nEnergyBins = h->GetBinContent(7);
   nMultBins   = h->GetBinContent(8);
   nptBins     = h->GetBinContent(9);
   netaBins    = h->GetBinContent(10);
}
*/

#endif /* Selection.h */
