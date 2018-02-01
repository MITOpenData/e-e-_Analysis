#ifndef SETTINGS
#define SETTINGS

#include <string>
#include <iostream>
#include "TH1F.h"

class Settings{
  public:

    //0 = ALEPH (defaults, rest are initialized in constructor)
    //1 - DELPHI
    //2 = Belle
    //3 = CMS

    int experiment = 0;
    bool isMC = false;
    int MCProcess = 5;
    //ALEPH data
    std::string inputFile = "/data/cmcginn/StudyMultSamples/ALEPH/LEP2/20180119/LEP2MC1996YDATAEMINI_recons_aftercut-MERGED.root";
    //std::string inputFile = "/data/cmcginn/GeneratorsHEPMC/PYTHIA8/0p0912/pthat1_Zee/20171208/outFile_MERGED_nEvt2000000_nMinChgPart0_RopeWalk0.root";//regular pythia 8
    //std::string inputFile = "/data/cmcginn/GeneratorsHEPMC/PYTHIA8/0p0912/pthat1_Zee_RopeWalk/20171208/outFile_MERGED_nEvt2000000_nMinChgPart0_RopeWalk1.root";//ropewalk pythia 8
    //std::string inputFile = "/data/cmcginn/GeneratorsHEPMC/PYTHIA8/0p0912/pthat1_Zee_minChgNPart30/20171208/outFile_MERGED__nEvt20000_nMinChgPart30_RopeWalk0.root";//high mult regular pythia
    //std::string inputFile = "/data/cmcginn/GeneratorsHEPMC/PYTHIA8/0p0912/pthat1_Zee_RopeWalk_minChgNPart30/20171208/outFile_MERGED_nEvt20000_nMinChgPart30_RopeWalk1.root";//high mult ropewalk


    //cuts
    bool doUseLeptons = false;

    //kinematics (if trig != assoc cuts, make sure doExcludeNTrigLT2 is set to false)
    float trigPt[2] = {0.4,100};
    float assocPt[2] = {0.4,100};
    float nTrkPt[2] = {0.4,100};
    
    //float etaCut = 1.8;//BEAM AXIS
    float etaCut = 5.0;//THRUST AXIS
    //beam  axis stuff
    //float etaPlotRange = 1.8;//this gets multiplied by 2
    //float dEtaBins = 36;//keep even
    //float dPhiBins = 36;//keep factor of 4
    //float etaPlotRange = 3.2;//BEAM AXIS this gets multiplied by 2 to give full deta range
    float etaPlotRange = 6.0;//THRUST AXIS this gets multiplied by 2
    float dEtaBins = 20;//keep even
    float dPhiBins = 20;//keep factor of 4

    float dEtaRangeToIntegrate[2] = {2.0,3.6};//try to make this correspond with bin edges based on above parameters

    //mixing
    int nMixedEvents = 1;
    int maxSkipSize  = 1;


    //plots
    bool useBeamMult = false;
    static const int nMultBins = 3;
    int multBinsLow[nMultBins]  = {0 , 20, 30};
    int multBinsHigh[nMultBins] = {20, 30, 999};

    bool calcKinematicsWrtThrust = false;

    //other
    bool doThrust = true;
    bool doChargedThrust = false;
    float thrustMatchWindow = 99.0;
    bool doMultMatch = true;
    bool doMissPCut = false;
    float MissPCut = 20;
    bool doExcludeNTrigLT2 = true;
    bool doAjCut = false;
    bool keep3jetEvts = false;
    float AjCut = 0.1;
    float thirdJetCut = 0.03;
    bool doAllData = true;
    int nEvts = 105;

    Settings();
    bool isInMultBin(int n, int bin);
    bool isInSameMultBin(int n1, int n2);

  private:

};

Settings::Settings()
{
  std::cout << "Getting settings.." << std::endl;
  if(experiment == 3){//CMS
    inputFile = "/data/flowex/CMSsample/TPCNtuple_MinBias_TuneCUETP8M1_5p02TeV-pythia8-HINppWinter16DR-NoPU_75X_mcRun2_asymptotic_ppAt5TeV_forest_v2_track.root";
    nMixedEvents = 5;
    trigPt[0] = 1; trigPt[1] = 3;
    assocPt[0] = 1; assocPt[1] = 3;
    nTrkPt[0] = 0.4; nTrkPt[1]=100;
    etaCut = 2.4;
    dEtaBins = 40;//keep even
    dPhiBins = 40;//keep factor of 4
    dEtaRangeToIntegrate[0] = 2; dEtaRangeToIntegrate[1] = 4;
    doThrust = false;
  }

  return;
}

bool Settings::isInMultBin(int n, int bin){
  if(bin >= nMultBins) std::cout << "Error in isInMultBin(): bin out of bounds!" << std::endl;
  if(n >= multBinsLow[bin] && n < multBinsHigh[bin]) return true;
  return false;
}

bool Settings::isInSameMultBin(int n1, int n2){
  bool isIt = false;
  for(int i = 0; i<nMultBins; i++){
    if(isInMultBin(n1,i) && isInMultBin(n2,i)) isIt = true;
  }
  return isIt;
}
#endif
