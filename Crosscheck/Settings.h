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
    //ALEPH data
    std::string inputFile = "/data/abaty/ALEPHTrees/cleaned_ALEPH_Data2-v3_Aug11_2017.root";

    //cuts
    bool doUseLeptons = false;

    //kinematics (if trig != assoc cuts, make sure doExcludeNTrigLT2 is set to false)
    float trigPt[2] = {0,999};
    float assocPt[2] = {0,999};

    float etaCut = 1.8;
    float dEtaBins = 36;//keep even
    float dPhiBins = 36;//keep factor of 4

    float dEtaRangeToIntegrate[2] = {2.0,3.6};//try to make this correspond with bin edges based on above parameters

    //mixing
    int nMixedEvents = 5;
    int maxSkipSize  = 3;


    //plots
    static const int nMultBins = 5;
    int multBinsLow[nMultBins]  = {0,   0 , 15, 25, 35};
    int multBinsHigh[nMultBins] = {999, 15, 25, 35, 999};

    //other
    bool doExcludeNTrigLT2 = true;
    bool doAllData = true;
    int nEvts = 1000000;

    Settings();
    bool isInMultBin(int n, int bin);

  private:

};

Settings::Settings()
{
  std::cout << "Getting settings.." << std::endl;
  if(experiment == 3){//CMS
    inputFile = "/data/flowex/CMSsample/TPCNtuple_MinBias_TuneCUETP8M1_5p02TeV-pythia8-HINppWinter16DR-NoPU_75X_mcRun2_asymptotic_ppAt5TeV_forest_v2_track.root";
    trigPt[0] = 1; trigPt[1] = 3;
    assocPt[0] = 1; assocPt[1] = 3;
    etaCut = 2.4;
    dEtaBins = 40;//keep even
    dPhiBins = 40;//keep factor of 4
    dEtaRangeToIntegrate[0] = 2; dEtaRangeToIntegrate[1] = 4;
  }

  return;
}

bool Settings::isInMultBin(int n, int bin){
  if(bin >= nMultBins) std::cout << "Error in isInMultBin(): bin out of bounds!" << std::endl;
  if(n >= multBinsLow[n] && n < multBinsHigh[bin]) return true;
  return false;
}

#endif
