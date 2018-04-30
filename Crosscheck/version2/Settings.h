#ifndef SETTINGS
#define SETTINGS

#include <string>
#include <iostream>
#include "TH1F.h"

class Settings{
  public:

    bool doParallel = false;

    //ALEPH data
    std::string inputFile = "/data/cmcginn/StudyMultSamples/ALEPH/LEP1/20180426/LEP1Data1992_recons_aftercut-MERGED.root";
    std::string inputFileMix = "/data/cmcginn/StudyMultSamples/ALEPH/LEP1/20180426/LEP1Data1992_recons_aftercut-MERGED_Mix.root";

    //kinematics (if trig != assoc cuts, make sure doExcludeNTrigLT2 is set to false)
    float trigPt[2] = {0.0,100};
    float assocPt[2] = {0.0,100};
    float nTrkPt[2] = {0.0,100};
    
    //beam  axis stuff
    float etaPlotRange = 3.2;//THRUST AXIS this gets multiplied by 2
    float dEtaBins = 20;//keep even
    float dPhiBins = 20;//keep factor of 4

    float dEtaRangeToIntegrate[2] = {2.0,3.6};//try to make this correspond with bin edges based on above parameters


    //plots
    static const int nMultBins = 3;
    int multBinsLow[nMultBins]  = {0 , 20, 30};
    int multBinsHigh[nMultBins] = {20, 30, 999};


    //other
    bool doThrust = false;
    bool doAllData = true;
    int nEvts = 50000;

    Settings();
    bool isInMultBin(int n, int bin);
    bool isInSameMultBin(int n1, int n2);

  private:

};

Settings::Settings()
{
  std::cout << "Getting settings.." << std::endl;
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
