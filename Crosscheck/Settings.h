#ifndef SETTINGS
#define SETTINGS

#include <string>
#include <iostream>
#include "TH1F.h"

class Settings{
  public:

    std::string inputFile = "/data/abaty/ALEPHTrees/cleaned_ALEPH_Data2-v3_Aug11_2017.root";

    //cuts
    bool doUseLeptons = true;

    float etaCut = 1.7;
    float dEtaBins = 32;//keep even
    float dPhiBins = 32;//keep factor of 4

    int nMixedEvents = 5;
    int maxSkipSize  = 3;

    static const int nMultBins = 5;
    int multBinsLow[nMultBins]  = {0,   0 , 15, 25, 35};
    int multBinsHigh[nMultBins] = {999, 15, 25, 35, 999};

    bool doAllData = true;
    int nEvts = 5000;

    Settings();
    bool isInMultBin(int n, int bin);

  private:

};

Settings::Settings()
{
  std::cout << "Getting settings.." << std::endl;
  return;
}

bool Settings::isInMultBin(int n, int bin){
  if(bin >= nMultBins) std::cout << "Error in isInMultBin(): bin out of bounds!" << std::endl;
  if(n >= multBinsLow[n] && n < multBinsHigh[bin]) return true;
  return false;
}

#endif
