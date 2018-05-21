#ifndef GLOBALJETVAR_H
#define GLOBALJETVAR_H

#include "TLorentzVector.h"
#include "DataProcessing/include/dataTools.h"

//extension of jetData found in DataProcessing directory - takes basic variables + adds a serious of standard inclusive jet measurements

class globalJetVar{
 public:
  float thrustMag;
  float thrustPx;
  float thrustPy;
  float thrustPz;
  float thrustMajorMag;
  float thrustMajorPx;
  float thrustMajorPy;
  float thrustMajorPz;
  float thrustMinorMag;
  float thrustMinorPx;
  float thrustMinorPy;
  float thrustMinorPz;

  static const int nVar = 12;
  std::string varStr[nVar] = {"thrustMag",
			      "thrustPx",
			      "thrustPy",
			      "thrustPz",
			      "thrustMajorMag",
			      "thrustMajorPx",
			      "thrustMajorPy",
			      "thrustMajorPz",
			      "thrustMinorMag",
			      "thrustMinorPx",
			      "thrustMinorPy",
			      "thrustMinorPz"};

  bool varIsGood[nVar];

  globalJetVar();
  void SetStatusAndAddressRead(TTree* inTree_p, std::vector<std::string> inList);
  void SetBranchWrite(TTree* inTree_p);
  void preFillClean();
};

globalJetVar::globalJetVar()
{
  thrustMag = -999;
  thrustPx = -999;
  thrustPy = -999;
  thrustPz = -999;
  thrustMajorMag = -999;
  thrustMajorPx = -999;
  thrustMajorPy = -999;
  thrustMajorPz = -999;
  thrustMinorMag = -999;
  thrustMinorPx = -999;
  thrustMinorPy = -999;
  thrustMinorPz = -999;

  for(int i = 0; i < nVar; ++i){varIsGood[i] = true;}
  return;
}

void globalJetVar::SetStatusAndAddressRead(TTree* inTree_p, std::vector<std::string> inList)
{
  if(inList.size() != 0){
    for(int i = 0; i < nVar; ++i){varIsGood[i] = false;}

    for(unsigned int i = 0; i < inList.size(); ++i){

      for(Int_t j = 0; j < nVar; ++j){
        if(inList.at(i).size() == varStr[j].size() && inList.at(i).find(varStr[j]) != std::string::npos){
          varIsGood[j] = true;
          break;
        }
      }
    }
  }

  for(Int_t i = 0; i < nVar; ++i){
    if(varIsGood[i]) inTree_p->SetBranchStatus(varStr[i].c_str(), 1);
  }

  if(varIsGood[0]) inTree_p->SetBranchAddress("thrustMag", &thrustMag);
  if(varIsGood[1]) inTree_p->SetBranchAddress("thrustPx", &thrustPx);
  if(varIsGood[2]) inTree_p->SetBranchAddress("thrustPy", &thrustPy);
  if(varIsGood[3]) inTree_p->SetBranchAddress("thrustPz", &thrustPz);
  if(varIsGood[4]) inTree_p->SetBranchAddress("thrustMajorMag", &thrustMajorMag);
  if(varIsGood[5]) inTree_p->SetBranchAddress("thrustMajorPx", &thrustMajorPx);
  if(varIsGood[6]) inTree_p->SetBranchAddress("thrustMajorPy", &thrustMajorPy);
  if(varIsGood[7]) inTree_p->SetBranchAddress("thrustMajorPz", &thrustMajorPz);
  if(varIsGood[8]) inTree_p->SetBranchAddress("thrustMinorMag", &thrustMinorMag);
  if(varIsGood[9]) inTree_p->SetBranchAddress("thrustMinorPx", &thrustMinorPx);
  if(varIsGood[10]) inTree_p->SetBranchAddress("thrustMinorPy", &thrustMinorPy);
  if(varIsGood[11]) inTree_p->SetBranchAddress("thrustMinorPz", &thrustMinorPz);

  return;
}

void globalJetVar::SetBranchWrite(TTree* inTree_p)
{
  inTree_p->Branch("thrustMag", &thrustMag, "thrustMag/F");
  inTree_p->Branch("thrustPx", &thrustPx, "thrustPx/F");
  inTree_p->Branch("thrustPy", &thrustPy, "thrustPy/F");
  inTree_p->Branch("thrustPz", &thrustPz, "thrustPz/F");
  inTree_p->Branch("thrustMajorMag", &thrustMajorMag, "thrustMajorMag/F");
  inTree_p->Branch("thrustMajorPx", &thrustMajorPx, "thrustMajorPx/F");
  inTree_p->Branch("thrustMajorPy", &thrustMajorPy, "thrustMajorPy/F");
  inTree_p->Branch("thrustMajorPz", &thrustMajorPz, "thrustMajorPz/F");
  inTree_p->Branch("thrustMinorMag", &thrustMinorMag, "thrustMinorMag/F");
  inTree_p->Branch("thrustMinorPx", &thrustMinorPx, "thrustMinorPx/F");
  inTree_p->Branch("thrustMinorPy", &thrustMinorPy, "thrustMinorPy/F");
  inTree_p->Branch("thrustMinorPz", &thrustMinorPz, "thrustMinorPz/F");

  return;
}

void globalJetVar::preFillClean()
{
  return;
} 

#endif
