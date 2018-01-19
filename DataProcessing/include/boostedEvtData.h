#ifndef BOOSTEDEVTDATA_H
#define BOOSTEDEVTDATA_H

#include "TTree.h"

class boostedEvtData{
 public:
  static const int nMaxBoostedPart = 10000;

  int nParticle;
  float WTAAxis_Theta;
  float WTAAxis_Phi;
  float pt[nMaxBoostedPart];
  float pmag[nMaxBoostedPart];
  float eta[nMaxBoostedPart];
  float theta[nMaxBoostedPart];
  float phi[nMaxBoostedPart];
  float boostx;
  float boosty;
  float boostz;
  float boost;

  static const int nVar = 12;
  std::string varStr[nVar] = {"nParticle",
                              "WTAAxis_Theta",
                              "WTAAxis_Phi",
			      "pt",
			      "pmag",
			      "eta",
			      "theta",
			      "phi",
                              "boostx",
                              "boosty",
                              "boostz",
                              "boost"
			     };

  bool varIsGood[nVar];
  
  boostedEvtData();
  void SetStatusAndAddressRead(TTree* inTree_p, std::vector<std::string> inList);
};

boostedEvtData::boostedEvtData()
{
  for(int i = 0; i < nVar; ++i){varIsGood[i] = true;}
  return;
}

void boostedEvtData::SetStatusAndAddressRead(TTree* inTree_p, std::vector<std::string> inList)
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

  if(varIsGood[0]) inTree_p->SetBranchAddress("nParticle", &nParticle);
  if(varIsGood[1]) inTree_p->SetBranchAddress("WTAAxis_Theta", &WTAAxis_Theta);
  if(varIsGood[2]) inTree_p->SetBranchAddress("WTAAxis_Phi", &WTAAxis_Phi);
  if(varIsGood[3]) inTree_p->SetBranchAddress("pt", pt);
  if(varIsGood[4]) inTree_p->SetBranchAddress("pmag", pmag);
  if(varIsGood[5]) inTree_p->SetBranchAddress("eta", eta);
  if(varIsGood[6]) inTree_p->SetBranchAddress("theta", theta);
  if(varIsGood[7]) inTree_p->SetBranchAddress("phi", phi);
  if(varIsGood[8]) inTree_p->SetBranchAddress("boostx", &boostx);
  if(varIsGood[9]) inTree_p->SetBranchAddress("boosty", &boosty);
  if(varIsGood[10]) inTree_p->SetBranchAddress("boostz", &boostz);
  if(varIsGood[11]) inTree_p->SetBranchAddress("boost", &boost);
  
  return;
}

#endif
