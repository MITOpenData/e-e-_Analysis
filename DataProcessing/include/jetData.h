#ifndef JETDATA_H
#define JETDATA_H

#include "TLorentzVector.h"

class jetData{
 public:
  static const Int_t nMaxJet = 500;
  static const int pwFlag = 6;

  int nref;
  Float_t jtpt[nMaxJet];
  Float_t jteta[nMaxJet];
  Float_t jtphi[nMaxJet];
  Int_t jtN[nMaxJet];
  Int_t jtNPW[nMaxJet][pwFlag];
  Float_t jtptFracPW[nMaxJet][pwFlag];

  static const int nVar = 7;
  std::string varStr[nVar] = {"nref",
			      "jtpt",
			      "jteta",
			      "jtphi",
			      "jtN",
			      "jtNPW",
			      "jtptFracPW"};


  TLorentzVector fourJet[nMaxJet];
  bool varIsGood[nVar];

  jetData();
  void SetStatusAndAddressRead(TTree* inTree_p, std::vector<std::string> inList);
  void SetBranchWrite(TTree* inTree_p);
};

jetData::jetData()
{
  nref = -999;
  for(Int_t i = 0; i < nMaxJet; ++i){
    jtpt[i] = -999;
    jtphi[i] = -999;
    jteta[i] = -999;
    jtN[i] = -999;

    for(Int_t j = 0; j < pwFlag; ++j){
      jtNPW[i][j] = -999;
      jtptFracPW[i][j] = -999;
    }
  }

  for(int i = 0; i < nVar; ++i){varIsGood[i] = true;}
  return;
}

void jetData::SetStatusAndAddressRead(TTree* inTree_p, std::vector<std::string> inList)
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

  if(varIsGood[0]) inTree_p->SetBranchAddress("nref", &nref);
  if(varIsGood[1]) inTree_p->SetBranchAddress("jtpt", jtpt);
  if(varIsGood[2]) inTree_p->SetBranchAddress("jteta", jteta);
  if(varIsGood[3]) inTree_p->SetBranchAddress("jtphi", jtphi);
  if(varIsGood[4]) inTree_p->SetBranchAddress("jtN", jtN);
  if(varIsGood[5]) inTree_p->SetBranchAddress("jtNPW", jtNPW);
  if(varIsGood[6]) inTree_p->SetBranchAddress("jtptFracPW", jtptFracPW);

  return;
}

void jetData::SetBranchWrite(TTree* inTree_p)
{
  inTree_p->Branch("nref", &nref, "nref/I");
  inTree_p->Branch("jtpt", jtpt, "jtpt[nref]/F");
  inTree_p->Branch("jteta", jteta, "jteta[nref]/F");
  inTree_p->Branch("jtphi", jtphi, "jtphi[nref]/F");
  inTree_p->Branch("jtN", jtN, "jtN[nref]/F");
  inTree_p->Branch("jtNPW", jtNPW, ("jtNPW[nref][" + std::to_string(pwFlag) + "]/F").c_str());
  inTree_p->Branch("jtptFracPW", jtptFracPW, ("jtptFracPW[nref][" + std::to_string(pwFlag) + "]/F").c_str());

  return;
}

#endif
