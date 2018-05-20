#ifndef GENERALJETVAR_H
#define GENERALJETVAR_H

#include "TLorentzVector.h"
#include "DataProcessing/include/dataTools.h"

//extension of jetData found in DataProcessing directory - takes basic variables + adds a serious of standard inclusive jet measurements

class generalJetVar{
 public:
  static const Int_t nMaxJet = 500;
  static const int pwFlag = 6;

  int nref;
  Float_t jtpt[nMaxJet];
  Float_t jteta[nMaxJet];
  Float_t jtphi[nMaxJet];
  Float_t jtm[nMaxJet];

  static const int nVar = 5;
  std::string varStr[nVar] = {"nref",
			      "jtpt",
			      "jteta",
			      "jtphi",
			      "jtm"};


  TLorentzVector fourJet[nMaxJet];
  bool varIsGood[nVar];

  generalJetVar();
  void SetStatusAndAddressRead(TTree* inTree_p, std::vector<std::string> inList);
  void SetBranchWrite(TTree* inTree_p);
  void preFillClean();
};

generalJetVar::generalJetVar()
{
  nref = -999;
  for(Int_t i = 0; i < nMaxJet; ++i){
    jtpt[i] = -999;
    jtphi[i] = -999;
    jtm[i] = -999;
    jteta[i] = -999;
  }

  for(int i = 0; i < nVar; ++i){varIsGood[i] = true;}
  return;
}

void generalJetVar::SetStatusAndAddressRead(TTree* inTree_p, std::vector<std::string> inList)
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
  if(varIsGood[4]) inTree_p->SetBranchAddress("jtm", jtm);

  return;
}

void generalJetVar::SetBranchWrite(TTree* inTree_p)
{
  inTree_p->Branch("nref", &nref, "nref/I");
  inTree_p->Branch("jtpt", jtpt, "jtpt[nref]/F");
  inTree_p->Branch("jteta", jteta, "jteta[nref]/F");
  inTree_p->Branch("jtphi", jtphi, "jtphi[nref]/F");
  inTree_p->Branch("jtm", jtm, "jtm[nref]/F");

  return;
}

void generalJetVar::preFillClean()
{
  for(Int_t i = 0; i < nref; ++i){
    jtpt[i] = reducedPrecision(jtpt[i]);
    jtphi[i] = reducedPrecision(jtphi[i]);
    jtm[i] = reducedPrecision(jtm[i]);
    jteta[i] = reducedPrecision(jteta[i]);
  }

  return;
} 

#endif
