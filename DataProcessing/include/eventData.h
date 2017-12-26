#ifndef EVENTDATA_H
#define EVENTDATA_H

#include <string>

class eventData{
 public:
  Float_t missP;
  Float_t missPt;
  Float_t missTheta;
  Float_t missPhi;
  Float_t missChargedP;
  Float_t missChargedPt;
  Float_t missChargedTheta;
  Float_t missChargedPhi;
  Int_t nChargedHadrons;
  Int_t nChargedHadrons_GT0p4;
  Int_t nChargedHadrons_GT0p4Thrust;

  //thrust axis variables
  Float_t Thrust;
  Float_t TTheta;
  Float_t TPhi;
  Float_t Thrust_charged;
  Float_t TTheta_charged;
  Float_t TPhi_charged;

  static const int nVar = 17;
  std::string varStr[nVar] = {"missP",
			      "missPt",
			      "missTheta",
			      "missPhi",
			      "missChargedP",
			      "missChargedPt",
			      "missChargedTheta",
			      "missChargedPhi",
			      "nChargedHadrons",
			      "nChargedHadrons_GT0p4",
			      "nChargedHadrons_GT0p4Thrust",
			      "Thrust",
			      "TTheta",
			      "TPhi",
			      "Thrust_charged",
			      "TTheta_charged",
			      "TPhi_charged"};

  bool varIsGood[nVar];

  eventData();
  void SetStatusAndAddressRead(TTree* inTree_p, std::vector<std::string> inList);
};

eventData::eventData()
{
  for(int i = 0; i < nVar; ++i){varIsGood[i] = true;}
  return;
}

void eventData::SetStatusAndAddressRead(TTree* inTree_p, std::vector<std::string> inList)
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

  if(varIsGood[0]) inTree_p->SetBranchAddress("missP", &missP);
  if(varIsGood[1]) inTree_p->SetBranchAddress("missPt", &missPt);
  if(varIsGood[2]) inTree_p->SetBranchAddress("missTheta", &missTheta);
  if(varIsGood[3]) inTree_p->SetBranchAddress("missPhi", &missPhi);
  if(varIsGood[4]) inTree_p->SetBranchAddress("missChargedP", &missChargedP);
  if(varIsGood[5]) inTree_p->SetBranchAddress("missChargedPt", &missChargedPt);
  if(varIsGood[6]) inTree_p->SetBranchAddress("missChargedTheta", &missChargedTheta);
  if(varIsGood[7]) inTree_p->SetBranchAddress("missChargedPhi", &missChargedPhi);
  if(varIsGood[8]) inTree_p->SetBranchAddress("nChargedHadrons", &nChargedHadrons);
  if(varIsGood[9]) inTree_p->SetBranchAddress("nChargedHadrons_GT0p4", &nChargedHadrons_GT0p4);
  if(varIsGood[10]) inTree_p->SetBranchAddress("nChargedHadrons_GT0p4Thrust", &nChargedHadrons_GT0p4Thrust);
  if(varIsGood[11]) inTree_p->SetBranchAddress("Thrust", &Thrust);
  if(varIsGood[12]) inTree_p->SetBranchAddress("TTheta", &TTheta);
  if(varIsGood[13]) inTree_p->SetBranchAddress("TPhi", &TPhi);
  if(varIsGood[14]) inTree_p->SetBranchAddress("Thrust_charged", &Thrust_charged);
  if(varIsGood[15]) inTree_p->SetBranchAddress("TTheta_charged", &TTheta_charged);
  if(varIsGood[16]) inTree_p->SetBranchAddress("TPhi_charged", &TPhi_charged);

  return;
}

#endif
