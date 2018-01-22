#ifndef EVENTDATA_H
#define EVENTDATA_H

#include <string>

class eventData{
 public:
  Bool_t passesWW;
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

  static const int nVar = 18;
  std::string varStr[nVar] = {"passesWW",
			      "missP",
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
  void SetBranchWrite(TTree* inTree_p);
};

eventData::eventData()
{
  passesWW = false;
  missP = -999;
  missPt = -999;
  missTheta = -999;
  missPhi = -999;
  missChargedP = -999;
  missChargedPt = -999;
  missChargedTheta = -999;
  missChargedPhi = -999;
  nChargedHadrons = -999;
  nChargedHadrons_GT0p4 = -999;
  nChargedHadrons_GT0p4Thrust = -999;
  Thrust = -999;
  TTheta = -999;
  TPhi = -999;
  Thrust_charged = -999;
  TTheta_charged = -999;
  TPhi_charged = -999;

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

  if(varIsGood[0]) inTree_p->SetBranchAddress("passesWW", &passesWW);
  if(varIsGood[1]) inTree_p->SetBranchAddress("missP", &missP);
  if(varIsGood[2]) inTree_p->SetBranchAddress("missPt", &missPt);
  if(varIsGood[3]) inTree_p->SetBranchAddress("missTheta", &missTheta);
  if(varIsGood[4]) inTree_p->SetBranchAddress("missPhi", &missPhi);
  if(varIsGood[5]) inTree_p->SetBranchAddress("missChargedP", &missChargedP);
  if(varIsGood[6]) inTree_p->SetBranchAddress("missChargedPt", &missChargedPt);
  if(varIsGood[7]) inTree_p->SetBranchAddress("missChargedTheta", &missChargedTheta);
  if(varIsGood[8]) inTree_p->SetBranchAddress("missChargedPhi", &missChargedPhi);
  if(varIsGood[9]) inTree_p->SetBranchAddress("nChargedHadrons", &nChargedHadrons);
  if(varIsGood[10]) inTree_p->SetBranchAddress("nChargedHadrons_GT0p4", &nChargedHadrons_GT0p4);
  if(varIsGood[11]) inTree_p->SetBranchAddress("nChargedHadrons_GT0p4Thrust", &nChargedHadrons_GT0p4Thrust);
  if(varIsGood[12]) inTree_p->SetBranchAddress("Thrust", &Thrust);
  if(varIsGood[13]) inTree_p->SetBranchAddress("TTheta", &TTheta);
  if(varIsGood[14]) inTree_p->SetBranchAddress("TPhi", &TPhi);
  if(varIsGood[15]) inTree_p->SetBranchAddress("Thrust_charged", &Thrust_charged);
  if(varIsGood[16]) inTree_p->SetBranchAddress("TTheta_charged", &TTheta_charged);
  if(varIsGood[17]) inTree_p->SetBranchAddress("TPhi_charged", &TPhi_charged);

  return;
}

void eventData::SetBranchWrite(TTree* inTree_p)
{
  inTree_p->Branch("passesWW", &passesWW, "passesWW/F");
  inTree_p->Branch("missP", &missP, "missP/F");
  inTree_p->Branch("missPt", &missPt, "missPt/F");
  inTree_p->Branch("missTheta", &missTheta, "missTheta/F");
  inTree_p->Branch("missPhi", &missPhi, "missPhi/F");
  inTree_p->Branch("missChargedP", &missChargedP, "missChargedP/F");
  inTree_p->Branch("missChargedPt", &missChargedPt, "missChargedPt/F");
  inTree_p->Branch("missChargedTheta", &missChargedTheta, "missChargedTheta/F");
  inTree_p->Branch("missChargedPhi", &missChargedPhi, "missChargedPhi/F");
  inTree_p->Branch("nChargedHadrons", &nChargedHadrons, "nChargedHadrons/I");
  inTree_p->Branch("nChargedHadrons_GT0p4", &nChargedHadrons_GT0p4, "nChargedHadrons_GT0p4/I");
  inTree_p->Branch("nChargedHadrons_GT0p4Thrust", &nChargedHadrons_GT0p4Thrust, "nChargedHadrons_GT0p4Thrust/I");
  inTree_p->Branch("Thrust", &Thrust, "Thrust/F");
  inTree_p->Branch("TTheta", &TTheta, "TTheta/F");
  inTree_p->Branch("TPhi", &TPhi, "TPhi/F");
  inTree_p->Branch("Thrust_charged", &Thrust_charged, "Thrust_charged/F");
  inTree_p->Branch("TTheta_charged", &TTheta_charged, "TTheta_charged/F");
  inTree_p->Branch("TPhi_charged", &TPhi_charged, "TPhi_charged/F");
  
  return;
}


#endif
