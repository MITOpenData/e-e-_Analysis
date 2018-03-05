#ifndef EVENTDATA_H
#define EVENTDATA_H

#include <string>

class eventData{
 public:
  Bool_t passesWW;
  Bool_t passesD2;
  Bool_t passesCW;
  Float_t d2;
  Float_t cW;

  Float_t wwInitFourJetSumPx;
  Float_t wwInitFourJetSumPy;
  Float_t wwInitFourJetSumPz;
  Float_t wwInitFourJetSumE;

  Float_t wwFinalFourJetSumPx;
  Float_t wwFinalFourJetSumPy;
  Float_t wwFinalFourJetSumPz;
  Float_t wwFinalFourJetSumE;

  Float_t wwInitFourJetPx[4];
  Float_t wwInitFourJetPy[4];
  Float_t wwInitFourJetPz[4];
  Float_t wwInitFourJetE[4];

  Float_t wwFinalFourJetPx[4];
  Float_t wwFinalFourJetPy[4];
  Float_t wwFinalFourJetPz[4];
  Float_t wwFinalFourJetE[4];

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

  //Sphericity variables
  Float_t Sphericity;
  Float_t STheta;
  Float_t SPhi;
  Float_t Aplanarity;
  Float_t Sphericity_linearized;
  Float_t STheta_linearized;
  Float_t SPhi_linearized;
  Float_t Aplanarity_linearized;
  Float_t C_linearized;
  Float_t D_linearized;

  static const int nVar = 48;
  std::string varStr[nVar] = {"passesWW",
			      "passesD2",
			      "passesCw",
			      "d2",
			      "cW",
			      "wwInitFourJetSumPx",
			      "wwInitFourJetSumPy",
			      "wwInitFourJetSumPz",
			      "wwInitFourJetSumE",
			      "wwFinalFourJetSumPx",
			      "wwFinalFourJetSumPy",
			      "wwFinalFourJetSumPz",
			      "wwFinalFourJetSumE",
			      "wwInitFourJetPx",
			      "wwInitFourJetPy",
			      "wwInitFourJetPz",
			      "wwInitFourJetE",
			      "wwFinalFourJetPx",
			      "wwFinalFourJetPy",
			      "wwFinalFourJetPz",
			      "wwFinalFourJetE",
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
			      "TPhi_charged",
                              "Sphericity",
                              "STheta",
                              "SPhi",
                              "Aplanarity",
                              "Sphericity_linearized",
                              "STheta_linearized",
                              "SPhi_linearized",
                              "Aplanarity_linearized",
                              "C_linearized",
                              "D_linearized"
                             };

  bool varIsGood[nVar];

  eventData();
  void SetStatusAndAddressRead(TTree* inTree_p, std::vector<std::string> inList);
  void SetBranchWrite(TTree* inTree_p);
};

eventData::eventData()
{
  passesWW = false;
  passesD2 = false;
  passesCW = false;
  d2 = -999.;
  cW = -999.;
  wwInitFourJetSumPx = -999.;
  wwInitFourJetSumPy = -999.;
  wwInitFourJetSumPz = -999.;
  wwInitFourJetSumE = -999.;
  wwFinalFourJetSumPx = -999.;
  wwFinalFourJetSumPy = -999.;
  wwFinalFourJetSumPz = -999.;
  wwFinalFourJetSumE = -999.;
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
  Sphericity = -999;
  STheta = -999;
  SPhi = -999;
  Aplanarity = -999;
  Sphericity_linearized = -999;
  STheta_linearized = -999;
  SPhi_linearized = -999;
  Aplanarity_linearized = -999;
  C_linearized = -999;
  D_linearized = -999;

  for(Int_t i = 0; i < 4; ++ i){
    wwInitFourJetPx[i] = -999.;
    wwInitFourJetPy[i] = -999.;
    wwInitFourJetPz[i] = -999.;
    wwInitFourJetE[i] = -999.;

    wwFinalFourJetPx[i] = -999.;
    wwFinalFourJetPy[i] = -999.;
    wwFinalFourJetPz[i] = -999.;
    wwFinalFourJetE[i] = -999.;
  }

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
  if(varIsGood[1]) inTree_p->SetBranchAddress("passesD2", &passesD2);
  if(varIsGood[2]) inTree_p->SetBranchAddress("passesCW", &passesCW);
  if(varIsGood[3]) inTree_p->SetBranchAddress("d2", &d2);
  if(varIsGood[4]) inTree_p->SetBranchAddress("cW", &cW);
  if(varIsGood[5]) inTree_p->SetBranchAddress("wwInitFourJetSumPx", &wwInitFourJetSumPx);
  if(varIsGood[6]) inTree_p->SetBranchAddress("wwInitFourJetSumPy", &wwInitFourJetSumPy);
  if(varIsGood[7]) inTree_p->SetBranchAddress("wwInitFourJetSumPz", &wwInitFourJetSumPz);
  if(varIsGood[8]) inTree_p->SetBranchAddress("wwInitFourJetSumE", &wwInitFourJetSumE);
  if(varIsGood[9]) inTree_p->SetBranchAddress("wwFinalFourJetSumPx", &wwFinalFourJetSumPx);
  if(varIsGood[10]) inTree_p->SetBranchAddress("wwFinalFourJetSumPy", &wwFinalFourJetSumPy);
  if(varIsGood[11]) inTree_p->SetBranchAddress("wwFinalFourJetSumPz", &wwFinalFourJetSumPz);
  if(varIsGood[12]) inTree_p->SetBranchAddress("wwFinalFourJetSumE", &wwFinalFourJetSumE);
  if(varIsGood[13]) inTree_p->SetBranchAddress("wwInitFourJetPx", wwInitFourJetPx);
  if(varIsGood[14]) inTree_p->SetBranchAddress("wwInitFourJetPy", wwInitFourJetPy);
  if(varIsGood[15]) inTree_p->SetBranchAddress("wwInitFourJetPz", wwInitFourJetPz);
  if(varIsGood[16]) inTree_p->SetBranchAddress("wwInitFourJetE", wwInitFourJetE);
  if(varIsGood[17]) inTree_p->SetBranchAddress("wwFinalFourJetPx", wwFinalFourJetPx);
  if(varIsGood[18]) inTree_p->SetBranchAddress("wwFinalFourJetPy", wwFinalFourJetPy);
  if(varIsGood[19]) inTree_p->SetBranchAddress("wwFinalFourJetPz", wwFinalFourJetPz);
  if(varIsGood[20]) inTree_p->SetBranchAddress("wwFinalFourJetE", wwFinalFourJetE);
  if(varIsGood[21]) inTree_p->SetBranchAddress("missP", &missP);
  if(varIsGood[22]) inTree_p->SetBranchAddress("missPt", &missPt);
  if(varIsGood[23]) inTree_p->SetBranchAddress("missTheta", &missTheta);
  if(varIsGood[24]) inTree_p->SetBranchAddress("missPhi", &missPhi);
  if(varIsGood[25]) inTree_p->SetBranchAddress("missChargedP", &missChargedP);
  if(varIsGood[26]) inTree_p->SetBranchAddress("missChargedPt", &missChargedPt);
  if(varIsGood[27]) inTree_p->SetBranchAddress("missChargedTheta", &missChargedTheta);
  if(varIsGood[28]) inTree_p->SetBranchAddress("missChargedPhi", &missChargedPhi);
  if(varIsGood[29]) inTree_p->SetBranchAddress("nChargedHadrons", &nChargedHadrons);
  if(varIsGood[30]) inTree_p->SetBranchAddress("nChargedHadrons_GT0p4", &nChargedHadrons_GT0p4);
  if(varIsGood[31]) inTree_p->SetBranchAddress("nChargedHadrons_GT0p4Thrust", &nChargedHadrons_GT0p4Thrust);
  if(varIsGood[32]) inTree_p->SetBranchAddress("Thrust", &Thrust);
  if(varIsGood[33]) inTree_p->SetBranchAddress("TTheta", &TTheta);
  if(varIsGood[34]) inTree_p->SetBranchAddress("TPhi", &TPhi);
  if(varIsGood[35]) inTree_p->SetBranchAddress("Thrust_charged", &Thrust_charged);
  if(varIsGood[36]) inTree_p->SetBranchAddress("TTheta_charged", &TTheta_charged);
  if(varIsGood[37]) inTree_p->SetBranchAddress("TPhi_charged", &TPhi_charged);
  if(varIsGood[38]) inTree_p->SetBranchAddress("Sphericity", &Sphericity);
  if(varIsGood[39]) inTree_p->SetBranchAddress("STheta", &STheta);
  if(varIsGood[40]) inTree_p->SetBranchAddress("SPhi", &SPhi);
  if(varIsGood[41]) inTree_p->SetBranchAddress("Aplanarity", &Aplanarity);
  if(varIsGood[42]) inTree_p->SetBranchAddress("Sphericity_linearized", &Sphericity_linearized);
  if(varIsGood[43]) inTree_p->SetBranchAddress("STheta_linearized", &STheta_linearized);
  if(varIsGood[44]) inTree_p->SetBranchAddress("SPhi_linearized", &SPhi_linearized);
  if(varIsGood[45]) inTree_p->SetBranchAddress("Aplanarity_linearized", &Aplanarity_linearized);
  if(varIsGood[46]) inTree_p->SetBranchAddress("C_linearized", &C_linearized);
  if(varIsGood[47]) inTree_p->SetBranchAddress("D_linearized", &D_linearized);

  return;
}

void eventData::SetBranchWrite(TTree* inTree_p)
{
  inTree_p->Branch("passesWW", &passesWW, "passesWW/O");
  inTree_p->Branch("passesD2", &passesD2, "passesD2/O");
  inTree_p->Branch("passesCW", &passesCW, "passesCW/O");
  inTree_p->Branch("d2", &d2, "d2/F");
  inTree_p->Branch("cW", &cW, "cW/F");
  inTree_p->Branch("wwInitFourJetSumPx", &wwInitFourJetSumPx, "wwInitFourJetSumPx/F");
  inTree_p->Branch("wwInitFourJetSumPy", &wwInitFourJetSumPy, "wwInitFourJetSumPy/F");
  inTree_p->Branch("wwInitFourJetSumPz", &wwInitFourJetSumPz, "wwInitFourJetSumPz/F");
  inTree_p->Branch("wwInitFourJetSumE", &wwInitFourJetSumE, "wwInitFourJetSumE/F");
  inTree_p->Branch("wwFinalFourJetSumPx", &wwFinalFourJetSumPx, "wwFinalFourJetSumPx/F");
  inTree_p->Branch("wwFinalFourJetSumPy", &wwFinalFourJetSumPy, "wwFinalFourJetSumPy/F");
  inTree_p->Branch("wwFinalFourJetSumPz", &wwFinalFourJetSumPz, "wwFinalFourJetSumPz/F");
  inTree_p->Branch("wwFinalFourJetSumE", &wwFinalFourJetSumE, "wwFinalFourJetSumE/F");
  inTree_p->Branch("wwInitFourJetPx", wwInitFourJetPx, "wwInitFourJetPx[4]/F");
  inTree_p->Branch("wwInitFourJetPy", wwInitFourJetPy, "wwInitFourJetPy[4]/F");
  inTree_p->Branch("wwInitFourJetPz", wwInitFourJetPz, "wwInitFourJetPz[4]/F");
  inTree_p->Branch("wwInitFourJetE", wwInitFourJetE, "wwInitFourJetE[4]/F");
  inTree_p->Branch("wwFinalFourJetPx", wwFinalFourJetPx, "wwFinalFourJetPx[4]/F");
  inTree_p->Branch("wwFinalFourJetPy", wwFinalFourJetPy, "wwFinalFourJetPy[4]/F");
  inTree_p->Branch("wwFinalFourJetPz", wwFinalFourJetPz, "wwFinalFourJetPz[4]/F");
  inTree_p->Branch("wwFinalFourJetE", wwFinalFourJetE, "wwFinalFourJetE[4]/F");
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
  inTree_p->Branch("Sphericity", &Sphericity,"Sphericity/F");
  inTree_p->Branch("STheta", &STheta,"STheta/F");
  inTree_p->Branch("SPhi", &SPhi,"SPhi/F");
  inTree_p->Branch("Aplanarity", &Aplanarity,"Aplanarity/F");
  inTree_p->Branch("Sphericity_linearized", &Sphericity_linearized,"Sphericity_linearized/F");
  inTree_p->Branch("STheta_linearized", &STheta_linearized,"STheta_linearized/F");
  inTree_p->Branch("SPhi_linearized", &SPhi_linearized,"SPhi_linearized/F");
  inTree_p->Branch("Aplanarity_linearized", &Aplanarity_linearized,"Aplanarity_linearized/F");
  inTree_p->Branch("C_linearized", &C_linearized,"C_linearized/F");
  inTree_p->Branch("D_linearized", &D_linearized,"D_linearized/F");
  
  return;
}


#endif
