#ifndef GLOBALJETVAR_H
#define GLOBALJETVAR_H

#include "TLorentzVector.h"
#include "DataProcessing/include/dataTools.h"

//extension of jetData found in DataProcessing directory - takes basic variables + adds a serious of standard inclusive jet measurements

class globalJetVar{
 public:
  int EventNo;
  int RunNo;
  int year;
  int subDir;
  int process;
  bool isMC;
  unsigned long long uniqueID;
  float energy;
  int nParticle;
  float thrustMag; // see definition from http://cds.cern.ch/record/690637/files/ep-2003-084.pdf, p.14
  float thrustPx;
  float thrustPy;
  float thrustPz;
  float thrustMajorMag; // see definition from http://cds.cern.ch/record/690637/files/ep-2003-084.pdf, p.14
  float thrustMajorPx;
  float thrustMajorPy;
  float thrustMajorPz;
  float thrustMinorMag; // see definition from http://cds.cern.ch/record/690637/files/ep-2003-084.pdf, p.15
  float thrustMinorPx;
  float thrustMinorPy;
  float thrustMinorPz;
  float oblateness; // see definition from http://cds.cern.ch/record/690637/files/ep-2003-084.pdf, p.15
  float sphericity; // see definition from http://cds.cern.ch/record/690637/files/ep-2003-084.pdf, p.15
  float aplanarity; // see definition from http://cds.cern.ch/record/690637/files/ep-2003-084.pdf, p.15
  float planarity; // see definition from http://cds.cern.ch/record/690637/files/ep-2003-084.pdf, p.15
  float eVis;
  float heavyJetMass; // see definition from http://cds.cern.ch/record/690637/files/ep-2003-084.pdf, p.15
  float lightJetMass; // see definition from http://cds.cern.ch/record/690637/files/ep-2003-084.pdf, p.15
  float jetMassDifference; // see definition from http://cds.cern.ch/record/690637/files/ep-2003-084.pdf, p.16
  float wideJetBroadening; // see definition from http://cds.cern.ch/record/690637/files/ep-2003-084.pdf, p.16
  float narrowJetBroadening; // see definition from http://cds.cern.ch/record/690637/files/ep-2003-084.pdf, p.16
  float totalJetBroadening; // see definition from http://cds.cern.ch/record/690637/files/ep-2003-084.pdf, p.16
  float cParam; // see definition from http://cds.cern.ch/record/690637/files/ep-2003-084.pdf, p.16
  float jetResParam4; // see definition from http://cds.cern.ch/record/690637/files/ep-2003-084.pdf, p.16
  float jetResParam4NegLog; // see definition from http://cds.cern.ch/record/690637/files/ep-2003-084.pdf, p.16

  static const int nVar = 35;
  std::string varStr[nVar] = {"EventNo",
			      "RunNo",
			      "year",
			      "subDir",
			      "process",
			      "isMC",
			      "uniqueID",
			      "energy",
			      "nParticle",
			      "thrustMag",
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
			      "thrustMinorPz",
			      "oblateness",
			      "sphericity",
			      "aplanarity",
			      "planarity",
			      "eVis",
			      "heavyJetMass",
			      "lightJetMass",
			      "jetMassDifference",
			      "wideJetBroadening",
			      "narrowJetBroadening",
			      "totalJetBroadening",
			      "cParam",
			      "jetResParam4",
			      "jetResParam4NegLog"};

  bool varIsGood[nVar];

  globalJetVar();
  void SetStatusAndAddressRead(TTree* inTree_p, std::vector<std::string> inList);
  void SetBranchWrite(TTree* inTree_p);
  void preFillClean();
};

globalJetVar::globalJetVar()
{
  EventNo = -999;
  RunNo = -999;
  year = -999;
  subDir = -999;
  process = -999;
  isMC = false;
  uniqueID = 0;
  energy = -999;
  nParticle = -999;
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
  oblateness = -999;
  sphericity = -999;
  aplanarity = -999;
  planarity = -999;
  eVis = -999;
  heavyJetMass = -999;
  lightJetMass = -999;
  jetMassDifference = -999;
  wideJetBroadening = -999;
  narrowJetBroadening = -999;
  totalJetBroadening = -999;
  cParam = -999;
  jetResParam4 = -999;
  jetResParam4NegLog = -999;

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


  if(varIsGood[0]) inTree_p->SetBranchAddress("EventNo", &EventNo);
  if(varIsGood[1]) inTree_p->SetBranchAddress("RunNo", &RunNo);
  if(varIsGood[2]) inTree_p->SetBranchAddress("year", &year);
  if(varIsGood[3]) inTree_p->SetBranchAddress("subDir", &subDir);
  if(varIsGood[4]) inTree_p->SetBranchAddress("process", &process);
  if(varIsGood[5]) inTree_p->SetBranchAddress("isMC", &isMC);
  if(varIsGood[6]) inTree_p->SetBranchAddress("uniqueID", &uniqueID);
  if(varIsGood[7]) inTree_p->SetBranchAddress("energy", &energy);
  if(varIsGood[8]) inTree_p->SetBranchAddress("nParticle", &nParticle);
  if(varIsGood[9]) inTree_p->SetBranchAddress("thrustMag", &thrustMag);
  if(varIsGood[10]) inTree_p->SetBranchAddress("thrustPx", &thrustPx);
  if(varIsGood[11]) inTree_p->SetBranchAddress("thrustPy", &thrustPy);
  if(varIsGood[12]) inTree_p->SetBranchAddress("thrustPz", &thrustPz);
  if(varIsGood[13]) inTree_p->SetBranchAddress("thrustMajorMag", &thrustMajorMag);
  if(varIsGood[14]) inTree_p->SetBranchAddress("thrustMajorPx", &thrustMajorPx);
  if(varIsGood[15]) inTree_p->SetBranchAddress("thrustMajorPy", &thrustMajorPy);
  if(varIsGood[16]) inTree_p->SetBranchAddress("thrustMajorPz", &thrustMajorPz);
  if(varIsGood[17]) inTree_p->SetBranchAddress("thrustMinorMag", &thrustMinorMag);
  if(varIsGood[18]) inTree_p->SetBranchAddress("thrustMinorPx", &thrustMinorPx);
  if(varIsGood[19]) inTree_p->SetBranchAddress("thrustMinorPy", &thrustMinorPy);
  if(varIsGood[20]) inTree_p->SetBranchAddress("thrustMinorPz", &thrustMinorPz);
  if(varIsGood[21]) inTree_p->SetBranchAddress("oblateness", &oblateness);
  if(varIsGood[22]) inTree_p->SetBranchAddress("sphericity", &sphericity);
  if(varIsGood[23]) inTree_p->SetBranchAddress("aplanarity", &aplanarity);
  if(varIsGood[24]) inTree_p->SetBranchAddress("planarity", &planarity);
  if(varIsGood[25]) inTree_p->SetBranchAddress("eVis", &eVis);  
  if(varIsGood[26]) inTree_p->SetBranchAddress("heavyJetMass", &heavyJetMass);
  if(varIsGood[27]) inTree_p->SetBranchAddress("lightJetMass", &lightJetMass);
  if(varIsGood[28]) inTree_p->SetBranchAddress("jetMassDifference", &jetMassDifference);
  if(varIsGood[29]) inTree_p->SetBranchAddress("wideJetBroadening", &wideJetBroadening);
  if(varIsGood[30]) inTree_p->SetBranchAddress("narrowJetBroadening", &narrowJetBroadening);
  if(varIsGood[31]) inTree_p->SetBranchAddress("totalJetBroadening", &totalJetBroadening);
  if(varIsGood[32]) inTree_p->SetBranchAddress("cParam", &cParam);
  if(varIsGood[33]) inTree_p->SetBranchAddress("jetResParam4", &jetResParam4);
  if(varIsGood[34]) inTree_p->SetBranchAddress("jetResParam4NegLog", &jetResParam4NegLog);

  return;
}

void globalJetVar::SetBranchWrite(TTree* inTree_p)
{
  inTree_p->Branch("EventNo", &EventNo, "EventNo/I");
  inTree_p->Branch("RunNo", &RunNo, "RunNo/I");
  inTree_p->Branch("year", &year, "year/I");
  inTree_p->Branch("subDir", &subDir, "subDir/I");
  inTree_p->Branch("process", &process, "process/I");
  inTree_p->Branch("isMC", &isMC, "isMC/O");
  inTree_p->Branch("uniqueID", &uniqueID, "uniqueID/l");
  inTree_p->Branch("energy", &energy, "energy/F");
  inTree_p->Branch("nParticle", &nParticle, "nParticle/F");
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
  inTree_p->Branch("oblateness", &oblateness, "oblateness/F");
  inTree_p->Branch("sphericity", &sphericity, "sphericity/F");
  inTree_p->Branch("aplanarity", &aplanarity, "aplanarity/F");
  inTree_p->Branch("planarity", &planarity, "planarity/F");
  inTree_p->Branch("eVis", &eVis, "eVis/F");
  inTree_p->Branch("heavyJetMass", &heavyJetMass, "heavyJetMass/F");
  inTree_p->Branch("lightJetMass", &lightJetMass, "lightJetMass/F");
  inTree_p->Branch("jetMassDifference", &jetMassDifference, "jetMassDifference/F");
  inTree_p->Branch("wideJetBroadening", &wideJetBroadening, "wideJetBroadening/F");
  inTree_p->Branch("narrowJetBroadening", &narrowJetBroadening, "narrowJetBroadening/F");
  inTree_p->Branch("totalJetBroadening", &totalJetBroadening, "totalJetBroadening/F");
  inTree_p->Branch("cParam", &cParam, "cParam/F");
  inTree_p->Branch("jetResParam4", &jetResParam4, "jetResParam4/F");
  inTree_p->Branch("jetResParam4NegLog", &jetResParam4NegLog, "jetResParam4NegLog/F");

  return;
}

void globalJetVar::preFillClean()
{
  return;
} 

#endif
