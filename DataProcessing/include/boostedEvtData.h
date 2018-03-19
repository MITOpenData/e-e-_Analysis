#ifndef BOOSTEDEVTDATA_H
#define BOOSTEDEVTDATA_H

#include "TTree.h"

#include "dataTools.h"

class boostedEvtData{
 public:
  static const int nMaxBoostedPart = 1000;

  int nParticle;
  float WTAAxis_Theta;
  float WTAAxis_ThetaPerp;
  float WTAAxis_Phi;
  float pt[nMaxBoostedPart];
  float pmag[nMaxBoostedPart];
  float eta[nMaxBoostedPart];
  float theta[nMaxBoostedPart];
  float phi[nMaxBoostedPart];
  float pt_Perp[nMaxBoostedPart];
  float pmag_Perp[nMaxBoostedPart];
  float eta_Perp[nMaxBoostedPart];
  float theta_Perp[nMaxBoostedPart];
  float phi_Perp[nMaxBoostedPart];
  float boostx;
  float boosty;
  float boostz;
  float boost;
  float particleWeight;

  static const int nVar = 19;
  std::string varStr[nVar] = {"nParticle",
                              "WTAAxis_Theta",
                              "WTAAxis_ThetaPerp",
                              "WTAAxis_Phi",
			      "pt",
			      "pmag",
			      "eta",
			      "theta",
			      "phi",
			      "pt_Perp",
			      "pmag_Perp",
			      "eta_Perp",
			      "theta_Perp",
			      "phi_Perp",
                              "boostx",
                              "boosty",
                              "boostz",
                              "boost",
                              "particleWeight",
};

  bool varIsGood[nVar];
  
  boostedEvtData();
  void SetStatusAndAddressRead(TTree* inTree_p, std::vector<std::string> inList);
  void SetBranchWrite(TTree* inTree_p);
  void preFillClean();
};

boostedEvtData::boostedEvtData()
{
  nParticle = -999;
  WTAAxis_Theta = -999;
  WTAAxis_ThetaPerp = -999;
  WTAAxis_Phi = -999;
  boostx = -999;
  boosty = -999;
  boostz = -999;
  boost = -999;
  particleWeight = 1; //should be initialized to 1 if not used in mixing

  for(Int_t pI = 0; pI < nMaxBoostedPart; ++pI){
    pt[pI] = -999;
    pmag[pI] = -999;
    eta[pI] = -999;
    theta[pI] = -999;
    phi[pI] = -999;

    pt_Perp[pI] = -999;
    pmag_Perp[pI] = -999;
    eta_Perp[pI] = -999;
    theta_Perp[pI] = -999;
    phi_Perp[pI] = -999;
  }

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
  if(varIsGood[2]) inTree_p->SetBranchAddress("WTAAxis_ThetaPerp", &WTAAxis_ThetaPerp);
  if(varIsGood[3]) inTree_p->SetBranchAddress("WTAAxis_Phi", &WTAAxis_Phi);
  if(varIsGood[4]) inTree_p->SetBranchAddress("pt", pt);
  if(varIsGood[5]) inTree_p->SetBranchAddress("pmag", pmag);
  if(varIsGood[6]) inTree_p->SetBranchAddress("eta", eta);
  if(varIsGood[7]) inTree_p->SetBranchAddress("theta", theta);
  if(varIsGood[8]) inTree_p->SetBranchAddress("phi", phi);
  if(varIsGood[9]) inTree_p->SetBranchAddress("pt_Perp", pt_Perp);
  if(varIsGood[10]) inTree_p->SetBranchAddress("pmag_Perp", pmag_Perp);
  if(varIsGood[11]) inTree_p->SetBranchAddress("eta_Perp", eta_Perp);
  if(varIsGood[12]) inTree_p->SetBranchAddress("theta_Perp", theta_Perp);
  if(varIsGood[13]) inTree_p->SetBranchAddress("phi_Perp", phi_Perp);
  if(varIsGood[14]) inTree_p->SetBranchAddress("boostx", &boostx);
  if(varIsGood[15]) inTree_p->SetBranchAddress("boosty", &boosty);
  if(varIsGood[16]) inTree_p->SetBranchAddress("boostz", &boostz);
  if(varIsGood[17]) inTree_p->SetBranchAddress("boost", &boost);
  if(varIsGood[18]) inTree_p->SetBranchAddress("particleWeight", &particleWeight);
  
  return;
}

void boostedEvtData::SetBranchWrite(TTree* inTree_p)
{
  inTree_p->Branch("nParticle", &nParticle, "nParticle/I");
  inTree_p->Branch("WTAAxis_Theta", &WTAAxis_Theta, "WTAAxis_Theta/F");
  inTree_p->Branch("WTAAxis_ThetaPerp", &WTAAxis_ThetaPerp, "WTAAxis_ThetaPerp/F");
  inTree_p->Branch("WTAAxis_Phi", &WTAAxis_Phi, "WTAAxis_Phi/F");
  inTree_p->Branch("pt", pt, "pt[nParticle]/F");
  inTree_p->Branch("pmag", pmag, "pmag[nParticle]/F");
  inTree_p->Branch("eta", eta, "eta[nParticle]/F");
  inTree_p->Branch("theta", theta, "theta[nParticle]/F");
  inTree_p->Branch("phi", phi, "phi[nParticle]/F");
  inTree_p->Branch("pt_Perp", pt_Perp, "pt_Perp[nParticle]/F");
  inTree_p->Branch("pmag_Perp", pmag_Perp, "pmag_Perp[nParticle]/F");
  inTree_p->Branch("eta_Perp", eta_Perp, "eta_Perp[nParticle]/F");
  inTree_p->Branch("theta_Perp", theta_Perp, "theta_Perp[nParticle]/F");
  inTree_p->Branch("phi_Perp", phi_Perp, "phi_Perp[nParticle]/F");
  inTree_p->Branch("boostx", &boostx, "boostx/F");
  inTree_p->Branch("boosty", &boosty, "boosty/F");
  inTree_p->Branch("boostz", &boostz, "boostz/F");
  inTree_p->Branch("boost", &boost, "boost/F");
  inTree_p->Branch("particleWeight", &particleWeight, "particleWeight/F");

  return;
}

void boostedEvtData::preFillClean()
{
  for(Int_t pI = 0; pI < nParticle; ++pI){
    pt[pI] = reducedPrecision(pt[pI]);
    pmag[pI] = reducedPrecision(pmag[pI]);
    eta[pI] = reducedPrecision(eta[pI]);
    theta[pI] = reducedPrecision(theta[pI]);
    phi[pI] = reducedPrecision(phi[pI]);

    pt_Perp[pI] = reducedPrecision(pt_Perp[pI]);
    pmag_Perp[pI] = reducedPrecision(pmag_Perp[pI]);
    eta_Perp[pI] = reducedPrecision(eta_Perp[pI]);
    theta_Perp[pI] = reducedPrecision(theta_Perp[pI]);
    phi_Perp[pI] = reducedPrecision(phi_Perp[pI]);
  }

  return;
}

#endif
