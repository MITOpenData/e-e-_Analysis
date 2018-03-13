#ifndef PARTICLEDATA_H
#define PARTICLEDATA_H

#include "TTree.h"

class particleData{
 public:
  static const int nMaxPart = 1000;

  int nParticle;
  int EventNo;
  int RunNo;
  int year;
  int subDir;
  int process;
  float Energy;
  int bFlag;
  float particleWeight;
  float bx;
  float by;
  float ebx;
  float eby;
  float px[nMaxPart];
  float py[nMaxPart];
  float pz[nMaxPart];
  float pt[nMaxPart];
  float pmag[nMaxPart];
  float rap[nMaxPart];
  float eta[nMaxPart];
  float theta[nMaxPart];
  float phi[nMaxPart];
  float mass[nMaxPart];
  float charge[nMaxPart];
  // Starting from 0, pwflag (via Marcello) - CHARGED_TRACK, CHARGED_LEPTONS1, CHARGED_LEPTONS2, V0, PHOTON, NEUTRAL_HADRON
  int pwflag[nMaxPart];
  int pid[nMaxPart];
  float d0[nMaxPart];
  float z0[nMaxPart];
  int ntpc[nMaxPart];
  int nitc[nMaxPart];
  int nvdet[nMaxPart];
  float vx[nMaxPart];
  float vy[nMaxPart];
  float vz[nMaxPart];

  //thrust axis variables
  float pt_wrtThr[nMaxPart];
  float eta_wrtThr[nMaxPart];
  float rap_wrtThr[nMaxPart];
  float theta_wrtThr[nMaxPart];
  float phi_wrtThr[nMaxPart];
  float pt_wrtThrPerp[nMaxPart];
  float eta_wrtThrPerp[nMaxPart];
  float rap_wrtThrPerp[nMaxPart];
  float theta_wrtThrPerp[nMaxPart];
  float phi_wrtThrPerp[nMaxPart];
  float pt_wrtChThr[nMaxPart];
  float eta_wrtChThr[nMaxPart];
  float rap_wrtChThr[nMaxPart];
  float theta_wrtChThr[nMaxPart];
  float phi_wrtChThr[nMaxPart];
  float pt_wrtChThrPerp[nMaxPart];
  float eta_wrtChThrPerp[nMaxPart];
  float rap_wrtChThrPerp[nMaxPart];
  float theta_wrtChThrPerp[nMaxPart];
  float phi_wrtChThrPerp[nMaxPart];

  static const int nVar = 54;
  std::string varStr[nVar] = {"nParticle",
			      "EventNo",
			      "RunNo",
			      "year",
			      "subDir",
			      "process",
			      "Energy",
			      "bFlag",
                              "particleWeight",
			      "bx",
			      "by",
			      "ebx",
			      "eby",
			      "px",
			      "py",
			      "pz",
			      "pt",
			      "pmag",
			      "rap",
			      "eta",
			      "theta",
			      "phi",
			      "mass",
			      "charge",
			      "pwflag",
			      "pid",
			      "d0",
			      "z0",
			      "ntpc",
			      "nitc",
			      "nvdet",
			      "vx",
			      "vy",
			      "vz",
			      "pt_wrtThr",
			      "eta_wrtThr",
			      "rap_wrtThr",
			      "theta_wrtThr",
			      "phi_wrtThr",
			      "pt_wrtThrPerp",
			      "eta_wrtThrPerp",
			      "rap_wrtThrPerp",
			      "theta_wrtThrPerp",
			      "phi_wrtThrPerp",
			      "pt_wrtChThr",
			      "eta_wrtChThr",
			      "rap_wrtChThr",
			      "theta_wrtChThr",
			      "phi_wrtChThr",
			      "pt_wrtChThrPerp",
			      "eta_wrtChThrPerp",
			      "rap_wrtChThrPerp",
			      "theta_wrtChThrPerp",
			      "phi_wrtChThrPerp"};

  bool varIsGood[nVar];
  
  particleData();
  void SetStatusAndAddressRead(TTree* inTree_p, std::vector<std::string> inList);
  void SetBranchWrite(TTree* inTree_p);
};

particleData::particleData()
{
  //init values in event they are written to tree but are meaningless - treat -999 as danger
  nParticle = -999;
  EventNo = -999;
  RunNo = -999;
  year = -999;
  subDir = -999;
  process = -999;
  Energy = -999;
  bFlag = -999;
  particleWeight = 1;//keep this default as 1 (only used for mixed events, but should be 1 otherwise)
  bx = -999;
  by = -999;
  ebx = -999;
  eby = -999;

  for(Int_t i = 0; i < nMaxPart; ++i){
    px[i] = -999;
    py[i] = -999;
    pz[i] = -999;
    pt[i] = -999;
    pmag[i] = -999;
    rap[i] = -999;
    eta[i] = -999;
    theta[i] = -999;
    phi[i] = -999;
    mass[i] = -999;
    charge[i] = -999;
    pwflag[i] = -999;
    pid[i] = -999;
    d0[i] = -999;
    z0[i] = -999;
    ntpc[i] = -999;
    nitc[i] = -999;
    nvdet[i] = -999;
    vx[i] = -999.;
    vy[i] = -999.;
    vz[i] = -999.;
    
    pt_wrtThr[i] = -999;
    eta_wrtThr[i] = -999;
    rap_wrtThr[i] = -999;
    theta_wrtThr[i] = -999;
    phi_wrtThr[i] = -999;
    pt_wrtThrPerp[i] = -999;
    eta_wrtThrPerp[i] = -999;
    rap_wrtThrPerp[i] = -999;
    theta_wrtThrPerp[i] = -999;
    phi_wrtThrPerp[i] = -999;
    pt_wrtChThr[i] = -999;
    eta_wrtChThr[i] = -999;
    rap_wrtChThr[i] = -999;
    theta_wrtChThr[i] = -999;
    phi_wrtChThr[i] = -999;
    pt_wrtChThrPerp[i] = -999;
    eta_wrtChThrPerp[i] = -999;
    rap_wrtChThrPerp[i] = -999;
    theta_wrtChThrPerp[i] = -999;
    phi_wrtChThrPerp[i] = -999;
  }
  
  for(int i = 0; i < nVar; ++i){varIsGood[i] = true;}
  return;
}

void particleData::SetStatusAndAddressRead(TTree* inTree_p, std::vector<std::string> inList)
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
  if(varIsGood[1]) inTree_p->SetBranchAddress("EventNo", &EventNo);
  if(varIsGood[2]) inTree_p->SetBranchAddress("RunNo", &RunNo);
  if(varIsGood[3]) inTree_p->SetBranchAddress("year", &year);
  if(varIsGood[4]) inTree_p->SetBranchAddress("subDir", &subDir);
  if(varIsGood[5]) inTree_p->SetBranchAddress("process", &process);
  if(varIsGood[6]) inTree_p->SetBranchAddress("Energy", &Energy);
  if(varIsGood[7]) inTree_p->SetBranchAddress("bFlag", &bFlag);
  if(varIsGood[8]) inTree_p->SetBranchAddress("particleWeight", &particleWeight);
  if(varIsGood[9]) inTree_p->SetBranchAddress("bx", &bx);
  if(varIsGood[10]) inTree_p->SetBranchAddress("by", &by);
  if(varIsGood[11]) inTree_p->SetBranchAddress("ebx", &ebx);
  if(varIsGood[12]) inTree_p->SetBranchAddress("eby", &eby);
  if(varIsGood[13]) inTree_p->SetBranchAddress("px", px);
  if(varIsGood[14]) inTree_p->SetBranchAddress("py", py);
  if(varIsGood[15]) inTree_p->SetBranchAddress("pz", pz);
  if(varIsGood[16]) inTree_p->SetBranchAddress("pt", pt);
  if(varIsGood[17]) inTree_p->SetBranchAddress("pmag", pmag);
  if(varIsGood[18]) inTree_p->SetBranchAddress("rap", rap);
  if(varIsGood[19]) inTree_p->SetBranchAddress("eta", eta);
  if(varIsGood[20]) inTree_p->SetBranchAddress("theta", theta);
  if(varIsGood[21]) inTree_p->SetBranchAddress("phi", phi);
  if(varIsGood[22]) inTree_p->SetBranchAddress("mass", mass);
  if(varIsGood[23]) inTree_p->SetBranchAddress("charge", charge);
  if(varIsGood[24]) inTree_p->SetBranchAddress("pwflag", pwflag);
  if(varIsGood[25]) inTree_p->SetBranchAddress("pid", pid);
  if(varIsGood[26]) inTree_p->SetBranchAddress("d0", d0);
  if(varIsGood[27]) inTree_p->SetBranchAddress("z0", z0);
  if(varIsGood[28]) inTree_p->SetBranchAddress("ntpc", ntpc);
  if(varIsGood[29]) inTree_p->SetBranchAddress("nitc", nitc);
  if(varIsGood[30]) inTree_p->SetBranchAddress("nvdet", nvdet);
  if(varIsGood[31]) inTree_p->SetBranchAddress("vx", vx);
  if(varIsGood[32]) inTree_p->SetBranchAddress("vy", vy);
  if(varIsGood[33]) inTree_p->SetBranchAddress("vz", vz);
  if(varIsGood[34]) inTree_p->SetBranchAddress("pt_wrtThr", pt_wrtThr);
  if(varIsGood[35]) inTree_p->SetBranchAddress("eta_wrtThr", eta_wrtThr);
  if(varIsGood[36]) inTree_p->SetBranchAddress("rap_wrtThr", rap_wrtThr);
  if(varIsGood[37]) inTree_p->SetBranchAddress("theta_wrtThr", theta_wrtThr);
  if(varIsGood[38]) inTree_p->SetBranchAddress("phi_wrtThr", phi_wrtThr);
  if(varIsGood[39]) inTree_p->SetBranchAddress("pt_wrtThrPerp", pt_wrtThrPerp);
  if(varIsGood[40]) inTree_p->SetBranchAddress("eta_wrtThrPerp", eta_wrtThrPerp);
  if(varIsGood[41]) inTree_p->SetBranchAddress("rap_wrtThrPerp", rap_wrtThrPerp);
  if(varIsGood[42]) inTree_p->SetBranchAddress("theta_wrtThrPerp", theta_wrtThrPerp);
  if(varIsGood[43]) inTree_p->SetBranchAddress("phi_wrtThrPerp", phi_wrtThrPerp);
  if(varIsGood[44]) inTree_p->SetBranchAddress("pt_wrtChThr", pt_wrtChThr);
  if(varIsGood[45]) inTree_p->SetBranchAddress("eta_wrtChThr", eta_wrtChThr);
  if(varIsGood[46]) inTree_p->SetBranchAddress("rap_wrtChThr", rap_wrtChThr);
  if(varIsGood[47]) inTree_p->SetBranchAddress("theta_wrtChThr", theta_wrtChThr);
  if(varIsGood[48]) inTree_p->SetBranchAddress("phi_wrtChThr", phi_wrtChThr);
  if(varIsGood[49]) inTree_p->SetBranchAddress("pt_wrtChThrPerp", pt_wrtChThrPerp);
  if(varIsGood[50]) inTree_p->SetBranchAddress("eta_wrtChThrPerp", eta_wrtChThrPerp);
  if(varIsGood[51]) inTree_p->SetBranchAddress("rap_wrtChThrPerp", rap_wrtChThrPerp);
  if(varIsGood[52]) inTree_p->SetBranchAddress("theta_wrtChThrPerp", theta_wrtChThrPerp);
  if(varIsGood[53]) inTree_p->SetBranchAddress("phi_wrtChThrPerp", phi_wrtChThrPerp);
  
  return;
}

void particleData::SetBranchWrite(TTree* inTree_p)
{
  inTree_p->Branch("EventNo", &EventNo, "EventNo/I");
  inTree_p->Branch("RunNo", &RunNo, "RunNo/I");
  inTree_p->Branch("year", &year, "year/I");
  inTree_p->Branch("subDir", &subDir, "subDir/I");
  inTree_p->Branch("process", &process, "process/I");
  inTree_p->Branch("Energy", &Energy, "Energy/F");
  inTree_p->Branch("bFlag", &bFlag, "bFlag/I");
  inTree_p->Branch("particleWeight", &particleWeight, "particleWeight/F");
  inTree_p->Branch("bx", &bx, "bx/F");
  inTree_p->Branch("by", &by, "by/F");
  inTree_p->Branch("ebx", &ebx, "ebx/F");
  inTree_p->Branch("eby", &eby, "eby/F");
  inTree_p->Branch("nParticle", &nParticle, "nParticle/I"); 
  inTree_p->Branch("px", px, "px[nParticle]/F");
  inTree_p->Branch("py", py, "py[nParticle]/F");
  inTree_p->Branch("pz", pz, "pz[nParticle]/F");
  inTree_p->Branch("pt", pt, "pt[nParticle]/F");
  inTree_p->Branch("pmag", pmag, "pmag[nParticle]/F");
  inTree_p->Branch("rap", rap, "rap[nParticle]/F");
  inTree_p->Branch("eta", eta, "eta[nParticle]/F");
  inTree_p->Branch("theta", theta, "theta[nParticle]/F");
  inTree_p->Branch("phi", phi, "phi[nParticle]/F");
  inTree_p->Branch("mass", mass, "mass[nParticle]/F");
  inTree_p->Branch("charge", charge, "charge[nParticle]/F");
  inTree_p->Branch("pwflag", pwflag, "pwflag[nParticle]/I");
  inTree_p->Branch("pid", pid, "pid[nParticle]/I");
  inTree_p->Branch("d0", d0, "d0[nParticle]/F");
  inTree_p->Branch("z0", z0, "z0[nParticle]/F");
  inTree_p->Branch("ntpc", ntpc, "ntpc[nParticle]/I");
  inTree_p->Branch("nitc", nitc, "nitc[nParticle]/I");
  inTree_p->Branch("nvdet", nvdet, "nvdet[nParticle]/I");
  inTree_p->Branch("vx", vx, "vx[nParticle]/F");
  inTree_p->Branch("vy", vy, "vy[nParticle]/F");
  inTree_p->Branch("vz", vz, "vz[nParticle]/F");
  inTree_p->Branch("pt_wrtThr", pt_wrtThr, "pt_wrtThr[nParticle]/F");
  inTree_p->Branch("eta_wrtThr", eta_wrtThr, "eta_wrtThr[nParticle]/F");
  inTree_p->Branch("rap_wrtThr", rap_wrtThr, "rap_wrtThr[nParticle]/F");
  inTree_p->Branch("theta_wrtThr", theta_wrtThr, "theta_wrtThr[nParticle]/F");
  inTree_p->Branch("phi_wrtThr", phi_wrtThr, "phi_wrtThr[nParticle]/F");
  inTree_p->Branch("pt_wrtThrPerp", pt_wrtThrPerp, "pt_wrtThrPerp[nParticle]/F");
  inTree_p->Branch("eta_wrtThrPerp", eta_wrtThrPerp, "eta_wrtThrPerp[nParticle]/F");
  inTree_p->Branch("rap_wrtThrPerp", rap_wrtThrPerp, "rap_wrtThrPerp[nParticle]/F");
  inTree_p->Branch("theta_wrtThrPerp", theta_wrtThrPerp, "theta_wrtThrPerp[nParticle]/F");
  inTree_p->Branch("phi_wrtThrPerp", phi_wrtThrPerp, "phi_wrtThrPerp[nParticle]/F");
  inTree_p->Branch("pt_wrtChThr", pt_wrtChThr, "pt_wrtChThr[nParticle]/F");
  inTree_p->Branch("eta_wrtChThr", eta_wrtChThr, "eta_wrtChThr[nParticle]/F");
  inTree_p->Branch("rap_wrtChThr", rap_wrtChThr, "rap_wrtChThr[nParticle]/F");
  inTree_p->Branch("theta_wrtChThr", theta_wrtChThr, "theta_wrtChThr[nParticle]/F");
  inTree_p->Branch("phi_wrtChThr", phi_wrtChThr, "phi_wrtChThr[nParticle]/F");
  inTree_p->Branch("pt_wrtChThrPerp", pt_wrtChThrPerp, "pt_wrtChThrPerp[nParticle]/F");
  inTree_p->Branch("eta_wrtChThrPerp", eta_wrtChThrPerp, "eta_wrtChThrPerp[nParticle]/F");
  inTree_p->Branch("rap_wrtChThrPerp", rap_wrtChThrPerp, "rap_wrtChThrPerp[nParticle]/F");
  inTree_p->Branch("theta_wrtChThrPerp", theta_wrtChThrPerp, "theta_wrtChThrPerp[nParticle]/F");
  inTree_p->Branch("phi_wrtChThrPerp", phi_wrtChThrPerp, "phi_wrtChThrPerp[nParticle]/F");
  
  return;
}

#endif
