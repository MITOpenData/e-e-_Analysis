//c and cpp dependencies
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

//root dependencies
#include <TFile.h>
#include <TTree.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TNamed.h"

//fastjet dependencies
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

//local dependencies
#include "include/doLocalDebug.h"
#include "include/checkMakeDir.h"
#include "include/particleData.h"
#include "include/eventSelection.h"
#include "include/eventData.h"
#include "include/jetData.h"
#include "include/boostedEvtData.h"
#include "include/thrustTools.h"
#include "include/boostTools.h"
#include "include/sphericityTools.h"
#include "include/processJets.h"
#include "include/uniqueIDTools.h"
#include "src/mixFile.cc"

bool getIsMC(std::string inStr)
{
  if(inStr.find("MC") != std::string::npos && inStr.find("Data") == std::string::npos) return true;
  else return false;
}

bool getIsRecons(std::string inStr)
{
  if(!getIsMC(inStr)) return true;
  if(inStr.find("_recons_") != std::string::npos) return true;
  else return false;
}

//currently compatible with aleph filepaths
int yearFromPath(std::string inStr)
{
  const int nY=10;
  int years[nY] = {1991,1992,1993,1994,1995,1996,1997,1998,1999,2000};
  for(int i = 0; i < nY; ++i){
    if(inStr.find("/"+ std::to_string(years[i]) + "/") != std::string::npos) return years[i];
    if(inStr.find("Data"+ std::to_string(years[i])) != std::string::npos) return years[i];
  }
  std::cout << "Given string \'" << inStr << "\' does not contain valid year. return -1" << std::endl;
  return -1;
}

int subDirFromPath(std::string inStr, bool isMC)
{
  if(!isMC) return -999;
  const int nSubDir = 8;
  const int subDir[nSubDir] = {183, 189, 192, 196, 200, 202, 205, 207};
  for(int i = 0; i < nSubDir; ++i){
    if(inStr.find("/"+ std::to_string(subDir[i]) + "/") != std::string::npos) return subDir[i];
  }
  std::cout << "Given string \'" << inStr << "\' does not contain valid subdir. return -1" << std::endl;
  return -1;
}

//implemented as according to here https://www.dropbox.com/s/zujvnpftoqb7w6k/MC_alephstatus.pdf?dl=0, s3, with mod for GGTT and GGUS
int processFromPath(std::string inStr)
{
  const int nProc = 13;
  std::string proc[nProc] = {"GGBB", "GGCC", "GGSS", "GGTT", "GGUD", "KQQ", "KWENU", "KWW4F", "PZEE", "PZZ", "TT", "ZNN", "GGUS"};
  for(int i = 0; i < nProc; ++i){
    if(inStr.find("/" + proc[i] + "/") != std::string::npos) return i;
  }
  return -1;
}

std::vector<std::string> processAlephString(std::string inStr)
{
  bool addDummy = false;
  bool addDummy2 = false;
  if(inStr.find("ALEPH_DATA") != std::string::npos) addDummy = true;
  if(inStr.find("MC_RECO") != std::string::npos) addDummy = true;
  if(inStr.find("MC_TRUE_AFTER_CUT") != std::string::npos) addDummy = true;
  if(inStr.find("MC_TRUE_BEFORE_CUT") != std::string::npos) addDummy = true;
  if(inStr.find("Primary vertex info") != std::string::npos) addDummy2 = true;

  const Int_t nSpecialRemoveStr = 2;
  std::string specialRemoveStr[nSpecialRemoveStr] = {"d0", "z0"};
  for(Int_t rI = 0; rI < nSpecialRemoveStr; ++rI){
    while(inStr.find(specialRemoveStr[rI]) != std::string::npos) inStr.replace(inStr.find(specialRemoveStr[rI]), specialRemoveStr[rI].size(), "");
  }

  const std::string goodChar = "0123456789-. ";
  unsigned int pos = 0;
  while(inStr.size() > pos){
    if(goodChar.find(inStr.substr(pos,1)) == std::string::npos) inStr.replace(pos,1,"");
    else ++pos;
  }

  std::vector<std::string> retV;
  if(addDummy){retV.push_back("-999."); retV.push_back("-999."); retV.push_back("-999."); retV.push_back("-999.");}
  if(addDummy2){retV.push_back("-998."); retV.push_back("-998."); retV.push_back("-998."); retV.push_back("-998.");}

  while(inStr.size() != 0){
    while(inStr.substr(0,1).find(" ") != std::string::npos){inStr.replace(0,1,"");}

    if(inStr.find(" ") != std::string::npos){
      retV.push_back(inStr.substr(0,inStr.find(" ")));
      inStr.replace(0,inStr.find(" ")+1,"");
    }
    else if(inStr.size() != 0){
      retV.push_back(inStr);
      inStr = "";
    }
  }

  return retV;
}


bool checkGeneral(std::string inStr, std::string checkStr)
{
  if(inStr.find(checkStr) != std::string::npos && inStr.size() == checkStr.size()) return true;
  return false;
}
bool check999(std::string inStr){return checkGeneral(inStr, "-999.");}
bool check998(std::string inStr){return checkGeneral(inStr, "-998.");}

//isNewInfo was originally for introduction of ntpc, etc.
//isNewInfo2 is for nitc, nvdet

void initVal(Int_t* in, Int_t set){(*in) = set; return;}
void initVal(Float_t* in, Float_t set){(*in) = set; return;}
void initVal(Double_t* in, Double_t set){(*in) = set; return;}
void initVal(std::vector<Int_t*> in, Int_t set)
{
  for(unsigned int i = 0; i < in.size(); ++i){
    if(in.at(i) == NULL) continue;
    initVal(in.at(i), set);
  }
  return;
}
void initVal(std::vector<Float_t*> in, Float_t set)
{
  for(unsigned int i = 0; i < in.size(); ++i){
    if(in.at(i) == NULL) continue;
    initVal(in.at(i), set);
  }
  return;
}
void initVal(std::vector<Double_t*> in, Float_t set)
{
  for(unsigned int i = 0; i < in.size(); ++i){
    if(in.at(i) == NULL) continue;
    initVal(in.at(i), set);
  }
  return;
}
void initTVector3(TVector3* in){(*in) = TVector3(0,0,0); return;}
void initTVector3(std::vector<TVector3*> in)
{
  for(unsigned int i = 0; i < in.size(); ++i){initTVector3(in.at(i));}
  return;
}

void doEndEvent(particleData* pData_p, eventData* eData_p, std::vector<boostedEvtData*> bData_p, Int_t counterParticles, Int_t nTrk, Int_t nTrk_GT0p4, Int_t nTrk_GT0p4Thrust, TVector3 netP, TVector3 netP_charged)
{
  pData_p->nParticle=counterParticles;
  
  for(unsigned int jI = 0; jI < bData_p.size(); ++jI){bData_p.at(jI)->nParticle = counterParticles;}

  TVector3 thrust, thrust_charged;
  initTVector3({&thrust, &thrust_charged});
  
  thrust = getThrust(pData_p->nParticle, pData_p->px, pData_p->py, pData_p->pz, THRUST::OPTIMAL);
  thrust_charged = getChargedThrust(pData_p->nParticle, pData_p->px, pData_p->py, pData_p->pz, pData_p->pwflag, THRUST::OPTIMAL);
  eData_p->nChargedHadrons = nTrk;
  eData_p->nChargedHadrons_GT0p4 = nTrk_GT0p4;
  eData_p->nChargedHadrons_GT0p4Thrust = nTrk_GT0p4Thrust;
  setThrustVariables(pData_p, eData_p, thrust, thrust_charged);
  eData_p->Thrust = thrust.Mag();
  eData_p->TTheta = thrust.Theta();
  eData_p->TPhi = thrust.Phi();
  eData_p->Thrust_charged = thrust_charged.Mag();
  eData_p->TTheta_charged = thrust_charged.Theta();
  eData_p->TPhi_charged = thrust_charged.Phi();
  eData_p->missP = netP.Mag();
  eData_p->missPt = netP.Perp();
  eData_p->missTheta = netP.Theta();
  eData_p->missPhi = netP.Phi();
  eData_p->missChargedP = netP_charged.Mag();
  eData_p->missChargedPt = netP_charged.Perp();
  eData_p->missChargedTheta = netP_charged.Theta();
  eData_p->missChargedPhi = netP_charged.Phi();

  Sphericity spher = Sphericity(pData_p->nParticle, pData_p->px, pData_p->py, pData_p->pz, pData_p->pwflag, true);
  spher.setTree(eData_p);

  eventSelection eSelection;
  eSelection.setEventSelection(pData_p, eData_p);

  eData_p->passesTotalChgEnergyMin = eSelection.getPassesTotalChgEnergyMin();
  eData_p->passesNTrkMin = eSelection.getPassesNTrkMin();
  eData_p->passesSTheta = eSelection.getPassesSTheta();
  eData_p->passesMissP = eSelection.getPassesMissP();

  eData_p->passesISR = eSelection.getPassesISR();
  eData_p->passesWW = eSelection.getPassesWW();

  eData_p->passesAll = (eData_p->passesNTupleAfterCut && eData_p->passesTotalChgEnergyMin && eData_p->passesNTrkMin && eData_p->passesSTheta && eData_p->passesMissP && eData_p->passesISR && eData_p->passesWW);

  return;
}

int scan(const std::string inFileName, const bool isNewInfo, const bool isNewInfo2, std::string outFileName="")
{
  std::cout << "Gets Here" << std::endl;

  if(!checkFile(inFileName)){
    std::cout << "Given inFileName \'" << inFileName << "\' is invalid, return 1." << std::endl;
    return 1;
  }

  uniqueIDTools idMaker;

  std::vector<std::string> fileList;
  if(inFileName.size() < 5){
    std::cout << "Given inFileName \'" << inFileName << "\' is invalid, return 1." << std::endl;
    return 1;
  }
  else if(inFileName.substr(inFileName.size()-4, 4).find(".txt") != std::string::npos){
    std::ifstream pathFile(inFileName.c_str());
    std::string paths;
    while(std::getline(pathFile, paths)){
      if(paths.size() == 0) continue;
      fileList.push_back(paths);
    }
    pathFile.close();
  }
  else if(inFileName.substr(inFileName.size()-5, 5).find(".list") != std::string::npos) fileList.push_back(inFileName);
  else if(inFileName.substr(inFileName.size()-6, 6).find(".aleph") != std::string::npos) fileList.push_back(inFileName);
  else{
    std::cout << "Given inFileName \'" << inFileName << "\' is invalid, return 1." << std::endl;
    return 1;
  }

  if(fileList.size() == 0){
    std::cout << "File list empty for input \'" << inFileName << "\'. return 1" << std::endl;
    return 1;
  }

  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  const bool isMC = getIsMC(fileList.at(0));
  const bool isRecons = getIsRecons(fileList.at(0));

  std::cout << "IsMC, recons: " << isMC << ", " << isRecons << std::endl;

  if(outFileName.size() == 0){
    outFileName = inFileName;
    if(outFileName.rfind(".") != std::string::npos) outFileName = outFileName.substr(0,outFileName.rfind("."));
    outFileName = outFileName + ".root";
    while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1,"");}
  }

  //define jetMaker here so rParam can be used in jetTree name for clarity
  const double jtPtCut = .01;
  const int nJtAlgo = 8;
  const double rParam[nJtAlgo] = {0.4, 0.4, 0.8, 0.8, -1., -1., -1., -1.};
  bool rParamIs8[nJtAlgo];
  for(int i = 0; i < nJtAlgo; ++i){rParamIs8[i] = rParam[i] > .799 && rParam[i] < .801;}
  const int nFinalClust[nJtAlgo] = {-1, -1, -1, -1, 2, 2, 3, 3};
  const double recombScheme[nJtAlgo] = {fastjet::E_scheme, fastjet::WTA_modp_scheme, fastjet::E_scheme, fastjet::WTA_modp_scheme, fastjet::E_scheme, fastjet::WTA_modp_scheme, fastjet::E_scheme, fastjet::WTA_modp_scheme};
  fastjet::JetDefinition jDef[nJtAlgo];
  fastjet::JetDefinition jDefReclust[nJtAlgo];
  
  for(int i = 0; i < nJtAlgo; ++i){
    if(rParam[i] > 0){
      jDef[i] = fastjet::JetDefinition(fastjet::ee_genkt_algorithm, rParam[i], -1, fastjet::RecombinationScheme(recombScheme[i]));
      jDefReclust[i] = fastjet::JetDefinition(fastjet::ee_genkt_algorithm, 5, -1, fastjet::RecombinationScheme(recombScheme[i]));
    }
    else{
      jDef[i] = fastjet::JetDefinition(fastjet::ee_kt_algorithm, fastjet::RecombinationScheme(recombScheme[i]));
      jDefReclust[i] = fastjet::JetDefinition(fastjet::ee_kt_algorithm, fastjet::RecombinationScheme(recombScheme[i]));
    }
  }

  const std::string partTreeName = "t";
  const std::string genPartTreeName = "tgen";
  const std::string genPartTreeBeforeName = "tgenBefore";
  std::string jetTreeName[nJtAlgo];
  std::string genJetTreeName[nJtAlgo];
  std::string genJetTreeBeforeName[nJtAlgo];
  std::string boostedTreeName[nJtAlgo];
  std::string genboostedTreeName[nJtAlgo];
  std::string genboostedTreeBeforeName[nJtAlgo];

  for(int i = 0; i < nJtAlgo; ++i){
    std::string recombSchemeStr = "EScheme";
    if(recombScheme[i] == fastjet::WTA_modp_scheme) recombSchemeStr = "WTAmodpScheme";

    std::string rOrJetStr = "akR" + std::to_string(int(rParam[i]*10));
    if(rParam[i] < 0) rOrJetStr = "ktN" + std::to_string(nFinalClust[i]);

    jetTreeName[i] = rOrJetStr + recombSchemeStr + "JetTree";
    genJetTreeName[i] = rOrJetStr + recombSchemeStr + "GenJetTree";
    genJetTreeBeforeName[i] = rOrJetStr + recombSchemeStr + "GenJetBeforeTree";

    boostedTreeName[i] = "Boosted" + recombSchemeStr + rOrJetStr + "Evt";
    genboostedTreeName[i] = "genBoosted" + recombSchemeStr + rOrJetStr + "Evt";
    genboostedTreeBeforeName[i] = "genBoosted" + recombSchemeStr + rOrJetStr + "BeforeEvt";
    if(recombSchemeStr.find("modpScheme") != std::string::npos){
      boostedTreeName[i].replace(boostedTreeName[i].find("modpScheme"), std::string("modpScheme").size(), "");
      genboostedTreeName[i].replace(genboostedTreeName[i].find("modpScheme"), std::string("modpScheme").size(), "");
      genboostedTreeBeforeName[i].replace(genboostedTreeBeforeName[i].find("modpScheme"), std::string("modpScheme").size(), "");
    }
    if(rOrJetStr.find("ak") != std::string::npos){
      boostedTreeName[i].replace(boostedTreeName[i].find("ak"), std::string("ak").size(), "");
      genboostedTreeName[i].replace(genboostedTreeName[i].find("ak"), std::string("ak").size(), "");
      genboostedTreeBeforeName[i].replace(genboostedTreeBeforeName[i].find("ak"), std::string("ak").size(), "");
    }
  }

  std::string finalPartTreeName = partTreeName;
  if(!isRecons) finalPartTreeName = genPartTreeName;

  std::string finalJetTreeName[nJtAlgo];
  for(int i = 0; i < nJtAlgo; ++i){
    finalJetTreeName[i] = jetTreeName[i];
    if(!isRecons) finalJetTreeName[i] = genJetTreeName[i];
  }

  TFile *hf = new TFile(outFileName.c_str(), "RECREATE");
  TTree *tout = new TTree(finalPartTreeName.c_str(), finalPartTreeName.c_str());
  TTree *bout[nJtAlgo];
  TTree *jout[nJtAlgo];
  for(int i = 0; i < nJtAlgo; ++i){
    jout[i] = NULL;
    jout[i] = new TTree(finalJetTreeName[i].c_str(), finalJetTreeName[i].c_str());
    
    bout[i] = NULL;
    if(recombScheme[i] == fastjet::WTA_modp_scheme){     
      if(rParamIs8[i] || nFinalClust[i] == 2) bout[i] = new TTree(boostedTreeName[i].c_str(), boostedTreeName[i].c_str());
    }
  }

  TTree *tgout = 0;
  TTree *bgout[nJtAlgo];
  TTree *jgout[nJtAlgo];

  TTree *tgBeforeout = 0;
  TTree *bgBeforeout[nJtAlgo];
  TTree *jgBeforeout[nJtAlgo];

  if(isRecons && isMC){
    tgout = new TTree(genPartTreeName.c_str(), genPartTreeName.c_str());
    tgBeforeout = new TTree(genPartTreeBeforeName.c_str(), genPartTreeBeforeName.c_str());

    for(int i = 0; i < nJtAlgo; ++i){
      jgout[i] = NULL;
      jgout[i] = new TTree(genJetTreeName[i].c_str(), genJetTreeName[i].c_str());

      jgBeforeout[i] = NULL;
      jgBeforeout[i] = new TTree(genJetTreeBeforeName[i].c_str(), genJetTreeBeforeName[i].c_str());

      bgout[i] = NULL;
      bgBeforeout[i] = NULL;
      if(recombScheme[i] == fastjet::WTA_modp_scheme){
	if(rParamIs8[i] || nFinalClust[i] == 2){
	  bgout[i] = new TTree(genboostedTreeName[i].c_str(), genboostedTreeName[i].c_str());
	  bgBeforeout[i] = new TTree(genboostedTreeBeforeName[i].c_str(), genboostedTreeBeforeName[i].c_str());
	}
      }
    }
  }

  eventSelection eSelection;

  particleData pData;
  jetData jData[nJtAlgo];
  eventData eData;
  boostedEvtData bData[nJtAlgo];

  eventSelection egSelection;
  particleData pgData;
  jetData jgData[nJtAlgo];
  eventData egData;
  boostedEvtData bgData[nJtAlgo];

  particleData pgBeforeData;
  jetData jgBeforeData[nJtAlgo];
  eventData egBeforeData;
  boostedEvtData bgBeforeData[nJtAlgo];

  pData.SetBranchWrite(tout);
  eData.SetBranchWrite(tout);
  //thrust quantities

  for(int i = 0; i < nJtAlgo; ++i){
    jData[i].SetBranchWrite(jout[i]);
    if(recombScheme[i] == fastjet::WTA_modp_scheme){
      if(rParamIs8[i] || nFinalClust[i] == 2) bData[i].SetBranchWrite(bout[i]);
    }
  }

  if(isRecons && isMC){
    pgData.SetBranchWrite(tgout);
    egData.SetBranchWrite(tgout);

    pgBeforeData.SetBranchWrite(tgBeforeout);
    egBeforeData.SetBranchWrite(tgBeforeout);

    for(int i = 0; i < nJtAlgo; ++i){
      jgData[i].SetBranchWrite(jgout[i]);
      jgBeforeData[i].SetBranchWrite(jgBeforeout[i]);
      if(recombScheme[i] == fastjet::WTA_modp_scheme){
	if(rParamIs8[i] || nFinalClust[i] == 2){bgData[i].SetBranchWrite(bgout[i]); bgBeforeData[i].SetBranchWrite(bgBeforeout[i]);}
      }
    }
  }

  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    std::cout << "Processing " << fI << "/" << fileList.size() << ", \'" << fileList.at(fI) << "\'" << std::endl;
    const int year = yearFromPath(fileList.at(fI));
    const int subDir = subDirFromPath(fileList.at(fI), isMC);
    const int process = processFromPath(fileList.at(fI));

    std::ifstream file(Form("%s",fileList.at(fI).c_str()));
    std::string getStr;
    
    int counterEntries=0;
    int counterParticles=0;
    TLorentzVector v;
    std::vector<fastjet::PseudoJet> particles;

    std::vector<int> runNo;
    std::vector<int> evtNo;
    TVector3 netP, netP_charged;
    initTVector3({&netP, &netP_charged});

    int nTrk=0;
    int nTrk_GT0p4=0;
    int nTrk_GT0p4Thrust=0;

    while(std::getline(file,getStr)){
      if(getStr.size() == 0) continue;
      else if(getStr.find("******") != std::string::npos) continue;
      else if(getStr.find("END_EVENT") != std::string::npos) continue; 
      else if(getStr.find("END_FILE") != std::string::npos) continue;
      std::vector<std::string> num = processAlephString(getStr);
      
      if(num.size() < 5){
	std::cout << "Number of columns for line \'" << getStr << "\' is invalid, size \'" << num.size() << "\'. return 1" << std::endl;
	//gotta cleanup before return
	delete tout;
	for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	  delete jout[jIter];
	  if(recombScheme[jIter] == fastjet::WTA_modp_scheme){
	    if(rParamIs8[jIter] || nFinalClust[jIter] == 2) delete bout[jIter];
	  }
	}
	
	if(isMC && isRecons){
	  delete tgout;
	  delete tgBeforeout;
	  for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	    delete jgout[jIter];
	    delete jgBeforeout[jIter];
	    if(recombScheme[jIter] == fastjet::WTA_modp_scheme){
	      if(rParamIs8[jIter] || nFinalClust[jIter] == 2){delete bgout[jIter]; delete bgBeforeout[jIter];}
	    }
	  }
	}
	
	hf->Close();
	delete hf;
	
	return 1;
      }
    
      if(check999(num.at(0)) && check999(num.at(1)) && check999(num.at(2))/* not all files do four 999 && check999(num.at(3)*/){
	std::vector<boostedEvtData*> bDataTempVect;
	for(Int_t jI = 0; jI < nJtAlgo; ++jI){bDataTempVect.push_back(&(bData[jI]));}
	doEndEvent(&pData, &eData, bDataTempVect, counterParticles, nTrk, nTrk_GT0p4, nTrk_GT0p4Thrust, netP, netP_charged);
	bDataTempVect.clear();

	if(counterEntries>0){pData.preFillClean(); tout->Fill();}
	//Processing particles->jets
	for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	  processJets(particles, jDef[jIter], jDefReclust[jIter], &(jData[jIter]), jtPtCut, rParam[jIter], nFinalClust[jIter]);
	}

	if(counterEntries>0){
	  for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	    jData[jIter].preFillClean();
	    jout[jIter]->Fill();
	    if(recombScheme[jIter]==fastjet::WTA_modp_scheme){
	      if(rParamIs8[jIter] || nFinalClust[jIter] == 2){
		if(jData[jIter].nref<2) setBoostedVariables(false, &pData, &(bData[jIter]));
		else{
		  TVector3 wtaBoost = findBack2BackBoost(jData[jIter].fourJet[0],jData[jIter].fourJet[1]);
		  setBoostedVariables(true, &pData, &(bData[jIter]), jData[jIter].fourJet[0], wtaBoost);
		}
		bData[jIter].preFillClean();
		bout[jIter]->Fill();
	      }
            }
	  }
	}      

	//clear particles for next iteration clustering
	particles.clear();

	unsigned int runPos = 3;
	unsigned int evtPos = 4;
	unsigned int ePos = 5;

	//again, 4 999 is max	if(check999(num.at(4))){runPos++; evtPos++; ePos++;}
	if(check999(num.at(3))){runPos++; evtPos++; ePos++;}

	pData.year = year;
	pData.subDir = subDir;
	pData.process = process;
	pData.isMC = isMC;

	pData.RunNo = std::stoi(num.at(runPos));
	pData.EventNo= std::stoi(num.at(evtPos));
	pData.Energy= std::stof(num.at(ePos));

	int tempSubDir = subDir;
	if(tempSubDir < 0) tempSubDir = 0;
	int tempProcess = process;
	if(tempProcess < 0) tempProcess = 15;
	 
	pData.uniqueID = idMaker.getID(isMC, tempProcess, tempSubDir, year, pData.RunNo, pData.EventNo);

	eData.passesNTupleAfterCut = true;

	initVal(&(pData.bFlag), -999.);
	initVal({&(pData.bx), &(pData.by), &(pData.ebx), &(pData.eby)}, -999.);

	runNo.push_back(pData.RunNo);
	evtNo.push_back(pData.EventNo);
      
	initTVector3({&netP, &netP_charged});
	initVal({&nTrk, &nTrk_GT0p4, &nTrk_GT0p4Thrust, &counterParticles}, 0);

	++counterEntries;	
	
	continue;
      }
      else if(check998(num.at(0)) && check998(num.at(1)) && check998(num.at(2))){
	pData.bFlag = std::stoi(num.at(4));
	pData.bx = std::stof(num.at(5));
	pData.by = std::stof(num.at(6));
	pData.ebx = std::stof(num.at(7));
	pData.eby = std::stof(num.at(8));
	
	continue;
      }

      
      // check the number of columns before assigning values
      bool assumePID = false;
      if(num.size() == 6 && !isNewInfo) assumePID = false; 
      else if(num.size() == 9 && isNewInfo) assumePID = false; 
      else if(num.size() == 11 && isNewInfo2) assumePID = false; 
      else if(!isNewInfo && (num.size() == 7 || num.size() == 8)) assumePID = true; 
      else if(isNewInfo && (num.size() == 10 || num.size() == 11)) assumePID = true; 
      else if(isNewInfo2 && (num.size() == 12 || num.size() == 13)) assumePID = true; 
      else{//return, this is an invalid format (or fix code here if format valid
	std::cout << "Number of columns for line \'" << getStr << "\' is invalid, size \'" << num.size() << "\'. return 1" << std::endl;
	//gotta cleanup before return
	delete tout;
	for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	  delete jout[jIter];
	  if(recombScheme[jIter]==fastjet::WTA_modp_scheme){
	    if(rParamIs8[jIter] || nFinalClust[jIter] == 2) delete bout[jIter];
	  }
	}

	if(isMC && isRecons){
	  delete tgout;
	  delete tgBeforeout;
	  for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	    delete jgout[jIter];
	    delete jgBeforeout[jIter];
	    if(recombScheme[jIter]==fastjet::WTA_modp_scheme){
	      if(rParamIs8[jIter] || nFinalClust[jIter] == 2){delete bgout[jIter]; delete bgBeforeout[jIter];}
	    }
	  }
	}
	
	hf->Close();
	delete hf;
	
	return 1;
      }
      
      
      float _px = std::stof(num.at(0));
      float _py = std::stof(num.at(1));
      float _pz = std::stof(num.at(2));
      float _m = std::stof(num.at(3));
      float _charge = std::stof(num.at(4));
      int _pwflag = std::stoi(num.at(5));

      if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
      
      pData.px[counterParticles]=_px;
      pData.py[counterParticles]=_py;
      pData.pz[counterParticles]=_pz;
      pData.mass[counterParticles]=_m;
      v.SetXYZM(_px,_py,_pz,_m);
      fastjet::PseudoJet particle(_px,_py,_pz,v.E());
      particle.set_user_index(_pwflag);
      particles.push_back(particle);
      pData.pt[counterParticles]=v.Pt();
      pData.pmag[counterParticles]=v.Rho(); //Added later on
      pData.rap[counterParticles]=v.Rapidity();
      pData.eta[counterParticles]=v.PseudoRapidity();
      pData.theta[counterParticles]=v.Theta();
      pData.phi[counterParticles]=v.Phi();
      
      pData.charge[counterParticles]=_charge;
      pData.pwflag[counterParticles]=_pwflag;
      //check before assigning PID
      if(assumePID) pData.pid[counterParticles]=std::stoi(num.at(6));
      else pData.pid[counterParticles]=-999;
     
      if((isNewInfo || isNewInfo2) && _pwflag < 4){
	int posOffset = 0;
	if(assumePID) posOffset = 1;

	pData.d0[counterParticles] = std::stof(num.at(6 + posOffset));
	pData.z0[counterParticles] = std::stof(num.at(7 + posOffset));
	pData.ntpc[counterParticles] =  std::stoi(num.at(8 + posOffset));	

	if(isNewInfo2){
	  pData.nitc[counterParticles] =  std::stoi(num.at(9 + posOffset));	
	  pData.nvdet[counterParticles] =  std::stoi(num.at(10 + posOffset));	
	}
	else initVal({&(pData.nitc[counterParticles]), &(pData.nvdet[counterParticles])}, -127);
      }
      else{
	initVal({&(pData.d0[counterParticles]), &(pData.z0[counterParticles]), NULL}, -999);
	initVal({&(pData.ntpc[counterParticles]), &(pData.nitc[counterParticles]), &(pData.nvdet[counterParticles])}, -127);
      }
      initVal({&(pData.vx[counterParticles]), &(pData.vy[counterParticles]), &(pData.vz[counterParticles])}, -999.);
      //missing momentum calculation and multiplicity calculation
      netP -= TVector3(_px,_py,_pz);
      if(_pwflag==0){
        netP_charged -= TVector3(_px,_py,_pz);
        nTrk++;
        if(v.Pt()>0.4) nTrk_GT0p4++; 
        if(v.Pt()>0.4) nTrk_GT0p4Thrust++; 
      } 
 
      ++counterParticles;	
    }

    if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  
    //Have to fill one last time since the condition for fill is dependent on NEXT EVENT existing, else we would lose last event per file
    std::vector<boostedEvtData*> bDataTempVect;
    for(Int_t jI = 0; jI < nJtAlgo; ++jI){bDataTempVect.push_back(&(bData[jI]));}
    doEndEvent(&pData, &eData, bDataTempVect, counterParticles, nTrk, nTrk_GT0p4, nTrk_GT0p4Thrust, netP, netP_charged);
    bDataTempVect.clear();
    
    if(counterEntries>0){pData.preFillClean(); tout->Fill();}
    for(int jIter = 0; jIter < nJtAlgo; ++jIter){
      processJets(particles, jDef[jIter], jDefReclust[jIter], &(jData[jIter]), jtPtCut, rParam[jIter], nFinalClust[jIter]);
    }
    if(counterEntries>0){
      for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	jData[jIter].preFillClean();
	jout[jIter]->Fill();
	if(recombScheme[jIter]==fastjet::WTA_modp_scheme){
	  if(rParamIs8[jIter] || nFinalClust[jIter] == 2){
	    if(jData[jIter].nref<2) setBoostedVariables(false, &pData, &(bData[jIter]));
	    else{
	      TVector3 wtaBoost = findBack2BackBoost(jData[jIter].fourJet[0],jData[jIter].fourJet[1]);
	      setBoostedVariables(true, &pData, &(bData[jIter]), jData[jIter].fourJet[0], wtaBoost);
	    }
	    bData[jIter].preFillClean();
	    bout[jIter]->Fill();
	  }
        }
      }
    }
    file.close();
  
    if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    if(isRecons && isMC){
      std::string genFileStr = fileList.at(fI);
      genFileStr.replace(genFileStr.find("_recons_"), 8, "_mctrue_");

      if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      std::ifstream fileGen(genFileStr.c_str());
      initTVector3({&netP, &netP_charged});
      initVal({&nTrk, &nTrk_GT0p4, &nTrk_GT0p4Thrust, &counterParticles, &counterEntries}, 0);

      while(std::getline(fileGen,getStr)){
	if(getStr.size() == 0) continue;
	else if(getStr.find("******") != std::string::npos) continue;
	else if(getStr.find("END_EVENT") != std::string::npos) continue;
	else if(getStr.find("END_FILE") != std::string::npos) continue;
	std::vector<std::string> num = processAlephString(getStr);

	if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	if(num.size() < 5){
	  std::cout << "Number of columns for line \'" << getStr << "\' is invalid, size \'" << num.size() << "\'. return 1" << std::endl;
	  //gotta cleanup before return
	  delete tout;
	  for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	    delete jout[jIter];
	    if(recombScheme[jIter]==fastjet::WTA_modp_scheme){
	      if(rParamIs8[jIter] || nFinalClust[jIter] == 2) delete bout[jIter];
	    }
	  }

	  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	  if(isMC && isRecons){
	    delete tgout;
	    delete tgBeforeout;
	    for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	      delete jgout[jIter];
	      if(recombScheme[jIter]==fastjet::WTA_modp_scheme){
		if(rParamIs8[jIter] || nFinalClust[jIter] == 2){delete bgout[jIter]; delete bgBeforeout[jIter];}
	      }
	    }
	  }

	  hf->Close();
	  delete hf;
	
	  return 1;
	}

	if(check999(num.at(0)) && check999(num.at(1)) && check999(num.at(2)) /* && check999(num.at(3))*/){ 
	  std::vector<boostedEvtData*> bgDataTempVect;
	  for(Int_t jI = 0; jI < nJtAlgo; ++jI){bgDataTempVect.push_back(&(bgData[jI]));}
	  doEndEvent(&pgData, &egData, bgDataTempVect, counterParticles, nTrk, nTrk_GT0p4, nTrk_GT0p4Thrust, netP, netP_charged);
	  bgDataTempVect.clear();

	  if(counterEntries>0){pgData.preFillClean(); tgout->Fill();}

	  //Processing particles->jets
	  for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	    processJets(particles, jDef[jIter], jDefReclust[jIter], &(jgData[jIter]), jtPtCut, rParam[jIter], nFinalClust[jIter]);
	    if(counterEntries>0){
	      jgData[jIter].preFillClean();
	      jgout[jIter]->Fill();
	      if(recombScheme[jIter]==fastjet::WTA_modp_scheme){
		if(rParamIs8[jIter] || nFinalClust[jIter] == 2){
		  if(jgData[jIter].nref<2) setBoostedVariables(false, &pgData, &(bgData[jIter]));
		  else{
		    TVector3 wtaBoost = findBack2BackBoost(jgData[jIter].fourJet[0],jgData[jIter].fourJet[1]);
		    setBoostedVariables(true, &pgData, &(bgData[jIter]), jgData[jIter].fourJet[0], wtaBoost);
		  }
		  bgData[jIter].preFillClean();
		  bgout[jIter]->Fill();
		}
	      }		
	    }
	  }
	  //clear particles for next iteration clustering
	  particles.clear();

	  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	  unsigned int runPos = 3;
	  unsigned int evtPos = 4;
	  unsigned int ePos = 5;
	  
	  if(check999(num.at(3))){runPos++; evtPos++; ePos++;}
	  //	  if(check999(num.at(4))){runPos++; evtPos++; ePos++;}

	  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	
	  pgData.year = year;
	  pgData.subDir = subDir;
	  pgData.process = process;
	  pgData.isMC = isMC;
	  pgData.RunNo = std::stoi(num.at(runPos));
	  pgData.EventNo= std::stoi(num.at(evtPos));
	  pgData.Energy= std::stof(num.at(ePos));
       
	  int tempSubDir = subDir;
	  if(tempSubDir < 0) tempSubDir = 0;
	  int tempProcess = process;
	  if(tempProcess < 0) tempProcess = 15;
	  
	  pgData.uniqueID = idMaker.getID(isMC, tempProcess, tempSubDir, year, pgData.RunNo, pgData.EventNo);

	  egData.passesNTupleAfterCut = true;


	  initVal(&(pData.bFlag), -999.);
	  initVal({&(pData.bx), &(pData.by), &(pData.ebx), &(pData.eby)}, -999.);	 
	  initTVector3({&netP, &netP_charged});
	  initVal({&nTrk, &nTrk_GT0p4, &nTrk_GT0p4Thrust, &counterParticles}, 0);
	  
	  if(pgData.RunNo != runNo.at(counterEntries) || pgData.EventNo != evtNo.at(counterEntries)){
	    std::cout << "Gen entries dont match reco for file \'" << genFileStr << "\'. return 1" << std::endl;
	    //gotta cleanup before return
	    delete tout;

	    for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	      delete jout[jIter];
	      
	      if(recombScheme[jIter]==fastjet::WTA_modp_scheme){
		if(rParamIs8[jIter] || nFinalClust[jIter] == 2) delete bout[jIter];
	      }
	    }

	    if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	    if(isMC && isRecons){
	      delete tgout;
	      delete tgBeforeout;
	      for(int jIter = 0; jIter < nJtAlgo; ++jIter){
		delete jgout[jIter];
		delete jgBeforeout[jIter];
		if(recombScheme[jIter]==fastjet::WTA_modp_scheme){
		  if(rParamIs8[jIter] || nFinalClust[jIter] == 2){delete bgout[jIter]; delete bgBeforeout[jIter];}
		}
	      }
	    }

	    if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	    hf->Close();
	    delete hf;
	    
	    return 1;
	  }
	  
	  ++counterEntries;	

	  continue;
	}
	else if(check998(num.at(0)) && check998(num.at(1)) && check998(num.at(2))){
	  pgData.bFlag = std::stoi(num.at(4));
	  pgData.bx = std::stof(num.at(5));
	  pgData.by = std::stof(num.at(6));
	  pgData.ebx = std::stof(num.at(7));
	  pgData.eby = std::stof(num.at(8));

	  continue;
	}
	
	// check the number of columns before assigning values
	bool assumePID = false;
	if(!isNewInfo && num.size() == 6) assumePID = false; 
	else if(isNewInfo && num.size() == 9) assumePID = false; 
	else if(isNewInfo2 && num.size() == 11) assumePID = false; 
	else if(!isNewInfo && (num.size() == 7 || num.size() == 8)) assumePID = true; 
	else if(isNewInfo && (num.size() == 10 || num.size() == 11)) assumePID = true; 
	else if(isNewInfo2 && (num.size() == 12 || num.size() == 13)) assumePID = true; 
	else{//return, this is an invalid format (or fix code here if format valid
	  std::cout << "Number of columns for line \'" << getStr << "\' is invalid, size \'" << num.size() << "\'. return 1" << std::endl;
	  //gotta cleanup before return
	  delete tout;
	  for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	    delete jout[jIter];
	    if(recombScheme[jIter]==fastjet::WTA_modp_scheme){
	      if(rParamIs8[jIter] || nFinalClust[jIter] == 2) delete bout[jIter];
	    }
	  }

	  if(isMC && isRecons){
	    delete tgout;
	    delete tgBeforeout;
	    for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	      delete jgout[jIter];
	      delete jgBeforeout[jIter];
              if(recombScheme[jIter]==fastjet::WTA_modp_scheme){
                if(rParamIs8[jIter] || nFinalClust[jIter] == 2){delete bgout[jIter]; delete bgBeforeout[jIter];}
              }
	    }
	  }

	  hf->Close();
	  delete hf;
	
	  return 1;
	}

	if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
      

	if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	float _px = std::stof(num.at(0));
	float _py = std::stof(num.at(1));
	float _pz = std::stof(num.at(2));
	float _m = std::stof(num.at(3));
	float _charge = std::stof(num.at(4));
	int _pwflag = std::stoi(num.at(5));

	if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	
	pgData.px[counterParticles]=_px;
	pgData.py[counterParticles]=_py;
	pgData.pz[counterParticles]=_pz;
	pgData.mass[counterParticles]=_m;
	v.SetXYZM(_px,_py,_pz,_m);
	fastjet::PseudoJet particle(_px,_py,_pz,v.E());
	particle.set_user_index(_pwflag);
	particles.push_back(particle);
	pgData.pt[counterParticles]=v.Pt();
	pgData.pmag[counterParticles]=v.Rho(); //Added later on
	pgData.rap[counterParticles]=v.Rapidity();
	pgData.eta[counterParticles]=v.PseudoRapidity();
	pgData.theta[counterParticles]=v.Theta();
	pgData.phi[counterParticles]=v.Phi();

	if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	
	pgData.charge[counterParticles]=_charge;
	pgData.pwflag[counterParticles]=_pwflag;
	//check before assigning PID
	if(assumePID) pgData.pid[counterParticles]=std::stoi(num.at(6));
	else pgData.pid[counterParticles]=-999;
      
	if((isNewInfo || isNewInfo2) && _pwflag < 4){
	  int posOffset = 0;
	  if(assumePID) posOffset += 1;
	  
	  pgData.vx[counterParticles] = std::stof(num.at(6+posOffset));
	  pgData.vy[counterParticles] = std::stof(num.at(7+posOffset));
	  pgData.vz[counterParticles] = std::stof(num.at(8+posOffset));

	  initVal({&(pgData.d0[counterParticles]), &(pgData.z0[counterParticles]), NULL}, -999);
	  initVal({&(pgData.ntpc[counterParticles]), &(pgData.nitc[counterParticles]), &(pgData.nvdet[counterParticles])}, -127);
	}
	else{
	  initVal({&(pgData.vx[counterParticles]), &(pgData.vy[counterParticles]), &(pgData.vz[counterParticles]), &(pgData.d0[counterParticles]), &(pgData.z0[counterParticles])}, -999);
	  initVal({&(pgData.ntpc[counterParticles]), &(pgData.nitc[counterParticles]), &(pgData.nvdet[counterParticles])}, -127);
	}
	
	if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

        //missing momentum calculation and multiplicity calculation
        netP -= TVector3(_px,_py,_pz);
        if(_pwflag==0){
          netP_charged -= TVector3(_px,_py,_pz);
          nTrk++;
          if(v.Pt()>0.4) nTrk_GT0p4++; 
          if(v.Pt()>0.4) nTrk_GT0p4Thrust++; 
        } 
      
	++counterParticles;	

	if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
      }
      //Have to fill one last time since the condition for fill is dependent on NEXT EVENT existing, else we would lose last event per file
      if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      std::vector<boostedEvtData*> bgDataTempVect;
      for(Int_t jI = 0; jI < nJtAlgo; ++jI){bgDataTempVect.push_back(&(bgData[jI]));}
      doEndEvent(&pgData, &egData, bgDataTempVect, counterParticles, nTrk, nTrk_GT0p4, nTrk_GT0p4Thrust, netP, netP_charged);
      bgDataTempVect.clear();
      
      if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      if(counterEntries>0){pgData.preFillClean(); tgout->Fill();}
      for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	processJets(particles, jDef[jIter], jDefReclust[jIter], &(jgData[jIter]), jtPtCut, rParam[jIter], nFinalClust[jIter]);
	if(counterEntries>0){
	  jgData[jIter].preFillClean();
	  jgout[jIter]->Fill();
	  if(recombScheme[jIter]==fastjet::WTA_modp_scheme){
	    if(rParamIs8[jIter] || nFinalClust[jIter] == 2){
	      if(jgData[jIter].nref<2) setBoostedVariables(false, &pgData, &(bgData[jIter]));
	      else{
		TVector3 wtaBoost = findBack2BackBoost(jgData[jIter].fourJet[0],jgData[jIter].fourJet[1]);
		setBoostedVariables(true, &pgData, &(bgData[jIter]), jgData[jIter].fourJet[0], wtaBoost);
	      }
	      bgData[jIter].preFillClean();
	      bgout[jIter]->Fill();
	    }
	  }
	}
      }
      if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      fileGen.close();
    }

    
    if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;


    if(isRecons && isMC){
      if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      std::string genFileStr = fileList.at(fI);
      genFileStr.replace(genFileStr.find("_recons_"), 8, "_mctrue_");
      genFileStr.replace(genFileStr.find("_aftercut-"), 10, "_beforecut-");

      if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      std::ifstream fileGen(genFileStr.c_str());
      initTVector3({&netP, &netP_charged});
      initVal({&nTrk, &nTrk_GT0p4, &nTrk_GT0p4Thrust, &counterParticles, &counterEntries}, 0);

      while(std::getline(fileGen,getStr)){
	if(getStr.size() == 0) continue;
	else if(getStr.find("******") != std::string::npos) continue;
	else if(getStr.find("END_EVENT") != std::string::npos) continue;
	else if(getStr.find("END_FILE") != std::string::npos) continue;
	std::vector<std::string> num = processAlephString(getStr);

	if(num.size() < 5){
	  std::cout << "Number of columns for line \'" << getStr << "\' is invalid, size \'" << num.size() << "\'. return 1" << std::endl;
	  //gotta cleanup before return
	  delete tout;
	  for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	    delete jout[jIter];
	    if(recombScheme[jIter]==fastjet::WTA_modp_scheme){
	      if(rParamIs8[jIter] || nFinalClust[jIter] == 2) delete bout[jIter];
	    }
	  }

	  if(isMC && isRecons){
	    delete tgout;
	    delete tgBeforeout;
	    for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	      delete jgout[jIter];
	      if(recombScheme[jIter]==fastjet::WTA_modp_scheme){
		if(rParamIs8[jIter] || nFinalClust[jIter] == 2){delete bgout[jIter]; delete bgBeforeout[jIter];}
	      }
	    }
	  }

	  hf->Close();
	  delete hf;
	
	  return 1;
	}

	if(check999(num.at(0)) && check999(num.at(1)) && check999(num.at(2)) /* && check999(num.at(3))*/){ 
	  std::vector<boostedEvtData*> bgBeforeDataTempVect;
	  for(Int_t jI = 0; jI < nJtAlgo; ++jI){bgBeforeDataTempVect.push_back(&(bgBeforeData[jI]));}
	  doEndEvent(&pgBeforeData, &egBeforeData, bgBeforeDataTempVect, counterParticles, nTrk, nTrk_GT0p4, nTrk_GT0p4Thrust, netP, netP_charged);
	  bgBeforeDataTempVect.clear();

	  if(counterEntries>0){pgBeforeData.preFillClean(); tgBeforeout->Fill();}

	  //Processing particles->jets
	  for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	    processJets(particles, jDef[jIter], jDefReclust[jIter], &(jgBeforeData[jIter]), jtPtCut, rParam[jIter], nFinalClust[jIter]);
	    if(counterEntries>0){
	      jgBeforeData[jIter].preFillClean();
	      jgBeforeout[jIter]->Fill();
	      if(recombScheme[jIter]==fastjet::WTA_modp_scheme){
		if(rParamIs8[jIter] || nFinalClust[jIter] == 2){
		  if(jgBeforeData[jIter].nref<2) setBoostedVariables(false, &pgBeforeData, &(bgBeforeData[jIter]));
		  else{
		    TVector3 wtaBoost = findBack2BackBoost(jgBeforeData[jIter].fourJet[0],jgBeforeData[jIter].fourJet[1]);
		    setBoostedVariables(true, &pgBeforeData, &(bgBeforeData[jIter]), jgBeforeData[jIter].fourJet[0], wtaBoost);
		  }
		  bgBeforeData[jIter].preFillClean();
		  bgBeforeout[jIter]->Fill();
		}
	      }		
	    }
	  }
	  //clear particles for next iteration clustering
	  particles.clear();

	  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	  unsigned int runPos = 3;
	  unsigned int evtPos = 4;
	  unsigned int ePos = 5;
	  
	  if(check999(num.at(3))){runPos++; evtPos++; ePos++;}
	  //	  if(check999(num.at(4))){runPos++; evtPos++; ePos++;}

	  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	  pgBeforeData.year = year;
	  pgBeforeData.subDir = subDir;
	  pgBeforeData.process = process;
	  pgBeforeData.isMC = isMC;
	  pgBeforeData.RunNo = std::stoi(num.at(runPos));
	  pgBeforeData.EventNo= std::stoi(num.at(evtPos));
	  pgBeforeData.Energy= std::stof(num.at(ePos));
	  
	  int tempSubDir = subDir;
	  if(tempSubDir < 0) tempSubDir = 0;
	  int tempProcess = process;
	  if(tempProcess < 0) tempProcess = 15;
	  
	  pgBeforeData.uniqueID = idMaker.getID(isMC, tempProcess, tempSubDir, year, pgBeforeData.RunNo, pgBeforeData.EventNo);

	  bool passesNTupleAfterCut = false;
	  for(unsigned int rI = 0; rI < runNo.size(); ++rI){
	    if(pgBeforeData.RunNo == runNo.at(rI) && pgBeforeData.EventNo == evtNo.at(rI)){
	      passesNTupleAfterCut = true;
	      break;
	    }
	  }

	  egBeforeData.passesNTupleAfterCut = passesNTupleAfterCut;

	  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	  initVal(&(pData.bFlag), -999.);
	  initVal({&(pData.bx), &(pData.by), &(pData.ebx), &(pData.eby)}, -999.);	 
	  initTVector3({&netP, &netP_charged});
	  initVal({&nTrk, &nTrk_GT0p4, &nTrk_GT0p4Thrust, &counterParticles}, 0);
	  
	  if(doLocalDebug){
	    std::cout << __FILE__ << ", " << __LINE__ << ", " << num.size() << ", " << counterEntries << std::endl;
	    std::cout << pgBeforeData.RunNo << ", " << pgBeforeData.EventNo << std::endl;
	    std::cout << runNo.size() << ", " << evtNo.size() << std::endl;
	  }

	  /*
	  if(pgBeforeData.RunNo != runNo.at(counterEntries) || pgBeforeData.EventNo != evtNo.at(counterEntries)){
	    std::cout << "Gen entries dont match reco for file \'" << genFileStr << "\'. return 1" << std::endl;
	    //gotta cleanup before return
	    delete tout;

	    for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	      delete jout[jIter];
	      
	      if(recombScheme[jIter]==fastjet::WTA_modp_scheme){
		if(rParamIs8[jIter] || nFinalClust[jIter] == 2) delete bout[jIter];
	      }
	    }

	    if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	    if(isMC && isRecons){
	      delete tgout;
	      delete tgBeforeout;
	      for(int jIter = 0; jIter < nJtAlgo; ++jIter){
		delete jgout[jIter];
		delete jgBeforeout[jIter];
		if(recombScheme[jIter]==fastjet::WTA_modp_scheme){
		  if(rParamIs8[jIter] || nFinalClust[jIter] == 2){delete bgout[jIter]; delete bgBeforeout[jIter];}
		}
	      }
	    }

	    if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	    hf->Close();
	    delete hf;
	    
	    return 1;
	  }
	  */
	  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	  
	  ++counterEntries;	

	  continue;
	}
	else if(check998(num.at(0)) && check998(num.at(1)) && check998(num.at(2))){
	  pgBeforeData.bFlag = std::stoi(num.at(4));
	  pgBeforeData.bx = std::stof(num.at(5));
	  pgBeforeData.by = std::stof(num.at(6));
	  pgBeforeData.ebx = std::stof(num.at(7));
	  pgBeforeData.eby = std::stof(num.at(8));

	  continue;
	}

	if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	
	// check the number of columns before assigning values
	bool assumePID = false;
	if(!isNewInfo && num.size() == 6) assumePID = false; 
	else if(isNewInfo && num.size() == 9) assumePID = false; 
	else if(isNewInfo2 && num.size() == 11) assumePID = false; 
	else if(!isNewInfo && (num.size() == 7 || num.size() == 8)) assumePID = true; 
	else if(isNewInfo && (num.size() == 10 || num.size() == 11)) assumePID = true; 
	else if(isNewInfo2 && (num.size() == 12 || num.size() == 13)) assumePID = true; 
	else{//return, this is an invalid format (or fix code here if format valid
	  std::cout << "Number of columns for line \'" << getStr << "\' is invalid, size \'" << num.size() << "\'. return 1" << std::endl;
	  //gotta cleanup before return
	  delete tout;
	  for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	    delete jout[jIter];
	    if(recombScheme[jIter]==fastjet::WTA_modp_scheme){
	      if(rParamIs8[jIter] || nFinalClust[jIter] == 2) delete bout[jIter];
	    }
	  }

	  if(isMC && isRecons){
	    delete tgout;
	    delete tgBeforeout;
	    for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	      delete jgout[jIter];
	      delete jgBeforeout[jIter];
              if(recombScheme[jIter]==fastjet::WTA_modp_scheme){
                if(rParamIs8[jIter] || nFinalClust[jIter] == 2){delete bgout[jIter]; delete bgBeforeout[jIter];}
              }
	    }
	  }

	  hf->Close();
	  delete hf;
	
	  return 1;
	}

	if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
      

	if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	float _px = std::stof(num.at(0));
	float _py = std::stof(num.at(1));
	float _pz = std::stof(num.at(2));
	float _m = std::stof(num.at(3));
	float _charge = std::stof(num.at(4));
	int _pwflag = std::stoi(num.at(5));

	if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	
	pgBeforeData.px[counterParticles]=_px;
	pgBeforeData.py[counterParticles]=_py;
	pgBeforeData.pz[counterParticles]=_pz;
	pgBeforeData.mass[counterParticles]=_m;
	v.SetXYZM(_px,_py,_pz,_m);
	fastjet::PseudoJet particle(_px,_py,_pz,v.E());
	particle.set_user_index(_pwflag);
	particles.push_back(particle);
	pgBeforeData.pt[counterParticles]=v.Pt();
	pgBeforeData.pmag[counterParticles]=v.Rho(); //Added later on
	pgBeforeData.rap[counterParticles]=v.Rapidity();
	pgBeforeData.eta[counterParticles]=v.PseudoRapidity();
	pgBeforeData.theta[counterParticles]=v.Theta();
	pgBeforeData.phi[counterParticles]=v.Phi();

	if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	
	pgBeforeData.charge[counterParticles]=_charge;
	pgBeforeData.pwflag[counterParticles]=_pwflag;
	//check before assigning PID
	if(assumePID) pgBeforeData.pid[counterParticles]=std::stoi(num.at(6));
	else pgBeforeData.pid[counterParticles]=-999;
      
	if((isNewInfo || isNewInfo2) && _pwflag < 4){
	  int posOffset = 0;
	  if(assumePID) posOffset += 1;
	  
	  pgBeforeData.vx[counterParticles] = std::stof(num.at(6+posOffset));
	  pgBeforeData.vy[counterParticles] = std::stof(num.at(7+posOffset));
	  pgBeforeData.vz[counterParticles] = std::stof(num.at(8+posOffset));

	  initVal({&(pgBeforeData.d0[counterParticles]), &(pgBeforeData.z0[counterParticles]), NULL}, -999);
	  initVal({&(pgBeforeData.ntpc[counterParticles]), &(pgBeforeData.nitc[counterParticles]), &(pgBeforeData.nvdet[counterParticles])}, -127);
	}
	else{
	  initVal({&(pgBeforeData.vx[counterParticles]), &(pgBeforeData.vy[counterParticles]), &(pgBeforeData.vz[counterParticles]), &(pgBeforeData.d0[counterParticles]), &(pgBeforeData.z0[counterParticles])}, -999);
	  initVal({&(pgBeforeData.ntpc[counterParticles]), &(pgBeforeData.nitc[counterParticles]), &(pgBeforeData.nvdet[counterParticles])}, -127);
	}
	
	if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

        //missing momentum calculation and multiplicity calculation
        netP -= TVector3(_px,_py,_pz);
        if(_pwflag==0){
          netP_charged -= TVector3(_px,_py,_pz);
          nTrk++;
          if(v.Pt()>0.4) nTrk_GT0p4++; 
          if(v.Pt()>0.4) nTrk_GT0p4Thrust++; 
        } 
      
	++counterParticles;	

	if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
      }
      //Have to fill one last time since the condition for fill is dependent on NEXT EVENT existing, else we would lose last event per file
      std::vector<boostedEvtData*> bgBeforeDataTempVect;
      for(Int_t jI = 0; jI < nJtAlgo; ++jI){bgBeforeDataTempVect.push_back(&(bgBeforeData[jI]));}
      doEndEvent(&pgBeforeData, &egBeforeData, bgBeforeDataTempVect, counterParticles, nTrk, nTrk_GT0p4, nTrk_GT0p4Thrust, netP, netP_charged);
      bgBeforeDataTempVect.clear();
      
      if(counterEntries>0){pgBeforeData.preFillClean(); tgBeforeout->Fill();}
      for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	processJets(particles, jDef[jIter], jDefReclust[jIter], &(jgBeforeData[jIter]), jtPtCut, rParam[jIter], nFinalClust[jIter]);
	if(counterEntries>0){
	  jgBeforeData[jIter].preFillClean();
	  jgBeforeout[jIter]->Fill();
	  if(recombScheme[jIter]==fastjet::WTA_modp_scheme){
	    if(rParamIs8[jIter] || nFinalClust[jIter] == 2){
	      if(jgBeforeData[jIter].nref<2) setBoostedVariables(false, &pgBeforeData, &(bgBeforeData[jIter]));
	      else{
		TVector3 wtaBoost = findBack2BackBoost(jgBeforeData[jIter].fourJet[0],jgBeforeData[jIter].fourJet[1]);
		setBoostedVariables(true, &pgBeforeData, &(bgBeforeData[jIter]), jgBeforeData[jIter].fourJet[0], wtaBoost);
	      }
	      bgBeforeData[jIter].preFillClean();
	      bgBeforeout[jIter]->Fill();
	    }
	  }
	}
      }


      fileGen.close();
    }

    runNo.clear();
    evtNo.clear();
  }

  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  hf->cd();
  
  tout->Write("", TObject::kOverwrite);
  delete tout;

  for(int i = 0; i < nJtAlgo; ++i){
    jout[i]->Write("", TObject::kOverwrite);
    delete jout[i];

    if(recombScheme[i]==fastjet::WTA_modp_scheme){
      if(rParamIs8[i] || nFinalClust[i] == 2){
	bout[i]->Write("", TObject::kOverwrite);
	delete bout[i];
      }
    }
  }

  if(isMC && isRecons){
    tgout->Write("", TObject::kOverwrite);
    delete tgout;

    tgBeforeout->Write("", TObject::kOverwrite);
    delete tgBeforeout;

    for(int jIter = 0; jIter < nJtAlgo; ++jIter){
      jgout[jIter]->Write("", TObject::kOverwrite);
      delete jgout[jIter];
      jgBeforeout[jIter]->Write("", TObject::kOverwrite);
      delete jgBeforeout[jIter];
      if(recombScheme[jIter]==fastjet::WTA_modp_scheme){
	if(rParamIs8[jIter] || nFinalClust[jIter] == 2){
	  bgout[jIter]->Write("", TObject::kOverwrite);
	  delete bgout[jIter];
	  bgBeforeout[jIter]->Write("", TObject::kOverwrite);
	  delete bgBeforeout[jIter];
	}
      }    
    }
  }
  
  hf->Close();
  delete hf;

  std::cout << "Job complete, return 0." << std::endl;

  return 0;
}

int main(int argc, char *argv[])
{
  if(argc != 4 && argc != 5){
    std::cout << "Usage: ./bin/scan.exe <inFileName> <isNewInfo> <isNewInfo2> <OPT-outFileName>" << std::endl;
    return 1;
  }

  std::vector<std::string> fileListTemp;
  std::string inFileName = argv[1];

  int nMix = -1;

  if(inFileName.substr(inFileName.size()-4, 4).find(".txt") != std::string::npos){
    std::ifstream pathFile(inFileName.c_str());
    std::string paths;
    while(std::getline(pathFile, paths)){
      if(paths.size() == 0) continue;
      fileListTemp.push_back(paths);
    }
    pathFile.close();
  }
  else if(inFileName.substr(inFileName.size()-5, 5).find(".list") != std::string::npos) fileListTemp.push_back(inFileName);
  else if(inFileName.substr(inFileName.size()-6, 6).find(".aleph") != std::string::npos) fileListTemp.push_back(inFileName);
  else{
    std::cout << "WARNING: inFileName \'" << inFileName << "\' cannot determine MC or Data, mixing 1 event only." << std::endl;
    nMix = 1;
  }

  if(nMix < 0){
    if(getIsMC(fileListTemp.at(0))) nMix = 1;
    else nMix = 3;
  }

  std::cout << "Begin processing..." << std::endl;
  int retVal = 0;
  if(argc == 4) retVal += scan(argv[1], argv[2], argv[3]);
  else if(argc == 5) retVal += scan(argv[1], argv[2], argv[3], argv[4]);
  
  std::cout << "Making mixing file..." << std::endl;
  if(argc == 5) retVal += makeMixFile(argv[4], argv[4], nMix);
  else          retVal += makeMixFile(argv[1], "", nMix);
  return retVal;
}
