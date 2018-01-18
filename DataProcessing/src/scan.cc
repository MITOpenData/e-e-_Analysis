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
#include "include/eventData.h"
#include "include/jetData.h"
#include "include/boostedEvtData.h"
#include "include/thrustTools.h"
#include "include/boostTools.h"

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

void processJets(std::vector<fastjet::PseudoJet> p, fastjet::JetDefinition jDef, fastjet::JetDefinition jDefReclust, jetData *d, const double ptCut = 0.1)
{
  d->nref = 0;
  if(p.size() > 0){
    fastjet::ClusterSequence cs(p, jDef);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
    for(unsigned int i = 0; i < jets.size(); ++i){
      if(jets.at(i).pt() < ptCut) break; //Arbitrarily low cut on jets, removes spike at phi zero when things become ill defined
      d->jtpt[d->nref] = jets.at(i).pt();
      d->jtphi[d->nref] = jets.at(i).phi_std();
      d->jteta[d->nref] = jets.at(i).eta();
      
      std::vector<fastjet::PseudoJet> jetConst = jets.at(i).constituents();
      d->jtN[d->nref] = jetConst.size();
      std::vector< std::vector<fastjet::PseudoJet> > subJets;

      for(int j = 0; j < 6; ++j){
	d->jtNPW[d->nref][j] = 0;
	d->jtptFracPW[d->nref][j] = 0;

	std::vector<fastjet::PseudoJet> tempSubJets;
	subJets.push_back(tempSubJets);
      }

      for(unsigned int k = 0; k < jetConst.size(); ++k){
	if(jetConst.at(k).user_index() < 0 || jetConst.at(k).user_index() > 5) continue;
	subJets.at(jetConst.at(k).user_index()).push_back(jetConst.at(k));
	d->jtNPW[d->nref][jetConst.at(k).user_index()]++;
      }
      
      for(int j = 0; j < 6; ++j){
	fastjet::ClusterSequence csSub(subJets.at(j), jDefReclust);

	std::vector<fastjet::PseudoJet> constTot = fastjet::sorted_by_pt(csSub.inclusive_jets());
	if(constTot.size() > 1){
	  std::cout << "WARNING - RECLUSTER OF CONSTITUENTS YIELDS GREATER THAN 1 JET" << std::endl;
	  std::cout << "Top jet: " << jets.at(i).pt() << ", " << jets.at(i).phi() << ", " << jets.at(i).eta() << std::endl;
	  std::cout << "Const: " << std::endl;
	  for(unsigned int k = 0; k < subJets.at(j).size(); ++k){
	    std::cout << " " << k << "/" << subJets.at(j).size() << ": " << subJets.at(j).at(k).pt() << ", " << subJets.at(j).at(k).phi() << ", " << subJets.at(j).at(k).eta() << std::endl;
	  }

	  std::cout << "Reclust: " << std::endl;
	  for(unsigned int k = 0; k < constTot.size(); ++k){
	    std::cout << " " << k << "/" << constTot.size() << ": " << constTot.at(k).pt() << ", " << constTot.at(k).phi() << ", " << constTot.at(k).eta() << std::endl;
	  }

	}
	else if(constTot.size() == 1) d->jtptFracPW[d->nref][j] = constTot.at(0).pt()/jets.at(i).pt();
	else d->jtptFracPW[d->nref][j] = 0.;
      }

      ++d->nref;
    }
    jets.clear();
    std::vector<fastjet::PseudoJet>().swap(jets);
  }
  return;
}

int scan(std::string inFileName, const bool isNewInfo, std::string outFileName="")
{
  if(!checkFile(inFileName)){
    std::cout << "Given inFileName \'" << inFileName << "\' is invalid, return 1." << std::endl;
    return 1;
  }

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
  const int nJtAlgo = 4;
  const double rParam[nJtAlgo] = {0.4, 0.4, 0.8, 0.8};
  const double recombScheme[nJtAlgo] = {fastjet::E_scheme, fastjet::WTA_modp_scheme, fastjet::E_scheme, fastjet::WTA_modp_scheme};
  fastjet::JetDefinition jDef[nJtAlgo];
  fastjet::JetDefinition jDefReclust[nJtAlgo];
  
  for(int i = 0; i < nJtAlgo; ++i){
    jDef[i] = fastjet::JetDefinition(fastjet::ee_genkt_algorithm, rParam[i], -1, fastjet::RecombinationScheme(recombScheme[i]));
    jDefReclust[i] = fastjet::JetDefinition(fastjet::ee_genkt_algorithm, 5, -1, fastjet::RecombinationScheme(recombScheme[i]));
  }

  const std::string partTreeName = "t";
  const std::string genPartTreeName = "tgen";
  std::string jetTreeName[nJtAlgo];
  std::string genJetTreeName[nJtAlgo];
  const std::string boostedTreeName = "BoostedWTAEvt";
  const std::string genboostedTreeName = "genBoostedWTAEvt";
  for(int i = 0; i < nJtAlgo; ++i){
    std::string recombSchemeStr = "EScheme";
    if(recombScheme[i] == fastjet::WTA_modp_scheme) recombSchemeStr = "WTAmodpScheme";

    jetTreeName[i] = "ak" + std::to_string(int(rParam[i]*10)) + recombSchemeStr + "JetTree";
    genJetTreeName[i] = "ak" + std::to_string(int(rParam[i]*10)) + recombSchemeStr + "GenJetTree";
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
  TTree *bout = new TTree(boostedTreeName.c_str(), boostedTreeName.c_str());
  TTree *jout[nJtAlgo];
  for(int i = 0; i < nJtAlgo; ++i){jout[i] = new TTree(finalJetTreeName[i].c_str(), finalJetTreeName[i].c_str());}


  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  TTree *tgout = 0;
  TTree *bgout = 0;
  TTree *jgout[nJtAlgo] = {0, 0, 0, 0};

  if(isRecons && isMC){
    tgout = new TTree(genPartTreeName.c_str(), genPartTreeName.c_str());
    bgout = new TTree(genboostedTreeName.c_str(), genboostedTreeName.c_str()); 

    for(int i = 0; i < nJtAlgo; ++i){jgout[i] = new TTree(genJetTreeName[i].c_str(), genJetTreeName[i].c_str());}
  }

  particleData pData;
  jetData jData[nJtAlgo];
  eventData eData;
  boostedEvtData bData;

  particleData pgData;
  jetData jgData[nJtAlgo];
  eventData egData;
  boostedEvtData bgData;

  tout->Branch("year", &pData.year, "year/I");
  tout->Branch("EventNo", &pData.EventNo,"EventNo/I");
  tout->Branch("RunNo", &pData.RunNo,"RunNo/I");
  tout->Branch("Energy", &pData.Energy,"Energy/F");
  tout->Branch("process", &pData.process, "process/I");
  tout->Branch("bFlag", &pData.bFlag, "bFlag/I");
  tout->Branch("bx", &pData.bx, "bx/F");
  tout->Branch("by", &pData.by, "by/F");
  tout->Branch("ebx", &pData.ebx, "ebx/F");
  tout->Branch("eby", &pData.eby, "eby/F");
  tout->Branch("nParticle", &pData.nParticle,"nParticle/I");
  tout->Branch("px", pData.px,"px[nParticle]/F");
  tout->Branch("py", pData.py,"py[nParticle]/F");
  tout->Branch("pz", pData.pz,"pz[nParticle]/F");
  tout->Branch("pt", pData.pt,"pt[nParticle]/F");
  tout->Branch("pmag", pData.pmag,"pmag[nParticle]/F");//Added later on
  tout->Branch("mass", pData.mass,"mass[nParticle]/F");
  tout->Branch("eta", pData.eta,"eta[nParticle]/F");
  tout->Branch("theta", pData.theta,"theta[nParticle]/F");
  tout->Branch("phi", pData.phi,"phi[nParticle]/F");
  tout->Branch("charge", pData.charge,"charge[nParticle]/F");
  tout->Branch("pwflag", pData.pwflag,"pwflag[nParticle]/I");
  tout->Branch("pid", pData.pid,"pid[nParticle]/I");
  tout->Branch("d0", pData.d0,"d0[nParticle]/F");
  tout->Branch("z0", pData.z0,"z0[nParticle]/F");
  tout->Branch("ntpc", pData.ntpc,"ntpc[nParticle]/I");
 
  //thrust quantities
  tout->Branch("Thrust",&eData.Thrust,"Thrust/F");
  tout->Branch("TTheta",&eData.TTheta,"TTheta/F");
  tout->Branch("TPhi",&eData.TPhi,"TPhi/F");
  tout->Branch("Thrust_charged",&eData.Thrust_charged,"Thrust_charged/F");
  tout->Branch("TTheta_charged",&eData.TTheta_charged,"TTheta_charged/F");
  tout->Branch("TPhi_charged",&eData.TPhi_charged,"TPhi_charged/F");
  tout->Branch("pt_wrtThr", pData.pt_wrtThr,"pt_wrtThr[nParticle]/F");
  tout->Branch("eta_wrtThr", pData.eta_wrtThr,"eta_wrtThr[nParticle]/F");
  tout->Branch("theta_wrtThr", pData.theta_wrtThr,"theta_wrtThr[nParticle]/F");
  tout->Branch("phi_wrtThr", pData.phi_wrtThr,"phi_wrtThr[nParticle]/F");
  tout->Branch("pt_wrtChThr", pData.pt_wrtChThr,"pt_wrtChThr[nParticle]/F");
  tout->Branch("eta_wrtChThr", pData.eta_wrtChThr,"eta_wrtChThr[nParticle]/F");
  tout->Branch("theta_wrtChThr", pData.theta_wrtChThr,"theta_wrtChThr[nParticle]/F");
  tout->Branch("phi_wrtChThr", pData.phi_wrtChThr,"phi_wrtChThr[nParticle]/F");

  //derived quantities
  tout->Branch("missP",&eData.missP,"missP/F");
  tout->Branch("missPt",&eData.missPt,"missPt/F");
  tout->Branch("missTheta",&eData.missTheta,"missTheta/F");
  tout->Branch("missPhi",&eData.missPhi,"missPhi/F");
  tout->Branch("missChargedP",&eData.missChargedP,"missP/F");
  tout->Branch("missChargedPt",&eData.missChargedPt,"missChargedPt/F");
  tout->Branch("missChargedTheta",&eData.missChargedTheta,"missChargedTheta/F");
  tout->Branch("missChargedPhi",&eData.missChargedPhi,"missChargedPhi/F");
  tout->Branch("nChargedHadrons",&eData.nChargedHadrons,"nChargedHadrons/I");
  tout->Branch("nChargedHadrons_GT0p4",&eData.nChargedHadrons_GT0p4,"nChargedHadrons_GT0p4/I");
  tout->Branch("nChargedHadrons_GT0p4Thrust",&eData.nChargedHadrons_GT0p4Thrust,"nChargedHadrons_GT0p4Thrust/I");

  bout->Branch("nParticle",&bData.nParticle,"nParticle/I");
  bout->Branch("WTAAxis_Theta",&bData.WTAAxis_Theta,"WTAAxis_Theta/F");
  bout->Branch("WTAAxis_Phi",&bData.WTAAxis_Phi,"WTAAxis_Phi/F");
  bout->Branch("pt", bData.pt,"pt[nParticle]/F");
  bout->Branch("pmag", bData.pmag,"pmag[nParticle]/F");//Added later on
  bout->Branch("eta", bData.eta,"eta[nParticle]/F");
  bout->Branch("theta", bData.theta,"theta[nParticle]/F");
  bout->Branch("phi", bData.phi,"phi[nParticle]/F");

  for(int i = 0; i < nJtAlgo; ++i){
    jout[i]->Branch("nref", &jData[i].nref,"nref/I");
    jout[i]->Branch("jtpt", jData[i].jtpt,"jtpt[nref]/F");
    jout[i]->Branch("jteta", jData[i].jteta,"jteta[nref]/F");
    jout[i]->Branch("jtphi", jData[i].jtphi,"jtphi[nref]/F");
    jout[i]->Branch("jtN", jData[i].jtN, "jtN[nref]/I");
    jout[i]->Branch("jtNPW", jData[i].jtNPW, "jtNPW[nref][6]/I");
    jout[i]->Branch("jtptFracPW", jData[i].jtptFracPW, "jtptFracPW[nref][6]/F");
  }

  if(isRecons && isMC){
    tgout->Branch("year", &pgData.year, "year/I");
    tgout->Branch("EventNo", &pgData.EventNo,"EventNo/I");
    tgout->Branch("RunNo", &pgData.RunNo,"RunNo/I");
    tgout->Branch("Energy", &pgData.Energy,"Energy/F");
    tgout->Branch("process", &pgData.process, "process/I");
    tgout->Branch("bFlag", &pgData.bFlag, "bFlag/I");
    tgout->Branch("bx", &pgData.bx, "bx/F");
    tgout->Branch("by", &pgData.by, "by/F");
    tgout->Branch("ebx", &pgData.ebx, "ebx/F");
    tgout->Branch("eby", &pgData.eby, "eby/F");
    tgout->Branch("nParticle", &pgData.nParticle,"nParticle/I");
    tgout->Branch("px", pgData.px,"px[nParticle]/F");
    tgout->Branch("py", pgData.py,"py[nParticle]/F");
    tgout->Branch("pz", pgData.pz,"pz[nParticle]/F");
    tgout->Branch("pt", pgData.pt,"pt[nParticle]/F");
    tgout->Branch("pmag", pgData.pmag,"pmag[nParticle]/F");//Added later on
    tgout->Branch("mass", pgData.mass,"mass[nParticle]/F");
    tgout->Branch("eta", pgData.eta,"eta[nParticle]/F");
    tgout->Branch("theta", pgData.theta,"theta[nParticle]/F");
    tgout->Branch("phi", pgData.phi,"phi[nParticle]/F");
    tgout->Branch("charge", pgData.charge,"charge[nParticle]/F");
    tgout->Branch("pwflag", pgData.pwflag,"pwflag[nParticle]/I");
    tgout->Branch("pid", pgData.pid,"pid[nParticle]/I");
    tgout->Branch("d0", pgData.d0,"d0[nParticle]/F");
    tgout->Branch("z0", pgData.z0,"z0[nParticle]/F");
    tgout->Branch("ntpc", pgData.ntpc,"ntpc[nParticle]/I");

    //thrust quantities
    tgout->Branch("Thrust",&egData.Thrust,"Thrust/F");
    tgout->Branch("TTheta",&egData.TTheta,"TTheta/F");
    tgout->Branch("TPhi",&egData.TPhi,"TPhi/F");
    tgout->Branch("Thrust_charged",&egData.Thrust_charged,"Thrust_charged/F");
    tgout->Branch("TTheta_charged",&egData.TTheta_charged,"TTheta_charged/F");
    tgout->Branch("TPhi_charged",&egData.TPhi_charged,"TPhi_charged/F");
    tgout->Branch("pt_wrtThr", pgData.pt_wrtThr,"pt_wrtThr[nParticle]/F");
    tgout->Branch("eta_wrtThr", pgData.eta_wrtThr,"eta_wrtThr[nParticle]/F");
    tgout->Branch("theta_wrtThr", pgData.theta_wrtThr,"theta_wrtThr[nParticle]/F");
    tgout->Branch("phi_wrtThr", pgData.phi_wrtThr,"phi_wrtThr[nParticle]/F");
    tgout->Branch("pt_wrtChThr", pgData.pt_wrtChThr,"pt_wrtChThr[nParticle]/F");
    tgout->Branch("eta_wrtChThr", pgData.eta_wrtChThr,"eta_wrtChThr[nParticle]/F");
    tgout->Branch("theta_wrtChThr", pgData.theta_wrtChThr,"theta_wrtChThr[nParticle]/F");
    tgout->Branch("phi_wrtChThr", pgData.phi_wrtChThr,"phi_wrtChThr[nParticle]/F");

    //derived quantities
    tgout->Branch("missP",&egData.missP,"missP/F");
    tgout->Branch("missPt",&egData.missPt,"missPt/F");
    tgout->Branch("missTheta",&egData.missTheta,"missTheta/F");
    tgout->Branch("missPhi",&egData.missPhi,"missPhi/F");
    tgout->Branch("missChargedP",&egData.missChargedP,"missP/F");
    tgout->Branch("missChargedPt",&egData.missChargedPt,"missChargedPt/F");
    tgout->Branch("missChargedTheta",&egData.missChargedTheta,"missChargedTheta/F");
    tgout->Branch("missChargedPhi",&egData.missChargedPhi,"missChargedPhi/F");
    tgout->Branch("nChargedHadrons",&egData.nChargedHadrons,"nChargedHadrons/I");
    tgout->Branch("nChargedHadrons_GT0p4",&egData.nChargedHadrons_GT0p4,"nChargedHadrons_GT0p4/I");
    tgout->Branch("nChargedHadrons_GT0p4Thrust",&egData.nChargedHadrons_GT0p4Thrust,"nChargedHadrons_GT0p4Thrust/I");

    bgout->Branch("nParticle",&bgData.nParticle,"nParticle/I");
    bgout->Branch("WTAAxis_Theta",&bgData.WTAAxis_Theta,"WTAAxis_Theta/F");
    bgout->Branch("WTAAxis_Phi",&bgData.WTAAxis_Phi,"WTAAxis_Phi/F");
    bgout->Branch("pt", bgData.pt,"pt[nParticle]/F");
    bgout->Branch("pmag", bgData.pmag,"pmag[nParticle]/F");//Added later on
    bgout->Branch("eta", bgData.eta,"eta[nParticle]/F");
    bgout->Branch("theta", bgData.theta,"theta[nParticle]/F");
    bgout->Branch("phi", bgData.phi,"phi[nParticle]/F");
    
    for(int i = 0; i < nJtAlgo; ++i){
      jgout[i]->Branch("nref", &jgData[i].nref,"nref/I");
      jgout[i]->Branch("jtpt", jgData[i].jtpt,"jtpt[nref]/F");
      jgout[i]->Branch("jteta", jgData[i].jteta,"jteta[nref]/F");
      jgout[i]->Branch("jtphi", jgData[i].jtphi,"jtphi[nref]/F");
      jgout[i]->Branch("jtN", jgData[i].jtN,"jtN[nref]/I");
      jgout[i]->Branch("jtNPW", jgData[i].jtNPW,"jtNPW[nref][6]/I");
      jgout[i]->Branch("jtptFracPW", jgData[i].jtptFracPW,"jtptFracPW[nref][6]/F");
    }
  }

  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    std::cout << "Processing " << fI << "/" << fileList.size() << ", \'" << fileList.at(fI) << "\'" << std::endl;
    const int year = yearFromPath(fileList.at(fI));
    const int process = processFromPath(fileList.at(fI));

    std::ifstream file(Form("%s",fileList.at(fI).c_str()));
    std::string getStr;
    
    int counterEntries=0;
    int counterParticles=0;
    TLorentzVector v;
    std::vector<fastjet::PseudoJet> particles;

    std::vector<int> runNo;
    std::vector<int> evtNo;
    TVector3 netP = TVector3(0,0,0);
    TVector3 thrust = TVector3(0,0,0);
    TVector3 netP_charged = TVector3(0,0,0);
    TVector3 thrust_charged = TVector3(0,0,0);
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
	delete bout;
	for(int jIter = 0; jIter < nJtAlgo; ++jIter){delete jout[jIter];}
	
	if(isMC && isRecons){
	  delete tgout;
	  delete bgout;
	  for(int jIter = 0; jIter < nJtAlgo; ++jIter){delete jgout[jIter];}
	}
	
	hf->Close();
	delete hf;
	
	return 1;
      }

      if(check999(num.at(0)) && check999(num.at(1)) && check999(num.at(2))/* not all files do four 999 && check999(num.at(3)*/){
	pData.nParticle=counterParticles;
        bData.nParticle=counterParticles;
        thrust = getThrust(pData.nParticle, pData.px, pData.py, pData.pz, THRUST::OPTIMAL); 
        thrust_charged = getChargedThrust(pData.nParticle, pData.px, pData.py, pData.pz, pData.pwflag, THRUST::OPTIMAL);
        setThrustVariables(&pData, &eData, thrust, thrust_charged);
        eData.Thrust = thrust.Mag();
        eData.TTheta = thrust.Theta();
        eData.TPhi = thrust.Phi();
        eData.Thrust_charged = thrust_charged.Mag();
        eData.TTheta_charged = thrust_charged.Theta();
        eData.TPhi_charged = thrust_charged.Phi();
      

	if(counterEntries>0) tout->Fill(); 
	
	//Processing particles->jets
	for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	  processJets(particles, jDef[jIter], jDefReclust[jIter], &(jData[jIter]), jtPtCut);
	}

	if(counterEntries>0){
	  for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	    jout[jIter]->Fill();
	  }
          bout->Fill();
	}        

	//clear particles for next iteration clustering
	particles.clear();

	unsigned int runPos = 3;
	unsigned int evtPos = 4;
	unsigned int ePos = 5;

	//again, 4 999 is max	if(check999(num.at(4))){runPos++; evtPos++; ePos++;}
	if(check999(num.at(3))){runPos++; evtPos++; ePos++;}

	pData.year = year;
	pData.process = process;
	pData.RunNo = std::stoi(num.at(runPos));
	pData.EventNo= std::stoi(num.at(evtPos));
	pData.Energy= std::stof(num.at(ePos));

	pData.bFlag = -999;
	pData.bx = -999.;
	pData.by = -999.;
	pData.ebx = -999.;
	pData.eby = -999.;
      
	runNo.push_back(pData.RunNo);
	evtNo.push_back(pData.EventNo);

        netP = TVector3(0,0,0);
        thrust = TVector3(0,0,0);
        netP_charged = TVector3(0,0,0);
        thrust_charged = TVector3(0,0,0);
        nTrk=0;
        nTrk_GT0p4=0;
        nTrk_GT0p4Thrust=0;

	counterParticles=0;   
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
      else if(!isNewInfo && (num.size() == 7 || num.size() == 8)) assumePID = true; 
      else if(isNewInfo && (num.size() == 10 || num.size() == 11)) assumePID = true; 
      else{//return, this is an invalid format (or fix code here if format valid
	std::cout << "Number of columns for line \'" << getStr << "\' is invalid, size \'" << num.size() << "\'. return 1" << std::endl;
	//gotta cleanup before return
	delete tout;
	delete bout;
	for(int jIter = 0; jIter < nJtAlgo; ++jIter){delete jout[jIter];}

	if(isMC && isRecons){
	  delete tgout;
	  delete bgout;
	  for(int jIter = 0; jIter < nJtAlgo; ++jIter){delete jgout[jIter];}
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
      pData.eta[counterParticles]=v.PseudoRapidity();
      pData.theta[counterParticles]=v.Theta();
      pData.phi[counterParticles]=v.Phi();
      
      pData.charge[counterParticles]=_charge;
      pData.pwflag[counterParticles]=_pwflag;
      //check before assigning PID
      if(assumePID) pData.pid[counterParticles]=std::stoi(num.at(6));
      else pData.pid[counterParticles]=-999;
     
      if(isNewInfo && _pwflag < 4){
	int posOffset = 0;
	if(assumePID) posOffset = 1;

	pData.d0[counterParticles] = std::stof(num.at(6 + posOffset));
	pData.z0[counterParticles] = std::stof(num.at(7 + posOffset));
	pData.ntpc[counterParticles] =  std::stoi(num.at(8 + posOffset));	
      }
      else{
	pData.d0[counterParticles] = -999.;
	pData.z0[counterParticles] = -999.;
	pData.ntpc[counterParticles] = -999;
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
      eData.missP = netP.Mag();
      eData.missPt = netP.Perp();
      eData.missTheta = netP.Theta();
      eData.missPhi = netP.Phi();
      eData.missChargedP = netP_charged.Mag();
      eData.missChargedPt = netP_charged.Perp();
      eData.missChargedTheta = netP_charged.Theta();
      eData.missChargedPhi = netP_charged.Phi();
      eData.nChargedHadrons = nTrk;
      eData.nChargedHadrons_GT0p4 = nTrk_GT0p4; 
      eData.nChargedHadrons_GT0p4Thrust = nTrk_GT0p4Thrust; 
 
      ++counterParticles;	
    }

    if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  
    //Have to fill one last time since the condition for fill is dependent on NEXT EVENT existing, else we would lose last event per file
    pData.nParticle=counterParticles;
    thrust = getThrust(pData.nParticle, pData.px, pData.py, pData.pz, THRUST::OPTIMAL); 
    thrust_charged = getChargedThrust(pData.nParticle, pData.px, pData.py, pData.pz, pData.pwflag, THRUST::OPTIMAL);
    setThrustVariables(&pData, &eData, thrust, thrust_charged);

    eData.Thrust = thrust.Mag();
    eData.TTheta = thrust.Theta();
    eData.TPhi = thrust.Phi();
    eData.Thrust_charged = thrust_charged.Mag();
    eData.TTheta_charged = thrust_charged.Theta();
    eData.TPhi_charged = thrust_charged.Phi();

    if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    
    if(counterEntries>0) tout->Fill(); 
    for(int jIter = 0; jIter < nJtAlgo; ++jIter){
      processJets(particles, jDef[jIter], jDefReclust[jIter], &(jData[jIter]), jtPtCut);
    }
    if(counterEntries>0){
      for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	jout[jIter]->Fill();
      }
      bout->Fill();
    }
    file.close();
  
    if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    if(isRecons && isMC){
      std::string genFileStr = fileList.at(fI);
      genFileStr.replace(genFileStr.find("_recons_"), 8, "_mctrue_");

      std::ifstream fileGen(genFileStr.c_str());
      counterEntries=0;
      counterParticles=0;
      netP = TVector3(0,0,0);
      thrust = TVector3(0,0,0);
      netP_charged = TVector3(0,0,0);
      thrust_charged = TVector3(0,0,0);
      nTrk=0;
      nTrk_GT0p4=0;
      nTrk_GT0p4Thrust=0;

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
          delete bout;
	  for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	    delete jout[jIter];
	  }

	  if(isMC && isRecons){
	    delete tgout;
	    delete bgout;
	    for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	      delete jgout[jIter];
	    }
	  }

	  hf->Close();
	  delete hf;
	
	  return 1;
	}

	if(check999(num.at(0)) && check999(num.at(1)) && check999(num.at(2)) /* && check999(num.at(3))*/){ 
	  pgData.nParticle=counterParticles;
	  bgData.nParticle=counterParticles;
          
          thrust = getThrust(pgData.nParticle, pgData.px, pgData.py, pgData.pz, THRUST::OPTIMAL); 
          thrust_charged = getChargedThrust(pgData.nParticle, pgData.px, pgData.py, pgData.pz, pgData.pwflag, THRUST::OPTIMAL);
          setThrustVariables(&pgData, &egData, thrust, thrust_charged);
          egData.Thrust = thrust.Mag();
          egData.TTheta = thrust.Theta();
          egData.TPhi = thrust.Phi();
          egData.Thrust_charged = thrust_charged.Mag();
          egData.TTheta_charged = thrust_charged.Theta();
          egData.TPhi_charged = thrust_charged.Phi();

	  if(counterEntries>0) tgout->Fill(); 

	  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;	

	  //Processing particles->jets
	  for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	    processJets(particles, jDef[jIter], jDefReclust[jIter], &(jgData[jIter]), jtPtCut);
	    if(counterEntries>0){
	      jgout[jIter]->Fill();
	    }
            bgout->Fill();
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
	  pgData.process = process;
	  pgData.RunNo = std::stoi(num.at(runPos));
	  pgData.EventNo= std::stoi(num.at(evtPos));
	  pgData.Energy= std::stof(num.at(ePos));

	  pgData.bFlag = -999;
	  pgData.bx = -999.;
	  pgData.by = -999.;
	  pgData.ebx = -999.;
	  pgData.eby = -999.;
	  
          netP = TVector3(0,0,0);
          thrust = TVector3(0,0,0);
          netP_charged = TVector3(0,0,0);
          thrust_charged = TVector3(0,0,0);
          nTrk=0;
          nTrk_GT0p4=0;
          nTrk_GT0p4Thrust=0;

	  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << ", " << genFileStr << ", " << counterEntries << ", " << evtNo.size() << std::endl;
	  if(doLocalDebug) std::cout << pgData.RunNo << ", " << pgData.EventNo << std::endl;
	  if(doLocalDebug) std::cout << runNo.at(counterEntries) << ", " << evtNo.at(counterEntries) << std::endl;
	  
	  if(pgData.RunNo != runNo.at(counterEntries) && pgData.EventNo != evtNo.at(counterEntries)){
	    std::cout << "Gen entries dont match reco for file \'" << genFileStr << "\'. return 1" << std::endl;
	    //gotta cleanup before return
	    delete tout;
	    delete bout;

	    for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	      delete jout[jIter];
	    }

	    if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	    if(isMC && isRecons){
	      delete tgout;
	      delete bgout;
	      for(int jIter = 0; jIter < nJtAlgo; ++jIter){
		delete jgout[jIter];
	      }
	    }

	    if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	    hf->Close();
	    delete hf;
	    
	    return 1;
	  }

	  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	  counterParticles=0;   

	  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	  
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
	else if(!isNewInfo && (num.size() == 7 || num.size() == 8)) assumePID = true; 
	else if(isNewInfo && (num.size() == 10 || num.size() == 11)) assumePID = true; 
	else{//return, this is an invalid format (or fix code here if format valid
	  std::cout << "Number of columns for line \'" << getStr << "\' is invalid, size \'" << num.size() << "\'. return 1" << std::endl;
	  //gotta cleanup before return
	  delete tout;
	  delete bout;
	  for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	    delete jout[jIter];
	  }

	  if(isMC && isRecons){
	    delete tgout;
	    delete bgout;
	    for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	      delete jgout[jIter];
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
	pgData.eta[counterParticles]=v.PseudoRapidity();
	pgData.theta[counterParticles]=v.Theta();
	pgData.phi[counterParticles]=v.Phi();

	if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	
	pgData.charge[counterParticles]=_charge;
	pgData.pwflag[counterParticles]=_pwflag;
	//check before assigning PID
	if(assumePID) pgData.pid[counterParticles]=std::stoi(num.at(6));
	else pgData.pid[counterParticles]=-999;

	if(isNewInfo && _pwflag < 4){
	  int posOffset = 0;
	  if(assumePID) posOffset += 1;
	  
	  pgData.d0[counterParticles] = std::stof(num.at(6 + posOffset));
	  pgData.z0[counterParticles] = std::stof(num.at(7 + posOffset));
	  pgData.ntpc[counterParticles] =  std::stoi(num.at(8 + posOffset));	
	}
	else{
	  pgData.d0[counterParticles] = -999.;
	  pgData.z0[counterParticles] = -999.;
	  pgData.ntpc[counterParticles] = -999;
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
        egData.missP = netP.Mag();
        egData.missPt = netP.Perp();
        egData.missTheta = netP.Theta();
        egData.missPhi = netP.Phi();
        egData.missChargedP = netP_charged.Mag();
        egData.missChargedPt = netP_charged.Perp();
        egData.missChargedTheta = netP_charged.Theta();
        egData.missChargedPhi = netP_charged.Phi();
        egData.nChargedHadrons = nTrk;
        egData.nChargedHadrons_GT0p4 = nTrk_GT0p4; 
        egData.nChargedHadrons_GT0p4Thrust = nTrk_GT0p4Thrust; 
      
	++counterParticles;	

	if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
      }
      //Have to fill one last time since the condition for fill is dependent on NEXT EVENT existing, else we would lose last event per file
      if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      pgData.nParticle=counterParticles;
      thrust = getThrust(pgData.nParticle, pgData.px, pgData.py, pgData.pz, THRUST::OPTIMAL); 
      thrust_charged = getChargedThrust(pgData.nParticle, pgData.px, pgData.py, pgData.pz, pgData.pwflag, THRUST::OPTIMAL );
      setThrustVariables(&pgData, &egData, thrust, thrust_charged);
      egData.Thrust = thrust.Mag();
      egData.TTheta = thrust.Theta();
      egData.TPhi = thrust.Phi();
      egData.Thrust_charged = thrust_charged.Mag();
      egData.TTheta_charged = thrust_charged.Theta();
      egData.TPhi_charged = thrust_charged.Phi();
      
      if(counterEntries>0) tgout->Fill(); 
      for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	processJets(particles, jDef[jIter], jDefReclust[jIter], &(jgData[jIter]), jtPtCut);
	if(counterEntries>0){
	  jgout[jIter]->Fill();
	}
        bgout->Fill();
      }
      
      if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      fileGen.close();
    }

    if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    
    runNo.clear();
    evtNo.clear();
  }

  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  hf->cd();
  
  tout->Write("", TObject::kOverwrite);
  delete tout;
  bout->Write("", TObject::kOverwrite);
  delete bout;

  for(int i = 0; i < nJtAlgo; ++i){
    jout[i]->Write("", TObject::kOverwrite);
    delete jout[i];
  }

  if(isMC && isRecons){
    tgout->Write("", TObject::kOverwrite);
    delete tgout;
    bgout->Write("", TObject::kOverwrite);
    delete bgout;

    for(int jIter = 0; jIter < nJtAlgo; ++jIter){
      jgout[jIter]->Write("", TObject::kOverwrite);
      delete jgout[jIter];
    }
  }
  
  hf->Close();
  delete hf;

  std::cout << "Job complete, return 0." << std::endl;

  return 0;
}

int main(int argc, char *argv[])
{
  if(argc != 3 && argc != 4){
    std::cout << "Usage: ./bin/scan.exe <inFileName> <isNewInfo> <OPT-outFileName>" << std::endl;
    return 1;
  }

  std::cout << "Begin processing..." << std::endl;
  int retVal = 0;
  if(argc == 3) retVal += scan(argv[1], std::stoi(argv[2]));
  else if(argc == 4) retVal += scan(argv[1], std::stoi(argv[2]), argv[3]);
  return retVal;
}
