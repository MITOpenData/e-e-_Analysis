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
#include "fastjet/PseudoJet.hh"

//local dependencies
#include "include/doLocalDebug.h"
#include "include/checkMakeDir.h"
#include "include/particleData.h"
#include "include/eventData.h"
#include "include/jetData.h"
#include "include/simpleJetMaker.h"
#include "include/thrustTools.h"

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
  const int nY=6;
  int years[nY] = {1995,1996,1997,1998,1999,2000};
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
  if(inStr.find("ALEPH_DATA") != std::string::npos) addDummy = true;

  const std::string goodChar = "0123456789-. ";
  unsigned int pos = 0;
  while(inStr.size() > pos){
    if(goodChar.find(inStr.substr(pos,1)) == std::string::npos) inStr.replace(pos,1,"");
    else ++pos;
  }

  std::vector<std::string> retV;
  if(addDummy){retV.push_back("-999."); retV.push_back("-999."); retV.push_back("-999."); retV.push_back("-999.");}

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

bool check999(std::string inStr)
{
  if(inStr.find("-999.") != std::string::npos && inStr.size() == 5) return true;
  return false;
}

void processJets(std::vector<fastjet::PseudoJet> p, simpleJetMaker j, jetData *d)
{
  d->nref = 0;
  if(p.size() > 0){
    std::vector<fastjet::PseudoJet> jets = j.getSimpleJets(p);
    for(unsigned int i = 0; i < jets.size(); ++i){
      if(jets.at(i).pt() < .1) break; //Arbitrarily low cut on jets, removes spike at phi zero when things become ill defined
      d->jtpt[d->nref] = jets.at(i).pt();
      d->jtphi[d->nref] = jets.at(i).phi_std();
      d->jteta[d->nref] = jets.at(i).eta();
      ++d->nref;
    }
    jets.clear();
    std::vector<fastjet::PseudoJet>().swap(jets);
  }
  return;
}

int scan(std::string inFileName, std::string outFileName="")
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
  const double rParam = 0.4;
  simpleJetMaker jMaker(rParam);

  const std::string partTreeName = "t";
  const std::string genPartTreeName = "tgen";
  const std::string jetTreeName = "ak" + std::to_string(int(rParam*10)) + "JetTree";
  const std::string genJetTreeName = "ak" + std::to_string(int(rParam*10)) + "GenJetTree";

  std::string finalPartTreeName = partTreeName;
  if(!isRecons) finalPartTreeName = genPartTreeName;
  std::string finalJetTreeName = jetTreeName;
  if(!isRecons) finalJetTreeName = genJetTreeName;

  TFile *hf = new TFile(outFileName.c_str(), "RECREATE");
  TTree *tout = new TTree(finalPartTreeName.c_str(), "");
  TTree *jout = new TTree(finalJetTreeName.c_str(), "");

  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  TTree *tgout=0;
  TTree *jgout=0;

  if(isRecons && isMC){
    tgout = new TTree(genPartTreeName.c_str(),"");
    jgout = new TTree(genJetTreeName.c_str(),"");
  }

  particleData pData;
  jetData jData;
  eventData eData;

  particleData pgData;
  jetData jgData;
  eventData egData;

  tout->Branch("year", &pData.year, "year/I");
  tout->Branch("EventNo", &pData.EventNo,"EventNo/I");
  tout->Branch("RunNo", &pData.RunNo,"RunNo/I");
  tout->Branch("Energy", &pData.Energy,"Energy/F");
  tout->Branch("process", &pData.process, "process/I");
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

  jout->Branch("nref", &jData.nref,"nref/I");
  jout->Branch("jtpt", jData.jtpt,"jtpt[nref]/F");
  jout->Branch("jteta", jData.jteta,"jteta[nref]/F");
  jout->Branch("jtphi", jData.jtphi,"jtphi[nref]/F");


  if(isRecons && isMC){
    tgout->Branch("year", &pgData.year, "year/I");
    tgout->Branch("EventNo", &pgData.EventNo,"EventNo/I");
    tgout->Branch("RunNo", &pgData.RunNo,"RunNo/I");
    tgout->Branch("Energy", &pgData.Energy,"Energy/F");
    tgout->Branch("process", &pgData.process, "process/I");
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
    
    jgout->Branch("nref", &jgData.nref,"nref/I");
    jgout->Branch("jtpt", jgData.jtpt,"jtpt[nref]/F");
    jgout->Branch("jteta", jgData.jteta,"jteta[nref]/F");
    jgout->Branch("jtphi", jgData.jtphi,"jtphi[nref]/F");
  }

  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    std::cout << "Processing \'" << fileList.at(fI) << "\'" << std::endl;
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
    TVector3 netP_charged = TVector3(0,0,0);
    int nTrk=0;
    int nTrk_GT0p4=0;
    int nTrk_GT0p4Thrust=0;

    while(std::getline(file,getStr)){
      if(getStr.size() == 0) continue;
      else if(getStr.find("******") != std::string::npos) continue;
      else if(getStr.find("END_EVENT") != std::string::npos) continue;
      else if(getStr.find("END_FILE") != std::string::npos) continue;
      std::vector<std::string> num = processAlephString(getStr);
      
      // check the number of columns before assigning values
      bool assumePID = false;
      if(num.size() == 6) assumePID = false; 
      else if(num.size() == 7 || num.size() == 8) assumePID = true; 
      else{//return, this is an invalid format (or fix code here if format valid
	std::cout << "Number of columns for line \'" << getStr << "\' is invalid, size \'" << num.size() << "\'. return 1" << std::endl;
	//gotta cleanup before return
	delete tout;
	delete jout;

	if(isMC && isRecons){
	  delete tgout;
	  delete jgout;
	}
	
	hf->Close();
	delete hf;
	
	return 1;
      }
      
      if(check999(num.at(0)) && check999(num.at(1)) && check999(num.at(2)) && check999(num.at(3))){ 
	pData.nParticle=counterParticles;
	if(counterEntries>0) tout->Fill(); 
	
	//Processing particles->jets
	processJets(particles, jMaker, &jData);
	if(counterEntries>0) jout->Fill();
	//clear particles for next iteration clustering
	particles.clear();

	unsigned int runPos = 4;
	unsigned int evtPos = 5;
	unsigned int ePos = 6;

	if(check999(num.at(4))){runPos++; evtPos++; ePos++;}

	pData.year = year;
	pData.process = process;
	pData.RunNo = std::stoi(num.at(runPos));
	pData.EventNo= std::stoi(num.at(evtPos));
	pData.Energy= std::stof(num.at(ePos));
      
	runNo.push_back(pData.RunNo);
	evtNo.push_back(pData.EventNo);

        netP = TVector3(0,0,0);
        netP_charged = TVector3(0,0,0);
        nTrk=0;
        nTrk_GT0p4=0;
        nTrk_GT0p4Thrust=0;

	counterParticles=0;   
	++counterEntries;	
	
	continue;
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
      particles.push_back(fastjet::PseudoJet(_px,_py,_pz,v.E()));
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
    //Have to fill one last time since the condition for fill is dependent on NEXT EVENT existing, else we would lose last event per file
    pData.nParticle=counterParticles;
    if(counterEntries>0) tout->Fill(); 
    processJets(particles, jMaker, &jData);
    if(counterEntries>0) jout->Fill();
    file.close();
  
    if(isRecons && isMC){
      std::string genFileStr = fileList.at(fI);
      genFileStr.replace(genFileStr.find("_recons_"), 8, "_mctrue_");

      std::ifstream fileGen(genFileStr.c_str());
      counterEntries=0;
      counterParticles=0;
      netP = TVector3(0,0,0);
      netP_charged = TVector3(0,0,0);
      nTrk=0;
      nTrk_GT0p4=0;
      nTrk_GT0p4Thrust=0;

      while(std::getline(fileGen,getStr)){
	if(getStr.size() == 0) continue;
	else if(getStr.find("******") != std::string::npos) continue;
	else if(getStr.find("END_EVENT") != std::string::npos) continue;
	else if(getStr.find("END_FILE") != std::string::npos) continue;
	std::vector<std::string> num = processAlephString(getStr);
      
	// check the number of columns before assigning values
	bool assumePID = false;
	if(num.size() == 6) assumePID = false; 
	else if(num.size() == 7 || num.size() == 8) assumePID = true; 
	else{//return, this is an invalid format (or fix code here if format valid
	  std::cout << "Number of columns for line \'" << getStr << "\' is invalid, size \'" << num.size() << "\'. return 1" << std::endl;
	  //gotta cleanup before return
	  delete tout;
	  delete jout;

	  if(isMC && isRecons){
	    delete tgout;
	    delete jgout;
	  }

	  hf->Close();
	  delete hf;
	
	  return 1;
	}

	if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
      
	if(check999(num.at(0)) && check999(num.at(1)) && check999(num.at(2)) && check999(num.at(3))){ 
	  pgData.nParticle=counterParticles;
	  if(counterEntries>0) tgout->Fill(); 
	  
	  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;	

	  //Processing particles->jets
	  processJets(particles, jMaker, &jgData);
	  if(counterEntries>0) jgout->Fill();
	  //clear particles for next iteration clustering
	  particles.clear();

	  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	  unsigned int runPos = 4;
	  unsigned int evtPos = 5;
	  unsigned int ePos = 6;
	  
	  if(check999(num.at(4))){runPos++; evtPos++; ePos++;}

	  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	  pgData.year = year;
	  pgData.process = process;
	  pgData.RunNo = std::stoi(num.at(runPos));
	  pgData.EventNo= std::stoi(num.at(evtPos));
	  pgData.Energy= std::stof(num.at(ePos));

          netP = TVector3(0,0,0);
          netP_charged = TVector3(0,0,0);
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
	    delete jout;

	    if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	    if(isMC && isRecons){
	      delete tgout;
	      delete jgout;
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
	particles.push_back(fastjet::PseudoJet(_px,_py,_pz,v.E()));
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
      if(counterEntries>0) tgout->Fill(); 
      processJets(particles, jMaker, &jgData);
      if(counterEntries>0) jgout->Fill();
      
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
  jout->Write("", TObject::kOverwrite);
  delete jout;

  if(isMC && isRecons){
    tgout->Write("", TObject::kOverwrite);
    delete tgout;
    jgout->Write("", TObject::kOverwrite);
    delete jgout;
  }
  
  hf->Close();
  delete hf;

  return 0;
}

int main(int argc, char *argv[])
{
  if(argc != 2 && argc != 3){
    std::cout << "Usage: ./bin/scan.exe <inFileName> <OPT-outFileName>" << std::endl;
    return 1;
  }

  std::cout << "Begin processing..." << std::endl;
  int retVal = 0;
  if(argc == 2) retVal += scan(argv[1]);
  else if(argc == 3) retVal += scan(argv[1], argv[2]);
  return retVal;
}
