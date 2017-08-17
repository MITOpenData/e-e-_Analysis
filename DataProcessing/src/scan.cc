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
#include "TNamed.h"

//fastjet dependencies
#include "fastjet/PseudoJet.hh"

//local dependencies
#include "include/checkMakeDir.h"
#include "include/particleData.h"
#include "include/jetData.h"
#include "include/simpleJetMaker.h"

std::vector<std::string> processAlephString(std::string inStr)
{
  std::vector<std::string> retV;

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

  // "/data/flowex/Datasamples/LEP2_MAIN/ROOTfiles/cleaned_ALEPH_Data-all.aleph"
  std::ifstream file(Form("%s",inFileName.c_str()));
  std::string getStr;
  
  if(outFileName.size() == 0){
    outFileName = inFileName + ".root";
    while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1,"");}
  }

  //define jetMaker here so rParam can be used in jetTree name for clarity
  const double rParam = 0.4;
  simpleJetMaker jMaker(rParam);
  const std::string jetTreeName = "ak" + std::to_string(int(rParam*10)) + "JetTree";

  TFile *hf = new TFile(outFileName.c_str(), "RECREATE");
  TTree *tout = new TTree("t","");
  TTree *jout = new TTree(jetTreeName.c_str(),"");

  particleData pData;
  jetData jData;

  tout->Branch("nParticle", &pData.nParticle,"nParticle/I");
  tout->Branch("EventNo", &pData.EventNo,"EventNo/I");
  tout->Branch("RunNo", &pData.RunNo,"RunNo/I");
  tout->Branch("Energy", &pData.Energy,"Energy/F");
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

  jout->Branch("nref", &jData.nref,"nref/I");
  jout->Branch("jtpt", jData.jtpt,"jtpt[nref]/F");
  jout->Branch("jteta", jData.jteta,"jteta[nref]/F");
  jout->Branch("jtphi", jData.jtphi,"jtphi[nref]/F");

  int counterEntries=0;
  int counterParticles=0;
  TLorentzVector v;

  std::vector<fastjet::PseudoJet> particles;

  while(std::getline(file,getStr)){
    if(getStr.size() == 0) continue;
    std::vector<std::string> num = processAlephString(getStr);
    if(check999(num.at(0)) && check999(num.at(1)) && check999(num.at(2))){ 
      pData.nParticle=counterParticles;

      if(counterEntries>0) tout->Fill(); 

      //Processing particles->jets
      processJets(particles, jMaker, &jData);
      if(counterEntries>0) jout->Fill();
      //clear particles for next iteration clustering
      particles.clear();

      pData.RunNo = std::stoi(num.at(4));
      pData.EventNo= std::stoi(num.at(5));
      pData.Energy= std::stof(num.at(6));

      counterParticles=0;   

      continue;
    }

    float _px = std::stof(num.at(0));
    float _py = std::stof(num.at(1));
    float _pz = std::stof(num.at(2));
    float _m = std::stof(num.at(3));
    float _charge = std::stof(num.at(4));
    int _pwflag = std::stoi(num.at(5));

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
    pData.pid[counterParticles]=0;
    ++counterParticles;	
    ++counterEntries;	
  }
  //Have to fill one last time since the condition for fill is dependent on NEXT EVENT existing, else we would lose last event per file
  pData.nParticle=counterParticles;
  if(counterEntries>0) tout->Fill(); 
  processJets(particles, jMaker, &jData);
  if(counterEntries>0) jout->Fill();

  hf->cd();

  tout->Write("", TObject::kOverwrite);
  delete tout;
  jout->Write("", TObject::kOverwrite);
  delete jout;

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
