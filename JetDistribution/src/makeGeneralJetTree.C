//A base for more generalized jet measurements (initial author Chris McGinn)
//Borrows heavily from base code found in DataProcessing (Most implementations of thrust by Austin Baty or taken from ALEPH/BELLE experiment code)
//Some general resources for what we measure:
//  http://cds.cern.ch/record/690637/files/ep-2003-084.pdf (ALEPH QCD at Ecm 91 and 209 GeV)
//  

//standard dependencies
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

//root dependencies
#include "TFile.h"
#include "TTree.h"
#include "TDatime.h"
#include "TVector3.h"
#include "TLorentzVector.h"

//fastjet dependencies
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

//non-local project dependencies
#include "DataProcessing/include/particleData.h"
#include "DataProcessing/include/checkMakeDir.h"
#include "DataProcessing/include/thrustTools.h"

//local project dependencies
#include "JetDistribution/include/globalJetVar.h"
#include "JetDistribution/include/generalJetVar.h"

void initTVector3(TVector3* in){(*in) = TVector3(0,0,0); return;}
void initTVector3(std::vector<TVector3*> in)
{
  for(unsigned int i = 0; i < in.size(); ++i){initTVector3(in.at(i));}
  return;
}

int makeGeneralJetTree(const std::string inName)
{
  std::vector<std::string> fileList;
  if(inName.find(".root") != std::string::npos) fileList.push_back(inName);
  else if(inName.find(".txt") != std::string::npos){
    std::ifstream file(inName.c_str());
    std::string tempStr;

    while(std::getline(file, tempStr)){
      if(tempStr.size() == 0) continue;
      
      if(tempStr.find(".root") != std::string::npos) fileList.push_back(inName);
    }

    file.close();
  }
  else{
    std::cout << "inName \'" << inName << "\' is not valid. return 1" << std::endl;
    return 1;
  }

  if(fileList.size() == 0){
    std::cout << "inName \'" << inName << "\' contains no valid root files. return 1" << std::endl;
    return 1;
  }

  std::string reducedInName = inName;
  while(reducedInName.find("/") != std::string::npos){reducedInName.replace(0, reducedInName.find("/")+1, "");}
  reducedInName.replace(reducedInName.find(".root"), 5, "");
  
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  //define jetMaker here so rParam can be used in jetTree name for clarity

  const double jtPtCut = .01;
  const int nJtAlgo = 4;
  const double rParam[nJtAlgo] = {0.4, 0.8, 0.4, 0.8};
  const int ktParam[nJtAlgo] = {-1, -1, 1, 1};
  //  bool rParamIs8[nJtAlgo];
  //  for(int i = 0; i < nJtAlgo; ++i){rParamIs8[i] = rParam[i] > .799 && rParam[i] < .801;}
  const double recombScheme[nJtAlgo] = {fastjet::E_scheme, fastjet::E_scheme, fastjet::E_scheme, fastjet::E_scheme};
  fastjet::JetDefinition jDef[nJtAlgo];
  fastjet::JetDefinition jDefReclust[nJtAlgo];
  const int nFinalClust[nJtAlgo] = {-1, -1, -1, -1};

  for(int i = 0; i < nJtAlgo; ++i){
    jDef[i] = fastjet::JetDefinition(fastjet::ee_genkt_algorithm, rParam[i], ktParam[i], fastjet::RecombinationScheme(recombScheme[i]));
    jDefReclust[i] = fastjet::JetDefinition(fastjet::ee_genkt_algorithm, 10, ktParam[i], fastjet::RecombinationScheme(recombScheme[i]));
  }

  checkMakeDir("output");
  const std::string outFileName = "output/" + reducedInName + "_GeneralJetTree_" + dateStr + ".root";
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TTree* globalOutTree_p = new TTree("generalJetTree_Global", "");
  TTree* outTree_p[nJtAlgo];
  globalJetVar gJetVar;
  generalJetVar jetVar[nJtAlgo];

  gJetVar.SetBranchWrite(globalOutTree_p);

  for(Int_t jI = 0; jI < nJtAlgo; ++jI){
    std::string recombSchemeStr = "EScheme";
    if(recombScheme[jI] == fastjet::WTA_modp_scheme) recombSchemeStr = "WTAmodpScheme";

    std::string ktOrAntiKt = "ak";
    if(ktParam[jI] == 1 || rParam[jI] < 0) ktOrAntiKt = "kt";

    std::string rOrJetStr = ktOrAntiKt + "R" + std::to_string(int(rParam[jI]*10));
    if(rParam[jI] < 0) rOrJetStr = ktOrAntiKt + "N" + std::to_string(nFinalClust[jI]);

    const std::string jetTreeName = "generalJetTree_" + rOrJetStr + recombSchemeStr; 
    outTree_p[jI] = new TTree(jetTreeName.c_str(), "");
    jetVar[jI].SetBranchWrite(outTree_p[jI]);
  }

  particleData pData;

  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    std::cout << "File " << fI << "/" << fileList.size() << ": " << fileList.at(fI) << std::endl;

    TFile* inFile_p = new TFile(fileList.at(fI).c_str(), "READ");
    TTree* inTree_p = (TTree*)inFile_p->Get("t");
    pData.SetStatusAndAddressRead(inTree_p, {"nParticle", "px", "py", "pz", "mass"});

    const Int_t nEntries = inTree_p->GetEntries();

    for(Int_t entry = 0; entry < nEntries; ++entry){
      if(entry%10000 == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;
      inTree_p->GetEntry(entry);

      std::vector<fastjet::PseudoJet> particles;

      TVector3 thrust, thrustMajor, thrustMinor;
      initTVector3({&thrust, &thrustMajor, &thrustMinor});
      thrust = getThrust(pData.nParticle, pData.px, pData.py, pData.pz, THRUST::OPTIMAL);      
      thrustMajor = getThrustMajor(thrust, pData.nParticle, pData.px, pData.py, pData.pz, NULL, THRUST::OPTIMAL, false);
      thrustMinor = getThrustMinor(thrust, thrustMajor, pData.nParticle, pData.px, pData.py, pData.pz, NULL, THRUST::OPTIMAL, false);

      gJetVar.thrustMag = thrust.Mag();
      gJetVar.thrustPx = thrust.Px();
      gJetVar.thrustPy = thrust.Py();
      gJetVar.thrustPz = thrust.Pz();
      gJetVar.thrustMajorMag = thrustMajor.Mag();
      gJetVar.thrustMajorPx = thrustMajor.Px();
      gJetVar.thrustMajorPy = thrustMajor.Py();
      gJetVar.thrustMajorPz = thrustMajor.Pz();
      gJetVar.thrustMinorMag = thrustMinor.Mag();
      gJetVar.thrustMinorPx = thrustMinor.Px();
      gJetVar.thrustMinorPy = thrustMinor.Py();
      gJetVar.thrustMinorPz = thrustMinor.Pz();
      
      for(Int_t pI = 0; pI < pData.nParticle; ++pI){
	TLorentzVector particle;
	particle.SetXYZM(pData.px[pI], pData.py[pI], pData.pz[pI], pData.mass[pI]);
	particles.push_back(fastjet::PseudoJet(pData.px[pI], pData.py[pI], pData.pz[pI], particle.E()));
      }

      if(particles.size() == 0) continue;

      for(Int_t algoI = 0; algoI < nJtAlgo; ++algoI){
	fastjet::ClusterSequence cs(particles, jDef[algoI]);
	std::vector<fastjet::PseudoJet> jets;
	jets = fastjet::sorted_by_pt(cs.inclusive_jets());

	jetVar[algoI].nref = 0;

	for(unsigned int i = 0; i < jets.size(); ++i){
	  if(jets.at(i).pt() < jtPtCut && rParam[algoI] > 0) break; //Arbitrarily low cut on jets, removes spike at phi zero when things become ill defined

	  jetVar[algoI].jtpt[jetVar[algoI].nref] = jets.at(i).pt();
	  jetVar[algoI].jtphi[jetVar[algoI].nref] = jets.at(i).phi_std();
	  jetVar[algoI].jtm[jetVar[algoI].nref] = jets.at(i).m();
	  jetVar[algoI].jteta[jetVar[algoI].nref] = jets.at(i).eta();
	  jetVar[algoI].fourJet[jetVar[algoI].nref] = TLorentzVector(jets.at(i).px(), jets.at(i).py(), jets.at(i).pz(), jets.at(i).E());
	  
	  ++(jetVar[algoI].nref);
	}
      }

      globalOutTree_p->Fill();
      for(Int_t jI = 0; jI < nJtAlgo; ++jI){
	jetVar[jI].preFillClean();
	outTree_p[jI]->Fill();
      }
    }
    
    inFile_p->Close();
    delete inFile_p;
  }

  outFile_p->cd();

  globalOutTree_p->Write("", TObject::kOverwrite);
  delete globalOutTree_p;
  for(Int_t jI = 0; jI < nJtAlgo; ++jI){
    outTree_p[jI]->Write("", TObject::kOverwrite);
    delete outTree_p[jI];
  }

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/makeGeneralJetTree.exe <inName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += makeGeneralJetTree(argv[1]);
  return retVal;
}
