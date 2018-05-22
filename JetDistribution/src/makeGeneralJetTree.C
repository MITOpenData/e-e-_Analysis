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
#include "TMatrixD.h"
#include "TMath.h"

//fastjet dependencies
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

//non-local project dependencies
#include "DataProcessing/include/particleData.h"
#include "DataProcessing/include/checkMakeDir.h"
#include "DataProcessing/include/thrustTools.h"
#include "DataProcessing/include/sphericityTools.h"

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
      
      if(tempStr.find(".root") != std::string::npos) fileList.push_back(tempStr);
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
  if(inName.find(".root") != std::string::npos) reducedInName.replace(reducedInName.find(".root"), 5, "");
  else if(inName.find(".txt") != std::string::npos) reducedInName.replace(reducedInName.find(".txt"), 4, "");
  
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  //define jetMaker here so rParam can be used in jetTree name for clarity

  const double jtPtCut = .01;
  const int nJtAlgo = 6;
  const double rParam[nJtAlgo] = {0.4, 0.8, 0.4, 0.8, -1, -1};
  const int ktParam[nJtAlgo] = {-1, -1, 1, 1, 100, 100};
  //  bool rParamIs8[nJtAlgo];
  //  for(int i = 0; i < nJtAlgo; ++i){rParamIs8[i] = rParam[i] > .799 && rParam[i] < .801;}
  const double recombScheme[nJtAlgo] = {fastjet::E_scheme, fastjet::E_scheme, fastjet::E_scheme, fastjet::E_scheme, fastjet::E_scheme, fastjet::E_scheme};
  fastjet::JetDefinition jDef[nJtAlgo];
  fastjet::JetDefinition jDefReclust[nJtAlgo];
  const int nFinalClust[nJtAlgo] = {-1, -1, -1, -1, 3, 4};

  for(int i = 0; i < nJtAlgo; ++i){
    if(rParam[i] > 0){
      jDef[i] = fastjet::JetDefinition(fastjet::ee_genkt_algorithm, rParam[i], ktParam[i], fastjet::RecombinationScheme(recombScheme[i]));
      jDefReclust[i] = fastjet::JetDefinition(fastjet::ee_genkt_algorithm, 10, ktParam[i], fastjet::RecombinationScheme(recombScheme[i]));
    }
    else{
      jDef[i] = fastjet::JetDefinition(fastjet::ee_kt_algorithm, fastjet::RecombinationScheme(recombScheme[i]));
      jDefReclust[i] = fastjet::JetDefinition(fastjet::ee_kt_algorithm, fastjet::RecombinationScheme(recombScheme[i]));
    }
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
    pData.SetStatusAndAddressRead(inTree_p, {"EventNo", "RunNo", "year", "subDir", "process", "isMC", "uniqueID", "nParticle", "px", "py", "pz", "mass"});

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

      gJetVar.EventNo = pData.EventNo;     
      gJetVar.RunNo = pData.RunNo;     
      gJetVar.year = pData.year;     
      gJetVar.subDir = pData.subDir;     
      gJetVar.process = pData.process;     
      gJetVar.isMC = pData.isMC;     
      gJetVar.uniqueID = pData.uniqueID;     
      gJetVar.nParticle = pData.nParticle;     
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
      gJetVar.oblateness = gJetVar.thrustMajorMag - gJetVar.thrustMinorMag;

      Sphericity spher(pData.nParticle, pData.px, pData.py, pData.pz, pData.pwflag, false);
      
      gJetVar.aplanarity = spher.aplanarity();
      gJetVar.planarity = spher.planarity();
      
      TLorentzVector hemiPos(0,0,0,0);
      TLorentzVector hemiNeg(0,0,0,0);
      double broadening1 = 0.;
      double broadening2 = 0.;
      double eVis = 0.;
      double pSum = 0.;
      double momTensor[3][3];

      for(Int_t i = 0; i < 3; ++i){
	for(Int_t j = 0; j < 3; ++j){
	  momTensor[i][j] = 0.;
	}
      }

      thrust.SetMag(1.);

      for(Int_t pI = 0; pI < pData.nParticle; ++pI){
	TLorentzVector particle;
	particle.SetXYZM(pData.px[pI], pData.py[pI], pData.pz[pI], pData.mass[pI]);
	particles.push_back(fastjet::PseudoJet(pData.px[pI], pData.py[pI], pData.pz[pI], particle.E()));
	TVector3 temp(pData.px[pI], pData.py[pI], pData.pz[pI]);     
	double tempArr[3] = {pData.px[pI], pData.py[pI], pData.pz[pI]};

	for(Int_t i = 0; i < 3; ++i){
	  for(Int_t j = 0; j < 3; ++j){
	    momTensor[i][j] += tempArr[i]*tempArr[j]/temp.Mag();
	  }
	}

	eVis += particle.E();
	if(temp.Dot(thrustMajor) > 0){
	  hemiPos += particle;
	  broadening1 += TMath::Abs(temp.Cross(thrust).Mag());
	}
	else{
	  hemiNeg += particle;
	  broadening2 += TMath::Abs(temp.Cross(thrust).Mag());
	}

	pSum += temp.Mag();
      }

      TMatrixD matrix(3,3);
      TVectorD eigenValues;
      for(Int_t i = 0; i < 3; ++i){
        for(Int_t j = 0; j < 3; ++j){
          matrix(i,j) = momTensor[i][j]/pSum;
        }
      }

      matrix = matrix.EigenVectors(eigenValues);


      gJetVar.eVis = eVis;
      gJetVar.heavyJetMass = hemiPos.E()*hemiPos.E() - (hemiPos.Px()*hemiPos.Px() + hemiPos.Py()*hemiPos.Py() + hemiPos.Pz()*hemiPos.Pz());
      gJetVar.heavyJetMass /= eVis*eVis;

      gJetVar.lightJetMass = hemiNeg.E()*hemiNeg.E() - (hemiNeg.Px()*hemiNeg.Px() + hemiNeg.Py()*hemiNeg.Py() + hemiNeg.Pz()*hemiNeg.Pz());
      gJetVar.lightJetMass /= eVis*eVis;

      if(gJetVar.heavyJetMass < gJetVar.lightJetMass){
	float tempVal = gJetVar.heavyJetMass;
	gJetVar.heavyJetMass = gJetVar.lightJetMass;
	gJetVar.lightJetMass = tempVal;
      }
      
      gJetVar.jetMassDifference = gJetVar.heavyJetMass - gJetVar.lightJetMass;

      broadening1 /= 2.*pSum;
      broadening2 /= 2.*pSum;

      gJetVar.wideJetBroadening = TMath::Max(broadening1, broadening2);
      gJetVar.narrowJetBroadening = TMath::Min(broadening1, broadening2);
      gJetVar.totalJetBroadening = gJetVar.wideJetBroadening + gJetVar.narrowJetBroadening;
      gJetVar.cParam = 3.*(eigenValues(0)*eigenValues(1) + eigenValues(1)*eigenValues(2) + eigenValues(2)*eigenValues(0));

      if(particles.size() == 0) continue;

      for(Int_t algoI = 0; algoI < nJtAlgo; ++algoI){
	fastjet::ClusterSequence cs(particles, jDef[algoI]);
	std::vector<fastjet::PseudoJet> jets;

	if(rParam[algoI] > 0) jets = fastjet::sorted_by_pt(cs.inclusive_jets());
	else jets = fastjet::sorted_by_pt(cs.exclusive_jets(nFinalClust[algoI]));

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

	//take exclusive clustering to four jets, and find the minimum distance parameters for jetResParam4
	//no procedure is listed in paper for this, but this seems correct - to be checked
	if(nFinalClust[algoI] == 4){
	  double minRes = 9999;
	  for(unsigned int i = 0; i < jets.size(); ++i){
	    double e1 = jets.at(i).E();
	    TVector3 v1(jets.at(i).px(), jets.at(i).py(), jets.at(i).pz());

	    for(unsigned int j = i+1; j < jets.size(); ++j){
	      double e2 = jets.at(j).E();
	      TVector3 v2(jets.at(j).px(), jets.at(j).py(), jets.at(j).pz());
	      double tempRes = 2.*TMath::Min(e1*e1,e2*e2)*(1. - v1.Dot(v2)/(v1.Mag()*v2.Mag()))/(eVis*eVis); //substituting dot product for cosine

	      if(tempRes < minRes) minRes = tempRes;
	    }
	  }
      
	  gJetVar.jetResParam4 = minRes;
	}     
      }

      gJetVar.jetResParam4NegLog = -TMath::Log(gJetVar.jetResParam4);

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
