#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TDatime.h"
#include "TVector3.h"
#include "TMath.h"

#include "DataProcessing/include/checkMakeDir.h"
#include "DataProcessing/include/returnRootFileContentsList.h"
#include "DataProcessing/include/particleData.h"
#include "DataProcessing/include/eventData.h"
#include "DataProcessing/include/removeVectorDuplicates.h"
#include "DataProcessing/include/thrustTools.h"
#include "DataProcessing/include/histDefUtility.h"

#include "TwoParticleCorrelation/include/getLogBins.h"
#include "TwoParticleCorrelation/include/getLinBins.h"

int sigMixThrustComp(const std::string inFileName1, const std::string inFileName2)
{
  if(!checkFile(inFileName1)){
    std::cout << "Given inFileName1 \'" << inFileName1 << "\' is invalid. return 1" << std::endl;
    return 1;
  }
  else if(!checkFile(inFileName2)){
    std::cout << "Given inFileName2 \'" << inFileName2 << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  std::vector<std::string> fileList1;
  std::vector<std::string> fileList2;

  if(inFileName1.find(".txt") != std::string::npos){
    std::ifstream file(inFileName1.c_str());
    std::string tempStr;
    while(std::getline(file, tempStr)){
      if(tempStr.size() == 0) continue;
      if(tempStr.find(".root") != std::string::npos && checkFile(tempStr)) fileList1.push_back(tempStr);
    }
    file.close();
  }
  else if(inFileName1.find(".root") != std::string::npos) fileList1.push_back(inFileName1);
  else{
    std::cout << "Given inFileName1 \'" << inFileName1 << "\' is not \'.txt\' or \'.root\'. return 1" << std::endl;
    return 1;
  }

  if(inFileName2.find(".txt") != std::string::npos){
    std::ifstream file(inFileName2.c_str());
    std::string tempStr;
    while(std::getline(file, tempStr)){
      if(tempStr.size() == 0) continue;
      if(tempStr.find(".root") != std::string::npos && checkFile(tempStr)) fileList2.push_back(tempStr);
    }
    file.close();
  }
  else if(inFileName2.find(".root") != std::string::npos) fileList2.push_back(inFileName2);
  else{
    std::cout << "Given inFileName2 \'" << inFileName2 << "\' is not \'.txt\' or \'.root\'. return 1" << std::endl;
    return 1;
  }



  TFile* inFile1_p = new TFile(fileList1.at(0).c_str(), "READ");
  std::vector<std::string> ttreeList1 = returnRootFileContentsList(inFile1_p, "TTree", "");
  inFile1_p->Close();
  delete inFile1_p;
  removeVectorDuplicates(&ttreeList1);

  TFile* inFile2_p = new TFile(fileList2.at(0).c_str(), "READ");
  std::vector<std::string> ttreeList2 = returnRootFileContentsList(inFile2_p, "TTree", "");
  inFile2_p->Close();
  delete inFile2_p;
  removeVectorDuplicates(&ttreeList2);

  bool containsTGen1 = false;
  bool containsTGen2 = false;

  for(unsigned int tI = 0; tI < ttreeList1.size(); ++tI){
    //    std::cout << "TTree " << tI << "/" << ttreeList1.size() << ": " << ttreeList1.at(tI) << std::endl;
    if(ttreeList1.at(tI).find("tgen") != std::string::npos){
      containsTGen1 = true;
      break;
    }
  }

  for(unsigned int tI = 0; tI < ttreeList2.size(); ++tI){
    //    std::cout << "TTree " << tI << "/" << ttreeList2.size() << ": " << ttreeList2.at(tI) << std::endl;
    if(ttreeList2.at(tI).find("tgen") != std::string::npos){
      containsTGen2 = true;
      break;
    }
  }

  /*
  if(!containsTGen1){
    std::cout << "inFile1 \'" << inFileName1 << "\' does not contain tgen TTree. return 1" << std::endl;

    inFile1_p->Close();
    delete inFile1_p;
    inFile2_p->Close();
    delete inFile2_p;
  }
  else if(!containsTGen2){
    std::cout << "inFile2 \'" << inFileName2 << "\' does not contain tgen TTree. return 1" << std::endl;

    inFile1_p->Close();
    delete inFile1_p;
    inFile2_p->Close();
    delete inFile2_p;
  }
  */

  bool containsTGen = containsTGen1 || containsTGen2;

  TDatime* date = new TDatime();
  std::string dataMCStr = "Data";
  if(containsTGen) dataMCStr = "MC";
 
  if(inFileName1.find("Pythia") != std::string::npos || inFileName1.find("PYTHIA") != std::string::npos || inFileName1.find("pythia") != std::string::npos) dataMCStr = "PYTHIA";

  TFile* outFile_p = new TFile(("output/sigMixThrustComp_" + dataMCStr + "_" + std::to_string(date->GetDate()) + ".root").c_str(), "RECREATE");
  delete date;

  const Int_t nMultBins = 3;
  const Int_t multBinsLow[nMultBins] = {0, 20, 30};
  const Int_t multBinsHi[nMultBins] = {20, 30, 200};  

  const Int_t nPtBins = 20;
  Float_t ptLow = 0.2;
  Float_t ptHi = 30.;
  Double_t ptBins[nPtBins+1];
  getLogBins(ptLow, ptHi, nPtBins, ptBins);

  const Int_t nEtaBins = 50;
  //  Float_t etaLow = -10.;
  //  Float_t etaHi = 10.;
  Double_t etaBins[nEtaBins+1] = {-10., -6.000, -5.5, -5.0, -4.5, -4.0, -3.800, -3.600, -3.400, -3.200, -3.000, -2.800, -2.600, -2.400, -2.200, -2.000, -1.800, -1.600, -1.400, -1.200, -1.000, -0.800, -0.600, -0.400, -0.200, 0.000, 0.200, 0.400, 0.600, 0.800, 1.000, 1.200, 1.400, 1.600, 1.800, 2.000, 2.200, 2.400, 2.600, 2.800, 3.000, 3.200, 3.400, 3.600, 3.800, 4.0, 4.5, 5.0, 5.5, 6.000, 10.};

  const Int_t nAbsEtaBins = 25;
  Double_t absEtaBins[nAbsEtaBins+1] = {0.0, 0.200, 0.400, 0.600, 0.800, 1.000, 1.200, 1.400, 1.600, 1.800, 2.000, 2.200, 2.400, 2.600, 2.800, 3.000, 3.200, 3.400, 3.600, 3.800, 4.0, 4.5, 5.0, 5.5, 6.000, 10.};

  const Int_t nRapBins = 100;
  Float_t rapLow = -10.;
  Float_t rapHi = 10.;

  const Int_t nThetaBins = 20;
  Float_t thetaLow = 0;
  Float_t thetaHi = TMath::Pi();

  const Int_t nPhiBins = 40;
  Float_t phiLow = -TMath::Pi()-.01;
  Float_t phiHi = TMath::Pi()+.01;
  Double_t phiBins[nPhiBins+1];
  getLinBins(phiLow, phiHi, nPhiBins, phiBins);

  TH1D* sigRecoPtWrtRecoThr_h[nMultBins+1];
  TH1D* sigGenPtWrtGenThr_h[nMultBins+1];
  TH1D* mixRecoPtWrtRecoThr_h[nMultBins+1];
  TH1D* mixGenPtWrtGenThr_h[nMultBins+1];

  TH1D* sigRecoEtaWrtRecoThr_h[nMultBins+1];
  TH1D* sigGenEtaWrtGenThr_h[nMultBins+1];
  TH1D* mixRecoEtaWrtRecoThr_h[nMultBins+1];
  TH1D* mixGenEtaWrtGenThr_h[nMultBins+1];

  TH1D* sigRecoAbsEtaWrtRecoThr_h[nMultBins+1];
  TH1D* sigGenAbsEtaWrtGenThr_h[nMultBins+1];
  TH1D* mixRecoAbsEtaWrtRecoThr_h[nMultBins+1];
  TH1D* mixGenAbsEtaWrtGenThr_h[nMultBins+1];

  TH1D* sigRecoRapWrtRecoThr_h[nMultBins+1];
  TH1D* sigGenRapWrtGenThr_h[nMultBins+1];
  TH1D* mixRecoRapWrtRecoThr_h[nMultBins+1];
  TH1D* mixGenRapWrtGenThr_h[nMultBins+1];

  TH1D* sigRecoThetaWrtRecoThr_h[nMultBins+1];
  TH1D* sigGenThetaWrtGenThr_h[nMultBins+1];
  TH1D* mixRecoThetaWrtRecoThr_h[nMultBins+1];
  TH1D* mixGenThetaWrtGenThr_h[nMultBins+1];

  TH1D* sigRecoPhiWrtRecoThr_h[nMultBins+1];
  TH1D* sigGenPhiWrtGenThr_h[nMultBins+1];
  TH1D* mixRecoPhiWrtRecoThr_h[nMultBins+1];
  TH1D* mixGenPhiWrtGenThr_h[nMultBins+1];   

  TH2D* sigRecoEtaPhiWrtRecoThr_h[nMultBins+1];
  TH2D* sigGenEtaPhiWrtGenThr_h[nMultBins+1];
  TH2D* mixRecoEtaPhiWrtRecoThr_h[nMultBins+1];
  TH2D* mixGenEtaPhiWrtGenThr_h[nMultBins+1];

  TH2D* sigRecoAbsEtaPhiWrtRecoThr_h[nMultBins+1];
  TH2D* sigGenAbsEtaPhiWrtGenThr_h[nMultBins+1];
  TH2D* mixRecoAbsEtaPhiWrtRecoThr_h[nMultBins+1];
  TH2D* mixGenAbsEtaPhiWrtGenThr_h[nMultBins+1];

  for(Int_t mI = 0; mI < nMultBins+1; ++mI){
    std::string multStr = "Mult" + std::to_string(multBinsLow[0]) + "to" + std::to_string(multBinsHi[nMultBins-1]);
    if(mI < nMultBins) multStr = "Mult" + std::to_string(multBinsLow[mI]) + "to" + std::to_string(multBinsHi[mI]);

    std::cout << multStr << ", " << mI << "/" << nMultBins << std::endl;

    sigRecoPtWrtRecoThr_h[mI] = new TH1D(("sigRecoPtWrtRecoThr_" + multStr + "_h").c_str(), ";Sig. reco. p_{T} w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{dp_{T}}", nPtBins, ptBins);
    sigGenPtWrtGenThr_h[mI] = new TH1D(("sigRecoPtWrtGenThr_" + multStr + "_h").c_str(), ";Sig. gen. p_{T} w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{dp_{T}}", nPtBins, ptBins);
    mixRecoPtWrtRecoThr_h[mI] = new TH1D(("mixRecoPtWrtRecoThr_" + multStr + "_h").c_str(), ";Mix. reco. p_{T} w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{dp_{T}}", nPtBins, ptBins);
    mixGenPtWrtGenThr_h[mI] = new TH1D(("mixRecoPtWrtGenThr_" + multStr + "_h").c_str(), ";Mix. gen. p_{T} w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{dp_{T}}", nPtBins, ptBins);
    
    sigRecoEtaWrtRecoThr_h[mI] = new TH1D(("sigRecoEtaWrtRecoThr_" + multStr + "_h").c_str(), ";Sig. reco. #eta w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#eta}", nEtaBins, etaBins);
    sigGenEtaWrtGenThr_h[mI] = new TH1D(("sigRecoEtaWrtGenThr_" + multStr + "_h").c_str(), ";Sig. gen. #eta w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#eta}", nEtaBins, etaBins);
    mixRecoEtaWrtRecoThr_h[mI] = new TH1D(("mixRecoEtaWrtRecoThr_" + multStr + "_h").c_str(), ";Mix. reco. #eta w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#eta}", nEtaBins, etaBins);
    mixGenEtaWrtGenThr_h[mI] = new TH1D(("mixRecoEtaWrtGenThr_" + multStr + "_h").c_str(), ";Mix. gen. #eta w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#eta}", nEtaBins, etaBins);

    sigRecoAbsEtaWrtRecoThr_h[mI] = new TH1D(("sigRecoAbsEtaWrtRecoThr_" + multStr + "_h").c_str(), ";Sig. reco. |#eta| w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d|#eta|}", nAbsEtaBins, absEtaBins);
    sigGenAbsEtaWrtGenThr_h[mI] = new TH1D(("sigRecoAbsEtaWrtGenThr_" + multStr + "_h").c_str(), ";Sig. gen. |#eta| w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d|#eta|}", nAbsEtaBins, absEtaBins);
    mixRecoAbsEtaWrtRecoThr_h[mI] = new TH1D(("mixRecoAbsEtaWrtRecoThr_" + multStr + "_h").c_str(), ";Mix. reco. |#eta| w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d|#eta|}", nAbsEtaBins, absEtaBins);
    mixGenAbsEtaWrtGenThr_h[mI] = new TH1D(("mixRecoAbsEtaWrtGenThr_" + multStr + "_h").c_str(), ";Mix. gen. |#eta| w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d|#eta|}", nAbsEtaBins, absEtaBins);
    
    sigRecoRapWrtRecoThr_h[mI] = new TH1D(("sigRecoRapWrtRecoThr_" + multStr + "_h").c_str(), ";Sig. reco. y w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{dy}", nRapBins, rapLow, rapHi);
    sigGenRapWrtGenThr_h[mI] = new TH1D(("sigRecoRapWrtGenThr_" + multStr + "_h").c_str(), ";Sig. gen. y w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{dy}", nRapBins, rapLow, rapHi);
    mixRecoRapWrtRecoThr_h[mI] = new TH1D(("mixRecoRapWrtRecoThr_" + multStr + "_h").c_str(), ";Mix. reco. y w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{dy}", nRapBins, rapLow, rapHi);
    mixGenRapWrtGenThr_h[mI] = new TH1D(("mixRecoRapWrtGenThr_" + multStr + "_h").c_str(), ";Mix. gen. y w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{dy}", nRapBins, rapLow, rapHi);
    
    sigRecoThetaWrtRecoThr_h[mI] = new TH1D(("sigRecoThetaWrtRecoThr_" + multStr + "_h").c_str(), ";Sig. reco. #theta w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#theta}", nThetaBins, thetaLow, thetaHi);
    sigGenThetaWrtGenThr_h[mI] = new TH1D(("sigRecoThetaWrtGenThr_" + multStr + "_h").c_str(), ";Sig. gen. #theta w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#theta}", nThetaBins, thetaLow, thetaHi);
    mixRecoThetaWrtRecoThr_h[mI] = new TH1D(("mixRecoThetaWrtRecoThr_" + multStr + "_h").c_str(), ";Mix. reco. #theta w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#theta}", nThetaBins, thetaLow, thetaHi);
    mixGenThetaWrtGenThr_h[mI] = new TH1D(("mixRecoThetaWrtGenThr_" + multStr + "_h").c_str(), ";Mix. gen. #theta w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#theta}", nThetaBins, thetaLow, thetaHi);
    
    sigRecoPhiWrtRecoThr_h[mI] = new TH1D(("sigRecoPhiWrtRecoThr_" + multStr + "_h").c_str(), ";Sig. reco. #phi w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#phi}", nPhiBins, phiLow, phiHi);
    sigGenPhiWrtGenThr_h[mI] = new TH1D(("sigRecoPhiWrtGenThr_" + multStr + "_h").c_str(), ";Sig. gen. #phi w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#phi}", nPhiBins, phiLow, phiHi);
    mixRecoPhiWrtRecoThr_h[mI] = new TH1D(("mixRecoPhiWrtRecoThr_" + multStr + "_h").c_str(), ";Mix. reco. #phi w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#phi}", nPhiBins, phiLow, phiHi);
    mixGenPhiWrtGenThr_h[mI] = new TH1D(("mixRecoPhiWrtGenThr_" + multStr + "_h").c_str(), ";Mix. gen. #phi w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#phi}", nPhiBins, phiLow, phiHi);


    sigRecoEtaPhiWrtRecoThr_h[mI] = new TH2D(("sigRecoEtaPhiWrtRecoThr_" + multStr + "_h").c_str(), ";Sig. reco. #eta w.r.t. Thrust;Sig. reco. #phi w.r.t. Thrust", nEtaBins, etaBins, nPhiBins, phiBins);
    sigGenEtaPhiWrtGenThr_h[mI] = new TH2D(("sigGenEtaPhiWrtGenThr_" + multStr + "_h").c_str(), ";Sig. reco. #eta w.r.t. Thrust;Sig. reco. #phi w.r.t. Thrust", nEtaBins, etaBins, nPhiBins, phiBins);
    mixRecoEtaPhiWrtRecoThr_h[mI] = new TH2D(("mixRecoEtaPhiWrtRecoThr_" + multStr + "_h").c_str(), ";Mix reco. #eta w.r.t. Thrust;Mix reco. #phi w.r.t. Thrust", nEtaBins, etaBins, nPhiBins, phiBins);
    mixGenEtaPhiWrtGenThr_h[mI] = new TH2D(("mixGenEtaPhiWrtGenThr_" + multStr + "_h").c_str(), ";Mix gen. #eta w.r.t. Thrust;Mix gen. #phi w.r.t. Thrust", nEtaBins, etaBins, nPhiBins, phiBins);

    sigRecoAbsEtaPhiWrtRecoThr_h[mI] = new TH2D(("sigRecoAbsEtaPhiWrtRecoThr_" + multStr + "_h").c_str(), ";Sig. reco. #eta w.r.t. Thrust;Sig. reco. #phi w.r.t. Thrust", nAbsEtaBins, absEtaBins, nPhiBins, phiBins);
    sigGenAbsEtaPhiWrtGenThr_h[mI] = new TH2D(("sigGenAbsEtaPhiWrtGenThr_" + multStr + "_h").c_str(), ";Sig. reco. #eta w.r.t. Thrust;Sig. reco. #phi w.r.t. Thrust", nAbsEtaBins, absEtaBins, nPhiBins, phiBins);
    mixRecoAbsEtaPhiWrtRecoThr_h[mI] = new TH2D(("mixRecoAbsEtaPhiWrtRecoThr_" + multStr + "_h").c_str(), ";Mix reco. #eta w.r.t. Thrust;Mix reco. #phi w.r.t. Thrust", nAbsEtaBins, absEtaBins, nPhiBins, phiBins);
    mixGenAbsEtaPhiWrtGenThr_h[mI] = new TH2D(("mixGenAbsEtaPhiWrtGenThr_" + multStr + "_h").c_str(), ";Mix gen. #eta w.r.t. Thrust;Mix gen. #phi w.r.t. Thrust", nAbsEtaBins, absEtaBins, nPhiBins, phiBins);
  }

  particleData pSigData;
  particleData pSigDataGen;
  particleData pMixData;
  particleData pMixDataGen;

  for(unsigned int fI = 0; fI < fileList1.size(); ++fI){
    std::cout << "File  " << fI << "/" << fileList1.size() << std::endl;

    inFile1_p = new TFile(fileList1.at(fI).c_str(), "READ");
    inFile2_p = new TFile(fileList2.at(fI).c_str(), "READ");

    TTree* tsig = (TTree*)inFile1_p->Get("t");
    TTree* tsiggen = NULL;
    if(containsTGen) tsiggen = (TTree*)inFile1_p->Get("tgen");
    
    TTree* tmix = (TTree*)inFile2_p->Get("t");
    TTree* tmixgen = NULL;
    if(containsTGen) tmixgen = (TTree*)inFile2_p->Get("tgen");
    
    pSigData.SetStatusAndAddressRead(tsig, {});
    if(containsTGen) pSigDataGen.SetStatusAndAddressRead(tsiggen, {});
    
    pMixData.SetStatusAndAddressRead(tmix, {});
    if(containsTGen) pMixDataGen.SetStatusAndAddressRead(tmixgen, {});
    
    const Int_t nEntries = tsig->GetEntries();
    //  const Double_t weight = 1./(Double_t)nEntries;
    
    for(Int_t entry = 0; entry < nEntries; ++entry){
      if(entry%10000 == 0) std::cout << " Entry: " << entry << "/" << nEntries << std::endl;
      //    if(entry == 100000) break;
      
      tsig->GetEntry(entry);
      if(containsTGen) tsiggen->GetEntry(entry);
      
      tmix->GetEntry(entry);
      if(containsTGen) tmixgen->GetEntry(entry);
      
      std::vector<int> pos;
      pos.push_back(nMultBins);
      for(Int_t mI = 0; mI < nMultBins; ++mI){
	if(pSigData.nParticle >= multBinsLow[mI] && pSigData.nParticle < multBinsHi[mI]){
	  pos.push_back(mI);
	  break;
	}
      }

      if(pos.size() != 2) std::cout << "ERROR" << std::endl;
    
      for(Int_t pI = 0; pI < pSigData.nParticle; ++pI){
	for(unsigned int mI = 0; mI < pos.size(); ++mI){
	  sigRecoPtWrtRecoThr_h[pos.at(mI)]->Fill(pSigData.pt_wrtThr[pI]);
	  sigRecoEtaWrtRecoThr_h[pos.at(mI)]->Fill(pSigData.eta_wrtThr[pI]);
	  sigRecoAbsEtaWrtRecoThr_h[pos.at(mI)]->Fill(TMath::Abs(pSigData.eta_wrtThr[pI]));
	  sigRecoRapWrtRecoThr_h[pos.at(mI)]->Fill(pSigData.rap_wrtThr[pI]);
	  sigRecoThetaWrtRecoThr_h[pos.at(mI)]->Fill(pSigData.theta_wrtThr[pI]);
	  sigRecoPhiWrtRecoThr_h[pos.at(mI)]->Fill(pSigData.phi_wrtThr[pI]);
	  sigRecoEtaPhiWrtRecoThr_h[pos.at(mI)]->Fill(pSigData.eta_wrtThr[pI], pSigData.phi_wrtThr[pI]);
	  sigRecoAbsEtaPhiWrtRecoThr_h[pos.at(mI)]->Fill(TMath::Abs(pSigData.eta_wrtThr[pI]), pSigData.phi_wrtThr[pI]);
	}
      }
      
      if(containsTGen){
	for(Int_t pI = 0; pI < pSigDataGen.nParticle; ++pI){
	  for(unsigned int mI = 0; mI < pos.size(); ++mI){
	    sigGenPtWrtGenThr_h[pos.at(mI)]->Fill(pSigDataGen.pt_wrtThr[pI]);
	    sigGenEtaWrtGenThr_h[pos.at(mI)]->Fill(pSigDataGen.eta_wrtThr[pI]);
	    sigGenAbsEtaWrtGenThr_h[pos.at(mI)]->Fill(TMath::Abs(pSigDataGen.eta_wrtThr[pI]));
	    sigGenRapWrtGenThr_h[pos.at(mI)]->Fill(pSigDataGen.rap_wrtThr[pI]);
	    sigGenThetaWrtGenThr_h[pos.at(mI)]->Fill(pSigDataGen.theta_wrtThr[pI]);
	    sigGenPhiWrtGenThr_h[pos.at(mI)]->Fill(pSigDataGen.phi_wrtThr[pI]);
	    sigGenEtaPhiWrtGenThr_h[pos.at(mI)]->Fill(pSigDataGen.eta_wrtThr[pI], pSigDataGen.phi_wrtThr[pI]);
	    sigGenAbsEtaPhiWrtGenThr_h[pos.at(mI)]->Fill(TMath::Abs(pSigDataGen.eta_wrtThr[pI]), pSigDataGen.phi_wrtThr[pI]);
	  }
	}
      }
          
      for(Int_t pI = 0; pI < pMixData.nParticle; ++pI){
	for(unsigned int mI = 0; mI < pos.size(); ++mI){
	  mixRecoPtWrtRecoThr_h[pos.at(mI)]->Fill(pMixData.pt_wrtThr[pI]);
	  mixRecoEtaWrtRecoThr_h[pos.at(mI)]->Fill(pMixData.eta_wrtThr[pI]);
	  mixRecoAbsEtaWrtRecoThr_h[pos.at(mI)]->Fill(TMath::Abs(pMixData.eta_wrtThr[pI]));
	  mixRecoRapWrtRecoThr_h[pos.at(mI)]->Fill(pMixData.rap_wrtThr[pI]);
	  mixRecoThetaWrtRecoThr_h[pos.at(mI)]->Fill(pMixData.theta_wrtThr[pI]);
	  mixRecoPhiWrtRecoThr_h[pos.at(mI)]->Fill(pMixData.phi_wrtThr[pI]);
	  mixRecoEtaPhiWrtRecoThr_h[pos.at(mI)]->Fill(pMixData.eta_wrtThr[pI], pMixData.phi_wrtThr[pI]);
	  mixRecoAbsEtaPhiWrtRecoThr_h[pos.at(mI)]->Fill(TMath::Abs(pMixData.eta_wrtThr[pI]), pMixData.phi_wrtThr[pI]);
	}
      }

      if(containsTGen){
	for(Int_t pI = 0; pI < pMixDataGen.nParticle; ++pI){
	  for(unsigned int mI = 0; mI < pos.size(); ++mI){
	    mixGenPtWrtGenThr_h[pos.at(mI)]->Fill(pMixDataGen.pt_wrtThr[pI]);
	    mixGenEtaWrtGenThr_h[pos.at(mI)]->Fill(pMixDataGen.eta_wrtThr[pI]);
	    mixGenAbsEtaWrtGenThr_h[pos.at(mI)]->Fill(TMath::Abs(pMixDataGen.eta_wrtThr[pI]));
	    mixGenRapWrtGenThr_h[pos.at(mI)]->Fill(pMixDataGen.rap_wrtThr[pI]);
	    mixGenThetaWrtGenThr_h[pos.at(mI)]->Fill(pMixDataGen.theta_wrtThr[pI]);
	    mixGenPhiWrtGenThr_h[pos.at(mI)]->Fill(pMixDataGen.phi_wrtThr[pI]);
	    mixGenEtaPhiWrtGenThr_h[pos.at(mI)]->Fill(pMixDataGen.eta_wrtThr[pI], pMixDataGen.phi_wrtThr[pI]);
	    mixGenAbsEtaPhiWrtGenThr_h[pos.at(mI)]->Fill(TMath::Abs(pMixDataGen.eta_wrtThr[pI]), pMixDataGen.phi_wrtThr[pI]);
	  }
	}
      }
    }
  
    inFile1_p->Close();
    delete inFile1_p;

    inFile2_p->Close();
    delete inFile2_p;
  }    
  
  outFile_p->cd();

  std::vector<TH1*> hists_p;
  std::vector<TH2*> hists2_p;
  
  for(Int_t mI = 0; mI < nMultBins+1; ++mI){
    hists_p.push_back(sigRecoPtWrtRecoThr_h[mI]);
    hists_p.push_back(sigRecoEtaWrtRecoThr_h[mI]);
    hists_p.push_back(sigRecoAbsEtaWrtRecoThr_h[mI]);
    hists_p.push_back(sigRecoRapWrtRecoThr_h[mI]);
    hists_p.push_back(sigRecoThetaWrtRecoThr_h[mI]);
    hists_p.push_back(sigRecoPhiWrtRecoThr_h[mI]);
    hists2_p.push_back(sigRecoEtaPhiWrtRecoThr_h[mI]);
    hists2_p.push_back(sigRecoAbsEtaPhiWrtRecoThr_h[mI]);
    
    if(containsTGen){
      hists_p.push_back(sigGenPtWrtGenThr_h[mI]);
      hists_p.push_back(sigGenEtaWrtGenThr_h[mI]);
      hists_p.push_back(sigGenAbsEtaWrtGenThr_h[mI]);
      hists_p.push_back(sigGenRapWrtGenThr_h[mI]);
      hists_p.push_back(sigGenThetaWrtGenThr_h[mI]);
      hists_p.push_back(sigGenPhiWrtGenThr_h[mI]);
      hists2_p.push_back(sigGenEtaPhiWrtGenThr_h[mI]);
      hists2_p.push_back(sigGenAbsEtaPhiWrtGenThr_h[mI]);
    }

    hists_p.push_back(mixRecoPtWrtRecoThr_h[mI]);
    hists_p.push_back(mixRecoEtaWrtRecoThr_h[mI]);
    hists_p.push_back(mixRecoAbsEtaWrtRecoThr_h[mI]);
    hists_p.push_back(mixRecoRapWrtRecoThr_h[mI]);
    hists_p.push_back(mixRecoThetaWrtRecoThr_h[mI]);
    hists_p.push_back(mixRecoPhiWrtRecoThr_h[mI]);
    hists2_p.push_back(mixRecoEtaPhiWrtRecoThr_h[mI]);
    hists2_p.push_back(mixRecoAbsEtaPhiWrtRecoThr_h[mI]);

    if(containsTGen){
      hists_p.push_back(mixGenPtWrtGenThr_h[mI]);
      hists_p.push_back(mixGenAbsEtaWrtGenThr_h[mI]);
      hists_p.push_back(mixGenEtaWrtGenThr_h[mI]);
      hists_p.push_back(mixGenRapWrtGenThr_h[mI]);
      hists_p.push_back(mixGenThetaWrtGenThr_h[mI]);
      hists_p.push_back(mixGenPhiWrtGenThr_h[mI]);
      hists2_p.push_back(mixGenEtaPhiWrtGenThr_h[mI]);
      hists2_p.push_back(mixGenAbsEtaPhiWrtGenThr_h[mI]);
    }
  }

  for(unsigned int hI = 0; hI < hists_p.size(); ++hI){
    /*    
    for(Int_t bI = 0; bI < hists_p.at(hI)->GetNbinsX(); ++bI){
      hists_p.at(hI)->SetBinContent(bI+1, hists_p.at(hI)->GetBinContent(bI+1)/hists_p.at(hI)->GetBinWidth(bI+1));
      hists_p.at(hI)->SetBinError(bI+1, hists_p.at(hI)->GetBinError(bI+1)/hists_p.at(hI)->GetBinWidth(bI+1));
    }
    */
    centerTitles(hists_p.at(hI));
    hists_p.at(hI)->Write("", TObject::kOverwrite);
    delete hists_p.at(hI);
  }

  for(unsigned int hI = 0; hI < hists2_p.size(); ++hI){
    centerTitles(hists2_p.at(hI));
    hists2_p.at(hI)->Write("", TObject::kOverwrite);
    delete hists2_p.at(hI);
  }

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/sigMixThrustComp.exe <inFileName1> <inFileName2>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += sigMixThrustComp(argv[1], argv[2]);
  return retVal;
}
