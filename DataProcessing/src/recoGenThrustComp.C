#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TDatime.h"
#include "TVector3.h"
#include "TMath.h"

#include "DataProcessing/include/checkMakeDir.h"
#include "DataProcessing/include/returnRootFileContentsList.h"
#include "DataProcessing/include/particleData.h"
#include "DataProcessing/include/eventData.h"
#include "DataProcessing/include/removeVectorDuplicates.h"
#include "DataProcessing/include/thrustTools.h"
#include "TwoParticleCorrelation/include/getLogBins.h"

int recoGenThrustComp(const std::string inFileName)
{
  if(!checkFile(inFileName)){
    std::cout << "Given inFileName \'" << inFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  std::vector<std::string> ttreeList = returnRootFileContentsList(inFile_p, "TTree", "");
  removeVectorDuplicates(&ttreeList);

  bool containsTGen = false;
  for(unsigned int tI = 0; tI < ttreeList.size(); ++tI){
    //    std::cout << "TTree " << tI << "/" << ttreeList.size() << ": " << ttreeList.at(tI) << std::endl;
    if(ttreeList.at(tI).find("tgen") != std::string::npos){
      containsTGen = true;
      break;
    }
  }

  if(!containsTGen){
    std::cout << "inFile \'" << inFileName << "\' does not contain tgen TTree. return 1" << std::endl;

    inFile_p->Close();
    delete inFile_p;
  }

  TDatime* date = new TDatime();
  TFile* outFile_p = new TFile(("output/recoGenThrustComp_" + std::to_string(date->GetDate()) + ".root").c_str(), "RECREATE");
  delete date;

  const Int_t nPtBins = 20;
  Float_t ptLow = 0.2;
  Float_t ptHi = 30.;
  Double_t ptBins[nPtBins+1];
  getLogBins(ptLow, ptHi, nPtBins, ptBins);

  const Int_t nEtaBins = 20;
  Float_t etaLow = -2.;
  Float_t etaHi = 2.;

  const Int_t nRapBins = 20;
  Float_t rapLow = -2.;
  Float_t rapHi = 2.;

  const Int_t nThetaBins = 20;
  Float_t thetaLow = 0;
  Float_t thetaHi = TMath::Pi();

  const Int_t nPhiBins = 20;
  Float_t phiLow = -TMath::Pi();
  Float_t phiHi = TMath::Pi();

  TH1F* recoPtWrtRecoThr_h = new TH1F("recoPtWrtRecoThr_h", ";p_{T} w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{dp_{T}}", nPtBins, ptBins);
  TH1F* recoPtWrtGenThr_h = new TH1F("recoPtWrtGenThr_h", ";p_{T} w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{dp_{T}}", nPtBins, ptBins);
  TH1F* recoEtaWrtRecoThr_h = new TH1F("recoEtaWrtRecoThr_h", ";#eta w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#eta}", nEtaBins, etaLow, etaHi);
  TH1F* recoEtaWrtGenThr_h = new TH1F("recoEtaWrtGenThr_h", ";#eta w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#eta}", nEtaBins, etaLow, etaHi);
  TH1F* recoRapWrtRecoThr_h = new TH1F("recoRapWrtRecoThr_h", ";y w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{dy}", nRapBins, rapLow, rapHi);
  TH1F* recoRapWrtGenThr_h = new TH1F("recoRapWrtGenThr_h", ";y w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{dy}", nRapBins, rapLow, rapHi);
  TH1F* recoThetaWrtRecoThr_h = new TH1F("recoThetaWrtRecoThr_h", ";#theta w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#theta}", nThetaBins, thetaLow, thetaHi);
  TH1F* recoThetaWrtGenThr_h = new TH1F("recoThetaWrtGenThr_h", ";#theta w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#theta}", nThetaBins, thetaLow, thetaHi);
  TH1F* recoPhiWrtRecoThr_h = new TH1F("recoPhiWrtRecoThr_h", ";#phi w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#phi}", nPhiBins, phiLow, phiHi);
  TH1F* recoPhiWrtGenThr_h = new TH1F("recoPhiWrtGenThr_h", ";#phi w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#phi}", nPhiBins, phiLow, phiHi);

  inFile_p->cd();
  particleData pData;
  particleData pDataGen;
  eventData eData;
  eventData eDataGen;

  TTree* t = (TTree*)inFile_p->Get("t");
  TTree* tgen = (TTree*)inFile_p->Get("tgen");

  pData.SetStatusAndAddressRead(t, {});
  pDataGen.SetStatusAndAddressRead(tgen, {});

  eData.SetStatusAndAddressRead(t, {});
  eDataGen.SetStatusAndAddressRead(tgen, {});
  
  const Int_t nEntries = t->GetEntries();
  const Double_t weight = 1./(Double_t)nEntries;

  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%10000 == 0) std::cout << "Entry: " << entry << "/" << nEntries << std::endl;
    //    if(entry == 100000) break;

    t->GetEntry(entry);
    tgen->GetEntry(entry);

    for(Int_t pI = 0; pI < pData.nParticle; ++pI){
      recoPtWrtRecoThr_h->Fill(pData.pt_wrtThr[pI], weight);
      recoEtaWrtRecoThr_h->Fill(pData.eta_wrtThr[pI], weight);
      recoRapWrtRecoThr_h->Fill(pData.rap_wrtThr[pI], weight);
      recoThetaWrtRecoThr_h->Fill(pData.theta_wrtThr[pI], weight);
      recoPhiWrtRecoThr_h->Fill(pData.phi_wrtThr[pI], weight);
    }

    TVector3 thrust = getThrust(pDataGen.nParticle, pDataGen.px, pDataGen.py, pDataGen.pz, THRUST::OPTIMAL);

    setThrustVariables(&pData, &eData, thrust, thrust);

    for(Int_t pI = 0; pI < pData.nParticle; ++pI){
      recoPtWrtGenThr_h->Fill(pData.pt_wrtThr[pI], weight);
      recoEtaWrtGenThr_h->Fill(pData.eta_wrtThr[pI], weight);
      recoRapWrtGenThr_h->Fill(pData.rap_wrtThr[pI], weight);
      recoThetaWrtGenThr_h->Fill(pData.theta_wrtThr[pI], weight);
      recoPhiWrtGenThr_h->Fill(pData.phi_wrtThr[pI], weight);
    }    
  }

  outFile_p->cd();

  std::vector<TH1*> hists_p;
  hists_p.push_back(recoPtWrtRecoThr_h);
  hists_p.push_back(recoEtaWrtRecoThr_h);
  hists_p.push_back(recoRapWrtRecoThr_h);
  hists_p.push_back(recoThetaWrtRecoThr_h);
  hists_p.push_back(recoPhiWrtRecoThr_h);

  hists_p.push_back(recoPtWrtGenThr_h);
  hists_p.push_back(recoEtaWrtGenThr_h);
  hists_p.push_back(recoRapWrtGenThr_h);
  hists_p.push_back(recoThetaWrtGenThr_h);
  hists_p.push_back(recoPhiWrtGenThr_h);

  for(unsigned int hI = 0; hI < hists_p.size(); ++hI){
    for(Int_t bI = 0; bI < hists_p.at(hI)->GetNbinsX(); ++bI){
      hists_p.at(hI)->SetBinContent(bI+1, hists_p.at(hI)->GetBinContent(bI+1)/hists_p.at(hI)->GetBinWidth(bI+1));
      hists_p.at(hI)->SetBinError(bI+1, hists_p.at(hI)->GetBinError(bI+1)/hists_p.at(hI)->GetBinWidth(bI+1));
    }

    hists_p.at(hI)->Write("", TObject::kOverwrite);
    delete hists_p.at(hI);
  }

  outFile_p->Close();
  delete outFile_p;

  inFile_p->Close();
  delete inFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/recoGenThrustComp.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += recoGenThrustComp(argv[1]);
  return retVal;
}
