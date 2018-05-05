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

  TFile* inFile1_p = new TFile(inFileName1.c_str(), "READ");
  std::vector<std::string> ttreeList1 = returnRootFileContentsList(inFile1_p, "TTree", "");
  removeVectorDuplicates(&ttreeList1);

  TFile* inFile2_p = new TFile(inFileName2.c_str(), "READ");
  std::vector<std::string> ttreeList2 = returnRootFileContentsList(inFile2_p, "TTree", "");
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

  TDatime* date = new TDatime();
  TFile* outFile_p = new TFile(("output/sigMixThrustComp_" + std::to_string(date->GetDate()) + ".root").c_str(), "RECREATE");
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

  TH1F* sigRecoPtWrtRecoThr_h = new TH1F("sigRecoPtWrtRecoThr_h", ";Sig. reco. p_{T} w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{dp_{T}}", nPtBins, ptBins);
  TH1F* sigGenPtWrtGenThr_h = new TH1F("sigRecoPtWrtGenThr_h", ";Sig. gen. p_{T} w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{dp_{T}}", nPtBins, ptBins);
  TH1F* mixRecoPtWrtRecoThr_h = new TH1F("mixRecoPtWrtRecoThr_h", ";Mix. reco. p_{T} w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{dp_{T}}", nPtBins, ptBins);
  TH1F* mixGenPtWrtGenThr_h = new TH1F("mixRecoPtWrtGenThr_h", ";Mix. gen. p_{T} w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{dp_{T}}", nPtBins, ptBins);

  TH1F* sigRecoEtaWrtRecoThr_h = new TH1F("sigRecoEtaWrtRecoThr_h", ";Sig. reco. #eta w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#eta}", nEtaBins, etaLow, etaHi);
  TH1F* sigGenEtaWrtGenThr_h = new TH1F("sigRecoEtaWrtGenThr_h", ";Sig. gen. #eta w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#eta}", nEtaBins, etaLow, etaHi);
  TH1F* mixRecoEtaWrtRecoThr_h = new TH1F("mixRecoEtaWrtRecoThr_h", ";Mix. reco. #eta w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#eta}", nEtaBins, etaLow, etaHi);
  TH1F* mixGenEtaWrtGenThr_h = new TH1F("mixRecoEtaWrtGenThr_h", ";Mix. gen. #eta w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#eta}", nEtaBins, etaLow, etaHi);

  TH1F* sigRecoRapWrtRecoThr_h = new TH1F("sigRecoRapWrtRecoThr_h", ";Sig. reco. y w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{dy}", nRapBins, rapLow, rapHi);
  TH1F* sigGenRapWrtGenThr_h = new TH1F("sigRecoRapWrtGenThr_h", ";Sig. gen. y w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{dy}", nRapBins, rapLow, rapHi);
  TH1F* mixRecoRapWrtRecoThr_h = new TH1F("mixRecoRapWrtRecoThr_h", ";Mix. reco. y w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{dy}", nRapBins, rapLow, rapHi);
  TH1F* mixGenRapWrtGenThr_h = new TH1F("mixRecoRapWrtGenThr_h", ";Mix. gen. y w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{dy}", nRapBins, rapLow, rapHi);

  TH1F* sigRecoThetaWrtRecoThr_h = new TH1F("sigRecoThetaWrtRecoThr_h", ";Sig. reco. #theta w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#theta}", nThetaBins, thetaLow, thetaHi);
  TH1F* sigGenThetaWrtGenThr_h = new TH1F("sigRecoThetaWrtGenThr_h", ";Sig. gen. #theta w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#theta}", nThetaBins, thetaLow, thetaHi);
  TH1F* mixRecoThetaWrtRecoThr_h = new TH1F("mixRecoThetaWrtRecoThr_h", ";Mix. reco. #theta w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#theta}", nThetaBins, thetaLow, thetaHi);
  TH1F* mixGenThetaWrtGenThr_h = new TH1F("mixRecoThetaWrtGenThr_h", ";Mix. gen. #theta w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#theta}", nThetaBins, thetaLow, thetaHi);

  TH1F* sigRecoPhiWrtRecoThr_h = new TH1F("sigRecoPhiWrtRecoThr_h", ";Sig. reco. #phi w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#phi}", nPhiBins, phiLow, phiHi);
  TH1F* sigGenPhiWrtGenThr_h = new TH1F("sigRecoPhiWrtGenThr_h", ";Sig. gen. #phi w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#phi}", nPhiBins, phiLow, phiHi);
  TH1F* mixRecoPhiWrtRecoThr_h = new TH1F("mixRecoPhiWrtRecoThr_h", ";Mix. reco. #phi w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#phi}", nPhiBins, phiLow, phiHi);
  TH1F* mixGenPhiWrtGenThr_h = new TH1F("mixRecoPhiWrtGenThr_h", ";Mix. gen. #phi w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{d#phi}", nPhiBins, phiLow, phiHi);

  particleData pSigData;
  particleData pSigDataGen;
  particleData pMixData;
  particleData pMixDataGen;

  TTree* tsig = (TTree*)inFile1_p->Get("t");
  TTree* tsiggen = (TTree*)inFile1_p->Get("tgen");

  TTree* tmix = (TTree*)inFile2_p->Get("t");
  TTree* tmixgen = (TTree*)inFile2_p->Get("tgen");

  pSigData.SetStatusAndAddressRead(tsig, {});
  pSigDataGen.SetStatusAndAddressRead(tsiggen, {});

  pMixData.SetStatusAndAddressRead(tmix, {});
  pMixDataGen.SetStatusAndAddressRead(tmixgen, {});
  
  const Int_t nEntries = tsig->GetEntries();
  const Double_t weight = 1./(Double_t)nEntries;

  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%10000 == 0) std::cout << "Entry: " << entry << "/" << nEntries << std::endl;
    //    if(entry == 100000) break;

    tsig->GetEntry(entry);
    tsiggen->GetEntry(entry);

    tmix->GetEntry(entry);
    tmixgen->GetEntry(entry);

    for(Int_t pI = 0; pI < pSigData.nParticle; ++pI){
      sigRecoPtWrtRecoThr_h->Fill(pSigData.pt_wrtThr[pI], weight);
      sigRecoEtaWrtRecoThr_h->Fill(pSigData.eta_wrtThr[pI], weight);
      sigRecoRapWrtRecoThr_h->Fill(pSigData.rap_wrtThr[pI], weight);
      sigRecoThetaWrtRecoThr_h->Fill(pSigData.theta_wrtThr[pI], weight);
      sigRecoPhiWrtRecoThr_h->Fill(pSigData.phi_wrtThr[pI], weight);
    }

    for(Int_t pI = 0; pI < pSigDataGen.nParticle; ++pI){
      sigGenPtWrtGenThr_h->Fill(pSigDataGen.pt_wrtThr[pI], weight);
      sigGenEtaWrtGenThr_h->Fill(pSigDataGen.eta_wrtThr[pI], weight);
      sigGenRapWrtGenThr_h->Fill(pSigDataGen.rap_wrtThr[pI], weight);
      sigGenThetaWrtGenThr_h->Fill(pSigDataGen.theta_wrtThr[pI], weight);
      sigGenPhiWrtGenThr_h->Fill(pSigDataGen.phi_wrtThr[pI], weight);
    }

    for(Int_t pI = 0; pI < pMixData.nParticle; ++pI){
      mixRecoPtWrtRecoThr_h->Fill(pMixData.pt_wrtThr[pI], weight);
      mixRecoEtaWrtRecoThr_h->Fill(pMixData.eta_wrtThr[pI], weight);
      mixRecoRapWrtRecoThr_h->Fill(pMixData.rap_wrtThr[pI], weight);
      mixRecoThetaWrtRecoThr_h->Fill(pMixData.theta_wrtThr[pI], weight);
      mixRecoPhiWrtRecoThr_h->Fill(pMixData.phi_wrtThr[pI], weight);
    }

    for(Int_t pI = 0; pI < pMixDataGen.nParticle; ++pI){
      mixGenPtWrtGenThr_h->Fill(pMixDataGen.pt_wrtThr[pI], weight);
      mixGenEtaWrtGenThr_h->Fill(pMixDataGen.eta_wrtThr[pI], weight);
      mixGenRapWrtGenThr_h->Fill(pMixDataGen.rap_wrtThr[pI], weight);
      mixGenThetaWrtGenThr_h->Fill(pMixDataGen.theta_wrtThr[pI], weight);
      mixGenPhiWrtGenThr_h->Fill(pMixDataGen.phi_wrtThr[pI], weight);
    }
  }

  outFile_p->cd();

  std::vector<TH1*> hists_p;
  hists_p.push_back(sigRecoPtWrtRecoThr_h);
  hists_p.push_back(sigRecoEtaWrtRecoThr_h);
  hists_p.push_back(sigRecoRapWrtRecoThr_h);
  hists_p.push_back(sigRecoThetaWrtRecoThr_h);
  hists_p.push_back(sigRecoPhiWrtRecoThr_h);

  hists_p.push_back(sigGenPtWrtGenThr_h);
  hists_p.push_back(sigGenEtaWrtGenThr_h);
  hists_p.push_back(sigGenRapWrtGenThr_h);
  hists_p.push_back(sigGenThetaWrtGenThr_h);
  hists_p.push_back(sigGenPhiWrtGenThr_h);

  hists_p.push_back(mixRecoPtWrtRecoThr_h);
  hists_p.push_back(mixRecoEtaWrtRecoThr_h);
  hists_p.push_back(mixRecoRapWrtRecoThr_h);
  hists_p.push_back(mixRecoThetaWrtRecoThr_h);
  hists_p.push_back(mixRecoPhiWrtRecoThr_h);

  hists_p.push_back(mixGenPtWrtGenThr_h);
  hists_p.push_back(mixGenEtaWrtGenThr_h);
  hists_p.push_back(mixGenRapWrtGenThr_h);
  hists_p.push_back(mixGenThetaWrtGenThr_h);
  hists_p.push_back(mixGenPhiWrtGenThr_h);

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

  inFile1_p->Close();
  delete inFile1_p;

  inFile2_p->Close();
  delete inFile2_p;

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
