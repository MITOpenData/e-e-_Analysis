#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TDatime.h"
#include "TVector3.h"

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

  TH1F* recoPtWrtRecoThr_h = new TH1F("recoPtWrtRecoThr_h", ";p_{T} w.r.t Reco. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{dp_{T}}", nPtBins, ptBins);
  TH1F* recoPtWrtGenThr_h = new TH1F("recoPtWrtGenThr_h", ";p_{T} w.r.t Gen. Thrust;#frac{1}{N_{evt}} #frac{dN_{particle}}{dp_{T}}", nPtBins, ptBins);

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
    }

    TVector3 thrust = getThrust(pDataGen.nParticle, pDataGen.px, pDataGen.py, pDataGen.pz, THRUST::OPTIMAL);

    setThrustVariables(&pData, &eData, thrust, thrust);

    for(Int_t pI = 0; pI < pData.nParticle; ++pI){
      recoPtWrtGenThr_h->Fill(pData.pt_wrtThr[pI], weight);
    }    
  }

  outFile_p->cd();

  for(Int_t bI = 0; bI < recoPtWrtRecoThr_h->GetNbinsX(); ++bI){
    recoPtWrtRecoThr_h->SetBinContent(bI+1, recoPtWrtRecoThr_h->GetBinContent(bI+1)/recoPtWrtRecoThr_h->GetBinWidth(bI+1));
    recoPtWrtRecoThr_h->SetBinError(bI+1, recoPtWrtRecoThr_h->GetBinError(bI+1)/recoPtWrtRecoThr_h->GetBinWidth(bI+1));

    recoPtWrtGenThr_h->SetBinContent(bI+1, recoPtWrtGenThr_h->GetBinContent(bI+1)/recoPtWrtGenThr_h->GetBinWidth(bI+1));
    recoPtWrtGenThr_h->SetBinError(bI+1, recoPtWrtGenThr_h->GetBinError(bI+1)/recoPtWrtGenThr_h->GetBinWidth(bI+1));
  }

  recoPtWrtRecoThr_h->Write("", TObject::kOverwrite);
  recoPtWrtGenThr_h->Write("", TObject::kOverwrite);

  delete recoPtWrtRecoThr_h;
  delete recoPtWrtGenThr_h;

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
