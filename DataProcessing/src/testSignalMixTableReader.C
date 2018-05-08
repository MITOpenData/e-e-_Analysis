#include <iostream>
#include <cstdlib>
#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TDatime.h"

#include "DataProcessing/include/particleData.h"
#include "DataProcessing/include/signalMixTableReader.h"

int testSignalMixTableReader(const std::string inSigFileName, const std::string inMixFileName, const bool isMC)
{
  const std::string fullPath =  std::getenv("STUDYMULTDIR");
  signalMixTableReader tableReader;
  if(isMC) tableReader.Init((fullPath + "/DataProcessing/tables/signalOverMixAbsEtaTable_MC_20180506.txt").c_str());
  else tableReader.Init((fullPath + "/DataProcessing/tables/signalOverMixAbsEtaTable_Data_20180508.txt").c_str());

  tableReader.Print();
  //  return 1;

  TDatime* date = new TDatime();

  std::string outFileName = inSigFileName;
  while(outFileName.find("/") != std::string::npos){
    outFileName.replace(0, outFileName.find("/")+1, "");
  }
  if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), 5, "");
  outFileName = "output/" + outFileName + "_SignalMixCheck_" + std::to_string(date->GetDate()) + ".root";

  const Int_t nEtaBins = 50;
  Float_t etaBins[nEtaBins+1] = {-10., -6.000, -5.5, -5.0, -4.5, -4.0, -3.800, -3.600, -3.400, -3.200, -3.000, -2.800, -2.600, -2.400, -2.200, -2.000, -1.800, -1.600, -1.400, -1.200, -1.000, -0.800, -0.600, -0.400, -0.200, 0.000, 0.200, 0.400, 0.600, 0.800, 1.000, 1.200, 1.400, 1.600, 1.800, 2.000, 2.200, 2.400, 2.600, 2.800, 3.000, 3.200, 3.400, 3.600, 3.800, 4.0, 4.5, 5.0, 5.5, 6.000, 10.};

  const Int_t nAbsEtaBins = 25;
  Float_t absEtaBins[nAbsEtaBins+1] = {0.000, 0.200, 0.400, 0.600, 0.800, 1.000, 1.200, 1.400, 1.600, 1.800, 2.000, 2.200, 2.400, 2.600, 2.800, 3.000, 3.200, 3.400, 3.600, 3.800, 4.0, 4.5, 5.0, 5.5, 6.000, 10.};

  const Int_t nMultBins = 3;
  const Int_t multBinsLow[nMultBins] = {0, 20, 30};
  const Int_t multBinsHi[nMultBins] = {20, 30, 200};


  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1F* signalEta_p[nMultBins+1];
  TH1F* mixEta_p[nMultBins+1];
  TH1F* mixEtaWeighted_p[nMultBins+1];

  TH1F* signalAbsEta_p[nMultBins+1];
  TH1F* mixAbsEta_p[nMultBins+1];
  TH1F* mixAbsEtaWeighted_p[nMultBins+1];

  for(Int_t mI = 0; mI < nMultBins+1; ++mI){
    std::string multStr = "Mult" + std::to_string(multBinsLow[0]) + "to" + std::to_string(multBinsHi[nMultBins-1]);
    if(mI < nMultBins) multStr = "Mult" + std::to_string(multBinsLow[mI]) + "to" + std::to_string(multBinsHi[mI]);

    signalEta_p[mI] = new TH1F(("signalEta_" + multStr + "_h").c_str(), ";#eta;Counts", nEtaBins, etaBins);
    mixEta_p[mI] = new TH1F(("mixEta_" + multStr + "_h").c_str(), ";#eta;Counts", nEtaBins, etaBins);
    mixEtaWeighted_p[mI] = new TH1F(("mixEtaWeighted_" + multStr + "_h").c_str(), ";#eta;Counts", nEtaBins, etaBins);    

    signalAbsEta_p[mI] = new TH1F(("signalAbsEta_" + multStr + "_h").c_str(), ";|#eta|;Counts", nAbsEtaBins, absEtaBins);
    mixAbsEta_p[mI] = new TH1F(("mixAbsEta_" + multStr + "_h").c_str(), ";|#eta|;Counts", nAbsEtaBins, absEtaBins);
    mixAbsEtaWeighted_p[mI] = new TH1F(("mixAbsEtaWeighted_" + multStr + "_h").c_str(), ";|#eta|;Counts", nAbsEtaBins, absEtaBins);    
  }

  particleData pDataSig;
  particleData pDataMix;

  TFile* inFileSig_p = new TFile(inSigFileName.c_str(), "READ");
  TTree* inTreeSig_p = (TTree*)inFileSig_p->Get("t");
  pDataSig.SetStatusAndAddressRead(inTreeSig_p, {"nParticle", "eta_wrtThr"});

  TFile* inFileMix_p = new TFile(inMixFileName.c_str(), "READ");
  TTree* inTreeMix_p = (TTree*)inFileMix_p->Get("t");
  pDataMix.SetStatusAndAddressRead(inTreeMix_p, {"nParticle", "eta_wrtThr"});

  const Int_t nEntriesSig = inTreeSig_p->GetEntries();

  std::cout << "Processing \'" << inSigFileName << "\'..." << std::endl;
  for(Int_t entry = 0; entry < nEntriesSig; ++entry){
    if(entry%10000 == 0) std::cout << " Entry: " << entry << "/" << nEntriesSig << std::endl;
    inTreeSig_p->GetEntry(entry);
    inTreeMix_p->GetEntry(entry);

    std::vector<int> pos;
    pos.push_back(nMultBins);

    for(Int_t pI = 0; pI < nMultBins; ++pI){
      if(multBinsLow[pI] <= pDataSig.nParticle && multBinsHi[pI] > pDataSig.nParticle){
	pos.push_back(pI);
	break;
      }
    }

    for(Int_t pI = 0; pI < pDataSig.nParticle; ++pI){
      Double_t weight = 1./(Double_t)nEntriesSig;

      for(unsigned int pI2 = 0; pI2 < pos.size(); ++pI2){
	signalEta_p[pos.at(pI2)]->Fill(pDataSig.eta_wrtThr[pI], weight);
	
	signalAbsEta_p[pos.at(pI2)]->Fill(TMath::Abs(pDataSig.eta_wrtThr[pI]), weight);
      }
    }

    for(Int_t pI = 0; pI < pDataMix.nParticle; ++pI){
      Double_t weight = 1./(Double_t)nEntriesSig;
      Double_t weight2 = tableReader.getSigOverMixFactor(pDataSig.nParticle, pDataMix.eta_wrtThr[pI])/(Double_t)nEntriesSig;

      for(unsigned int pI2 = 0; pI2 < pos.size(); ++pI2){
	mixEta_p[pos.at(pI2)]->Fill(pDataMix.eta_wrtThr[pI], weight);
	mixEtaWeighted_p[pos.at(pI2)]->Fill(pDataMix.eta_wrtThr[pI], weight2);
	
	mixAbsEta_p[pos.at(pI2)]->Fill(TMath::Abs(pDataMix.eta_wrtThr[pI]), weight);
	mixAbsEtaWeighted_p[pos.at(pI2)]->Fill(TMath::Abs(pDataMix.eta_wrtThr[pI]), weight2);
      }
    }
  }

  inFileSig_p->Close();
  delete inFileSig_p;

  inFileMix_p->Close();
  delete inFileMix_p;

  outFile_p->cd();

  
  for(Int_t mI = 0; mI < nMultBins+1; ++mI){
    //     signalEta_p[mI]->Scale(1./signalEta_p[mI]->Integral());
    //    mixEta_p[mI]->Scale(1./mixEta_p[mI]->Integral());
    //    mixEtaWeighted_p[mI]->Scale(1./mixEtaWeighted_p[mI]->Integral());
    
    for(Int_t bI = 0; bI < signalEta_p[mI]->GetNbinsX(); ++bI){
      signalEta_p[mI]->SetBinContent(bI+1, signalEta_p[mI]->GetBinContent(bI+1)/signalEta_p[mI]->GetBinWidth(bI+1));
      signalEta_p[mI]->SetBinError(bI+1, signalEta_p[mI]->GetBinError(bI+1)/signalEta_p[mI]->GetBinWidth(bI+1));

      mixEta_p[mI]->SetBinContent(bI+1, mixEta_p[mI]->GetBinContent(bI+1)/mixEta_p[mI]->GetBinWidth(bI+1));
      mixEta_p[mI]->SetBinError(bI+1, mixEta_p[mI]->GetBinError(bI+1)/mixEta_p[mI]->GetBinWidth(bI+1));

      mixEtaWeighted_p[mI]->SetBinContent(bI+1, mixEtaWeighted_p[mI]->GetBinContent(bI+1)/mixEtaWeighted_p[mI]->GetBinWidth(bI+1));
      mixEtaWeighted_p[mI]->SetBinError(bI+1, mixEtaWeighted_p[mI]->GetBinError(bI+1)/mixEtaWeighted_p[mI]->GetBinWidth(bI+1));
    }

    signalEta_p[mI]->Write("", TObject::kOverwrite);
    mixEta_p[mI]->Write("", TObject::kOverwrite);
    mixEtaWeighted_p[mI]->Write("", TObject::kOverwrite);
  
    delete signalEta_p[mI];
    delete mixEta_p[mI];
    delete mixEtaWeighted_p[mI];

    //     signalAbsEta_p[mI]->Scale(1./signalAbsEta_p[mI]->Integral());
    //    mixAbsEta_p[mI]->Scale(1./mixAbsEta_p[mI]->Integral());
    //    mixAbsEtaWeighted_p[mI]->Scale(1./mixAbsEtaWeighted_p[mI]->Integral());

 
    for(Int_t bI = 0; bI < signalAbsEta_p[mI]->GetNbinsX(); ++bI){
      signalAbsEta_p[mI]->SetBinContent(bI+1, signalAbsEta_p[mI]->GetBinContent(bI+1)/signalAbsEta_p[mI]->GetBinWidth(bI+1));
      signalAbsEta_p[mI]->SetBinError(bI+1, signalAbsEta_p[mI]->GetBinError(bI+1)/signalAbsEta_p[mI]->GetBinWidth(bI+1));

      mixAbsEta_p[mI]->SetBinContent(bI+1, mixAbsEta_p[mI]->GetBinContent(bI+1)/mixAbsEta_p[mI]->GetBinWidth(bI+1));
      mixAbsEta_p[mI]->SetBinError(bI+1, mixAbsEta_p[mI]->GetBinError(bI+1)/mixAbsEta_p[mI]->GetBinWidth(bI+1));

      mixAbsEtaWeighted_p[mI]->SetBinContent(bI+1, mixAbsEtaWeighted_p[mI]->GetBinContent(bI+1)/mixAbsEtaWeighted_p[mI]->GetBinWidth(bI+1));
      mixAbsEtaWeighted_p[mI]->SetBinError(bI+1, mixAbsEtaWeighted_p[mI]->GetBinError(bI+1)/mixAbsEtaWeighted_p[mI]->GetBinWidth(bI+1));
    }
   
    signalAbsEta_p[mI]->Write("", TObject::kOverwrite);
    mixAbsEta_p[mI]->Write("", TObject::kOverwrite);
    mixAbsEtaWeighted_p[mI]->Write("", TObject::kOverwrite);

  
    delete signalAbsEta_p[mI];
    delete mixAbsEta_p[mI];
    delete mixAbsEtaWeighted_p[mI];
  }

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 4){
    std::cout << "Usage: ./bin/testSignalMixTableReader.exe <inSigFileName> <inMixFileName> <isMC>" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += testSignalMixTableReader(argv[1], argv[2], (bool)std::stoi(argv[3]));
  return retVal;
}
