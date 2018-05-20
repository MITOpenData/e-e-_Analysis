#include <vector>
#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>

#include "TFile.h"
#include "TH1D.h"
#include "TDatime.h"
#include "TMath.h"

#include "DataProcessing/include/particleData.h"
#include "DataProcessing/include/signalMixTableReader.h"

int testSignalMixTableReader(const std::string inSigFileName, const std::string inMixFileName, const bool isMC, const bool isPYTHIA)
{
  const std::string fullPath =  std::getenv("STUDYMULTDIR");
  signalMixTableReader tableReader;
  if(isPYTHIA) tableReader.Init((fullPath + "/DataProcessing/tables/signalOverMixAbsEtaTable_PYTHIA_20180510.txt").c_str());
  else if(isMC) tableReader.Init((fullPath + "/DataProcessing/tables/signalOverMixAbsEtaTable_MC_20180506.txt").c_str());
  else tableReader.Init((fullPath + "/DataProcessing/tables/signalOverMixAbsEtaPhiTable_Data_20180520.txt").c_str());

  tableReader.Print();
  //  return 1;

  std::vector<std::string> fileListSig;
  std::vector<std::string> fileListMix;

  if(!checkFile(inSigFileName)){
    std::cout << "inSigFileName \'" << inSigFileName << "\' is not a valid file. return 1" << std::endl;
    return 1;
  }
  if(!checkFile(inMixFileName)){
    std::cout << "inMixFileName \'" << inMixFileName << "\' is not a valid file. return 1" << std::endl;
    return 1;
  }
  
  if(inSigFileName.find(".root") != std::string::npos) fileListSig.push_back(inSigFileName);
  else if(inSigFileName.find(".txt") != std::string::npos){
    std::ifstream file(inSigFileName.c_str());
    std::string tempStr;
    while(std::getline(file, tempStr)){
      if(tempStr.size() == 0) continue;
      if(tempStr.find(".root") != std::string::npos) fileListSig.push_back(tempStr);
    }
    file.close();
  }
  else{
    std::cout << "inSigFileName \'" << inSigFileName << "\' is not .root or .txt file. return 1" << std::endl;
    return 1;
  }

  if(inMixFileName.find(".root") != std::string::npos) fileListMix.push_back(inMixFileName);
  else if(inMixFileName.find(".txt") != std::string::npos){
    std::ifstream file(inMixFileName.c_str());
    std::string tempStr;
    while(std::getline(file, tempStr)){
      if(tempStr.size() == 0) continue;
      if(tempStr.find(".root") != std::string::npos) fileListMix.push_back(tempStr);
    }
    file.close();
  }
  else{
    std::cout << "inMixFileName \'" << inMixFileName << "\' is not .root or .txt file. return 1" << std::endl;
    return 1;
  }
  
  TDatime* date = new TDatime();

  std::string outFileName = fileListSig.at(0);
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
  TH1D* signalEta_p[nMultBins+1];
  TH1D* mixEta_p[nMultBins+1];
  TH1D* mixEtaWeighted_p[nMultBins+1];
  TH1D* mixEtaWeighted2_p[nMultBins+1];

  TH1D* signalAbsEta_p[nMultBins+1];
  TH1D* mixAbsEta_p[nMultBins+1];
  TH1D* mixAbsEtaWeighted_p[nMultBins+1];
  TH1D* mixAbsEtaWeighted2_p[nMultBins+1];

  for(Int_t mI = 0; mI < nMultBins+1; ++mI){
    std::string multStr = "Mult" + std::to_string(multBinsLow[0]) + "to" + std::to_string(multBinsHi[nMultBins-1]);
    if(mI < nMultBins) multStr = "Mult" + std::to_string(multBinsLow[mI]) + "to" + std::to_string(multBinsHi[mI]);

    signalEta_p[mI] = new TH1D(("signalEta_" + multStr + "_h").c_str(), ";#eta;Counts", nEtaBins, etaBins);
    mixEta_p[mI] = new TH1D(("mixEta_" + multStr + "_h").c_str(), ";#eta;Counts", nEtaBins, etaBins);
    mixEtaWeighted_p[mI] = new TH1D(("mixEtaWeighted_" + multStr + "_h").c_str(), ";#eta;Counts", nEtaBins, etaBins);    
    mixEtaWeighted2_p[mI] = new TH1D(("mixEtaWeighted2_" + multStr + "_h").c_str(), ";#eta;Counts", nEtaBins, etaBins);    

    signalAbsEta_p[mI] = new TH1D(("signalAbsEta_" + multStr + "_h").c_str(), ";|#eta|;Counts", nAbsEtaBins, absEtaBins);
    mixAbsEta_p[mI] = new TH1D(("mixAbsEta_" + multStr + "_h").c_str(), ";|#eta|;Counts", nAbsEtaBins, absEtaBins);
    mixAbsEtaWeighted_p[mI] = new TH1D(("mixAbsEtaWeighted_" + multStr + "_h").c_str(), ";|#eta|;Counts", nAbsEtaBins, absEtaBins);    
    mixAbsEtaWeighted2_p[mI] = new TH1D(("mixAbsEtaWeighted2_" + multStr + "_h").c_str(), ";|#eta|;Counts", nAbsEtaBins, absEtaBins);    
  }

  particleData pDataSig;
  particleData pDataMix;


  for(unsigned int fI = 0; fI < fileListSig.size(); ++fI){

    TFile* inFileSig_p = new TFile(fileListSig.at(fI).c_str(), "READ");
    TTree* inTreeSig_p = (TTree*)inFileSig_p->Get("t");
    pDataSig.SetStatusAndAddressRead(inTreeSig_p, {"nParticle", "eta_wrtThr"});
    
    TFile* inFileMix_p = new TFile(fileListMix.at(fI).c_str(), "READ");
    TTree* inTreeMix_p = (TTree*)inFileMix_p->Get("t");
    pDataMix.SetStatusAndAddressRead(inTreeMix_p, {"nParticle", "eta_wrtThr"});
    
    const Int_t nEntriesSig = inTreeSig_p->GetEntries();
    //  const Double_t weight = 1./(Double_t)nEntriesSig;
    
    std::cout << "Processing \'" << fileListSig.at(fI) << ", " << fileListMix.at(fI) << "\'..." << std::endl;

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
	for(unsigned int pI2 = 0; pI2 < pos.size(); ++pI2){
	  signalEta_p[pos.at(pI2)]->Fill(pDataSig.eta_wrtThr[pI]);
	  
	  signalAbsEta_p[pos.at(pI2)]->Fill(TMath::Abs(pDataSig.eta_wrtThr[pI]));
	}
      }
      
      for(Int_t pI = 0; pI < pDataMix.nParticle; ++pI){
	for(unsigned int pI2 = 0; pI2 < pos.size(); ++pI2){
	  bool weight = 1.;
	  if(!tableReader.is2D()) weight = tableReader.getSigOverMixFactor(pDataSig.nParticle, pDataMix.eta_wrtThr[pI]);
	  else weight = tableReader.getSigOverMixFactor2D(pDataSig.nParticle, pDataMix.eta_wrtThr[pI], pDataMix.phi_wrtThr[pI]);

	  mixEta_p[pos.at(pI2)]->Fill(pDataMix.eta_wrtThr[pI]);
	  mixEtaWeighted_p[pos.at(pI2)]->Fill(pDataMix.eta_wrtThr[pI], weight);
	  
	  mixAbsEta_p[pos.at(pI2)]->Fill(TMath::Abs(pDataMix.eta_wrtThr[pI]));
	  mixAbsEtaWeighted_p[pos.at(pI2)]->Fill(TMath::Abs(pDataMix.eta_wrtThr[pI]), weight);
	}
      }
    }
    
    inFileSig_p->Close();
    delete inFileSig_p;
    
    inFileMix_p->Close();
    delete inFileMix_p;
  }


  outFile_p->cd();

  
  for(Int_t mI = 0; mI < nMultBins+1; ++mI){
    //     signalEta_p[mI]->Scale(1./signalEta_p[mI]->Integral());
    //    mixEta_p[mI]->Scale(1./mixEta_p[mI]->Integral());
    //    mixEtaWeighted_p[mI]->Scale(1./mixEtaWeighted_p[mI]->Integral());
    
    for(Int_t bI = 0; bI < signalEta_p[mI]->GetNbinsX(); ++bI){
      /*
      signalEta_p[mI]->SetBinContent(bI+1, signalEta_p[mI]->GetBinContent(bI+1)/signalEta_p[mI]->GetBinWidth(bI+1));
      signalEta_p[mI]->SetBinError(bI+1, signalEta_p[mI]->GetBinError(bI+1)/signalEta_p[mI]->GetBinWidth(bI+1));

      mixEta_p[mI]->SetBinContent(bI+1, mixEta_p[mI]->GetBinContent(bI+1)/mixEta_p[mI]->GetBinWidth(bI+1));
      mixEta_p[mI]->SetBinError(bI+1, mixEta_p[mI]->GetBinError(bI+1)/mixEta_p[mI]->GetBinWidth(bI+1));

      mixEtaWeighted_p[mI]->SetBinContent(bI+1, mixEtaWeighted_p[mI]->GetBinContent(bI+1)/mixEtaWeighted_p[mI]->GetBinWidth(bI+1));
      mixEtaWeighted_p[mI]->SetBinError(bI+1, mixEtaWeighted_p[mI]->GetBinError(bI+1)/mixEtaWeighted_p[mI]->GetBinWidth(bI+1));
      */
    }

    for(Int_t bI = 0; bI < signalEta_p[mI]->GetNbinsX(); ++bI){
      mixEtaWeighted2_p[mI]->SetBinContent(bI+1, mixEta_p[mI]->GetBinContent(bI+1)*tableReader.getSigOverMixFactor(multBinsLow[mI], mixEtaWeighted2_p[mI]->GetBinCenter(bI+1)));
      mixEtaWeighted2_p[mI]->SetBinError(bI+1, 0);
    }

    signalEta_p[mI]->Write("", TObject::kOverwrite);
    mixEta_p[mI]->Write("", TObject::kOverwrite);
    mixEtaWeighted_p[mI]->Write("", TObject::kOverwrite);
    mixEtaWeighted2_p[mI]->Write("", TObject::kOverwrite);
  
    delete signalEta_p[mI];
    delete mixEta_p[mI];
    delete mixEtaWeighted_p[mI];
    delete mixEtaWeighted2_p[mI];

    //     signalAbsEta_p[mI]->Scale(1./signalAbsEta_p[mI]->Integral());
    //    mixAbsEta_p[mI]->Scale(1./mixAbsEta_p[mI]->Integral());
    //    mixAbsEtaWeighted_p[mI]->Scale(1./mixAbsEtaWeighted_p[mI]->Integral());

 
    for(Int_t bI = 0; bI < signalAbsEta_p[mI]->GetNbinsX(); ++bI){
      /*
      signalAbsEta_p[mI]->SetBinContent(bI+1, signalAbsEta_p[mI]->GetBinContent(bI+1)/signalAbsEta_p[mI]->GetBinWidth(bI+1));
      signalAbsEta_p[mI]->SetBinError(bI+1, signalAbsEta_p[mI]->GetBinError(bI+1)/signalAbsEta_p[mI]->GetBinWidth(bI+1));

      mixAbsEta_p[mI]->SetBinContent(bI+1, mixAbsEta_p[mI]->GetBinContent(bI+1)/mixAbsEta_p[mI]->GetBinWidth(bI+1));
      mixAbsEta_p[mI]->SetBinError(bI+1, mixAbsEta_p[mI]->GetBinError(bI+1)/mixAbsEta_p[mI]->GetBinWidth(bI+1));

      mixAbsEtaWeighted_p[mI]->SetBinContent(bI+1, mixAbsEtaWeighted_p[mI]->GetBinContent(bI+1)/mixAbsEtaWeighted_p[mI]->GetBinWidth(bI+1));
      mixAbsEtaWeighted_p[mI]->SetBinError(bI+1, mixAbsEtaWeighted_p[mI]->GetBinError(bI+1)/mixAbsEtaWeighted_p[mI]->GetBinWidth(bI+1));
      */
    }


    for(Int_t bI = 0; bI < signalAbsEta_p[mI]->GetNbinsX(); ++bI){
      mixAbsEtaWeighted2_p[mI]->SetBinContent(bI+1, mixAbsEta_p[mI]->GetBinContent(bI+1)*tableReader.getSigOverMixFactor(multBinsLow[mI], mixAbsEtaWeighted2_p[mI]->GetBinCenter(bI+1)));
      mixAbsEtaWeighted2_p[mI]->SetBinError(bI+1, 0);
    }
   
    signalAbsEta_p[mI]->Write("", TObject::kOverwrite);
    mixAbsEta_p[mI]->Write("", TObject::kOverwrite);
    mixAbsEtaWeighted_p[mI]->Write("", TObject::kOverwrite);
    mixAbsEtaWeighted2_p[mI]->Write("", TObject::kOverwrite);

  
    delete signalAbsEta_p[mI];
    delete mixAbsEta_p[mI];
    delete mixAbsEtaWeighted_p[mI];
    delete mixAbsEtaWeighted2_p[mI];
  }

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 5){
    std::cout << "Usage: ./bin/testSignalMixTableReader.exe <inSigFileName> <inMixFileName> <isMC> <isPYTHIA>" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += testSignalMixTableReader(argv[1], argv[2], (bool)std::stoi(argv[3]), (bool)std::stoi(argv[4]));
  return retVal;
}
