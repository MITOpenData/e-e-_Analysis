//c and c++ dependencies
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

//root dependencies
#include "TFile.h"
#include "TTree.h"
#include "TDatime.h"
#include "TH1F.h"
#include "TMath.h"
#include "TNamed.h"

//include depdencies
#include "include/eventData.h"
#include "include/particleData.h"
#include "include/ep2003_084_Fig1L.h"

int multDist(const std::string inFileName, std::string outFileName = "")
{
  std::vector<std::string> fileList;
  if(inFileName.find(".root") != std::string::npos) fileList.push_back(inFileName);
  else if(inFileName.find(".txt") != std::string::npos){
    std::ifstream file(inFileName.c_str());
    std::string tempFileStr;

    while(std::getline(file, tempFileStr)){
      if(tempFileStr.size() == 0) continue;
      if(tempFileStr.find(".root") == std::string::npos) continue;

      fileList.push_back(tempFileStr);
    }
  }

  if(fileList.size() == 0){
    std::cout << "Warning: inFileName \'" << inFileName << "\' is not a valid input. please check and try again. return 1" << std::endl;
    return 1;
  }

  TDatime* date = new TDatime();
  if(outFileName.size() == 0){
    outFileName = inFileName;
    while(outFileName.find("/") != std::string::npos) outFileName.replace(0,outFileName.find("/")+1,"");
  }
  if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), 5,"");
  else if(outFileName.find(".txt") != std::string::npos) outFileName.replace(outFileName.find(".txt"), 4,"");
  outFileName = outFileName + "_MultDist_" + std::to_string(date->GetDate()) + ".root";

  ep2003_084_Fig1L f1L;
  const Int_t f1LNPts = f1L.getNPoints();
  std::vector<Double_t> f1LBinsV = f1L.getBins();
  Double_t f1LBins[f1LNPts+1];

  for(Int_t i = 0; i < f1LNPts+1; ++i){
    f1LBins[i] = f1LBinsV.at(i);
    std::cout << "i: " << f1LBins[i] << std::endl;

  }

  const Double_t energyLow = 203.;
  const Double_t energyHigh = 209.;
  const Double_t piPlusMass = .1395706;

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1F* energy_p = new TH1F("energy_h", ";Energy;Counts", 50, energyLow, energyHigh);
  TH1F* year_p = new TH1F("year_h", ";year;Counts", 10, 1990.5, 2000.5);
  TH1F* fig1L_h = new TH1F("fig1L_h", ";;", f1LNPts, f1LBins);
  TH1F* fig1L_Corr_h = new TH1F("fig1L_Corr_h", ";;", f1LNPts, f1LBins);
  TH1F* fig1L_Paper_h = new TH1F("fig1L_Paper_h", ";;", f1LNPts, f1LBins);

  for(Int_t i = 0; i < f1LNPts; ++i){
    fig1L_Paper_h->SetBinContent(i+1, f1L.getPointYVal(i));
    fig1L_Paper_h->SetBinError(i+1, 0.);
  }

  Int_t goodEventf1L = 0;
  Int_t badEventf1LEnergy = 0;
  Int_t badEventf1LOther = 0;

  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    TFile* inFile_p = new TFile(fileList.at(fI).c_str(), "READ");
    TTree* inTree_p = (TTree*)inFile_p->Get("t");
    
    eventData eData;
    particleData pData;
    
    inTree_p->SetBranchStatus("*", 0);
    inTree_p->SetBranchStatus("Energy", 1);
    inTree_p->SetBranchStatus("year", 1);
    inTree_p->SetBranchStatus("nParticle", 1);
    inTree_p->SetBranchStatus("px", 1);
    inTree_p->SetBranchStatus("py", 1);
    inTree_p->SetBranchStatus("pz", 1);
    inTree_p->SetBranchStatus("pwflag", 1);
    inTree_p->SetBranchStatus("charge", 1);
    
    inTree_p->SetBranchAddress("Energy", &(pData.Energy));
    inTree_p->SetBranchAddress("year", &(pData.year));
    inTree_p->SetBranchAddress("nParticle", &(pData.nParticle));
    inTree_p->SetBranchAddress("px", pData.px);
    inTree_p->SetBranchAddress("py", pData.py);
    inTree_p->SetBranchAddress("pz", pData.pz);
    inTree_p->SetBranchAddress("pwflag", pData.pwflag);
    inTree_p->SetBranchAddress("charge", pData.charge);
    
    const Int_t nEntries = inTree_p->GetEntries();
    
    for(Int_t entry = 0; entry < nEntries; ++entry){
      if(entry%10000 == 0) std::cout << "Entry: " << entry << "/" << nEntries << std::endl;
      
      inTree_p->GetEntry(entry);

      if(pData.Energy < energyLow || pData.Energy > energyHigh){
	++badEventf1LEnergy;
	continue;
      }

      Int_t nChg = 0.;
      Double_t eChg = 0.;
      
      for(Int_t pI = 0; pI < pData.nParticle; ++pI){
	if(pData.charge[pI] == 0) continue;

	++nChg;
	Double_t mom = TMath::Sqrt(pData.px[pI]*pData.px[pI] + pData.py[pI]*pData.py[pI] + pData.pz[pI]*pData.pz[pI]);
	eChg += TMath::Sqrt(mom*mom + piPlusMass*piPlusMass);
      }

      //For these cuts see here: http://www.sciencedirect.com/science/article/pii/S0370157397000458
      if(eChg < 15. || nChg < 5){
	++badEventf1LOther;
	continue;
      }

      ++goodEventf1L;

      year_p->Fill(pData.year);
      energy_p->Fill(pData.Energy);
      
      for(Int_t pI = 0; pI < pData.nParticle; ++pI){
	if(pData.charge[pI] == 0) continue;
	
	Double_t mom = TMath::Sqrt(pData.px[pI]*pData.px[pI] + pData.py[pI]*pData.py[pI] + pData.pz[pI]*pData.pz[pI]);
	fig1L_h->Fill(2.*mom/pData.Energy);
	fig1L_Corr_h->Fill(2.*mom/pData.Energy);
      }
    }
    
    inFile_p->Close();
    delete inFile_p;
  }

  fig1L_h->Scale(1./(Double_t)goodEventf1L);
  fig1L_Corr_h->Scale(1./(Double_t)goodEventf1L);

  for(Int_t bI = 0; bI < fig1L_h->GetNbinsX(); ++bI){
    Double_t binCent = fig1L_h->GetBinCenter(bI+1);
    Double_t corrVal = f1L.getCorrVal(binCent);

    fig1L_h->SetBinContent(bI+1, fig1L_h->GetBinContent(bI+1)/fig1L_h->GetBinWidth(bI+1));
    fig1L_h->SetBinError(bI+1, fig1L_h->GetBinError(bI+1)/fig1L_h->GetBinWidth(bI+1));

    fig1L_Corr_h->SetBinContent(bI+1, corrVal*fig1L_Corr_h->GetBinContent(bI+1)/fig1L_Corr_h->GetBinWidth(bI+1));
    fig1L_Corr_h->SetBinError(bI+1, corrVal*fig1L_Corr_h->GetBinError(bI+1)/fig1L_Corr_h->GetBinWidth(bI+1));
  }

  outFile_p->cd();
  energy_p->Write("", TObject::kOverwrite);
  year_p->Write("", TObject::kOverwrite);
  fig1L_h->Write("", TObject::kOverwrite);
  fig1L_Corr_h->Write("", TObject::kOverwrite);
  fig1L_Paper_h->Write("", TObject::kOverwrite);

  TH1F* corr = f1L.getHist();
  corr->Write("", TObject::kOverwrite);

  TNamed goodEventf1LName("goodEventf1L", std::to_string(goodEventf1L).c_str());
  TNamed badEventf1LEnergyName("badEventf1LEnergy", std::to_string(badEventf1LEnergy).c_str());
  TNamed badEventf1LOtherName("badEventf1LOther", std::to_string(badEventf1LOther).c_str());

  goodEventf1LName.Write("", TObject::kOverwrite);
  badEventf1LEnergyName.Write("", TObject::kOverwrite);
  badEventf1LOtherName.Write("", TObject::kOverwrite);

  delete energy_p;
  delete year_p;
  delete fig1L_h;
  delete fig1L_Corr_h;
  delete fig1L_Paper_h;

  outFile_p->Close();
  delete outFile_p;

  delete date;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2 && argc != 3){
    std::cout << "Usage ./multDist.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 2) retVal += multDist(argv[1]);
  else if(argc == 3) retVal += multDist(argv[1], argv[2]);
  return retVal;
}
