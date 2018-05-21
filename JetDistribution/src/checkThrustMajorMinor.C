//Checking that thrust calculations give reasonable answers (dot products are all zero, magnitudes ordered, etc.)
//Original author Chris McGinn

//cpp dependencies
#include <iostream>
#include <string>
#include <iostream>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TDatime.h"
#include "TVector3.h"

//Non-Local StudyMult (DataProcessing, TwoParticleCorrelation, etc.) dependencies
#include "DataProcessing/include/checkMakeDir.h"
#include "DataProcessing/include/histDefUtility.h"
#include "TwoParticleCorrelation/include/getLinBins.h"
#include "TwoParticleCorrelation/include/getLogBins.h"

//Local StudyMult (JetDistribution) dependencies
#include "JetDistribution/include/globalJetVar.h"

int checkThrustMajorMinor(const std::string inFileName)
{
  if(!checkFile(inFileName)){
    std::cout << "Given inFileName \'" << inFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }
  else if(inFileName.find(".root") == std::string::npos){
    std::cout << "Given inFileName \'" << inFileName << "\' is not a root file. return 1" << std::endl;
    return 1;
  }

  checkMakeDir("output");

  TDatime* date = new TDatime();
  std::string outFileName = inFileName;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  outFileName.replace(outFileName.find(".root"), 5, "");
  outFileName = "output/" + outFileName + "_CheckThrustMajorMinor_" + std::to_string(date->GetDate()) + ".root";
  delete date;

  const Int_t nBins = 80;
  const Float_t thrustLow = 0.;
  const Float_t thrustHi = 1.;
  Double_t bins[nBins];
  getLinBins(thrustLow, thrustHi, nBins, bins);

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1F* thrust_h = new TH1F("thrust_h", ";Thrust;Counts", nBins, bins);
  TH1F* thrustMajor_h = new TH1F("thrustMajor_h", ";Thrust (Major);Counts", nBins, bins);
  TH1F* thrustMinor_h = new TH1F("thrustMinor_h", ";Thrust (Minor);Counts", nBins, bins);
  TH1F* thrustDotMajor_h = new TH1F("thrustDotMajor_h", ";T#circT_{MA};Counts", 100, 0, .1);
  TH1F* thrustDotMinor_h = new TH1F("thrustDotMinor_h", ";T#circT_{MI};Counts", 100, 0, .1);
  TH1F* thrustMajorDotMinor_h = new TH1F("thrustMajorDotMinor_h", ";T_{MA}#circT_{MI};Counts", 100, 0, .1);
  TH1F* thrustMagMinusThrustMajorMag_h = new TH1F("thrustMagMinusThrustMajorMag_h", ";T - T_{MA};Counts", 200, -1, 1);
  TH1F* thrustMajorMagMinusThrustMinorMag_h = new TH1F("thrustMajorMagMinusThrustMinorMag_h", ";T_{MA} - T_{MI};Counts", 200, -1, 1);

  centerTitles({thrust_h, thrustMajor_h, thrustMinor_h, thrustDotMajor_h, thrustDotMinor_h, thrustMajorDotMinor_h, thrustMagMinusThrustMajorMag_h, thrustMajorMagMinusThrustMinorMag_h});
  setSumW2({thrust_h, thrustMajor_h, thrustMinor_h, thrustDotMajor_h, thrustDotMinor_h, thrustMajorDotMinor_h, thrustMagMinusThrustMajorMag_h, thrustMajorMagMinusThrustMinorMag_h});

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* inTree_p = (TTree*)inFile_p->Get("generalJetTree_Global");
  globalJetVar jetVar;
  jetVar.SetStatusAndAddressRead(inTree_p, {});

  const Int_t nEntries = inTree_p->GetEntries();

  std::cout << "Processing " << nEntries << " events..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%10000 == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;
    inTree_p->GetEntry(entry);

    thrust_h->Fill(jetVar.thrustMag);
    thrustMajor_h->Fill(jetVar.thrustMajorMag);
    thrustMinor_h->Fill(jetVar.thrustMinorMag);

    TVector3 thrust(jetVar.thrustPx, jetVar.thrustPy, jetVar.thrustPz);
    TVector3 thrustMajor(jetVar.thrustMajorPx, jetVar.thrustMajorPy, jetVar.thrustMajorPz);
    TVector3 thrustMinor(jetVar.thrustMinorPx, jetVar.thrustMinorPy, jetVar.thrustMinorPz);
    
    double thrustDotMajor = thrust.Dot(thrustMajor);
    double thrustDotMinor = thrust.Dot(thrustMinor);
    double thrustMajorDotMinor = thrustMajor.Dot(thrustMinor);

    if(thrustDotMajor > thrustDotMajor_h->GetBinLowEdge(thrustDotMajor_h->GetNbinsX()+1)) thrustDotMajor = thrustDotMajor_h->GetBinCenter(thrustDotMajor_h->GetNbinsX());
    if(thrustDotMinor > thrustDotMinor_h->GetBinLowEdge(thrustDotMinor_h->GetNbinsX()+1)) thrustDotMinor = thrustDotMinor_h->GetBinCenter(thrustDotMinor_h->GetNbinsX());
    if(thrustMajorDotMinor > thrustMajorDotMinor_h->GetBinLowEdge(thrustMajorDotMinor_h->GetNbinsX()+1)) thrustMajorDotMinor = thrustMajorDotMinor_h->GetBinCenter(thrustMajorDotMinor_h->GetNbinsX());

    thrustDotMajor_h->Fill(thrustDotMajor);
    thrustDotMinor_h->Fill(thrustDotMinor);
    thrustMajorDotMinor_h->Fill(thrustMajorDotMinor);

    thrustMagMinusThrustMajorMag_h->Fill(jetVar.thrustMag - jetVar.thrustMajorMag);
    thrustMajorMagMinusThrustMinorMag_h->Fill(jetVar.thrustMajorMag - jetVar.thrustMinorMag);
  }

  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  thrust_h->Write("", TObject::kOverwrite);
  thrustMajor_h->Write("", TObject::kOverwrite);
  thrustMinor_h->Write("", TObject::kOverwrite);
  thrustDotMajor_h->Write("", TObject::kOverwrite);
  thrustDotMinor_h->Write("", TObject::kOverwrite);
  thrustMajorDotMinor_h->Write("", TObject::kOverwrite);
  thrustMagMinusThrustMajorMag_h->Write("", TObject::kOverwrite);
  thrustMajorMagMinusThrustMinorMag_h->Write("", TObject::kOverwrite);

  delete thrust_h;
  delete thrustMajor_h;
  delete thrustMinor_h;
  delete thrustDotMajor_h;
  delete thrustDotMinor_h;
  delete thrustMajorDotMinor_h;
  delete thrustMagMinusThrustMajorMag_h;
  delete thrustMajorMagMinusThrustMinorMag_h;
  
  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/checkThrustMajorMinor.exe <inName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += checkThrustMajorMinor(argv[1]);
  return retVal;
}
