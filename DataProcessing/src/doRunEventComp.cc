#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TDatime.h"
#include "TObjArray.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMath.h"

void doFill(TH1F* hist1_p, TH1F* hist2_p, TH1F* histDelta_p, Float_t val1, Float_t val2)
{
  hist1_p->Fill(val1);
  hist2_p->Fill(val2);
  histDelta_p->Fill(val1 - val2);
  return;
}

void doFillArr(TH1F* hist1_p, TH1F* hist2_p, TH1F* histDelta_p, Int_t size1, Int_t size2, Float_t val1[], Float_t val2[])
{
  for(Int_t pI = 0; pI < size1; ++pI){hist1_p->Fill(val1[pI]);}
  for(Int_t pI = 0; pI < size2; ++pI){hist2_p->Fill(val2[pI]);}
  if(size1 == size2){
    for(Int_t pI = 0; pI < size1; ++pI){histDelta_p->Fill(val1[pI] - val2[pI]);}
  }
  return;
}

void doFillArr(TH1F* hist1_p, TH1F* hist2_p, TH1F* histDelta_p, Int_t size1, Int_t size2, Int_t val1[], Int_t val2[])
{
  for(Int_t pI = 0; pI < size1; ++pI){hist1_p->Fill(val1[pI]);}
  for(Int_t pI = 0; pI < size2; ++pI){hist2_p->Fill(val2[pI]);}
  if(size1 == size2){
    for(Int_t pI = 0; pI < size1; ++pI){histDelta_p->Fill(val1[pI] - val2[pI]);}
  }
  return;
}

int doRunEventComp(const std::string inFileName1, const std::string inFileName2, std::string outFileName = "")
{
  TDatime* date = new TDatime();
  std::string inFileNameCombo = inFileName1;
  while(inFileNameCombo.find("/") != std::string::npos) inFileNameCombo.replace(0, inFileNameCombo.find("/")+1, "");
  while(inFileNameCombo.find(".root") != std::string::npos) inFileNameCombo.replace(inFileNameCombo.find(".root"), 5, "");
  
  std::string appOutFileName = inFileName2;
  while(appOutFileName.find("/") != std::string::npos) appOutFileName.replace(0, appOutFileName.find("/")+1, "");
  inFileNameCombo = inFileNameCombo + "_" + appOutFileName;

  while(inFileNameCombo.find(".root") != std::string::npos) inFileNameCombo.replace(inFileNameCombo.find(".root"), 5, "");    

  if(outFileName.size() == 0) outFileName = inFileNameCombo;
  if(outFileName.find(".root") != std::string::npos){outFileName.replace(outFileName.find(".root"), 5, "");}

  outFileName = outFileName + "_DORUNLUMIEVTCOMP_" + std::to_string(date->GetDate()) + ".root";

  Int_t RunNo1_;
  Int_t EventNo1_;
  std::map<Int_t, Int_t> runEvtMap1;
  std::map<Int_t, Int_t> runMinEvtMap1;
  std::map<Int_t, Int_t> runMaxEvtMap1;

  Int_t RunNo2_;
  Int_t EventNo2_;
  std::map<Int_t, Int_t> runEvtMap2;
  std::map<Int_t, Int_t> runMinEvtMap2;
  std::map<Int_t, Int_t> runMaxEvtMap2;

  std::map<ULong64_t, Int_t> f1RunEvtToEntry;

  TFile* inFile1_p = new TFile(inFileName1.c_str(), "READ");
  TTree* inTree1_p = (TTree*)inFile1_p->Get("t");
  inTree1_p->SetBranchStatus("*", 0);
  inTree1_p->SetBranchStatus("RunNo", 1);
  inTree1_p->SetBranchStatus("EventNo", 1);

  inTree1_p->SetBranchAddress("RunNo", &RunNo1_);
  inTree1_p->SetBranchAddress("EventNo", &EventNo1_);

  Int_t runNoMin = inTree1_p->GetMinimum("RunNo");
  Int_t runNoMax = inTree1_p->GetMaximum("RunNo");

  Int_t eventNoMin = inTree1_p->GetMinimum("EventNo");
  Int_t eventNoMax = inTree1_p->GetMaximum("EventNo");

  for(Int_t entry = 0; entry < inTree1_p->GetEntries(); ++entry){
    inTree1_p->GetEntry(entry);

    ULong64_t key = RunNo1_*10000000 + EventNo1_;
    if(f1RunEvtToEntry.find(key) != f1RunEvtToEntry.end()){
      std::cout << "Uhoh key duplication \'" << key << "\'..." << std::endl;
    }
    f1RunEvtToEntry[key] = entry;
  }

  TFile* inFile2_p = new TFile(inFileName2.c_str(), "READ");
  TTree* inTree2_p = (TTree*)inFile2_p->Get("t");
  inTree2_p->SetBranchStatus("*", 0);
  inTree2_p->SetBranchStatus("RunNo", 1);
  inTree2_p->SetBranchStatus("EventNo", 1);

  inTree2_p->SetBranchAddress("RunNo", &RunNo2_);
  inTree2_p->SetBranchAddress("EventNo", &EventNo2_);

  runNoMin = TMath::Min(runNoMin, (Int_t)inTree2_p->GetMinimum("RunNo"));
  runNoMax = TMath::Max(runNoMax, (Int_t)inTree2_p->GetMaximum("RunNo"));

  eventNoMin = TMath::Min(eventNoMin, (Int_t)inTree1_p->GetMinimum("EventNo"));
  eventNoMax = TMath::Max(eventNoMax, (Int_t)inTree1_p->GetMaximum("EventNo"));

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");

  const Int_t nRunNoBins = runNoMax - runNoMin + 1;
  TH1F* histRunNo1_p = new TH1F("histRunNo1_h", "RunNo 1", nRunNoBins, runNoMin - 0.5, runNoMax + 0.5); 
  TH1F* histRunNo2_p = new TH1F("histRunNo2_h", "RunNo 2", nRunNoBins, runNoMin - 0.5, runNoMax + 0.5);

  TH1F* histRunNo1_ReducedBins_p = new TH1F("histRunNo1_ReducedBins_h", "RunNo 1", 100, runNoMin - 0.5, runNoMax + 0.5); 
  TH1F* histRunNo2_ReducedBins_p = new TH1F("histRunNo2_ReducedBins_h", "RunNo 2", 100, runNoMin - 0.5, runNoMax + 0.5);


  for(Int_t entry = 0; entry < inTree1_p->GetEntries(); ++entry){
    inTree1_p->GetEntry(entry);

    runEvtMap1.find(RunNo1_) == runEvtMap1.end() ? runEvtMap1[RunNo1_] = 0 : ++runEvtMap1[RunNo1_];
    runMinEvtMap1.find(RunNo1_) == runMinEvtMap1.end() ? runMinEvtMap1[RunNo1_] = EventNo1_ : runMinEvtMap1[RunNo1_] = TMath::Min(EventNo1_, runMinEvtMap1[RunNo1_]);
    runMaxEvtMap1.find(RunNo1_) == runMaxEvtMap1.end() ? runMaxEvtMap1[RunNo1_] = EventNo1_ : runMaxEvtMap1[RunNo1_] = TMath::Max(EventNo1_, runMaxEvtMap1[RunNo1_]);

    histRunNo1_p->Fill(RunNo1_);
    histRunNo1_ReducedBins_p->Fill(RunNo1_);
  }

  for(Int_t entry = 0; entry < inTree2_p->GetEntries(); ++entry){
    inTree2_p->GetEntry(entry);

    runEvtMap2.find(RunNo2_) == runEvtMap2.end() ? runEvtMap2[RunNo2_] = 0 : ++runEvtMap2[RunNo2_];
    runMinEvtMap2.find(RunNo2_) == runMinEvtMap2.end() ? runMinEvtMap2[RunNo2_] = EventNo2_ : runMinEvtMap2[RunNo2_] = TMath::Min(EventNo2_, runMinEvtMap2[RunNo2_]);
    runMaxEvtMap2.find(RunNo2_) == runMaxEvtMap2.end() ? runMaxEvtMap2[RunNo2_] = EventNo2_ : runMaxEvtMap2[RunNo2_] = TMath::Max(EventNo2_, runMaxEvtMap2[RunNo2_]);

    histRunNo2_p->Fill(RunNo2_);
    histRunNo2_ReducedBins_p->Fill(RunNo2_);
  }

  outFile_p->cd();
  const Int_t nRun1Hist_ = runEvtMap1.size(); 
  TH1F* histEventNo1_p[nRun1Hist_];
  std::map<Int_t, Int_t> runHistPos1;

  const Int_t nRun2Hist_ = runEvtMap2.size(); 
  TH1F* histEventNo2_p[nRun2Hist_];
  std::map<Int_t, Int_t> runHistPos2;

  Int_t rI = 0;
  for(std::map<Int_t, Int_t>::iterator it = runEvtMap1.begin(); it != runEvtMap1.end(); ++it){
    histEventNo1_p[rI] = new TH1F(("histEventNo1_Run" + std::to_string(it->first) + "_h").c_str(), ";Events;Counts", 100, runMinEvtMap1[it->first]-0.5, runMaxEvtMap1[it->first]+0.5);
    runHistPos1[it->first] = rI;
    ++rI;
  }

  rI = 0;
  for(std::map<Int_t, Int_t>::iterator it = runEvtMap2.begin(); it != runEvtMap2.end(); ++it){
    histEventNo2_p[rI] = new TH1F(("histEventNo2_Run" + std::to_string(it->first) + "_h").c_str(), ";Events;Counts", 100, runMinEvtMap2[it->first]-0.5, runMaxEvtMap2[it->first]+0.5);
    runHistPos2[it->first] = rI;
    ++rI;
  }

  for(Int_t entry = 0; entry < inTree1_p->GetEntries(); ++entry){
    inTree1_p->GetEntry(entry);
    histEventNo1_p[runHistPos1[RunNo1_]]->Fill(EventNo1_);
  }

  for(Int_t entry = 0; entry < inTree2_p->GetEntries(); ++entry){
    inTree2_p->GetEntry(entry);
    histEventNo2_p[runHistPos2[RunNo2_]]->Fill(EventNo2_);
  }
  


  const Int_t nEntries = inTree2_p->GetEntries();
  std::cout << "Doing full processing..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%1000 == 0) std::cout << " Entry: " << entry << "/" << nEntries << std::endl;
    inTree2_p->GetEntry(entry);

    ULong64_t key = RunNo2_*10000000 + EventNo2_;
    if(f1RunEvtToEntry.find(key) == f1RunEvtToEntry.end()){
      std::cout << "Uhoh missing key \'" << key << "\'..." << std::endl;
      continue;
    }
    
    inTree1_p->GetEntry(f1RunEvtToEntry[key]);
  }

  inFile2_p->Close();
  delete inFile2_p;

  inFile1_p->Close();
  delete inFile1_p;

  outFile_p->cd();

  histRunNo1_p->Write("", TObject::kOverwrite);
  histRunNo2_p->Write("", TObject::kOverwrite);

  histRunNo1_ReducedBins_p->Write("", TObject::kOverwrite);
  histRunNo2_ReducedBins_p->Write("", TObject::kOverwrite);

  outFile_p->mkdir("evtFiles1");
  outFile_p->cd("evtFiles1");

  for(Int_t i = 0; i < nRun1Hist_; ++i){
    histEventNo1_p[i]->Write("", TObject::kOverwrite);
    delete histEventNo1_p[i];
  }

  outFile_p->cd();
  outFile_p->mkdir("evtFiles2");
  outFile_p->cd("evtFiles2");

  for(Int_t i = 0; i < nRun2Hist_; ++i){
    histEventNo2_p[i]->Write("", TObject::kOverwrite);
    delete histEventNo2_p[i];
  }  

  outFile_p->cd();

  delete histRunNo1_p;
  delete histRunNo2_p;

  delete histRunNo1_ReducedBins_p;
  delete histRunNo2_ReducedBins_p;

  TNamed nameFile1("nameFile1", inFileName1.c_str());
  TNamed nameFile2("nameFile2", inFileName2.c_str());

  nameFile1.Write("", TObject::kOverwrite);
  nameFile2.Write("", TObject::kOverwrite);

  outFile_p->Close();
  delete outFile_p;

  delete date;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3 && argc != 4){
    std::cout << "Usage ./doRunEventComp.exe <inFileName1> <inFileName2> <outFileName-optional>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 3) retVal += doRunEventComp(argv[1], argv[2]);
  else if(argc == 4) retVal += doRunEventComp(argv[1], argv[2], argv[3]);
  return retVal;
}
