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

#include "include/particleData.h"
#include "include/eventData.h"
#include "include/jetData.h"
#include "include/returnRootFileContentsList.h"
#include "include/removeVectorDuplicates.h"
#include "include/histDefUtility.h"
#include "include/smartJetName.h"
#include "include/doLocalDebug.h"

void doFill(TH1F* hist1_p, TH1F* hist2_p, TH1F* histDelta_p, Bool_t val1, Bool_t val2)
{
  Float_t val1Out, val2Out;
  val1 ? val1Out = 1 : val1Out = 0;
  val2 ? val2Out = 2 : val1Out = 0;

  hist1_p->Fill(val1Out);
  hist2_p->Fill(val2Out);
  if(val1-val2 > histDelta_p->GetBinLowEdge(histDelta_p->GetNbinsX()+1)) histDelta_p->Fill(histDelta_p->GetBinCenter(histDelta_p->GetNbinsX()));
  else if(val1-val2 < histDelta_p->GetBinLowEdge(1)) histDelta_p->Fill(histDelta_p->GetBinCenter(1));
  else histDelta_p->Fill(val1-val2);
	  
  return;
}

void doFill(TH1F* hist1_p, TH1F* hist2_p, TH1F* histDelta_p, Int_t val1, Int_t val2)
{
  hist1_p->Fill(val1);
  hist2_p->Fill(val2);
  if(val1-val2 > histDelta_p->GetBinLowEdge(histDelta_p->GetNbinsX()+1)) histDelta_p->Fill(histDelta_p->GetBinCenter(histDelta_p->GetNbinsX()));
  else if(val1-val2 < histDelta_p->GetBinLowEdge(1)) histDelta_p->Fill(histDelta_p->GetBinCenter(1));
  else histDelta_p->Fill(val1-val2);
	  
  return;
}

void doFill(TH1F* hist1_p, TH1F* hist2_p, TH1F* histDelta_p, Float_t val1, Float_t val2)
{
  hist1_p->Fill(val1);
  hist2_p->Fill(val2);
  if(val1-val2 > histDelta_p->GetBinLowEdge(histDelta_p->GetNbinsX()+1)) histDelta_p->Fill(histDelta_p->GetBinCenter(histDelta_p->GetNbinsX()));
  else if(val1-val2 < histDelta_p->GetBinLowEdge(1)) histDelta_p->Fill(histDelta_p->GetBinCenter(1));
  else histDelta_p->Fill(val1-val2);
	  
  return;
}

void doFillArr(TH1F* hist1_p, TH1F* hist2_p, TH1F* histDelta_p, Int_t size1, Int_t size2, Float_t val1[], Float_t val2[])
{
  for(Int_t pI = 0; pI < size1; ++pI){hist1_p->Fill(val1[pI]);}
  for(Int_t pI = 0; pI < size2; ++pI){hist2_p->Fill(val2[pI]);}
  if(size1 == size2){
    for(Int_t pI = 0; pI < size1; ++pI){
      if(val1[pI]-val2[pI] > histDelta_p->GetBinLowEdge(histDelta_p->GetNbinsX()+1)) histDelta_p->Fill(histDelta_p->GetBinCenter(histDelta_p->GetNbinsX()));
      else if(val1[pI]-val2[pI] < histDelta_p->GetBinLowEdge(1)) histDelta_p->Fill(histDelta_p->GetBinCenter(1));
      else histDelta_p->Fill(val1[pI]-val2[pI]);
    }
  }
  return;
}

void doFillArr(TH1F* hist1_p, TH1F* hist2_p, TH1F* histDelta_p, Int_t size1, Int_t size2, Int_t val1[], Int_t val2[])
{
  for(Int_t pI = 0; pI < size1; ++pI){hist1_p->Fill(val1[pI]);}
  for(Int_t pI = 0; pI < size2; ++pI){hist2_p->Fill(val2[pI]);}
  if(size1 == size2){
    for(Int_t pI = 0; pI < size1; ++pI){
      if(val1[pI]-val2[pI] > histDelta_p->GetBinLowEdge(histDelta_p->GetNbinsX()+1)) histDelta_p->Fill(histDelta_p->GetBinCenter(histDelta_p->GetNbinsX()));
      else if(val1[pI]-val2[pI] < histDelta_p->GetBinLowEdge(1)) histDelta_p->Fill(histDelta_p->GetBinCenter(1));
      else histDelta_p->Fill(val1[pI]-val2[pI]);      
    }
  }
  return;
}


int doComparison(const std::string inFileName1, const std::string inFileName2, std::string outFileName = "")
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

  outFileName = outFileName + "_DOCOMP_" + std::to_string(date->GetDate()) + ".root";

  std::vector<std::string> listOfCompBranches;
  std::vector<std::string> listExcluded;
  std::vector<double> branchMins;
  std::vector<double> branchMaxs;

  particleData pData1;
  particleData pData2;

  eventData eData1;
  eventData eData2;

  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::map<ULong64_t, Int_t> f1RunEvtToEntry;

  TFile* inFile1_p = new TFile(inFileName1.c_str(), "READ");
  TTree* inTree1_p = (TTree*)inFile1_p->Get("t");
  inTree1_p->SetBranchStatus("*", 0);
  inTree1_p->SetBranchStatus("RunNo", 1);
  inTree1_p->SetBranchStatus("EventNo", 1);

  inTree1_p->SetBranchAddress("RunNo", &pData1.RunNo);
  inTree1_p->SetBranchAddress("EventNo", &pData1.EventNo);

  for(Int_t entry = 0; entry < inTree1_p->GetEntries(); ++entry){
    inTree1_p->GetEntry(entry);

    ULong64_t key = pData1.RunNo*10000000 + pData1.EventNo;
    if(f1RunEvtToEntry.find(key) != f1RunEvtToEntry.end()){
      std::cout << "Uhoh key duplication \'" << key << "\'..." << std::endl;
    }
    f1RunEvtToEntry[key] = entry;
  }

  TObjArray* list1_p = (TObjArray*)inTree1_p->GetListOfBranches();
  for(Int_t i = 0; i < list1_p->GetEntries(); ++i){
    listOfCompBranches.push_back(list1_p->At(i)->GetName());
    branchMins.push_back(inTree1_p->GetMinimum(list1_p->At(i)->GetName()));
    branchMaxs.push_back(inTree1_p->GetMaximum(list1_p->At(i)->GetName()));
  }

  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::vector<std::string> listOfJetTrees1 = returnRootFileContentsList(inFile1_p, "TTree", "JetTree");
  removeVectorDuplicates(&listOfJetTrees1);
  std::vector< std::vector<std::string> > listOfJetTreeBranches;
  std::vector< std::vector<double> > listOfJetTreeMins;
  std::vector< std::vector<double> > listOfJetTreeMaxs;
  for(unsigned int jI = 0; jI < listOfJetTrees1.size(); ++jI){
    TTree* tempTree_p = (TTree*)inFile1_p->Get(listOfJetTrees1.at(jI).c_str());
    TObjArray* tempList = (TObjArray*)tempTree_p->GetListOfBranches();
    std::vector<std::string> tempV;
    std::vector<double> tempMins;
    std::vector<double> tempMaxs;

    for(Int_t i = 0; i < tempList->GetEntries(); ++i){
      tempV.push_back(tempList->At(i)->GetName());
      tempMins.push_back(tempTree_p->GetMinimum(tempList->At(i)->GetName()));
      tempMaxs.push_back(tempTree_p->GetMaximum(tempList->At(i)->GetName()));
    }

    listOfJetTreeBranches.push_back(tempV);
    listOfJetTreeMins.push_back(tempMins);
    listOfJetTreeMaxs.push_back(tempMaxs);
  }

  inFile1_p->Close();
  delete inFile1_p;

  TFile* inFile2_p = new TFile(inFileName2.c_str(), "READ");
  TTree* inTree2_p = (TTree*)inFile2_p->Get("t");
  TObjArray* list2_p = (TObjArray*)inTree2_p->GetListOfBranches();
  std::vector<bool> inTreeBools;

  for(Int_t i = 0; i < list2_p->GetEntries(); ++i){
    inTreeBools.push_back(false);
  }

  unsigned int pos = 0;
  while(pos < listOfCompBranches.size()){
    std::string tempStr = listOfCompBranches.at(pos);

    bool isInArr = false;
    for(Int_t i = 0; i < list2_p->GetEntries(); ++i){
      std::string tempStr2 = list2_p->At(i)->GetName();

      if(tempStr.size() == tempStr2.size() && tempStr.find(tempStr2) != std::string::npos){
	isInArr = true;
	inTreeBools.at(i) = true;
	break;
      }
    }

    if(!isInArr){
      listExcluded.push_back(listOfCompBranches.at(pos));
      listOfCompBranches.erase(listOfCompBranches.begin()+pos);
      branchMins.erase(branchMins.begin()+pos);
      branchMaxs.erase(branchMaxs.begin()+pos);      
    }
    else{
      Double_t tempMin = inTree2_p->GetMinimum(listOfCompBranches.at(pos).c_str());
      Double_t tempMax = inTree2_p->GetMaximum(listOfCompBranches.at(pos).c_str());

      if(branchMins.at(pos) > tempMin) branchMins.at(pos) = tempMin;
      if(branchMaxs.at(pos) < tempMax) branchMaxs.at(pos) = tempMax;
      ++pos;
    }
  }

  for(Int_t i = 0; i < list2_p->GetEntries(); ++i){
    std::string listName = list2_p->At(i)->GetName();
    if(!inTreeBools.at(i)) listExcluded.push_back(listName);
  }

  std::cout << "Set of var not in both trees: " << std::endl;
  for(unsigned int i = 0; i < listExcluded.size(); ++i){
    std::cout << " " << i << "/" << listExcluded.size() << ": " << listExcluded.at(i) << std::endl;
  }

  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::vector<std::string> listOfJetTrees2 = returnRootFileContentsList(inFile2_p, "TTree", "JetTree");
  removeVectorDuplicates(&listOfJetTrees2);

  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  
  unsigned int lI = 0;
  while(lI < listOfJetTrees1.size()){
    bool isGood = false;

  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    for(unsigned int lI2 = 0; lI2 < listOfJetTrees2.size(); ++lI2){
      if(listOfJetTrees1.at(lI).size() == listOfJetTrees2.at(lI2).size() && listOfJetTrees1.at(lI).find(listOfJetTrees2.at(lI2)) != std::string::npos){
	isGood = true;
	std::string tempLOJT2 = listOfJetTrees2.at(lI);
	listOfJetTrees2.at(lI) = listOfJetTrees2.at(lI2);
	listOfJetTrees2.at(lI2) = tempLOJT2;
	break;
      }
    }

    if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    if(!isGood){
      std::string tempLOJT2 = smartJetName(listOfJetTrees1.at(lI), listOfJetTrees2);

      if(tempLOJT2.size() != 0){
	Int_t prevPos = -1;
	for(unsigned int i = 0; i < listOfJetTrees2.size(); ++i){
	  if(sameString(tempLOJT2, listOfJetTrees2.at(i))){
	    prevPos = i;
	    break;
	  }
	}

	if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	isGood = true;
	listOfJetTrees2.at(prevPos) = listOfJetTrees2.at(lI);
	listOfJetTrees2.at(lI) = tempLOJT2;
      }
    }

    if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    if(!isGood){
      listOfJetTrees1.erase(listOfJetTrees1.begin()+lI);
      listOfJetTreeBranches.erase(listOfJetTreeBranches.begin()+lI);
      listOfJetTreeMins.erase(listOfJetTreeMins.begin()+lI);
      listOfJetTreeMaxs.erase(listOfJetTreeMaxs.begin()+lI);
    }
    else{
      TTree* tempTree_p = (TTree*)inFile2_p->Get(listOfJetTrees2.at(lI).c_str());
      TObjArray* tempList_p = (TObjArray*)tempTree_p->GetListOfBranches();

      unsigned int j = 0;
      while(j < listOfJetTreeBranches.at(lI).size()){
	bool branchIsGood = false;

	for(Int_t i = 0; i < tempList_p->GetEntries(); ++i){
	  std::string tempBranchName = tempList_p->At(i)->GetName();
	  if(tempBranchName.size() == listOfJetTreeBranches.at(lI).at(j).size() && tempBranchName.find(listOfJetTreeBranches.at(lI).at(j)) != std::string::npos){
	    branchIsGood = true;
	    break;
	  }
	}

	if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	if(!branchIsGood){
	  listOfJetTreeBranches.at(lI).erase(listOfJetTreeBranches.at(lI).begin() + j);
	  listOfJetTreeMins.at(lI).erase(listOfJetTreeMins.at(lI).begin() + j);
	  listOfJetTreeMaxs.at(lI).erase(listOfJetTreeMaxs.at(lI).begin() + j);
	}
	else{
	  if(tempTree_p->GetMinimum(listOfJetTreeBranches.at(lI).at(j).c_str()) < listOfJetTreeMins.at(lI).at(j))
	    listOfJetTreeMins.at(lI).at(j) = tempTree_p->GetMinimum(listOfJetTreeBranches.at(lI).at(j).c_str());
	  if(tempTree_p->GetMaximum(listOfJetTreeBranches.at(lI).at(j).c_str()) > listOfJetTreeMaxs.at(lI).at(j))
	    listOfJetTreeMaxs.at(lI).at(j) = tempTree_p->GetMaximum(listOfJetTreeBranches.at(lI).at(j).c_str());
	  
	  ++j;
	}
      }

      ++lI;
    }
  }

  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  inFile2_p->Close();
  delete inFile2_p;

  const Int_t nVarToComp = listOfCompBranches.size();
  const Int_t nJetTrees = listOfJetTrees1.size();
  jetData jData1[nJetTrees];
  jetData jData2[nJetTrees];

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1F* hist1_p[nVarToComp];
  TH1F* hist2_p[nVarToComp];
  TH1F* hist_Rat1Over2_p[nVarToComp];
  TH1F* hist_Delta1From2_EvtByEvt_p[nVarToComp];

  std::vector< std::vector< TH1F* > > hist1_Jet_p;
  std::vector< std::vector< TH1F* > > hist2_Jet_p;
  std::vector< std::vector< TH1F* > > hist_Rat1Over2_Jet_p;
  std::vector< std::vector< TH1F* > > hist_Delta1From2_EvtByEvt_Jet_p;


  for(Int_t i = 0; i < nVarToComp; ++i){
    double tempInterval = branchMaxs.at(i) - branchMins.at(i);
    branchMaxs.at(i) += tempInterval/10.;
    branchMins.at(i) -= tempInterval/10.;

    if(tempInterval == 0){
      branchMaxs.at(i) += 1;
      branchMins.at(i) -= 1;
    }
    std::cout << branchMins.at(i) << ", " << branchMaxs.at(i) << std::endl;

    hist1_p[i] = new TH1F((listOfCompBranches.at(i) + "_file1_h").c_str(), (";" + listOfCompBranches.at(i) + ";Counts").c_str(), 100, branchMins.at(i), branchMaxs.at(i));
    hist2_p[i] = new TH1F((listOfCompBranches.at(i) + "_file2_h").c_str(), (";" + listOfCompBranches.at(i) + ";Counts").c_str(), 100, branchMins.at(i), branchMaxs.at(i));
    hist_Rat1Over2_p[i] = new TH1F((listOfCompBranches.at(i) + "_file_Rat1Over2_h").c_str(), (";" + listOfCompBranches.at(i) + ";File1/File2").c_str(), 100, branchMins.at(i), branchMaxs.at(i));
    hist_Delta1From2_EvtByEvt_p[i] = new TH1F((listOfCompBranches.at(i) + "_Delta1From2_EvtByEvt_h").c_str(), (";" + listOfCompBranches.at(i) + "_{File 1}" + "-" + listOfCompBranches.at(i) + "_{File 2};Counts").c_str(), 100, -1, 1);
    centerTitles({hist1_p[i], hist2_p[i], hist_Rat1Over2_p[i], hist_Delta1From2_EvtByEvt_p[i]});
  }

  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(Int_t jI = 0; jI < nJetTrees; ++jI){
    hist1_Jet_p.push_back( {} );
    hist2_Jet_p.push_back( {} );
    hist_Rat1Over2_Jet_p.push_back( {} );
    hist_Delta1From2_EvtByEvt_Jet_p.push_back( {} );

    std::string jetStr = listOfJetTrees1.at(jI);

    for(unsigned int bI = 0; bI < listOfJetTreeBranches.at(jI).size(); ++bI){
      std::string brStr = listOfJetTreeBranches.at(jI).at(bI);
      std::string nameStr = brStr + "_" + jetStr;

      Double_t tempMin = listOfJetTreeMins.at(jI).at(bI);
      Double_t tempMax = listOfJetTreeMaxs.at(jI).at(bI);
      Double_t interval = (tempMax - tempMin)/10.;

      if(interval < 0.0001){
	tempMax += 1;
	tempMin -= 1;
      }
      else{
	tempMax += interval;
	tempMin -= interval;
      }

      hist1_Jet_p.at(jI).push_back(0);
      hist2_Jet_p.at(jI).push_back(0);
      hist_Rat1Over2_Jet_p.at(jI).push_back(0);
      hist_Delta1From2_EvtByEvt_Jet_p.at(jI).push_back(0);

      hist1_Jet_p.at(jI).at(bI) = new TH1F((nameStr + "_file1_h").c_str(), (";" + brStr + " (" + jetStr + ");Counts").c_str(), 100, tempMin, tempMax);
      hist2_Jet_p.at(jI).at(bI) = new TH1F((nameStr + "_file2_h").c_str(), (";" + brStr + " (" + jetStr + ");Counts").c_str(), 100, tempMin, tempMax);
      hist_Rat1Over2_Jet_p.at(jI).at(bI) = new TH1F((nameStr + "_file_Rat1Over2_h").c_str(), (";" + brStr + " (" + jetStr + ");File1/File2").c_str(), 100, tempMin, tempMax);
      hist_Delta1From2_EvtByEvt_Jet_p.at(jI).at(bI) = new TH1F((nameStr + "_Delta1From2_EvtByEvt_h").c_str(), (";" + brStr +"_{File 1} - " + brStr + "_{File 2}" + " (" + jetStr + ");Counts").c_str(), 100, -1, 1);
      centerTitles({hist1_Jet_p.at(jI).at(bI), hist2_Jet_p.at(jI).at(bI), hist_Rat1Over2_Jet_p.at(jI).at(bI), hist_Delta1From2_EvtByEvt_Jet_p.at(jI).at(bI)});
    }
  }

  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  inFile1_p = new TFile(inFileName1.c_str(), "READ");
  inTree1_p = (TTree*)inFile1_p->Get("t");
  inTree1_p->SetBranchStatus("*", 0);
  pData1.SetStatusAndAddressRead(inTree1_p, listOfCompBranches);
  eData1.SetStatusAndAddressRead(inTree1_p, listOfCompBranches);

  TTree* jetTree1_p[nJetTrees];
  for(Int_t jI = 0; jI < nJetTrees; ++jI){
    jetTree1_p[jI] = (TTree*)inFile1_p->Get(listOfJetTrees1.at(jI).c_str());
    jetTree1_p[jI]->SetBranchStatus("*", 0);
    jData1[jI].SetStatusAndAddressRead(jetTree1_p[jI], listOfJetTreeBranches.at(jI));
  }

  inFile2_p = new TFile(inFileName2.c_str(), "READ");
  inTree2_p = (TTree*)inFile2_p->Get("t");
  inTree2_p->SetBranchStatus("*", 0);
  pData2.SetStatusAndAddressRead(inTree2_p, listOfCompBranches);
  eData2.SetStatusAndAddressRead(inTree2_p, listOfCompBranches);

  TTree* jetTree2_p[nJetTrees];
  for(Int_t jI = 0; jI < nJetTrees; ++jI){
    jetTree2_p[jI] = (TTree*)inFile2_p->Get(listOfJetTrees2.at(jI).c_str());
    jetTree2_p[jI]->SetBranchStatus("*", 0);
    jData2[jI].SetStatusAndAddressRead(jetTree2_p[jI], listOfJetTreeBranches.at(jI));
  }

  const Int_t nEntries = inTree2_p->GetEntries();

  std::cout << "Doing full processing..." << std::endl;

  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%1000 == 0) std::cout << " Entry: " << entry << "/" << nEntries << std::endl;
    inTree2_p->GetEntry(entry);

    for(Int_t jI = 0; jI < nJetTrees; ++jI){jetTree2_p[jI]->GetEntry(entry);}

    ULong64_t key = pData2.RunNo*10000000 + pData2.EventNo;
    if(f1RunEvtToEntry.find(key) == f1RunEvtToEntry.end()){
      std::cout << "Uhoh missing key \'" << key << "\'..." << std::endl;
      continue;
    }
    
    inTree1_p->GetEntry(f1RunEvtToEntry[key]);

    for(Int_t jI = 0; jI < nJetTrees; ++jI){jetTree1_p[jI]->GetEntry(f1RunEvtToEntry[key]);}

    for(unsigned int lI = 0; lI < listOfCompBranches.size(); ++lI){
      std::string tempS = listOfCompBranches.at(lI);
      if(tempS.find("nParticle") != std::string::npos && tempS.size() == std::string("nParticle").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle);
      else if(tempS.find("EventNo") != std::string::npos && tempS.size() == std::string("EventNo").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.EventNo, pData2.EventNo);
      else if(tempS.find("RunNo") != std::string::npos && tempS.size() == std::string("RunNo").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.RunNo, pData2.RunNo);
      else if(tempS.find("year") != std::string::npos && tempS.size() == std::string("year").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.year, pData2.year);
      else if(tempS.find("process") != std::string::npos && tempS.size() == std::string("process").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.process, pData2.process);
      else if(tempS.find("Energy") != std::string::npos && tempS.size() == std::string("Energy").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.Energy, pData2.Energy);
      else if(tempS.find("bFlag") != std::string::npos && tempS.size() == std::string("bFlag").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.bFlag, pData2.bFlag);
      else if(tempS.find("bx") != std::string::npos && tempS.size() == std::string("bx").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.bx, pData2.bx);
      else if(tempS.find("by") != std::string::npos && tempS.size() == std::string("by").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.by, pData2.by);
      else if(tempS.find("ebx") != std::string::npos && tempS.size() == std::string("ebx").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.ebx, pData2.ebx);
      else if(tempS.find("eby") != std::string::npos && tempS.size() == std::string("eby").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.eby, pData2.eby);
      else if(tempS.find("px") != std::string::npos && tempS.size() == std::string("px").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.px, pData2.px);
      else if(tempS.find("py") != std::string::npos && tempS.size() == std::string("py").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.py, pData2.py);
      else if(tempS.find("pz") != std::string::npos && tempS.size() == std::string("pz").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.pz, pData2.pz);
      else if(tempS.find("pt") != std::string::npos && tempS.size() == std::string("pt").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.pt, pData2.pt);
      else if(tempS.find("pmag") != std::string::npos && tempS.size() == std::string("pmag").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.pmag, pData2.pmag);
      else if(tempS.find("eta") != std::string::npos && tempS.size() == std::string("eta").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.eta, pData2.eta);
      else if(tempS.find("rap") != std::string::npos && tempS.size() == std::string("rap").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.rap, pData2.rap);
      else if(tempS.find("theta") != std::string::npos && tempS.size() == std::string("theta").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.theta, pData2.theta);
      else if(tempS.find("phi") != std::string::npos && tempS.size() == std::string("phi").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.phi, pData2.phi);
      else if(tempS.find("mass") != std::string::npos && tempS.size() == std::string("mass").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.mass, pData2.mass);
      else if(tempS.find("charge") != std::string::npos && tempS.size() == std::string("charge").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.charge, pData2.charge);
      else if(tempS.find("pwflag") != std::string::npos && tempS.size() == std::string("pwflag").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.pwflag, pData2.pwflag);
      else if(tempS.find("pid") != std::string::npos && tempS.size() == std::string("pid").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.pid, pData2.pid);
      else if(tempS.find("d0") != std::string::npos && tempS.size() == std::string("d0").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.d0, pData2.d0);
      else if(tempS.find("z0") != std::string::npos && tempS.size() == std::string("z0").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.z0, pData2.z0);
      else if(tempS.find("ntpc") != std::string::npos && tempS.size() == std::string("ntpc").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.ntpc, pData2.ntpc);
      else if(tempS.find("pt_wrtThr") != std::string::npos && tempS.size() == std::string("pt_wrtThr").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.pt_wrtThr, pData2.pt_wrtThr);
      else if(tempS.find("eta_wrtThr") != std::string::npos && tempS.size() == std::string("eta_wrtThr").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.eta_wrtThr, pData2.eta_wrtThr);
      else if(tempS.find("rap_wrtThr") != std::string::npos && tempS.size() == std::string("rap_wrtThr").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.rap_wrtThr, pData2.rap_wrtThr);
      else if(tempS.find("theta_wrtThr") != std::string::npos && tempS.size() == std::string("theta_wrtThr").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.theta_wrtThr, pData2.theta_wrtThr);
      else if(tempS.find("phi_wrtThr") != std::string::npos && tempS.size() == std::string("phi_wrtThr").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.phi_wrtThr, pData2.phi_wrtThr);
      else if(tempS.find("pt_wrtChThr") != std::string::npos && tempS.size() == std::string("pt_wrtChThr").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.pt_wrtChThr, pData2.pt_wrtChThr);
      else if(tempS.find("eta_wrtChThr") != std::string::npos && tempS.size() == std::string("eta_wrtChThr").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.eta_wrtChThr, pData2.eta_wrtChThr);
      else if(tempS.find("rap_wrtChThr") != std::string::npos && tempS.size() == std::string("rap_wrtChThr").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.rap_wrtChThr, pData2.rap_wrtChThr);
      else if(tempS.find("theta_wrtChThr") != std::string::npos && tempS.size() == std::string("theta_wrtChThr").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.theta_wrtChThr, pData2.theta_wrtChThr);
      else if(tempS.find("phi_wrtChThr") != std::string::npos && tempS.size() == std::string("phi_wrtChThr").size()) doFillArr(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], pData1.nParticle, pData2.nParticle, pData1.phi_wrtChThr, pData2.phi_wrtChThr);
      else if(tempS.find("passesWW") != std::string::npos && tempS.size() == std::string("passesWW").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], eData1.passesWW, eData2.passesWW);
      else if(tempS.find("missP") != std::string::npos && tempS.size() == std::string("missP").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], eData1.missP, eData2.missP);
      else if(tempS.find("missPt") != std::string::npos && tempS.size() == std::string("missPt").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], eData1.missPt, eData2.missPt);
      else if(tempS.find("missTheta") != std::string::npos && tempS.size() == std::string("missTheta").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], eData1.missTheta, eData2.missTheta);
      else if(tempS.find("missPhi") != std::string::npos && tempS.size() == std::string("missPhi").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], eData1.missPhi, eData2.missPhi);
      else if(tempS.find("missChargedP") != std::string::npos && tempS.size() == std::string("missChargedP").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], eData1.missChargedP, eData2.missChargedP);
      else if(tempS.find("missChargedPt") != std::string::npos && tempS.size() == std::string("missChargedPt").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], eData1.missChargedPt, eData2.missChargedPt);
      else if(tempS.find("missChargedTheta") != std::string::npos && tempS.size() == std::string("missChargedTheta").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], eData1.missChargedTheta, eData2.missChargedTheta);
      else if(tempS.find("missChargedPhi") != std::string::npos && tempS.size() == std::string("missChargedPhi").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], eData1.missChargedPhi, eData2.missChargedPhi);
      else if(tempS.find("nChargedHadrons") != std::string::npos && tempS.size() == std::string("nChargedHadrons").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], eData1.nChargedHadrons, eData2.nChargedHadrons);
      else if(tempS.find("nChargedHadrons_GT0p4") != std::string::npos && tempS.size() == std::string("nChargedHadrons_GT0p4").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], eData1.nChargedHadrons_GT0p4, eData2.nChargedHadrons_GT0p4);
      else if(tempS.find("nChargedHadrons_GT0p4Thrust") != std::string::npos && tempS.size() == std::string("nChargedHadrons_GT0p4Thrust").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], eData1.nChargedHadrons_GT0p4Thrust, eData2.nChargedHadrons_GT0p4Thrust);
      else if(tempS.find("Thrust") != std::string::npos && tempS.size() == std::string("Thrust").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], eData1.Thrust, eData2.Thrust);
      else if(tempS.find("TTheta") != std::string::npos && tempS.size() == std::string("TTheta").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], eData1.TTheta, eData2.TTheta);
      else if(tempS.find("TPhi") != std::string::npos && tempS.size() == std::string("TPhi").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], eData1.TPhi, eData2.TPhi);
      else if(tempS.find("Thrust_charged") != std::string::npos && tempS.size() == std::string("Thrust_charged").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], eData1.Thrust_charged, eData2.Thrust_charged);
      else if(tempS.find("TTheta_charged") != std::string::npos && tempS.size() == std::string("TTheta_charged").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], eData1.TTheta_charged, eData2.TTheta_charged);
      else if(tempS.find("TPhi_charged") != std::string::npos && tempS.size() == std::string("TPhi_charged").size()) doFill(hist1_p[lI], hist2_p[lI], hist_Delta1From2_EvtByEvt_p[lI], eData1.TPhi_charged, eData2.TPhi_charged);
    }

    for(Int_t jI = 0; jI < nJetTrees; ++jI){
      for(unsigned bI = 0; bI < listOfJetTreeBranches.at(jI).size(); ++bI){
	std::string tempS = listOfJetTreeBranches.at(jI).at(bI);
	if(tempS.find("nref") != std::string::npos && tempS.size() == std::string("nref").size()) doFill(hist1_Jet_p.at(jI).at(bI), hist2_Jet_p.at(jI).at(bI), hist_Delta1From2_EvtByEvt_Jet_p.at(jI).at(bI), jData1[jI].nref, jData2[jI].nref);
	else if(tempS.find("jtpt") != std::string::npos && tempS.size() == std::string("jtpt").size()) doFillArr(hist1_Jet_p.at(jI).at(bI), hist2_Jet_p.at(jI).at(bI), hist_Delta1From2_EvtByEvt_Jet_p.at(jI).at(bI), jData1[jI].nref, jData2[jI].nref, jData1[jI].jtpt, jData2[jI].jtpt);
	else if(tempS.find("jteta") != std::string::npos && tempS.size() == std::string("jteta").size()) doFillArr(hist1_Jet_p.at(jI).at(bI), hist2_Jet_p.at(jI).at(bI), hist_Delta1From2_EvtByEvt_Jet_p.at(jI).at(bI), jData1[jI].nref, jData2[jI].nref, jData1[jI].jteta, jData2[jI].jteta);
	else if(tempS.find("jtphi") != std::string::npos && tempS.size() == std::string("jtphi").size()) doFillArr(hist1_Jet_p.at(jI).at(bI), hist2_Jet_p.at(jI).at(bI), hist_Delta1From2_EvtByEvt_Jet_p.at(jI).at(bI), jData1[jI].nref, jData2[jI].nref, jData1[jI].jtphi, jData2[jI].jtphi);
	else if(tempS.find("jtm") != std::string::npos && tempS.size() == std::string("jtm").size()) doFillArr(hist1_Jet_p.at(jI).at(bI), hist2_Jet_p.at(jI).at(bI), hist_Delta1From2_EvtByEvt_Jet_p.at(jI).at(bI), jData1[jI].nref, jData2[jI].nref, jData1[jI].jtm, jData2[jI].jtm);
	else if(tempS.find("jtN") != std::string::npos && tempS.size() == std::string("jtN").size()) doFillArr(hist1_Jet_p.at(jI).at(bI), hist2_Jet_p.at(jI).at(bI), hist_Delta1From2_EvtByEvt_Jet_p.at(jI).at(bI), jData1[jI].nref, jData2[jI].nref, jData1[jI].jtN, jData2[jI].jtN);
	//Note cant handle 2d arrays proper will add later
      }
    }
  }

  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  inFile2_p->Close();
  delete inFile2_p;

  inFile1_p->Close();
  delete inFile1_p;

  outFile_p->cd();

  const Double_t splitPoint = 0.35;

  std::vector<std::string> listOfPdf;
  std::vector<std::string> listOfVar;

  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(Int_t i = 0; i < nVarToComp; ++i){
    hist1_p[i]->Write("", TObject::kOverwrite);
    hist2_p[i]->Write("", TObject::kOverwrite);

    hist1_p[i]->Sumw2();
    hist2_p[i]->Sumw2();
    hist_Rat1Over2_p[i]->Sumw2();

    hist_Rat1Over2_p[i]->Divide(hist1_p[i], hist2_p[i]);
    hist_Rat1Over2_p[i]->Write("", TObject::kOverwrite);

    hist_Delta1From2_EvtByEvt_p[i]->Write("", TObject::kOverwrite);

    TCanvas* canv_p = new TCanvas("canv_c", "canv_c", 1000, 500);
    canv_p->SetTopMargin(0.01);
    canv_p->SetRightMargin(0.01);
    canv_p->SetLeftMargin(0.01);
    canv_p->SetBottomMargin(0.01);

    TPad* pad1_p = new TPad("pad1", "pad1", 0.0, splitPoint, 0.5, 1.0);
    pad1_p->Draw();
    pad1_p->SetTopMargin(0.01);
    pad1_p->SetRightMargin(0.01);
    pad1_p->SetBottomMargin(0.01);
    pad1_p->SetLeftMargin(pad1_p->GetLeftMargin()*1.3);

    TPad* pad2_p = new TPad("pad2", "pad2", 0.0, 0.0, 0.5, splitPoint);
    pad2_p->Draw();
    pad2_p->SetTopMargin(0.01);
    pad2_p->SetRightMargin(0.01);
    pad2_p->SetBottomMargin(pad1_p->GetLeftMargin()*1./splitPoint);
    pad2_p->SetLeftMargin(pad1_p->GetLeftMargin());

    TPad* pad3_p = new TPad("pad3", "pad3", 0.5, 0.0, 1.0, 1.0);
    pad3_p->Draw();
    pad3_p->SetTopMargin(0.01);
    pad3_p->SetRightMargin(0.01);
    pad3_p->SetLeftMargin(pad1_p->GetLeftMargin());
    pad3_p->SetBottomMargin(pad1_p->GetLeftMargin());

    pad1_p->cd();

    hist1_p[i]->GetXaxis()->SetTitleFont(43);
    hist2_p[i]->GetXaxis()->SetTitleFont(43);
    hist_Rat1Over2_p[i]->GetXaxis()->SetTitleFont(43);
    hist_Delta1From2_EvtByEvt_p[i]->GetXaxis()->SetTitleFont(43);
    hist1_p[i]->GetXaxis()->SetTitleSize(20);
    hist2_p[i]->GetXaxis()->SetTitleSize(20);
    hist_Rat1Over2_p[i]->GetXaxis()->SetTitleSize(20);
    hist_Delta1From2_EvtByEvt_p[i]->GetXaxis()->SetTitleSize(20);

    hist1_p[i]->GetYaxis()->SetTitleFont(43);
    hist2_p[i]->GetYaxis()->SetTitleFont(43);
    hist_Rat1Over2_p[i]->GetYaxis()->SetTitleFont(43);
    hist_Delta1From2_EvtByEvt_p[i]->GetYaxis()->SetTitleFont(43);
    hist1_p[i]->GetYaxis()->SetTitleSize(20);
    hist2_p[i]->GetYaxis()->SetTitleSize(20);
    hist_Rat1Over2_p[i]->GetYaxis()->SetTitleSize(20);
    hist_Delta1From2_EvtByEvt_p[i]->GetYaxis()->SetTitleSize(20);


    hist1_p[i]->GetXaxis()->SetLabelFont(43);
    hist2_p[i]->GetXaxis()->SetLabelFont(43);
    hist_Rat1Over2_p[i]->GetXaxis()->SetLabelFont(43);
    hist_Delta1From2_EvtByEvt_p[i]->GetXaxis()->SetLabelFont(43);
    hist1_p[i]->GetXaxis()->SetLabelSize(20);
    hist2_p[i]->GetXaxis()->SetLabelSize(20);
    hist_Rat1Over2_p[i]->GetXaxis()->SetLabelSize(20);
    hist_Delta1From2_EvtByEvt_p[i]->GetXaxis()->SetLabelSize(20);

    hist1_p[i]->GetYaxis()->SetLabelFont(43);
    hist2_p[i]->GetYaxis()->SetLabelFont(43);
    hist_Rat1Over2_p[i]->GetYaxis()->SetLabelFont(43);
    hist_Delta1From2_EvtByEvt_p[i]->GetYaxis()->SetLabelFont(43);
    hist1_p[i]->GetYaxis()->SetLabelSize(20);
    hist2_p[i]->GetYaxis()->SetLabelSize(20);
    hist_Rat1Over2_p[i]->GetYaxis()->SetLabelSize(20);
    hist_Delta1From2_EvtByEvt_p[i]->GetYaxis()->SetLabelSize(20);

    hist1_p[i]->GetYaxis()->SetTitleOffset(hist1_p[i]->GetYaxis()->GetTitleOffset()*3.);
    hist_Rat1Over2_p[i]->GetYaxis()->SetTitleOffset(hist1_p[i]->GetYaxis()->GetTitleOffset());
    hist_Rat1Over2_p[i]->GetXaxis()->SetTitleOffset(hist_Rat1Over2_p[i]->GetXaxis()->GetTitleOffset()*3.);

    hist1_p[i]->DrawCopy("HIST E1");
    hist2_p[i]->DrawCopy("SAME *HIST E1");

    pad2_p->cd();
    hist_Rat1Over2_p[i]->SetMaximum(1.3);
    hist_Rat1Over2_p[i]->SetMinimum(0.7);
    hist_Rat1Over2_p[i]->DrawCopy("P E1");

    pad3_p->cd();
    hist_Delta1From2_EvtByEvt_p[i]->DrawCopy("HIST E1");

    std::string pdfStr = listOfCompBranches.at(i) + "_" + inFileNameCombo + "_" + std::to_string(date->GetDate()) + ".pdf";
    canv_p->SaveAs(("pdfDir/" + pdfStr).c_str());
    listOfPdf.push_back(pdfStr);
    listOfVar.push_back(listOfCompBranches.at(i));

    delete pad1_p;
    delete pad2_p;
    delete pad3_p;

    delete canv_p;

    delete hist1_p[i];
    delete hist2_p[i];
    delete hist_Rat1Over2_p[i];

    delete hist_Delta1From2_EvtByEvt_p[i];
  }

  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(Int_t jI = 0; jI < nJetTrees; ++jI){
    const std::string jetStr = listOfJetTrees1.at(jI);
    
    for(unsigned int i = 0; i < listOfJetTreeBranches.at(jI).size(); ++i){
      hist1_Jet_p.at(jI).at(i)->Write("", TObject::kOverwrite);
      hist2_Jet_p.at(jI).at(i)->Write("", TObject::kOverwrite);
      
      hist1_Jet_p.at(jI).at(i)->Sumw2();
      hist2_Jet_p.at(jI).at(i)->Sumw2();
      hist_Rat1Over2_Jet_p.at(jI).at(i)->Sumw2();
      
      hist_Rat1Over2_Jet_p.at(jI).at(i)->Divide(hist1_Jet_p.at(jI).at(i), hist2_Jet_p.at(jI).at(i));
      hist_Rat1Over2_Jet_p.at(jI).at(i)->Write("", TObject::kOverwrite);
      
      hist_Delta1From2_EvtByEvt_Jet_p.at(jI).at(i)->Write("", TObject::kOverwrite);
      
      TCanvas* canv_p = new TCanvas("canv_c", "canv_c", 1000, 500);
      canv_p->SetTopMargin(0.01);
      canv_p->SetRightMargin(0.01);
      canv_p->SetLeftMargin(0.01);
      canv_p->SetBottomMargin(0.01);
      
      TPad* pad1_p = new TPad("pad1", "pad1", 0.0, splitPoint, 0.5, 1.0);
      pad1_p->Draw();
      pad1_p->SetTopMargin(0.01);
      pad1_p->SetRightMargin(0.01);
      pad1_p->SetBottomMargin(0.01);
      pad1_p->SetLeftMargin(pad1_p->GetLeftMargin()*1.3);
      
      TPad* pad2_p = new TPad("pad2", "pad2", 0.0, 0.0, 0.5, splitPoint);
      pad2_p->Draw();
      pad2_p->SetTopMargin(0.01);
      pad2_p->SetRightMargin(0.01);
      pad2_p->SetBottomMargin(pad1_p->GetLeftMargin()*1./splitPoint);
      pad2_p->SetLeftMargin(pad1_p->GetLeftMargin());
      
      TPad* pad3_p = new TPad("pad3", "pad3", 0.5, 0.0, 1.0, 1.0);
      pad3_p->Draw();
      pad3_p->SetTopMargin(0.01);
      pad3_p->SetRightMargin(0.01);
      pad3_p->SetLeftMargin(pad1_p->GetLeftMargin());
      pad3_p->SetBottomMargin(pad1_p->GetLeftMargin());
      
      pad1_p->cd();
      
      hist1_Jet_p.at(jI).at(i)->GetXaxis()->SetTitleFont(43);
      hist2_Jet_p.at(jI).at(i)->GetXaxis()->SetTitleFont(43);
      hist_Rat1Over2_Jet_p.at(jI).at(i)->GetXaxis()->SetTitleFont(43);
      hist_Delta1From2_EvtByEvt_Jet_p.at(jI).at(i)->GetXaxis()->SetTitleFont(43);
      hist1_Jet_p.at(jI).at(i)->GetXaxis()->SetTitleSize(20);
      hist2_Jet_p.at(jI).at(i)->GetXaxis()->SetTitleSize(20);
      hist_Rat1Over2_Jet_p.at(jI).at(i)->GetXaxis()->SetTitleSize(20);
      hist_Delta1From2_EvtByEvt_Jet_p.at(jI).at(i)->GetXaxis()->SetTitleSize(20);
      
      hist1_Jet_p.at(jI).at(i)->GetYaxis()->SetTitleFont(43);
      hist2_Jet_p.at(jI).at(i)->GetYaxis()->SetTitleFont(43);
      hist_Rat1Over2_Jet_p.at(jI).at(i)->GetYaxis()->SetTitleFont(43);
      hist_Delta1From2_EvtByEvt_Jet_p.at(jI).at(i)->GetYaxis()->SetTitleFont(43);
      hist1_Jet_p.at(jI).at(i)->GetYaxis()->SetTitleSize(20);
      hist2_Jet_p.at(jI).at(i)->GetYaxis()->SetTitleSize(20);
      hist_Rat1Over2_Jet_p.at(jI).at(i)->GetYaxis()->SetTitleSize(20);
      hist_Delta1From2_EvtByEvt_Jet_p.at(jI).at(i)->GetYaxis()->SetTitleSize(20);
      
      
      hist1_Jet_p.at(jI).at(i)->GetXaxis()->SetLabelFont(43);
      hist2_Jet_p.at(jI).at(i)->GetXaxis()->SetLabelFont(43);
      hist_Rat1Over2_Jet_p.at(jI).at(i)->GetXaxis()->SetLabelFont(43);
      hist_Delta1From2_EvtByEvt_Jet_p.at(jI).at(i)->GetXaxis()->SetLabelFont(43);
      hist1_Jet_p.at(jI).at(i)->GetXaxis()->SetLabelSize(20);
      hist2_Jet_p.at(jI).at(i)->GetXaxis()->SetLabelSize(20);
      hist_Rat1Over2_Jet_p.at(jI).at(i)->GetXaxis()->SetLabelSize(20);
      hist_Delta1From2_EvtByEvt_Jet_p.at(jI).at(i)->GetXaxis()->SetLabelSize(20);
      
      hist1_Jet_p.at(jI).at(i)->GetYaxis()->SetLabelFont(43);
      hist2_Jet_p.at(jI).at(i)->GetYaxis()->SetLabelFont(43);
      hist_Rat1Over2_Jet_p.at(jI).at(i)->GetYaxis()->SetLabelFont(43);
      hist_Delta1From2_EvtByEvt_Jet_p.at(jI).at(i)->GetYaxis()->SetLabelFont(43);
      hist1_Jet_p.at(jI).at(i)->GetYaxis()->SetLabelSize(20);
      hist2_Jet_p.at(jI).at(i)->GetYaxis()->SetLabelSize(20);
      hist_Rat1Over2_Jet_p.at(jI).at(i)->GetYaxis()->SetLabelSize(20);
      hist_Delta1From2_EvtByEvt_Jet_p.at(jI).at(i)->GetYaxis()->SetLabelSize(20);
      
      hist1_Jet_p.at(jI).at(i)->GetYaxis()->SetTitleOffset(hist1_Jet_p.at(jI).at(i)->GetYaxis()->GetTitleOffset()*3.);
      hist_Rat1Over2_Jet_p.at(jI).at(i)->GetYaxis()->SetTitleOffset(hist1_Jet_p.at(jI).at(i)->GetYaxis()->GetTitleOffset());
      hist_Rat1Over2_Jet_p.at(jI).at(i)->GetXaxis()->SetTitleOffset(hist_Rat1Over2_Jet_p.at(jI).at(i)->GetXaxis()->GetTitleOffset()*3.);
      
      hist1_Jet_p.at(jI).at(i)->DrawCopy("HIST E1");
      hist2_Jet_p.at(jI).at(i)->DrawCopy("SAME *HIST E1");
      
      pad2_p->cd();
      hist_Rat1Over2_Jet_p.at(jI).at(i)->SetMaximum(1.3);
      hist_Rat1Over2_Jet_p.at(jI).at(i)->SetMinimum(0.7);
      hist_Rat1Over2_Jet_p.at(jI).at(i)->DrawCopy("P E1");
      
      pad3_p->cd();
      hist_Delta1From2_EvtByEvt_Jet_p.at(jI).at(i)->DrawCopy("HIST E1");
      
      std::string pdfStr = listOfJetTreeBranches.at(jI).at(i) + "_" + jetStr+ "_" + inFileNameCombo + "_" + std::to_string(date->GetDate()) + ".pdf";
      canv_p->SaveAs(("pdfDir/" + pdfStr).c_str());
      listOfPdf.push_back(pdfStr);
      listOfVar.push_back((jetStr + "_" + listOfJetTreeBranches.at(jI).at(i)).c_str());
      
      delete pad1_p;
      delete pad2_p;
      delete pad3_p;
      
      delete canv_p;

      delete hist1_Jet_p.at(jI).at(i);
      delete hist2_Jet_p.at(jI).at(i);
      delete hist_Rat1Over2_Jet_p.at(jI).at(i);
      
      delete hist_Delta1From2_EvtByEvt_Jet_p.at(jI).at(i);
    }
  }

  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  TNamed nameFile1("nameFile1", inFileName1.c_str());
  TNamed nameFile2("nameFile2", inFileName2.c_str());

  nameFile1.Write("", TObject::kOverwrite);
  nameFile2.Write("", TObject::kOverwrite);

  outFile_p->Close();
  delete outFile_p;

  std::string inFile1TexStr = inFileName1;
  std::string inFile2TexStr = inFileName2;

  std::string tempStr;
  while(inFile1TexStr.find("_") != std::string::npos){
    tempStr = tempStr + inFile1TexStr.substr(0, inFile1TexStr.find("_"));
    tempStr = tempStr + "\\_";
    inFile1TexStr.replace(0, inFile1TexStr.find("_")+1, "");
  }
  tempStr = tempStr + inFile1TexStr;
  inFile1TexStr = tempStr;

  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  tempStr = "";
  while(inFile2TexStr.find("_") != std::string::npos){
    tempStr = tempStr + inFile2TexStr.substr(0, inFile2TexStr.find("_"));
    tempStr = tempStr + "\\_";
    inFile2TexStr.replace(0, inFile2TexStr.find("_")+1, "");
  }
  tempStr = tempStr + inFile2TexStr;
  inFile2TexStr = tempStr;


  const std::string outTexFileName = "pdfDir/" + inFileNameCombo + "_" + std::to_string(date->GetDate()) + ".tex";
  std::ofstream fileTex(outTexFileName.c_str());

  fileTex << "\\RequirePackage{xspace}" << std::endl;
  fileTex << "\\RequirePackage{amsmath}" << std::endl;

  fileTex << std::endl;

  fileTex << "\\documentclass[xcolor=dvipsnames]{beamer}" << std::endl;
  fileTex << "\\usetheme{Warsaw}" << std::endl;
  fileTex << "\\setbeamercolor{structure}{fg=NavyBlue!90!NavyBlue}" << std::endl;
  fileTex << "\\setbeamercolor{footlinecolor}{fg=white,bg=lightgray}" << std::endl;
  fileTex << "\\newcommand{\\pt}{\\ensuremath{p_{\\mathrm{T}}}\\xspace}" << std::endl;
  fileTex << "\\setbeamersize{text margin left=5pt,text margin right=5pt}" << std::endl;

  fileTex << std::endl;

  fileTex << "\\setbeamertemplate{frametitle}" << std::endl;
  fileTex << "{" << std::endl;
  fileTex << "  \\nointerlineskip" << std::endl;
  fileTex << "  \\begin{beamercolorbox}[sep=0.3cm, ht=1.8em, wd=\\paperwidth]{frametitle}" << std::endl;
  fileTex << "    \\vbox{}\\vskip-2ex%" << std::endl;
  fileTex << "    \\strut\\insertframetitle\\strut" << std::endl;
  fileTex << "    \\vskip-0.8ex%" << std::endl;
  fileTex << "  \\end{beamercolorbox}" << std::endl;
  fileTex << "}" << std::endl;

  fileTex << std::endl;

  fileTex << "\\setbeamertemplate{footline}{%" << std::endl;
  fileTex << "  \\begin{beamercolorbox}[sep=.8em,wd=\\paperwidth,leftskip=0.5cm,rightskip=0.5cm]{footlinecolor}" << std::endl;
  fileTex << "    \\hspace{0.3cm}%" << std::endl;
  fileTex << "    \\hfill\\insertauthor \\hfill\\insertpagenumber" << std::endl;
  fileTex << "  \\end{beamercolorbox}%" << std::endl;
  fileTex << "}" << std::endl;

  fileTex << std::endl;
  
  fileTex << "\\setbeamertemplate{navigation symbols}{}" << std::endl;
  fileTex << "\\setbeamertemplate{itemize item}[circle]" << std::endl;
  fileTex << "\\setbeamertemplate{itemize subitem}[circle]" << std::endl;
  fileTex << "\\setbeamertemplate{itemize subsubitem}[circle]" << std::endl;
  fileTex << "\\setbeamercolor{itemize item}{fg=black}" << std::endl;
  fileTex << "\\setbeamercolor{itemize subitem}{fg=black}" << std::endl;
  fileTex << "\\setbeamercolor{itemize subsubitem}{fg=black}" << std::endl;

  fileTex << std::endl;

  fileTex << "\\definecolor{links}{HTML}{00BFFF}" << std::endl;
  fileTex << "\\hypersetup{colorlinks,linkcolor=,urlcolor=links}" << std::endl;

  fileTex << std::endl;

  fileTex << "\\author[CM]{Placeholder}" << std::endl;
  fileTex << "\\begin{document}" << std::endl;

  fileTex << std::endl;

  fileTex << "\\begin{frame}" << std::endl;
  fileTex << "\\frametitle{\\centerline{Sample Validation (" << date->GetYear() << "." << date->GetMonth() << "." << date->GetDay() << ")}}" << std::endl;
  fileTex << " \\begin{itemize}" << std::endl;
  fileTex << "  \\fontsize{8}{8}\\selectfont" << std::endl;
  fileTex << "  \\item{" << inFile1TexStr << "}" << std::endl;
  fileTex << "  \\item{" << inFile2TexStr << "}" << std::endl;
  fileTex << " \\end{itemize}" << std::endl;
  fileTex << "\\end{frame}" << std::endl;

  for(unsigned int i = 0; i < listOfPdf.size(); ++i){
    std::string varStr = listOfVar.at(i);
    std::string newVarStr;
    while(varStr.find("_") != std::string::npos){
      newVarStr = newVarStr + varStr.substr(0, varStr.find("_"));
      newVarStr = newVarStr + "\\_";
      varStr.replace(0, varStr.find("_")+1, "");
    }
    newVarStr = newVarStr + varStr;

    fileTex << std::endl;
    fileTex << "\\begin{frame}" << std::endl;
    fileTex << "\\frametitle{\\centerline{" << newVarStr << "}}" << std::endl;
    fileTex << "\\begin{center}" << std::endl;
    fileTex << "\\includegraphics[width=0.8\\textwidth]{" << listOfPdf.at(i) << "}" << std::endl;
    fileTex << "\\end{center}" << std::endl;
    fileTex << "\\begin{itemize}" << std::endl;
    fileTex << "\\fontsize{8}{8}\\selectfont" << std::endl;
    fileTex << "\\item{" << newVarStr << "}" << std::endl;
    fileTex << "\\end{itemize}" << std::endl;
    fileTex << "\\end{frame}" << std::endl;
  }

  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  
  fileTex << std::endl;
  fileTex << "\\end{document}" << std::endl;

  fileTex.close();

  delete date;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3 && argc != 4){
    std::cout << "Usage ./doComparison.exe <inFileName1> <inFileName2> <outFileName-optional>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 3) retVal += doComparison(argv[1], argv[2]);
  else if(argc == 4) retVal += doComparison(argv[1], argv[2], argv[3]);
  return retVal;
}
