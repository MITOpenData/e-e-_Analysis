//c and cpp dependencies
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

//root dependencies
#include <TFile.h>
#include <TTree.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TNamed.h"

//local dependencies
#include "include/doLocalDebug.h"
#include "include/checkMakeDir.h"
#include "include/particleData.h"
#include "include/eventSelection.h"
#include "include/eventData.h"
#include "include/jetData.h"
#include "include/boostedEvtData.h"
#include "include/thrustTools.h"
#include "include/mixTools.h"
#include "include/boostTools.h"
#include "include/sphericityTools.h"
#include "include/returnRootFileContentsList.h"
#include "include/removeVectorDuplicates.h"

int makeMixFile(std::string inputFile, std::string outputFile = "", const int nEvts2Mix = 3){

  //open file and get main tree
  //if given a .aleph file, look for .root file in local directory instead (for scan.cc implementation)
  if(inputFile.find(".aleph") != std::string::npos){
    while(inputFile.find("/") != std::string::npos){inputFile.replace(0, inputFile.find("/")+1,"");}
    size_t lastindex = inputFile.find_last_of(".");
    inputFile = inputFile.substr(0, lastindex)+".root";
  }
  TFile * input = TFile::Open(inputFile.c_str(),"read");
  TTree * inTree1_p = (TTree*)input->Get("t");
  
  //compile list of branches and set them
  std::vector<std::string> listOfBranches;
  TObjArray* list1_p = (TObjArray*)inTree1_p->GetListOfBranches();
  for(Int_t i = 0; i < list1_p->GetEntries(); ++i)  listOfBranches.push_back(list1_p->At(i)->GetName());
  
  //particle and event info
  particleData pdataSig;
  pdataSig.SetStatusAndAddressRead(inTree1_p, listOfBranches);
  eventData edataSig;  
  edataSig.SetStatusAndAddressRead(inTree1_p, listOfBranches);

  //wta branches
  std::vector<std::string> listOfBoostedTrees1 = returnRootFileContentsList(input, "TTree", "Boosted");
  removeVectorDuplicates(&listOfBoostedTrees1);
  std::vector< std::vector<std::string> > listOfBoostedTreeBranches;
  for(unsigned int jI = 0; jI < listOfBoostedTrees1.size(); ++jI){
    TTree* tempTree_p = (TTree*)input->Get(listOfBoostedTrees1.at(jI).c_str());
    TObjArray* tempList = (TObjArray*)tempTree_p->GetListOfBranches();
    std::vector<std::string> tempV;

    for(Int_t i = 0; i < tempList->GetEntries(); ++i) tempV.push_back(tempList->At(i)->GetName());

    listOfBoostedTreeBranches.push_back(tempV);
  }
  const Int_t nBoostedTrees = listOfBoostedTrees1.size();
  boostedEvtData bData1[nBoostedTrees];
  TTree* bTree1_p[nBoostedTrees];
  for(Int_t jI = 0; jI < nBoostedTrees; ++jI){
    bTree1_p[jI] = (TTree*)input->Get(listOfBoostedTrees1.at(jI).c_str());
    bTree1_p[jI]->SetBranchStatus("*", 0);
    bData1[jI].SetStatusAndAddressRead(bTree1_p[jI], listOfBoostedTreeBranches.at(jI));
  }


  //open output file
  TFile * output;
  if(outputFile.compare("")==0){
    while(inputFile.find("/") != std::string::npos){inputFile.replace(0, inputFile.find("/")+1,"");}
    size_t lastindex = inputFile.find_last_of(".");
    output = TFile::Open(Form("%s_Mix.root",(inputFile.substr(0, lastindex)).c_str()),"recreate");
  }
  else{  
    output = TFile::Open(outputFile.c_str(),"recreate");
  }
  
  //particle info
  particleData pdataMix;
  TTree * outTree1_p = new TTree("t","t");
  outTree1_p->SetDirectory(output);
  pdataMix.SetBranchWrite(outTree1_p);

  //wta trees
  boostedEvtData bDataMix[nBoostedTrees];
  TTree * outTree1_b[nBoostedTrees];
  for(Int_t jI = 0; jI < nBoostedTrees; ++jI){
    outTree1_b[jI] = new TTree(listOfBoostedTrees1.at(jI).c_str(), listOfBoostedTrees1.at(jI).c_str());
    outTree1_b[jI]->SetDirectory(output);
    bDataMix[jI].SetBranchWrite(outTree1_b[jI]);
  }
  
 
  //loop through input file for the signal events
  for(int i = 0; i<inTree1_p->GetEntries(); i++){
    //for testing
    if(i>100) break;
    
    inTree1_p->GetEntry(i);
    for(Int_t jI = 0; jI < nBoostedTrees; ++jI) bTree1_p[jI]->GetEntry(i);
    //set temp variables for signal event if needed
    //get thrust axes
    TVector3 thrustAxis, thrustAxis_ch;
    thrustAxis.SetMagThetaPhi(1, edataSig.TTheta, edataSig.TPhi);  
    thrustAxis_ch.SetMagThetaPhi(1, edataSig.TTheta_charged, edataSig.TPhi_charged);  

    //get WTA and boosts
    TVector3 WTAAxis[nBoostedTrees], WTABoost[nBoostedTrees];
    for(Int_t jI = 0; jI < nBoostedTrees; ++jI){
      WTAAxis[jI].SetMagThetaPhi(1,bData1[jI].WTAAxis_Theta,bData1[jI].WTAAxis_Phi);
      WTABoost[jI].SetXYZ(bData1[jI].boostx, bData1[jI].boosty, bData1[jI].boostz);
    }

    //setup loop over other events in file to get mixed events
    resetMixEvt(&pdataMix);
    for(Int_t jI = 0; jI < nBoostedTrees; ++jI) resetMixEvtBoosted(&bDataMix[jI]);

    int mixedEventsFound = 0;
    for(int j = i+1; mixedEventsFound < nEvts2Mix; j++){
      //wrap around if you hit end of file
      if(j==inTree1_p->GetEntries()) j=0;

      //if we check the entire file and haven't found enough mixed events, give up 
      if(i==j){
        std::cout << "Warning! Only " << mixedEventsFound << "/" << nEvts2Mix << " found in file for event index i!  Giving up with less than the requested number of mixed events..." << std::endl;
        break;
      }
  
      inTree1_p->GetEntry(j);

      //make sure mixed event is a good one
      bool isGoodEvent = true;
      if(isGoodEvent) mixedEventsFound++;

      //particle loop is in here
      appendMixEvt(&pdataMix, &pdataSig, thrustAxis, thrustAxis_ch, (float)mixedEventsFound);
      for(Int_t jI = 0; jI < nBoostedTrees; ++jI) appendMixEvtBoosted(&bDataMix[jI], &pdataSig, WTAAxis[jI], WTABoost[jI], (float)mixedEventsFound);
    }
    outTree1_p->Fill();
    for(Int_t jI = 0; jI < nBoostedTrees; ++jI) outTree1_b[jI]->Fill();
  }
  output->Write();
  
  for(Int_t jI = 0; jI < nBoostedTrees; ++jI) delete outTree1_b[jI];
  delete outTree1_p;
  output->Close();


  std::cout << "Job Complete!" << std::endl;
  return 0;
}

