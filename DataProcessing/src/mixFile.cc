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

#include <time.h>
#include <sys/time.h>

//local dependencies
#include "DataProcessing/include/doLocalDebug.h"
#include "DataProcessing/include/checkMakeDir.h"
#include "DataProcessing/include/particleData.h"
#include "DataProcessing/include/eventData.h"
#include "DataProcessing/include/jetData.h"
#include "DataProcessing/include/boostedEvtData.h"
#include "DataProcessing/include/thrustTools.h"
#include "DataProcessing/include/mixTools.h"
#include "DataProcessing/include/mixMap.h"
#include "DataProcessing/include/boostTools.h"
#include "DataProcessing/include/sphericityTools.h"
#include "DataProcessing/include/returnRootFileContentsList.h"
#include "DataProcessing/include/removeVectorDuplicates.h"

double get_wall_time(){
  struct timeval time;
  if (gettimeofday(&time,NULL)){
    //  Handle error
    return 0;
  }
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

int makeMixFile(std::string inputFile, const bool isMC, std::string outputFile = "", const int nEvts2Mix = 3, const int maxMult = 35){
  double startTime = get_wall_time();
  double getTime = 0;

  //open file and get main tree
  //if given a .aleph file, look for .root file in local directory instead (for scan.cc implementation)
  if(inputFile.find(".aleph") != std::string::npos){
    while(inputFile.find("/") != std::string::npos){inputFile.replace(0, inputFile.find("/")+1,"");}
    size_t lastindex = inputFile.find_last_of(".");
    inputFile = inputFile.substr(0, lastindex)+".root";
  }
  TFile * input = TFile::Open(inputFile.c_str(),"read");
  TTree * inTree1_p = (TTree*)input->Get("t");
  TTree* inMCTree1_p = NULL;
  if(isMC) inMCTree1_p = (TTree*)input->Get("tgen");
  
  //compile list of branches and set them
  //std::vector<std::string> listOfBranches;
  //TObjArray* list1_p = (TObjArray*)inTree1_p->GetListOfBranches();
  //for(Int_t i = 0; i < list1_p->GetEntries(); ++i)  listOfBranches.push_back(list1_p->At(i)->GetName());
  
  
  //particle and event info
  inTree1_p->SetBranchStatus("*",0);
  if(isMC) inMCTree1_p->SetBranchStatus("*",0);

  particleData pdataSig;
  eventData edataSig;  
  particleData pdataMCSig;
  eventData edataMCSig;  

  std::vector<std::string> listOfBranches2 = {"passesAll", "nChargedHadronsHP","passesLEP1TwoPC"};
  edataSig.SetStatusAndAddressRead(inTree1_p, listOfBranches2); 

  //set up mixing multiplicity map
  MixMap m = MixMap(maxMult);
  for(int i = 0; i<inTree1_p->GetEntries(); i++){
    inTree1_p->GetEntry(i);
    if(!(edataSig.passesLEP1TwoPC)) continue;
    m.addElement(edataSig.nChargedHadronsHP,i);
  }
  
  std::vector<std::string> listOfBranches = {"nParticle", "process", "Energy", "px", "py", "pz", "pt", "pmag", "rap", "eta", "theta", "phi", "highPurity", "pwflag", "passesAll","passesLEP1TwoPC", "nChargedHadronsHP", "Thrust", "TTheta", "TPhi", "Thrust_charged", "TTheta_charged", "TPhi_charged"};
  pdataSig.SetStatusAndAddressRead(inTree1_p, listOfBranches);
  edataSig.SetStatusAndAddressRead(inTree1_p, listOfBranches);

  if(isMC){
    pdataMCSig.SetStatusAndAddressRead(inMCTree1_p, listOfBranches);
    edataMCSig.SetStatusAndAddressRead(inMCTree1_p, listOfBranches);
  }

  //wta branches
  std::vector<std::string> listOfBoostedTrees1 = returnRootFileContentsList(input, "TTree", "Boosted");
  removeVectorDuplicates(&listOfBoostedTrees1);
  std::vector< std::vector<std::string> > listOfBoostedTreeBranches;
  for(unsigned int jI = 0; jI < listOfBoostedTrees1.size(); ++jI){
   // TTree* tempTree_p = (TTree*)input->Get(listOfBoostedTrees1.at(jI).c_str());
   // TObjArray* tempList = (TObjArray*)tempTree_p->GetListOfBranches();
   // std::vector<std::string> tempV;
    std::vector<std::string> tempV = {"WTAAxis_Theta","WTAAxis_ThetaPerp","WTAAxis_Phi","boostx","boosty","boostz"};

   // for(Int_t i = 0; i < tempList->GetEntries(); ++i) tempV.push_back(tempList->At(i)->GetName());

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
    size_t lastindex = outputFile.find_last_of(".");
    output = TFile::Open(Form("%s_Mix.root",(outputFile.substr(0, lastindex)).c_str()),"recreate");
  }
  
  //particle info
  particleData pdataMix;
  TTree * outTree1_p = new TTree("t","t");
  outTree1_p->SetDirectory(output);
  pdataMix.SetBranchWrite(outTree1_p);

  particleData pdataMCMix;
  TTree * outMCTree1_p = NULL;
  if(isMC){
    outMCTree1_p = new TTree("tgen","tgen");
    outMCTree1_p->SetDirectory(output);
    pdataMCMix.SetBranchWrite(outMCTree1_p);
  }

  //evt info
  eventData edataMix;
  eventData edataMCMix;
  edataMix.SetBranchWrite(outTree1_p);
  if(isMC) edataMCMix.SetBranchWrite(outMCTree1_p);

  //wta trees
  boostedEvtData bDataMix[nBoostedTrees];
  TTree * outTree1_b[nBoostedTrees];
  for(Int_t jI = 0; jI < nBoostedTrees; ++jI){
    outTree1_b[jI] = new TTree(listOfBoostedTrees1.at(jI).c_str(), listOfBoostedTrees1.at(jI).c_str());
    outTree1_b[jI]->SetDirectory(output);
    bDataMix[jI].SetBranchWrite(outTree1_b[jI]);
  }

  const Int_t printInterval = 20;
  const Int_t printNEntries = TMath::Max(1,(int)(inTree1_p->GetEntries()/printInterval));

  double startTimeLoop = get_wall_time();

  //loop through input file for the signal events
  for(int i = 0; i<inTree1_p->GetEntries(); i++){
    if(i%printNEntries == 0) std::cout << "Mixing entry: " << i << "/" << inTree1_p->GetEntries() << std::endl;
    //for testing
    //if(i>1000) break;
    double tempTime = get_wall_time(); 
    inTree1_p->GetEntry(i);
    if(isMC) inMCTree1_p->GetEntry(i);
    for(Int_t jI = 0; jI < nBoostedTrees; ++jI) bTree1_p[jI]->GetEntry(i);
    getTime += get_wall_time()-tempTime; 

    //set temp variables for signal event if needed
    //get thrust axes
    TVector3 thrustAxis, thrustAxis_ch;
    thrustAxis.SetMagThetaPhi(edataSig.Thrust, edataSig.TTheta, edataSig.TPhi);  
    thrustAxis_ch.SetMagThetaPhi(edataSig.Thrust_charged, edataSig.TTheta_charged, edataSig.TPhi_charged);  
   
    int signalMultiplicity = edataSig.nChargedHadronsHP; 
    int signalProcess = pdataSig.process;
    int signalEnergy = (int)(pdataSig.Energy < 100);

    //get WTA and boosts
    TVector3 WTAAxis[nBoostedTrees], WTABoost[nBoostedTrees];
    for(Int_t jI = 0; jI < nBoostedTrees; ++jI){
      WTAAxis[jI].SetMagThetaPhi(1,bData1[jI].WTAAxis_Theta,bData1[jI].WTAAxis_Phi);
      WTABoost[jI].SetXYZ(bData1[jI].boostx, bData1[jI].boosty, bData1[jI].boostz);
    }

    //setup loop over other events in file to get mixed events
    resetMixEvt(&pdataMix);
    if(isMC){
      resetMixEvt(&pdataMCMix);
    }
    for(Int_t jI = 0; jI < nBoostedTrees; ++jI) resetMixEvtBoosted(&bDataMix[jI]);
   
    m.setCurrentElement(signalMultiplicity,i);//makes sure you find 'close' events in the file to speed things up
  
    int mixedEventsFound = 0;
    while(mixedEventsFound < nEvts2Mix){
      int j = m.getNextElement((unsigned char)signalMultiplicity, i);
      tempTime = get_wall_time();
      inTree1_p->GetEntry(j);
      if(isMC) inMCTree1_p->GetEntry(j);

      getTime += get_wall_time()-tempTime;     
  
      //if we check the entire file and haven't found enough mixed events, give up and mix it with itself.  Warning if it passes evt sel
      if(i==j){
        if(edataSig.passesLEP1TwoPC) std::cout << "Warning! Only " << mixedEventsFound << "/" << nEvts2Mix << " found in file for event index " <<  i << "!  Giving up with less than the requested number of mixed events..." << std::endl;
        appendMixEvt(&edataMix, &pdataMix, &pdataSig, thrustAxis, thrustAxis_ch, (float)mixedEventsFound);
        if(isMC) appendMixEvt(&edataMCMix, &pdataMCMix, &pdataMCSig, thrustAxis, thrustAxis_ch, (float)mixedEventsFound);
        for(Int_t jI = 0; jI < nBoostedTrees; ++jI) appendMixEvtBoosted(&bDataMix[jI], &pdataSig, WTAAxis[jI], WTABoost[jI], (float)mixedEventsFound);
        break;
      }
  

      //make sure mixed event is a good one (using Sig tree here because it is the input tree)
      bool processesMatch = (signalProcess == pdataSig.process);
      bool energiesMatch = (signalEnergy == (int)(pdataSig.Energy < 100));

      if(processesMatch && energiesMatch) mixedEventsFound++;
      else continue;

      //particle loop is in here
      appendMixEvt(&edataMix, &pdataMix, &pdataSig, thrustAxis, thrustAxis_ch, (float)mixedEventsFound);
      if(isMC) appendMixEvt(&edataMCMix, &pdataMCMix, &pdataMCSig, thrustAxis, thrustAxis_ch, (float)mixedEventsFound);
      for(Int_t jI = 0; jI < nBoostedTrees; ++jI) appendMixEvtBoosted(&bDataMix[jI], &pdataSig, WTAAxis[jI], WTABoost[jI], (float)mixedEventsFound);
    }
    pdataMix.preFillClean();
    if(isMC) pdataMCMix.preFillClean();
    outTree1_p->Fill();
    if(isMC) outMCTree1_p->Fill();
    for(Int_t jI = 0; jI < nBoostedTrees; ++jI){bDataMix[jI].preFillClean(); outTree1_b[jI]->Fill();}
  }
  output->Write("", TObject::kOverwrite);
  
  for(Int_t jI = 0; jI < nBoostedTrees; ++jI) delete outTree1_b[jI];
  delete outTree1_p;
  if(isMC) delete outMCTree1_p;
  output->Close();


  std::cout << "Job Complete!" << std::endl;
  double endTime = get_wall_time();
  std::cout << "Mix wall time: " << endTime-startTime << std::endl;
  std::cout << "Looping wall time: " << endTime-startTimeLoop << std::endl;
  std::cout << "getEntry time: " << getTime << " " << getTime/(endTime-startTimeLoop)*100 << "\%" << std::endl;
  return 0;
}

