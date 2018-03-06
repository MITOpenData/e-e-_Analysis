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

#include "DataProcessing/include/particleData.h"
#include "DataProcessing/include/eventData.h"
#include "DataProcessing/include/jetData.h"
#include "DataProcessing/include/returnRootFileContentsList.h"
#include "DataProcessing/include/removeVectorDuplicates.h"
#include "DataProcessing/include/histDefUtility.h"
#include "include/handleMultipleFiles.h"
#include "include/Selection.h"
// hist_p = {hist1_p,hist2_p, ... , histN_p}
// histDelta_p[] = {hist_Delta1From2,hist_Delta1From3, ... , hist_Delta1FromN, hist_Delta2From3, ... , hist_Delta2FromN, ... , hist_DeltaN-1FromN}
// val[] = {val1, val2, ... , valN}

void doFill(std::vector< TH1F* > hist_p, std::vector< TH1F* > histDelta_p, std::vector< Float_t > val)
{
    if(hist_p.size() != val.size()){std::cout<<"Not the same number of histograms as values"<<std::endl; return;}
    if(histDelta_p.size() != (val.size()*(val.size()-1)/2)) {std::cout<<"Not the N(N-1)/2 of delta histograms"<<std::endl; return;}
    
    unsigned int loc = 0;
    for(unsigned int hI = 0; hI< hist_p.size();++hI)
    {
        hist_p.at(hI)->Fill(val.at(hI));
        for(unsigned int hII = hI+1; hII< hist_p.size();++hII)
        {
            if(val.at(hI)-val.at(hII) > histDelta_p.at(loc)->GetBinLowEdge(histDelta_p.at(loc)->GetNbinsX()+1)) histDelta_p.at(loc)->Fill(histDelta_p.at(loc)->GetBinCenter(histDelta_p.at(loc)->GetNbinsX()));
            else if(val.at(hI)-val.at(hII) < histDelta_p.at(loc)->GetBinLowEdge(1)) histDelta_p.at(loc)->Fill(histDelta_p.at(loc)->GetBinCenter(1));
            else histDelta_p.at(loc)->Fill(val.at(hI)-val.at(hII));
            ++loc;
        }
    }
    if ( (loc) != (unsigned int) histDelta_p.size() ) std::cout<<"Incorrect number of delta histograms filled in doFill"<<std::endl;
    
    return;
}

void doFill(std::vector< TH1F* > hist_p, std::vector< TH1F* > histDelta_p, std::vector< Int_t > val)
{
    if(hist_p.size() != val.size()){std::cout<<"Not the same number of histograms as values"<<std::endl; return;}
    if(histDelta_p.size() != (val.size()*(val.size()-1)/2)) {std::cout<<"Not the N(N-1)/2 of delta histograms"<<std::endl; return;}
    
    unsigned int loc = 0;
    for(unsigned int hI = 0; hI< hist_p.size();++hI)
    {
        hist_p.at(hI)->Fill(val.at(hI));
        for(unsigned int hII = hI+1; hII< hist_p.size();++hII)
        {
            if(val.at(hI)-val.at(hII) > histDelta_p.at(loc)->GetBinLowEdge(histDelta_p.at(loc)->GetNbinsX()+1)) histDelta_p.at(loc)->Fill(histDelta_p.at(loc)->GetBinCenter(histDelta_p.at(loc)->GetNbinsX()));
            else if(val.at(hI)-val.at(hII) < histDelta_p.at(loc)->GetBinLowEdge(1)) histDelta_p.at(loc)->Fill(histDelta_p.at(loc)->GetBinCenter(1));
            else histDelta_p.at(loc)->Fill(val.at(hI)-val.at(hII));
            ++loc;
        }
    }
    if ( (loc) != (unsigned int) histDelta_p.size() ) std::cout<<"Incorrect number of delta histograms filled in doFill"<<std::endl;
    
    return;
}

void doFillArr(std::vector< TH1F* > hist_p, std::vector< TH1F* > histDelta_p, std::vector<std::vector<Float_t> > val)
{
    if(hist_p.size() != val.size()){std::cout<<"Not the same number of histograms as values"<<std::endl; return;}
    if(histDelta_p.size() != (val.size()*(val.size()-1)/2)) {std::cout<<"Not the N(N-1)/2 of delta histograms"<<std::endl; return;}
    
    unsigned int loc = 0;
    for (unsigned int hI = 0; hI < hist_p.size(); ++hI)
    {
        for(unsigned int pI = 0; pI < val.at(hI).size(); ++pI) hist_p.at(hI)->Fill(val.at(hI).at(pI));
        for(unsigned int hII = hI+1; hII< hist_p.size(); ++hII)
        {
            if(val.at(hI).size() == val.at(hII).size())
            {
                for(unsigned int pI = 0; pI < val.at(hI).size(); ++pI)
                {
                    if(val.at(hI).at(pI)-val.at(hII).at(pI) > histDelta_p.at(loc)->GetBinLowEdge(histDelta_p.at(loc)->GetNbinsX()+1)) histDelta_p.at(loc)->Fill(histDelta_p.at(loc)->GetBinCenter(histDelta_p.at(loc)->GetNbinsX()));
                    else if(val.at(hI).at(pI)-val.at(hII).at(pI) < histDelta_p.at(loc)->GetBinLowEdge(1)) histDelta_p.at(loc)->Fill(histDelta_p.at(loc)->GetBinCenter(1));
                    else histDelta_p.at(loc)->Fill(val.at(hI).at(pI)-val.at(hII).at(pI));
                }
            }
            else std::cout<<"Incorrect number of entries in val array passed into doFillArr"<<std::endl;
            ++loc;
        }
    }
    if ( (loc) != (unsigned int) histDelta_p.size() ) std::cout<<"Incorrect number of delta histograms filled in doFillArr"<<std::endl;
   return;
}

void doFillArr(std::vector< TH1F* > hist_p, std::vector< TH1F* > histDelta_p, std::vector<std::vector<Int_t> > val)
{
    if(hist_p.size() != val.size()){std::cout<<"Not the same number of histograms as values"<<std::endl; return;}
    if(histDelta_p.size() != (val.size()*(val.size()-1)/2)) {std::cout<<"Not the N(N-1)/2 of delta histograms"<<std::endl; return;}
    
    unsigned int loc = 0;
    for (unsigned int hI = 0; hI < hist_p.size(); ++hI)
    {
        for(unsigned int pI = 0; pI < val.at(hI).size(); ++pI) hist_p.at(hI)->Fill(val.at(hI).at(pI));
        for(unsigned int hII = hI+1; hII< hist_p.size(); ++hII)
        {
            if(val.at(hI).size() == val.at(hII).size())
            {
                for(unsigned int pI = 0; pI < val.at(hI).size(); ++pI)
                {
                    if(val.at(hI).at(pI)-val.at(hII).at(pI) > histDelta_p.at(loc)->GetBinLowEdge(histDelta_p.at(loc)->GetNbinsX()+1)) histDelta_p.at(loc)->Fill(histDelta_p.at(loc)->GetBinCenter(histDelta_p.at(loc)->GetNbinsX()));
                    else if(val.at(hI).at(pI)-val.at(hII).at(pI) < histDelta_p.at(loc)->GetBinLowEdge(1)) histDelta_p.at(loc)->Fill(histDelta_p.at(loc)->GetBinCenter(1));
                    else histDelta_p.at(loc)->Fill(val.at(hI).at(pI)-val.at(hII).at(pI));
                }
            }
            else std::cout<<"Incorrect number of entries in val array passed into doFillArr"<<std::endl;
            ++loc;
        }
    }
    if ( (loc) != (unsigned int) histDelta_p.size() ) std::cout<<"Incorrect number of delta histograms filled in doFillArr"<<std::endl;
    return;
}

void formatTH1F(TH1F * h)
{
    h->SetStats(0);
    h->GetXaxis()->SetTitleFont(43);
    h->GetXaxis()->SetTitleSize(20);
    h->GetYaxis()->SetTitleFont(43);
    h->GetYaxis()->SetTitleSize(20);
    h->GetXaxis()->SetLabelFont(43);
    h->GetXaxis()->SetLabelSize(20);
    h->GetYaxis()->SetLabelFont(43);
    h->GetYaxis()->SetLabelSize(20);
    h->GetYaxis()->SetTitleOffset(h->GetYaxis()->GetTitleOffset()*3.);
    h->GetXaxis()->SetTitleOffset(h->GetXaxis()->GetTitleOffset()*3.);
}

int plotDQC(const std::string inFileName, std::string outFileName = "", int doECut = 0, int drawHist = 0)
{
    
  ///// check if there are multiple files /////
  std::vector<std::string> fileList = handleMultipleFiles(inFileName);
  if(fileList.size() == 0){
      std::cout << "Given input \'" << inFileName << "\' doesn't produce valid input. return 1" << std::endl;
      return 1;
  }
    
  //////////////////// CREATE OUTFILE NAME ////////////////
    // create a outFileName if the one supplied is nothing
    // load the date
  TDatime* date = new TDatime();
  if(outFileName.size() == 0) outFileName = inFileName;
  if(outFileName.find(".root") != std::string::npos){outFileName.replace(outFileName.find(".root"), 5, "");}
  outFileName = outFileName + "_PLOTS_" + std::to_string(date->GetDate()) + ".root";
    
    
  //////////////////// LOADING ALL OF THE BRANCHES TO BE PLOTTED ////////////////
  // these are the master list of branches to plot and associated ranges
    
  std::cout << "Determining branches to be plotted..." << std::endl;
  std::vector<std::string> listOfCompBranches;
  std::vector<double> branchMins;
  std::vector<double> branchMaxs;

  std::vector<particleData> pData;
  std::vector<eventData> eData;
    
  std::vector< TFile* > inFile;
  std::vector< TTree* > inTree;
    
  for (unsigned int fI = 0; fI < fileList.size() ; ++fI)
  {
      particleData pData_p;
      pData.push_back(pData_p);
      eventData eData_p;
      eData.push_back(eData_p);
      TFile* inFile_p;
      inFile.push_back(inFile_p);
      TTree* inTree_p;
      inTree.push_back(inTree_p);
  }
    

  ///////// LOADING THE VARIABLES TO PLOT FROM THE FILES /////////
    
  /////// LOADING THE INITIAL BRANCHES/RANGES FROM FILE 1 ///////
  for(unsigned int i = 0; i<fileList.size();++i)std::cout<<fileList.at(i).c_str()<<std::endl;
  
  inFile.at(0) = new TFile(fileList.at(0).c_str(), "READ");
  inTree.at(0) = (TTree*)inFile.at(0)->Get("t");
  // turn on only the runNo and EventNo branches
  //inTree.at(0)->SetBranchStatus("*", 0);
  //inTree.at(0)->SetBranchStatus("RunNo", 1);
  //inTree.at(0)->SetBranchStatus("EventNo", 1);
  inTree.at(0)->SetBranchAddress("RunNo", &pData.at(0).RunNo);
  inTree.at(0)->SetBranchAddress("EventNo", &pData.at(0).EventNo);
  
  // load the master branch list from the first tree from file 1
  std::cout<<"Initializing Mins and Maxs..."<<std::endl;
  TObjArray* list1_p = (TObjArray*)inTree.at(0)->GetListOfBranches();
  for(Int_t i = 0; i < list1_p->GetEntries(); ++i){
    listOfCompBranches.push_back(list1_p->At(i)->GetName());
    // load in the minimum and maximums from the branches
    branchMins.push_back(inTree.at(0)->GetMinimum(list1_p->At(i)->GetName()));
    branchMaxs.push_back(inTree.at(0)->GetMaximum(list1_p->At(i)->GetName()));
    //std::cout << listOfCompBranches.at(i)<<" "<< branchMins.at(i) << " "<< inTree.at(0)->GetMinimum(list1_p->At(i)->GetName()) << ", " << branchMaxs.at(i) << " "<<inTree.at(0)->GetMaximum(list1_p->At(i)->GetName()) <<std::endl;
  }
  
  // get the list of jet trees in the file
  
  std::vector<std::string> listOfJetTrees1 = returnRootFileContentsList(inFile.at(0), "TTree", "JetTree");
  removeVectorDuplicates(&listOfJetTrees1);
  std::vector< std::vector<std::string> > listOfJetTreeBranches;
  // for the minimum and maximum of the jet trees
  std::vector< std::vector<double> > listOfJetTreeMins;
  std::vector< std::vector<double> > listOfJetTreeMaxs;
  for(unsigned int jI = 0; jI < listOfJetTrees1.size(); ++jI){
    // load the temporary jet tree
    TTree* tempTree_p = (TTree*)inFile.at(0)->Get(listOfJetTrees1.at(jI).c_str());
    TObjArray* tempList = (TObjArray*)tempTree_p->GetListOfBranches();
    std::vector<std::string> tempV;
    std::vector<double> tempMins;
    std::vector<double> tempMaxs;
    
    // load the mins and maxs from that tree
    for(Int_t i = 0; i < tempList->GetEntries(); ++i){
      tempV.push_back(tempList->At(i)->GetName());
      tempMins.push_back(tempTree_p->GetMinimum(tempList->At(i)->GetName()));
      tempMaxs.push_back(tempTree_p->GetMaximum(tempList->At(i)->GetName()));
    }
    
    listOfJetTreeBranches.push_back(tempV);
    listOfJetTreeMins.push_back(tempMins);
    listOfJetTreeMaxs.push_back(tempMaxs);
  }
  
  // clean up
  inFile.at(0)->Close();
  delete inFile.at(0);
  
  ////////////////////// DONE LOADING FROM FILE 1 //////////////////////
  
  
  ////////////////////// LOOP OVER OTHER FILES //////////////////////
  
  for (unsigned int fI = 1; fI < fileList.size(); ++fI)
    {
      inFile.at(fI) = new TFile(fileList.at(fI).c_str(), "READ");
      inTree.at(fI) = (TTree*)inFile.at(fI)->Get("t");
      TObjArray* list2_p = (TObjArray*)inTree.at(fI)->GetListOfBranches();
      
      inTree.at(fI)->SetBranchStatus("*", 0);
      inTree.at(fI)->SetBranchStatus("RunNo", 1);
      inTree.at(fI)->SetBranchStatus("EventNo", 1);
      inTree.at(fI)->SetBranchAddress("RunNo", &pData.at(fI).RunNo);
      inTree.at(fI)->SetBranchAddress("EventNo", &pData.at(fI).EventNo);
      
      unsigned int pos = 0;
      while(pos < listOfCompBranches.size()){
          std::string tempStr = listOfCompBranches.at(pos);
          
          bool isInArr = false;
          for(Int_t i = 0; i < list2_p->GetEntries(); ++i){
              
              // check if the branch is in both files
              std::string tempStr2 = list2_p->At(i)->GetName();
              
              if(tempStr.size() == tempStr2.size() && tempStr.find(tempStr2) != std::string::npos){
                  isInArr = true;
                  break;
              }
          }
          
          // if the branch is not in both branches then delete it from the list of things to be plotted
          if(!isInArr){
              listOfCompBranches.erase(listOfCompBranches.begin()+pos);
              branchMins.erase(branchMins.begin()+pos);
              branchMaxs.erase(branchMaxs.begin()+pos);
          }
          else{
              // if not then update the min and max for that variable
              Double_t tempMin = inTree.at(fI)->GetMinimum(listOfCompBranches.at(pos).c_str());
              Double_t tempMax = inTree.at(fI)->GetMaximum(listOfCompBranches.at(pos).c_str());
              
              if(branchMins.at(pos) > tempMin) branchMins.at(pos) = tempMin;
              if(branchMaxs.at(pos) < tempMax) branchMaxs.at(pos) = tempMax;
              ++pos;
          }
      }
      
      // load the list of jet trees from file
      std::vector<std::string> listOfJetTrees2 = returnRootFileContentsList(inFile.at(fI), "TTree", "JetTree");
      removeVectorDuplicates(&listOfJetTrees2);
      
      unsigned int lI = 0;
      while(lI < listOfJetTrees1.size()){
          bool isGood = false;
          
          // check if the jet tree is in both file 1 and file 2
          for(unsigned int lI2 = 0; lI2 < listOfJetTrees2.size(); ++lI2){
              if(listOfJetTrees1.at(lI).size() == listOfJetTrees2.at(lI2).size() && listOfJetTrees1.at(lI).find(listOfJetTrees2.at(lI2)) != std::string::npos){
                  isGood = true;
                  break;
              }
          }
          
          // if the jet tree is NOT in both files then remove it
          if(!isGood){
              listOfJetTrees1.erase(listOfJetTrees1.begin()+lI);
              listOfJetTreeBranches.erase(listOfJetTreeBranches.begin()+lI);
              listOfJetTreeMins.erase(listOfJetTreeMins.begin()+lI);
              listOfJetTreeMaxs.erase(listOfJetTreeMaxs.begin()+lI);
          }
          // if the jet tree is in both files then we go into each tree and check that all of the branches are the same
          else{
              TTree* tempTree_p = (TTree*)inFile.at(fI)->Get(listOfJetTrees1.at(lI).c_str());
              TObjArray* tempList_p = (TObjArray*)tempTree_p->GetListOfBranches();
              
              unsigned int j = 0;
              while(j < listOfJetTreeBranches.at(lI).size()){
                  bool branchIsGood = false;
                  
                  // check if branch is in both lists
                  for(Int_t i = 0; i < tempList_p->GetEntries(); ++i){
                      std::string tempBranchName = tempList_p->At(i)->GetName();
                      if(tempBranchName.size() == listOfJetTreeBranches.at(lI).at(j).size() && tempBranchName.find(listOfJetTreeBranches.at(lI).at(j)) != std::string::npos){
                          branchIsGood = true;
                          break;
                      }
                  }
                  
                  // if branch isnt in the second one then remove it
                  if(!branchIsGood){
                      listOfJetTreeBranches.at(lI).erase(listOfJetTreeBranches.at(lI).begin() + j);
                      listOfJetTreeMins.at(lI).erase(listOfJetTreeMins.at(lI).begin() + j);
                      listOfJetTreeMaxs.at(lI).erase(listOfJetTreeMaxs.at(lI).begin() + j);
                  }
                  // if it is then update the min and max
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
      
      inFile.at(fI)->Close();
      delete inFile.at(fI);
  }
  ////////////////////// DONE LOOPING OVER OTHER FILES //////////////////////

  /////////////////// DONE LOADING THE VARIABLES TO PLOT FROM THE FILES ///////////////////
    
 /////////// THIS IS WHERE PLOTTING BEGINS ///////////
    // create a histogram for everything found in both files
    // individual plots are created
    // ratio plots are created
    // delta plots are created as well
 /////////// /////////// /////////// /////////// ////
    
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  const Int_t nVarToComp = listOfCompBranches.size();
  const Int_t nJetTrees = listOfJetTrees1.size();
  // .at(file).at(jet tree name)
  std::vector< std::vector<jetData> > jData;
    
    
  ////////// INITIALIZING HISTOGRAMS //////////
  std::cout << "Initializing histograms..." << std::endl;
  std::vector< std::vector<TH1F*> > hist;
  std::vector< std::vector<TH1F*> > histRat;
  std::vector< std::vector<TH1F*> > histDelta;
    
  std::vector< std::vector< std::vector< TH1F* > > > hist_Jet;
  std::vector< std::vector< std::vector< TH1F* > > > histRat_Jet;
  std::vector< std::vector< std::vector< TH1F* > > > histDelta_Jet;
    
  for (unsigned int fI = 0; fI < fileList.size(); ++fI)
  {
      jetData jData1_p[nJetTrees];
      std::vector<jetData> jData1;
      jData1.assign(jData1_p,jData1_p + nJetTrees);
      jData.push_back(jData1);
      
      TH1F* temphist1_p[nVarToComp];
      std::vector<TH1F*> temphist1;
      temphist1.assign(temphist1_p, temphist1_p + nVarToComp);
      hist.push_back(temphist1);
      
      std::vector< std::vector< TH1F* > > hist1_Jet_p;
      hist_Jet.push_back(hist1_Jet_p);
      
      for(Int_t jI = 0; jI < nJetTrees; ++jI)
      {
          const int nJetBranches = listOfJetTreeBranches.at(jI).size();
          TH1F* temphist_Jet_p[nJetBranches];
          std::vector<TH1F*> temphist_Jet;
          temphist_Jet.assign(temphist_Jet_p, temphist_Jet_p + listOfJetTreeBranches.at(jI).size());
          hist_Jet.at(fI).push_back(temphist_Jet);
      }
                           
      for (unsigned int fII = fI+1; fII < fileList.size(); ++fII)
      {
          TH1F* temphist_Rat_p[nVarToComp];
          std::vector<TH1F*> temphist_Rat;
          temphist_Rat.assign(temphist_Rat_p, temphist_Rat_p + nVarToComp);
          histRat.push_back(temphist_Rat);
          
          TH1F* temphist_Delt_p[nVarToComp];
          std::vector<TH1F*> temphist_Delt;
          temphist_Delt.assign(temphist_Delt_p, temphist_Delt_p + nVarToComp);
          histDelta.push_back(temphist_Delt);
          
          std::vector< std::vector< TH1F* > > hist_Rat_Jet_p;
          histRat_Jet.push_back(hist_Rat_Jet_p);
          
          std::vector< std::vector< TH1F* > > hist_Delta_Jet_p;
          histDelta_Jet.push_back(hist_Delta_Jet_p);
          
          for(Int_t jI = 0; jI < nJetTrees; ++jI)
          {
              TH1F* temphistRat_Jet_p[nVarToComp];
              std::vector<TH1F*> temphistRat_Jet;
              temphistRat_Jet.assign(temphistRat_Jet_p, temphistRat_Jet_p + listOfJetTreeBranches.at(jI).size());
              histRat_Jet.at(fI).push_back(temphistRat_Jet);
              
              TH1F* temphistDelta_Jet_p[nVarToComp];
              std::vector<TH1F*> temphistDelta_Jet;
              temphistDelta_Jet.assign(temphistDelta_Jet_p, temphistDelta_Jet_p + listOfJetTreeBranches.at(jI).size());
              histDelta_Jet.at(fI).push_back(temphistDelta_Jet);
          }
          
      }
  }
  
  unsigned int nf = fileList.size();
  unsigned int nTriangNum = nf*(nf-1)/2.;
  // check to make sure that the correct number of histograms were created
  if (hist.size() != nf) std::cout<<"Incorrect number of file lists created for hist"<<std::endl;
  if (hist.at(0).size() != (unsigned int) nVarToComp) std::cout<<"Incorrect number of histograms created"<<std::endl;
  
  if (hist_Jet.size() != nf) std::cout<<"Incorrect number of file lists created for hist_Jet"<<std::endl;
  if (hist_Jet.at(0).size() != (unsigned int) nJetTrees) std::cout<<"Incorrect number of jet tree lists created"<<std::endl;
  if (hist_Jet.at(0).at(0).size() != (unsigned int) listOfJetTreeBranches.at(0).size() ) std::cout<<"Incorrect number of jet tree histograms created"<<std::endl;
  
  if(nf>=2) // if there are 2 or more files
  {
    if (histRat.size()!= nTriangNum) std::cout<<"Incorrect number of file lists created for histRat"<<std::endl;
    if (histRat.at(0).size()!= (unsigned int) nVarToComp) std::cout<<"Incorrect number of ratio histograms created"<<std::endl;
    
    if (histDelta.size()!= nTriangNum) std::cout<<"Incorrect number of file lists created for histDelta"<<std::endl;
    if (histDelta.at(0).size()!= (unsigned int) nVarToComp) std::cout<<"Incorrect number of delta histograms created"<<std::endl;

    if (histRat_Jet.size() != nTriangNum) std::cout<<"Incorrect number of file lists created for histRat_Jet"<<std::endl;
    if (histRat_Jet.at(0).size() != (unsigned int) nJetTrees) std::cout<<"Incorrect number of jet tree ratio lists created"<<std::endl;
    if (histRat_Jet.at(0).at(0).size() != (unsigned int) listOfJetTreeBranches.at(0).size() ) std::cout<<"Incorrect number of jet tree ratio histograms created"<<std::endl;
      
    if (histDelta_Jet.size()!= nTriangNum) std::cout<<"Incorrect number of jet delta histograms created"<<std::endl;
    if (histDelta_Jet.at(0).size() != (unsigned int) nJetTrees) std::cout<<"Incorrect number of jet tree delta lists created"<<std::endl;
    if (histDelta_Jet.at(0).at(0).size() != (unsigned int) listOfJetTreeBranches.at(0).size() ) std::cout<<"Incorrect number of jet tree delta histograms created"<<std::endl;
  }
  
    
    
  // check that the brachMax and brachMin are not equal and add a little bit of space to make the plot look nicer
  for(Int_t i = 0; i < nVarToComp; ++i){
    double tempInterval = branchMaxs.at(i) - branchMins.at(i);
    branchMaxs.at(i) += tempInterval/10.;
    branchMins.at(i) -= tempInterval/10.;

    if(tempInterval == 0){
      branchMaxs.at(i) += 1;
      branchMins.at(i) -= 1;
    }
    //std::cout << listOfCompBranches.at(i)<<" "<< branchMins.at(i) << ", " << branchMaxs.at(i) << std::endl;
  // initialize the histograms for the particle data
    unsigned int loc = 0;
    for(unsigned int hI = 0; hI < fileList.size(); ++hI)
    {
        hist.at(hI).at(i) = new TH1F(Form((listOfCompBranches.at(i) + "_file%d_h").c_str(),hI), (";" + listOfCompBranches.at(i) + ";Counts").c_str(), 100, branchMins.at(i), branchMaxs.at(i));
        for(unsigned int hII = hI+1; hII < fileList.size(); ++hII)
        {
            histRat.at(loc).at(i) = new TH1F(Form((listOfCompBranches.at(i) + "_file_Rat%dOver%d_h").c_str(),hI,hII), (";" + listOfCompBranches.at(i) + Form(";File%d/File%d",hI,hII)).c_str(), 100, branchMins.at(i), branchMaxs.at(i));
            histDelta.at(loc).at(i) = new TH1F(Form((listOfCompBranches.at(i) + "_Delta%dFrom%d_EvtByEvt_h").c_str(),hI,hII), (";" + listOfCompBranches.at(i) + Form("_{File %d}",hI) + "-" + listOfCompBranches.at(i) + Form("_{File %d};Counts",hII)).c_str(), 100, -1, 1);
            centerTitles({histRat.at(loc).at(i), histDelta.at(loc).at(i)});
            ++loc;
        }
        centerTitles(hist.at(hI).at(i));
    }
    if ( (loc) != nTriangNum ) std::cout<<"Incorrect number of delta and ratio histograms initialized"<<std::endl;
    //hist1_p[i] = new TH1F((listOfCompBranches.at(i) + "_file1_h").c_str(), (";" + listOfCompBranches.at(i) + ";Counts").c_str(), 100, branchMins.at(i), branchMaxs.at(i));
  }
  
  for(Int_t jI = 0; jI < nJetTrees; ++jI)
  {
    std::string jetStr = listOfJetTrees1.at(jI);
      // fix the interval for the jet max and mins
    for(unsigned int bI = 0; bI < listOfJetTreeBranches.at(jI).size(); ++bI)
    {
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
        
      unsigned int loc = 0;
      // initialize the histograms for the jet data
      for (unsigned int hI = 0; hI < fileList.size(); ++hI)
      {
          // at(file).at(jet Tree).at(branch)
          hist_Jet.at(hI).at(jI).at(bI) = new TH1F(Form((nameStr + "_file%d_h").c_str(),hI), (";" + brStr + " (" + jetStr + ");Counts").c_str(), 100, tempMin, tempMax);
          for (unsigned int hII = hI+1; hII < fileList.size(); ++hII)
          {
              histRat_Jet.at(loc).at(jI).at(bI) = new TH1F(Form((nameStr + "_file_Rat%dOver%d_h").c_str(),hI,hII), (";" + brStr + " (" + jetStr + Form(");File%d/File%d",hI,hII)).c_str(), 100, tempMin, tempMax);
              histDelta_Jet.at(loc).at(jI).at(bI) = new TH1F(Form((nameStr + "_Delta%dFrom%d_EvtByEvt_h").c_str(),hI,hII), (";" + brStr + Form("_{File %d}",hI) + "-" + brStr + Form("_{File %d}",hII) + " (" + jetStr + ");Counts").c_str(), 100, -1, 1);
              centerTitles({histRat_Jet.at(loc).at(jI).at(bI), histDelta_Jet.at(loc).at(jI).at(bI)});
              ++loc;
          }
          centerTitles(hist_Jet.at(hI).at(jI).at(bI));
      }
      if ( (loc) != nTriangNum ) std::cout<<"Incorrect number of jet delta and ratio histograms initialized"<<std::endl;
    }
  }
  ////////// DONE INITIALIZING HISTOGRAMS //////////
    
    
  //// TURN ON THE BRANCHES THAT WE WILL BE PLOTTING ////
         
  //.at(file).at(tree)
  std::vector< std::vector<TTree*> > jetTree;
  for (unsigned int fI = 0; fI < fileList.size(); ++fI)
  {
      inFile.at(fI) = new TFile(fileList.at(fI).c_str(), "READ");
      inTree.at(fI) = (TTree*)inFile.at(fI)->Get("t");
      inTree.at(fI)->SetBranchStatus("*", 0);
      // turn on the necessary branches
      pData.at(fI).SetStatusAndAddressRead(inTree.at(fI), listOfCompBranches);
      eData.at(fI).SetStatusAndAddressRead(inTree.at(fI), listOfCompBranches);
      // load all of the necessary jet trees and turn on the necessary branches
      TTree* jetTree_p[nJetTrees];
      for(Int_t jI = 0; jI < nJetTrees; ++jI){
          jetTree_p[jI] = (TTree*)inFile.at(fI)->Get(listOfJetTrees1.at(jI).c_str());
          jetTree_p[jI]->SetBranchStatus("*", 0);
          jData.at(fI).at(jI).SetStatusAndAddressRead(jetTree_p[jI], listOfJetTreeBranches.at(jI));
      }
      std::vector<TTree*> jettree;
      jettree.assign(jetTree_p, jetTree_p + nJetTrees);
      jetTree.push_back(jettree);
  }
        //// DONE TURNING ON THE BRANCHES THAT WE WILL BE PLOTTING ////

  Int_t nEntries = inTree.at(0)->GetEntries();
  if (nEntries > 1000000) nEntries = 1000000;
  std::cout << "Doing full processing..." << std::endl;
  Selection s = Selection();

  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%1000 == 0) std::cout << " Entry: " << entry << "/" << nEntries << std::endl;
    
    inTree.at(0)->GetEntry(entry);
    
    // look at around the 91.5 GeV range for LEP2
    if (doECut & (pData.at(0).Energy > 100.0)) continue; 

    // apply the event selection criteria
    Int_t nTrk = 0;
    if(!s.donTrkThrust) nTrk = s.ridge_eventSelection(eData.at(0).passesWW, eData.at(0).missP, pData.at(0).nParticle, jData.at(0).at(0).nref, jData.at(0).at(0).jtpt, jData.at(0).at(0).jteta, eData.at(0).STheta, pData.at(0).mass, pData.at(0).ntpc, pData.at(0).theta, pData.at(0).pmag, pData.at(0).d0, pData.at(0).z0, pData.at(0).pwflag);
    if(s.donTrkThrust) nTrk = s.ridge_eventSelection(eData.at(0).passesWW, eData.at(0).missP, pData.at(0).nParticle, jData.at(0).at(0).nref, jData.at(0).at(0).jtpt, jData.at(0).at(0).jteta, eData.at(0).STheta, pData.at(0).mass, pData.at(0).ntpc, pData.at(0).theta, pData.at(0).pmag, pData.at(0).d0, pData.at(0).z0, pData.at(0).pwflag);
    if( nTrk < 0) continue;


    if(inFileName.find("LEP2") != std::string::npos && eData.at(0).passesWW == 0) continue;
    for(Int_t jI = 0; jI < nJetTrees; ++jI){jetTree.at(0)[jI]->GetEntry(entry);}
    // check if the key is in the list

    // load the remaining trees and apply the event selection
    for( unsigned int fI = 1; fI < fileList.size(); ++fI)
    {
        inTree.at(fI)->GetEntry(entry);
    }
    /////// FILLING THE HISTOGRAMS ////////
    for(unsigned int lI = 0; lI < listOfCompBranches.size(); ++lI){
      std::string tempS = listOfCompBranches.at(lI);
        
      std::vector< TH1F* > hist_p = {};
      for(unsigned int hI = 0; hI < hist.size(); ++hI){hist_p.push_back(hist.at(hI).at(lI));}
      std::vector< TH1F* >  histDelta_p = {};
      for(unsigned int hI = 0; hI < histDelta.size(); ++hI){histDelta_p.push_back(histDelta.at(hI).at(lI));}
      
      if(tempS.find("nParticle") != std::string::npos && tempS.size() == std::string("nParticle").size())
      {
          std::vector< Int_t > pData_nParticle;
          for(unsigned int pI = 0; pI < pData.size(); ++pI){pData_nParticle.push_back(pData.at(pI).nParticle);}
          doFill(hist_p,histDelta_p,pData_nParticle);
      }
      else if(tempS.find("EventNo") != std::string::npos && tempS.size() == std::string("EventNo").size())
      {
          std::vector< Int_t > pData_EventNo;
          for(unsigned int pI = 0; pI < pData.size(); ++pI){pData_EventNo.push_back(pData.at(pI).EventNo);}
          doFill(hist_p,histDelta_p,pData_EventNo);
      }
      else if(tempS.find("RunNo") != std::string::npos && tempS.size() == std::string("RunNo").size())
      {
          std::vector< Int_t > pData_RunNo;
          for(unsigned int pI = 0; pI < pData.size(); ++pI){pData_RunNo.push_back(pData.at(pI).RunNo);}
          doFill(hist_p,histDelta_p,pData_RunNo);
      }
      else if(tempS.find("year") != std::string::npos && tempS.size() == std::string("year").size())
      {
          std::vector< Int_t > pData_year;
          for(unsigned int pI = 0; pI < pData.size(); ++pI){pData_year.push_back(pData.at(pI).year);}
          doFill(hist_p,histDelta_p,pData_year);
      }
      else if(tempS.find("process") != std::string::npos && tempS.size() == std::string("process").size())
      {
          std::vector< Int_t > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI){pData_process.push_back(pData.at(pI).process);}
          doFill(hist_p,histDelta_p,pData_process);
      }
      else if(tempS.find("Energy") != std::string::npos && tempS.size() == std::string("Energy").size())
      {
          std::vector< Float_t > pData_temp;
          for(unsigned int pI = 0; pI < pData.size(); ++pI){pData_temp.push_back(pData.at(pI).Energy);}
          doFill(hist_p,histDelta_p,pData_temp);
      }
      else if(tempS.find("bFlag") != std::string::npos && tempS.size() == std::string("bFlag").size())
      {
          std::vector< Float_t > pData_temp;
          for(unsigned int pI = 0; pI < pData.size(); ++pI){pData_temp.push_back(pData.at(pI).bFlag);}
          doFill(hist_p,histDelta_p,pData_temp);
      }
      else if(tempS.find("bx") != std::string::npos && tempS.size() == std::string("bx").size())
      {
          std::vector< Float_t > pData_temp;
          for(unsigned int pI = 0; pI < pData.size(); ++pI){pData_temp.push_back(pData.at(pI).bx);}
          doFill(hist_p,histDelta_p,pData_temp);
      }
      else if(tempS.find("by") != std::string::npos && tempS.size() == std::string("by").size())
      {
          std::vector< Float_t > pData_temp;
          for(unsigned int pI = 0; pI < pData.size(); ++pI){pData_temp.push_back(pData.at(pI).by);}
          doFill(hist_p,histDelta_p,pData_temp);
      }
      else if(tempS.find("ebx") != std::string::npos && tempS.size() == std::string("ebx").size())
      {
          std::vector< Float_t > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI){pData_process.push_back(pData.at(pI).ebx);}
          doFill(hist_p,histDelta_p,pData_process);
      }
      else if(tempS.find("eby") != std::string::npos && tempS.size() == std::string("eby").size())
      {
          std::vector< Float_t > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI){pData_process.push_back(pData.at(pI).eby);}
          doFill(hist_p,histDelta_p,pData_process);
      }
          
      else if(tempS.find("px") != std::string::npos && tempS.size() == std::string("px").size())
      {
          std::vector< std::vector<Float_t> > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).px,pData.at(pI).px+pData.at(pI).nParticle);
              pData_process.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_process);
      }
        
      else if(tempS.find("py") != std::string::npos && tempS.size() == std::string("py").size())
      {
          std::vector< std::vector<Float_t> > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).py,pData.at(pI).py+pData.at(pI).nParticle);
              pData_process.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_process);
      }
      else if(tempS.find("pz") != std::string::npos && tempS.size() == std::string("pz").size())
      {
          std::vector< std::vector<Float_t> > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).pz,pData.at(pI).pz+pData.at(pI).nParticle);
              pData_process.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_process);
      }
      else if(tempS.find("pt") != std::string::npos && tempS.size() == std::string("pt").size())
      {
          std::vector< std::vector<Float_t> > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).pt,pData.at(pI).pt+pData.at(pI).nParticle);
              pData_process.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_process);
      }
      else if(tempS.find("pmag") != std::string::npos && tempS.size() == std::string("pmag").size())
      {
          std::vector< std::vector<Float_t> > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).pmag,pData.at(pI).pmag+pData.at(pI).nParticle);
              pData_process.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_process);
      }
      else if(tempS.find("rap") != std::string::npos && tempS.size() == std::string("rap").size())
      {
          std::vector< std::vector<Float_t> > pData_temp;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).rap,pData.at(pI).rap+pData.at(pI).nParticle);
              pData_temp.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_temp);
      }
      else if(tempS.find("eta") != std::string::npos && tempS.size() == std::string("eta").size())
      {
          std::vector< std::vector<Float_t> > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).eta,pData.at(pI).eta+pData.at(pI).nParticle);
              pData_process.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_process);
      }
      else if(tempS.find("theta") != std::string::npos && tempS.size() == std::string("theta").size())
      {
          std::vector< std::vector<Float_t> > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).theta,pData.at(pI).theta+pData.at(pI).nParticle);
              pData_process.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_process);
      }
      else if(tempS.find("phi") != std::string::npos && tempS.size() == std::string("phi").size())
      {
          std::vector< std::vector<Float_t> > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).phi,pData.at(pI).phi+pData.at(pI).nParticle);
              pData_process.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_process);
      }
      else if(tempS.find("mass") != std::string::npos && tempS.size() == std::string("mass").size())
      {
          std::vector< std::vector<Float_t> > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).mass,pData.at(pI).mass+pData.at(pI).nParticle);
              pData_process.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_process);
      }
      else if(tempS.find("charge") != std::string::npos && tempS.size() == std::string("charge").size())
      {
          std::vector< std::vector<Float_t> > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).charge,pData.at(pI).charge+pData.at(pI).nParticle);
              pData_process.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_process);
      }
      else if(tempS.find("pwflag") != std::string::npos && tempS.size() == std::string("pwflag").size())
      {
          std::vector< std::vector<Float_t> > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).pwflag,pData.at(pI).pwflag+pData.at(pI).nParticle);
              pData_process.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_process);
      }
      else if(tempS.find("pid") != std::string::npos && tempS.size() == std::string("pid").size())
      {
          std::vector< std::vector<Float_t> > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).pid,pData.at(pI).pid+pData.at(pI).nParticle);
              pData_process.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_process);
      }
      else if(tempS.find("d0") != std::string::npos && tempS.size() == std::string("d0").size())
      {
          std::vector< std::vector<Float_t> > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).d0,pData.at(pI).d0+pData.at(pI).nParticle);
              pData_process.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_process);
      }
      else if(tempS.find("z0") != std::string::npos && tempS.size() == std::string("z0").size())
      {
          std::vector< std::vector<Float_t> > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).z0,pData.at(pI).z0+pData.at(pI).nParticle);
              pData_process.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_process);
      }
      else if(tempS.find("ntpc") != std::string::npos && tempS.size() == std::string("ntpc").size())
      {
          std::vector< std::vector<Int_t> > pData_temp;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Int_t> temp;
              temp.assign(pData.at(pI).ntpc,pData.at(pI).ntpc+pData.at(pI).nParticle);
              pData_temp.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_temp);
      }
      else if(tempS.find("nitc") != std::string::npos && tempS.size() == std::string("nitc").size())
      {
          std::vector< std::vector<Int_t> > pData_temp;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Int_t> temp;
              temp.assign(pData.at(pI).nitc,pData.at(pI).nitc+pData.at(pI).nParticle);
              pData_temp.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_temp);
      }
      else if(tempS.find("nvdet") != std::string::npos && tempS.size() == std::string("nvdet").size())
      {
          std::vector< std::vector<Int_t> > pData_temp;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Int_t> temp;
              temp.assign(pData.at(pI).nvdet,pData.at(pI).nvdet+pData.at(pI).nParticle);
              pData_temp.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_temp);
      }
      else if(tempS.find("vx") != std::string::npos && tempS.size() == std::string("vx").size())
      {
          std::vector< std::vector<Int_t> > pData_temp;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Int_t> temp;
              temp.assign(pData.at(pI).vx,pData.at(pI).vx+pData.at(pI).nParticle);
              pData_temp.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_temp);
      }
      else if(tempS.find("vy") != std::string::npos && tempS.size() == std::string("vy").size())
      {
          std::vector< std::vector<Int_t> > pData_temp;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Int_t> temp;
              temp.assign(pData.at(pI).vy,pData.at(pI).vy+pData.at(pI).nParticle);
              pData_temp.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_temp);
      }
      else if(tempS.find("vz") != std::string::npos && tempS.size() == std::string("vz").size())
      {
          std::vector< std::vector<Int_t> > pData_temp;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Int_t> temp;
              temp.assign(pData.at(pI).vz,pData.at(pI).vz+pData.at(pI).nParticle);
              pData_temp.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_temp);
      }
      else if(tempS.find("pt_wrtThr") != std::string::npos && tempS.size() == std::string("pt_wrtThr").size())
      {
          std::vector< std::vector<Float_t> > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).pt_wrtThr,pData.at(pI).pt_wrtThr+pData.at(pI).nParticle);
              pData_process.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_process);
      }
      else if(tempS.find("eta_wrtThr") != std::string::npos && tempS.size() == std::string("eta_wrtThr").size())
      {
          std::vector< std::vector<Float_t> > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).eta_wrtThr,pData.at(pI).eta_wrtThr+pData.at(pI).nParticle);
              pData_process.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_process);
      }
      else if(tempS.find("rap_wrtThr") != std::string::npos && tempS.size() == std::string("rap_wrtThr").size())
      {
          std::vector< std::vector<Float_t> > pData_temp;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).rap_wrtThr,pData.at(pI).rap_wrtThr+pData.at(pI).nParticle);
              pData_temp.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_temp);
      }
      else if(tempS.find("theta_wrtThr") != std::string::npos && tempS.size() == std::string("theta_wrtThr").size())
      {
          std::vector< std::vector<Float_t> > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).theta_wrtThr,pData.at(pI).theta_wrtThr+pData.at(pI).nParticle);
              pData_process.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_process);
      }
      else if(tempS.find("phi_wrtThr") != std::string::npos && tempS.size() == std::string("phi_wrtThr").size())
      {
          std::vector< std::vector<Float_t> > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).phi_wrtThr,pData.at(pI).phi_wrtThr+pData.at(pI).nParticle);
              pData_process.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_process);
      }
      else if(tempS.find("pt_wrtChThr") != std::string::npos && tempS.size() == std::string("pt_wrtChThr").size())
      {
          std::vector< std::vector<Float_t> > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).pt_wrtChThr,pData.at(pI).pt_wrtChThr+pData.at(pI).nParticle);
              pData_process.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_process);
      }
      else if(tempS.find("eta_wrtChThr") != std::string::npos && tempS.size() == std::string("eta_wrtChThr").size())
      {
          std::vector< std::vector<Float_t> > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).eta_wrtChThr,pData.at(pI).eta_wrtChThr+pData.at(pI).nParticle);
              pData_process.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_process);
      }
      else if(tempS.find("rap_wrtChThr") != std::string::npos && tempS.size() == std::string("rap_wrtChThr").size())
      {
          std::vector< std::vector<Float_t> > pData_temp;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).rap_wrtChThr,pData.at(pI).rap_wrtChThr+pData.at(pI).nParticle);
              pData_temp.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_temp);
      }
      else if(tempS.find("theta_wrtChThr") != std::string::npos && tempS.size() == std::string("theta_wrtChThr").size())
      {
          std::vector< std::vector<Float_t> > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).theta_wrtChThr,pData.at(pI).theta_wrtChThr+pData.at(pI).nParticle);
              pData_process.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_process);
      }
      else if(tempS.find("phi_wrtChThr") != std::string::npos && tempS.size() == std::string("phi_wrtChThr").size())
      {
          std::vector< std::vector<Float_t> > pData_process;
          for(unsigned int pI = 0; pI < pData.size(); ++pI)
          {
              std::vector<Float_t> temp;
              temp.assign(pData.at(pI).phi_wrtChThr,pData.at(pI).phi_wrtChThr+pData.at(pI).nParticle);
              pData_process.push_back(temp);
          }
          doFillArr(hist_p,histDelta_p,pData_process);
      }

      // start of eData plotting
      else if(tempS.find("passesWW") != std::string::npos && tempS.size() == std::string("passesWW").size())
      {
          std::vector< Int_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).passesWW);}
          doFill(hist_p,histDelta_p,eData_temp);
      }  
      else if(tempS.find("missP") != std::string::npos && tempS.size() == std::string("missP").size())
      {
          std::vector< Float_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).missP);}
          doFill(hist_p,histDelta_p,eData_temp);
      }
      else if(tempS.find("missPt") != std::string::npos && tempS.size() == std::string("missPt").size())
      {
          std::vector< Float_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).missPt);}
          doFill(hist_p,histDelta_p,eData_temp);
      }
      else if(tempS.find("missTheta") != std::string::npos && tempS.size() == std::string("missTheta").size())
      {
          std::vector< Float_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).missTheta);}
          doFill(hist_p,histDelta_p,eData_temp);
      }
      else if(tempS.find("missPhi") != std::string::npos && tempS.size() == std::string("missPhi").size())
      {
          std::vector< Float_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).missPhi);}
          doFill(hist_p,histDelta_p,eData_temp);
      }
      else if(tempS.find("missChargedP") != std::string::npos && tempS.size() == std::string("missChargedP").size())
      {
          std::vector< Float_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).missChargedP);}
          doFill(hist_p,histDelta_p,eData_temp);
      }
      else if(tempS.find("missChargedPt") != std::string::npos && tempS.size() == std::string("missChargedPt").size())
      {
          std::vector< Float_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).missChargedPt);}
          doFill(hist_p,histDelta_p,eData_temp);
      }
      else if(tempS.find("missChargedTheta") != std::string::npos && tempS.size() == std::string("missChargedTheta").size())
      {
          std::vector< Float_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).missChargedTheta);}
          doFill(hist_p,histDelta_p,eData_temp);
      }
      else if(tempS.find("missChargedPhi") != std::string::npos && tempS.size() == std::string("missChargedPhi").size())
      {
          std::vector< Float_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).missChargedPhi);}
          doFill(hist_p,histDelta_p,eData_temp);
      }
      else if(tempS.find("nChargedHadrons") != std::string::npos && tempS.size() == std::string("nChargedHadrons").size())
      {
          std::vector< Int_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).nChargedHadrons);}
          doFill(hist_p,histDelta_p,eData_temp);
      }
      else if(tempS.find("nChargedHadrons_GT0p4") != std::string::npos && tempS.size() == std::string("nChargedHadrons_GT0p4").size())
      {
          std::vector< Int_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).nChargedHadrons_GT0p4);}
          doFill(hist_p,histDelta_p,eData_temp);
      }
      else if(tempS.find("nChargedHadrons_GT0p4Thrust") != std::string::npos && tempS.size() == std::string("nChargedHadrons_GT0p4Thrust").size())
      {
          std::vector< Int_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).nChargedHadrons_GT0p4Thrust);}
          doFill(hist_p,histDelta_p,eData_temp);
      }
      else if(tempS.find("Thrust") != std::string::npos && tempS.size() == std::string("Thrust").size())
      {
          std::vector< Float_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).Thrust);}
          doFill(hist_p,histDelta_p,eData_temp);
      }
      else if(tempS.find("TTheta") != std::string::npos && tempS.size() == std::string("TTheta").size())
      {
          std::vector< Float_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).TTheta);}
          doFill(hist_p,histDelta_p,eData_temp);
      }
      else if(tempS.find("TPhi") != std::string::npos && tempS.size() == std::string("TPhi").size())
      {
          std::vector< Float_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).TPhi);}
          doFill(hist_p,histDelta_p,eData_temp);
      }
      else if(tempS.find("Thrust_charged") != std::string::npos && tempS.size() == std::string("Thrust_charged").size())
      {
          std::vector< Float_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).Thrust_charged);}
          doFill(hist_p,histDelta_p,eData_temp);
      }
      else if(tempS.find("TTheta_charged") != std::string::npos && tempS.size() == std::string("TTheta_charged").size())
      {
          std::vector< Float_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).TTheta_charged);}
          doFill(hist_p,histDelta_p,eData_temp);
      }
      else if(tempS.find("TPhi_charged") != std::string::npos && tempS.size() == std::string("TPhi_charged").size())
      {
          std::vector< Float_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).TPhi_charged);}
          doFill(hist_p,histDelta_p,eData_temp);
      }
      else if(tempS.find("Sphericity") != std::string::npos && tempS.size() == std::string("Sphericity").size())
      {
          std::vector< Float_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).Sphericity);}
          doFill(hist_p,histDelta_p,eData_temp);
      }   
      else if(tempS.find("STheta") != std::string::npos && tempS.size() == std::string("STheta").size())
      {
          std::vector< Float_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).STheta);}
          doFill(hist_p,histDelta_p,eData_temp);
      }  
      else if(tempS.find("SPhi") != std::string::npos && tempS.size() == std::string("SPhi").size())
      {
          std::vector< Float_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).SPhi);}
          doFill(hist_p,histDelta_p,eData_temp);
      } 
      else if(tempS.find("Aplanarity") != std::string::npos && tempS.size() == std::string("Aplanarity").size())
      {
          std::vector< Float_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).Aplanarity);}
          doFill(hist_p,histDelta_p,eData_temp);
      }    
      else if(tempS.find("Sphericity_linearized") != std::string::npos && tempS.size() == std::string("Sphericity_linearized").size())
      {
          std::vector< Float_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).Sphericity_linearized);}
          doFill(hist_p,histDelta_p,eData_temp);
      } 
      else if(tempS.find("SPhi_linearized") != std::string::npos && tempS.size() == std::string("SPhi_linearized").size())
      {
          std::vector< Float_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).SPhi_linearized);}
          doFill(hist_p,histDelta_p,eData_temp);
      } 
      else if(tempS.find("Aplanarity_linearized") != std::string::npos && tempS.size() == std::string("Aplanarity_linearized").size())
      {
          std::vector< Float_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).Aplanarity_linearized);}
          doFill(hist_p,histDelta_p,eData_temp);
      }
      else if(tempS.find("C_linearized") != std::string::npos && tempS.size() == std::string("C_linearized").size())
      {
          std::vector< Float_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).C_linearized);}
          doFill(hist_p,histDelta_p,eData_temp);
      }
      else if(tempS.find("D_linearized") != std::string::npos && tempS.size() == std::string("D_linearized").size())
      {
          std::vector< Float_t > eData_temp;
          for(unsigned int eI = 0; eI < eData.size(); ++eI){eData_temp.push_back(eData.at(eI).D_linearized);}
          doFill(hist_p,histDelta_p,eData_temp);
      }

    }
      
    
    ///////////////////////////// ANTHONY YOU ARE HERE AND NEED TO FINISH SETTING UP THE FILLS FOR THE JETS AND THEN ADD IN THE WRITING/DRAWING FOR JETS /////////////////////////////
    for(Int_t jI = 0; jI < nJetTrees; ++jI)
    {
      for(unsigned bI = 0; bI < listOfJetTreeBranches.at(jI).size(); ++bI)
      {
          std::string tempS = listOfJetTreeBranches.at(jI).at(bI);

          std::vector< TH1F* > hist_Jet_p = {};
          for(unsigned int hI = 0; hI < hist_Jet.size(); ++hI){hist_Jet_p.push_back(hist_Jet.at(hI).at(jI).at(bI));}
          std::vector< TH1F* >  histDelta_Jet_p = {};
          for(unsigned int hI = 0; hI < histDelta_Jet.size(); ++hI){histDelta_Jet_p.push_back(histDelta_Jet.at(hI).at(jI).at(bI));}
          // .at(file).at(jet tree name)

          if(tempS.find("nref") != std::string::npos && tempS.size() == std::string("nref").size())
          {
              std::vector< Int_t > jData_temp;
              for(unsigned int jDI = 0; jDI < jData.size(); ++jDI){jData_temp.push_back(jData.at(jDI).at(jI).nref);}
              doFill(hist_Jet_p,histDelta_Jet_p,jData_temp);
          }
          else if(tempS.find("jtpt") != std::string::npos && tempS.size() == std::string("jtpt").size())
          {
              std::vector< std::vector<Float_t> > jData_temp;
              for(unsigned int jDI = 0; jDI < pData.size(); ++jDI)
              {
                  std::vector<Float_t> temp;
                  temp.assign(jData.at(jDI).at(jI).jtpt,jData.at(jDI).at(jI).jtpt + jData.at(jDI).at(jI).nref);
                  jData_temp.push_back(temp);
              }
              doFillArr(hist_Jet_p,histDelta_Jet_p,jData_temp);
          }
          else if(tempS.find("jteta") != std::string::npos && tempS.size() == std::string("jteta").size())
          {
              std::vector< std::vector<Float_t> > jData_temp;
              for(unsigned int jDI = 0; jDI < pData.size(); ++jDI)
              {
                  std::vector<Float_t> temp;
                  temp.assign(jData.at(jDI).at(jI).jteta,jData.at(jDI).at(jI).jteta + jData.at(jDI).at(jI).nref);
                  jData_temp.push_back(temp);
              }
              doFillArr(hist_Jet_p,histDelta_Jet_p,jData_temp);
          }
          else if(tempS.find("jtphi") != std::string::npos && tempS.size() == std::string("jtphi").size())
          {
              std::vector< std::vector<Float_t> > jData_temp;
              for(unsigned int jDI = 0; jDI < pData.size(); ++jDI)
              {
                  std::vector<Float_t> temp;
                  temp.assign(jData.at(jDI).at(jI).jtphi,jData.at(jDI).at(jI).jtphi + jData.at(jDI).at(jI).nref);
                  jData_temp.push_back(temp);
              }
              doFillArr(hist_Jet_p,histDelta_Jet_p,jData_temp);
          }
          else if(tempS.find("jtm") != std::string::npos && tempS.size() == std::string("jtm").size())
          {
              std::vector< std::vector<Float_t> > jData_temp;
              for(unsigned int jDI = 0; jDI < pData.size(); ++jDI)
              {
                  std::vector<Float_t> temp;
                  temp.assign(jData.at(jDI).at(jI).jtm,jData.at(jDI).at(jI).jtm + jData.at(jDI).at(jI).nref);
                  jData_temp.push_back(temp);
              }
              doFillArr(hist_Jet_p,histDelta_Jet_p,jData_temp);
          }
          else if(tempS.find("jtN") != std::string::npos && tempS.size() == std::string("jtN").size())
          {
              std::vector< std::vector<Int_t> > jData_temp;
              for(unsigned int jDI = 0; jDI < pData.size(); ++jDI)
              {
                  std::vector<Int_t> temp;
                  temp.assign(jData.at(jDI).at(jI).jtN,jData.at(jDI).at(jI).jtN + jData.at(jDI).at(jI).nref);
                  jData_temp.push_back(temp);
              }
              doFillArr(hist_Jet_p,histDelta_Jet_p,jData_temp);
          }
          // not yet plotting jtNPW, jtptFracPW
          //Note cant handle 2d arrays proper will add later
      }
    }
    
  }
  for (int fI = fileList.size()-1; fI >= 0; --fI)
  {
      inFile.at(fI)->Close();
      delete inFile.at(fI);
  }
  outFile_p->cd();
                           
  const Double_t splitPoint = 0.35;
                           
  std::vector<std::string> listOfPdf;
  std::vector<std::string> listOfVar;
                           
  //if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
                           
  // Writing and plotting all of the histograms
  TCanvas* canv_p;
  TPad* pad1; TPad* pad2_p; TPad* pad3_p;
  std::vector<TPad*> pad2;
  std::vector<TPad*> pad3;
  for(int nI = 0; nI < nVarToComp; ++nI)
  {
      if(drawHist)
      {
        canv_p = new TCanvas("canv_c", "canv_c", 1000, 500);
        canv_p->SetTopMargin(0.01);
        canv_p->SetRightMargin(0.01);
        canv_p->SetLeftMargin(0.01);
        canv_p->SetBottomMargin(0.01);
        
        TPad* pad1 = new TPad("pad1", "pad1", 0.0, splitPoint, 0.5, 1.0);
        pad1->Draw();
        pad1->SetTopMargin(0.01);
        pad1->SetRightMargin(0.01);
        pad1->SetBottomMargin(0.01);
        pad1->SetLeftMargin(pad1->GetLeftMargin()*1.3);
        
        std::vector<TPad*> pad2;
        std::vector<TPad*> pad3;
        
        for (unsigned int pI = 0; pI < histRat.size();pI++)
        {
            pad2_p = new TPad("pad2", "pad2", 0.0, 0.0, 0.5, splitPoint);
            pad2_p->Draw();
            pad2_p->SetTopMargin(0.01);
            pad2_p->SetRightMargin(0.01);
            pad2_p->SetBottomMargin(pad2_p->GetLeftMargin()*1./splitPoint);
            pad2_p->SetLeftMargin(pad2_p->GetLeftMargin());
            pad2.push_back(pad2_p);
            
            pad3_p = new TPad("pad3", "pad3", 0.5, 0.0, 1.0, 1.0);
            pad3_p->Draw();
            pad3_p->SetTopMargin(0.01);
            pad3_p->SetRightMargin(0.01);
            pad3_p->SetLeftMargin(pad3_p->GetLeftMargin());
            pad3_p->SetBottomMargin(pad3_p->GetLeftMargin());
            pad3.push_back(pad3_p);
        }
      }
      unsigned int loc = 0;
      for(unsigned int hI = 0; hI < hist.size(); ++hI)
      {
          hist.at(hI).at(nI)->Write("", TObject::kOverwrite);
          hist.at(hI).at(nI)->Sumw2();
          if(drawHist)
          {
            pad1->cd();
            formatTH1F(hist.at(hI).at(nI));
            if(hI == 0) hist.at(hI).at(nI)->DrawCopy("HIST E1");
            else hist.at(hI).at(nI)->DrawCopy("SAME *HIST E1");
          }
          for(unsigned int hII = hI+1; hII < hist.size(); ++hII)
          {
              histRat.at(loc).at(nI)->Sumw2();
              histRat.at(loc).at(nI)->Divide(hist.at(hI).at(nI),hist.at(hII).at(nI));
              histRat.at(loc).at(nI)->Write("", TObject::kOverwrite);
              if(drawHist)
              {
                pad2.at(loc)->cd();
                formatTH1F(histRat.at(loc).at(nI));
                histRat.at(loc).at(nI)->SetMaximum(1.3);
                histRat.at(loc).at(nI)->SetMinimum(0.7);
                histRat.at(loc).at(nI)->DrawCopy("P E1");
              }
              histDelta.at(loc).at(nI)->Write("", TObject::kOverwrite);
              if(drawHist)
              {
                pad3.at(loc)->cd();
                formatTH1F(histDelta.at(loc).at(nI));
                histDelta.at(loc).at(nI)->DrawCopy("HIST E1");
              }  
              
              ++loc;
          }
      }
      if ( (loc) != nTriangNum ) std::cout<<"Incorrect number of jet delta and ratio histograms initialized"<<std::endl;
      
      std::string pdfStr = listOfCompBranches.at(nI) + "_" + outFileName + "_" + std::to_string(date->GetDate()) + ".pdf";
      if(drawHist) canv_p->SaveAs(("../pdfDir/" + pdfStr).c_str());
      listOfPdf.push_back(pdfStr);
      listOfVar.push_back(listOfCompBranches.at(nI));
      
      // clean up before next variable
      if(drawHist)
      {
        delete pad1;
        for (unsigned int pI = 0; pI < histRat.size(); pI++)
        {
            delete pad2.at(pI);
            delete pad3.at(pI);
        }
        delete canv_p;
      }
      for (unsigned int hI = 0; hI < hist.size(); hI++) delete hist.at(hI).at(nI);
      for (unsigned int hI = 0; hI < histRat.size(); hI++){ delete histRat.at(hI).at(nI);delete histDelta.at(hI).at(nI); }
  }

  // ADD JET PLOTTING IN AFTER THE ABOVE WORKS
    std::cout<<"nJetTrees "<<nJetTrees<<" "<< listOfJetTreeBranches.at(0).size() << " "<<listOfJetTreeBranches.at(1).size() << " "<<listOfJetTreeBranches.at(2).size() << " "<<listOfJetTreeBranches.at(3).size() << " "<<std::endl;
    for(int jI = 0; jI < nJetTrees; ++jI)
    {
        const std::string jetStr = listOfJetTrees1.at(jI);
        for(unsigned int bI = 0; bI < listOfJetTreeBranches.at(jI).size(); ++bI)
        {
          if(drawHist)
          {
            canv_p = new TCanvas("canv_c", "canv_c", 1000, 500);
            canv_p->SetTopMargin(0.01);
            canv_p->SetRightMargin(0.01);
            canv_p->SetLeftMargin(0.01);
            canv_p->SetBottomMargin(0.01);
            
            pad1 = new TPad("pad1", "pad1", 0.0, splitPoint, 0.5, 1.0);
            pad1->Draw();
            pad1->SetTopMargin(0.01);
            pad1->SetRightMargin(0.01);
            pad1->SetBottomMargin(0.01);
            pad1->SetLeftMargin(pad1->GetLeftMargin()*1.3);
            
            std::vector<TPad*> pad2;
            std::vector<TPad*> pad3;

            for (unsigned int fI = 0; fI < histRat_Jet.size(); fI++)
            {
                pad2_p = new TPad("pad2", "pad2", 0.0, 0.0, 0.5, splitPoint);
                pad2_p->Draw();
                pad2_p->SetTopMargin(0.01);
                pad2_p->SetRightMargin(0.01);
                pad2_p->SetBottomMargin(pad2_p->GetLeftMargin()*1./splitPoint);
                pad2_p->SetLeftMargin(pad2_p->GetLeftMargin());
                pad2.push_back(pad2_p);
                
                pad3_p = new TPad("pad3", "pad3", 0.5, 0.0, 1.0, 1.0);
                pad3_p->Draw();
                pad3_p->SetTopMargin(0.01);
                pad3_p->SetRightMargin(0.01);
                pad3_p->SetLeftMargin(pad3_p->GetLeftMargin());
                pad3_p->SetBottomMargin(pad3_p->GetLeftMargin());
                pad3.push_back(pad3_p);
             
            }
          }
            unsigned int loc = 0;
            for(unsigned int hI = 0; hI < hist_Jet.size(); ++hI)
            {
                hist_Jet.at(hI).at(jI).at(bI)->Write("", TObject::kOverwrite);
                hist_Jet.at(hI).at(jI).at(bI)->Sumw2();
                if(drawHist)
                {
                  pad1->cd();
                  formatTH1F(hist_Jet.at(hI).at(jI).at(bI));
                  if(hI == 0) hist_Jet.at(hI).at(jI).at(bI)->DrawCopy("HIST E1");
                  else hist_Jet.at(hI).at(jI).at(bI)->DrawCopy("SAME *HIST E1");
                }

                for(unsigned int hII = hI+1; hII < hist_Jet.size(); ++hII)
                {
                    histRat_Jet.at(loc).at(jI).at(bI)->Sumw2();
                    histRat_Jet.at(loc).at(jI).at(bI)->Divide(hist_Jet.at(hI).at(jI).at(bI),hist_Jet.at(hII).at(jI).at(bI));
                    histRat_Jet.at(loc).at(jI).at(bI)->Write("", TObject::kOverwrite);
                    if(drawHist)
                    {
                      pad2.at(loc)->cd();
                      formatTH1F(histRat_Jet.at(loc).at(jI).at(bI));
                      histRat_Jet.at(loc).at(jI).at(bI)->SetMaximum(1.3);
                      histRat_Jet.at(loc).at(jI).at(bI)->SetMinimum(0.7);
                      histRat_Jet.at(loc).at(jI).at(bI)->DrawCopy("P E1");
                    }
                    histDelta_Jet.at(loc).at(jI).at(bI)->Write("", TObject::kOverwrite);
                    if(drawHist)
                    {
                      pad3.at(loc)->cd();
                      formatTH1F(histDelta_Jet.at(loc).at(jI).at(bI));
                      histDelta_Jet.at(loc).at(jI).at(bI)->DrawCopy("HIST E1");
                    }
                    ++loc;
                }
            }

            if ( (loc) != nTriangNum ) std::cout<<"Incorrect number of jet delta and ratio histograms written"<<std::endl;
            
            std::string pdfStr = listOfJetTreeBranches.at(jI).at(bI) + "_" + jetStr+ "_" + outFileName + "_" + std::to_string(date->GetDate()) + ".pdf";
            if(drawHist)canv_p->SaveAs(("../pdfDir/" + pdfStr).c_str());
            listOfPdf.push_back(pdfStr);
            listOfVar.push_back((jetStr + "_" + listOfJetTreeBranches.at(jI).at(bI)).c_str());
            // clean up before next variable
            if(drawHist)
            {
              delete pad1;
              for (unsigned int pI = 0; pI < histRat_Jet.size(); pI++)
              {
                  delete pad2.at(pI);
                  delete pad3.at(pI);
              }
              delete canv_p;
            }
            for(unsigned int hI = 0; hI < hist_Jet.size(); ++hI) delete hist_Jet.at(hI).at(jI).at(bI);
            for(unsigned int hI = 0; hI < histRat_Jet.size(); ++hI) { delete histRat_Jet.at(hI).at(jI).at(bI); histDelta_Jet.at(hI).at(jI).at(bI); }
        }
    }
    
  std::vector<TNamed> nameFile;
  for(unsigned int fI = 0; fI < fileList.size(); ++fI)
  {
      TNamed nameFile1(Form("nameFile%d",fI), fileList.at(fI).c_str());
      nameFile1.Write("", TObject::kOverwrite);
  }
                    
  outFile_p->Close();
  delete outFile_p;

  // deleted all latex stuff after this add all of that in later after debugging

  delete date;

  return 0;
}

                           
int main(int argc, char* argv[])
{
  if(argc != 2 && argc != 3 && argc != 4){
    std::cout << "Usage ./plotDQC.exe <inFileName1> <outFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 2) retVal += plotDQC(argv[1]);
  else if(argc == 3) retVal += plotDQC(argv[1], argv[2]);
  else if(argc == 4) retVal += plotDQC(argv[1], argv[2], atoi(argv[3]));
  else if(argc == 5) retVal += plotDQC(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]));
  return retVal;
}
