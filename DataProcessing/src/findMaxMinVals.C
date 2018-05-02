#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "TFile.h"
#include "TTree.h"

#include "DataProcessing/include/particleData.h"

int findMaxMinVals(const std::string inFileName)
{
  if(inFileName.find(".txt") == std::string::npos){
    std::cout << "inFileName \'" << inFileName << "\' is not a .txt file. return 1" << std::endl;
    return 1;
  }

  std::vector<std::string> fileList;
  std::ifstream pathFile(inFileName.c_str());
  std::string paths;
  while(std::getline(pathFile, paths)){
    if(paths.size() == 0) continue;
    fileList.push_back(paths);
  }
  pathFile.close();

  particleData pData;

  Int_t maxRunNo = -1;
  Int_t maxEventNo = -1;
  Int_t maxSubDir = -1;
  Int_t maxYear = -1;
  Int_t maxProcess = -1;

  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    std::cout << "File " << fI << "/" << fileList.size() << std::endl;

    TFile* inFile_p = new TFile(fileList.at(fI).c_str(), "READ");
    TTree* inTree_p = (TTree*)inFile_p->Get("t");

    std::vector<std::string> listOfCompBranches = {"RunNo", "EventNo", "subDir", "year", "process"};

    pData.SetStatusAndAddressRead(inTree_p, listOfCompBranches);

    const Int_t nEntries = inTree_p->GetEntries();

    for(Int_t entry = 0; entry < nEntries; ++entry){
      inTree_p->GetEntry(entry);

      if(maxRunNo < pData.RunNo) maxRunNo = pData.RunNo;
      if(maxEventNo < pData.EventNo) maxEventNo = pData.EventNo;
      if(maxSubDir < pData.subDir) maxSubDir = pData.subDir;
      if(maxYear < pData.year) maxYear = pData.year;
      if(maxProcess < pData.process) maxProcess = pData.process;
    }

    inFile_p->Close();
    delete inFile_p;
  }

  std::cout << "Max RunNo: " << maxRunNo << std::endl;
  std::cout << "Max EventNo: " << maxEventNo << std::endl;
  std::cout << "Max SubDir: " << maxSubDir << std::endl;
  std::cout << "Max Year: " << maxYear << std::endl;
  std::cout << "Max Process: " << maxProcess << std::endl;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){ 
    std::cout << "Usage: ./findMaxMinVals.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += findMaxMinVals(argv[1]);
  return retVal;
}
