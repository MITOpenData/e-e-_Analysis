#ifndef HEPPLAINTXTREADER
#define HEPPLAINTXTREADER

//Original Author: Chris McGinn
//basic reader for comparisons
//cpp dependencies
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

//ROOT dependencies
#include "TH1D.h"

class hepPlainTxtReader{
 public:
  bool isInit;
  std::vector<std::string> histStrings;
  std::vector<double> xValLow;
  std::vector<double> xValHi;
  std::vector<double> yVal;
  std::vector<double> yErrStat;

  hepPlainTxtReader();
  hepPlainTxtReader(const std::string inFileName, const std::string histName);

  void Print();
  void Reset();
  void Init(const std::string inFileName, const std::string histName);
  void GetHistogram(TH1D** inHist_p, const std::string histName);
};

hepPlainTxtReader::hepPlainTxtReader()
{
  isInit = false;
  return;
}

hepPlainTxtReader::hepPlainTxtReader(const std::string inFileName, const std::string histName)
{
  Reset();
  Init(inFileName, histName);
  return;
}

void hepPlainTxtReader::Print()
{
  if(!isInit){
    std::cout << "hepPlainTxtReader::Print called, but not initialized. return" << std::endl;
    return;
  }

  std::cout << "Printing..." << std::endl;
  for(unsigned int i = 0; i < histStrings.size(); ++i){
    std::cout << " " << histStrings.at(i) << std::endl;
  }
  std::cout << "Done printing." << std::endl;
  
  return;
}

void hepPlainTxtReader::Reset()
{
  isInit = false;
  histStrings.clear();
  xValLow.clear();
  xValHi.clear();
  yVal.clear();
  yErrStat.clear();

  return;
}

void hepPlainTxtReader::Init(const std::string inFileName, const std::string histName)
{
  Reset();
  isInit = false;

  std::ifstream file(inFileName.c_str());
  std::string tempStr;

  std::vector<std::string> validNames;

  bool isFound = false;
  while(std::getline(file, tempStr)){
    if(tempStr.size() == 0) continue;
    while(tempStr.substr(0, 1).find(" ") != std::string::npos){tempStr.replace(0,1,"");}
    if(tempStr.size() == 0) continue;
    if(tempStr.find("Path") == std::string::npos) continue;

    tempStr.replace(0, tempStr.find("Path:")+5, "");
    while(tempStr.substr(0,1).find(" ") != std::string::npos){tempStr.replace(0,1,"");}
    validNames.push_back(tempStr);

    if(tempStr.find(histName.c_str()) != std::string::npos){
      histStrings.push_back(tempStr);

      bool startX = false;

      while(std::getline(file, tempStr)){
	if(tempStr.size() == 0) break;
	while(tempStr.substr(0, 1).find(" ") != std::string::npos){tempStr.replace(0,1,"");}
	while(tempStr.substr(0, 1).find(" ") != std::string::npos){tempStr.replace(0,1,"");}
	if(tempStr.size() == 0) break;
	if(tempStr.find("Path") != std::string::npos) break;

	histStrings.push_back(tempStr);

	if(startX){
	  tempStr.replace(0,1,"");
	  //	  std::cout << "\'" << tempStr << "\'" << std::endl;
	  
	  std::vector<std::string> stringBreakDown;

	  while(tempStr.size() != 0 && tempStr.find("\t") != std::string::npos){
	    stringBreakDown.push_back(tempStr.substr(0, tempStr.find("\t")));
	    tempStr.replace(0, tempStr.find("\t")+1, "");
	    while(tempStr.substr(0,1).find("\t") != std::string::npos){tempStr.replace(0,1,"");}
	  }

	  if(tempStr.size() != 0) stringBreakDown.push_back(tempStr);
	
	  xValLow.push_back(std::stod(stringBreakDown.at(1)));
	  xValHi.push_back(std::stod(stringBreakDown.at(2)));
	  yVal.push_back(std::stod(stringBreakDown.at(3)));
	  yErrStat.push_back(std::stod(stringBreakDown.at(4)));
	}

	if(tempStr.find("xdesc") != std::string::npos) startX = true;
      }

      isFound = true;
    }

    if(isFound) break;
  }

  file.close();

  if(!isFound){
    std::cout << "WARNING IN HEPPLAINTXTREADER::INIT: Initialization did not find \'" << histName << "\' did not find \'" << inFileName << "\'." << std::endl;
    
    std::cout << " Please pick from valid names: " << std::endl;
    for(unsigned int fI = 0; fI < validNames.size(); ++fI){
      if(fI%5 == 0) std::cout << "  ";
      std::cout << "\'" <<  validNames.at(fI) << "\', " ;
      if(fI%5 == 4) std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << "return uninitialized" << std::endl;
    return;
  }

  validNames.clear();

  isInit = true;

  return;
}

void hepPlainTxtReader::GetHistogram(TH1D** inHist_p, const std::string histName)
{
  if(!isInit){
    std::cout << "hepPlainTxtReader::GETHISTOGRAM: called, but reader not initialized. return without histogram set" << std::endl;
    return;
  }

  if((*inHist_p) != NULL){
    std::cout << "WARNING IN HEPPLAINTXTREADER::GETHISTOGRAM: inHist_p must be NULL! return without histogram set." << std::endl;
    return;
  }

  Int_t nBins = xValLow.size();
  Double_t bins[nBins+1];
  for(unsigned int i = 0; i < xValLow.size(); ++i){
    bins[i] = xValLow.at(i);
  }
  bins[nBins] = xValHi.at(nBins-1);

  //  std::cout << "BINS: " << std::endl;
  //  for(Int_t bI = 0; bI < nBins+1; ++bI){
  //    std::cout << " " << bins[bI] << ",";
  //  }
  //  std::cout << std::endl;
  
  (*inHist_p) = new TH1D(histName.c_str(), "", nBins, bins);
  for(unsigned int i = 0; i < yVal.size(); ++i){
    (*inHist_p)->SetBinContent(i+1, yVal.at(i));
    (*inHist_p)->SetBinError(i+1, yErrStat.at(i));
  }

  return;
}

#endif
