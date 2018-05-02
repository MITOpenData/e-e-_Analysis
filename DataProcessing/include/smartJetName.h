#ifndef SMARTJETNAME_H
#define SMARTJETNAME_H

#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "TFile.h"

#include "DataProcessing/include/returnRootFileContentsList.h"

bool sameString(const std::string inStr1, const std::string inStr2){return (inStr1.size() == inStr2.size() && inStr1.find(inStr2) != std::string::npos);}

bool stringInVect(const std::string inStr, std::vector<std::string> inVect)
{
  bool isIn = false;
  for(unsigned int i = 0; i < inVect.size(); ++i){
    if(sameString(inVect.at(i), inStr)){
      isIn = true;
      break;
    }
  }
  return isIn;
}


std::string smartJetName(const std::string inTreeGuess, std::vector<std::string> fileStr)
{
  std::string outStr = "";

  for(unsigned int i = 0; i < fileStr.size(); ++i){
    if(sameString(fileStr.at(i), inTreeGuess)){
      outStr = inTreeGuess;
      break;
    }
  }

  if(outStr.size() == 0){
    std::map<std::string, std::vector<std::string>> guessToAvailable;
    guessToAvailable["ak4JetTree"] = {"ak4ESchemeJetTree", "akR4ESchemeJetTree"};
    guessToAvailable["ak8JetTree"] = {"ak8ESchemeJetTree", "akR8ESchemeJetTree"};
    guessToAvailable["ak4ESchemeJetTree"] = {"akR4ESchemeJetTree", "ak4JetTree"};
    guessToAvailable["ak4WTAmodpSchemeJetTree"] = {"akR4WTAmodpSchemeJetTree"};
    guessToAvailable["ak8ESchemeJetTree"] = {"akR8ESchemeJetTree", "ak8JetTree"};
    guessToAvailable["ak8WTAmodpSchemeJetTree"] = {"akR8WTAmodpSchemeJetTree"};
    guessToAvailable["akR4ESchemeJetTree"] = {"ak4ESchemeJetTree", "ak4JetTree"};
    guessToAvailable["akR4WTAmodpSchemeJetTree"] = {"ak4WTAmodpSchemeJetTree"};
    guessToAvailable["akR8ESchemeJetTree"] = {"ak8ESchemeJetTree", "ak8JetTree"};
    guessToAvailable["akR8WTAmodpSchemeJetTree"] = {"ak8WTAmodpSchemeJetTree"};

    if(guessToAvailable.count(inTreeGuess) != 0){
      std::vector<std::string> available = guessToAvailable[inTreeGuess];

      for(unsigned int i = 0; i < available.size(); ++i){
	if(stringInVect(available.at(i), fileStr)){
	  outStr = available.at(i);
	  break;
	}
      }
    }
  }

  if(outStr.size() == 0){
    std::cout << "WARNING: Input jet tree guess \'" << inTreeGuess << "\' is not valid for given file, return empty string." << std::endl;
  }
  else if(outStr.size() != inTreeGuess.size() || outStr.find(inTreeGuess) == std::string::npos){
    std::cout << "WARNNING: Input jet tree guess \'" << inTreeGuess << "\' is not found in file, substituting \'" << outStr << "\'. Please fix if this substitution is not wanted." << std::endl;
  }

  return outStr;
}

std::string smartJetName(const std::string inTreeGuess, TFile* inFile_p)
{
  std::vector<std::string> fileStr = returnRootFileContentsList(inFile_p, "TTree", "Jet");
  unsigned int pos = 0;
  while(pos < fileStr.size()){
    bool isGood = true;
    for(unsigned int i = pos+1; i < fileStr.size(); ++i){
      if(sameString(fileStr.at(i), fileStr.at(pos))){
        fileStr.erase(fileStr.begin()+i);
        isGood = false;
        break;
      }
    }
    if(isGood) pos++;
  }
  
  return smartJetName(inTreeGuess, fileStr);
}

#endif
