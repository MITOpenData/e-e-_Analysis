#ifndef SIGNALMIXTABLEREADER_h
#define SIGNALMIXTABLEREADER_h

#include <string>
#include <map>
#include <fstream>
#include <vector>

#include "TMath.h"

#include "DataProcessing/include/checkMakeDir.h"

class signalMixTableReader{
 public:
  bool isInit = false;

  signalMixTableReader(){};
  signalMixTableReader(const std::string inFileName);

  bool isNum(const std::string inStr);
  bool Init(const std::string inFileName);
  unsigned int getKey(int mult, double inAbsEta);
  double getSigOverMixFactor(int mult, double inAbsEta);
  double getSigOverMixNum(int mult, double inAbsEta);
  double getSigOverMixDenom(int mult, double inAbsEta);

  void Print();

  std::map<unsigned int, double> globalMap;
  std::map<unsigned int, double> globalMapNum;
  std::map<unsigned int, double> globalMapDenom;
  std::vector<int> multLowBounds;
  std::vector<int> multHiBounds;
  std::vector<double> absEtaLowBounds;
  std::vector<double> absEtaHiBounds;
};

signalMixTableReader::signalMixTableReader(const std::string inFileName)
{
  Init(inFileName);
  return;
}

bool signalMixTableReader::isNum(const std::string inStr)
{
  if(std::string("0123456789").find(inStr.substr(0,1)) != std::string::npos) return true;
  return false;
}

bool signalMixTableReader::Init(const std::string inFileName)
{
  if(!checkFile(inFileName)){
    std::cout << "SIGNALMIXTABLEREADER: Given inFileName \'" << inFileName << "\' is invalid. not initialized" << std::endl;
    return false;
  }

  std::ifstream file(inFileName.c_str());
  std::string tempStr;

  unsigned int pos = 0;

  while(std::getline(file, tempStr)){
    if(tempStr.size() == 0) continue;


    if(tempStr.find("TABLE") != std::string::npos){
      pos = 0;
      multLowBounds.push_back(std::stoi(tempStr.substr(tempStr.find(":")+1, tempStr.find("<"))));
      tempStr.replace(0, tempStr.find("<")+1, "");
      multHiBounds.push_back(std::stoi(tempStr.substr(tempStr.find("<")+1, tempStr.size())));
    }
    else if(isNum(tempStr)){
      std::string firstNumStr = tempStr.substr(0, tempStr.find(","));
      tempStr.replace(0, tempStr.find(",")+1, "");
      std::string secondNumStr = tempStr.substr(0, tempStr.find(",")); 
      tempStr.replace(0, tempStr.find(",")+1, "");
      std::string thirdNumStr = tempStr.substr(0, tempStr.find(","));
      tempStr.replace(0, tempStr.find(",")+1, "");
      std::string fourthNumStr = tempStr.substr(0, tempStr.find(","));

      if(multLowBounds.size() == 1){
	absEtaLowBounds.push_back(std::stod(firstNumStr));
	absEtaHiBounds.push_back(std::stod(secondNumStr));
      }

      unsigned int key = (multHiBounds.size()-1) + 10*(pos);
      ++pos;

      std::cout << " Setting " << thirdNumStr << "/" << fourthNumStr << std::endl;
      globalMap[key] = std::stod(thirdNumStr)/std::stod(fourthNumStr);
      globalMapNum[key] = std::stod(thirdNumStr);
      globalMapDenom[key] = std::stod(fourthNumStr);
    }
  }

  file.close();

  isInit = true;
  return true;
}


unsigned int signalMixTableReader::getKey(int mult, double absEta)
{
  if(!isInit){
    std::cout << "SIGNALMIXTABLEREADER: Error, not initialized properly, return 0" << std::endl;
    return 0;
  }

  absEta = TMath::Abs(absEta);

  unsigned int mpos = 0;
  for(unsigned int mI = 0; mI < multHiBounds.size(); ++mI){
    if(mult >= multLowBounds.at(mI) && mult < multHiBounds.at(mI)){
      mpos = mI;
      break;
    }
  }

  unsigned int aepos = 0;
  for(unsigned int aeI = 0; aeI < absEtaHiBounds.size(); ++aeI){
    if(absEta >= absEtaLowBounds.at(aeI) && absEta < absEtaHiBounds.at(aeI)){
      aepos = aeI;
      break;
    }
  }

  unsigned int key = mpos + 10*aepos;
  return key;
}

double signalMixTableReader::getSigOverMixFactor(int mult, double absEta)
{
  if(!isInit){
    std::cout << "SIGNALMIXTABLEREADER: Error, not initialized properly, return 0" << std::endl;
    return 0;
  }

  unsigned int key = getKey(mult, absEta);
  return globalMap[key];
}

double signalMixTableReader::getSigOverMixNum(int mult, double absEta)
{
  if(!isInit){
    std::cout << "SIGNALMIXTABLEREADER: Error, not initialized properly, return 0" << std::endl;
    return 0;
  }

  unsigned int key = getKey(mult, absEta);
  return globalMapNum[key];
}

double signalMixTableReader::getSigOverMixDenom(int mult, double absEta)
{
  if(!isInit){
    std::cout << "SIGNALMIXTABLEREADER: Error, not initialized properly, return 0" << std::endl;
    return 0;
  }

  unsigned int key = getKey(mult, absEta);
  return globalMapDenom[key];
}


void signalMixTableReader::Print()
{
  std::cout << "Print table..." << std::endl;
  for(std::map<unsigned int,double>::iterator it = globalMap.begin(); it != globalMap.end(); ++it){
    std::cout << " " << it->first << ", " << it->second << std::endl;
  }

  return;
}

#endif
