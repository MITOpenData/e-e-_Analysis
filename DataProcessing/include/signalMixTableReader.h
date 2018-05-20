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
  bool isInit_ = false;
  bool is2D_ = false;

  signalMixTableReader(){};
  signalMixTableReader(const std::string inFileName);

  bool isNum(const std::string inStr);
  bool Init(const std::string inFileName);
  unsigned int getKey(int mult, double inAbsEta, double phi);
  double getSigOverMixFactor(int mult, double inAbsEta);
  double getSigOverMixNum(int mult, double inAbsEta);
  double getSigOverMixDenom(int mult, double inAbsEta);
  double getSigOverMixFactor2D(int mult, double inAbsEta, double phi);
  double getSigOverMixNum2D(int mult, double inAbsEta, double phi);
  double getSigOverMixDenom2D(int mult, double inAbsEta, double phi);

  bool is2D(){return is2D_;}

  void Print();

  std::map<unsigned int, double> globalMap;
  std::map<unsigned int, double> globalMapNum;
  std::map<unsigned int, double> globalMapDenom;
  std::vector<int> multLowBounds;
  std::vector<int> multHiBounds;
  std::vector<double> absEtaLowBounds;
  std::vector<double> absEtaHiBounds;
  std::vector<double> phiLowBounds;
  std::vector<double> phiHiBounds;
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

  std::getline(file, tempStr);
  std::getline(file, tempStr);

  int lines = 0;
  while(tempStr.find(",") != std::string::npos){
    tempStr.replace(0, tempStr.find(",") + 1, "");
    ++lines;
  }

  if(lines == 8) is2D_ = true;
  
  file.close();
  file.open(inFileName.c_str());

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

      std::string fifthNumStr = "";
      std::string sixthNumStr = "";      
      if(is2D_){
	fifthNumStr = thirdNumStr;
	sixthNumStr = fourthNumStr;

	tempStr.replace(0, tempStr.find(",")+1, "");
	thirdNumStr = tempStr.substr(0, tempStr.find(","));
	tempStr.replace(0, tempStr.find(",")+1, "");
	fourthNumStr = tempStr.substr(0, tempStr.find(","));
      }

      if(multLowBounds.size() == 1){
	absEtaLowBounds.push_back(std::stod(firstNumStr));
	absEtaHiBounds.push_back(std::stod(secondNumStr));

	if(is2D_){
	  phiLowBounds.push_back(std::stod(fifthNumStr));
	  phiHiBounds.push_back(std::stod(sixthNumStr));
	}
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

  isInit_ = true;
  return true;
}


unsigned int signalMixTableReader::getKey(int mult, double absEta, double phi = 0.)
{
  if(!isInit_){
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
    bool isAE = absEta >= absEtaLowBounds.at(aeI) && absEta < absEtaHiBounds.at(aeI);
    bool isPhi = true;
    if(is2D_) isPhi = phi >= phiLowBounds.at(aeI) && phi < phiHiBounds.at(aeI);

    if(isAE && isPhi){
      aepos = aeI;
      break;
    }
  }

  unsigned int key = mpos + 10*aepos;
  return key;
}

double signalMixTableReader::getSigOverMixFactor(int mult, double absEta)
{
  if(!isInit_){
    std::cout << "SIGNALMIXTABLEREADER: Error, not initialized properly, return 0" << std::endl;
    return 0;
  }
  if(is2D_){
    std::cout << "SIGNALMIXTABLEREADER: Error, call of getSigOverMixFactor() will not work when input table is 2d, please use getSigOverMixFactor2D(), return 0" << std::endl;
    return 0;
  }

  unsigned int key = getKey(mult, absEta);
  return globalMap[key];
}

double signalMixTableReader::getSigOverMixNum(int mult, double absEta)
{
  if(!isInit_){
    std::cout << "SIGNALMIXTABLEREADER: Error, not initialized properly, return 0" << std::endl;
    return 0;
  }
  if(is2D_){
    std::cout << "SIGNALMIXTABLEREADER: Error, call of getSigOverMixNum() will not work when input table is 2d, please use getSigOverMixNum2D(), return 0" << std::endl;
    return 0;
  }

  unsigned int key = getKey(mult, absEta);
  return globalMapNum[key];
}

double signalMixTableReader::getSigOverMixDenom(int mult, double absEta)
{
  if(!isInit_){
    std::cout << "SIGNALMIXTABLEREADER: Error, not initialized properly, return 0" << std::endl;
    return 0;
  }
  if(is2D_){
    std::cout << "SIGNALMIXTABLEREADER: Error, call of getSigOverMixDenom() will not work when input table is 2d, please use getSigOverMixDenom2D(), return 0" << std::endl;
    return 0;
  }

  unsigned int key = getKey(mult, absEta);
  return globalMapDenom[key];
}


double signalMixTableReader::getSigOverMixFactor2D(int mult, double absEta, double phi)
{
  if(!isInit_){
    std::cout << "SIGNALMIXTABLEREADER: Error, not initialized properly, return 0" << std::endl;
    return 0;
  }
  if(!is2D_){
    std::cout << "SIGNALMIXTABLEREADER: Error, call of getSigOverMixFactor2D() will not work when input table is 1d, please use getSigOverMixFactor(), return 0" << std::endl;
    return 0;
  }


  unsigned int key = getKey(mult, absEta, phi);
  return globalMap[key];
}

double signalMixTableReader::getSigOverMixNum2D(int mult, double absEta, double phi)
{
  if(!isInit_){
    std::cout << "SIGNALMIXTABLEREADER: Error, not initialized properly, return 0" << std::endl;
    return 0;
  }
  if(!is2D_){
    std::cout << "SIGNALMIXTABLEREADER: Error, call of getSigOverMixNum2D() will not work when input table is 1d, please use getSigOverMixNum(), return 0" << std::endl;
    return 0;
  }

  unsigned int key = getKey(mult, absEta, phi);
  return globalMapNum[key];
}

double signalMixTableReader::getSigOverMixDenom2D(int mult, double absEta, double phi)
{
  if(!isInit_){
    std::cout << "SIGNALMIXTABLEREADER: Error, not initialized properly, return 0" << std::endl;
    return 0;
  }
  if(!is2D_){
    std::cout << "SIGNALMIXTABLEREADER: Error, call of getSigOverMixDenom2D() will not work when input table is 1d, please use getSigOverMixDenom(), return 0" << std::endl;
    return 0;
  }


  unsigned int key = getKey(mult, absEta, phi);
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
