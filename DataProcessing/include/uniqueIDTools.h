#ifndef UNIQUEIDTOOLS_H
#define UNIQUEIDTOOLS_H

#include <math.h>
#include <iostream>

class uniqueIDTools{
 public:
  static const unsigned long long eventPow = 0;
  static const unsigned long long runPow = 5;
  static const unsigned long long yearPow = 10;
  static const unsigned long long subDirPow = 12;
  static const unsigned long long processPow = 15;
  static const unsigned long long isMCPow = 17;

  uniqueIDTools(){};
  unsigned long long getID(bool isMC, int process,  int subDir, int year, int runNo, int eventNo);
};

unsigned long long uniqueIDTools::getID(bool isMC, int process, int subDir, int year, int runNo, int eventNo)
{
  unsigned long long outID;
  outID = (unsigned long long)eventNo*(unsigned long long)std::pow(10, eventPow);
  outID += (unsigned long long)runNo*(unsigned long long)std::pow(10, runPow);
  year == 2000 ? year = 0 : year -= 1900;
  outID += (unsigned long long)year*(unsigned long long)std::pow(10, yearPow);
  outID += (unsigned long long)subDir*(unsigned long long)std::pow(10, subDirPow);
  outID += (unsigned long long)process*(unsigned long long)std::pow(10, processPow);
  if(isMC) outID += (unsigned long long)std::pow(10, isMCPow);

  return outID;
}

#endif
