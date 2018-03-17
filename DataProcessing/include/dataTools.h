#ifndef DATATOOLS_H
#define DATATOOLS_H

#include <math.h>

const int alephFloatDecSize = 3;

float reducedPrecision(float inVal, int newPrecision = alephFloatDecSize)
{
  int newValInt = inVal*std::pow(10, newPrecision);
  float outVal = newValInt;
  return outVal/std::pow(10, newPrecision);
}

#endif
