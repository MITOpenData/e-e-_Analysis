#ifndef getLinBins_h
#define getLinBins_h

#include "TMath.h"
void getLinBins(const float lower, const float higher, const int nBins, double bins[])
{
    bins[0] = lower;
    bins[nBins] = higher;
    
    float interval = (bins[nBins] - bins[0])/nBins;
    
    for(int iter = 1; iter < nBins; iter++){
        bins[iter] = bins[0] + iter*interval;
    }
    
    return;
}


#endif
