//
//  DQC_Settings
//  
//
//  Created by Anthony Badea on 1/22/18.
//

#ifndef DQC_SETTINGS_H
#define DQC_SETTINGS_H

//c and c++ dependencies
#include <sstream>
#include <string>
#include <iostream>
#include <stdlib.h>

//root dependencies
#include <TMath.h>
#include <TString.h>
#define PI 3.14159265358979

//local headers
#include "getLogBins.h"
#include "getLinBins.h"

class DQC_Settings
{
    public:
        // Data sample independent binning
    
        // declare log binning for probability of pt/multiplicity
        static const int nBinsP = 80;
        Double_t binsPLog[nBinsP+1]; // 80 is the number of bins
        Double_t logLow = .0000005;
        Double_t logHi = .3;
    
        // declare linear binning for probability of eta and phi
        static const int nBinsL = 80;
        Double_t binsPLin[nBinsL+1]; // 80 is the number of bins
        Double_t linLow = 0;
        Double_t linHi = 0.15;
    
        // declare binning for ratio plots
        static const int nBinsR = 80;
        Double_t binsRLin[nBinsR+1]; // 80 is the number of bins
        Double_t rLow = 0;
        Double_t rHi = 2.0;
    
        // Data sample dependent binning
    
        // declare binning for x-axis
        Int_t nBinsMult;
        Double_t multHi;
        Double_t multLow;
        std::vector<Double_t> binsMult;
        // to convert the vector to a array do std::vector<double> v; double* a = &v[0];
        Int_t nBinsPt;
        Double_t ptHi;
        Double_t ptLow;
        std::vector<Double_t> binsPt;
    
        Int_t nBinsEta;
        Double_t etaHi;
        Double_t etaLow;
        std::vector<Double_t> binsEta;
    
        Int_t nBinsPhi;
        Double_t phiHi;
        Double_t phiLow;
        std::vector<Double_t> binsPhi;
    
        Int_t nBinsY;
        Double_t yHi;
        Double_t yLow;
        std::vector<Double_t> binsY;
    
        DQC_Settings(const TString datalabel);
        int eeplot_Settings( const TString datalabel );
        int genBinning();
    
    private:
};

//int DQC_Settings::nBinsP = 80;

DQC_Settings::DQC_Settings(const TString datalabel)
{
    std::cout << "Getting settings.." << std::endl;
    
    getLogBins(logLow, logHi, nBinsP, binsPLog);
    getLinBins(linLow,linHi,nBinsL,binsPLin);
    getLinBins(rLow, rHi, nBinsR, binsRLin);
    
    eeplot_Settings(datalabel);
    genBinning();
    return;
}

int DQC_Settings::eeplot_Settings( const TString datalabel )
{
    if(datalabel == "PYTHIA8")
    {
        nBinsMult = 100; multLow = 0; multHi = 100;
        nBinsPt = 100; ptLow = 0; ptHi = 60;
        nBinsEta = 100; etaLow = -3.0; etaHi = 3.0;
        nBinsPhi = 100; phiLow = -TMath::Pi()/2.; phiHi = TMath::Pi();
        nBinsY = 100; yLow = 0; yHi = 8;
    }
    if(datalabel == "LEP1")
    {
        nBinsMult = 80; multLow = 0; multHi = 80;
        nBinsPt = 100; ptLow = 0; ptHi = 60;
        nBinsEta = 100; etaLow = -2.5; etaHi = 2.5;
        nBinsPhi = 100; phiLow = -TMath::Pi()/2.; phiHi = TMath::Pi();
        nBinsY = 100; yLow = 0; yHi = 8;
    }
    if(datalabel == "LEP2")
    {
        nBinsMult = 90; multLow = 0; multHi = 90;
        nBinsPt = 100; ptLow = 0; ptHi = 60;
        nBinsEta = 100; etaLow = -1.75; etaHi = 1.75;
        nBinsPhi = 100; phiLow = -TMath::Pi()/2.; phiHi = TMath::Pi();
        nBinsY = 100; yLow = 0; yHi = 8;
    }
    
    return 0;
}

// generate the binning based on data sample
int DQC_Settings::genBinning()
{
    Double_t binsmult[nBinsMult+1];
    Double_t binspt[nBinsPt+1];
    Double_t binseta[nBinsEta+1];
    Double_t binsphi[nBinsPhi+1];
    Double_t binsy[nBinsY+1];
    
    getLinBins(multLow, multHi, nBinsMult, binsmult);
    getLinBins(ptLow, ptHi, nBinsPt, binspt);
    getLinBins(etaLow, etaHi, nBinsEta, binseta);
    getLinBins(phiLow, phiHi, nBinsPhi, binsphi);
    getLinBins(yLow, yHi, nBinsY, binsy);
    
    // fill the vector bins
    binsMult.assign(binsmult, binsmult + nBinsMult + 1);
    binsPt.assign(binspt, binspt + nBinsPt + 1);
    binsEta.assign(binseta, binseta + nBinsEta + 1);
    binsPhi.assign(binsphi, binsphi + nBinsPhi + 1);
    binsY.assign(binsy, binsy + nBinsY + 1);
    
    return 1;
}
#endif /* DQC_Settings.h */
