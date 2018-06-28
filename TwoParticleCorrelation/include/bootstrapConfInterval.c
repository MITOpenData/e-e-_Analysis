//
//  bootstrapConfInterval.c
//  
//
//  Created by Anthony Badea on 6/25/18.
//

// C dependencies
#include <stdio.h>
#include <vector>
#include <array>


// Root Dependencies
#include "TH1F.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TFile.h"

// Local Dependencies
#include "RNGFromDist.C"


// Create nHists histograms with identical binning and boundaries as h
// Sample N points = same number of entries in h (numEntries)
// Perform fit on each histogram using functionalForm
// Create distribution of each parameter
// returns [95% above zero, Low,High] where Low,High are the lower and upper values of the 95% confidence interval about the fit parameter
float ** bootstrapConfInterval(TH1F *h, TF1 *f, int nHists)
{
    using namespace std;
    TH1::SetDefaultSumw2();
    // return between -pi/2 and 3pi/2
    RNGFromDist r = RNGFromDist(f);
    
    // Initialize Variables
    static const int numHists = nHists;
    static const int numEntries = h->GetEntries();
    static const float histMin = h->GetXaxis()->GetBinLowEdge(1);
    static const float histMax = h->GetXaxis()->GetBinLowEdge(h->GetNbinsX()+1);
    std::cout<<histMin<<" "<<histMax<<std::endl;
    static const int numBins = h->GetNbinsX();
    static const int numParams = f->GetNpar();
    std::cout<<"Maximum of input fit: "<<f->GetMaximum(histMin,histMax)<<std::endl;
    for (int nP = 0; nP < f->GetNpar(); nP++)
    {
        std::cout<<Form("Input fit parameters %d: ",nP)<<f->GetParameter(nP)<<std::endl;
    }
    static const int paramBins = 100;
    
    // Initialize bootstrap parameters
    // Percentiles to be used for boostrap method to give 95% confidence interval
    // compute the entry in the array of the sorted boostrap parameters to find the percentile
    Double_t xq[2] = {0.025,0.975};
    static const int deltaLow = xq[0]*nHists;
    static const int deltaHigh = xq[1]*nHists;
    // will store the intervals and be returned
    float** intervals = new float*[numParams];
    vector<float> bootstrapDiffsAboutZero[numParams];
    vector<float> bootstrapDiffsAboutFitParam[numParams];
    
    // Initialize Histograms
    TH1F *sampleHists[numHists];
    TF1 *sampleFits[numHists];
    TH1F *paramDists[numParams];
    TH1F *bootstrapDiffsAboutZero_h[numParams];
    TH1F *bootstrapDiffsAboutFitParam_h[numParams];
    
    // initialize parameter distributions +- 10 from parameter of original fit
    for (int nP = 0; nP<numParams; nP++)
    {
        paramDists[nP] = new TH1F(Form("par_%d",nP),Form("par_%d",nP), paramBins, f->GetParameter(nP)-10, f->GetParameter(nP)+10);
        bootstrapDiffsAboutZero_h[nP] = new TH1F(Form("bootstrapDiffsAboutZero_h_%d",nP),Form("bootstrapDiffsAboutZero_h_%d",nP), 1000,-0.5,0.5);
        bootstrapDiffsAboutFitParam_h[nP] = new TH1F(Form("bootstrapDiffsAboutParam_h_%d",nP),Form("bootstrapDiffsAboutParam_h_%d",nP), 1000,-0.5,0.5);
        
    }
    
    
    for (int nH = 0; nH<numHists; nH++)
    {
        // initialize sample histograms and fits
        sampleHists[nH] = new TH1F(Form("clone_%d",nH), "", numBins, histMin,histMax); // histMin, histMax);
        sampleFits[nH] = (TF1*) f->Clone(Form("fit_%d",nH));
        sampleFits[nH]->SetParameter(0,100);
        
        std::cout<<Form("Generating dataset %d",nH)<<std::endl;
        // fill same histograms
        for (int M = 0; M<numEntries; M++)
        {
            if(M%5000 == 0) std::cout<<Form("%d / %d",M,numEntries)<<std::endl;
            // generate random number between histMin and histMax following distribution given by input function
            sampleHists[nH]->Fill(r.getRand());
        }
        
        // perform fit
        sampleHists[nH]->Fit(Form("fit_%d",nH),"Q N 0");
        sampleHists[nH]->SetStats(0);
        
        std::cout<<"filling param hists"<<std::endl;
        // fill parameter distributions
        
        for (int nP = 0; nP<numParams; nP++)
        {
            paramDists[nP]->Fill(sampleFits[nH]->GetParameter(nP));
            bootstrapDiffsAboutZero_h[nP]->Fill(sampleFits[nH]->GetParameter(nP) - 0.0);
            bootstrapDiffsAboutZero[nP].push_back(sampleFits[nH]->GetParameter(nP) - 0.0);
            bootstrapDiffsAboutFitParam_h[nP]->Fill(sampleFits[nH]->GetParameter(nP) - f->GetParameter(nP));
            bootstrapDiffsAboutFitParam[nP].push_back(sampleFits[nH]->GetParameter(nP) - f->GetParameter(nP));
        }
    }
    
    // Find 95% confidence interval of each of the paramDists about above fit parameter. Find 95% upper bound above zero.
    //cout<<"Confidence intervals before passing to test()"<<endl;
    for (int nP = 0; nP<numParams; nP++)
    {
        sort(bootstrapDiffsAboutZero[nP].begin(),bootstrapDiffsAboutZero[nP].end());
        sort(bootstrapDiffsAboutFitParam[nP].begin(),bootstrapDiffsAboutFitParam[nP].end());
        //cout<<Form("v_%d: [%f,%f]",nP,bootstrapDiffsAboutFitParam[nP][deltaLow],bootstrapDiffsAboutFitParam[nP][deltaHigh])<<endl;
        intervals[nP] = new float[3];
        intervals[nP][0] = bootstrapDiffsAboutZero[nP][deltaHigh];
        intervals[nP][1] = bootstrapDiffsAboutFitParam[nP][deltaLow];
        intervals[nP][2] = bootstrapDiffsAboutFitParam[nP][deltaHigh];
    }
    
    TFile *fout = new TFile("testBOOTSTRAPCONF.root","recreate");
    for (int nH = 0; nH < numHists; nH++){
        sampleHists[nH]->Write("",TObject::kOverwrite);
        sampleFits[nH]->Write("",TObject::kOverwrite);
    }
    for (int nP = 0; nP<numParams; nP++){
        paramDists[nP]->Write("",TObject::kOverwrite);
        bootstrapDiffsAboutZero_h[nP]->Write("",TObject::kOverwrite);
        bootstrapDiffsAboutFitParam_h[nP]->Write("",TObject::kOverwrite);
    }
    
    // Clean up
    for (int nP = 0; nP<numParams; nP++){
        delete bootstrapDiffsAboutFitParam_h[nP];
        delete bootstrapDiffsAboutZero_h[nP];
        delete paramDists[nP];
    }
    for (int nH = 0; nH < numHists; nH++){ delete sampleHists[nH]; delete sampleFits[nH];}
    
    
    fout->Close();
    delete fout;
    
    return intervals;
}

// create gaussian distribution with 100 data points to test if sampleFit working properly
int test()
{
    using namespace std;
    
    TH1F *gausHist = new TH1F("gausHist","gausHist",100,-10,10);
    static const float mean = 0;
    static const float sigma = 5;
    static const int numPoints = 10000;
    
    TRandom *r0 = new TRandom();
    for (int i = 0; i < numPoints; i++)
    {
        gausHist->Fill(r0->Gaus(mean,sigma));
    }
    
    string functionalForm = "gaus(0)"; // [0]*exp(-0.5*((x-[1])/[2])**2)
    TF1 *f1 = new TF1("gaus",functionalForm.c_str(),-10,10);
    gausHist->Fit("gaus","Q N 0");
    f1->GetRandom(-10,10);
    float** confidenceIntervals = bootstrapConfInterval(gausHist,f1,10);
    
    cout<<"Confidence intervals"<<endl;
    
    for (int nP = 0; nP < f1->GetNpar(); nP++)
    {
        float lowConf = f1->GetParameter(nP) + confidenceIntervals[nP][0]; // plus is correct here because if the delta is negative then the number will be shifted down
        float highConf = f1->GetParameter(nP) + confidenceIntervals[nP][1]; // plus is also correct here because if the delta is positive then the number will be shifted up
        cout<<Form("v_%d: [%f,%f]",nP,f1->GetParameter(nP) + confidenceIntervals[nP][0],f1->GetParameter(nP)+ confidenceIntervals[nP][1])<<endl;
        // important: clean up memory
        delete [] confidenceIntervals[nP];
    }
    
    delete [] confidenceIntervals;
    
    return 0;
}













