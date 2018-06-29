//
//  templateFit.c
//  
//
//  Created by Anthony Badea on 6/29/18.
//

// C dependencies
#include <stdio.h>

// Root Dependencies
#include "TH1F.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TCanvas.h"

TF1* templateFit(TH1F *highMult, TH1F *periphLowMult)
{
    // define lambda for the template distribution = peripheral (low multiplicity) histogram + flow contribution. pass in parameters by reference
    auto templateDist = [periphLowMult](double* x, double* pars){
        //const double* xx = x[0];
        
        // the fit parameters, i.e. the histogram weight
        const double F = pars[0];
        const double ridgePars[2] = {pars[1], pars[2]};
        // get content of the histogram for this point
        const double yPeriph = periphLowMult->GetBinContent(periphLowMult->GetXaxis()->FindFixBin(x[0])); // use fix bin because we do not want to extend the axis if the given point belongs to an under-/overflow bin. the input x's should be the proper ones and we want to know if they are not
        TF1* yRidge = new TF1("yRidge","[0]*(1+2*[1]*cos(2*x))");
        return F*yPeriph + yRidge->EvalPar(&x[0],ridgePars);
    };
    
    // create function for the peripheral and ridge distributions. create sum of them using lambda
    TF1* yTempl = new TF1("yTempl",templateDist,-TMath::Pi()/2.0,3*TMath::Pi()/2.0,3); // last param is number of free param used by function
    highMult->Fit("yTempl","L Q N 0"); // Poisson (not Gaussian) statistics so use L for loglikelihood
    
    return yTempl;
}

int test(const std::string inFileName)
{
    TFile *f = new TFile(inFileName.c_str(),"read");
    TH1F *periphLowMult = (TH1F*)f->Get("h_deltaphi0_0_4_10_0_0");
    TH1F *highMult_10_20 = (TH1F*)f->Get("h_deltaphi0_0_10_20_0_0");
    TH1F *highMult_20_30 = (TH1F*)f->Get("h_deltaphi0_0_20_30_0_0");
    TH1F *highMult_30_999 = (TH1F*)f->Get("h_deltaphi0_0_30_999_0_0");
    TH1F *highMult_35_999 = (TH1F*)f->Get("h_deltaphi0_0_35_999_0_0");
    
    // Finish here // 
    TF1 *periphLowMult_fit = new TF1("periphLowMult_fit","gaus(0)",-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
    periphLowMult->Fit("periphLowMult_fit");
    
    TF1 *yTempl_10_20 = templateFit(highMult_10_20,periphLowMult);
    TF1 *yTempl_20_30 = templateFit(highMult_20_30,periphLowMult);
    TF1 *yTempl_30_999 = templateFit(highMult_30_999,periphLowMult);
    TF1 *yTempl_35_999 = templateFit(highMult_35_999,periphLowMult);
    
    TCanvas *c = new TCanvas("c","c",500,500);
    c->Divide(2,2);
    c->cd(1); highMult_10_20->Draw(); yTempl_10_20->Draw("C same");
    c->cd(2); highMult_20_30->Draw(); yTempl_20_30->Draw("C same");
    c->cd(3); highMult_30_999->Draw(); yTempl_30_999->Draw("C same");
    c->cd(4); highMult_35_999->Draw(); yTempl_35_999->Draw("C same");
    
    std::cout<<Form("(10,20): F = %f , G = %f, v2 = %f",yTempl_10_20->GetParameter(0), yTempl_10_20->GetParameter(1), yTempl_10_20->GetParameter(2))<<std::endl;
    std::cout<<Form("(20,30): F = %f , G = %f, v2 = %f",yTempl_20_30->GetParameter(0), yTempl_20_30->GetParameter(1), yTempl_20_30->GetParameter(2))<<std::endl;
    std::cout<<Form("(30,999): F = %f , G = %f, v2 = %f",yTempl_30_999->GetParameter(0), yTempl_30_999->GetParameter(1), yTempl_30_999->GetParameter(2))<<std::endl;
    std::cout<<Form("(35,999): F = %f , G = %f, v2 = %f",yTempl_35_999->GetParameter(0), yTempl_35_999->GetParameter(1), yTempl_35_999->GetParameter(2))<<std::endl;
    
    return 0;
}


// Basic test example to show the process for a double gaussian
// Create a fit function from a C function
double double_gaussian(double* x, double* pars) {
    // Here I just create the histograms I would like to fit with weights.
    // You don't have to do this here, all that is needed is that the histograms
    // are visible in this function.
    static bool initialized = false;
    static TH1D h1("h1", "", 100, -10, 10);
    static TH1D h2("h2", "", 100, -10, 10);
    if (!initialized) {
        int n = 100000;
        for (int i=0; i<n; ++i) h1.Fill(gRandom->Gaus(-3));
        for (int i=0; i<n; ++i) h2.Fill(gRandom->Gaus(+3));
        h1.Scale(1./n);
        h2.Scale(1./n);
        initialized = true;
    }
    
    const double xx = x[0]; // use x[1] to get 2nd dimension, x[2] for 3rd ...
    // the fit parameters, i.e. the histogram weights
    const double w1 = pars[0];
    const double w2 = pars[1];
    
    // get content of the histograms for this point
    const double y1 = h1.GetBinContent(h1.GetXaxis()->FindFixBin(xx));
    const double y2 = h2.GetBinContent(h2.GetXaxis()->FindFixBin(xx));
    
    return w1*y1 + w2*y2;
}

int fit() {
    TH1D *h = new TH1D("h", "data", 100, -10, 10);
    for (int i=0; i<100; ++i) h->Fill(gRandom->Gaus(-3));
    for (int i=0; i<100; ++i) h->Fill(gRandom->Gaus(+3));
    
    // create a fit function from a C function
    TF1 *f = new TF1("double_gaussian", double_gaussian, -10, 10, 2);
    
    h->Fit("double_gaussian", "L Q N 0");
    TCanvas *c = new TCanvas("c","c",500,500);
    h->Draw("E1");
    f->Draw("same");
    
    return 0;
}
