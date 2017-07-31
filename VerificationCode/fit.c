#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>
#include <stdlib.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TFormula.h>
//#include "fourier.h"
using namespace std;

#define PI 3.1415926
//enum SIMPLEPID {PHOTON, ELECTRON, PION, MUON, KAON, PROTON};
enum SIMPLEPWLAG {CHARGED_TRACK, CHARGED_LEPTONS1, CHARGED_LEPTONS2, V0, PHOTON, NEUTRAL_HADRON};

// Try to make the template composed of histograms. Fit the two distributions first. Fit the background first. Instead of putting the background elements in ftotal. first fit the background and then use the function generated from the fit in the ftotal.

// Check to make sure that the cosine fit is properly done because the blue curve has a minimum below the max of the data. if there was a signal then we should see a maximum there

// Try to superimpose the background plot to the original plot and see if the final result makes sense. Plot the background, cosine, fit, and actual to see if the actual - background = cosine and so on. See if the plots add up to the correct thing.

// Think about how to improve template overall. Talk to camelia. 
TH1D *background;// = analysis(0, 0,  0, 30, 40, 50, 0, 1);

// = analysis(0, 0,  0, 0, 10, 50, 0, 1);

Double_t ftotal(Double_t *x, Double_t *par) {
    Double_t xx = x[0];
    Int_t bin = background->GetXaxis()->FindBin(xx);
    //Int_t bin2 = background->GetXaxis()->FindBin(0.);
    Double_t br = par[0]*background->GetBinContent(bin);
    //Double_t arg = (2*xx);
    Double_t sr = par[1] + 2*par[2]*TMath::Cos(2*xx);
    return br + sr;
}

Double_t cosine(Double_t *x1, Double_t *par1){
    Double_t xx = x1[0];
    Double_t y1 = par1[0] + 2*par1[1]*TMath::Cos(2*xx);
    return y1;
}

void subtract()
{
    // Get the background and the result data
    TFile *file = new TFile("correlation_0_10.root");
    background = (TH1D*)file->Get("h_deltaphi_0");
    
    TFile *file2 = new TFile("correlation_30_40.root");
    TH1D* result = (TH1D*)file2->Get("h_deltaphi_0");
    
    // Prepare the function as described in the paper to be fitted
    TF1 *ftot = new TF1("ftot",ftotal,-3.1416/2., 3.1416*1.5, 3); //3 specifies the number of parameters we have to fit.
    Double_t norm = result->GetMaximum();
    //Set the initial parameters
    ftot->SetParameters(0.001, 0.5*norm, 0.5*norm);
    
    // Perform the fit
    result->Fit("ftot", "b");
    
    // Get the parameters so that we can plot the background and cosine parts seperately from the fit
    double p0 = ftot->GetParameter(0);
    double p1 = ftot->GetParameter(1);
    double p2 = ftot->GetParameter(2);
    
    double chi2 = ftot->GetChisquare();
    
    // Normalize so that the fit is of the form F*Background + G(1+v_2 cos(2x)) --> Find G here
    double Y_area = result->Integral(-3.1416/2., 3.1416*1.5);
    double Y_temp_area = ftot->Integral(-3.1416/2., 3.1416*1.5);
    double G_dash = Y_area/Y_temp_area;
    double G = G_dash*p1;
    double v_22 = p2/p1;
    
    // Plot the background and cosine portions seperately from the fit
    // MAKE A COPY OF THE BACKGROUND AND THEN DO IT
    //background->Scale(G_dash*p0);

    TF1 *cosine_effect = new TF1("cosine_effect", cosine, -3.1416/2., 3.1416*1.5, 2);
    cosine_effect->SetParameter(0, p1);
    cosine_effect->SetParameter(1, G*p2/p1);
    
    gROOT->ForceStyle();
    TCanvas *c1 = new TCanvas("c1","",600,600);
    
    result->SetMarkerStyle(4);
    background->SetMarkerStyle(20);
    result->SetMarkerColor(4);
    background->SetMarkerColor(2);
    
    
    
    c1->cd();
    
    result->Draw("C");
    
    //background->Draw("same");
    
    //cosine_effect->Draw("same");
    
    c1->SaveAs("fit.pdf");
    
    cout<<"Y_area: "<<Y_area<<endl;
    cout<<"Y_temp_area: "<<Y_temp_area<<endl;
    cout<<"ChiSquare: "<<chi2<<endl;
    cout<<"G: "<<G<<endl;
    cout<<"G_dash: "<<G_dash<<endl;
    
    
    file->Close();
    file2->Close();
}
