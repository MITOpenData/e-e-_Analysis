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

// add a linear component p[0] +p[1]**xx
// try the unbinned fit and see what happens (usually safer)
// get uncertainity of the coefficients and the error and the relative error
TH1D *background;// = analysis(0, 0,  0, 30, 40, 50, 0, 1);
TF1 *total;
// = analysis(0, 0,  0, 0, 10, 50, 0, 1);

Double_t ftotal(Double_t *x, Double_t *par) {
    Double_t xx = x[0];
    //Int_t bin = background->GetXaxis()->FindBin(xx);
    //Int_t bin2 = background->GetXaxis()->FindBin(0.);
    Double_t br = par[0]*total->Eval(xx);
    //Double_t arg = (2*xx);
    Double_t sr = par[1] + 2*par[2]*TMath::Cos(2*xx);
    return br + sr;
}

Double_t cosine(Double_t *x1, Double_t *par1){
    Double_t xx = x1[0];
    Double_t y1 = par1[0] + 2*par1[1]*TMath::Cos(2*xx);
    return y1;
}

Double_t background_fit(Double_t *x, Double_t *par){
    Double_t xx = x[0];
    Double_t y = par[0] + par[1]*TMath::Cos(par[2]*xx);
    return y;
}

Double_t double_gauss_poly3(Double_t *x, Double_t *par)
{
    Double_t g1 = par[0]*TMath::Exp(-0.5*TMath::Power(((x[0]-PI)/par[1]),2));
    Double_t g2 = par[2]*TMath::Exp(-0.5*TMath::Power(((x[0]-PI)/par[3]),2));
    //Double_t poly3 = par[5]*TMath::Power(x[0],3) + par[6]*TMath::Power(x[0],2) + par[7]*x[0] + par[8];
    return g1 + g2+par[4]*x[0]+par[5];
}

void subtract()
{
    TH1::SetDefaultSumw2();
    //analysis(0, 0, 0, 0, 30, 50, 0, 1);
    TFile *file = new TFile("correlation_1_0_8.root");
    background = (TH1D*)file->Get("h_deltaphi_0");
    /*
    TF1 *g1 = new TF1("b1","gaus",-3.1416/2.,0.8);
    TF1 *g2 = new TF1("b2","gaus",0.8,2.);
    TF1 *g3 = new TF1("b3","gaus",2.,3.1416*1.5);

    // The total is the sum of the three, each has 3 parameters
    total = new TF1("mstotal","gaus(0)+gaus(3)+gaus(6)",-3.1416/2., 3.1416*1.5);
    
    Double_t par[9];
    background->Fit(g1,"R,0");
    background->Fit(g2,"R+,0");
    background->Fit(g3,"R+,0");
    g1->GetParameters(&par[0]);
    g2->GetParameters(&par[3]);
    g3->GetParameters(&par[6]);
    total->SetParameters(par);
    background->Fit(total,"R+");
    */
    total = new TF1("g2p3",double_gauss_poly3, -3.1416/2., 3.1416*1.5, 6);
    total->SetParameters(0,0.0003);
    total->SetParameters(1,0.5);
    total->SetParLimits(1,0.3,1);
    total->SetParameters(2,0.0003);
    total->SetParameters(3,0.1);
    total->SetParLimits(3,0.1,0.15);
    total->SetParameters(4,0);
    total->SetParameters(5,0);
    /*
    total->SetParameters(5,0);
    total->SetParameters(6,0);
    total->SetParameters(7,0);
    total->SetParameters(8,0);*/
    background->Fit("g2p3","MR");
    
    //analysis(0, 0, 0, 140, 170, 50, 0, 1);
    TFile *file2 = new TFile("correlation_1_8_20.root");
    TH1D* result = (TH1D*)file2->Get("h_deltaphi_0");
    
    
    
    TF1 *ftot = new TF1("ftot",ftotal,-3.1416/2., 3.1416*1.5, 3); //3 specifies the number of parameters we have to fit.
    Double_t norm = result->GetMaximum();
    
    ftot->SetParameter(0,0.1);
    ftot->SetParameter(1,0.5*norm);//
    ftot->SetParameter(2, 0.5*norm);
 
    
    
    result->Fit("ftot", "b");
    
    double p0 = ftot->GetParameter(0);
    double p1 = ftot->GetParameter(1);
    double p2 = ftot->GetParameter(2);
    //double p3 = ftot->GetParameter(3);
    //Int_t b = background->GetXaxis()->FindBin(0.);
    //double p = background->GetBinContent(b);
    
    double chi2 = ftot->GetChisquare();
    
    double Y_area = result->Integral(-3.1416/2., 3.1416*1.5);
    double Y_temp_area = ftot->Integral(-3.1416/2., 3.1416*1.5);
    double G_dash = Y_area/Y_temp_area;
    double G = G_dash*p1;
    double v_22 = p2/p1;
    
    //background->Scale(G_dash*p0);
    /*
     for (int i = 0; i<50; i++)
     {
     background->AddBinContent(i,p1);
     }
     */
    TF1 *cosine_effect = new TF1("cosine_effect", cosine, -3.1416/2., 3.1416*1.5, 2);
    cosine_effect->SetParameter(0, p1);
    cosine_effect->SetParameter(1, G*p2/p1);
    
    //gROOT->ForceStyle();
    TCanvas *c1 = new TCanvas("c1","",600,600);
    TCanvas *c2 = new TCanvas("c2","",600,600);
    
    result->SetMarkerStyle(4);
    background->SetMarkerStyle(20);
    result->SetMarkerColor(4);
    background->SetMarkerColor(2);
    cosine_effect->SetMarkerStyle(4);
    cosine_effect->SetMarkerColor(5);
    
    
    c1->cd();
    
    result->Draw();
    
    //background->Draw("same");
    
    cosine_effect->Draw("same");
    
    c1->SaveAs("fit_nbin50.pdf");
    
    cout<<"Y_area: "<<Y_area<<endl;
    cout<<"Y_temp_area: "<<Y_temp_area<<endl;
    cout<<"ChiSquare: "<<chi2<<endl;
    cout<<"G: "<<G<<endl;
    cout<<"G_dash: "<<G_dash<<endl;
    
    c2->cd();
    background->Draw();
    c2->SaveAs("fit_background_nbin50.pdf");
    
    file->Close();
    file2->Close();
}
