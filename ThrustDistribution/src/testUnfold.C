// Standard library
#include <math.h>
#include <iostream>
#include <fstream>

// ROOT Library
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TLine.h>
#include <TF1.h>
#include <TCut.h>
#include <TPad.h>
#include "TRandom.h"
#include <TLegend.h>
// Unfolding
#include "bayesianUnfold.h"

// Constants
#define CanvasSizeX 400
#define CanvasSizeY 400



Double_t smear (Double_t xt)
{
  Double_t xeff= 0.3 + (1.0-0.3)/20*(xt+10.0);  // efficiency
  Double_t x= gRandom->Rndm();
//  if (x>xeff) return -10;
  Double_t xsmear= gRandom->Gaus(-5,1);     // bias and smear
  return xt+xsmear;
}



void testUnfold()
{
    TH1D *hM = new TH1D("hM","",200,-50,50);
    TH1D *hMMC = new TH1D("hMMC","",200,-50,50);
    TH1D *hT = new TH1D("hT","",200,-50,50);
    TH2D *hR = new TH2D("hR","",200,-50,50,200,-50,50);
    
   
  // Train with a Breit-Wigner, mean 0.3 and width 2.5.
  for (Int_t i= 0; i<10000000; i++) {
    Double_t xt= gRandom->Rndm()*90-45;
//    Double_t xt= gRandom->Gaus(0.0, 8.0);
    Double_t x= smear (xt);
    hR->Fill(xt,x);
    hMMC->Fill(x);
  }

  for (Int_t i=0; i<1000000; i++) {
    Double_t xt= gRandom->Gaus (0.0, 8.0), x= smear (xt);
    hT->Fill(xt);
    hM->Fill(x);
  }
  
  
  
  
    // rescale response function
  for (Int_t m= 1; m<=hR->GetNbinsY(); m++)
  {
     double Measured = hM->GetBinContent(m);
     double MeasuredMC = hMMC->GetBinContent(m);
	
     double Ratio=0;
     if (MeasuredMC) {
        Ratio = Measured / MeasuredMC;
     }
     
     for (Int_t t= 1; t<=hR->GetNbinsX(); t++)
     {
	double var = hR->GetBinContent(t,m);
	double varErr = hR->GetBinError(t,m);
	hR->SetBinContent(t,m,var*Ratio);
	hR->SetBinError(t,m,varErr*Ratio);
     }  
  }
  
  for (int i=0;i<100;i++) {
     bayesianUnfold myUnfolding0(hR,hT,0);
     myUnfolding0.unfold(hM,5);
  
     TH1D *hU0 = (TH1D*)myUnfolding0.hPrior->Clone();
     TH1D *hT2 = (TH1D*) hR->ProjectionX();

     for (Int_t t= 1; t<=hR->GetNbinsX(); t++)
     {
        double Measured = hU0->GetBinContent(t);
        double MeasuredMC = hT2->GetBinContent(t);
	
        double Ratio=0;
        if (MeasuredMC) {
           Ratio = Measured / MeasuredMC;
        }
       
        for (Int_t m= 1; m<=hR->GetNbinsY(); m++)
        {
 	  double var = hR->GetBinContent(t,m);
   	   double varErr = hR->GetBinError(t,m);
   	   hR->SetBinContent(t,m,var*Ratio);
	   hR->SetBinError(t,m,varErr*Ratio);
        }  
     }
  }
  
  
  bayesianUnfold myUnfolding(hR,hT,0);
  myUnfolding.unfold(hM,5);
  
  TCanvas *c = new TCanvas("c","",600,600);

  myUnfolding.hResMatrix->Draw("col");
  
  TCanvas *c2 = new TCanvas("c2","",600,600);
//  c2->SetLogy();
  
  TH1D *hU = (TH1D*)myUnfolding.hPrior->Clone();
  hU->SetName("hU");
  TH1D *hRepM = (TH1D*)myUnfolding.hReproduced->Clone();
  hRepM->SetName("hRepM");
  hRepM->SetLineColor(6);
  hRepM->Draw();
  hU->Draw("same");
  hT->SetLineColor(4);
  hT->Draw("same");
  hM->SetLineColor(2);
  hM->Draw("same");
  
  TLegend *leg = new TLegend(0.68,0.7,0.93,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry("hT","Truth","l");
  leg->AddEntry("hM","Measured","l");
  leg->AddEntry("hU","Unfolded","l");
  leg->AddEntry("hRepM","Reproduced","l");
  leg->Draw();
  
  TCanvas *c4 = new TCanvas("c4","",CanvasSizeX,CanvasSizeY);
  TH1F *hRatioM = (TH1F*)hM->Clone();
  hRatioM->SetName("hRatioM");
  hRatioM->SetXTitle("Measured / Reproduced Measured");
  hRatioM->Divide(hRepM);
  hRatioM->Draw();

  TCanvas *c5 = new TCanvas("c5","",CanvasSizeX,CanvasSizeY);
  TH1F *hRatioU = (TH1F*)hU->Clone();
  hRatioU->Divide(hT);
  hRatioU->SetXTitle("Unfolded / Truth");
  hRatioU->Draw();
  
}
