#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <iostream>
#include <stdlib.h>
#include <TStyle.h>

void plot(int isBelle=1,int mult=40,int cuthigh=1){
  
    TFile*finput=new TFile(Form("ROOTfiles/myoutput_isBelle%d_minMult%d.root",isBelle,mult));
    finput->cd();
    TH2F*h_2D=(TH2F*)finput->Get("h_2D");
    TH2F*h_2Dmix=(TH2F*)finput->Get("h_2Dmix");
    TH2F*h_ratio=(TH2F*)finput->Get("h_ratio");
 
    TCanvas * c1 = new  TCanvas("c1","c1",600,600);
    c1->SetPhi(40);
    c1->SetTheta(60);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    c1->Range(0,0,1,1);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(0);
    c1->SetTickx(1);
    c1->SetTicky(1);
    c1->SetLeftMargin(0.17);
    c1->SetRightMargin(0.08);
    c1->SetTopMargin(0.1);
    c1->SetBottomMargin(0.15);
    c1->SetFrameLineColor(0);
    c1->SetFrameBorderMode(0);
 
    h_ratio->SetTitle(Form("BELLE e^{+}e^{-}, 0.1 < p_{T}< 4 GeV/c, N >= %d",mult));
    h_ratio->GetXaxis()->CenterTitle();
    h_ratio->GetYaxis()->CenterTitle();
    h_ratio->GetXaxis()->CenterTitle();
    h_ratio->GetZaxis()->SetTitle("R(#Delta#eta, #Delta#phi)");
    
   
    h_ratio->GetXaxis()->SetTitleOffset(1.4);
    h_ratio->GetXaxis()->SetTitleSize(0.06);
    h_ratio->GetXaxis()->SetLabelSize(0.06);

    h_ratio->GetYaxis()->SetTitleOffset(1.4);
    h_ratio->GetYaxis()->SetTitleSize(0.06);
    h_ratio->GetYaxis()->SetLabelSize(0.06);

    h_ratio->GetZaxis()->SetTitleOffset(1.2);
    h_ratio->GetZaxis()->SetTitleSize(0.06);
    h_ratio->GetZaxis()->SetLabelSize(0.04);
    if(cuthigh) h_ratio->GetZaxis()->SetRangeUser(-6,-3);
    //h_ratio->SetNdivisions(1,"Z");
    h_ratio->Draw("surf1FB");
    if (isBelle) c1->SaveAs(Form("Plots/canvasRidgeBelleMult%dCutHigh%d.pdf",mult,cuthigh));
    
    
    TH1D*h_deltaphi[3];
    for (int i=0;i<3;i++)h_deltaphi[i]=(TH1D*)finput->Get(Form("h_deltaphi%d",i));
    
    TCanvas * c2 = new TCanvas("c2","c2",500,500);
    c2->cd(1);
    h_deltaphi[2]->Draw("ep");
    c2->SaveAs(Form("Plots/canvasProjection_isBelle%d_mult%d.pdf",isBelle,mult));

    
    

}

/*

    
*/