#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <iostream>
#include <stdlib.h>
#include <TStyle.h>

void plot(int isBelle=1,int mult=30,int cuthigh=0){
  
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
    c1->SetLeftMargin(0.22);
    c1->SetRightMargin(0.08);
    c1->SetTopMargin(0.1);
    c1->SetBottomMargin(0.15);
    c1->SetFrameLineColor(0);
    c1->SetFrameBorderMode(0);
 
    h_ratio->SetTitle(Form("BELLE e^{+}e^{-}, 0.1 < p_{T}< 4 GeV/c, N >= %d",mult));
    h_ratio->GetXaxis()->CenterTitle();
    h_ratio->GetYaxis()->CenterTitle();
    h_ratio->GetXaxis()->CenterTitle();
    h_ratio->GetZaxis()->SetTitle(" #frac{1}{N_{trig}} #frac{d^{2}N^{pair}}{d#Delta#eta#Delta#phi}");
    h_ratio->GetXaxis()->SetTitleOffset(1.4);
    h_ratio->GetXaxis()->SetTitleSize(0.06);
    h_ratio->GetXaxis()->SetLabelSize(0.06);
    h_ratio->GetYaxis()->SetTitleOffset(1.4);
    h_ratio->GetYaxis()->SetTitleSize(0.06);
    h_ratio->GetYaxis()->SetLabelSize(0.06);
    h_ratio->GetZaxis()->SetTitleOffset(2.3);
    h_ratio->GetZaxis()->SetTitleSize(0.038);
    h_ratio->GetZaxis()->SetLabelSize(0.04);
    if(cuthigh) h_ratio->GetZaxis()->SetRangeUser(-6,-3);
    //h_ratio->SetNdivisions(1,"Z");
    h_ratio->Draw("surf1FB");
    if (isBelle) c1->SaveAs(Form("Plots/canvasRidgeBelleMult%dCutHigh%d.pdf",mult,cuthigh));
    
    
    TH1D*h_deltaphi[3];
    for (int i=0;i<3;i++)h_deltaphi[i]=(TH1D*)finput->Get(Form("h_deltaphi%d",i));
    
    
    
    h_deltaphi[0]->GetYaxis()->SetTitle(" #frac{1}{N_{trig}} #frac{dN^{pair}}{d#Delta#phi} (0<#Delta#eta<1)");
    h_deltaphi[1]->GetYaxis()->SetTitle(" #frac{1}{N_{trig}} #frac{dN^{pair}}{d#Delta#phi} (1<#Delta#eta<2)");
    h_deltaphi[2]->GetYaxis()->SetTitle(" #frac{1}{N_{trig}} #frac{dN^{pair}}{d#Delta#phi} (2<#Delta#eta<3)");
    
    for (int i=0;i<3;i++){
      h_deltaphi[i]->SetTitle("");
      h_deltaphi[i]->GetYaxis()->SetTitleOffset(1.9);
      h_deltaphi[i]->GetXaxis()->SetTitleOffset(1.);
      h_deltaphi[i]->GetYaxis()->SetTitleSize(0.05);
      h_deltaphi[i]->GetXaxis()->SetTitleSize(0.06);
      h_deltaphi[i]->GetXaxis()->SetLabelSize(0.055);
      h_deltaphi[i]->GetYaxis()->SetLabelSize(0.055);
      //h_deltaphi[i]->SetMinimum(-6.);
      //h_deltaphi[i]->SetMaximum(-3.);      
      h_deltaphi[i]->GetYaxis()->SetLabelSize(0.055);
      h_deltaphi[i]->Sumw2();
    }


    TCanvas * c0_1 = new TCanvas("c0_1","c0_1",600,500);
    c0_1->SetLeftMargin(0.23);
    c0_1->SetRightMargin(0.05);
    h_deltaphi[0]->Draw("pe");
    c0_1->SaveAs(Form("Plots/canvasProjection_isBelle%d_mult%d_eta01.pdf",isBelle,mult));
    TCanvas * c1_2 = new TCanvas("c1_2","c1_2",600,500);
    c1_2->SetLeftMargin(0.23);
    c1_2->SetRightMargin(0.05);
    h_deltaphi[1]->Draw("pe");
    c1_2->SaveAs(Form("Plots/canvasProjection_isBelle%d_mult%d_eta12.pdf",isBelle,mult));
    TCanvas * c2_3 = new TCanvas("c2_3","c2_3",600,500);
    c2_3->SetLeftMargin(0.23);
    c2_3->SetRightMargin(0.05);
    h_deltaphi[2]->Draw("pe");
    c2_3->SaveAs(Form("Plots/canvasProjection_isBelle%d_mult%d_eta23.pdf",isBelle,mult));

}

/*

    
*/