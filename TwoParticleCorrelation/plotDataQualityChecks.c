//
//  plotDataQualityChecks.c
//  
//
//  Created by Anthony Badea on 1/11/18.
//

#include <iostream>
#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TDatime.h"
#include "TLine.h"
#include "TLegend.h"
#include "xjjrootuti.h"
int plotDataQualityChecks()
{
    xjjroot::setgstyle();
    
    // load data files
    TFile *LEP1 = new TFile("qualityCheck/outFile_LEP1.root");
    TH1F *LEP1_cmult = (TH1F*)gDirectory->Get("cmult");
    
    TH1F *LEP1_cmom = (TH1F*)gDirectory->Get("cmom");
    TH1F *LEP1_ceta = (TH1F*)gDirectory->Get("ceta");
    
    TFile *LEP2 = new TFile("qualityCheck/outFile_LEP2.root");
    TH1F *LEP2_cmult = (TH1F*)gDirectory->Get("cmult");
    TH1F *LEP2_cmom = (TH1F*)gDirectory->Get("cmom");
    TH1F *LEP2_ceta = (TH1F*)gDirectory->Get("ceta");
    
    TFile *PYTHIA8 = new TFile("qualityCheck/outFile_PYTHIA8.root");
    TH1F *PYTHIA8_cmult = (TH1F*)gDirectory->Get("cmult");
    TH1F *PYTHIA8_cmom = (TH1F*)gDirectory->Get("cmom");
    TH1F *PYTHIA8_ceta = (TH1F*)gDirectory->Get("ceta");
    
    // multiplicity spectrum
    TCanvas *mult_LEP1_PYTHIA8 = new TCanvas("mult","",600,600);
    TH2F *hempty_mult1 = new TH2F("Multiplicity Distribution",";Multiplicity;Probability",1,0,95,1,0,0.15);
    xjjroot::sethempty(hempty_mult1,0,0.3);
    hempty_mult1->Draw();
    xjjroot::setthgrstyle(LEP1_cmult, kRed, 21, 1.2, kRed, 1, 1, -1, -1, -1);
    LEP1_cmult->Draw("pe same");
    xjjroot::setthgrstyle(PYTHIA8_cmult, kBlue, 21, 1.2, kBlue, 1, 1, -1, -1, -1);
    PYTHIA8_cmult->Draw("pe same");
    TLegend *leg_mult_LEP1_PYTHIA8 = new TLegend(0.67,0.7,1.1,0.88);
    xjjroot::setleg(leg_mult_LEP1_PYTHIA8);
    leg_mult_LEP1_PYTHIA8->AddEntry(LEP1_cmult,"LEP1","p");
    leg_mult_LEP1_PYTHIA8->AddEntry(PYTHIA8_cmult,"PYTHIA8","p");
    leg_mult_LEP1_PYTHIA8->Draw();
    xjjroot::drawtex(0.2,0.876,"LEP1 vs PYTHIA8");
    mult_LEP1_PYTHIA8->SaveAs("qualityCheck/mult_LEP1_PYTHIA8.pdf");
    
    TCanvas *mult_LEP2 = new TCanvas("mult","",600,600);
    TH2F *hempty_mult2 = new TH2F("Multiplicity Distribution",";Multiplicity;Probability",1,0,95,1,0,0.15);
    xjjroot::sethempty(hempty_mult2,0,0.3);
    hempty_mult2->Draw();
    xjjroot::setthgrstyle(LEP2_cmult, kBlack, 21, 1.2, kBlack, 1, 1, -1, -1, -1);
    LEP2_cmult->Draw("pe same");
    TLegend *leg_mult_LEP2 = new TLegend(0.67,0.7,1.1,0.88);
    xjjroot::setleg(leg_mult_LEP2);
    leg_mult_LEP2->AddEntry(LEP2_cmult,"LEP2","p");
    leg_mult_LEP2->Draw();
    xjjroot::drawtex(0.2,0.876,"LEP2");
    mult_LEP2->SaveAs("qualityCheck/mult_LEP2.pdf");
    
    // pt spectra
    TCanvas *pt_LEP1_PYTHIA8 = new TCanvas("pt","",600,600);
    TH2F *hempty_pt1 = new TH2F("PT Spectra",";PT;Probability",1,0,95,1,0,0.15);
    xjjroot::sethempty(hempty_pt1,0,0.3);
    hempty_pt1->Draw();
    xjjroot::setthgrstyle(LEP1_cmom, kRed, 21, 1.2, kRed, 1, 1, -1, -1, -1);
    LEP1_cmom->Draw("pe same");
    xjjroot::setthgrstyle(PYTHIA8_cmom, kBlue, 21, 1.2, kBlue, 1, 1, -1, -1, -1);
    PYTHIA8_cmom->Draw("pe same");
    TLegend *leg_pt_LEP1_PYTHIA8 = new TLegend(0.67,0.7,1.1,0.88);
    xjjroot::setleg(leg_pt_LEP1_PYTHIA8);
    leg_pt_LEP1_PYTHIA8->AddEntry(LEP1_cmom,"LEP1","p");
    leg_pt_LEP1_PYTHIA8->AddEntry(PYTHIA8_cmom,"PYTHIA8","p");
    leg_pt_LEP1_PYTHIA8->Draw();
    xjjroot::drawtex(0.2,0.876,"LEP1 vs PYTHIA8");
    pt_LEP1_PYTHIA8->SaveAs("qualityCheck/pt_LEP1_PYTHIA8.pdf");
    
    TCanvas *pt_LEP2 = new TCanvas("pt","",600,600);
    TH2F *hempty_pt2 = new TH2F("PT Spectra",";PT;Probability",1,0,95,1,0,0.15);
    xjjroot::sethempty(hempty_pt2,0,0.3);
    hempty_pt2->Draw();
    xjjroot::setthgrstyle(LEP2_cmom, kBlack, 21, 1.2, kBlack, 1, 1, -1, -1, -1);
    LEP2_cmom->Draw("pe same");
    TLegend *leg_pt_LEP2 = new TLegend(0.67,0.7,1.1,0.88);
    xjjroot::setleg(leg_pt_LEP2);
    leg_pt_LEP2->AddEntry(LEP2_cmom,"LEP2","p");
    leg_pt_LEP2->Draw();
    xjjroot::drawtex(0.2,0.876,"LEP2");
    pt_LEP2->SaveAs("qualityCheck/pt_LEP2.pdf");
    
    // eta spectra
    TCanvas *eta_LEP1_PYTHIA8 = new TCanvas("eta","",600,600);
    TH2F *hempty_eta1 = new TH2F("Eta Spectra",";Eta;Probability",1,-5,5,1,0,0.06);
    xjjroot::sethempty(hempty_eta1,0,0.3);
    hempty_eta1->Draw();
    xjjroot::setthgrstyle(LEP1_ceta, kRed, 21, 1.2, kRed, 1, 1, -1, -1, -1);
    LEP1_ceta->Draw("pe same");
    xjjroot::setthgrstyle(PYTHIA8_ceta, kBlue, 21, 1.2, kBlue, 1, 1, -1, -1, -1);
    PYTHIA8_ceta->Draw("pe same");
    TLegend *leg_eta_LEP1_PYTHIA8 = new TLegend(0.67,0.7,1.1,0.88);
    xjjroot::setleg(leg_eta_LEP1_PYTHIA8);
    leg_eta_LEP1_PYTHIA8->AddEntry(LEP1_ceta,"LEP1","p");
    leg_eta_LEP1_PYTHIA8->AddEntry(PYTHIA8_ceta,"PYTHIA8","p");
    leg_eta_LEP1_PYTHIA8->Draw();
    xjjroot::drawtex(0.2,0.876,"LEP1 vs PYTHIA8");
    eta_LEP1_PYTHIA8->SaveAs("qualityCheck/eta_LEP1_PYTHIA8.pdf");
    
    TCanvas *eta_LEP2 = new TCanvas("eta","",600,600);
    TH2F *hempty_eta2 = new TH2F("Eta Spectra",";Eta;Probability",1,-5,5,1,0,0.05);
    xjjroot::sethempty(hempty_eta2,0,0.3);
    hempty_eta2->Draw();
    xjjroot::setthgrstyle(LEP2_ceta, kBlack, 21, 1.2, kBlack, 1, 1, -1, -1, -1);
    LEP2_ceta->Draw("pe same");
    TLegend *leg_eta_LEP2 = new TLegend(0.67,0.7,1.1,0.88);
    xjjroot::setleg(leg_eta_LEP2);
    leg_eta_LEP2->AddEntry(LEP2_ceta,"LEP2","p");
    leg_eta_LEP2->Draw();
    xjjroot::drawtex(0.2,0.876,"LEP2");
    eta_LEP2->SaveAs("qualityCheck/eta_LEP2.pdf");
    
    return 0;
}

int drawEta(TString filename = "qualityCheck/outFile_LEP1.root",  // input file
            TString datalabel = "LEP1"  // Text in upper-left corner
           )
{
    TFile *file = new TFile(filename.Data());
    TH1F *aeta = (TH1F*)gDirectory->Get("aeta");
    TH1F *neta = (TH1F*)gDirectory->Get("neta");
    TH1F *ceta = (TH1F*)gDirectory->Get("ceta");
    
    xjjroot::setgstyle();
    // Plot all, neutral, charged eta
    TCanvas *alleta = new TCanvas("alleta","alleta",600,600);
    TH2F *hemptyeta = new TH2F("",";eta;Probability",1,-5,5,1,0,0.06);
    xjjroot::sethempty(hemptyeta,0,0.3);
    hemptyeta->Draw();
    xjjroot::setthgrstyle(aeta, kBlack, 20, 1.2, kRed, 1, 1, -1, -1, -1);
    aeta->Draw("pe same");
    xjjroot::setthgrstyle(ceta, kRed, 21, 1.2, kRed, 1, 1, -1, -1, -1);
    ceta->Draw("pe same");
    xjjroot::setthgrstyle(neta, kBlue, 22, 1.2, kBlue, 1, 1, -1, -1, -1);
    neta->Draw("pe same");
    TLegend* etaleg = new TLegend(0.67,0.7,1.1,0.88);
    xjjroot::setleg(etaleg);
    etaleg->AddEntry(aeta,"All","p");
    etaleg->AddEntry(ceta,"Charged Hadrons","p");
    etaleg->AddEntry(neta,"Neutral Hadrons + Photons","p");
    etaleg->Draw();
    xjjroot::drawtex(0.2,0.876,datalabel);
    alleta->SaveAs(Form("qualityCheck/%s_eta.pdf",datalabel.Data()));
    return 0;
}



