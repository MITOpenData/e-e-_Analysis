//
//  plotDataQualityChecks.c
//  
//
//  Created by Anthony Badea on 1/11/18.
//

//c and c++ dependencies
#include <iostream>
#include <string>

//root dependencies
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TDatime.h"
#include "TLine.h"
#include "TLegend.h"
#include "TH2F.h"

//local headers
#include "../include/xjjrootuti.h"
#include "../include/getLogBins.h"
#include "../include/getLinBins.h"
int plotDataQualityChecks()
{
    // setting ROOT defaults
    xjjroot::setgstyle();
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    // declare log binning for probability of pt/multiplicity
    const int nBinsY = 80;
    Double_t binsYLog[nBinsY+1];
    const Double_t logLow = .0000005;
    const Double_t logHi = .3;
    getLogBins(logLow, logHi, nBinsY, binsYLog);
    
    // declare linear binning for probability of eta
    Double_t binsYLin[nBinsY+1];
    const Double_t linLow = 0;
    const Double_t linHi = 0.1;
    getLinBins(linLow,linHi,nBinsY,binsYLin);
    
    // declare binning for x-axis
    int nBinsMult = 100;
    Double_t binsMult[nBinsMult+1];
    Double_t multHi = 100;
    Double_t multLow = 0;
    getLinBins(multLow, multHi, nBinsMult, binsMult);
    
    int nBinsPt = 60;
    Double_t binsPt[nBinsPt+1];
    Double_t ptHi = 60;
    Double_t ptLow = 0;
    getLinBins(ptLow, ptHi, nBinsPt, binsPt);
    
    int nBinsEta = 60;
    Double_t binsEta[nBinsEta+1];
    Double_t etaHi = 3.0;
    Double_t etaLow = -etaHi;
    getLinBins(etaLow, etaHi, nBinsEta, binsEta);
    
    // load data files
    TFile *LEP1 = new TFile("inputs/qualityCheck/outFile_LEP1.root");
    TH1F *LEP1_cmult = (TH1F*)gDirectory->Get("cmult");
    TH1F *LEP1_cpt = (TH1F*)gDirectory->Get("cpt");
    TH1F *LEP1_ceta = (TH1F*)gDirectory->Get("ceta");
    
    /*
    TFile *LEP2 = new TFile("inputs/qualityCheck/outFile_LEP2.root");
    TH1F *LEP2_cmult = (TH1F*)gDirectory->Get("cmult");
    TH1F *LEP2_cpt = (TH1F*)gDirectory->Get("cpt");
    TH1F *LEP2_ceta = (TH1F*)gDirectory->Get("ceta");
    */
    
    TFile *PYTHIA8 = new TFile("inputs/qualityCheck/outFile_PYTHIA8.root");
    TH1F *PYTHIA8_cmult = (TH1F*)gDirectory->Get("cmult");
    TH1F *PYTHIA8_cpt = (TH1F*)gDirectory->Get("cpt");
    TH1F *PYTHIA8_ceta = (TH1F*)gDirectory->Get("ceta");
    
    // multiplicity spectrum
    TCanvas *mult_LEP1_PYTHIA8 = new TCanvas("mult","",600,600);
    gPad->SetLogy();
    TH2F *hempty_mult1 = new TH2F("Multiplicity Distribution",";Multiplicity;Probability",nBinsMult,binsMult,nBinsY,binsYLog);
    //TH2F *hempty_mult1 = new TH2F("Multiplicity Distribution",";Multiplicity;Probability",1,0,95,1,0,0.15);
    xjjroot::sethempty(hempty_mult1,0,0.3);
    hempty_mult1->DrawCopy();
    xjjroot::setthgrstyle(LEP1_cmult, kRed, 21, 1.2, kRed, 1, 1, -1, -1, -1);
    LEP1_cmult->DrawCopy("pe same");
    xjjroot::setthgrstyle(PYTHIA8_cmult, kBlue, 21, 1.2, kBlue, 1, 1, -1, -1, -1);
    PYTHIA8_cmult->DrawCopy("pe same");
    TLegend *leg_mult_LEP1_PYTHIA8 = new TLegend(0.67,0.7,1.1,0.88);
    xjjroot::setleg(leg_mult_LEP1_PYTHIA8);
    leg_mult_LEP1_PYTHIA8->AddEntry(LEP1_cmult,"LEP1","p");
    leg_mult_LEP1_PYTHIA8->AddEntry(PYTHIA8_cmult,"PYTHIA8","p");
    leg_mult_LEP1_PYTHIA8->Draw("SAME");
    xjjroot::drawtex(0.2,0.876,"LEP1 vs PYTHIA8");
    mult_LEP1_PYTHIA8->SaveAs("pdfDir/mult_LEP1_PYTHIA8.pdf");
    
    delete leg_mult_LEP1_PYTHIA8;
    delete hempty_mult1;
    delete mult_LEP1_PYTHIA8;

    /*
    TCanvas *mult_LEP2 = new TCanvas("mult","",600,600);
    TH2F *hempty_mult2 = new TH2F("Multiplicity Distribution",";Multiplicity;Probability",1,0,95,1,0,0.15);
    xjjroot::sethempty(hempty_mult2,0,0.3);
    hempty_mult2->DrawCopy();
    xjjroot::setthgrstyle(LEP2_cmult, kBlack, 21, 1.2, kBlack, 1, 1, -1, -1, -1);
    LEP2_cmult->DrawCopy("pe same");
    TLegend *leg_mult_LEP2 = new TLegend(0.67,0.7,1.1,0.88);
    xjjroot::setleg(leg_mult_LEP2);
    leg_mult_LEP2->AddEntry(LEP2_cmult,"LEP2","p");
    leg_mult_LEP2->Draw("SAME");
    xjjroot::drawtex(0.2,0.876,"LEP2");
    mult_LEP2->SaveAs("pdfDir/mult_LEP2.pdf");
    
    delete leg_mult_LEP2;
    delete hempty_mult2;
    delete mult_LEP2;
     */
    
    // pt spectra
    TCanvas *pt_LEP1_PYTHIA8 = new TCanvas("pt","",600,600);
    gPad->SetLogy();
    TH2F *hempty_pt1 = new TH2F("PT Spectra",";pt;Probability",nBinsPt,binsPt,nBinsY,binsYLog);
    //TH2F *hempty_pt1 = new TH2F("PT Spectra",";PT;Probability",1,0,95,1,0,0.15);
    xjjroot::sethempty(hempty_pt1,0,0.3);
    hempty_pt1->DrawCopy();
    xjjroot::setthgrstyle(LEP1_cpt, kRed, 21, 1.2, kRed, 1, 1, -1, -1, -1);
    LEP1_cpt->DrawCopy("pe same");
    xjjroot::setthgrstyle(PYTHIA8_cpt, kBlue, 21, 1.2, kBlue, 1, 1, -1, -1, -1);
    PYTHIA8_cpt->DrawCopy("pe same");
    TLegend *leg_pt_LEP1_PYTHIA8 = new TLegend(0.67,0.7,1.1,0.88);
    xjjroot::setleg(leg_pt_LEP1_PYTHIA8);
    leg_pt_LEP1_PYTHIA8->AddEntry(LEP1_cpt,"LEP1","p");
    leg_pt_LEP1_PYTHIA8->AddEntry(PYTHIA8_cpt,"PYTHIA8","p");
    leg_pt_LEP1_PYTHIA8->Draw("SAME");
    xjjroot::drawtex(0.35,0.876,"LEP1 vs PYTHIA8");
    pt_LEP1_PYTHIA8->SaveAs("pdfDir/pt_LEP1_PYTHIA8.pdf");
    
    delete leg_pt_LEP1_PYTHIA8;
    delete hempty_pt1;
    delete pt_LEP1_PYTHIA8;

    /*
    TCanvas *pt_LEP2 = new TCanvas("pt","",600,600);
    TH2F *hempty_pt2 = new TH2F("PT Spectra",";PT;Probability",1,0,95,1,0,0.15);
    xjjroot::sethempty(hempty_pt2,0,0.3);
    hempty_pt2->DrawCopy();
    xjjroot::setthgrstyle(LEP2_cpt, kBlack, 21, 1.2, kBlack, 1, 1, -1, -1, -1);
    LEP2_cpt->DrawCopy("pe same");
    TLegend *leg_pt_LEP2 = new TLegend(0.67,0.7,1.1,0.88);
    xjjroot::setleg(leg_pt_LEP2);
    leg_pt_LEP2->AddEntry(LEP2_cpt,"LEP2","p");
    leg_pt_LEP2->Draw("SAME");
    xjjroot::drawtex(0.2,0.876,"LEP2");
    pt_LEP2->SaveAs("pdfDir/pt_LEP2.pdf");

    delete hempty_pt2;
    delete leg_pt_LEP2;
    delete pt_LEP2;
    */
    
    // eta spectra
    TCanvas *eta_LEP1_PYTHIA8 = new TCanvas("eta","",600,600);
    TH2F *hempty_eta1 = new TH2F("Eta Spectra",";eta;Probability",nBinsEta,binsEta,nBinsEta,binsYLin);
    //TH2F *hempty_eta1 = new TH2F("Eta Spectra",";Eta;Probability",1,-5,5,1,0,0.06);
    xjjroot::sethempty(hempty_eta1,0,0.3);
    hempty_eta1->DrawCopy();
    xjjroot::setthgrstyle(LEP1_ceta, kRed, 21, 1.2, kRed, 1, 1, -1, -1, -1);
    LEP1_ceta->DrawCopy("pe same");
    xjjroot::setthgrstyle(PYTHIA8_ceta, kBlue, 21, 1.2, kBlue, 1, 1, -1, -1, -1);
    PYTHIA8_ceta->DrawCopy("pe same");
    TLegend *leg_eta_LEP1_PYTHIA8 = new TLegend(0.67,0.7,1.1,0.88);
    xjjroot::setleg(leg_eta_LEP1_PYTHIA8);
    leg_eta_LEP1_PYTHIA8->AddEntry(LEP1_ceta,"LEP1","p");
    leg_eta_LEP1_PYTHIA8->AddEntry(PYTHIA8_ceta,"PYTHIA8","p");
    leg_eta_LEP1_PYTHIA8->Draw("SAME");
    xjjroot::drawtex(0.2,0.876,"LEP1 vs PYTHIA8");
    eta_LEP1_PYTHIA8->SaveAs("pdfDir/eta_LEP1_PYTHIA8.pdf");
    
    delete leg_eta_LEP1_PYTHIA8;
    delete hempty_eta1;
    delete eta_LEP1_PYTHIA8;

    /*
    TCanvas *eta_LEP2 = new TCanvas("eta","",600,600);
    TH2F *hempty_eta2 = new TH2F("Eta Spectra",";Eta;Probability",1,-5,5,1,0,0.05);
    xjjroot::sethempty(hempty_eta2,0,0.3);
    hempty_eta2->DrawCopy();
    xjjroot::setthgrstyle(LEP2_ceta, kBlack, 21, 1.2, kBlack, 1, 1, -1, -1, -1);
    LEP2_ceta->DrawCopy("pe same");
    TLegend *leg_eta_LEP2 = new TLegend(0.67,0.7,1.1,0.88);
    xjjroot::setleg(leg_eta_LEP2);
    leg_eta_LEP2->AddEntry(LEP2_ceta,"LEP2","p");
    leg_eta_LEP2->Draw("SAME");
    xjjroot::drawtex(0.2,0.876,"LEP2");
    eta_LEP2->SaveAs("pdfDir/eta_LEP2.pdf");

    delete hempty_eta2;
    delete leg_eta_LEP2;
    delete eta_LEP2;
     */
    
    PYTHIA8->Close();
    delete PYTHIA8;
    //LEP2->Close();
    //delete LEP2;
    LEP1->Close();
    delete LEP1;
    
    return 0;
}

int drawEta
        (
            const std::string inFileName, // input file
            const std::string datalabel // Text in upper-left corner
        )
{
    TFile *file = new TFile(inFileName.c_str());
    
    xjjroot::setgstyle();
    // Plot all, neutral, charged eta
    TCanvas *alleta = new TCanvas("alleta","alleta",600,600);
    TH2F *hemptyeta = new TH2F("",";eta;Probability",1,-5,5,1,0,0.06);
    xjjroot::sethempty(hemptyeta,0,0.3);
    hemptyeta->DrawCopy();
    TH1F *aeta = (TH1F*)gDirectory->Get("aeta");
    xjjroot::setthgrstyle(aeta, kBlack, 20, 1.2, kRed, 1, 1, -1, -1, -1);
    aeta->DrawCopy("pe same");
    TH1F *ceta = (TH1F*)gDirectory->Get("ceta");
    xjjroot::setthgrstyle(ceta, kRed, 21, 1.2, kRed, 1, 1, -1, -1, -1);
    ceta->DrawCopy("pe same");
    TH1F *neta = (TH1F*)gDirectory->Get("neta");
    xjjroot::setthgrstyle(neta, kBlue, 22, 1.2, kBlue, 1, 1, -1, -1, -1);
    neta->DrawCopy("pe same");
    TLegend* etaleg = new TLegend(0.67,0.7,1.1,0.88);
    xjjroot::setleg(etaleg);
    etaleg->AddEntry(aeta,"All","p");
    etaleg->AddEntry(ceta,"Charged Hadrons","p");
    etaleg->AddEntry(neta,"Neutral Hadrons + Photons","p");
    etaleg->Draw("SAME");
    xjjroot::drawtex(0.2,0.876,datalabel.c_str());
    alleta->SaveAs(Form("pdfDir/%s_eta.pdf",datalabel.c_str()));
    
    delete etaleg;
    delete neta;
    delete ceta;
    delete aeta;
    delete alleta;
    
    file->Close();
    delete file;
    
    return 0;
}


int main(int argc, char* argv[])
{
    if(argc > 1)
    {
        std::cout << "Usage: ./" << argv[0] << std::endl;
        return 1;
    }
    
    int retVal = 0;
    retVal += plotDataQualityChecks();
    return retVal;
}
