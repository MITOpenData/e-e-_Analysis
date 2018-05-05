//
//  nTrkOfflineContributions.c
//  
//
//  Created by Anthony Badea on 4/30/18.
//
//

//c and c++ dependencies
#include <sstream>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <fstream>

//root dependencies
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TFormula.h>
#include <TNtuple.h>
#include "TDatime.h"
#include <TLegend.h>

//local headers
#include "DataProcessing/include/smartJetName.h"
#include "DataProcessing/include/trackSelection.h"
#include "TwoParticleCorrelation/include/fourier.h"
#include "TwoParticleCorrelation/include/Selection.h"
#include "TwoParticleCorrelation/include/ProgressBar.h"
#include "TwoParticleCorrelation/include/formatHists.h"
#include "TwoParticleCorrelation/include/TPCNtupleData.h"
#include "ThrustDistribution/include/getLogBins.h"
#include "Utilities/include/plotLogo.h"
#include "include/HttStyles.cc"
#include "include/vanGoghPalette.h"
#include "include/formatMultipanelViewer.h"

int nTrkOfflineContributions(const std::string inFileName)
{
    Selection s = Selection();
    vanGoghPalette vGP = vanGoghPalette();
    
    //declare some binnings first for the thrust distribution
    const int nBins = 45;
    Double_t bins[nBins+1];
    const Double_t logLow = 0.5501;
    const Double_t logHi = 1.0499;
    getLogBins(logLow, logHi, nBins, bins);
    
    // full distributions
    TH1D * nTrkOffline = new TH1D("nTrkOffline","nTrkOffline;N_{Trk}^{Offline};Entries",54,0,54);
    InitHist(nTrkOffline, "N_{Trk}^{Offline}", "Entries", kBlack, 1001,0);
    TH1F *thrust = new TH1F("Thrust","Thrust;Thrust;Entries",nBins,bins);
    InitHist(thrust, "Thrust", "#frac{1}{#sigma} #frac{d#sigma}{dT}", kBlack, 1001,0);
    thrust->GetXaxis()->SetTitleSize(18);
    thrust->GetXaxis()->SetTitleFont(43);
    thrust->SetTitleOffset(2.000,"X");
    thrust->GetXaxis()->SetLabelSize(14);
    thrust->GetXaxis()->SetLabelFont(43);
    thrust->GetYaxis()->SetTitleSize(18);
    thrust->GetYaxis()->SetTitleFont(43);
    thrust->SetTitleOffset(3.000,"Y");
    
    TH1F *nJets = new TH1F("nJets","nJets;N_{Jets}^{Offline};Entries",25,0,25);
    InitHist(nJets, "N_{Jets}^{Offline}", "Entries", kBlack, 1001,0);
    nJets->GetXaxis()->SetTitleSize(18);
    nJets->GetXaxis()->SetTitleFont(43);
    nJets->SetTitleOffset(2.000,"X");
    nJets->GetXaxis()->SetLabelSize(14);
    nJets->GetXaxis()->SetLabelFont(43);
    nJets->GetYaxis()->SetTitleSize(18);
    nJets->GetYaxis()->SetTitleFont(43);
    nJets->SetTitleOffset(3.000,"Y");
    
    // Subranges of NtrkOffline
    static const int numRanges = 5;
    int nTrkMax[numRanges+1] = {0,8,16,24,32,999};
    //static const int numRanges = 6;
    //int nTrkMax[numRanges+1] = {0,8,16,24,32,40,999};
    int colors[numRanges] = {0,1,2,3,4}; // gold,blue,red,orange,green,peach
    TH1D * nTrkOffline_Cut[numRanges];
    TH1F * thrust_Cut[numRanges];
    TH1F * nJets_Cut[numRanges];
    
    for(unsigned int i = 0; i < numRanges; ++i)
    {
        nTrkOffline_Cut[i] = new TH1D(Form("nTrkOffline_%d",nTrkMax[i+1]), Form("nTrkOffline_%d;N_{Trk}^{Offline};Entries",nTrkMax[i+1]), 54,0,54); //10, nTrkMax[i],nTrkMax[i+1]
        InitHist(nTrkOffline_Cut[i], "N_{Trk}^{Offline}", "Entries", kBlack, 1001,0);
        thrust_Cut[i] = new TH1F(Form("Thrust_%d",nTrkMax[i+1]), Form("Thrust_%d;Thrust;Entries",nTrkMax[i+1]),nBins,bins);
        InitHist(thrust_Cut[i], "", "Entries", kBlack, 1001,0);
        nJets_Cut[i] = new TH1F(Form("nJets_%d",nTrkMax[i+1]), Form("nJets_%d;N_{Jets}^{Offline};Entries",nTrkMax[i+1]),25,0,25);
    }
    
    /// Initialize the tree and variables for use
    TFile * f = TFile::Open(inFileName.c_str(),"read");
    TTree * t = (TTree*)f->Get("t");
    TTree * jt = (TTree*)f->Get(smartJetName("ak4ESchemeJetTree", f).c_str());
    
    Int_t nevent = (Int_t)t->GetEntries();
    int nChargedHadronsHP;      t->SetBranchAddress("nChargedHadronsHP",&nChargedHadronsHP);
    //float nChargedHadronsHP_Corrected;      t->SetBranchAddress("nChargedHadronsHP_Corrected",&nChargedHadronsHP_Corrected);
    float Thrust;       t->SetBranchAddress("Thrust",&Thrust);
    int nref;       jt->SetBranchAddress("nref",&nref);
    
    
    for (Int_t i=0;i<nevent;i++)
    {
        t->GetEntry(i);
        jt->GetEntry(i);
        
        // no cut
        nTrkOffline->Fill(nChargedHadronsHP);
        thrust->Fill(Thrust);
        nJets->Fill(nref);
        
        // 0-10
        if(nChargedHadronsHP >= nTrkMax[0] && nChargedHadronsHP < nTrkMax[1])
        {
            nTrkOffline_Cut[0]->Fill(nChargedHadronsHP);
            thrust_Cut[0]->Fill(Thrust);
            nJets_Cut[0]->Fill(nref);
        }
        // 10-20
        else if(nChargedHadronsHP >= nTrkMax[1] && nChargedHadronsHP < nTrkMax[2])
        {
            nTrkOffline_Cut[1]->Fill(nChargedHadronsHP);
            thrust_Cut[1]->Fill(Thrust);
            nJets_Cut[1]->Fill(nref);
        }
        // 20-30
        else if(nChargedHadronsHP >= nTrkMax[2] && nChargedHadronsHP < nTrkMax[3])
        {
            nTrkOffline_Cut[2]->Fill(nChargedHadronsHP);
            thrust_Cut[2]->Fill(Thrust);
            nJets_Cut[2]->Fill(nref);
        }
        // 30-40
        else if(nChargedHadronsHP >= nTrkMax[3] && nChargedHadronsHP < nTrkMax[4])
        {
            nTrkOffline_Cut[3]->Fill(nChargedHadronsHP);
            thrust_Cut[3]->Fill(Thrust);
            nJets_Cut[3]->Fill(nref);
        }
        //40-Inf
        else
        {
            nTrkOffline_Cut[4]->Fill(nChargedHadronsHP);
            thrust_Cut[4]->Fill(Thrust);
            nJets_Cut[4]->Fill(nref);
        }
    }
    
    std::string collisionLabel = "e^{+}e^{-} #rightarrow hadrons, #sqrt{s} = 91 GeV";
    std::string ALEPH = "ALEPH Archived Data";
    /************************************************/
    // nTrkOffline                                   /
    /************************************************/
    
    TCanvas *c_nTrkOffline_Cut = new TCanvas("c_nTrkOffline_Cut", "c_nTrkOffline_Cut", 600, 600);
    c_nTrkOffline_Cut->SetLeftMargin     (0.18);
    c_nTrkOffline_Cut->SetRightMargin    (0.05);
    //c_nTrkOffline_Cut->SetTopMargin      (0.08);
    c_nTrkOffline_Cut->SetBottomMargin   (0.15);
    gPad->SetLogy();
    c_nTrkOffline_Cut->SetTickx(1);
    c_nTrkOffline_Cut->SetTicky(1);
    nTrkOffline->Draw();

    TLegend *leg_nTrk = new TLegend(0.46,0.72,0.96,0.8);
    leg_nTrk->SetBorderSize(0);
    leg_nTrk->SetFillStyle(0);
    leg_nTrk->AddEntry(nTrkOffline,collisionLabel.c_str(),"t");
    leg_nTrk->AddEntry(nTrkOffline,ALEPH.c_str(),"t");
    leg_nTrk->Draw();
    
    for(unsigned int i = 0; i < numRanges; ++i)
    {
        //nTrkOffline_Cut[i]->SetFillStyle(3000+i);
        nTrkOffline_Cut[i]->SetFillColorAlpha(vGP.getColor(colors[i]),0.85);
        nTrkOffline_Cut[i]->Draw("same");
    }
    plotLogo(1,1,1);
    

    /************************************************/
    // Thrust                                        /
    /************************************************/
    
    static const int rows = 2;
    static const int columns = 3;
    static const int nPads = rows * columns;
    
    TCanvas *multiPanel_thrust = new TCanvas("multiPanel_thrust","multiPanel_thrust",1000,1000);
    TPad* pads_thrust[nPads];
    formatMultipanelViewer(multiPanel_thrust, pads_thrust, rows, columns);

    TLatex latex;
    latex.SetTextSize(18);
    latex.SetTextAlign(13);
    latex.SetTextFont(43);
    
    multiPanel_thrust->cd(1);
    gPad->SetLogy();
    thrust->Draw();
    latex.DrawLatex(0.59,1000,"Inclusive N_{Trk}^{Offline}");
    
    for(unsigned int i = 0; i < numRanges; ++i)
    {
        multiPanel_thrust->cd(i+2);
        gPad->SetLogy();
        thrust->Draw();
        thrust_Cut[i]->SetFillColorAlpha(vGP.getColor(colors[i]),0.85);
        thrust_Cut[i]->Draw("same");
        if(i<numRanges-1) latex.DrawLatex(0.59,1000,Form("%d #leq N_{Trk}^{Offline} < %d",nTrkMax[i],nTrkMax[i+1]));
        else latex.DrawLatex(0.59,1000,Form("N_{Trk}^{Offline} #geq %d",nTrkMax[i]));
    }
    
    multiPanel_thrust->cd(1);
    TLegend *leg_ThrLabel = new TLegend(0.2,0.65,0.5,0.75);
    leg_ThrLabel->SetBorderSize(0);
    leg_ThrLabel->SetFillStyle(0);
    leg_ThrLabel->SetTextSize(0.04);
    //leg_ThrLabel->SetTextSize(42);
    leg_ThrLabel->AddEntry(nTrkOffline,collisionLabel.c_str(),"t");
    leg_ThrLabel->AddEntry(nTrkOffline,ALEPH.c_str(),"t");
    //leg_ThrLabel->AddEntry(thrust_Cut[0],Form("%d #leq N_{Trk}^{Offline} < %d",nTrkMax[0],nTrkMax[0+1]),"t");
    leg_ThrLabel->Draw();
    multiPanel_thrust->cd();
    plotLogo(1,1,1);
    
    /************************************************/
    // nJet                                          /
    /************************************************/
    
    TCanvas *multiPanel_nJets = new TCanvas("multiPanel_nJets","multiPanel_nJets",1000,1000);
    TPad* pads_nJets[nPads];
    formatMultipanelViewer(multiPanel_nJets, pads_nJets, rows, columns);
    
    TLatex latex_nJet;
    latex_nJet.SetTextSize(13);
    latex_nJet.SetTextAlign(13);
    latex_nJet.SetTextFont(43);
    
    multiPanel_nJets->cd(1);
    gPad->SetLogy();
    nJets->Draw();
    latex_nJet.DrawLatex(16,950,"Inclusive N_{Trk}^{Offline}");
    
    for(unsigned int i = 0; i < numRanges; ++i)
    {
        multiPanel_nJets->cd(i+2);
        gPad->SetLogy();
        nJets->Draw();
        nJets_Cut[i]->SetFillColorAlpha(vGP.getColor(colors[i]),0.85);
        nJets_Cut[i]->Draw("same");
        if(i<numRanges-1) latex_nJet.DrawLatex(16,950,Form("%d #leq N_{Trk}^{Offline} < %d",nTrkMax[i],nTrkMax[i+1]));
        else latex_nJet.DrawLatex(16,950,Form("N_{Trk}^{Offline} #geq %d",nTrkMax[i]));
    }
    
    multiPanel_nJets->cd(1);
    TLegend *leg_nJetsLabel = new TLegend(0.47,0.67,0.94,0.75);
    leg_nJetsLabel->SetBorderSize(0);
    leg_nJetsLabel->SetFillStyle(0);
    leg_nJetsLabel->SetTextSize(0.031);
    //leg_ThrLabel->SetTextSize(42);
    leg_nJetsLabel->AddEntry(nJets,collisionLabel.c_str(),"t");
    leg_nJetsLabel->AddEntry(nJets,ALEPH.c_str(),"t");
    //leg_ThrLabel->AddEntry(thrust_Cut[0],Form("%d #leq N_{Trk}^{Offline} < %d",nTrkMax[0],nTrkMax[0+1]),"t");
    leg_nJetsLabel->Draw();
    multiPanel_nJets->cd();
    plotLogo(1,1,1);
    
    
    /************************************************/
    // Save                                          /
    /************************************************/
    
    c_nTrkOffline_Cut->RedrawAxis();
    c_nTrkOffline_Cut->SaveAs("plots/nTrkOffline.eps");
    for(unsigned int i = 0; i < numRanges+1; ++i)
    {
        multiPanel_thrust->cd(i+1)->RedrawAxis();
        multiPanel_nJets->cd(i+1)->RedrawAxis();
    }
    multiPanel_thrust->SaveAs("plots/multiPanel_thrust_nTrk.eps");
    multiPanel_nJets->SaveAs("plots/multiPanel_nJets_nTrk.eps");
    
    return 0;
}
