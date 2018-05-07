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

void formatMultiPanelHist(TH1F *hist)
{
    hist->GetXaxis()->SetTitleSize(24);
    hist->GetXaxis()->SetTitleFont(43);
    hist->SetTitleOffset(2.000,"X");
    hist->GetXaxis()->SetLabelSize(20);
    hist->GetXaxis()->SetLabelFont(43);
    hist->GetXaxis()->SetNdivisions(5);
    hist->GetYaxis()->SetTitleSize(24);
    hist->GetYaxis()->SetTitleFont(43);
    hist->GetYaxis()->SetLabelSize(20);
    hist->GetYaxis()->SetLabelFont(43);
    hist->SetTitleOffset(3.000,"Y");
}

int nTrkOfflineContributions(const std::string inFileNameData,  //  Data
                             const std::string inFileNameMC,    // monte carlo
                             int data,  // 0 full directory, 1 local run on single file
                             int monteCarlo     // 1 - with Monte Carlo
                             )
{
    Selection s = Selection();
    vanGoghPalette vGP = vanGoghPalette();
    
    std::cout<<"Initializing histograms and loading data..."<<std::endl;
    /************************************************/
    // Initialize Histograms                         /
    /************************************************/
    
    const int nBins = 45;
    Double_t bins[nBins+1];
    const Double_t logLow = 0.5501;
    const Double_t logHi = 1.0499;
    getLogBins(logLow, logHi, nBins, bins);
    
    // Subranges of NtrkOffline
    static const int numRanges = 5;
    int nTrkMin[numRanges] = {4,10,20,30,35};
    int nTrkMax[numRanges] = {10,20,30,999,999};
    int colors[numRanges] = {0,1,2,3,4}; // gold,blue,red,orange,green,peach
    
        // Data //
    
    // Full Distributions //
    TH1D * nTrkOffline = new TH1D("nTrkOffline","nTrkOffline;N_{Trk}^{Offline};Entries",54,0,54);
    InitHist(nTrkOffline, "N_{Trk}^{Offline}", "Entries", kBlack, 1001,0);
    
    TH1F *thrust = new TH1F("Thrust","Thrust;Thrust;Entries",nBins,bins);
    InitHist(thrust, "Thrust", "#frac{1}{#sigma} #frac{d#sigma}{dT}", kBlack, 1001,0);
    formatMultiPanelHist(thrust);
    
    TH1F *nJets = new TH1F("nJets","nJets;N_{Jets}^{Offline};Entries",24,0.0001,24.9999);
    InitHist(nJets, "N_{Jets}^{Offline}", "Entries", kBlack, 1001,0);
    formatMultiPanelHist(nJets);
    
    // nTrk cuts //
    TH1D * nTrkOffline_Cut[numRanges];
    TH1F * thrust_Cut[numRanges];
    TH1F * nJets_Cut[numRanges];
    
    for(unsigned int i = 0; i < numRanges; ++i)
    {
        nTrkOffline_Cut[i] = new TH1D(Form("nTrkOffline_%d",nTrkMin[i]), Form("nTrkOffline_%d;N_{Trk}^{Offline};Entries",nTrkMin[i]), 54,0,54); //10, nTrkMax[i],nTrkMax[i+1]
        InitHist(nTrkOffline_Cut[i], "N_{Trk}^{Offline}", "Entries", kBlack, 1001,0);
        thrust_Cut[i] = new TH1F(Form("Thrust_%d",nTrkMin[i]), Form("Thrust_%d;Thrust;Entries",nTrkMin[i]),nBins,bins);
        InitHist(thrust_Cut[i], "", "Entries", kBlack, 1001,0);
        nJets_Cut[i] = new TH1F(Form("nJets_%d",nTrkMin[i]), Form("nJets_%d;N_{Jets}^{Offline};Entries",nTrkMin[i]),24,0.0001,24.9999);
        InitHist(nJets_Cut[i], "", "Entries", kBlack, 1001,0);
    }
    
    
        // Monte Carlo //
    
    // Full Distributions //
    TH1D * nTrkOffline_MC = new TH1D("nTrkOffline_MC","nTrkOffline_MC;N_{Trk}^{Offline};Entries",54,0,54);
    InitHist(nTrkOffline_MC, "N_{Trk}^{Offline}", "Entries", kBlack, 1001,0);
    
    TH1F *thrust_MC = new TH1F("Thrust_MC","Thrust_MC;Thrust;Entries",nBins,bins);
    InitHist(thrust_MC, "Thrust", "#frac{1}{#sigma} #frac{d#sigma}{dT}", kBlack, 1001,0);
    formatMultiPanelHist(thrust_MC);

    TH1F *nJets_MC = new TH1F("nJets_MC","nJets_MC;N_{Jets}^{Offline};Entries",24,0.0001,24.9999);
    InitHist(nJets_MC, "N_{Jets}^{Offline}", "Entries", kBlack, 1001,0);
    formatMultiPanelHist(nJets_MC);

    // nTrk cuts //
    TH1D * nTrkOffline_MC_Cut[numRanges];
    TH1F * thrust_MC_Cut[numRanges];
    TH1F * nJets_MC_Cut[numRanges];
    
    for(unsigned int i = 0; i < numRanges; ++i)
    {
        nTrkOffline_MC_Cut[i] = new TH1D(Form("nTrkOffline_MC_%d",nTrkMin[i]), Form("nTrkOffline_MC_%d;N_{Trk}^{Offline};Entries",nTrkMin[i]), 54,0,54); //10, nTrkMax[i],nTrkMax[i+1]
        InitHist(nTrkOffline_MC_Cut[i], "N_{Trk}^{Offline}", "Entries", kBlack, 1001,0);
        thrust_MC_Cut[i] = new TH1F(Form("Thrust_MC_%d",nTrkMin[i]), Form("Thrust_MC_%d;Thrust;Entries",nTrkMin[i]),nBins,bins);
        InitHist(thrust_MC_Cut[i], "", "Entries", kBlack, 1001,0);
        nJets_MC_Cut[i] = new TH1F(Form("nJets_MC_%d",nTrkMin[i]), Form("nJets_MC_%d;N_{Jets}^{Offline};Entries",nTrkMin[i]),24,0.0001,24.9999);
    }

    /************************************************/
    // Initialize Tree and Variables                 /
    /************************************************/
    
    TChain chain_t("t");
    TChain chain_jt("akR4ESchemeJetTree"); // should use the smartJetName method but not sure how to implement it here so hardcoded name for now
    
    TChain chain_MC_t("t");
    TChain chain_MC_jt("akR4ESchemeJetTree");
    
    // Data
    std::string LEP1Data1992_recons_aftercut = "/LEP1Data1992_recons_aftercut-MERGED.root";
    std::string LEP1Data1993_recons_aftercut = "/LEP1Data1993_recons_aftercut-MERGED.root";
    std::string LEP1Data1994P1_recons_aftercut = "/LEP1Data1994P1_recons_aftercut-MERGED.root";
    std::string LEP1Data1994P2_recons_aftercut = "/LEP1Data1994P2_recons_aftercut-MERGED.root";
    std::string LEP1Data1994P3_recons_aftercut = "/LEP1Data1994P3_recons_aftercut-MERGED.root";
    std::string LEP1Data1995_recons_aftercut = "/LEP1Data1995_recons_aftercut-MERGED.root";
    
    std::string LEP1Data1992_recons_aftercut_Mix = "/LEP1Data1992_recons_aftercut-MERGED_Mix.root";
    std::string LEP1Data1993_recons_aftercut_Mix = "/LEP1Data1993_recons_aftercut-MERGED_Mix.root";
    std::string LEP1Data1994P1_recons_aftercut_Mix = "/LEP1Data1994P1_recons_aftercut-MERGED_Mix.root";
    std::string LEP1Data1994P2_recons_aftercut_Mix = "/LEP1Data1994P2_recons_aftercut-MERGED_Mix.root";
    std::string LEP1Data1994P3_recons_aftercut_Mix = "/LEP1Data1994P3_recons_aftercut-MERGED_Mix.root";
    std::string LEP1Data1995_recons_aftercut_Mix = "/LEP1Data1995_recons_aftercut-MERGED_Mix.root";
    
    
    // Monte Carlo
    std::string alephMCRecoAfterCutPaths_1994 = "/alephMCRecoAfterCutPaths_1994.root";
    
    std::string alephMCRecoAfterCutPaths_1994_Mix = "/alephMCRecoAfterCutPaths_1994_Mix.root";

    if (data == 0) // full data set
    {
        chain_t.Add( (inFileNameData+LEP1Data1992_recons_aftercut).c_str() );
        chain_t.Add( (inFileNameData+LEP1Data1993_recons_aftercut).c_str() );
        chain_t.Add( (inFileNameData+LEP1Data1994P1_recons_aftercut).c_str() );
        chain_t.Add( (inFileNameData+LEP1Data1994P2_recons_aftercut).c_str() );
        chain_t.Add( (inFileNameData+LEP1Data1994P3_recons_aftercut).c_str() );
        chain_t.Add( (inFileNameData+LEP1Data1995_recons_aftercut).c_str() );
        
        chain_jt.Add( (inFileNameData+LEP1Data1992_recons_aftercut).c_str() );
        chain_jt.Add( (inFileNameData+LEP1Data1993_recons_aftercut).c_str() );
        chain_jt.Add( (inFileNameData+LEP1Data1994P1_recons_aftercut).c_str() );
        chain_jt.Add( (inFileNameData+LEP1Data1994P2_recons_aftercut).c_str() );
        chain_jt.Add( (inFileNameData+LEP1Data1994P3_recons_aftercut).c_str() );
        chain_jt.Add( (inFileNameData+LEP1Data1995_recons_aftercut).c_str() );
    }
    
    if (data == 1) // single file run
    {
        chain_t.Add(inFileNameData.c_str());
        
        chain_jt.Add(inFileNameData.c_str());
    }
    
    if (monteCarlo == 1) // monte-carlo
    {
        chain_MC_t.Add( (inFileNameMC+alephMCRecoAfterCutPaths_1994).c_str() );
        
        chain_MC_jt.Add( (inFileNameMC+alephMCRecoAfterCutPaths_1994).c_str() );
    }

    std::cout<<"Looping over data..."<<std::endl;
    /************************************************/
    // Loop over data                                /
    /************************************************/
    
    //Int_t nevent = (Int_t)chain_t->GetEntries();
    Int_t nevent = (Int_t)chain_t.GetEntries();
    
    int nChargedHadronsHP;      chain_t.SetBranchAddress("nChargedHadronsHP",&nChargedHadronsHP); //t->SetBranchAddress("nChargedHadronsHP",&nChargedHadronsHP);
    float Thrust;       chain_t.SetBranchAddress("Thrust",&Thrust); //t->SetBranchAddress("Thrust",&Thrust);
    int nref;       chain_jt.SetBranchAddress("nref",&nref); //jt->SetBranchAddress("nref",&nref);
    
    // Setup Progress bar
    ProgressBar Bar(cout, nevent);
    Bar.SetStyle(1);
    unsigned int entryDiv = (nevent > 200) ? nevent / 200 : 1;
    
    for (Int_t i=0;i<nevent;i++)
    {
        //t->GetEntry(i);
        //jt->GetEntry(i);
        
        chain_t.GetEntry(i);
        chain_jt.GetEntry(i);
        
        Bar.Update(i);
        Bar.PrintWithMod(entryDiv);
        
        // no cut
        nTrkOffline->Fill(nChargedHadronsHP);
        thrust->Fill(Thrust);
        nJets->Fill(nref);
        
        // 4-10
        if(nChargedHadronsHP >= nTrkMin[0] && nChargedHadronsHP < nTrkMax[0])
        {
            nTrkOffline_Cut[0]->Fill(nChargedHadronsHP);
            thrust_Cut[0]->Fill(Thrust);
            nJets_Cut[0]->Fill(nref);
        }
        // 10-20
        else if(nChargedHadronsHP >= nTrkMin[1] && nChargedHadronsHP < nTrkMax[1])
        {
            nTrkOffline_Cut[1]->Fill(nChargedHadronsHP);
            thrust_Cut[1]->Fill(Thrust);
            nJets_Cut[1]->Fill(nref);
        }
        // 20-30
        else if(nChargedHadronsHP >= nTrkMin[2] && nChargedHadronsHP < nTrkMax[2])
        {
            nTrkOffline_Cut[2]->Fill(nChargedHadronsHP);
            thrust_Cut[2]->Fill(Thrust);
            nJets_Cut[2]->Fill(nref);
        }
        // 30-999
        else if(nChargedHadronsHP >= nTrkMin[3] && nChargedHadronsHP < nTrkMax[3])
        {
            nTrkOffline_Cut[3]->Fill(nChargedHadronsHP);
            thrust_Cut[3]->Fill(Thrust);
            nJets_Cut[3]->Fill(nref);
            
             // 35-999
            if ( nChargedHadronsHP >= nTrkMin[4] && nChargedHadronsHP < nTrkMax[4] )
            {
                nTrkOffline_Cut[4]->Fill(nChargedHadronsHP);
                thrust_Cut[4]->Fill(Thrust);
                nJets_Cut[4]->Fill(nref);
            }
        }
    }
    
    std::cout<<"Looping over monte carlo..."<<std::endl;
    /************************************************/
    // Loop over Monte Carlo                         /
    /************************************************/
    
    //Int_t nevent = (Int_t)chain_t->GetEntries();
    Int_t nevent_MC = (Int_t)chain_MC_t.GetEntries();
    
    int nChargedHadronsHP_MC;      chain_MC_t.SetBranchAddress("nChargedHadronsHP",&nChargedHadronsHP_MC); //t->SetBranchAddress("nChargedHadronsHP",&nChargedHadronsHP);
    float Thrust_MC;       chain_MC_t.SetBranchAddress("Thrust",&Thrust_MC); //t->SetBranchAddress("Thrust",&Thrust);
    int nref_MC;       chain_MC_jt.SetBranchAddress("nref",&nref_MC); //jt->SetBranchAddress("nref",&nref);
    
    // Setup Progress bar
    ProgressBar Bar_MC(cout, nevent_MC);
    Bar_MC.SetStyle(1);
    unsigned int entryDiv_MC = (nevent_MC > 200) ? nevent_MC / 200 : 1;
    
    for (Int_t i=0;i<nevent_MC;i++)
    {
        //t->GetEntry(i);
        //jt->GetEntry(i);
        
        chain_MC_t.GetEntry(i);
        chain_MC_jt.GetEntry(i);
        
        Bar_MC.Update(i);
        Bar_MC.PrintWithMod(entryDiv);
        
        // no cut
        nTrkOffline_MC->Fill(nChargedHadronsHP_MC);
        thrust_MC->Fill(Thrust_MC);
        nJets_MC->Fill(nref_MC);
        
        // 4-10
        if(nChargedHadronsHP_MC >= nTrkMin[0] && nChargedHadronsHP_MC < nTrkMax[0])
        {
            nTrkOffline_MC_Cut[0]->Fill(nChargedHadronsHP_MC);
            thrust_MC_Cut[0]->Fill(Thrust_MC);
            nJets_MC_Cut[0]->Fill(nref_MC);
        }
        // 10-20
        else if(nChargedHadronsHP_MC >= nTrkMin[1] && nChargedHadronsHP_MC < nTrkMax[1])
        {
            nTrkOffline_MC_Cut[1]->Fill(nChargedHadronsHP_MC);
            thrust_MC_Cut[1]->Fill(Thrust_MC);
            nJets_MC_Cut[1]->Fill(nref_MC);
        }
        // 20-30
        else if(nChargedHadronsHP_MC >= nTrkMin[2] && nChargedHadronsHP_MC < nTrkMax[2])
        {
            nTrkOffline_MC_Cut[2]->Fill(nChargedHadronsHP_MC);
            thrust_MC_Cut[2]->Fill(Thrust_MC);
            nJets_MC_Cut[2]->Fill(nref_MC);
        }
        // 30-999
        else if(nChargedHadronsHP_MC >= nTrkMin[3] && nChargedHadronsHP_MC < nTrkMax[3])
        {
            nTrkOffline_MC_Cut[3]->Fill(nChargedHadronsHP_MC);
            thrust_MC_Cut[3]->Fill(Thrust_MC);
            nJets_MC_Cut[3]->Fill(nref_MC);
            
            // 35-999
            if ( nChargedHadronsHP_MC >= nTrkMin[4] && nChargedHadronsHP_MC < nTrkMax[4] )
            {
                nTrkOffline_MC_Cut[4]->Fill(nChargedHadronsHP_MC);
                thrust_MC_Cut[4]->Fill(Thrust_MC);
                nJets_MC_Cut[4]->Fill(nref_MC);
            }
        }
    }
    
    std::cout<<"Plotting..."<<std::endl;
    /************************************************/
    // Data Plotting                                 /
    /************************************************/
    
    std::string collisionLabel = "e^{+}e^{-} #rightarrow hadrons, #sqrt{s} = 91 GeV";
    std::string ALEPH = "ALEPH Archived Data";
    std::string PYTHIA = "Archived PYTHIA 6.1 MC";
    
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
    latex.DrawLatex(0.624,3000,"Inclusive N_{Trk}^{Offline}"); // 0.585
    
    for(unsigned int i = 0; i < numRanges; ++i)
    {
        multiPanel_thrust->cd(i+2);
        gPad->SetLogy();
        thrust->Draw();
        thrust_Cut[i]->SetFillColorAlpha(vGP.getColor(colors[i]),0.85);
        thrust_Cut[i]->Draw("same");
        if(i<numRanges-2) latex.DrawLatex(0.624,3000,Form("%d #leq N_{Trk}^{Offline} < %d",nTrkMin[i],nTrkMax[i]));
        else latex.DrawLatex(0.624,3000,Form("N_{Trk}^{Offline} #geq %d",nTrkMin[i]));
    }
    
    multiPanel_thrust->cd(1);
    TLegend *leg_ThrLabel = new TLegend(0.22,0.62,0.67,0.75); //0.2,0.65,0.5,0.75
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
    latex_nJet.SetTextSize(18);
    latex_nJet.SetTextAlign(13);
    latex_nJet.SetTextFont(43);
    
    multiPanel_nJets->cd(1);
    gPad->SetLogy();
    nJets->Draw();
    latex_nJet.DrawLatex(14,5500,"Inclusive N_{Trk}^{Offline}");
    
    for(unsigned int i = 0; i < numRanges; ++i)
    {
        multiPanel_nJets->cd(i+2);
        gPad->SetLogy();
        nJets->Draw();
        nJets_Cut[i]->SetFillColorAlpha(vGP.getColor(colors[i]),0.85);
        nJets_Cut[i]->Draw("same");
        if(i<numRanges-2) latex_nJet.DrawLatex(14,5500,Form("%d #leq N_{Trk}^{Offline} < %d",nTrkMin[i],nTrkMax[i]));
        else latex_nJet.DrawLatex(16,5500,Form("N_{Trk}^{Offline} #geq %d",nTrkMin[i]));
    }
    
    multiPanel_nJets->cd(1);
    TLegend *leg_nJetsLabel = new TLegend(0.23,0.12,0.68,0.25); //0.47,0.67,0.94,0.75
    leg_nJetsLabel->SetBorderSize(0);
    leg_nJetsLabel->SetFillStyle(0);
    leg_nJetsLabel->SetTextSize(0.04);
    leg_nJetsLabel->AddEntry(nJets,collisionLabel.c_str(),"t");
    leg_nJetsLabel->AddEntry(nJets,ALEPH.c_str(),"t");
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
    
    /************************************************/
    // Monte Carlo Plotting                          /
    /************************************************/
    
    // data has open markers, monte carlo has full markers
    if (monteCarlo == 1)
    {
        // Normalize all histograms to allow comparison between Monte Carlo and Data (clone data histograms) //
        float nTrkOffline_Integral = nTrkOffline->Integral();
        float thrust_Integral = thrust->Integral();
        float nJets_Integral = nJets->Integral();
        
        float nTrkOffline_MC_Integral = nTrkOffline_MC->Integral();
        float thrust_MC_Integral = thrust_MC->Integral();
        float nJets_MC_Integral = nJets_MC->Integral();
        
        TH1D *nTrkOffline_Clone = (TH1D*)nTrkOffline->Clone("nTrkOffline_Clone");
        nTrkOffline_Clone->SetMarkerStyle(24); //open circle
        nTrkOffline_Clone->SetMarkerSize(2);
        TH1F *thrust_Clone = (TH1F*)thrust->Clone("thrust_Clone");
        thrust_Clone->SetMarkerStyle(24); //open circle
        thrust_Clone->SetMarkerSize(1);
        TH1F *nJets_Clone = (TH1F*)nJets->Clone("nJets_Clone");
        nJets_Clone->SetMarkerStyle(24); //open circle
        nJets_Clone->SetMarkerSize(1);
        
        TH1D * nTrkOffline_Clone_Cut[numRanges];
        TH1F * thrust_Clone_Cut[numRanges];
        TH1F * nJets_Clone_Cut[numRanges];
        
        nTrkOffline_Clone->Scale(1./nTrkOffline_Integral);
        thrust_Clone->Scale(1./thrust_Integral);
        nJets_Clone->Scale(1./nJets_Integral);
        
        nTrkOffline_MC->Scale(1./nTrkOffline_MC_Integral);
        nTrkOffline_MC->SetMarkerStyle(33); //filled diamond
        nTrkOffline_MC->SetMarkerSize(2);
        thrust_MC->Scale(1./thrust_MC_Integral);
        thrust_MC->SetMarkerStyle(33); //filled diamond
        thrust_MC->SetMarkerSize(1);
        nJets_MC->Scale(1./nJets_MC_Integral);
        nJets_MC->SetMarkerStyle(33); //filled diamond
        nJets_MC->SetMarkerSize(1);
        
        for(unsigned int i = 0; i < numRanges; ++i)
        {
            // normalize these to the full distributions
            nTrkOffline_Clone_Cut[i] = (TH1D*)nTrkOffline_Cut[i]->Clone(Form("nTrkOffline_Clone_%d",nTrkMin[i]));
            nTrkOffline_Clone_Cut[i]->Scale(1./nTrkOffline_Integral);
            nTrkOffline_Clone_Cut[i]->SetMarkerStyle(25); //open box
            nTrkOffline_Clone_Cut[i]->SetMarkerSize(1);
            nTrkOffline_Clone_Cut[i]->SetMarkerColorAlpha(vGP.getColor(colors[i]),0.85);
            
            thrust_Clone_Cut[i] = (TH1F*)thrust_Cut[i]->Clone(Form("thrust_Clone_%d",nTrkMin[i]));
            thrust_Clone_Cut[i]->Scale(1./thrust_Integral);
            thrust_Clone_Cut[i]->SetMarkerStyle(25); //open box
            thrust_Clone_Cut[i]->SetMarkerSize(1);
            thrust_Clone_Cut[i]->SetMarkerColorAlpha(vGP.getColor(colors[i]),0.85);
            
            nJets_Clone_Cut[i] = (TH1F*)nJets_Cut[i]->Clone(Form("nJets_Clone_%d",nTrkMin[i]));
            nJets_Clone_Cut[i]->Scale(1./nJets_Integral);
            nJets_Clone_Cut[i]->SetMarkerStyle(25); //open box
            nJets_Clone_Cut[i]->SetMarkerSize(1);
            nJets_Clone_Cut[i]->SetMarkerColorAlpha(vGP.getColor(colors[i]),0.85);
        
            nTrkOffline_MC_Cut[i]->Scale(1./nTrkOffline_MC_Integral);
            nTrkOffline_MC_Cut[i]->SetMarkerStyle(34); //filled cross
            nTrkOffline_MC_Cut[i]->SetMarkerSize(1);
            nTrkOffline_MC_Cut[i]->SetMarkerColorAlpha(vGP.getColor(colors[i]),0.85);
            
            thrust_MC_Cut[i]->Scale(1./thrust_MC_Integral);
            thrust_MC_Cut[i]->SetMarkerStyle(34); //filled cross
            thrust_MC_Cut[i]->SetMarkerSize(1);
            thrust_MC_Cut[i]->SetMarkerColorAlpha(vGP.getColor(colors[i]),0.85);
            
            nJets_MC_Cut[i]->Scale(1./nJets_MC_Integral);
            nJets_MC_Cut[i]->SetMarkerStyle(34); //filled cross
            nJets_MC_Cut[i]->SetMarkerSize(1);
            nJets_MC_Cut[i]->SetMarkerColorAlpha(vGP.getColor(colors[i]),0.85);
        }
        
        /************************************************/
        // nTrkOffline                                   /
        /************************************************/
        
        TCanvas *c_nTrkOffline_MC_Cut = new TCanvas("c_nTrkOffline_MC_Cut", "c_nTrkOffline_MC_Cut", 600, 600);
        c_nTrkOffline_MC_Cut->SetLeftMargin     (0.18);
        c_nTrkOffline_MC_Cut->SetRightMargin    (0.05);
        c_nTrkOffline_MC_Cut->SetBottomMargin   (0.15);
        gPad->SetLogy();
        c_nTrkOffline_MC_Cut->SetTickx(1);
        c_nTrkOffline_MC_Cut->SetTicky(1);
        nTrkOffline_Clone->Draw();
        nTrkOffline_MC->Draw("same");
        
        TLegend *leg_nTrk_MC = new TLegend(0.22,0.22,0.63,0.32);
        leg_nTrk_MC->SetBorderSize(0);
        leg_nTrk_MC->SetFillStyle(0);
        leg_nTrk_MC->SetTextSize(0.035);
        leg_nTrk_MC->AddEntry(nTrkOffline_MC,collisionLabel.c_str(),"t");
        leg_nTrk_MC->AddEntry(nTrkOffline_Clone,ALEPH.c_str(),"p");
        leg_nTrk_MC->AddEntry(nTrkOffline_MC,PYTHIA.c_str(),"p");
        leg_nTrk_MC->Draw();
        
        /*
        for(unsigned int i = 0; i < numRanges; ++i)
        {
            nTrkOffline_Cut[i]->SetFillColorAlpha(vGP.getColor(colors[i]),0.85);
            nTrkOffline_Cut[i]->Draw("same");
        }
         */
        plotLogo(1,1,1);
        
        /************************************************/
        // Thrust                                        /
        /************************************************/
        
        TCanvas *multiPanel_thrust_MC = new TCanvas("multiPanel_thrust_MC","multiPanel_thrust_MC",1000,1000);
        TPad* pads_thrust_MC[nPads];
        formatMultipanelViewer(multiPanel_thrust_MC, pads_thrust_MC, rows, columns);
        
        TLatex latex_thrust_MC;
        latex_thrust_MC.SetTextSize(18);
        latex_thrust_MC.SetTextAlign(13);
        latex_thrust_MC.SetTextFont(43);
        
        multiPanel_thrust_MC->cd(1);
        gPad->SetLogy();
        thrust_Clone->Draw();
        thrust_MC->Draw("same");
        latex_thrust_MC.DrawLatex(0.624,0.05,"Inclusive N_{Trk}^{Offline}");
        
        for(unsigned int i = 0; i < numRanges; ++i)
        {
            multiPanel_thrust_MC->cd(i+2);
            gPad->SetLogy();
            thrust_Clone->Draw();
            thrust_MC->Draw("same");
            thrust_Clone_Cut[i]->Draw("same");
            thrust_MC_Cut[i]->Draw("same");
            
            if(i<numRanges-2) latex_thrust_MC.DrawLatex(0.624,0.05,Form("%d #leq N_{Trk}^{Offline} < %d",nTrkMin[i],nTrkMax[i]));
            else latex_thrust_MC.DrawLatex(0.624,0.05,Form("N_{Trk}^{Offline} #geq %d",nTrkMin[i]));
        }
        
        multiPanel_thrust_MC->cd(1);
        TLegend *leg_ThrLabel_MC = new TLegend(0.22,0.62,0.67,0.75);
        leg_ThrLabel_MC->SetBorderSize(0);
        leg_ThrLabel_MC->SetFillStyle(0);
        leg_ThrLabel_MC->SetTextSize(0.04);
        leg_ThrLabel_MC->AddEntry(thrust_MC,collisionLabel.c_str(),"t");
        leg_ThrLabel_MC->AddEntry(thrust_Clone,ALEPH.c_str(),"p");
        leg_ThrLabel_MC->AddEntry(thrust_MC,PYTHIA.c_str(),"p");
        leg_ThrLabel_MC->Draw();
        multiPanel_thrust_MC->cd();
        plotLogo(1,1,1);
        
        /************************************************/
        // nJet                                          /
        /************************************************/
        
        TCanvas *multiPanel_nJets_MC = new TCanvas("multiPanel_nJets_MC","multiPanel_nJets_MC",1000,1000);
        TPad* pads_nJets_MC[nPads];
        formatMultipanelViewer(multiPanel_nJets_MC, pads_nJets_MC, rows, columns);
        
        TLatex latex_nJet_MC;
        latex_nJet_MC.SetTextSize(18);
        latex_nJet_MC.SetTextAlign(13);
        latex_nJet_MC.SetTextFont(43);
        
        multiPanel_nJets_MC->cd(1);
        gPad->SetLogy();
        nJets_Clone->Draw();
        nJets_MC->Draw("same");
        latex_nJet_MC.DrawLatex(14,0.08,"Inclusive N_{Trk}^{Offline}");
        
        for(unsigned int i = 0; i < numRanges; ++i)
        {
            multiPanel_nJets_MC->cd(i+2);
            gPad->SetLogy();
            nJets_Clone->Draw();
            nJets_MC->Draw("same");
            nJets_Clone_Cut[i]->Draw("same");
            nJets_MC_Cut[i]->Draw("same");
            if(i<numRanges-2) latex_nJet_MC.DrawLatex(14,0.08,Form("%d #leq N_{Trk}^{Offline} < %d",nTrkMin[i],nTrkMax[i]));
            else latex_nJet_MC.DrawLatex(16,0.08,Form("N_{Trk}^{Offline} #geq %d",nTrkMin[i]));
        }
        
        multiPanel_nJets_MC->cd(1);
        TLegend *leg_nJetsLabel_MC = new TLegend(0.23,0.12,0.68,0.25);
        leg_nJetsLabel_MC->SetBorderSize(0);
        leg_nJetsLabel_MC->SetFillStyle(0);
        leg_nJetsLabel_MC->SetTextSize(0.04);
        leg_nJetsLabel_MC->AddEntry(nJets_MC,collisionLabel.c_str(),"t");
        leg_nJetsLabel_MC->AddEntry(nJets_Clone,ALEPH.c_str(),"p");
        leg_nJetsLabel_MC->AddEntry(nJets_MC,PYTHIA.c_str(),"p");
        leg_nJetsLabel_MC->Draw();
        multiPanel_nJets_MC->cd();
        plotLogo(1,1,1);
        
        /************************************************/
        // Save                                          /
        /************************************************/
        
        c_nTrkOffline_MC_Cut->RedrawAxis();
        c_nTrkOffline_MC_Cut->SaveAs("plots/nTrkOffline_MC.eps");
        for(unsigned int i = 0; i < numRanges+1; ++i)
        {
            multiPanel_thrust_MC->cd(i+1)->RedrawAxis();
            multiPanel_nJets_MC->cd(i+1)->RedrawAxis();
        }
        multiPanel_thrust_MC->SaveAs("plots/multiPanel_thrust_nTrk_MC.eps");
        multiPanel_nJets_MC->SaveAs("plots/multiPanel_nJets_nTrk_MC.eps");
        
    }
    
    return 0;
}
