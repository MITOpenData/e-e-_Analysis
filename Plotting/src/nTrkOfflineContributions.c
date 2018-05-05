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

int nTrkOfflineContributions(const std::string inFileName, 		// Input file
                             std::string outFileName    	// Output file
                            )
{
    
    Selection s = Selection();
    
    /********************************************************************************************************************/
    // Initialize histograms
    /********************************************************************************************************************/
    
    //declare some binnings first for the thrust distribution
    const int nBins = 45;
    Double_t bins[nBins+1];
    const Double_t logLow = .55;
    const Double_t logHi = 1.05;
    getLogBins(logLow, logHi, nBins, bins);
    
    // full distributions
    TH1D * nTrkOffline = new TH1D("nTrkOffline","nTrkOffline;N_{Trk}^{Offline};Entries",54,0,54);
    formatTH1(nTrkOffline,.8,1.5);
    TH1F *Thrust = new TH1F("Thrust","Thrust;Thrust;Entries",nBins,bins);
    formatTH1(Thrust,.8,1.5);
    TH1F *nJets = new TH1F("nJets","nJets;N_{Jets}^{Offline};Entries",25,0,25);
    formatTH1(nJets,.8,1.5);
    
    // Subranges of NtrkOffline
    static const int numRanges = 5;
    int nTrkMax[numRanges+1] = {0,10,20,30,40,999};
    TH1D * nTrkOffline_Cut[numRanges];
    TH1F * Thrust_Cut[numRanges];
    TH1F * nJets_Cut[numRanges];
    
    for(unsigned int i = 0; i < numRanges; ++i)
    {
        nTrkOffline_Cut[i] = new TH1D(Form("nTrkOffline_%d",nTrkMax[i+1]), Form("nTrkOffline_%d;N_{Trk}^{Offline};Entries",nTrkMax[i+1]), 10, nTrkMax[i],nTrkMax[i+1]);
        Thrust_Cut[i] = new TH1F(Form("Thrust_%d",nTrkMax[i+1]), Form("Thrust_%d;Thrust;Entries",nTrkMax[i+1]),nBins,bins);
        nJets_Cut[i] = new TH1F(Form("nJets_%d",nTrkMax[i+1]), Form("nJets_%d;N_{Jets}^{Offline};Entries",nTrkMax[i+1]),25,0,25);
    }

    /********************************************************************************************************************/
    // Define the output file
    /********************************************************************************************************************/
    
    /// Initialize the trees for use
    std::cout<<"Initializing trees for use..."<<std::endl;
    std::string jtTreeName = "";
    if (s.jttree == 0) jtTreeName = "akR4ESchemeJetTree";
    if (s.jttree == 1) jtTreeName = "akR4WTAmodpSchemeJetTree";
    if (s.jttree == 2) jtTreeName = "akR8ESchemeJetTree";
    if (s.jttree == 3) jtTreeName = "akR8WTAmodpSchemeJetTree";
    if (s.jttree == 4) jtTreeName = "ktN2WTAmodpSchemeJetTree";
    
    if(inFileName.find(".root") != std::string::npos){
        TFile* temp_p = new TFile(inFileName.c_str(), "READ");
        jtTreeName = smartJetName(jtTreeName, temp_p);
        temp_p->Close();
        delete temp_p;
    }
    
    // files and variables for input
    // get the correct tree for the analysis
    std::string boostTree = "BoostedWTAR8Evt";
    
    TChain * t = new TChain("t");       			t->Add(inFileName.c_str());
    TChain * boost_t = new TChain(boostTree.c_str());       	boost_t->Add(inFileName.c_str());
    TChain * jt = new TChain(jtTreeName.c_str());       	jt->Add(inFileName.c_str());
    
    TPCNtupleData data(s.doBelle, s.doThrust, s.doWTA);      	data.setupTPCTree(t,boost_t,jt);
    
    Int_t nevent = (Int_t)t->GetEntries();
    
    // Setup Progress bar
    ProgressBar Bar(cout, nevent);
    Bar.SetStyle(1);
    unsigned int entryDiv = (nevent > 200) ? nevent / 200 : 1;

    for (Int_t i=0;i<nevent;i++)
    {
        t->GetEntry(i);
        boost_t->GetEntry(i);
        jt->GetEntry(i);
        Bar.Update(i);
        Bar.PrintWithMod(entryDiv);
        
        // nTrk calculation
        Int_t nTrk = s.ridge_eventSelection(&data.event, &data.jet, &data.particle);
        
        nTrkOffline->Fill(nTrk);
        Thrust->Fill(data.event.Thrust);
        nJets->Fill(data.jet.nref);
        
        if(nTrk >= nTrkMax[0] && nTrk < nTrkMax[1]) {       nTrkOffline_Cut[0]->Fill(nTrk);         Thrust_Cut[0]->Fill(data.event.Thrust);         nJets_Cut[0]->Fill(data.jet.nref);} // 0-10
        else if(nTrk >= nTrkMax[1] && nTrk < nTrkMax[2]) {       nTrkOffline_Cut[1]->Fill(nTrk);         Thrust_Cut[1]->Fill(data.event.Thrust);         nJets_Cut[1]->Fill(data.jet.nref);} // 10-20
        else if(nTrk >= nTrkMax[2] && nTrk < nTrkMax[3]) {       nTrkOffline_Cut[2]->Fill(nTrk);         Thrust_Cut[2]->Fill(data.event.Thrust);         nJets_Cut[2]->Fill(data.jet.nref);} // 20-30
        else if(nTrk >= nTrkMax[3] && nTrk < nTrkMax[4]) {       nTrkOffline_Cut[3]->Fill(nTrk);         Thrust_Cut[3]->Fill(data.event.Thrust);         nJets_Cut[3]->Fill(data.jet.nref);} // 30-40
        else {       nTrkOffline_Cut[4]->Fill(nTrk);         Thrust_Cut[4]->Fill(data.event.Thrust);         nJets_Cut[4]->Fill(data.jet.nref);} //40-Inf
        
    }
    
    TCanvas *c_nTrkOffline_Cut = new TCanvas("c_nTrkOffline_Cut","c_nTrkOffline_Cut",1000,1000);
    nTrkOffline->Draw();
    for(unsigned int i = 0; i < numRanges; ++i)
    {
        nTrkOffline_Cut[i]->SetFillStyle(3000+i);
        nTrkOffline_Cut[i]->SetFillColor(i+1);
        nTrkOffline_Cut[i]->Draw("same");
    }
    
    /*
    TCanvas *c_Thrust_Cut = new TCanvas("c_Thrust_Cut","c_Thrust_Cut",1000,1000);
    c_Thrust_Cut->SetLogy();
    Thrust->Draw();
    for(unsigned int i = 0; i < numRanges; ++i)
    {
        Thrust_Cut[i]->SetFillStyle(3000+i);
        Thrust_Cut[i]->SetFillColor(i+1);
        Thrust_Cut[i]->Draw("same");
    }
    */
    TCanvas *c_Thrust_Cut_summary = new TCanvas("c_Thrust_Cut_summary","c_Thrust_Cut_summary",1000,1000);
    c_Thrust_Cut_summary->Divide(2,3);
      //align at top
    for(unsigned int i = 0; i < numRanges; ++i)
    {
        c_Thrust_Cut_summary->cd(i+1);
        gPad->SetLogy();
        Thrust->Draw();
        Thrust_Cut[i]->SetFillStyle(3000+i);
        Thrust_Cut[i]->SetFillColor(i+1);
        Thrust_Cut[i]->Draw("same");
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.SetTextAlign(13);
        latex.DrawLatex(.2,.8,Form("%d <= Ntrk < %d",nTrkMax[i],nTrkMax[i+1]));
    }
    
    TCanvas *c_nJets_summary = new TCanvas("c_nJets_summary","c_nJets_summary",1000,1000);
    c_nJets_summary->Divide(2,3);
    //align at top
    for(unsigned int i = 0; i < numRanges; ++i)
    {
        c_nJets_summary->cd(i+1);
        nJets->Draw();
        nJets_Cut[i]->SetFillStyle(3000+i);
        nJets_Cut[i]->SetFillColor(i+1);
        nJets_Cut[i]->Draw("same");
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.SetTextAlign(13);
        latex.DrawLatex(.7,.8,Form("%d <= Ntrk < %d",nTrkMax[i],nTrkMax[i+1]));
    }
    
    
    // Save the canvas
    c_nTrkOffline_Cut->SaveAs("plots/nTrkOffline.pdf");
    c_Thrust_Cut_summary->SaveAs("plots/Thrust_nTrkOffline_summary.pdf");
    c_nJets_summary->SaveAs("plots/nJets_nTrkOffline_Cut.pdf");
    
    return 0;
}
