//
//  paperPlot.c
//  
//
//  Created by Anthony Badea on 4/27/18.
//
//

// C headers
#include <stdio.h>

// root headers
#include <TPad.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLatex.h>

// local headers
#include "include/styleUtil.h"
#include "TwoParticleCorrelation/include/Selection.h"
#include "TwoParticleCorrelation/include/utilities.h"
#include "TwoParticleCorrelation/include/xjjrootuti.h"
#include "Utilities/include/plotLogo.h"
#include "include/formatMultipanelViewer.h"


/**************************************************************************************/
// Handsome ROOT Style for Correlation Functions
/**************************************************************************************/

// Handsome Canvas for 2D histogram
void formatCFViewer(TCanvas *c)
{
    // ****** The MAGIC angle ****** //
    c->SetTheta(60.839);
    c->SetPhi(38.0172);
}

// Handsome Correlation Function Histogram
void formatCFTH2F(TH2F *h,float offsetx, float offsety, float offsetz)
{
    h->Sumw2();
    h->SetTitleOffset(offsetx,"X");
    h->SetTitleOffset(offsety,"Y");
    h->SetTitleOffset(offsetz,"Z");
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
    h->GetZaxis()->CenterTitle();
    h->GetXaxis()->SetTitle("#Delta#eta");
    h->GetYaxis()->SetTitle("#Delta#phi");
    h->GetZaxis()->SetTitle("#frac{1}{N^{trig}} #frac{d^{2}N^{pair}}{d#Delta#etad#Delta#phi}");
    h->SetNdivisions(505,"X");
    h->SetNdivisions(505,"Y");
    h->SetNdivisions(505,"Z");
    h->SetStats(0);
}

// Handsome Delta Phi Projection Histogram
void formatProjTH1F(TH1F * h, float offsetx, float offsety)
{
    h->SetTitleOffset(offsetx,"X");
    h->SetTitleOffset(offsety,"Y");
    h->GetXaxis()->CenterTitle();
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->CenterTitle();
    h->SetStats(0);
    
    h->GetXaxis()->SetTitle("#Delta#eta");
    h->GetYaxis()->SetTitle("#frac{1}{N^{trig}} #frac{dN^{pair}}{#delta#phi}");
}

/**************************************************************************************/
// Final Paper plotting
/**************************************************************************************/
int paperPlot(
              const std::string inFileName1,
              const int rows,
              const int columns
              )
{
    
    // ROOT Global setting
    TH1::SetDefaultSumw2();    TH2::SetDefaultSumw2();
    
    Selection s;
    TFile * f1 = TFile::Open(inFileName1.c_str(),"read");
    
    // ***** List of CFs to plot ***** //
    
    static const int num2PCs = 1;   // number of correlation functions
    static const int numProj = 7;   // number of projections
    double etaranges[numProj+1]={1.5,2.5,2.5,5,1.5,5,2.6,10};
    Int_t minbin,maxbin; // used for projections
    static const int nPads = rows * columns;
    // nPads should equal numProj but currently it does not
    
    int E[num2PCs] = {0};   // energy
    int m[num2PCs] = {2};  // multiplicity
    int pt[num2PCs] = {0};  // pt // 0.4 - 100
    int et[num2PCs] = {0};  // eta // Beam <1.6 // Thrust <4.5
    
    // ***** Load the histograms ***** //
    
    TH2F * signal2PC[num2PCs];
    TH2F * bkgrnd2PC[num2PCs];
    TH2F * ratio2PC[num2PCs];
    
    TH1D * nEvtSigHist1 = (TH1D*)f1->Get("nEvtSigHisto");
    TH1D * nEvtBkgHist1 = (TH1D*)f1->Get("nEvtSigHisto");
    
    TH1F * h_deltaphi[num2PCs][numProj];
    
    for (unsigned int i = 0; i < num2PCs; i++)
    {
        signal2PC[i] = (TH2F*)f1->Get(Form("signal2PC_%d_%d_%d_%d_%d",E[i],s.multBinsLow[m[i]],s.multBinsHigh[m[i]],pt[i],et[i]));
        std::cout<<signal2PC[i]->GetEntries()<<std::endl;
        if (signal2PC[i]->GetEntries()==0)
        {
            std::cout<<Form("BIG PROBLEM SIGNAL DISTRIBUTION %d HAS 0 ENTRIES",i)<<std::endl;
            break;
        }
        signal2PC[i]->Scale(1./nEvtSigHist1->GetBinContent(m[i]+1)); // plus 1 because 0 is the underflow bin
        
        bkgrnd2PC[i] = (TH2F*)f1->Get(Form("bkgrnd2PC_%d_%d_%d_%d_%d",E[i],s.multBinsLow[m[i]],s.multBinsHigh[m[i]],pt[i],et[i]));
        bkgrnd2PC[i]->Scale(1./nEvtBkgHist1->GetBinContent(m[i]+1));
        
        ratio2PC[i] = (TH2F*) signal2PC[i]->Clone(Form("ratio2PC_%d_%d_%d_%d_%d",E[i],s.multBinsLow[m[i]],s.multBinsHigh[m[i]],pt[i],et[i]));
        ratio2PC[i]->Reset();
        calculateRatio(signal2PC[i],bkgrnd2PC[i],ratio2PC[i]);
        formatCFTH2F(ratio2PC[i],1.5,1.5,2);    // format the plot
        
        // Perform 1D projections
        for (unsigned int j = 0; j < numProj; j++)
        {
            minbin =  ratio2PC[i]->GetXaxis()->FindBin(etaranges[j]);
            maxbin =  ratio2PC[i]->GetXaxis()->FindBin(etaranges[j+1]);
            
            h_deltaphi[i][j]  = (TH1F*) ratio2PC[i]->ProjectionY(Form("h_deltaphi%d_%d_%d_%d_%d_%d",j,E[i],s.multBinsLow[m[i]],s.multBinsHigh[m[i]],pt[i],et[i]),minbin,maxbin);
            h_deltaphi[i][j]->SetName(Form("h_deltaphi%d_%d_%d_%d_%d_%d",j,E[i],s.multBinsLow[m[i]],s.multBinsHigh[m[i]],pt[i],et[i]));
            h_deltaphi[i][j]->GetXaxis()->SetTitle("#Delta#phi");
            h_deltaphi[i][j]->SetTitle(Form("#Delta#phi, #Delta#eta (%1.1f, %1.1f), Multipliplicity (%d, %d)",etaranges[j],etaranges[j+1],s.multBinsLow[m[i]],s.multBinsHigh[m[i]]));
            h_deltaphi[i][j]->GetYaxis()->SetTitle("Y(#Delta#phi)");
            h_deltaphi[i][j]->Scale(1./(maxbin-minbin+1));
            formatProjTH1F(h_deltaphi[i][j],.8,1.5); // format the plot
        }
    }
    
    // ***** Plot to beautiful canvas ***** //
    for (unsigned int i = 0; i < num2PCs; i++)
    {
        TCanvas *CFViewer = new TCanvas("CFViewer","CFViewer",1000,1000);
        formatCFViewer(CFViewer);
        ratio2PC[i]->Draw("surf1 fb");
        plotLogo(1,1.6,1);
        CFViewer->SaveAs(Form("ratio2PC_%d_%d_%d_%d_%d.pdf",E[i],s.multBinsLow[m[i]],s.multBinsHigh[m[i]],pt[i],et[i]));
        
        TCanvas *multiPanel = new TCanvas("multiPanel","multiPanel",1000,1000);
        TPad* pads[nPads];
        formatMultipanelViewer(multiPanel, pads, rows, columns);
        
        TLatex* lat = new TLatex (0, 0, "hi");
        setTextAbovePad(lat,pads[0],0.38,0.2);
        lat->SetTextSize(0.05);
        
        for (unsigned int j = 0; j < numProj; j++)
        {
            multiPanel->cd(j+1);
            h_deltaphi[i][j]->Draw();
            lat->Draw();
        }
        multiPanel->SaveAs(Form("proj_%d_%d_%d_%d_%d.pdf",E[i],s.multBinsLow[m[i]],s.multBinsHigh[m[i]],pt[i],et[i]));
    }
    
    return 0;
}

/*
TH1F* test[nPads];



//TLatex* lat1 = new TLatex (0, 0, "hi");
//setTextBelowPad(lat1,pads[0],0.38,0.13);
//lat1->SetTextSize(0.05);

for (int i=0;i<nPads;i++)
{
    test[i] = new TH1F(Form("test%d",i),"",100,0.0,1.0);
    formatTH1F(test[i],.8,1.5);
    multiPanel->cd(i+1);
    test[i]->Draw();
    lat->Draw();
    //lat1->Draw();
}
*/








