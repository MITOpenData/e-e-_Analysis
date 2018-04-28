//
//  multipanelPlotter.c
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
#include "../include/styleUtil.h"

void formatTH1F(TH1F * h, float offsetx, float offsety)
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


int multipanelPlotter(const int rows, const int columns)
{
    // Parameters for the canvas division
    static const int nPads = rows * columns;
    static const float leftMargin = 0.1;
    static const float rightMargin = 0.1;
    static const float bottomMargin = 0.1;
    static const float topMargin = 0.1;
    static const float xMargin = 0.0;
    static const float yMargin = 0.0;
    static const float frameWidth = 0.8;
    static const float frameHeight = 0.8;
    static const float yMinOffset = 0.05;
    

    TCanvas* multPanel = new TCanvas("multPanel","multPanel",1200,1200);
    TPad* pads[nPads];
    divideCanvas(multPanel,pads,rows,columns,leftMargin,rightMargin,bottomMargin,topMargin,xMargin,yMargin,frameWidth,frameHeight,yMinOffset);
    TH1F* test[nPads];
    
    TLatex* lat = new TLatex (0, 0, "hi");
    setTextAbovePad(lat,pads[0],0.38,0.2);
    lat->SetTextSize(0.05);
    
    //TLatex* lat1 = new TLatex (0, 0, "hi");
    //setTextBelowPad(lat1,pads[0],0.38,0.13);
    //lat1->SetTextSize(0.05);
    
    for (int i=0;i<nPads;i++)
    {
        test[i] = new TH1F(Form("test%d",i),"",100,0.0,1.0);
        formatTH1F(test[i],.8,1.5);
        multPanel->cd(i+1);
        test[i]->Draw();
        lat->Draw();
        //lat1->Draw();
    }
    return 0;
}
