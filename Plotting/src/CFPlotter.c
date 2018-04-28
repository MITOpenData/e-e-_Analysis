//
//  CFPlotter.c
//  
//
//  Created by Anthony Badea on 4/27/18.
//
//

// C dependencies
#include <stdio.h>

// root dependencies
#include <TCanvas.h>
#include <TVector3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>
#include <TMath.h>

// local dependencies
/**************************************************************************************/
// Handsome ROOT Style
/**************************************************************************************/

// Handsome Canvas for 2D histogram
TCanvas *CFViewer(const char *canvasName, const char *title, double x, double y)
{
    TCanvas *c = new TCanvas(canvasName,title,x,y);
    // ****** The MAGIC angle ****** //
    c->SetTheta(60.839);
    c->SetPhi(38.0172);
    return c;
}

// Handsome Correlation Function Histogram
void formatCFTH2F(TH2F *h,float offsetx, float offsety, float offsetz)
{
    h->Sumw2();
    h->SetTitleOffset(offset,"X");
    h->SetTitleOffset(offset,"Y");
    h->SetTitleOffset(offsetZ,"Z");
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
    h->GetZaxis()->CenterTitle();
    h->SetNdivisions(505,"X");
    h->SetNdivisions(505,"Y");
    h->SetNdivisions(505,"Z");
}
