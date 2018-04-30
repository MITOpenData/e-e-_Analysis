//
//  formatHists.h
//  
//
//  Created by Anthony Badea on 4/30/18.
//
//

#ifndef formatHists_h
#define formatHists_h



// C dependencies
#include <stdio.h>

// root dependencies
#include <TCanvas.h>
#include <TVector3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>
#include <TMath.h>

/**************************************************************************************/
// General purpose histogram formatting
/**************************************************************************************/

// nice formatTH1(h,.8,1.5);
void formatTH1(TH1F * h, float offsetx, float offsety)
{
    //h->Sumw2();
    h->SetTitleOffset(offsetx,"X");
    h->SetTitleOffset(offsety,"Y");
    h->GetXaxis()->CenterTitle();
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->CenterTitle();
    h->SetStats(0);
}

void formatTH1(TH1D * h, float offsetx, float offsety)
{
    //h->Sumw2();
    h->SetTitleOffset(offsetx,"X");
    h->SetTitleOffset(offsety,"Y");
    h->GetXaxis()->CenterTitle();
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->CenterTitle();
    h->SetStats(0);
}

#endif /* formatHists_h */
