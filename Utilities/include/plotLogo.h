#include <TPad.h>
#include <TImage.h>
#include <TStyle.h>
#include <iostream>

void plotLogo( Float_t v_scale, Float_t skew )
 {
 
    TImage *img = TImage::Open("../../AnalysisNote/ALEPH/images/Logo/MODLogo.jpg");
    if (!img) {
       cout << "+++ Could not open image MODLogo.gif" << endl;
       return;
    }
       
    img->SetConstRatio(kFALSE);
    UInt_t h_ = img->GetHeight();
    UInt_t w_ = img->GetWidth();
 
    Float_t r = w_/h_;
    gPad->Update();
    Float_t rpad = Double_t(gPad->VtoAbsPixel(0) - gPad->VtoAbsPixel(1))/(gPad->UtoAbsPixel(1) - gPad->UtoAbsPixel(0));
    r *= rpad;
 
    Float_t d = 0.055;
    // absolute coordinates
    Float_t x1R = 1 - gStyle->GetPadRightMargin(); 
    Float_t y1B = 1 - gStyle->GetPadTopMargin()+.01; // we like the logo to sit a bit above the histo 
 
    Float_t x1L = x1R - d*r/skew;
    Float_t y1T = y1B + d*v_scale*skew;
    if (y1T>0.99) y1T = 0.99;
 
    TPad *p1 = new TPad("imgpad", "imgpad", x1L, y1B, x1R, y1T );
    p1->SetRightMargin(0);
    p1->SetBottomMargin(0);
    p1->SetLeftMargin(0);
    p1->SetTopMargin(0);
    p1->Draw();
 
    Int_t xSizeInPixel = p1->UtoAbsPixel(1) - p1->UtoAbsPixel(0);
    Int_t ySizeInPixel = p1->VtoAbsPixel(0) - p1->VtoAbsPixel(1);
    cout <<x1L<<" "<<y1B<<" "<<x1R<<" "<<y1T<<" "<<endl;
    if (xSizeInPixel<=5 || ySizeInPixel<=5) {
       cout <<xSizeInPixel<<" "<<ySizeInPixel<<endl;
       delete p1;
       return; // ROOT doesn't draw smaller than this
    }
    p1->cd();
    img->Draw();
 } 
