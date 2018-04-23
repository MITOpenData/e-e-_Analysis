#include "Plotting.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

Plotting::Plotting() {
}


TH2F*Plotting::GetEmpty(TString histoname,TString AxisTitleX,TString AxisTitleY,double min_X, double max_X,double min_Y, double max_Y){

  TH2F* hempty=new TH2F(histoname.Data(),histoname.Data(),50,min_X,max_X,10,min_Y,max_Y);  
  hempty->GetXaxis()->CenterTitle();
  hempty->GetYaxis()->CenterTitle();
  hempty->GetXaxis()->SetTitle(AxisTitleX.Data());
  hempty->GetYaxis()->SetTitle(AxisTitleY.Data());
  hempty->GetXaxis()->SetTitleOffset(0.9);
  hempty->GetYaxis()->SetTitleOffset(1.4);
  hempty->GetXaxis()->SetTitleSize(0.05);
  hempty->GetYaxis()->SetTitleSize(0.05);
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelSize(0.035);
  hempty->GetYaxis()->SetLabelSize(0.035);  
  return hempty;
}

TLegend* Plotting::GetLegend(double xmin,double xmax,double ymin,double ymax){

   TLegend *leg = new TLegend(xmin,xmax,ymin,ymax);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   return leg;
   
}
