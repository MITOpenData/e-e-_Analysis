#include <iostream>
#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TDatime.h"
#include "TLine.h"
#include "TLegend.h"

#include "include/kirchnerPalette.h"

void prettyHist(TH1* hist_p)
{
  hist_p->GetXaxis()->SetTitleFont(43);
  hist_p->GetYaxis()->SetTitleFont(43);
  hist_p->GetXaxis()->SetLabelFont(43);
  hist_p->GetYaxis()->SetLabelFont(43);

  hist_p->GetXaxis()->SetTitleSize(14);
  hist_p->GetYaxis()->SetTitleSize(14);
  hist_p->GetXaxis()->SetLabelSize(10);
  hist_p->GetYaxis()->SetLabelSize(10);

  hist_p->GetXaxis()->CenterTitle();
  hist_p->GetYaxis()->CenterTitle();

  hist_p->GetXaxis()->SetTitleOffset(hist_p->GetXaxis()->GetTitleOffset()*2.5);
  hist_p->GetYaxis()->SetTitleOffset(TMath::Max(1., (Double_t)hist_p->GetYaxis()->GetTitleOffset())*2.5);

  std::cout << hist_p->GetXaxis()->GetTitleOffset() << ", " << hist_p->GetYaxis()->GetTitleOffset() << std::endl;

  return;
}

int plotMultDist(const std::string inFileName)
{
  kirchnerPalette col;

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  //  TH1F* fig1L_p = (TH1F*)inFile_p->Get("fig1L_h");
  TH1F* fig1L_Corr_p = (TH1F*)inFile_p->Get("fig1L_Corr_h");
  TH1F* fig1L_Paper_p = (TH1F*)inFile_p->Get("fig1L_Paper_h");
  prettyHist(fig1L_Corr_p);
  prettyHist(fig1L_Paper_p);

  fig1L_Corr_p->SetMarkerStyle(20);
  fig1L_Paper_p->SetMarkerStyle(20);

  fig1L_Corr_p->SetMarkerSize(1);
  fig1L_Paper_p->SetMarkerSize(1);

  fig1L_Corr_p->SetMarkerColor(col.getColor(0));
  fig1L_Paper_p->SetMarkerColor(col.getColor(2));

  fig1L_Corr_p->SetLineColor(col.getColor(0));
  fig1L_Paper_p->SetLineColor(col.getColor(2));
  
  TCanvas* canv_p = new TCanvas("canv_c", "canv_c", 400, 500);
  canv_p->SetTopMargin(0.01);
  canv_p->SetLeftMargin(0.01);
  canv_p->SetBottomMargin(0.01);
  canv_p->SetRightMargin(0.01);

  canv_p->cd();

  const Double_t splitPoint = .33;

  TPad* pad1_p = new TPad("pad1", "pad1", 0.0, splitPoint, 1.0, 1.0);
  pad1_p->Draw();
  pad1_p->SetRightMargin(0.01);
  pad1_p->SetTopMargin(0.01);
  pad1_p->SetBottomMargin(0.0);
  pad1_p->SetLeftMargin(pad1_p->GetLeftMargin()*1.5);

  TPad* pad2_p = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, splitPoint);
  pad2_p->Draw();
  pad2_p->SetRightMargin(0.01);
  pad2_p->SetTopMargin(0.0);
  pad2_p->SetLeftMargin(pad1_p->GetLeftMargin());
  pad2_p->SetBottomMargin(pad1_p->GetLeftMargin()*(1-splitPoint)/splitPoint);

  canv_p->cd();
  pad1_p->cd();

  fig1L_Corr_p->SetMaximum(1000.);
  fig1L_Corr_p->SetMinimum(.001);
  fig1L_Corr_p->DrawCopy("HIST E1");
  fig1L_Paper_p->DrawCopy("HIST E1 SAME");

  gStyle->SetOptStat(0);
  gPad->SetLogy();

  TLegend* leg_p = new TLegend(0.6, 0.9, 0.6, 0.9);
  leg_p->SetBorderSize(0);
  leg_p->SetFillStyle(0);
  leg_p->SetFillColor(0);
  leg_p->AddEntry(fig1L_Paper_p, "Paper", "P L");
  leg_p->AddEntry(fig1L_Corr_p, "Reprocess", "P L");

  leg_p->Draw("SAME");

  canv_p->cd();
  pad2_p->cd();

  fig1L_Corr_p->Divide(fig1L_Paper_p);
  fig1L_Corr_p->GetYaxis()->SetTitle("Reprocess/Paper");
  
  Double_t min = fig1L_Corr_p->GetMaximum();
  Double_t max = fig1L_Corr_p->GetMinimum();

  for(Int_t i = 0; i < fig1L_Corr_p->GetNbinsX(); ++i){
    if(min > fig1L_Corr_p->GetBinContent(i+1)) min = fig1L_Corr_p->GetBinContent(i+1);
    if(max < fig1L_Corr_p->GetBinContent(i+1)) max = fig1L_Corr_p->GetBinContent(i+1);
  }

  Double_t interval = max - min;

  max += interval/10;
  min - interval/10 < 0 ? min = 0: min = min - interval/10;

  fig1L_Corr_p->SetMaximum(TMath::Max(1.5, max));
  fig1L_Corr_p->SetMinimum(TMath::Min(0.5, min));
  fig1L_Corr_p->DrawCopy("E1 P");
  gStyle->SetOptStat(0);

  TLine* line_p = new TLine();
  line_p->SetLineStyle(2);
  line_p->DrawLine(fig1L_Corr_p->GetBinLowEdge(1), 1, fig1L_Corr_p->GetBinLowEdge(fig1L_Corr_p->GetNbinsX()+1), 1);
  delete line_p;

  TDatime* date = new TDatime();
  canv_p->SaveAs(("pdfDir/f1L_multDist_" + std::to_string(date->GetDate()) + ".pdf").c_str());
  delete date;

  delete pad1_p;
  delete pad2_p;

  delete canv_p;

  inFile_p->Close();
  delete inFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage ./plotMultDist.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 2) retVal += plotMultDist(argv[1]);
  return retVal;
}
