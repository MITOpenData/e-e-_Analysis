#include <iostream>
#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TDatime.h"
#include "TStyle.h"
#include "TLegend.h"

#include "DataProcessing/include/checkMakeDir.h"
#include "DataProcessing/include/kirchnerPalette.h"
#include "DataProcessing/include/histDefUtility.h"

int plotRecoGenThrustComp(const std::string inFileName)
{
  const int nHist = 5;
  const std::string hist[nHist] = {"Pt", "Eta", "Rap", "Theta", "Phi"};
  const float yPadFrac = .35;

  const Int_t nStyles = 2;
  const Int_t styles[nStyles] = {20, 21};

  kirchnerPalette col;

  if(!checkFile(inFileName)){
    std::cout << "Given inFileName \'" << inFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  TDatime* date = new TDatime();
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");

  for(Int_t hI = 0; hI < nHist; ++hI){
    TH1F* hist1_p = (TH1F*)inFile_p->Get(("sigReco" + hist[hI] + "WrtRecoThr_h").c_str());
    TH1F* hist2_p = (TH1F*)inFile_p->Get(("mixReco" + hist[hI] + "WrtRecoThr_h").c_str());

    centerTitles({hist1_p, hist2_p});
    setSumW2({hist1_p, hist2_p});

    hist1_p->SetMarkerStyle(styles[0]); 
    hist1_p->SetMarkerSize(0.6);
    hist1_p->SetMarkerColor(col.getColor(0));
    hist1_p->SetLineColor(col.getColor(0));

    hist1_p->GetXaxis()->SetTitleFont(43);
    hist1_p->GetYaxis()->SetTitleFont(43);
    hist1_p->GetXaxis()->SetLabelFont(43);
    hist1_p->GetYaxis()->SetLabelFont(43);
    hist1_p->GetXaxis()->SetTitleSize(14);
    hist1_p->GetYaxis()->SetTitleSize(14);
    hist1_p->GetXaxis()->SetLabelSize(12);
    hist1_p->GetYaxis()->SetLabelSize(12);

    hist1_p->GetXaxis()->SetTitleOffset(hist1_p->GetXaxis()->GetTitleOffset()*1.5);
    hist1_p->GetYaxis()->SetTitleOffset(hist1_p->GetXaxis()->GetTitleOffset());
    hist1_p->GetXaxis()->SetTitleOffset(hist1_p->GetXaxis()->GetTitleOffset()*2.);

    hist2_p->SetMarkerStyle(styles[1]); 
    hist2_p->SetMarkerSize(0.6);
    hist2_p->SetMarkerColor(col.getColor(2));
    hist2_p->SetLineColor(col.getColor(2));

    TLegend* leg_p = new TLegend(0.65, 0.7, 0.9, 0.9);
    leg_p->SetBorderSize(0);
    leg_p->SetFillColor(0);
    leg_p->SetFillStyle(0);
    leg_p->SetTextFont(43);
    leg_p->SetTextSize(14);

    leg_p->AddEntry(hist1_p, "Reco. Signal", "P L");
    leg_p->AddEntry(hist2_p, "Reco. Mix", "P L");

    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
    canv_p->cd();
    TPad* pads_p[2];
    pads_p[0] = new TPad("pad0", "", 0.0, yPadFrac, 1.0, 1.0);
    pads_p[0]->SetTopMargin(0.01);
    pads_p[0]->SetRightMargin(0.01);
    pads_p[0]->SetBottomMargin(0.001);
    pads_p[0]->SetLeftMargin(pads_p[0]->GetLeftMargin()*1.3);
    pads_p[0]->Draw("SAME");
    pads_p[0]->cd();

    canv_p->cd();
    pads_p[1] = new TPad("pad1", "", 0.0, 0.0, 1.0, yPadFrac);
    pads_p[1]->SetTopMargin(0.001);
    pads_p[1]->SetRightMargin(0.01);
    pads_p[1]->SetBottomMargin(pads_p[0]->GetLeftMargin()/yPadFrac);
    pads_p[1]->SetLeftMargin(pads_p[0]->GetLeftMargin());
    pads_p[1]->Draw("SAME");
    pads_p[1]->cd();

    canv_p->cd();
    pads_p[0]->cd();

    hist1_p->DrawCopy("HIST E1 P");
    hist2_p->DrawCopy("HIST E1 P SAME");
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    if(hist[hI].find("Pt") != std::string::npos) gPad->SetLogx();

    leg_p->Draw("SAME");

    canv_p->cd();
    pads_p[1]->cd();

    hist1_p->Divide(hist2_p);
    hist1_p->SetMarkerColor(1);
    hist1_p->SetLineColor(1);
    
    hist1_p->GetYaxis()->SetTitle("#frac{Reco.}{Gen.} Thrust");
    hist1_p->SetMaximum(1.15);
    hist1_p->SetMinimum(0.85);
    hist1_p->GetYaxis()->SetNdivisions(505);

    hist1_p->DrawCopy("HIST E1 P");
    if(hist[hI].find("Pt") != std::string::npos) gPad->SetLogx();
    gStyle->SetOptStat(0);

    hist1_p->SetMarkerColor(col.getColor(0));
    hist1_p->SetLineColor(col.getColor(0));

    canv_p->SaveAs(("pdfDir/sigMixRecoThrustComp_" + hist[hI] + "_" + std::to_string(date->GetDate()) + ".pdf").c_str());

    delete leg_p;
    delete pads_p[1];
    delete pads_p[0];
    delete canv_p;
  }

  inFile_p->Close();
  delete inFile_p;

  delete date;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/plotRecoGenThrustComp.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += plotRecoGenThrustComp(argv[1]);
  return retVal;
}
