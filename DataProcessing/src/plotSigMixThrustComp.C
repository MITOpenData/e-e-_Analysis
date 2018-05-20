#include <iostream>
#include <string>
#include <fstream>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TDatime.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"

#include "DataProcessing/include/checkMakeDir.h"
#include "DataProcessing/include/kirchnerPalette.h"
#include "DataProcessing/include/histDefUtility.h"

int plotSigMixThrustComp(const std::string inFileName)
{
  std::string dataMCStr = "Data";
  if(inFileName.find("MC") != std::string::npos) dataMCStr = "MC";

  if(inFileName.find("Pythia") != std::string::npos || inFileName.find("PYTHIA") != std::string::npos || inFileName.find("pythia") != std::string::npos) dataMCStr = "PYTHIA";

  const int nHist = 6;
  const std::string hist[nHist] = {"Pt", "Eta", "AbsEta", "Rap", "Theta", "Phi"};
  const float yPadFrac = .35;
  const float canvMarginFrac = .1;
  const float xVal2 = 0.5;
  const float yVal2 = 0.5/(1. - canvMarginFrac);
  const float padMarginFrac = .1*1./xVal2;

  const float xVal2Bottom = xVal2*(xVal2 + canvMarginFrac)/(xVal2);
  //  const float xVal2Bottom = xVal2/(1. - canvMarginFrac);
  const float padMarginFracBottom = .1*1./xVal2Bottom;

  const int nHist2 = 2;
  const std::string hist2[nHist2] = {"EtaPhi", "AbsEtaPhi"};

  const Int_t nStyles = 2;
  const Int_t styles[nStyles] = {20, 21};

  const Int_t nMultBins = 3;
  const Int_t multBinsLow[nMultBins] = {0, 20, 30};
  const Int_t multBinsHi[nMultBins] = {20, 30, 200};

  kirchnerPalette col;

  if(!checkFile(inFileName)){
    std::cout << "Given inFileName \'" << inFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  TDatime* date = new TDatime();
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");

  for(Int_t mI = 0; mI < nMultBins+1; ++mI){
    std::string multStr = "Mult" + std::to_string(multBinsLow[0]) + "to" + std::to_string(multBinsHi[nMultBins-1]);
    std::string multStr2 = std::to_string(multBinsLow[0]) + "< N_{particle} <" + std::to_string(multBinsHi[nMultBins-1]);
    if(mI < nMultBins){
      multStr = "Mult" + std::to_string(multBinsLow[mI]) + "to" + std::to_string(multBinsHi[mI]);    
      multStr2 = std::to_string(multBinsLow[mI]) + "< N_{particle} <" + std::to_string(multBinsHi[mI]);
    }    

    for(Int_t hI = 0; hI < nHist; ++hI){
      TH1D* hist1_p = (TH1D*)inFile_p->Get(("sigReco" + hist[hI] + "WrtRecoThr_" + multStr + "_h").c_str());
      TH1D* hist2_p = (TH1D*)inFile_p->Get(("mixReco" + hist[hI] + "WrtRecoThr_" + multStr + "_h").c_str());
      
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
    
      TLatex* label_p = new TLatex();
      label_p->SetTextFont(43);
      label_p->SetTextSize(12);
      label_p->SetNDC();



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
      
      Double_t overMin = 100;
      for(Int_t hI = 0; hI < hist1_p->GetNbinsX(); ++hI){
	if(hist1_p->GetBinContent(hI+1) > 0 && hist1_p->GetBinContent(hI+1) < overMin) overMin = hist1_p->GetBinContent(hI+1);
	if(hist2_p->GetBinContent(hI+1) > 0 && hist2_p->GetBinContent(hI+1) < overMin) overMin = hist2_p->GetBinContent(hI+1);
      }

      hist1_p->SetMaximum(TMath::Max(hist1_p->GetMaximum(), hist2_p->GetMaximum())*2.);
      hist1_p->SetMinimum(overMin/5.);


      hist1_p->DrawCopy("HIST E1 P");
      hist2_p->DrawCopy("HIST E1 P SAME");
      gStyle->SetOptStat(0);
      gPad->SetLogy();
      if(hist[hI].find("Pt") != std::string::npos) gPad->SetLogx();


      label_p->DrawLatex(0.3, 0.9, multStr2.c_str());
      leg_p->Draw("SAME");
      delete label_p;
      
      canv_p->cd();
      pads_p[1]->cd();

      if(hist[hI].find("AbsEta") != std::string::npos && mI < nMultBins){
	checkMakeDir("output");
	std::ofstream file;
	if(mI != 0) file.open(("output/signalOverMixAbsEtaTable_" + dataMCStr + "_" + std::to_string(date->GetDate()) + ".txt").c_str(), std::ios::app);
	else file.open(("output/signalOverMixAbsEtaTable_" + dataMCStr + "_" + std::to_string(date->GetDate()) + ".txt").c_str());

	file << "TABLE " << mI << "/" << nMultBins << ": " << multStr2 << std::endl;
	file << "AbsEtaLow,AbsEtaHigh,Numerator,Denominator,Error,RelativeError" << std::endl;
	for(Int_t bI = 0; bI < hist1_p->GetNbinsX(); ++bI){
	  Double_t val = hist1_p->GetBinContent(bI+1)/hist2_p->GetBinContent(bI+1);
	  Double_t err1 = hist1_p->GetBinError(bI+1)/hist1_p->GetBinContent(bI+1);
	  Double_t err2 = hist2_p->GetBinError(bI+1)/hist2_p->GetBinContent(bI+1);
	  Double_t err = TMath::Sqrt(err1*err1 + err2*err2);

	  file << hist1_p->GetBinLowEdge(bI+1) << ", " << hist1_p->GetBinLowEdge(bI+2) << ", " << hist1_p->GetBinContent(bI+1) << ", " << hist2_p->GetBinContent(bI+1) << ", " << val*err << ", "  << err << std::endl;
	}

	file.close();
      }
      
      hist1_p->Divide(hist2_p);

      hist1_p->SetMarkerColor(1);
      hist1_p->SetLineColor(1);
      
      hist1_p->GetYaxis()->SetTitle("#frac{Signal}{Mix}");
      hist1_p->SetMaximum(2.0);
      hist1_p->SetMinimum(0.0);
      hist1_p->GetYaxis()->SetNdivisions(505);

      if(hist[hI].find("Eta") != std::string::npos || hist[hI].find("AbsEta") != std::string::npos){
	Double_t max = 0;
	Double_t min = 2;
	for(Int_t hI = 0; hI < hist1_p->GetNbinsX(); ++hI){
	  if(hist1_p->GetBinContent(hI+1) > 0 && hist1_p->GetBinContent(hI+1) < min) min = hist1_p->GetBinContent(hI+1);
	  if(hist1_p->GetBinContent(hI+1) > max) max = hist1_p->GetBinContent(hI+1);
	}

	hist1_p->SetMinimum(min/5.);
	hist1_p->SetMaximum(max*5.);
      }

      hist1_p->DrawCopy("HIST E1 P");
      if(hist[hI].find("Pt") != std::string::npos) gPad->SetLogx();
      gStyle->SetOptStat(0);
      if(hist[hI].find("Eta") != std::string::npos) gPad->SetLogy();
      
      hist1_p->SetMarkerColor(col.getColor(0));
      hist1_p->SetLineColor(col.getColor(0));
      
      canv_p->SaveAs(("pdfDir/sigMixRecoThrustComp_" + hist[hI] + "_" + multStr + "_" + dataMCStr + "_" + std::to_string(date->GetDate()) + ".pdf").c_str());
      
      delete leg_p;
      delete pads_p[1];
      delete pads_p[0];
      delete canv_p;
    }

    for(Int_t hI = 0; hI < nHist2; ++hI){
      TH2D* hist1_p = (TH2D*)inFile_p->Get(("sigReco" + hist2[hI] + "WrtRecoThr_" + multStr + "_h").c_str());
      TH2D* hist2_p = (TH2D*)inFile_p->Get(("mixReco" + hist2[hI] + "WrtRecoThr_" + multStr + "_h").c_str());

      centerTitles({hist1_p, hist2_p});
      setSumW2({hist1_p, hist2_p});

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
      
      TLatex* label_p = new TLatex();
      label_p->SetTextFont(43);
      label_p->SetTextSize(12);
      label_p->SetNDC();

            
      TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
      canv_p->cd();
      TPad* pads_p[4];
      pads_p[0] = new TPad("pad0", "", 0.0, yVal2, xVal2, 1.0);
      pads_p[0]->SetTopMargin(0.01);
      pads_p[0]->SetRightMargin(0.001);
      pads_p[0]->SetBottomMargin(0.001);
      pads_p[0]->SetLeftMargin(padMarginFrac);
      pads_p[0]->Draw("SAME");
      pads_p[0]->cd();
      
      canv_p->cd();
      pads_p[1] = new TPad("pad1", "", xVal2, yVal2, 1.0, 1.0);
      pads_p[1]->SetTopMargin(0.01);
      pads_p[1]->SetRightMargin(padMarginFrac);
      pads_p[1]->SetBottomMargin(0.001);
      pads_p[1]->SetLeftMargin(0.001);
      pads_p[1]->Draw("SAME");
      pads_p[1]->cd();

      canv_p->cd();
      pads_p[2] = new TPad("pad2", "", 0.0, 0.0, xVal2Bottom, yVal2);
      pads_p[2]->SetTopMargin(0.001);
      pads_p[2]->SetRightMargin(padMarginFracBottom);
      pads_p[2]->SetBottomMargin(padMarginFrac);
      pads_p[2]->SetLeftMargin(padMarginFracBottom);
      pads_p[2]->Draw("SAME");
      pads_p[2]->cd();

      canv_p->cd();
      pads_p[3] = new TPad("pad3", "", xVal2Bottom, 0.0, 1.0, yVal2);

      pads_p[3]->SetTopMargin(0.001);
      pads_p[3]->SetRightMargin(0.01);
      pads_p[3]->SetBottomMargin(padMarginFrac);
      pads_p[3]->SetLeftMargin(0.001);
      pads_p[3]->Draw("SAME");
      pads_p[3]->cd();
      
      canv_p->cd();
      pads_p[0]->cd();
      
      Double_t overMin = 100;
      for(Int_t hIX = 0; hIX < hist1_p->GetNbinsX(); ++hIX){
	for(Int_t hIY = 0; hIY < hist1_p->GetNbinsX(); ++hIY){
	  if(hist1_p->GetBinContent(hIX+1, hIY+1) > 0 && hist1_p->GetBinContent(hIX+1, hIY+1) < overMin) overMin = hist1_p->GetBinContent(hIX+1, hIY+1);
	  if(hist2_p->GetBinContent(hIX+1, hIY+1) > 0 && hist2_p->GetBinContent(hIX+1, hIY+1) < overMin) overMin = hist2_p->GetBinContent(hIX+1, hIY+1);
	}
      }

      hist1_p->SetMaximum(TMath::Max(hist1_p->GetMaximum(), hist2_p->GetMaximum())*2.);
      hist1_p->SetMinimum(overMin/5.);
      hist1_p->GetYaxis()->SetNdivisions(505);
      hist2_p->GetYaxis()->SetNdivisions(505);

      hist1_p->GetYaxis()->SetTitle("#phi w.r.t. Thrust");

      hist1_p->DrawCopy("COL");
      gStyle->SetOptStat(0);
      gPad->SetLogz();

      label_p->DrawLatex(0.4, 0.8, "Signal");

      canv_p->cd();
      pads_p[1]->cd();

      hist2_p->SetMaximum(hist1_p->GetMaximum());
      hist2_p->SetMinimum(overMin/5.);

      hist2_p->DrawCopy("COLZ");
      gStyle->SetOptStat(0);
      gPad->SetLogz();

      label_p->DrawLatex(0.2, 0.8, "Mixed Event");
    
      if(hist2[hI].find("AbsEtaPhi") != std::string::npos && mI < nMultBins){
	checkMakeDir("output");
	std::ofstream file;
	if(mI != 0) file.open(("output/signalOverMixAbsEtaPhiTable_" + dataMCStr + "_" + std::to_string(date->GetDate()) + ".txt").c_str(), std::ios::app);
	else file.open(("output/signalOverMixAbsEtaPhiTable_" + dataMCStr + "_" + std::to_string(date->GetDate()) + ".txt").c_str());

	file << "TABLE " << mI << "/" << nMultBins << ": " << multStr2 << std::endl;
	file << "AbsEtaLow,AbsEtaHigh,PhiLow,PhiHigh,Numerator,Denominator,Error,RelativeError" << std::endl;
	for(Int_t bIX = 0; bIX < hist1_p->GetNbinsX(); ++bIX){
	  for(Int_t bIY = 0; bIY < hist1_p->GetNbinsY(); ++bIY){
	    Double_t val = hist1_p->GetBinContent(bIX+1, bIY+1)/hist2_p->GetBinContent(bIX+1, bIY+1);
	    Double_t err1 = hist1_p->GetBinError(bIX+1, bIY+1)/hist1_p->GetBinContent(bIX+1, bIY+1);
	    Double_t err2 = hist2_p->GetBinError(bIX+1, bIY+1)/hist2_p->GetBinContent(bIX+1, bIY+1);
	    Double_t err = TMath::Sqrt(err1*err1 + err2*err2);
	    
	    file << hist1_p->GetXaxis()->GetBinLowEdge(bIX+1) << ", " << hist1_p->GetXaxis()->GetBinLowEdge(bIX+2) << ", " << hist1_p->GetYaxis()->GetBinLowEdge(bIY+1) << ", " << hist1_p->GetYaxis()->GetBinLowEdge(bIY+2) << ", " << hist1_p->GetBinContent(bIX+1, bIY+1) << ", " << hist2_p->GetBinContent(bIX+1, bIY+1) << ", " << val*err << ", "  << err << std::endl;
	  }
	}

	file.close();
      }

      canv_p->cd();
      pads_p[2]->cd();
      
      hist1_p->Divide(hist2_p);
      
      hist1_p->GetYaxis()->SetTitle("#frac{Signal}{Mix}");
      hist1_p->SetMaximum(2.0);
      hist1_p->SetMinimum(0.0);
      hist1_p->GetYaxis()->SetNdivisions(505);

      /*
      if(hist[hI].find("Eta") != std::string::npos || hist[hI].find("AbsEta") != std::string::npos){
	Double_t max = 0;
	Double_t min = 2;
	for(Int_t hI = 0; hI < hist1_p->GetNbinsX(); ++hI){
	  if(hist1_p->GetBinContent(hI+1) > 0 && hist1_p->GetBinContent(hI+1) < min) min = hist1_p->GetBinContent(hI+1);
	  if(hist1_p->GetBinContent(hI+1) > max) max = hist1_p->GetBinContent(hI+1);
	}

	std::cout << "Min, max: " << min << ", " << max << std::endl;

	hist1_p->SetMinimum(min/5.);
	hist1_p->SetMaximum(max*5.);
      }
      */

      hist1_p->DrawCopy("COLZ");
      gStyle->SetOptStat(0);
      gPad->SetLogz();

      canv_p->SaveAs(("pdfDir/sigMixRecoThrustComp_" + hist2[hI] + "_" + multStr + "_" + dataMCStr + "_" + std::to_string(date->GetDate()) + ".pdf").c_str());

      delete pads_p[3];
      delete pads_p[2];
      delete pads_p[1];
      delete pads_p[0];

      delete canv_p;
    }
  }

  inFile_p->Close();
  delete inFile_p;

  delete date;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/plotSigMixThrustComp.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += plotSigMixThrustComp(argv[1]);
  return retVal;
}
