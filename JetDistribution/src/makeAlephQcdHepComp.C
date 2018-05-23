//cpp dependencies
#include <iostream>
#include <string>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TDatime.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"

//Non-Local StudyMult (DataProcessing, etc.) dependencies
#include "DataProcessing/include/checkMakeDir.h"
#include "DataProcessing/include/histDefUtility.h"
#include "Plotting/include/vanGoghPalette.h"

//Local StudyMult (JetDistribution) dependencies
#include "JetDistribution/include/hepPlainTxtReader.h"
#include "JetDistribution/include/globalJetVar.h"

int makeAlephQcdHepComp(const std::string inFileName)
{
  if(!checkFile(inFileName)){
    std::cout << "Given inFileName \'" << inFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }
  else if(inFileName.find(".root") == std::string::npos){
    std::cout << "Given inFileName \'" << inFileName << "\' is not a root file. return 1" << std::endl;
    return 1;
  }

  vanGoghPalette vg;

  std::string outFileName = inFileName;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  outFileName.replace(outFileName.find(".root"), 5, "");
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  outFileName = "output/" + outFileName + "_AlephQcdHepComp_" + dateStr + ".root";
  delete date;

  const Int_t nHist = 8;
  const std::string histStr[nHist] = {"majorThrust", "minorThrust", "oblateness", "sphericity", "aplanarity", "planarity", "massDiff", "jetResParam4NegLog"};
  const std::string tableStr[nHist] = {"/HepData/6834/d101-x1-y1", "/HepData/6834/d109-x1-y1", "/HepData/6834/d140-x1-y1", "/HepData/6834/d148-x1-y1", "/HepData/6834/d125-x1-y1", "/HepData/6834/d132-x1-y1", "/HepData/6834/d117-x1-y1", "/HepData/6834/d172-x1-y1"};
  const std::string xStr[nHist] = {"Thrust Major", "Thrust Minor", "Oblateness", "Sphericity", "Aplanarity", "Planarity", "jet mass difference M_{D}", "4-jet resolution parameter -ln y_{4}"};
  const std::string yStr[nHist] = {"1/#sigma d#sigma/dTmajor", "1/#sigma d#sigma/dTminor", "1/#sigma d#sigma/dO", "1/#sigma d#sigma/dS", "1/#sigma d#sigma/dA", "1/#sigma d#sigma/dP", "1/#sigma d#sigma/dM_{D}" "1/#sigma d#sigma/d-ln y_{4}"};

  TH1D* alephQCDHist_p[nHist];
  TH1D* studyMultHist_p[nHist];
  
  hepPlainTxtReader hep;
  const std::string fullPath =  std::getenv("STUDYMULTDIR");

  for(Int_t hI = 0; hI < nHist; ++hI){
    alephQCDHist_p[hI] = NULL;

    hep.Init(std::string(fullPath + "/JetDistribution/tables/hepDataALEPHStudiesOfQCD2004_10_1140_epjc_s2004_01891_4.txt"), tableStr[hI]);
    hep.Print();

    hep.GetHistogram((&(alephQCDHist_p[hI])), std::string(histStr[hI] + "_AlephQcdHist_h"));

    const Int_t nBins = alephQCDHist_p[hI]->GetNbinsX();
    Double_t bins[nBins+1];
    for(Int_t bI = 0; bI < nBins+1; ++bI){
      bins[bI] = alephQCDHist_p[hI]->GetBinLowEdge(bI+1);
    }

    studyMultHist_p[hI] = new TH1D((histStr[hI] + "_StudyMultHist_h").c_str(), ";;", nBins, bins);

    alephQCDHist_p[hI]->GetXaxis()->SetTitle(xStr[hI].c_str());
    alephQCDHist_p[hI]->GetYaxis()->SetTitle(yStr[hI].c_str());
    alephQCDHist_p[hI]->GetXaxis()->SetTitleFont(43);
    alephQCDHist_p[hI]->GetXaxis()->SetLabelFont(43);
    alephQCDHist_p[hI]->GetYaxis()->SetTitleFont(43);
    alephQCDHist_p[hI]->GetYaxis()->SetLabelFont(43);
    alephQCDHist_p[hI]->GetXaxis()->SetTitleSize(16);
    alephQCDHist_p[hI]->GetXaxis()->SetLabelSize(16);
    alephQCDHist_p[hI]->GetYaxis()->SetTitleSize(16);
    alephQCDHist_p[hI]->GetYaxis()->SetLabelSize(16);

    studyMultHist_p[hI]->GetXaxis()->SetTitle(xStr[hI].c_str());
    studyMultHist_p[hI]->GetYaxis()->SetTitle(yStr[hI].c_str());
    studyMultHist_p[hI]->GetXaxis()->SetTitleFont(43);
    studyMultHist_p[hI]->GetXaxis()->SetLabelFont(43);
    studyMultHist_p[hI]->GetYaxis()->SetTitleFont(43);
    studyMultHist_p[hI]->GetYaxis()->SetLabelFont(43);
    studyMultHist_p[hI]->GetXaxis()->SetTitleSize(16);
    studyMultHist_p[hI]->GetXaxis()->SetLabelSize(16);
    studyMultHist_p[hI]->GetYaxis()->SetTitleSize(16);
    studyMultHist_p[hI]->GetYaxis()->SetLabelSize(16);

    alephQCDHist_p[hI]->SetMarkerColor(vg.getColor(0));
    alephQCDHist_p[hI]->SetMarkerSize(0.8);
    alephQCDHist_p[hI]->SetMarkerStyle(20);
    alephQCDHist_p[hI]->SetLineColor(vg.getColor(0));

    studyMultHist_p[hI]->SetMarkerColor(vg.getColor(1));
    studyMultHist_p[hI]->SetMarkerSize(0.8);
    studyMultHist_p[hI]->SetMarkerStyle(47);
    studyMultHist_p[hI]->SetLineColor(vg.getColor(1));
    
    centerTitles({alephQCDHist_p[hI], studyMultHist_p[hI]});
    setSumW2({alephQCDHist_p[hI], studyMultHist_p[hI]});
  }

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* inTree_p = (TTree*)inFile_p->Get("generalJetTree_Global");
  globalJetVar jetVar;
  jetVar.SetStatusAndAddressRead(inTree_p, {});

  const Int_t nEntries = inTree_p->GetEntries();

  Int_t validEvents = 0;

  std::cout << "Processing " << nEntries << " entries..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%10000 == 0) std::cout << " Entry: " << entry << "/" << nEntries << std::endl;
    inTree_p->GetEntry(entry);
    
    if(jetVar.energy < 203.) continue;
    if(jetVar.energy > 209.) continue;

    ++validEvents;

    Double_t varVals[nHist] = {jetVar.thrustMajorMag, jetVar.thrustMinorMag, jetVar.oblateness, jetVar.sphericity, jetVar.aplanarity, jetVar.planarity, jetVar.jetMassDifference, jetVar.jetResParam4NegLog};
    
    for(Int_t hI = 0; hI < nHist; ++hI){
      if(studyMultHist_p[hI]->GetBinLowEdge(1) > varVals[hI] || varVals[hI] > studyMultHist_p[hI]->GetBinLowEdge(studyMultHist_p[hI]->GetNbinsX()+1)) continue;
      studyMultHist_p[hI]->Fill(varVals[hI]);
    }
  }

  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  for(Int_t hI = 0; hI < nHist; ++hI){
    alephQCDHist_p[hI]->Write("", TObject::kOverwrite);

    studyMultHist_p[hI]->Scale(1./validEvents);
    for(Int_t bI = 0; bI < studyMultHist_p[hI]->GetNbinsX(); ++bI){
      studyMultHist_p[hI]->SetBinContent(bI+1, studyMultHist_p[hI]->GetBinContent(bI+1)/studyMultHist_p[hI]->GetBinWidth(bI+1));
      studyMultHist_p[hI]->SetBinError(bI+1, studyMultHist_p[hI]->GetBinError(bI+1)/studyMultHist_p[hI]->GetBinWidth(bI+1));
    }

    studyMultHist_p[hI]->Write("", TObject::kOverwrite);
  }

  TLatex* label_p = new TLatex();
  label_p->SetTextFont(43);
  label_p->SetTextSize(16);
  label_p->SetNDC();

  TLegend* leg_p = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg_p->SetTextFont(43);
  leg_p->SetTextSize(16);
  leg_p->SetBorderSize(0);
  leg_p->SetFillStyle(0);
  leg_p->SetFillColor(0);
  leg_p->AddEntry(studyMultHist_p[0], "StudyMult", "P L");
  leg_p->AddEntry(alephQCDHist_p[0], "ALEPH QCD", "P L");

  checkMakeDir("pdfDir");
  for(Int_t hI = 0; hI < nHist; ++hI){
    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
    canv_p->SetTopMargin(0.01);
    canv_p->SetRightMargin(0.01);
    canv_p->SetLeftMargin(0.12);
    canv_p->SetBottomMargin(canv_p->GetLeftMargin());

    canv_p->cd();

    studyMultHist_p[hI]->GetXaxis()->SetNdivisions(505);
    studyMultHist_p[hI]->GetYaxis()->SetNdivisions(505);

    double minVal = 9999999;
    double maxVal = 0;
    for(Int_t bI = 0; bI < studyMultHist_p[hI]->GetNbinsX(); ++bI){
      if(studyMultHist_p[hI]->GetBinContent(bI+1) > maxVal) maxVal = studyMultHist_p[hI]->GetBinContent(bI+1);
      if(alephQCDHist_p[hI]->GetBinContent(bI+1) > maxVal) maxVal = alephQCDHist_p[hI]->GetBinContent(bI+1);

      if(studyMultHist_p[hI]->GetBinContent(bI+1) < minVal && studyMultHist_p[hI]->GetBinContent(bI+1) > 0.) minVal = studyMultHist_p[hI]->GetBinContent(bI+1);
      if(alephQCDHist_p[hI]->GetBinContent(bI+1) < minVal && alephQCDHist_p[hI]->GetBinContent(bI+1) > 0.) minVal = alephQCDHist_p[hI]->GetBinContent(bI+1);
    }

    studyMultHist_p[hI]->SetMaximum(maxVal*5.);
    studyMultHist_p[hI]->SetMinimum(minVal/5.);

    studyMultHist_p[hI]->DrawCopy("HIST E1 P");
    alephQCDHist_p[hI]->DrawCopy("HIST E1 P SAME");

    leg_p->Draw("SAME");
    label_p->DrawLatex(.3, .9, "203 < E_{CM} < 209");

    gPad->RedrawAxis();
    gPad->SetLogy();
    gStyle->SetOptStat(0);

    const std::string saveName = "pdfDir/" + histStr[hI] + "_AlephComp_" + dateStr + ".pdf";
    canv_p->SaveAs(saveName.c_str());
    delete canv_p;
  }

  delete label_p;
  delete leg_p;

  for(Int_t hI = 0; hI < nHist; ++hI){
    delete alephQCDHist_p[hI];
    delete studyMultHist_p[hI];
  }

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/makeAlephQcdHepComp.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += makeAlephQcdHepComp(argv[1]);
  return retVal;
}
