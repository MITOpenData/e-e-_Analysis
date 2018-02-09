//root dependencies
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TFile.h"
#include "TAxis.h"
#include "TAttMarker.h"
#include "TAttLine.h"
#include "TAttAxis.h"
#include "TColor.h"
#include "TF1.h"

//local headers
#include "include/Selection.h"
#include "include/utilities.h"

void formatTPCAxes(TH2F * h, float offsetx, float offsety, float offsetz){
   h->SetTitleOffset(offsetx,"X");
   h->SetTitleOffset(offsety,"Y");
   h->SetTitleOffset(offsetz,"Z");
   h->GetXaxis()->CenterTitle();
   h->GetYaxis()->CenterTitle();
   h->GetZaxis()->CenterTitle();
   h->SetNdivisions(505,"X");
   h->SetNdivisions(505,"Y");
   h->SetNdivisions(505,"Z");
   h->GetXaxis()->SetTitle("#Delta#eta");
   h->GetYaxis()->SetTitle("#Delta#phi");
   h->GetZaxis()->SetTitle("#frac{1}{N^{trig}} #frac{d^{2}N^{pair}}{d#Delta#etad#Delta#phi}");
   h->SetStats(0);
}
void formatTH1F(TH1F * h, float offsetx, float offsety){
   h->SetTitleOffset(offsetx,"X");
   h->SetTitleOffset(offsety,"Y");
   h->GetXaxis()->CenterTitle();
   h->GetYaxis()->CenterTitle();
   h->SetStats(0);
}

void formatZaxis(TH2F * h, bool minIsZero = true){
  float maximum = -1;
  float minimum = 99999;
  float averageSum = 0;
  float averageCounter = 0;
    
  // exclude the main peak that comes from self correlation of jets
  float dphi_peakMin = -1;
  float dphi_peakMax = 1;
  float deta_peakMin = -2;
  float deta_peakMax = 2;
    
  for(int i = 1; i<h->GetXaxis()->GetNbins()+1; i++){
    for(int j = 1; j<h->GetYaxis()->GetNbins()+1; j++){
        if( h->GetYaxis()->GetBinCenter(j) >= dphi_peakMin && h->GetYaxis()->GetBinCenter(j) <= dphi_peakMax && h->GetXaxis()->GetBinCenter(i) >= deta_peakMin && h->GetXaxis()->GetBinCenter(i) <= deta_peakMax) continue;
      if(h->GetBinContent(i,j)>=maximum) maximum = h->GetBinContent(i,j);
      if(h->GetBinContent(i,j)<=minimum) minimum = h->GetBinContent(i,j);
      averageSum += h->GetBinContent(i,j);
      averageCounter++;   
    }
  }
  float mean = averageSum/averageCounter;
  std::cout << maximum << " " << minimum << " " << averageSum/averageCounter << std::endl;

  h->GetZaxis()->SetRangeUser(minimum-0.2*(mean-minimum),maximum);
  //if(minIsZero==true) h->GetZaxis()->SetRangeUser(0,maximum);
  //else                h->GetZaxis()->SetRangeUser(minimum-0.2*(mean-minimum),mean+2.0*(mean-minimum));
}

void TPCPlots(const std::string inFileName1, const std::string inFileName2, const std::string dataName)
{
  Selection s;
    
  gStyle->SetLegendBorderSize(0);

  TH2F * sig1[s.nMultBins], *bkg1[s.nMultBins], *ratio1[s.nMultBins];
  TH2F * sig2[s.nMultBins], *bkg2[s.nMultBins], *ratio2[s.nMultBins];
  TH2F * r_sig[s.nMultBins], *r_bkg[s.nMultBins], *r_ratio[s.nMultBins];
  TH1F * longRangeYield1[s.nMultBins], *longRangeYield2[s.nMultBins];
  TH1F * r_longRangeYield[s.nMultBins];

  TFile * f1 = TFile::Open(inFileName1.c_str(),"read");
  TFile * f2 = TFile::Open(inFileName2.c_str(),"read");
   
  TH1D * nEvtSigHist1 = (TH1D*)f1->Get("nEvtSigHisto"); 
  TH1D * nEvtBkgHist1 = (TH1D*)f1->Get("nEvtSigHisto"); 

  TH1D * nEvtSigHist2 = (TH1D*)f2->Get("nEvtSigHisto"); 
  TH1D * nEvtBkgHist2 = (TH1D*)f2->Get("nEvtSigHisto"); 

  Float_t etaPlotRange = s.getEtaPlotRange();
  for(int i = 0; i<s.nMultBins; i++)
  {
    //if(i>1) continue;
    sig1[i] = (TH2F*)f1->Get(Form("signal2PC_%d_%d",s.multBinsLow[i],s.multBinsHigh[i]));
    sig1[i]->Scale(1./nEvtSigHist1->GetBinContent(i));
    bkg1[i] = (TH2F*)f1->Get(Form("bkgrnd2PC_%d_%d",s.multBinsLow[i],s.multBinsHigh[i]));
    bkg1[i]->Scale(1./nEvtBkgHist1->GetBinContent(i));
    ratio1[i] = new TH2F(Form("ratio2PC_%d_%d",s.multBinsLow[i],s.multBinsHigh[i]),";#Delta#eta;#Delta#Phi",s.dEtaBins,-etaPlotRange,etaPlotRange,s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
    ratio1[i]->Sumw2();
    longRangeYield1[i] = new TH1F(Form("longRangeYield_%d_%d",s.multBinsLow[i],s.multBinsHigh[i]),";#Delta#phi;Y(#Delta#Phi)",s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
    longRangeYield1[i]->Sumw2();
    calculateRatio(sig1[i],bkg1[i],ratio1[i]);
    getLongRangeYield(s,ratio1[i],longRangeYield1[i]);

    sig2[i] = (TH2F*)f2->Get(Form("signal2PC_%d_%d",s.multBinsLow[i],s.multBinsHigh[i]));
    sig2[i]->Scale(1./nEvtSigHist2->GetBinContent(i));
    bkg2[i] = (TH2F*)f2->Get(Form("bkgrnd2PC_%d_%d",s.multBinsLow[i],s.multBinsHigh[i]));
    bkg2[i]->Scale(1./nEvtBkgHist2->GetBinContent(i));
    ratio2[i] = new TH2F(Form("ratio2PC_%d_%d",s.multBinsLow[i],s.multBinsHigh[i]),";#Delta#eta;#Delta#Phi",s.dEtaBins,-etaPlotRange,etaPlotRange,s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
    ratio2[i]->Sumw2();
    longRangeYield2[i] = new TH1F(Form("longRangeYield_%d_%d",s.multBinsLow[i],s.multBinsHigh[i]),";#Delta#phi;Y(#Delta#Phi)",s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
    longRangeYield2[i]->Sumw2();
    calculateRatio(sig2[i],bkg2[i],ratio2[i]);
    getLongRangeYield(s,ratio2[i],longRangeYield2[i]);

    r_sig[i] = (TH2F*)sig1[i]->Clone(Form("%s_r_signal2PC_%d_%d",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));    r_sig[i]->Divide(sig2[i]);
    r_bkg[i] = (TH2F*)bkg1[i]->Clone(Form("%s_r_bkgrnd2PC_%d_%d",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));    r_bkg[i]->Divide(bkg2[i]);
    r_ratio[i] = (TH2F*)ratio1[i]->Clone(Form("%s_r_ratio2PC_%d_%d",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));   r_ratio[i]->Divide(ratio2[i]);
    r_longRangeYield[i] = (TH1F*)longRangeYield1[i]->Clone(Form("%s_r_longRangeYield_%d_%d",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i])); r_longRangeYield[i]->Divide(longRangeYield2[i]);

    TLegend * l = new TLegend(0.6,0.8,0.95,0.95);
    l->SetFillStyle(0);
    if(s.experiment==0)  l->AddEntry((TObject*)0,"ALEPH e^{+}e^{-}, #sqrt{s}= XXX GeV","");
    if(s.experiment==1)  l->AddEntry((TObject*)0,"DELPHI e^{+}e^{-}","");
    if(s.experiment==2)  l->AddEntry((TObject*)0,"Belle e^{+}e^{-}","");
    if(s.experiment==3)  l->AddEntry((TObject*)0,"CMS pp","");
    l->AddEntry((TObject*)0,Form("%d<=N^{trk}<%d",s.multBinsLow[i],s.multBinsHigh[i]),"");
    l->AddEntry((TObject*)0,Form("|#eta|<%1.1f",s.etaCut),"");
    l->AddEntry((TObject*)0,Form("%1.1f<p_{T}<%1.1f",s.trigPt[0],s.trigPt[1]),"");

    TCanvas * c1 = new TCanvas("c1","c1",800,800);
    c1->SetLeftMargin(0.2);
    c1->SetTheta(60.839);
    c1->SetPhi(38.0172);
    formatTPCAxes(sig1[i],1.5,1.5,2);
    formatZaxis(sig1[i],1);
    sig1[i]->Draw("surf1 fb");
    l->Draw("same");
    //sig[i]->GetZaxis()->SetRangeUser(0.9*sig[i]->GetMinimum(),0.4*sig[i]->GetMaximum());
    c1->SaveAs(Form("../pdfDir/%s_signal1_%d_%d.png",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    c1->SaveAs(Form("../pdfDir/%s_signal1_%d_%d.pdf",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    c1->SaveAs(Form("../pdfDir/%s_signal1_%d_%d.C",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
      
    formatTPCAxes(sig2[i],1.5,1.5,2);
    formatZaxis(sig2[i],1);
    sig2[i]->Draw("surf1 fb");
    l->Draw("same");
    //sig[i]->GetZaxis()->SetRangeUser(0.9*sig[i]->GetMinimum(),0.4*sig[i]->GetMaximum());
    c1->SaveAs(Form("../pdfDir/%s_signal2_%d_%d.png",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    c1->SaveAs(Form("../pdfDir/%s_signal2_%d_%d.pdf",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    c1->SaveAs(Form("../pdfDir/%s_signal2_%d_%d.C",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
      
    formatTPCAxes(r_sig[i],1.5,1.5,2);
    formatZaxis(r_sig[i],1);
    r_sig[i]->Draw("surf1 fb");
    l->Draw("same");
    //sig[i]->GetZaxis()->SetRangeUser(0.9*sig[i]->GetMinimum(),0.4*sig[i]->GetMaximum());
    c1->SaveAs(Form("../pdfDir/%s_r_signal_%d_%d.png",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    c1->SaveAs(Form("../pdfDir/%s_r_signal_%d_%d.pdf",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    c1->SaveAs(Form("../pdfDir/%s_r_signal_%d_%d.C",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    
    formatTPCAxes(bkg1[i],1.5,1.5,2);
    formatZaxis(bkg1[i],1);
    bkg1[i]->Draw("surf1 fb");
    l->Draw("same");
    //bkg[i]->GetZaxis()->SetRangeUser(0.9*bkg[i]->GetMinimum(),1.1*bkg[i]->GetMaximum());
    c1->SaveAs(Form("../pdfDir/%s_background1_%d_%d.png",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    c1->SaveAs(Form("../pdfDir/%s_background1_%d_%d.pdf",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    c1->SaveAs(Form("../pdfDir/%s_background1_%d_%d.C",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    
    formatTPCAxes(bkg2[i],1.5,1.5,2);
    formatZaxis(bkg2[i],1);
    bkg2[i]->Draw("surf1 fb");
    l->Draw("same");
    //bkg[i]->GetZaxis()->SetRangeUser(0.9*bkg[i]->GetMinimum(),1.1*bkg[i]->GetMaximum());
    c1->SaveAs(Form("../pdfDir/%s_background2_%d_%d.png",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    c1->SaveAs(Form("../pdfDir/%s_background2_%d_%d.pdf",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    c1->SaveAs(Form("../pdfDir/%s_background2_%d_%d.C",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    
    formatTPCAxes(r_bkg[i],1.5,1.5,2);
    formatZaxis(r_bkg[i],1);
    r_bkg[i]->Draw("surf1 fb");
    l->Draw("same");
    //bkg[i]->GetZaxis()->SetRangeUser(0.9*bkg[i]->GetMinimum(),1.1*bkg[i]->GetMaximum());
    c1->SaveAs(Form("../pdfDir/%s_r_background_%d_%d.png",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    c1->SaveAs(Form("../pdfDir/%s_r_background_%d_%d.pdf",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    c1->SaveAs(Form("../pdfDir/%s_r_background_%d_%d.C",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
      
    formatTPCAxes(ratio1[i],1.5,1.5,2);
    formatZaxis(ratio1[i],0);
    ratio1[i]->Draw("surf1 fb");
    l->Draw("same");
    //ratio[i]->GetZaxis()->SetRangeUser(0.9*ratio[i]->GetMinimum(),0.4*ratio[i]->GetMaximum());
    c1->SaveAs(Form("../pdfDir/%s_ratio1_%d_%d.png",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    c1->SaveAs(Form("../pdfDir/%s_ratio1_%d_%d.pdf",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    c1->SaveAs(Form("../pdfDir/%s_ratio1_%d_%d.C",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
 
    formatTPCAxes(ratio2[i],1.5,1.5,2);
    formatZaxis(ratio2[i],0);
    ratio2[i]->Draw("surf1 fb");
    l->Draw("same");
    //ratio[i]->GetZaxis()->SetRangeUser(0.9*ratio[i]->GetMinimum(),0.4*ratio[i]->GetMaximum());
    c1->SaveAs(Form("../pdfDir/%s_ratio2_%d_%d.png",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    c1->SaveAs(Form("../pdfDir/%s_ratio2_%d_%d.pdf",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    c1->SaveAs(Form("../pdfDir/%s_ratio2_%d_%d.C",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
      
    formatTPCAxes(r_ratio[i],1.5,1.5,2);
    formatZaxis(r_ratio[i],0);
    r_ratio[i]->Draw("surf1 fb");
    l->Draw("same");
    //ratio[i]->GetZaxis()->SetRangeUser(0.9*ratio[i]->GetMinimum(),0.4*ratio[i]->GetMaximum());
    c1->SaveAs(Form("../pdfDir/%s_r_ratio_%d_%d.png",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    c1->SaveAs(Form("../pdfDir/%s_r_ratio_%d_%d.pdf",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    c1->SaveAs(Form("../pdfDir/%s_r_ratio_%d_%d.C",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
      
    formatTH1F(longRangeYield1[i],1.5,1.5);
    longRangeYield1[i]->Draw();
    l->Draw("same");
    c1->SaveAs(Form("../pdfDir/%s_r_longRangeYield1_%d_%d.png",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    c1->SaveAs(Form("../pdfDir/%s_r_longRangeYield1_%d_%d.pdf",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    c1->SaveAs(Form("../pdfDir/%s_r_longRangeYield1_%d_%d.C",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));

    formatTH1F(longRangeYield2[i],1.5,1.5);
    longRangeYield2[i]->Draw();
    l->Draw("same");
    c1->SaveAs(Form("../pdfDir/%s_r_longRangeYield2_%d_%d.png",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    c1->SaveAs(Form("../pdfDir/%s_r_longRangeYield2_%d_%d.pdf",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));
    c1->SaveAs(Form("../pdfDir/%s_r_longRangeYield2_%d_%d.C",dataName.c_str(),s.multBinsLow[i],s.multBinsHigh[i]));



    delete c1;
    delete l;
  }

  TH1F *eta1, *theta1, *pt1, *phi1, *TTheta1, *TPhi1;
  pt1 = (TH1F*)f1->Get("pt");
  eta1 = (TH1F*)f1->Get("eta");
  theta1 = (TH1F*)f1->Get("theta");
  phi1 = (TH1F*)f1->Get("phi");
  TTheta1 = (TH1F*)f1->Get("T_theta");
  TPhi1 = (TH1F*)f1->Get("T_phi");
    
  TH1F *eta2, *theta2, *pt2, *phi2, *TTheta2, *TPhi2;
  pt2 = (TH1F*)f2->Get("pt");
  eta2 = (TH1F*)f2->Get("eta");
  theta2 = (TH1F*)f2->Get("theta");
  phi2 = (TH1F*)f2->Get("phi");
  TTheta2 = (TH1F*)f2->Get("T_theta");
  TPhi2 = (TH1F*)f2->Get("T_phi");

  TCanvas * c2 = new TCanvas("c2","c2",800,800);
    
  eta1->Draw("p");
  eta1->SetTitle(";#eta;N");
  eta1->SetStats(0);
  c2->SaveAs(Form("../pdfDir/%s_eta1.png",dataName.c_str()));
  c2->SaveAs(Form("../pdfDir/%s_eta1.pdf",dataName.c_str()));
  c2->SaveAs(Form("../pdfDir/%s_eta1.C",dataName.c_str()));
  eta2->Draw("p");
  eta2->SetTitle(";#eta;N");
  eta2->SetStats(0);
  c2->SaveAs(Form("../pdfDir/%s_eta2.png",dataName.c_str()));
  c2->SaveAs(Form("../pdfDir/%s_eta2.pdf",dataName.c_str()));
  c2->SaveAs(Form("../pdfDir/%s_eta2.C",dataName.c_str()));
    
  phi1->Draw("p");
  phi1->SetTitle(";#phi;N");
  phi1->SetStats(0);
  c2->SaveAs(Form("../pdfDir/%s_phi1.png",dataName.c_str()));
  c2->SaveAs(Form("../pdfDir/%s_phi1.pdf",dataName.c_str()));
  c2->SaveAs(Form("../pdfDir/%s_phi1.C",dataName.c_str()));
  phi2->Draw("p");
  phi2->SetTitle(";#phi;N");
  phi2->SetStats(0);
  c2->SaveAs(Form("../pdfDir/%s_phi2.png",dataName.c_str()));
  c2->SaveAs(Form("../pdfDir/%s_phi2.pdf",dataName.c_str()));
  c2->SaveAs(Form("../pdfDir/%s_phi2.C",dataName.c_str()));
    
  theta1->Draw("p");
  theta1->SetTitle(";#theta;N");
  theta1->SetStats(0);
  c2->SaveAs(Form("../pdfDir/%s_theta1.png",dataName.c_str()));
  c2->SaveAs(Form("../pdfDir/%s_theta1.pdf",dataName.c_str()));
  c2->SaveAs(Form("../pdfDir/%s_theta1.C",dataName.c_str()));
  theta2->Draw("p");
  theta2->SetTitle(";#theta;N");
  theta2->SetStats(0);
  c2->SaveAs(Form("../pdfDir/%s_theta2.png",dataName.c_str()));
  c2->SaveAs(Form("../pdfDir/%s_theta2.pdf",dataName.c_str()));
  c2->SaveAs(Form("../pdfDir/%s_theta2.C",dataName.c_str()));
    
  TTheta1->Draw("p");
  TTheta1->SetTitle(";#theta_{thrust};N");
  TTheta1->SetStats(0);
  c2->SaveAs(Form("../pdfDir/%s_TTheta1.png",dataName.c_str()));
  c2->SaveAs(Form("../pdfDir/%s_TTheta1.pdf",dataName.c_str()));
  c2->SaveAs(Form("../pdfDir/%s_TTheta1.C",dataName.c_str()));
  TTheta2->Draw("p");
  TTheta2->SetTitle(";#theta_{thrust};N");
  TTheta2->SetStats(0);
  c2->SaveAs(Form("../pdfDir/%s_TTheta2.png",dataName.c_str()));
  c2->SaveAs(Form("../pdfDir/%s_TTheta2.pdf",dataName.c_str()));
  c2->SaveAs(Form("../pdfDir/%s_TTheta2.C",dataName.c_str()));
    
  TPhi1->Draw("p");
  TPhi1->SetTitle(";#phi_{thrust};N");
  TPhi1->SetStats(0);
  c2->SaveAs(Form("../pdfDir/%s_TPhi1.png",dataName.c_str()));
  c2->SaveAs(Form("../pdfDir/%s_TPhi1.pdf",dataName.c_str()));
  c2->SaveAs(Form("../pdfDir/%s_TPhi1.C",dataName.c_str()));
  TPhi2->Draw("p");
  TPhi2->SetTitle(";#phi_{thrust};N");
  TPhi2->SetStats(0);
  c2->SaveAs(Form("../pdfDir/%s_TPhi2.png",dataName.c_str()));
  c2->SaveAs(Form("../pdfDir/%s_TPhi2.pdf",dataName.c_str()));
  c2->SaveAs(Form("../pdfDir/%s_TPhi2.C",dataName.c_str()));

  c2->SetLogy();
  pt1->Draw("p");
  pt1->SetTitle(";p_{T};N");
  pt1->SetStats(0);
  c2->SaveAs(Form("../pdfDir/%s_pt1.png",dataName.c_str()));
  c2->SaveAs(Form("../pdfDir/%s_pt1.pdf",dataName.c_str()));
  c2->SaveAs(Form("../pdfDir/%s_pt1.C",dataName.c_str()));
  c2->SetLogy();
  pt2->Draw("p");
  pt2->SetTitle(";p_{T};N");
  pt2->SetStats(0);
  c2->SaveAs(Form("../pdfDir/%s_pt2.png",dataName.c_str()));
  c2->SaveAs(Form("../pdfDir/%s_pt2.pdf",dataName.c_str()));
  c2->SaveAs(Form("../pdfDir/%s_pt2.C",dataName.c_str()));
    
  delete c2;
    
  // Cross Checks
}
