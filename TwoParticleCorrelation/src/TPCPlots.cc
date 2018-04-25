//root dependencies
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TVirtualPad.h"
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
#include "../include/xjjrootuti.h"

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
   h->GetXaxis()->SetTitleSize(0.05);
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

void setupCanvas(TVirtualPad *c1)
{
          c1->SetLeftMargin(0.2);
          c1->SetTheta(60.839);
          c1->SetPhi(38.0172);
}

int TPCPlots(const std::string inFileName1, const std::string outfilename, const std::string plotsname, int doOneBin=0)
{
  // ROOT Global setting
  TH1::SetDefaultSumw2();    TH2::SetDefaultSumw2();

  Selection s;
    
  gStyle->SetLegendBorderSize(0);

  /// Initialize the histograms
  std::cout<<"Initializing histograms..."<<std::endl;
  Int_t nEnergyBins = s.nEnergyBins;
  Int_t nMultBins = s.nMultBins;
  Int_t nptBins = s.nptBins;
  Int_t netaBins = s.netaBins;
  Int_t experiment = s.experiment;
  bool doThrust = s.doThrust;
  bool doWTA = s.doWTA;
  bool doPerp = s.doPerp;
  bool doGen = s.doGen;
  bool getThetaAngle = s.getThetaAngle;
  

  TH2F * signal2PC[nEnergyBins][nMultBins][nptBins][netaBins];
  TH2F * bkgrnd2PC[nEnergyBins][nMultBins][nptBins][netaBins];
  TH2F * ratio2PC[nEnergyBins][nMultBins][nptBins][netaBins];

  TFile * f1 = TFile::Open(inFileName1.c_str(),"read");
  TFile * fout = TFile::Open(outfilename.c_str(),"recreate");
  
  TH1F *hMetaData = (TH1F*) f1->Get("hMetaData");
  if (hMetaData!=0){
     cout <<"Read configuration from data file:"<<endl;
     experiment = hMetaData->GetBinContent(1);
     doThrust   = hMetaData->GetBinContent(2);
     doWTA      = hMetaData->GetBinContent(3);
     doPerp     = hMetaData->GetBinContent(4);
     doGen      = hMetaData->GetBinContent(5);
     getThetaAngle = hMetaData->GetBinContent(6);
     nEnergyBins = hMetaData->GetBinContent(7);
     nMultBins   = hMetaData->GetBinContent(8);
     nptBins     = hMetaData->GetBinContent(9);
     netaBins    = hMetaData->GetBinContent(10);
  }

  TH1D * nEvtSigHist1 = (TH1D*)f1->Get("nEvtSigHisto");
  TH1D * nEvtBkgHist1 = (TH1D*)f1->Get("nEvtSigHisto"); 

  Float_t etaPlotRange = s.getEtaPlotRange();
  double etaranges[8]={1.5,2.5,2.2,10,2.4,10,2.6,10};
  Int_t minbin,maxbin;
  Int_t nCanvas=0;
  for(int e = 0; e<nEnergyBins; e++)
  {
    if (doOneBin&&e!=0) continue;
    for(int m = 0; m<nMultBins; m++)
    {
      //if (doOneBin&&m!=0) continue;
      for(int p = 0; p<nptBins; p++)
      {
        if (doOneBin&&p!=0) continue;
        for(int et = 0; et<netaBins; et++)
        {
          if (doOneBin&&et!=0) continue;
          // preparing variables for the legend and checking if valid axis
          std::string sqrts = "#sqrt{s}";
          if(s.energyBinsHigh[e] <= 100) sqrts = sqrts + "=91GeV";
          else sqrts = sqrts + ">100GeV";
          Float_t etaCut = 0;
          Float_t ptLow = 0;
          Float_t ptHigh = 0;
	  Int_t ana=0;
          if (hMetaData==0&&inFileName1.find("perp1") != std::string::npos ) doPerp=1;
	      if ((hMetaData==0&&inFileName1.find("thrust1") != std::string::npos) || doThrust) {ana=1; etaCut = s.etaBinsLow_wrtThr[et]; ptLow = s.ptBinsLow_wrtThr[p]; ptHigh = s.ptBinsHigh_wrtThr[p];}
          else if((hMetaData==0&&inFileName1.find("wta1") != std::string::npos) || doWTA) {ana=2;etaCut = s.etaBinsLow_wrtWTA[et]; ptLow = s.ptBinsLow_wrtWTA[p]; ptHigh = s.ptBinsHigh_wrtWTA[p];}
          else if((hMetaData==0&&inFileName1.find("wta1") != std::string::npos&&inFileName1.find("perp1") != std::string::npos&&inFileName1.find("thrust1") != std::string::npos) || (doThrust==0&&doWTA==0)) {ana=3;etaCut = s.etaBinsLow_wrtBeam[et]; ptLow = s.ptBinsLow_wrtBeam[p]; ptHigh = s.ptBinsHigh_wrtBeam[p];}
          else {std::cout<<"Unknown axis type...breaking"<<std::endl; break;}
          TLegend * l = new TLegend(0.6,0.8,0.95,0.95);
          l->SetFillStyle(0);
          if(experiment==0)  l->AddEntry((TObject*)0,Form("ALEPH e^{+}e^{-}, %s",sqrts.c_str()),"");
          if(experiment==1)  l->AddEntry((TObject*)0,"DELPHI e^{+}e^{-}","");
          if(experiment==2)  l->AddEntry((TObject*)0,"Belle e^{+}e^{-}","");
          if(experiment==3)  l->AddEntry((TObject*)0,"CMS pp","");
          l->AddEntry((TObject*)0,Form("%d<=N^{trk}<%d",s.multBinsLow[m],s.multBinsHigh[m]),"");
          l->AddEntry((TObject*)0,Form("|#eta|<%1.1f",etaCut),"");
          l->AddEntry((TObject*)0,Form("%1.1f<p_{T}<%1.1f GeV",ptLow,ptHigh),"");
	  if (ana==1) {
	     if (doPerp) { 
	        l->AddEntry((TObject*)0,"Thrust Perp Axis Analysis","");
	     } else {
	        l->AddEntry((TObject*)0,"Thrust Axis Analysis","");
	     }	
	  } else if (ana==2) {
	     if (doPerp) {
   	        l->AddEntry((TObject*)0,"WTA Axis Analysis","");
             } else {
	        l->AddEntry((TObject*)0,"WTA Perp Axis Analysis","");
             }
          } else if (ana==3) {
	     l->AddEntry((TObject*)0,"Beam Axis Analysis","");
          }
          // loading histograms
          signal2PC[e][m][p][et] = (TH2F*)f1->Get(Form("signal2PC_%d_%d_%d_%d_%d",e,s.multBinsLow[m],s.multBinsHigh[m],p,et));
          signal2PC[e][m][p][et]->Scale(1./nEvtSigHist1->GetBinContent(m+1)); // plus 1 because 0 is the underflow bin
          bkgrnd2PC[e][m][p][et] = (TH2F*)f1->Get(Form("bkgrnd2PC_%d_%d_%d_%d_%d",e,s.multBinsLow[m],s.multBinsHigh[m],p,et));
          bkgrnd2PC[e][m][p][et]->Scale(1./nEvtBkgHist1->GetBinContent(m+1));
          ratio2PC[e][m][p][et] = (TH2F*) signal2PC[e][m][p][et]->Clone(Form("ratio2PC_%d_%d_%d_%d_%d",e,s.multBinsLow[m],s.multBinsHigh[m],p,et));
	  ratio2PC[e][m][p][et]->Reset();
          calculateRatio(signal2PC[e][m][p][et],bkgrnd2PC[e][m][p][et],ratio2PC[e][m][p][et]);
          if (signal2PC[e][m][p][et]->GetEntries()==0) continue;
	  
          // For performing 1D projection
          TH1F*h_deltaphi[7];
        
          for (Int_t j=0;j<7;j =j+2)
          {
              // if (i==0) h_deltaphi[m]->SetFillColor(kRed);
              minbin =  ratio2PC[e][m][p][et]->GetXaxis()->FindBin(etaranges[j]);
              maxbin =  ratio2PC[e][m][p][et]->GetXaxis()->FindBin(etaranges[j+1]);
              h_deltaphi[j]  = (TH1F*) ratio2PC[e][m][p][et]->ProjectionY(Form("h_deltaphi%d_%d_%d_%d_%d_%d",j,e,s.multBinsLow[m],s.multBinsHigh[m],p,et),minbin,maxbin);
              h_deltaphi[j]->SetName(Form("h_deltaphi%d_%d_%d_%d_%d_%d",j,e,s.multBinsLow[m],s.multBinsHigh[m],p,et));
              h_deltaphi[j]->GetXaxis()->SetTitle("#Delta#phi");
              if (getThetaAngle)  h_deltaphi[j]->SetTitle(Form("#Delta#phi, #Delta#theta (%1.1f, %1.1f), Multipliplicity (%1.1d, %1.1d)",etaranges[j],etaranges[j+1], s.multBinsLow[m],s.multBinsHigh[m]));
              else            h_deltaphi[j]->SetTitle(Form("#Delta#phi, #Delta#eta (%1.1f, %1.1f), Multipliplicity (%d, %d)",etaranges[j],etaranges[j+1], s.multBinsLow[m],s.multBinsHigh[m]));
              h_deltaphi[j]->GetYaxis()->SetTitle("Y(#Delta#phi)");
              h_deltaphi[j]->Scale(1./(maxbin-minbin+1));
          }

          // drawing plots
	  nCanvas++;
	  
          TCanvas * cAll = new TCanvas(Form("c_%d",nCanvas),Form("Summary%d",nCanvas),800,800);
	  cAll->Divide(2,2);
	  for (Int_t iCanvas = 1; iCanvas<4; iCanvas++) setupCanvas(cAll->GetPad(iCanvas));
          
	  TCanvas * c1 = new TCanvas("c1","c1",800,800);
          setupCanvas(c1);
          formatTPCAxes(signal2PC[e][m][p][et],1.5,1.5,2);
          formatZaxis(signal2PC[e][m][p][et],1);
          c1->cd();
	  signal2PC[e][m][p][et]->Draw("surf1 fb");
          l->Draw("same");
          //sig[m]->GetZaxis()->SetRangeUser(0.9*sig[m]->GetMinimum(),0.4*sig[m]->GetMaximum());
          c1->SaveAs(Form("%s_signal_%d_%d_%d_%d_%d.png",plotsname.c_str(),e,s.multBinsLow[m],s.multBinsHigh[m],p,et));
          c1->SaveAs(Form("%s_signal_%d_%d_%d_%d_%d.pdf",plotsname.c_str(),e,s.multBinsLow[m],s.multBinsHigh[m],p,et));
          
	  cAll->cd(1);
	  signal2PC[e][m][p][et]->Draw("surf1 fb");
          l->Draw("same");
          
	  //c1->SaveAs(Form("%s_signal1_%d_%d_%d_%d_%d.C",plotsname.c_str(),e,s.multBinsLow[m],s.multBinsHigh[m],p,et));
            
          formatTPCAxes(bkgrnd2PC[e][m][p][et],1.5,1.5,2);
          formatZaxis(bkgrnd2PC[e][m][p][et],1);
          c1->cd();
          bkgrnd2PC[e][m][p][et]->Draw("surf1 fb");
          l->Draw("same");
          c1->SaveAs(Form("%s_background_%d_%d_%d_%d_%d.png",plotsname.c_str(),e,s.multBinsLow[m],s.multBinsHigh[m],p,et));
          c1->SaveAs(Form("%s_background_%d_%d_%d_%d_%d.pdf",plotsname.c_str(),e,s.multBinsLow[m],s.multBinsHigh[m],p,et));
          //c1->SaveAs(Form("%s_background1_%d_%d_%d_%d_%d.C",plotsname.c_str(),s.multBinsLow[m],s.multBinsHigh[m]));
          
	  cAll->cd(2);
	  bkgrnd2PC[e][m][p][et]->Draw("surf1 fb");
          l->Draw("same");
          
	    
          formatTPCAxes(ratio2PC[e][m][p][et],1.5,1.5,2);
          formatZaxis(ratio2PC[e][m][p][et],0);
          c1->cd();
          ratio2PC[e][m][p][et]->Draw("surf1 fb");
          l->Draw("same");
          c1->SaveAs(Form("%s_ratio2PC_%d_%d_%d_%d_%d.png",plotsname.c_str(),e,s.multBinsLow[m],s.multBinsHigh[m],p,et));
          c1->SaveAs(Form("%s_ratio2PC_%d_%d_%d_%d_%d.pdf",plotsname.c_str(),e,s.multBinsLow[m],s.multBinsHigh[m],p,et));
          //c1->SaveAs(Form("%s_ratio2PC_%d_%d_%d_%d_%d.C",plotsname.c_str(),s.multBinsLow[m],s.multBinsHigh[m]));
          
	  cAll->cd(3);
	  ratio2PC[e][m][p][et]->Draw("surf1 fb");
	  ratio2PC[e][m][p][et]->Write();
          l->Draw("same");
            
          TCanvas * c2 = new TCanvas("c2","dphi",600,600);
          //c2->Divide(2,2);
          
          // Fourier decomposition
          TLatex* texv2;
          TF1 *f1 = new TF1("f1","[0]*(1+2*([1]*cos(1*x)+[2]*cos(2*x)+[3]*cos(3*x)+[4]*cos(4*x)+[5]*cos(5*x)+[6]*cos(6*x)))");
          for (Int_t j=0;j<4;j++) {
            //c2->cd(j+1);
            //if (i ==2) for(Int_t k = 0; k<h_deltaphi[j*2]->GetNbinsX();k++)std::cout<<h_deltaphi[j*2]->GetBinContent(k)<<std::endl;
            for(Int_t k = 0; k<h_deltaphi[j*2]->GetNbinsX();k++)if(std::isnan(h_deltaphi[j*2]->GetBinError(k+1))) {h_deltaphi[j*2]->SetBinError(k+1,0);}
            formatTH1F(h_deltaphi[j*2],.8,1.5);
            h_deltaphi[j*2]->Draw();
            xjjroot::drawtex(0.15,0.876,Form("%s",plotsname.c_str()));
            h_deltaphi[j*2]->Fit("f1");
            h_deltaphi[j*2]->Fit("f1");
            h_deltaphi[j*2]->SetStats(0);
            /*
            if (j==0){
               h_deltaphi[0]->Fit("f1");
               h_deltaphi[0]->Fit("f1");
               h_deltaphi[0]->SetStats(0);
               texv2 = new TLatex(1, h_deltaphi[0]->GetMaximum(), Form("v_{2}=%.3f #pm %.3f",f1->GetParameter(2),f1->GetParError(2)));
               std::cout<<"ELLIPTIC FLOW"<<f1->GetParameter(1)<<std::endl;
              }
              */
            c2->SaveAs(Form("%s_deltaphi_%d_%d_%d_%d_%d_%d.png",plotsname.c_str(),e,s.multBinsLow[m],s.multBinsHigh[m],p,et,j));
            c2->SaveAs(Form("%s_deltaphi_%d_%d_%d_%d_%d_%d.pdf",plotsname.c_str(),e,s.multBinsLow[m],s.multBinsHigh[m],p,et,j));
            }
            //c2->SaveAs(Form("%s_deltaphi_%d_%d_%d_%d_%d.png",plotsname.c_str(),e,s.multBinsLow[m],s.multBinsHigh[m],p,et));
            //c2->SaveAs(Form("%s_deltaphi_%d_%d_%d_%d_%d.pdf",plotsname.c_str(),e,s.multBinsLow[m],s.multBinsHigh[m],p,et));
            //c2->SaveAs(Form("%s_longRangeYield1_%d_%d.C",plotsname.c_str(),s.multBinsLow[m],s.multBinsHigh[m]));

          cAll->cd(4);
          h_deltaphi[0]->Draw();
          texv2->Draw();
//	  h_deltaphi[0]->Write();
	  delete c1;
          /*
	  delete l;
          delete signal2PC[e][m][p][et];
          delete bkgrnd2PC[e][m][p][et];
          delete ratio2PC[e][m][p][et];
          for (Int_t i=0;i<4;i++)  delete h_deltaphi[i*2];
	  */
	      cAll->SaveAs(Form("%s_ASummary_%d_%d_%d_%d_%d.png",plotsname.c_str(),e,s.multBinsLow[m],s.multBinsHigh[m],p,et));
	      cAll->SaveAs(Form("%s_ASummary_%d_%d_%d_%d_%d.pdf",plotsname.c_str(),e,s.multBinsLow[m],s.multBinsHigh[m],p,et));
            
        }
      }
    }
  }
  
  if (doOneBin) return 0;

  TH1F *eta1, *theta1, *pt1, *phi1, *TTheta1, *TPhi1;
  pt1 = (TH1F*)f1->Get("pt");
  eta1 = (TH1F*)f1->Get("eta");
  theta1 = (TH1F*)f1->Get("theta");
  phi1 = (TH1F*)f1->Get("phi");
  TTheta1 = (TH1F*)f1->Get("T_theta");
  TPhi1 = (TH1F*)f1->Get("T_phi");

  TCanvas * c2 = new TCanvas("c2","c2",800,800);
    
  eta1->Draw("p");
  eta1->SetTitle(";#eta;N");
  eta1->SetStats(0);
  c2->SaveAs(Form("%s_eta1.png",plotsname.c_str()));
  c2->SaveAs(Form("%s_eta1.pdf",plotsname.c_str()));
  //c2->SaveAs(Form("%s_eta1.C",plotsname.c_str()));
    
  phi1->Draw("p");
  phi1->SetTitle(";#phi;N");
  phi1->SetStats(0);
  c2->SaveAs(Form("%s_phi1.png",plotsname.c_str()));
  c2->SaveAs(Form("%s_phi1.pdf",plotsname.c_str()));
  //c2->SaveAs(Form("%s_phi1.C",plotsname.c_str()));
    
  theta1->Draw("p");
  theta1->SetTitle(";#theta;N");
  theta1->SetStats(0);
  c2->SaveAs(Form("%s_theta1.png",plotsname.c_str()));
  c2->SaveAs(Form("%s_theta1.pdf",plotsname.c_str()));
  //c2->SaveAs(Form("%s_theta1.C",plotsname.c_str()));
    
  TTheta1->Draw("p");
  TTheta1->SetTitle(";#theta_{thrust};N");
  TTheta1->SetStats(0);
  c2->SaveAs(Form("%s_TTheta1.png",plotsname.c_str()));
  c2->SaveAs(Form("%s_TTheta1.pdf",plotsname.c_str()));
  //c2->SaveAs(Form("%s_TTheta1.C",plotsname.c_str()));
    
  TPhi1->Draw("p");
  TPhi1->SetTitle(";#phi_{thrust};N");
  TPhi1->SetStats(0);
  c2->SaveAs(Form("%s_TPhi1.png",plotsname.c_str()));
  c2->SaveAs(Form("%s_TPhi1.pdf",plotsname.c_str()));
  //c2->SaveAs(Form("%s_TPhi1.C",plotsname.c_str()));

  c2->SetLogy();
  pt1->Draw("p");
  pt1->SetTitle(";p_{T};N");
  pt1->SetStats(0);
  c2->SaveAs(Form("%s_pt1.png",plotsname.c_str()));
  c2->SaveAs(Form("%s_pt1.pdf",plotsname.c_str()));
  //c2->SaveAs(Form("../pdfDir/%s_pt1.C",plotsname.c_str()));
  c2->SetLogy();
    

  fout->Write();
  delete c2;
  return 0;

}
