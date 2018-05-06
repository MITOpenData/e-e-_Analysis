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
#include "TBox.h"

//local headers
#include "include/Selection.h"
#include "include/utilities.h"
#include "../include/xjjrootuti.h"
#include "../../Utilities/include/plotLogo.h"

void setSys(TH1F *h, float sys=0.03)
{
   for (int i=1;i<=h->GetNbinsX();i++){
      h->SetBinError(i,h->GetBinContent(i)*0.03);
   }
}

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

int compareProjection(const std::string inFileName1, const std::string mcFileName, const std::string outfilename, const std::string plotsname, int doOneBin=0)
{
  // ROOT Global setting
  TH1::SetDefaultSumw2();    TH2::SetDefaultSumw2();

  Selection s;
    
  gStyle->SetLegendBorderSize(0);
  gStyle->SetPadTopMargin(0.08);
      
  
  //gPad->SetTopMargin(0.08);

  /// Initialize the histograms
  std::cout<<"Initializing histograms..."<<std::endl;

  bool getThetaAngle = s.getThetaAngle;
  
  double sys = 0.03;
  


  TFile * f1 = TFile::Open(inFileName1.c_str(),"read");
  TFile * fMC = TFile::Open(mcFileName.c_str(),"read");
  TFile * fout = TFile::Open(outfilename.c_str(),"recreate");
  
  TH1F *hMetaData = (TH1F*) f1->Get("hMetaData");
  
  if (hMetaData!=0){
     cout <<"Read configuration from data file:"<<endl;
     s.experiment = hMetaData->GetBinContent(1);
     s.doThrust   = hMetaData->GetBinContent(2);
     s.doWTA      = hMetaData->GetBinContent(3);
     s.doPerp     = hMetaData->GetBinContent(4);
     s.doGen      = hMetaData->GetBinContent(5);
     /*
     s.getThetaAngle = hMetaData->GetBinContent(6);
     s.nEnergyBins = hMetaData->GetBinContent(7);
     s.nMultBins   = hMetaData->GetBinContent(8);
     s.nptBins     = hMetaData->GetBinContent(9);
     s.netaBins    = hMetaData->GetBinContent(10);
     */
  }

  s.Init();
  Int_t nEnergyBins = s.nEnergyBins;
  Int_t nMultBins = s.nMultBins;
  Int_t nptBins = s.nptBins;
  Int_t netaBins = s.netaBins;
  Int_t experiment = s.experiment;
  bool doThrust = s.doThrust;
  bool doWTA = s.doWTA;
  bool doPerp = s.doPerp;
  bool doGen = s.doGen;  

  int cCount=0;
  TH2F * dphi[nEnergyBins][nMultBins][nptBins][netaBins];
  TH2F * dphiMC[nEnergyBins][nMultBins][nptBins][netaBins];

  Float_t etaPlotRange = s.getEtaPlotRange();
  double etaranges[8]={0,1.6,1.6,3.0,2,3.0};
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
          for (int np = 0; np<3; np++) {
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
             TLegend * l = new TLegend(0.29,0.59,0.73,0.91);
             l->SetFillStyle(0);
             if(experiment==0)  l->AddEntry((TObject*)0,Form("ALEPH e^{+}e^{-}, %s",sqrts.c_str()),"");
             if(experiment==1)  l->AddEntry((TObject*)0,"DELPHI e^{+}e^{-}","");
             if(experiment==2)  l->AddEntry((TObject*)0,"Belle e^{+}e^{-}","");
             if(experiment==3)  l->AddEntry((TObject*)0,"CMS pp","");
             l->AddEntry((TObject*)0,Form("%d#leqN_{Trk}^{Offline}<%d",s.multBinsLow[m],s.multBinsHigh[m]),"");
             l->AddEntry((TObject*)0,Form("|#eta|<%1.1f, %1.1f<p_{T}<%1.1f GeV",etaCut,ptLow,ptHigh),"");
             l->AddEntry((TObject*)0,Form("%1.1f<|#Delta#eta|<%1.1f",etaranges[np*2],etaranges[np*2+1]),"");
	  /*
	  if (ana==1) {
	     if (doPerp) { 
	        l->AddEntry((TObject*)0,"Thrust Perp Axis","");
	     } else {
	        l->AddEntry((TObject*)0,"Thrust Axis","");
	     }	
	  } else if (ana==2) {
	     if (doPerp) {
   	        l->AddEntry((TObject*)0,"WTA Axis","");
             } else {
	        l->AddEntry((TObject*)0,"WTA Perp Axis","");
             }
          } else if (ana==3) {
	     l->AddEntry((TObject*)0,"Beam Axis","");
          }
	  */
          // loading histograms
          dphi[e][m][p][et] = (TH2F*)f1->Get(Form("h_deltaphi%d_%d_%d_%d_%d_%d",np*2,e,s.multBinsLow[m],s.multBinsHigh[m],p,et));
          dphi[e][m][p][et]->SetName(Form("hData_deltaphi%d_%d_%d_%d_%d_%d",np*2,e,s.multBinsLow[m],s.multBinsHigh[m],p,et));
	  dphiMC[e][m][p][et] = (TH2F*)fMC->Get(Form("h_deltaphi%d_%d_%d_%d_%d_%d",np*2,e,s.multBinsLow[m],s.multBinsHigh[m],p,et));
          dphiMC[e][m][p][et]->SetName(Form("hMC_deltaphi%d_%d_%d_%d_%d_%d",np*2,e,s.multBinsLow[m],s.multBinsHigh[m],p,et));
          if (dphi[e][m][p][et]->GetFunction("f1")!=0) {
	     dphi[e][m][p][et]->GetFunction("f1")->Delete();
	  }
	  dphi[e][m][p][et]->SetAxisRange(0,3.14,"x");
	  double minY=dphi[e][m][p][et]->GetMinimum();
	  double maxY=dphi[e][m][p][et]->GetMaximum();
	  if (dphiMC[e][m][p][et]->GetMaximum()>maxY) maxY=dphiMC[e][m][p][et]->GetMaximum();
	  if (dphiMC[e][m][p][et]->GetMinimum()<minY) minY=dphiMC[e][m][p][et]->GetMinimum();
	  dphi[e][m][p][et]->SetAxisRange(minY-0.1*(maxY-minY),maxY+0.1*(maxY-minY),"Y");
	  TCanvas *c = new TCanvas(Form("c%d",cCount),"",1000,1000);
	  gPad->SetTopMargin(0.08);
	  cCount++;
	  dphi[e][m][p][et]->SetMarkerSize(2);
	  dphi[e][m][p][et]->SetLineWidth(3);
	  dphiMC[e][m][p][et]->SetLineWidth(2);
	  dphi[e][m][p][et]->SetTitle("");
	  dphi[e][m][p][et]->SetTitleOffset(1.3,"Y");
	  dphi[e][m][p][et]->SetTitleOffset(1.3,"X");
	  
	  dphi[e][m][p][et]->Draw();
	  TH1F *hSys = (TH1F*)dphi[e][m][p][et]->Clone(Form("hSys_deltaphi%d_%d_%d_%d_%d_%d",np*2,e,s.multBinsLow[m],s.multBinsHigh[m],p,et));
	  setSys(hSys,sys);
	  hSys->SetLineWidth(15);
	  hSys->SetLineColor(kGray+1);
	  hSys->SetMarkerStyle(0);
	  hSys->Draw("same");

	  dphiMC[e][m][p][et]->SetLineColor(4);
	  dphiMC[e][m][p][et]->SetMarkerStyle(0);
	  dphiMC[e][m][p][et]->Draw("hist e c same");
	  dphiMC[e][m][p][et]->Draw("hist c same");
  	  dphi[e][m][p][et]->Draw("same");
	  
	  l->AddEntry(dphi[e][m][p][et],"Archived Data","pl");
	  l->AddEntry(hSys,"Systematical Uncertainty","l");
	  l->AddEntry(dphiMC[e][m][p][et],"Archived PYTHIA 6.1 MC","pl");
	  
	  l->Draw();
 	  plotLogo(1,1,1);
          c->SaveAs(Form("%s_dphiComparison_%d_%d_%d_%d_%d_%d.eps",plotsname.c_str(),e,s.multBinsLow[m],s.multBinsHigh[m],p,et,np));
          c->SaveAs(Form("%s_dphiComparison_%d_%d_%d_%d_%d_%d.png",plotsname.c_str(),e,s.multBinsLow[m],s.multBinsHigh[m],p,et,np));
	  
	  }
        }
      }
    }
  }
  
  return 0;

}
