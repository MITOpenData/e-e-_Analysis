#include "DataQuality.h"

void CompareMCData(){


  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(0.043);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.145);
  gStyle->SetTitleX(.0f);
  gStyle->SetTitleY(.0f);
  gStyle->SetOptStat(0);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(0.8);
  
  const int samples=2;
  TString namefiles[samples]={"samplevalidationDataThrust.root","samplevalidationMCThrust.root"};
  TFile* files[samples];

  for (int index=0;index<samples;index++){
    files[index]=new TFile(namefiles[index].Data());
  }
   gethists(files[0], files[1], "hist");

  Plotting*myplotLeg=new Plotting();
  TH2F*hempty1D[nvariables1D] [nsteps];


    for(int i=0;i<nvariables1D;i++) {
        for(int j=0;j<nsteps;j++) {
           hempty1D[i][j]=(TH2F*)myplot->GetEmpty(namehistoempty1D[i][j].Data(),xnamelabelhisto1D[i].Data(),ynamelabelhisto1D[i].Data(),xlowerbin1D[i],xupperbin1D[i],ylowerbin1D[i],yupperbin1D[i]);
           canvasVariables[i][j]=new TCanvas(canvasname[i][j].Data(),canvasname[i][j].Data(),500,500);
           hempty1D[i][j]->Draw();
           hVariable1DData[i][j]->SetMaximum(hVariable1DData[i][j]->GetMaximum()*2.);
           hVariable1DData[i][j]->SetLineColor(1);
           hVariable1DMC[i][j]->SetLineColor(2);
           hVariable1DMC[i][j]->SetMarkerColor(2);
           hVariable1DData[i][j]->SetLineWidth(3);
           hVariable1DMC[i][j]->SetLineWidth(2);
           hVariable1DMC[i][j]->SetMarkerColor(2);
           hVariable1DData[i][j]->Draw("same");
           hVariable1DMC[i][j]->Draw("same");
           leg[i][j]=new TLegend(0.4205154,0.7232817,0.8200301,0.8721362);
           leg[i][j]->AddEntry(hVariable1DData[i][j],"DATA","pl");
           leg[i][j]->AddEntry(hVariable1DMC[i][j],"MC","pl");
           leg[i][j]->Draw();
           canvasVariables[i][j]->SaveAs(Form("datamc/%s.pdf",canvasname[i][j].Data()));
        }
     }
  }
  