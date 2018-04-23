#include "Plotting.cpp"
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

  Plotting*myplot=new Plotting();

   hVariable1DData[0][0]->Draw();
   hVariable1DMC[0][0]->Draw("same");

   TCanvas*c=new TCanvas("c","c",500,400);
    for(int i=0;i<nvariables2D;i++) {
        for(int j=0;j<nsteps;j++) {
        
        
        }
     }
  }
  
 /*
  
  
  
  
  
  
  
  
  
  
  
  TString namentupleHF[samples]={"ntDkpi","ntKp"};
  TString namentupleHi[samples]={"ntHi","ntHi"};
  TString selections[samples]={"Dgen==23333","Bgen==23333"};
  TString DBlabel[samples]={"D^{0} meson","B^{+} meson"};
  int colours[samples]={1,2};
  int markerstyle[samples]={21,22};
  int width[samples]={2,2};
  TString label[2]={"D pp 5 TeV","B pp 5 TeV"};
  
  Plotting*myplot=new Plotting();

  const int nvariables=3;
  TString string_xaxis[nvariables]={"Decay Length (cm)","Decay Length Error (cm)","p_{T} (GeV)"};
  TString namehempty[nvariables]={"hemptyDecay","hemptyDecayErr","hPt"};
  TString namehisto[samples][nvariables];
  namehisto[0][0]="hDecayD";
  namehisto[0][1]="hDecayErrD";
  namehisto[0][2]="hPtCandidateD";
  namehisto[1][0]="hDecayB";
  namehisto[1][1]="hDecayErrB";
  namehisto[1][2]="hPtCandidateB";
  TString variables[samples][nvariables];
  variables[0][0]="DsvpvDistance";
  variables[0][1]="DsvpvDisErr";
  variables[0][2]="Dpt";
  variables[1][0]="BsvpvDistance";
  variables[1][1]="BsvpvDisErr";
  variables[1][2]="Bpt";
  
  double nbins[nvariables]={2000,2000,2000};
  double lowerrangex[nvariables]={0.001,0.001,1.};
  double upperrangex[nvariables]={5,0.3,100};
  double lowerrangey[nvariables]={0.00001,0.00001,0.00001};
  double upperrangey[nvariables]={1e10,1e10,1e10};

  TFile *files[samples];
  TTree *ntupleHF[samples];
  TTree *ntupleHi[samples];
  TH1F* hVariable[samples][nvariables];
  TH2F* hempty[nvariables];
  TLegend*leg[nvariables];

  for (int index=0;index<samples;index++){
    files[index]=new TFile(namefiles[index].Data());
    ntupleHF[index]=(TTree*)files[index]->Get(namentupleHF[index].Data());
    ntupleHi[index]=(TTree*)files[index]->Get(namentupleHi[index].Data());
    ntupleHF[index]->AddFriend(ntupleHi[index]);
    
   for (int indexv=0;indexv<nvariables;indexv++){
      hVariable[index][indexv]=new TH1F(namehisto[index][indexv].Data(),"",nbins[indexv],lowerrangex[indexv],upperrangex[indexv]);
      ntupleHF[index]->Project(namehisto[index][indexv].Data(),variables[index][indexv].Data(),TCut("pthatweight")*(TCut(selections[index].Data())));
      hVariable[index][indexv]->SetLineColor(colours[index]);
      hVariable[index][indexv]->SetMarkerColor(colours[index]);
      hVariable[index][indexv]->SetLineWidth(width[index]);
    }
  }

  for (int indexv=0;indexv<nvariables;indexv++){
    hempty[indexv]=(TH2F*)myplot->GetEmpty(namehempty[indexv],string_xaxis[indexv].Data(),"Entries",lowerrangex[indexv],upperrangex[indexv],lowerrangey[indexv],upperrangey[indexv]);
    leg[indexv]=(TLegend*)myplot->GetLegend(0.4205154,0.7232817,0.8200301,0.8721362);
     for (int index=0;index<samples;index++){
       leg[indexv]->AddEntry(hVariable[index][indexv],DBlabel[index].Data(),"pl");
     }
  }
  
   TCanvas*c=new TCanvas("c","c",1200,400);
    c->Divide(3,1);
    TFile*foutput=new TFile("output.root","recreate");
    foutput->cd();

    for (int indexv=0;indexv<nvariables;indexv++){
      c->cd(indexv+1);
      gPad->SetLogx();
      gPad->SetLogy();
      hempty[indexv]->Draw();       
       for (int index=0;index<samples;index++){
         hVariable[index][indexv]->Draw("same");
         hVariable[index][indexv]->Write();  
         cout<<"average of variable ="<<variables[index][indexv]<<", is="<<hVariable[index][indexv]->GetMean()<<endl;
       }
       leg[indexv]->Draw();
    }
    c->SaveAs("StudyDBVariables.pdf");
}
*/