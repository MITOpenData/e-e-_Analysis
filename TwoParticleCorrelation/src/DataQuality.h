#include <iostream>
#include <iomanip>
#include <TGraphAsymmErrors.h> 
#include <TFile.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TF1.h>
#include <TString.h>
#include "include/TPCNtupleData.h"

int const nsteps=5;
TString namestep[nsteps]={"PreES","AfterESQCD","AfterESQCDBarrelCut","AfterESQCDBarrelCutHighPur","AfterESQCDBarrelCutHighPurEtaCut"};

int const nvariables1D=3;
TString namevariable1D[nvariables1D]={"Phi","eta","pt"};
TString namehisto1D[nvariables1D] [nsteps];
TString namelabelhisto1D[nvariables1D]={";#phi;Entries",";#eta;Entries",";p_{T};Entries"};
TH1F*hVariable1D[nvariables1D] [nsteps];
int nbins1D[nvariables1D]={100,100,100};
double lowerbin1D[nvariables1D]={-3.14,-3.14,0.};
double upperbin1D[nvariables1D]={3.14,3.14,10};


int const nvariables2D=1;
TString namevariable2D[nvariables2D]={"Phi_eta"};
TString namehisto2D[nvariables2D] [nsteps];
TString namelabelhisto2D[nvariables2D]={";#phi;#eta"};
TH2F*hVariable2D[nvariables2D] [nsteps];

int nbins2D_x[nvariables2D]={100};
double lowerbin2D_x[nvariables2D]={-3.14};
double upperbin2D_x[nvariables2D]={3.14};
int nbins2D_y[nvariables2D]={100};
double lowerbin2D_y[nvariables2D]={-3.14};
double upperbin2D_y[nvariables2D]={3.14};

int createhists(Option_t* option)
{
  TString opt  = option;
  opt.ToLower();
  if(opt=="savehist")
    {
    std::cout<<"step0"<<std::endl;
    for(int i=0;i<nvariables1D;i++) { 
        for(int j=0;j<nsteps;j++) {
          namehisto1D[i][j]=Form("histo%s%s",namevariable1D[i].Data(),namestep[j].Data());
          hVariable1D[i][j] = new TH1F(namehisto1D[i][j].Data(), namelabelhisto1D[i].Data(), nbins1D[i], lowerbin1D[i], upperbin1D[i]);
          hVariable1D[i][j]->Sumw2();
      }
    }
    
   for(int i=0;i<nvariables2D;i++) { 
        for(int j=0;j<nsteps;j++) {
          namehisto2D[i][j]=Form("histo%s%s",namevariable2D[i].Data(),namestep[j].Data());
          hVariable2D[i][j] = new TH2F(namehisto2D[i][j].Data(), namelabelhisto2D[i].Data(), nbins2D_x[i], lowerbin2D_x[i], upperbin2D_x[i],nbins2D_y[i], lowerbin2D_y[i], upperbin2D_y[i]);
          hVariable2D[i][j]->Sumw2();
      }
    }

   std::cout<<"step1"<<std::endl;
   return 0;
  }
  return 1;
}

int writehists(Option_t* option)
{
  TString opt  = option;
  opt.ToLower();
  if(opt=="savehist")
    {
    for(int i=0;i<nvariables1D;i++) {
        for(int j=0;j<nsteps;j++) {
           hVariable1D[i][j]->Write();
        }
     }
    for(int i=0;i<nvariables2D;i++) {
        for(int j=0;j<nsteps;j++) {
           hVariable2D[i][j]->Write();
        }
     }
      return 0;
    }
    return 1;
  }

void fillHisto(TPCNtupleData data, int indexcandidate, int step=0){

      hVariable1D[0][step]->Fill(data.getPhi(indexcandidate));
      hVariable1D[1][step]->Fill(data.getEta(indexcandidate));
      hVariable1D[2][step]->Fill(data.getPt(indexcandidate));
      hVariable2D[0][step]->Fill(data.getEta(indexcandidate),data.getPhi(indexcandidate));
}