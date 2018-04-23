#include <iostream>
#include <iomanip>
#include <TGraphAsymmErrors.h> 
#include <TFile.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TF1.h>
#include <TString.h>
#include "Plotting.h"
#include "include/TPCNtupleData.h"

int const nsteps=5;
TString namestep[nsteps]={"PreES","AfterESQCD","AfterESQCDBarrelCut","AfterESQCDBarrelCutHighPur","AfterESQCDBarrelCutHighPurEtaCut"};

int const nvariables1D=3;
TString namevariable1D[nvariables1D]={"Phi","eta","pt"};
TString namehisto1D[nvariables1D] [nsteps];
TString namelabelhisto1D[nvariables1D]={";#phi;Entries",";#eta;Entries",";p_{T};Entries"};
TString xnamelabelhisto1D[nvariables1D]={"#phi","#eta","p_{T}"};
TString ynamelabelhisto1D[nvariables1D]={"Entries","Entries","Entries"};

TH1F*hVariable1D[nvariables1D] [nsteps];
int nbins1D[nvariables1D]={100,100,100};
double xlowerbin1D[nvariables1D]={-3.14,-3.14,0.};
double xupperbin1D[nvariables1D]={3.14,3.14,10};
double ylowerbin1D[nvariables1D]={0., 0.,0.};
double yupperbin1D[nvariables1D]={10000, 10000,10000};



// CHANGE UP TO HERE ONLY.



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

TH1F*hVariable1DData[nvariables1D] [nsteps];
TH1F*hVariable1DMC[nvariables1D] [nsteps];
TH2F*hVariable2DData[nvariables2D] [nsteps];
TH2F*hVariable2DMC[nvariables2D] [nsteps];


TH2F*hempty1D[nvariables1D] [nsteps];
TString namehistoempty1D[nvariables1D][nsteps];
TCanvas*canvasVariables[nvariables1D] [nsteps];
TLegend*leg[nvariables1D] [nsteps];
TString canvasname[nvariables1D][nsteps];


int createhists(Option_t* option)
{
  TString opt  = option;
  opt.ToLower();
  if(opt=="savehist")
    {
    for(int i=0;i<nvariables1D;i++) { 
        for(int j=0;j<nsteps;j++) {
          namehisto1D[i][j]=Form("histo%s%s",namevariable1D[i].Data(),namestep[j].Data());
          hVariable1D[i][j] = new TH1F(namehisto1D[i][j].Data(), namelabelhisto1D[i].Data(), nbins1D[i], xlowerbin1D[i], xupperbin1D[i]);
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

int gethists(TFile* infData, TFile* infMC, Option_t* option)
{
  Plotting *myplot=new Plotting();
  TString opt  = option;
  opt.ToLower();
  if(opt=="hist"){
    for(int i=0;i<nvariables1D;i++) {
        for(int j=0;j<nsteps;j++) {
           hVariable1DData[i][j] = (TH1F*)infData->Get(Form("histo%s%s",namevariable1D[i].Data(),namestep[j].Data()));
           hVariable1DMC[i][j] = (TH1F*)infMC->Get(Form("histo%s%s",namevariable1D[i].Data(),namestep[j].Data()));
           namehistoempty1D[i][j]=Form("histoempty%s%s",namevariable1D[i].Data(),namestep[j].Data());
           hempty1D[i][j]=(TH2F*)myplot->GetEmpty(namehistoempty1D[i][j].Data(),xnamelabelhisto1D[i].Data(),ynamelabelhisto1D[i].Data(),xlowerbin1D[i],xupperbin1D[i],ylowerbin1D[i],yupperbin1D[i]);
           canvasname[i][j]=Form("canvas%s%s",namevariable1D[i].Data(),namestep[j].Data());
         } 
      }
    for(int i=0;i<nvariables2D;i++) {
        for(int j=0;j<nsteps;j++) {
           hVariable2DData[i][j] = (TH2F*)infData->Get(Form("histo%s%s",namevariable2D[i].Data(),namestep[j].Data()));
           hVariable2DMC[i][j] = (TH2F*)infMC->Get(Form("histo%s%s",namevariable2D[i].Data(),namestep[j].Data()));
         } 
      }
    return 0;
    }
  std::cout<<"error: invalid option for gethists()"<<std::endl;
  return 1;
}



void fillHisto(TPCNtupleData data, int indexcandidate, int step=0){
      hVariable1D[0][step]->Fill(data.getPhi(indexcandidate));
      hVariable1D[1][step]->Fill(data.getEta(indexcandidate));
      hVariable1D[2][step]->Fill(data.getPt(indexcandidate));
      
           
      hVariable2D[0][step]->Fill(data.getEta(indexcandidate),data.getPhi(indexcandidate));
}