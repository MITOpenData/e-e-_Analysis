//
//  thrust_distribution.c
//
//
//  Created by Anthony Badea on 10/24/17.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TStyle.h>
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "getLogBins.h"


////////////////////////////////////////////////////////////
//////////////////// Thrust Distribution  //////////////////
// Select events by the following criteria
// Compile/Run with .x thrust_distribution.c+()
// 1. Events accepted with at least 5 good tracks  Done
// 2. Total charged energy > 15GeV  Not done
// 3. Energy (91,91.5) to focus on Z pole
// 3. |cos(TTheta)|<0.9  Done
// 4. missP/Energy<0.3  Done
////////////////////////////////////////////////////////////
void thrust_distribution(TString filename = "/home/abadea/Documents/20171022/alephDataPaths_LEP1.root", // file used
                         Float_t cut_missP = 0.3,   // upper bound on missP/energy
                         Float_t min_TTheta = 0.0, // lower cut on TTheta
                         Float_t max_TTheta = 3.5, // upper cut on TTheta --> currently full range of TTheta
                         Float_t min_Energy = 91, // lower cut on Energy   i.e. Energy>91 to take event
                         Float_t max_Energy = 91.5, // upper cut on Energy    i.e. Energy<91.5 to take event
                         Int_t isCharged = 0, // 0 to include all particle, 1 to only look at charged particles
                         Int_t isGen = 0    // 1 to use gen level
)
{
    TFile *f = new TFile(filename);
    TTree *t1;
    t1 = (TTree*)f->Get("t");
    if(isGen){t1 = (TTree*)f->Get("tgen");}
    
    
    Int_t nParticle;
    Float_t px[5000];
    Float_t py[5000];
    Float_t pz[5000];
    Float_t theta[5000];
    Float_t TTheta;
    Float_t TPhi;
    Int_t pwflag[5000];
    Float_t Energy;
    Float_t missP;
    
    t1->SetBranchAddress("nParticle",&nParticle);
    t1->SetBranchAddress("px",px);
    t1->SetBranchAddress("py",py);
    t1->SetBranchAddress("pz",pz);
    t1->SetBranchAddress("theta",theta);
    t1->SetBranchAddress("TTheta", &TTheta);
    t1->SetBranchAddress("TPhi", &TPhi);
    t1->SetBranchAddress("pwflag",pwflag);
    t1->SetBranchAddress("Energy", &Energy);
    t1->SetBranchAddress("missP", &missP);
    
    Int_t nevent = (Int_t)t1->GetEntries();
    Int_t nevent_process = nevent;
    std::vector<float> T;
    
    const int nBins = 100;
    Double_t bins[nBins+1];
    getLogBins(.05, 1, nBins, bins);
    TH1D *h_thrust = new TH1D("h_thrust",";Thrust;#frac{1}{#sigma} #frac{d#sigma}{dT}",nBins,bins);
    h_thrust->Sumw2();
    for(Int_t i=0;i<nevent_process;i++)
    {
        //cout<<"EVENT "<<i<<endl;
        
        t1->GetEntry(i);
        if (i%10000==0) cout <<i<<"/"<<nevent_process<<endl;

        // NUMBER 2: Cut on Energy
        if(Energy<min_Energy || Energy>max_Energy){continue;}
        // Cut on TTheta
        if(TTheta<min_TTheta || TTheta>max_TTheta)continue;
        // NUMBER 3: Cut of TTheta
        if(fabs(cos(TTheta))>0.9)continue;
        // NUMBER 4: Cut on missP/Energy
        if(missP/Energy>cut_missP)continue;
        
        
        TVector3 v1(1,0,0);
        v1.SetPhi(TPhi);   // keeping rho and theta
        v1.SetTheta(TTheta); // keeping rho and phi
        v1.SetMag(1.0);  // keeping theta and phi
        
        
        Float_t T_sum = 0.0;
        Float_t T_mag = 0.0;
        Int_t pwflag0 = 0;
        
        for (Int_t j=0;j<nParticle;j++)
        {
            //if(isCharged == 1 && pwflag[j]!=0){continue;}
            if(pwflag[j]==0)pwflag0+=1;
            TVector3 v2(px[j],py[j],pz[j]);
            if(fabs(v1.Dot(v2))>v2.Mag())cout<<"DOT = "<<abs(v1.Dot(v2))<<" V2 Magnitude = "<<v2.Mag()<<endl;
            T_sum += fabs(v1.Dot(v2));
            T_mag += v2.Mag();
        }
        
        // NUMBER 1: if we get 5 charged particles
        if(pwflag0>=5)h_thrust->Fill(T_sum/T_mag);
    }

    
   // Double_t scale = 1.0/( h_thrust->GetXaxis()->GetBinWidth(1)*h_thrust->Integral());
    h_thrust->Scale(1.0/h_thrust->Integral());
    for(int i = 0; i < nBins; ++i)
    {
        h_thrust->SetBinContent(i+1, h_thrust->GetBinContent(i+1)/h_thrust->GetBinWidth(i+1));
        h_thrust->SetBinError(i+1, h_thrust->GetBinError(i+1)/h_thrust->GetBinWidth(i+1));
    }
   
    TCanvas *c2 = new TCanvas("c2","",200,10,500,500);
    gStyle->SetOptStat(0);
    h_thrust->Draw();
    c2->SaveAs(Form("tDist_%.2f_%.2f_%d.pdf",min_TTheta,max_TTheta,isCharged));
}


void run()
{
    Float_t min_TTheta[3] = {0.0,0.75,2.4};
    Float_t max_TTheta[3] = {0.75,2.4,3.5};
    for(int i=0;i<2;i++)
    {
        for(int j = 0;j<3;j++)
        {
            thrust_distribution("/home/abadea/Documents/20171022/alephDataPaths_LEP2_1995to2000.root",min_TTheta[j],max_TTheta[j],i);
        }
    }
}


