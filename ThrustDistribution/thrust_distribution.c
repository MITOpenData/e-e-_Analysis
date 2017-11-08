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
#include <TLine.h>

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

std::vector<float> compute_errors();

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
    
    const int nBins = 50;
    Double_t bins[nBins+1];
    getLogBins(.05, 0.5, nBins, bins);
    TH1D *h_one_minus_thrust = new TH1D("h_one_minus_thrust",";1-Thrust;#frac{1}{#sigma} #frac{d#sigma}{dT}",nBins,bins);
    h_one_minus_thrust->Sumw2();
    TH1D *h_thrust = new TH1D("h_thrust",";Thrust;#frac{1}{#sigma} #frac{d#sigma}{dT}",43,0.57,1);
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
        if(pwflag0>=5)
        {
            Float_t Thrust = T_sum/T_mag;
            h_thrust->Fill(Thrust);
            if((1.0-Thrust)<.05){h_one_minus_thrust->Fill(h_one_minus_thrust->GetBinCenter(1));}
            else{h_one_minus_thrust->Fill(1-Thrust);}
        }
    }

    
    Double_t scale = 1.0/( h_thrust->GetXaxis()->GetBinWidth(1)*h_thrust->Integral());
    h_thrust->Scale(scale);
    h_one_minus_thrust->Scale(1.0/h_one_minus_thrust->Integral());
    
    
    TCanvas *c2 = new TCanvas("T","",200,10,500,500);
    gStyle->SetOptStat(0);
    h_thrust->GetXaxis()->CenterTitle();
    h_thrust->GetYaxis()->CenterTitle();
    h_thrust->SetMarkerStyle(4);
    h_thrust->Draw();
    
    TFile *hdata = new TFile("HEPData-ins636645-v1-Table54.root");
    TH1F *hep;
    hdata->cd("Table 54");
    hep = (TH1F*)gDirectory->Get("Hist1D_y1");
    hep->SetMarkerStyle(5);
    hep->Draw("hist p SAME");
    
    // drawing approximate errors from HEP data
    TLine* line_thrust = new TLine();
    Double_t binLowEdge = 0;
    Double_t binHiEdge = 0;
    Double_t binval = 0;
    Float_t sys = 0;
    std::vector<float> quad_errors = compute_errors();
    for (int i = 0;i<h_thrust->GetNbinsX();i++)
    {
        binLowEdge = h_thrust->GetBinLowEdge(i);
        binHiEdge = binLowEdge + h_thrust->GetBinWidth(i);
        binval = h_thrust->GetBinContent(i);
        sys = quad_errors[i];
        line_thrust->DrawLine(binLowEdge, binval - sys, binHiEdge, binval - sys);
        line_thrust->DrawLine(binLowEdge, binval + sys, binHiEdge, binval + sys);
    }
    
    
    for(int i = 0; i < nBins; ++i)
    {
        Float_t binWidth = h_one_minus_thrust->GetBinWidth(i+1);
        if(i==0)binWidth+=0.05;
        h_one_minus_thrust->SetBinContent(i+1, h_one_minus_thrust->GetBinContent(i+1)/binWidth);
        h_one_minus_thrust->SetBinError(i+1, h_one_minus_thrust->GetBinError(i+1)/binWidth);
    }
    
    TCanvas *c1 = new TCanvas("log(1-T)","",200,10,500,500);
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    gPad->SetLogx();
    h_one_minus_thrust->GetXaxis()->CenterTitle();
    h_one_minus_thrust->GetYaxis()->CenterTitle();
    h_one_minus_thrust->SetMarkerStyle(4);
    h_one_minus_thrust->Draw();
    
    TLine* line_one_minus_thrust = new TLine();
    Double_t binLowEdge_T = 0;
    Double_t binHiEdge_T = 0;
    Double_t j;
    float max_error;
    float error_temp;
    Int_t binx;
    //std::vector<float> one_minus_T_errors;
    
    // ratio of bin value to error
    TH1D *h_ratio_one_minus_thrust = (TH1D*) h_one_minus_thrust->Clone();
    for (int i = 0;i<h_one_minus_thrust->GetNbinsX();i++)
    {
        // get low and high bins for 1-T
        binLowEdge_T = h_one_minus_thrust->GetBinLowEdge(i);
        binHiEdge_T = binLowEdge_T + h_one_minus_thrust->GetBinWidth(i);
        
        // get low and high bins for corresponding T values
        binHiEdge = 1 - binLowEdge_T;
        binLowEdge = 1 - binHiEdge_T;
        j = binLowEdge;
        max_error = 0.0;
        error_temp = 0.0;
        while(j<binHiEdge)
        {
            cout<<"j "<<j<<endl;
            binx = h_thrust->GetXaxis()->FindBin(j);
            error_temp = quad_errors[binx]/h_thrust->GetBinContent(binx) * h_one_minus_thrust->GetBinContent(i);
            cout<<"i "<<i<<" error_temp "<<error_temp<<endl;
            if(error_temp>max_error)max_error = error_temp;
            j+=h_thrust->GetBinWidth(binx);
        }
        cout<<"MAX ERROR"<<endl;
        // Now we have the maximum error value
        //one_minus_T_errors.push_back(max_error);
        binval = h_one_minus_thrust->GetBinContent(i);
        line_one_minus_thrust->DrawLine(binLowEdge_T, binval - max_error, binHiEdge_T, binval - max_error);
        line_one_minus_thrust->DrawLine(binLowEdge_T, binval + max_error, binHiEdge_T, binval + max_error);
        if(max_error>0)
        {
            h_ratio_one_minus_thrust->SetBinContent(i,max_error/h_one_minus_thrust->GetBinContent(i));
            //cout<<"central value = "<<h_one_minus_thrust->GetBinContent(i)<<" error = "<<max_error<<endl;
        }
    }
    
    //c2->SaveAs(Form("tDist_%.2f_%.2f_%d.pdf",min_TTheta,max_TTheta,isCharged));

    TCanvas *c3 = new TCanvas("ratio of central value and error for 1-T","",200,10,500,500);
    gStyle->SetOptStat(0);
    gPad->SetLogx();
    h_ratio_one_minus_thrust->GetYaxis()->SetTitle("Systematic/Central Value");
    h_ratio_one_minus_thrust->GetYaxis()->SetTitleOffset(1.5);
    h_ratio_one_minus_thrust->GetYaxis()->CenterTitle();
    h_ratio_one_minus_thrust->GetXaxis()->CenterTitle();
    h_ratio_one_minus_thrust->Draw();
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


std::vector<float> compute_errors()
{
    TString hep_file = "HEPData-ins636645-v1-Table54.root";
    TFile *hdata = new TFile(hep_file);
    hdata->cd("Table 54");
    // e1 statistical error
    // e2 systematic error 1
    // e3 systematic error 2
    TH1F *e2 = (TH1F*)gDirectory->Get("Hist1D_y1_e2");
    TH1F *e3 = (TH1F*)gDirectory->Get("Hist1D_y1_e3");
    
    int nBins = 0;
    if(e2->GetNbinsX() == e3->GetNbinsX()){nBins = e2->GetNbinsX(); cout<<"Number of Bins"<<e2->GetNbinsX()<<endl;}
    else{cout<<"Different number of bins"<<endl;}
    
    // compute Quadrature of errors of e2,e3
    std::vector<float> quad_errors;
    
    float error;
    float e2_error;
    float e3_error;
    float sum;
    
    for (int i=0;i<=e2->GetNbinsX();i++)
    {
        e2_error = e2->GetBinError(i);
        e3_error = e3->GetBinError(i);
        //cout<<"e2 error "<<e2->GetBinError(i)<<endl;
        //cout<<"e3 error "<<e3->GetBinError(i)<<endl;
        sum  = pow(e2_error,2)+pow(e3_error,2);
        error = pow(sum,0.5);
        //cout<<"quad error "<<error<<endl;
        quad_errors.push_back(error);
    }
    return quad_errors;
}

