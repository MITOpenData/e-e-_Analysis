//
//  thrust_distribution.c
//
//
//  Created by Anthony Badea on 10/24/17.
//  Edit by Yen-Jie Lee


//c and c++ dependencies
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

//root dependencies
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLegend.h>
#include <TNtuple.h>

//local headers
#include "../include/getLogBins.h"
#include "../include/doGlobalDebug.h"

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

std::vector<float> compute_errors(std::string dataname);
void relative_error();
int thrust_distribution(TString filename = "~/lxplus/LEP2_all_merged_20180122.root", //file used
                         std::string dataname = "LEP2",
                         Float_t cut_missP = 0.3,   // upper bound on missP/energy
                         Float_t min_TTheta = 0.0, // lower cut on TTheta
                         Float_t max_TTheta = 3.5, // upper cut on TTheta --> currently full range of TTheta
                         Int_t min_nParticle = 0, // lower cut on multiplicity
                         Int_t max_nParticle = 9999 // upper cut on multiplicity
)
{
  //declare some binnings first - just get it out of the way
  const int nBins = 80;
  const int nBinsSys = 8000;
  Double_t bins[nBins+1];
  Double_t binsSys[nBinsSys+1];
  const Double_t logLow = .005;
  const Double_t logHi = .4;
  getLogBins(logLow, logHi, nBins, bins);
  getLogBins(logLow, logHi, nBinsSys, binsSys);

  std::string sysFileName;

  Float_t min_Energy; // lower cut on Energy   i.e. Energy>91 to take event
  Float_t max_Energy; // upper cut on Energy   i.e. Energy<91.5 to take event
  Int_t isGen = 0;    // 1 to use gen level
  Int_t doCorr = 1;
  
  if(dataname == "LEP1") {
     sysFileName = "../inputs/HEPData-ins636645-v1-Table54.root";
     min_Energy = 91;
     max_Energy = 91.5;
  }

  if(dataname == "LEP1MC") {
     sysFileName = "../inputs/HEPData-ins636645-v1-Table54.root";
     min_Energy = 91;
     max_Energy = 91.5;
     doCorr = 0;
  }

  if(dataname == "LEP1MCGen") {
     sysFileName = "../inputs/HEPData-ins636645-v1-Table54.root";
     min_Energy = 91;
     max_Energy = 91.5;
     isGen=1;
     doCorr =0;
  }

  if(dataname == "LEP2") {
     sysFileName = "../inputs/HEPData-ins636645-v1-Table61.root";
     min_Energy = 203;
     max_Energy = 209;
  }
  
  // Published Results
  const int nbin96 = 21;
  Float_t bins96[nbin96+1]={0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.05,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4};
  Float_t bins96i[nbin96+1];
  for (int i=0;i<nbin96+1;i++){
     bins96i[nbin96-i]=1-bins96[i];
  }

    
  TH1F *hppe96 = new TH1F("hppe96","",nbin96,bins96i);
  hppe96->SetBinContent(21,1.017);
  hppe96->SetBinContent(20,6.035);
  hppe96->SetBinContent(19,12.44);
  hppe96->SetBinContent(18,16.07);
  hppe96->SetBinContent(17,16.45);
  hppe96->SetBinContent(16,15.25);
  hppe96->SetBinContent(15,13.38);
  hppe96->SetBinContent(14,11.58);
  hppe96->SetBinContent(13,9.346);
  hppe96->SetBinContent(12,7.159);
  hppe96->SetBinContent(11,5.088);
  hppe96->SetBinContent(10,3.427);
  hppe96->SetBinContent(9,2.482);
  hppe96->SetBinContent(8,1.847);
  hppe96->SetBinContent(7,1.390);
  hppe96->SetBinContent(6,1.072);
  hppe96->SetBinContent(5,0.8465);
  hppe96->SetBinContent(4,0.5661);
  hppe96->SetBinContent(3,0.3065);
  hppe96->SetBinContent(2,0.1248);
  hppe96->SetBinContent(1,0.0184);
  
  hppe96->SetLineColor(kGreen+1);
       
  TFile* outFile_p = new TFile(Form("outFile_%s_%d_%d.root",dataname.c_str(),min_nParticle,max_nParticle), "RECREATE"); // we will write our th1 to a file for debugging purposes
  TNtuple *nt = new TNtuple("nt","","Thrust:GenThrust:T:sumP:nch");
  TNtuple *ntGen = new TNtuple("ntGen","","GenThrust");

  TH1D *h_thrust = new TH1D("h_thrust",";Thrust;#frac{1}{#sigma} #frac{d#sigma}{dT}",42,0.58,1); //moving declarations here for clarity
  h_thrust->Sumw2();
  TH1D *h_one_minus_thrust = new
  TH1D("h_one_minus_thrust",";1-Thrust;#frac{1}{#sigma} #frac{d#sigma}{dT}",20,0.0,0.4);
  h_one_minus_thrust->Sumw2();
  TH1D *h_ratio_one_minus_thrust = (TH1D*)h_one_minus_thrust->Clone();
  h_one_minus_thrust->Sumw2();
  TH1D *h_one_minus_thrust_log = new TH1D("h_one_minus_thrust_log",";1-Thrust log scale;#frac{1}{#sigma} #frac{d#sigma}{dT}",nBins,bins);
  h_one_minus_thrust_log->Sumw2();
  
  TH1D* h_one_minus_thrust_log_SysUp = new TH1D("h_one_minus_thrust_log_SysUp", "", nBins, bins);
  TH1D* h_one_minus_thrust_log_SysDown = new TH1D("h_one_minus_thrust_log_SysDown", "", nBins, bins);

  TH1D* h_one_minus_thrust_log_SysUpGraph = new TH1D("h_one_minus_thrust_log_SysUpGraph", "", nBinsSys, binsSys);
  TH1D* h_one_minus_thrust_log_SysDownGraph = new TH1D("h_one_minus_thrust_log_SysDownGraph", "", nBinsSys, binsSys);

  h_one_minus_thrust_log_SysUpGraph->SetMarkerStyle(20);
  h_one_minus_thrust_log_SysUpGraph->SetMarkerSize(0.01);
  h_one_minus_thrust_log_SysUpGraph->SetMarkerColor(kRed);
  h_one_minus_thrust_log_SysUpGraph->SetLineColor(kRed);

  h_one_minus_thrust_log_SysDownGraph->SetMarkerStyle(20);
  h_one_minus_thrust_log_SysDownGraph->SetMarkerSize(0.01);
  h_one_minus_thrust_log_SysDownGraph->SetMarkerColor(kRed);
  h_one_minus_thrust_log_SysDownGraph->SetLineColor(kRed);

  h_one_minus_thrust_log_SysUp->SetMarkerStyle(20);
  h_one_minus_thrust_log_SysUp->SetMarkerSize(0.01);
  h_one_minus_thrust_log_SysUp->SetMarkerColor(kRed);
  h_one_minus_thrust_log_SysUp->SetLineColor(kRed);

  h_one_minus_thrust_log_SysDown->SetMarkerStyle(20);
  h_one_minus_thrust_log_SysDown->SetMarkerSize(0.01);
  h_one_minus_thrust_log_SysDown->SetMarkerColor(kRed);
  h_one_minus_thrust_log_SysDown->SetLineColor(kRed);

  TH1D *h_thrustNoCorr = new TH1D("h_thrustNoCorr",";ThrustNoCorr;#frac{1}{#sigma} #frac{d#sigma}{dT}",42,0.58,1); //moving declarations here for clarity
  h_thrustNoCorr->Sumw2();
  h_thrustNoCorr->SetMarkerStyle(24);
  h_thrustNoCorr->SetMarkerColor(2);

  TH1D *h_one_minus_thrustNoCorr = new TH1D("h_one_minus_thrustNoCorr",";1-ThrustNoCorr;#frac{1}{#sigma}  #frac{d#sigma}{dT}",20,0.0,0.4);
  h_one_minus_thrustNoCorr->Sumw2();
  h_one_minus_thrustNoCorr->SetMarkerStyle(24);
  h_one_minus_thrustNoCorr->SetMarkerColor(2);
  
  TH1D *h_one_minus_thrustNoCorr_log = new TH1D("h_one_minus_thrustNoCorr_log",";1-ThrustNoCorr log scale;#frac{1}{#sigma} #frac{d#sigma}{dT}",nBins,bins);
  h_one_minus_thrustNoCorr_log->Sumw2();

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  TFile *f = new TFile(filename, "READ"); //specify reading only - dont want to mod input
  if(f->IsZombie())
  {
    std::cout << "BAD FILE" << std::endl;
    return -1;
  }
    TTree *t1;
    TTree *tgen;
    if(isGen) t1 = (TTree*)f->Get("tgen"); //a little strange to do both Get()'s in event isGen -> here only do one
    else t1 = (TTree*)f->Get("t");
            
    Int_t nParticle;
    Float_t px[5000];
    Float_t pmag[5000];
    Float_t mass[5000];
    Float_t py[5000];
    Float_t pz[5000];
    Float_t pt[5000];
    Float_t theta[5000];
    Float_t TTheta;
    Float_t STheta;
    Float_t TPhi;
    Float_t T_Thrust;
    Float_t GenThrust;
    Int_t pwflag[5000];
    Float_t Energy;
    Float_t missP;
    Float_t missTheta;
    bool passesWW;
    Int_t nChargedHadrons;
    
    t1->SetBranchAddress("nParticle",&nParticle);
    t1->SetBranchAddress("px",px);
    t1->SetBranchAddress("py",py);
    t1->SetBranchAddress("pz",pz);
    t1->SetBranchAddress("pt",pt);
    t1->SetBranchAddress("pmag",pmag);
    t1->SetBranchAddress("mass",mass);
    t1->SetBranchAddress("theta",theta);
    t1->SetBranchAddress("Thrust",&T_Thrust);
    t1->SetBranchAddress("nChargedHadrons",&nChargedHadrons);
    t1->SetBranchAddress("TTheta", &TTheta);
    t1->SetBranchAddress("STheta", &STheta);
    t1->SetBranchAddress("TPhi", &TPhi);
    t1->SetBranchAddress("pwflag",pwflag);
    t1->SetBranchAddress("Energy", &Energy);
    t1->SetBranchAddress("missP", &missP);
    t1->SetBranchAddress("missTheta", &missTheta);
    t1->SetBranchAddress("passesWW", &passesWW);

    if(dataname == "LEP1MC") {
       tgen = (TTree*) f->Get("tgen");
       tgen->SetBranchAddress("Thrust",&GenThrust);
    }

    
    Int_t nevent = (Int_t)t1->GetEntries();
    Int_t nevent_process = nevent;

    //if(doGlobalDebug) nevent_process = 1000;
    std::vector<float> T;
    
    
    //for(int i = 0;i<h_one_minus_thrust->GetNbinsX();i++){std::cout<<h_one_minus_thrust->GetBinLowEdge(i)<<std::endl;}
    for(Int_t i=0;i<nevent_process;i++)
    {
        t1->GetEntry(i);
        if(dataname == "LEP1MC") tgen->GetEntry(i);
     
        if (i%100000==0) std::cout <<i<<"/"<<nevent_process<<std::endl;

        if(Energy<min_Energy || Energy>max_Energy){continue;}
        if(fabs(cos(STheta))>0.82)continue;
        if(missP/Energy>cut_missP)continue;
//	if(fabs(cos(missTheta))>0.8) continue;
        // cut on multiplicity
        //if(nParticle < min_nParticle || nParticle > max_nParticle) continue;
        
        TVector3 v1(1,0,0);
        v1.SetPhi(TPhi);   // keeping rho and theta
        v1.SetTheta(TTheta); // keeping rho and phi
        v1.SetMag(1.0);  // keeping theta and phi
        
        Float_t T_sum = 0.0;
        Float_t T_mag = 0.0;
        Int_t pwflag0 = 0;
        Float_t sumP=0;
	Int_t nch=0;
	bool veto=0;
	Int_t nobj=0;
        for (Int_t j=0;j<nParticle;j++)
        {
	    //if (pmag[j]/Energy>0.5) veto=1;  
	    if (pwflag[j]==0&&fabs(cos(theta[j]))<0.94&&pt[j]>0.2){
	       sumP+=fabs(sqrt(pmag[j]*pmag[j]+mass[j]*mass[j]));
	       nch++;
	       nobj++;
               TVector3 v2(px[j],py[j],pz[j]);
               if(fabs(v1.Dot(v2))>v2.Mag())std::cout<<"DOT = "<<abs(v1.Dot(v2))<<" V2 Magnitude = "<<v2.Mag()<<std::endl;
               T_sum += fabs(v1.Dot(v2));
               T_mag += v2.Mag();

	    } else if (pwflag[j]!=0&&fabs(cos(theta[j]))<0.98&&pmag[j]>0.4) {
	       nobj++;
               TVector3 v2(px[j],py[j],pz[j]);
               if(fabs(v1.Dot(v2))>v2.Mag())std::cout<<"DOT = "<<abs(v1.Dot(v2))<<" V2 Magnitude = "<<v2.Mag()<<std::endl;
               T_sum += fabs(v1.Dot(v2));
               T_mag += v2.Mag();
	    }   

            //if(isCharged == 1 && pwflag[j]!=0){continue;}
            if(pwflag[j]==0)pwflag0+=1;
            
        }
        ntGen->Fill(T_Thrust);	

        if (veto) continue;
	if (sumP<15) continue;
	if (nobj<13) continue;
	pwflag0=nChargedHadrons;
        // NUMBER 1: if we get 5 charged particles
	if (nch>=5)
        {
	    Float_t Thrust = T_sum/T_mag;
            h_thrust->Fill(Thrust, 1);
	        h_one_minus_thrust->Fill(1-Thrust, 1); //You have unbinned data! Fill here
	        h_one_minus_thrust_log->Fill(1-Thrust, 1); //You have unbinned data! Fill here

            h_thrustNoCorr->Fill(Thrust);
	        h_one_minus_thrustNoCorr->Fill(1-Thrust); //You have unbinned data! Fill here
	        h_one_minus_thrustNoCorr_log->Fill(1-Thrust); //You have unbinned data! Fill here
	    nt->Fill(Thrust,GenThrust,T_Thrust,sumP,nch);	
        }
    }

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    
    Double_t scale = 1.0/( h_thrust->GetXaxis()->GetBinWidth(1)*h_thrust->Integral());
    h_thrust->Scale(scale);    

    scale = 1.0/( h_thrustNoCorr->GetXaxis()->GetBinWidth(1)*h_thrustNoCorr->Integral());
    h_thrustNoCorr->Scale(scale);    
    
    TCanvas *c2 = new TCanvas("T","",200,10,500,500);
    gStyle->SetOptStat(0);
    h_thrust->GetXaxis()->CenterTitle();
    h_thrust->GetYaxis()->CenterTitle();
    h_thrust->SetTitleOffset(1.5,"Y");
    h_thrust->GetYaxis()->SetRangeUser(0,20);
    h_thrust->SetMarkerStyle(20);
    h_thrust->SetMarkerSize(1);
    
    TFile *corrFileGen;
    TFile *corrFileReco;
    
    if (doCorr){
       corrFileGen = new TFile("outFile_LEP1MCGen_0_9999.root");
       corrFileReco = new TFile("outFile_LEP1MC_0_9999.root");
    
       TH1D *hGen = (TH1D*)corrFileGen->Get("h_thrust");
       hGen->SetName("hGenMC");
       TH1D *hReco = (TH1D*)corrFileReco->Get("h_thrust");
       hGen->SetName("hRecoMC");
       hGen->Divide(hReco);
       TCanvas *c_cor = new TCanvas("c_cor","");
       hGen->Draw();
       h_thrust->Multiply(hGen);
    }    
    c2->cd();
    h_thrust->Draw();
    outFile_p->cd();
    nt->Write();
    ntGen->Write();
    h_thrust->Write();
    //h_thrustNoCorr->Draw("same");
    
    
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    std::vector<float>quad_errors;
    Double_t binLowEdge = 0;
    Double_t binHiEdge = 0;
    Double_t binval = 0;
    // drawing approximate errors from HEP Table 54 data
    if(dataname=="LEP1")
    {
        TFile *hdata = new TFile("../inputs/HEPData-ins636645-v1-Table54.root");
        TH1F *hep;
        hdata->cd("Table 54");
        hep = (TH1F*)gDirectory->Get("Hist1D_y1");
        hep->SetMarkerStyle(5);
	hep->SetLineColor(kBlue);
	hep->SetLineWidth(2);
        hep->Draw("hist SAME");
        hppe96->Draw("hist same");
        TLine* line_thrust = new TLine();
        Float_t sys = 0;
        quad_errors = compute_errors(dataname);
	TLegend *leg = new TLegend(0.3,0.5,0.7,0.9);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->AddEntry(h_thrust,"ALEPH archived data (LEP1)","pl");
	leg->AddEntry(hppe96,"PR294(1998)1","l");
	leg->AddEntry(hep,"EPJC35(2004)457","l");
	leg->Draw();
    }

    c2->SaveAs(Form("mithig_T_%s_%d_%d.pdf",dataname.c_str(),min_nParticle,max_nParticle));
    
    // Fill the 1-T distribution after the correction to T
    //Deleted version that was just rebin of binned data - when possible rebin starting from unbinned!
    Double_t scale_one_minus_T = 1.0/( h_one_minus_thrust->GetXaxis()->GetBinWidth(1)*h_one_minus_thrust->Integral());
    h_one_minus_thrust->Scale(scale_one_minus_T);

    scale_one_minus_T = 1.0/( h_one_minus_thrustNoCorr->GetXaxis()->GetBinWidth(1)*h_one_minus_thrustNoCorr->Integral());
    h_one_minus_thrustNoCorr->Scale(scale_one_minus_T);
    

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    
    TCanvas *c1 = new TCanvas("1-T","",200,10,500,500);
    c1->cd();
    gStyle->SetOptStat(0);
    //gPad->SetLogy();
    //gPad->SetLogx();
    h_one_minus_thrust->GetXaxis()->CenterTitle();
    h_one_minus_thrust->GetYaxis()->CenterTitle();
    h_one_minus_thrust->GetYaxis()->SetRangeUser(0,20);
    h_one_minus_thrust->SetMarkerStyle(4);
    h_one_minus_thrust->Draw();
    h_one_minus_thrustNoCorr->Draw("same");

    //Drawing approximate errors from HEP Table 61 data
    if(dataname=="LEP2")
    {
        TFile *hdata = new TFile("../inputs/HEPData-ins636645-v1-Table61.root");
        TH1F *hep;
        hdata->cd("Table 61");
        hep = (TH1F*)gDirectory->Get("Hist1D_y1");
        hep->SetMarkerStyle(5);
        hep->Draw("hist SAME");

        TLine* line_thrust = new TLine();
        Float_t sys = 0;
        quad_errors = compute_errors(dataname);
        for (int i = 1;i<=h_thrust->GetNbinsX();i++)
        {
            binLowEdge = h_thrust->GetBinLowEdge(i);
            binHiEdge = binLowEdge + h_thrust->GetBinWidth(i);
            binval = h_thrust->GetBinContent(i);
            sys = quad_errors[i+2]; // +2 since thrust distribution starts at 0.6 while HEP errors start at 0.57
	        //std::cout<<"binLow = "<<binLowEdge<<std::endl;
	        //std::cout<<"quad_errors["<<i<<"] = "<<sys<<std::endl;
            line_thrust->DrawLine(binLowEdge, binval - sys, binHiEdge, binval - sys);
            line_thrust->DrawLine(binLowEdge, binval + sys, binHiEdge, binval + sys);
        }
    }

    c1->SaveAs(Form("mithig_1-T_%s_%d_%d.pdf",dataname.c_str(),min_nParticle,max_nParticle));

    TLine* line_one_minus_thrust = new TLine();
    Double_t binLowEdge_T = 0;
    Double_t binHiEdge_T = 0;
    Double_t j;
    float max_error;
    float error_temp;
    Int_t binx;
    return 0; //even if void - useful for searching the control path
}


void run()
{
    Float_t min_TTheta[3] = {0.0,0.75,2.4};
    Float_t max_TTheta[3] = {0.75,2.4,3.5};
    for(int i=0;i<2;i++)
    {
        for(int j = 0;j<3;j++)
        {
            thrust_distribution("/home/abadea/Documents/20171022/alephDataPaths_LEP2_1995to2000.root","LEP1",min_TTheta[j],max_TTheta[j],i);
        }
    }
}


std::vector<float> compute_errors(std::string dataname)
{
  TString hep_file;
  TFile* hdata;
  TH1F* e2 = 0;
  TH1F* e3 = 0;
  if(dataname=="LEP1") 
  {
    hep_file = "../inputs/HEPData-ins636645-v1-Table54.root";
    hdata = new TFile(hep_file);
    hdata->cd("Table 54");
  
    e2 = (TH1F*)gDirectory->Get("Hist1D_y1_e2");
    e3 = (TH1F*)gDirectory->Get("Hist1D_y1_e3");
  }
  if(dataname=="LEP2")
  {
    hep_file = "../inputs/HEPData-ins636645-v1-Table61.root";
    hdata = new TFile(hep_file);
    hdata->cd("Table 61");
  
    e2 = (TH1F*)gDirectory->Get("Hist1D_y1_e1");
    e3 = (TH1F*)gDirectory->Get("Hist1D_y1_e2");
  }
  
  if(e2->GetNbinsX() == e3->GetNbinsX()) std::cout<<"Number of Bins"<<e2->GetNbinsX()<<std::endl;
  else std::cout<<"Different number of bins"<<std::endl;
  
  // compute Quadrature of errors of e2,e3
  std::vector<float> quad_errors;
  
  float error;
  float e2_error;
  float e3_error;
  float sum;
  
  for(int i = 0; i <= e2->GetNbinsX(); i++){
    e2_error = e2->GetBinContent(i+1);
    e3_error = e3->GetBinContent(i+1);
    sum  = pow(e2_error,2)+pow(e3_error,2);
    error = pow(sum,0.5);
    quad_errors.push_back(error);
  }

  return quad_errors;
}

void relative_error()
{
    TFile *outFile = new TFile("outFile.root");
    TH1D *h_one_minus_thrust_log;
    TH1D *h_one_minus_thrust_log_SysDown;
    TH1D *h_one_minus_thrust_log_SysUp;
    h_one_minus_thrust_log = (TH1D*)outFile->Get("h_one_minus_thrust_log");
    h_one_minus_thrust_log_SysDown = (TH1D*)outFile->Get("h_one_minus_thrust_log_SysDown");
    h_one_minus_thrust_log_SysUp = (TH1D*)outFile->Get("h_one_minus_thrust_log_SysUp");

    // Divide
    h_one_minus_thrust_log_SysUp->Divide(h_one_minus_thrust_log);
    h_one_minus_thrust_log_SysDown->Divide(h_one_minus_thrust_log);
    
    // Set to same range
    h_one_minus_thrust_log_SysUp->SetMinimum(0.65);
    h_one_minus_thrust_log_SysUp->SetMaximum(1.35);
    
    h_one_minus_thrust_log_SysDown->SetMinimum(0.65);
    h_one_minus_thrust_log_SysDown->SetMaximum(1.35);
    
    TCanvas *c1 = new TCanvas("log(1-T)","",200,10,500,500);
    gStyle->SetOptStat(0);
    gPad->SetLogx();
    
    h_one_minus_thrust_log_SysUp->GetXaxis()->SetTitle("1-Thrust");
    h_one_minus_thrust_log_SysUp->GetXaxis()->CenterTitle();
    h_one_minus_thrust_log_SysUp->GetYaxis()->SetTitle("Systematic Error");
    h_one_minus_thrust_log_SysUp->GetYaxis()->CenterTitle();
    
    h_one_minus_thrust_log_SysUp->Draw("HIST");
    h_one_minus_thrust_log_SysDown->Draw("HIST SAME");
    
    gPad->RedrawAxis();
    c1->SaveAs("SystematicError.pdf");
    return;
}

