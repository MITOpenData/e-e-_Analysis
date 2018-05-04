#include <iostream>
#include <cstdlib>

#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TH3.h"
#include "TH2.h"
#include "TH1.h"
#include "TVector.h"
#include "getLogBins.h"
#include "getLinBins.h"

int GetNaiveEff()
{
	int ptbins = 30;
	int thetabins = 20;
	int phibins = 20;
	Float_t minpt = 0.2;
	Float_t maxpt = 30;
	int maxmult = 999;
	TFile* f = TFile::Open("/afs/cern.ch/work/m/mipeters/alephMCRecoAfterCutPaths_1994.root");

	Double_t ptlogbins[ptbins+1];
	getLogBins(minpt,maxpt,ptbins,ptlogbins);
	//for(int j=0;j<thetabins;j++) std::cout << ptlogbins[j] << std::endl;
	Double_t thetalinbins[thetabins+1];
	getLinBins(0,3.15,thetabins,thetalinbins);
	//for(int j=0;j<thetabins;j++) std::cout << thetalinbins[j] << std::endl;
	Double_t philinbins[phibins+1];
	getLinBins(-3.15,3.15,phibins,philinbins);
	//for(int j=0;j<thetabins;j++) std::cout << philinbins[j] << std::endl;

	TTree* t = (TTree*)f->Get("t");
	TTree* tgen = (TTree*)f->Get("tgen");
	ULong64_t nevents = t->GetEntriesFast();
	ULong64_t neventsgen = tgen->GetEntriesFast();
	ULong64_t nevt = nevents<neventsgen? nevents:neventsgen;

	TH3F* eff3D = new TH3F("eff3D","eff3D",ptbins,ptlogbins,thetabins,thetalinbins,phibins,philinbins);
	TH3F* effgen3D = new TH3F("effgen3D","effgen3D",ptbins,ptlogbins,thetabins,thetalinbins,phibins,philinbins);
	TH2F* eff = new TH2F("eff","eff",ptbins,ptlogbins,thetabins,thetalinbins);
	TH2F* effgen = new TH2F("effgen","effgen",ptbins,ptlogbins,thetabins,thetalinbins);
	TH1F* efftheta_clean = new TH1F("efftheta_clean","efftheta_clean",thetabins,thetalinbins);
	TH1F* effthetagen_clean = new TH1F("effthetagen_clean","effthetagen_clean",thetabins,thetalinbins);
	TH1F* effpt_clean = new TH1F("effpt_clean","effpt_clean",ptbins,ptlogbins);
	TH1F* effptgen_clean = new TH1F("effptgen_clean","effptgen_clean",ptbins,ptlogbins);
	TH1F* effphi_clean = new TH1F("effphi_clean","effphi_clean",phibins,philinbins);
	TH1F* effphigen_clean = new TH1F("effphigen_clean","effphigen_clean",phibins,philinbins);

	Float_t pt[maxmult];
	Float_t ptgen[maxmult];
	Float_t theta[maxmult];
	Float_t thetagen[maxmult];
	Float_t phi[maxmult];
	Float_t phigen[maxmult];
	Int_t mult = 0;
	Int_t multgen = 0;
	Bool_t highPurity[maxmult];
	Bool_t highPuritygen[maxmult];
	Short_t pwflag[maxmult];
	Short_t pwflaggen[maxmult];
	t->SetBranchAddress("pt",pt);
	tgen->SetBranchAddress("pt",ptgen);
	t->SetBranchAddress("theta",theta);
	tgen->SetBranchAddress("theta",thetagen);
	t->SetBranchAddress("phi",phi);
	tgen->SetBranchAddress("phi",phigen);
	t->SetBranchAddress("nParticle",&mult);
	tgen->SetBranchAddress("nParticle",&multgen);
	t->SetBranchAddress("highPurity",highPurity);
	tgen->SetBranchAddress("highPurity",highPuritygen);
	t->SetBranchAddress("pwflag",pwflag);
	tgen->SetBranchAddress("pwflag",pwflaggen);
	for(ULong64_t i=0;i<nevt;i++)
	{
		if(i % 1000 == 0) std::cout << i << "/" << nevt << std::endl;
		t->GetEntry(i);
		tgen->GetEntry(i);
		for(int j=0;j<mult;j++) 
		{
			if(!highPurity[j] || pwflag[j]>2) continue;
			eff3D->Fill(pt[j],theta[j],phi[j]);
			eff->Fill(pt[j],theta[j]);
			if(pt[j]>1) efftheta_clean->Fill(theta[j]);
			if(theta[j]>0.35 && theta[j]<2.8) effpt_clean->Fill(pt[j]);
			if(pt[j]>1 && theta[j]>0.35 && theta[j]<2.8) effphi_clean->Fill(phi[j]);
		}
		for(int j=0;j<multgen;j++) 
		{
			if(pwflaggen[j]>2) continue;
			effgen3D->Fill(ptgen[j],thetagen[j],phigen[j]);
			effgen->Fill(ptgen[j],thetagen[j]);
			if(ptgen[j]>1) effthetagen_clean->Fill(thetagen[j]);
			if(thetagen[j]>0.35 && thetagen[j]<2.8) effptgen_clean->Fill(ptgen[j]);
			if(ptgen[j]>1 && thetagen[j]>0.35 && thetagen[j]<2.8) effphigen_clean->Fill(phigen[j]);
		}
	}
	eff->Sumw2();
	effgen->Sumw2();
	eff3D->Sumw2();
	effgen3D->Sumw2();
	efftheta_clean->Sumw2();
	effthetagen_clean->Sumw2();
	effpt_clean->Sumw2();
	effptgen_clean->Sumw2();
	effphi_clean->Sumw2();
	effphigen_clean->Sumw2();

	efftheta_clean->Divide(effthetagen_clean);
	effpt_clean->Divide(effptgen_clean);
	effphi_clean->Divide(effphigen_clean);
	eff3D->Divide(effgen3D);
	eff->Divide(effgen);
	//print out code for efficiency table
	/*
	std::cout << "Float_t eff[ptbins][thetabins][phibins] = ";
	std::cout << "{";
	for(int i=0;i<ptbins;i++)
	{
		std::cout << "{";
		for(int j=0;j<thetabins;j++)
		{
			std::cout << "{";
			for(int k=0;k<phibins;k++)
			{
				std::cout << eff->GetBinContent(i+1,j+1,k+1);
				if(k!=phibins-1) std::cout << ", ";
			}
			std::cout << "}";
			if(j!=thetabins-1) std::cout << "," << std::endl;
		}
		std::cout << "}";
		if(i!=ptbins-1) std::cout << "," << std::endl;
	}
	std::cout << "};";
	*/
	std::string fullPath = std::getenv("STUDYMULTDIR");
	std::string outpath = fullPath.append("/DataProcessing/tables/efficiency_hist.root");
	TFile* out = new TFile(outpath.c_str(),"recreate");
	eff->Write("eff");
	eff3D->Write("eff3D");
	efftheta_clean->Write("efftheta_clean");
	effpt_clean->Write("effpt_clean");
	effphi_clean->Write("effphi_clean");
	out->Close();
	f->Close();
	return 0;
}
