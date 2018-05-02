#include <iostream>

#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TH3.h"
#include "TVector.h"

int GetNaiveEff()
{
	int ptbins = 50;
	int thetabins = 20;
	int phibins = 20;
	Float_t maxpt = 30;
	int maxmult = 999;
	TFile* f = TFile::Open("/afs/cern.ch/work/m/mipeters/alephMCRecoAfterCutPaths_1994.root");
	TTree* t = (TTree*)f->Get("t");
	TTree* tgen = (TTree*)f->Get("tgen");
	ULong64_t nevents = t->GetEntriesFast();
	ULong64_t neventsgen = tgen->GetEntriesFast();
	ULong64_t nevt = nevents<neventsgen?nevents:neventsgen;
	TH3F* eff = new TH3F("eff","eff",ptbins,0,maxpt,thetabins,0,TMath::Pi(),phibins,-TMath::Pi(),TMath::Pi());
	TH3F* effgen = new TH3F("effgen","effgen",ptbins,0,maxpt,thetabins,0,TMath::Pi(),phibins,-TMath::Pi(),TMath::Pi());
	TH1F* efftheta_clean = new TH1F("efftheta_clean","efftheta_clean",thetabins,0,TMath::Pi());
	TH1F* effthetagen_clean = new TH1F("effthetagen_clean","effthetagen_clean",thetabins,0,TMath::Pi());
	TH1F* effpt_clean = new TH1F("effpt_clean","effpt_clean",ptbins,0,maxpt);
	TH1F* effptgen_clean = new TH1F("effptgen_clean","effptgen_clean",ptbins,0,maxpt);
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
			eff->Fill(pt[j],theta[j],phi[j]);
			if(pt[j]>1) efftheta_clean->Fill(theta[j]);
			if(theta[j]>0.35 && theta[j]<2.8) effpt_clean->Fill(pt[j]);
		}
		for(int j=0;j<multgen;j++) 
		{
			if(pwflaggen[j]>2) continue;
			effgen->Fill(ptgen[j],thetagen[j],phigen[j]);
			if(ptgen[j]>1) effthetagen_clean->Fill(thetagen[j]);
			if(thetagen[j]>0.35 && thetagen[j]<2.8) effptgen_clean->Fill(ptgen[j]);
		}
	}
	eff->Sumw2();
	effgen->Sumw2();
	efftheta_clean->Sumw2();
	effthetagen_clean->Sumw2();
	effpt_clean->Sumw2();
	effptgen_clean->Sumw2();
	efftheta_clean->Divide(effthetagen_clean);
	effpt_clean->Divide(effptgen_clean);
	eff->Divide(effgen);
	//print out code for efficiency table
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
	TFile* out = new TFile("tables/efficiency_hist.root","recreate");
	eff->Write("eff");
	efftheta_clean->Write("efftheta_clean");
	effpt_clean->Write("effpt_clean");
	out->Close();
	f->Close();
	return 0;
}
