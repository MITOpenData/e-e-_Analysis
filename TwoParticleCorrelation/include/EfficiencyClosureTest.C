#include <iostream>

#include "TH3.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"

int EfficiencyClosureTest()
{
	int ptbins = 50;
	int thetabins = 20;
	int phibins = 20;
	Float_t maxpt = 30;
	Int_t maxmult = 999;

	TFile* f = new TFile("/afs/cern.ch/work/m/mipeters/alephMCRecoAfterCutPaths_1994.root","read");
	TTree* t = (TTree*)f->Get("t");
	TTree* tgen = (TTree*)f->Get("tgen");

	TFile* e = new TFile("DataProcessing/tables/efficiency_hist.root","read");
	TH3F* eff = (TH3F*)e->Get("eff");

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

	TH1F* ptratio = new TH1F("ptratio","pt reco/gen efficiency-corrected",ptbins,0,maxpt);
	TH1F* ptgenfactor = new TH1F("ptgenfactor","ptgenfactor",ptbins,0,maxpt);
	TH1F* thetaratio = new TH1F("thetaratio","#theta reco/gen efficiency-corrected",thetabins,0,TMath::Pi());
	TH1F* thetagenfactor = new TH1F("thetagenfactor","thetagenfactor",thetabins,0,TMath::Pi());
	TH1F* phiratio = new TH1F("phiratio","#phi reco/gen efficiency-corrected",phibins,-TMath::Pi(),TMath::Pi());
	TH1F* phigenfactor = new TH1F("phigenfactor","ptgenfactor",phibins,-TMath::Pi(),TMath::Pi());

	ULong64_t nevents = t->GetEntriesFast();
	ULong64_t neventsgen = tgen->GetEntriesFast();
	ULong64_t nevt = nevents<neventsgen?nevents:neventsgen;

	for(ULong64_t i=0;i<nevt;i++)
	{
		if(i % 1000 == 0) std::cout << i << "/" << nevt << std::endl;
		t->GetEntry(i);
		tgen->GetEntry(i);
		for(int j=0;j<mult;j++)
		{
			if(!highPurity[j] || theta[j]<0.35 || theta[j]>2.8 || pt[j]<1 || pwflag[j]>2) continue;
			if(eff->GetBinContent(eff->FindBin(pt[j],theta[j],phi[j]))!=0)
			{
				ptratio->Fill(pt[j],1./(eff->GetBinContent(eff->FindBin(pt[j],theta[j],phi[j]))));
				thetaratio->Fill(theta[j],1./(eff->GetBinContent(eff->FindBin(pt[j],theta[j],phi[j]))));
				phiratio->Fill(phi[j],1./(eff->GetBinContent(eff->FindBin(pt[j],theta[j],phi[j]))));
			}
		}
		for(int j=0;j<multgen;j++)
		{
			if(thetagen[j]<0.35 || thetagen[j]>2.8 || ptgen[j]<1 || pwflaggen[j]>2) continue;
			ptgenfactor->Fill(ptgen[j]);
			thetagenfactor->Fill(thetagen[j]);
			phigenfactor->Fill(phigen[j]);
		}
	}
	ptratio->Sumw2();
	thetaratio->Sumw2();
	phiratio->Sumw2();
	ptgenfactor->Sumw2();
	thetagenfactor->Sumw2();
	phigenfactor->Sumw2();
	ptratio->Divide(ptgenfactor);
	thetaratio->Divide(thetagenfactor);
	phiratio->Divide(phigenfactor);

	gStyle->SetOptStat(0);

	TCanvas* c = new TCanvas("c","c",800,800);
	ptratio->Draw("pe");
	c->SaveAs("ptratio.png");
	thetaratio->Draw("pe");
	c->SaveAs("thetaratio.png");
	phiratio->Draw("pe");
	c->SaveAs("phiratio.png");
	
	f->Close();
	return 0;
}
