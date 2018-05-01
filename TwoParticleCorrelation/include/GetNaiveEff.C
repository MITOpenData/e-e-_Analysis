#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TH3.h"
#include "TVector.h"

int GetNaiveEff()
{
	int ptbins = 50;
	int thetabins = 10;
	int phibins = 10;
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
	Int_t pwflag[maxmult];
	Int_t pwflaggen[maxmult];
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
		for(int j=0;j<mult;j++) {if(!highPurity[j] || pwflag[j]>2) continue; eff->Fill(pt[j],theta[j],phi[j]);}
		for(int j=0;j<multgen;j++) {if(pwflaggen[j]>2) continue; effgen->Fill(ptgen[j],thetagen[j],phigen[j]);}
	}
	eff->Sumw2();
	effgen->Sumw2();
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
	TFile* out = new TFile("efficiency_hist.root","recreate");
	eff->Write("eff");
	TVector* total = new TVector(1);
	total[0] = (int)nevt;
	total->Write("totalevts");
	out->Close();
	f->Close();
	return 0;
}