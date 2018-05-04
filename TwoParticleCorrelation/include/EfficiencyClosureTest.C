#include <iostream>

#include "TH3.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "getLogBins.h"
#include "getLinBins.h"
#include "../../DataProcessing/include/alephTrkEfficiency.h"
#include "../../DataProcessing/include/xjjrootuti.h"

int EfficiencyClosureTest()
{
	int ptbins = 50;
	int thetabins = 20;
	int phibins = 20;
	Float_t minpt = 0.2;
	Float_t maxpt = 30;
	Int_t maxmult = 999;

	TFile* f = new TFile("/afs/cern.ch/work/m/mipeters/alephMCRecoAfterCutPaths_1994.root","read");
	TTree* t = (TTree*)f->Get("t");
	TTree* tgen = (TTree*)f->Get("tgen");

	alephTrkEfficiency* eff = new alephTrkEfficiency();

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

	Double_t ptlogbins[ptbins+1];
	getLogBins(minpt,maxpt,ptbins,ptlogbins);

	TH1F* ptuncorr = new TH1F("ptuncorr","pt uncorrected",ptbins,ptlogbins);
	TH1F* ptratio = new TH1F("ptratio","pt reco/gen efficiency-corrected",ptbins,ptlogbins);
	TH1F* ptgenfactor = new TH1F("ptgenfactor","ptgenfactor",ptbins,ptlogbins);
	TH1F* thetauncorr = new TH1F("thetauncorr","thetauncorr",thetabins,0,TMath::Pi());
	TH1F* thetaratio = new TH1F("thetaratio","#theta reco/gen efficiency-corrected",thetabins,0,TMath::Pi());
	TH1F* thetagenfactor = new TH1F("thetagenfactor","thetagenfactor",thetabins,0,TMath::Pi());
	TH1F* phiuncorr = new TH1F("phiuncorr","phiuncorr",phibins,-TMath::Pi(),TMath::Pi());
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
			ptuncorr->Fill(pt[j]);
			thetauncorr->Fill(theta[j]);
			phiuncorr->Fill(phi[j]);
			if(eff->efficiency(theta[j],phi[j],pt[j])!=0)
			{
				ptratio->Fill(pt[j],1./eff->efficiency(theta[j],phi[j],pt[j]));
				thetaratio->Fill(theta[j],1./eff->efficiency(theta[j],phi[j],pt[j]));
				phiratio->Fill(phi[j],1./eff->efficiency(theta[j],phi[j],pt[j]));
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
	TH1F* ptreco = (TH1F*)ptratio->Clone();
	ptratio->Divide(ptgenfactor);
	TH1F* thetareco = (TH1F*)thetaratio->Clone();
	thetaratio->Divide(thetagenfactor);
	TH1F* phireco = (TH1F*)phiratio->Clone();
	phiratio->Divide(phigenfactor);

	xjjroot::setgstyle();

	TCanvas* c = new TCanvas("c","c",800,800);
	TPad* top = new TPad("top","top",0,0.25,1.,1.);
	top->SetMargin(xjjroot::margin_pad_left,xjjroot::margin_pad_right,xjjroot::margin_pad_bottom,xjjroot::margin_pad_top);
	TPad* bottom = new TPad("bottom","bottom",0.,0.,1.,0.3);
	bottom->SetMargin(xjjroot::margin_pad_left,xjjroot::margin_pad_right,xjjroot::margin_pad_bottom+.12,xjjroot::margin_pad_top);
	TLegend* lpt = new TLegend(0.7,0.75,0.9,0.9);
	TLegend* ltheta = new TLegend(0.7,0.75,0.9,0.9);
	TLegend* lphi = new TLegend(0.7,0.5,0.9,0.65);
	TH2F* hemptypt_top = new TH2F("pt_top",";;#frac{dN}{dp_{t}}",1,0.,30.,1,0.,ptreco->GetMaximum()*1.1);
	TH2F* hemptypt_ratio = new TH2F("pt_r",";p_{t};corrected/gen",1,0.,30.,1,0.96,1.04);
	TH2F* hemptytheta_top = new TH2F("theta_top",";;#frac{dN}{d#theta}",1,0.,TMath::Pi(),1,0.,thetareco->GetMaximum()*1.1);
	TH2F* hemptytheta_ratio = new TH2F("theta_r",";#theta;corrected/gen",1,0.,TMath::Pi(),1,0.7,1.3);
	TH2F* hemptyphi_top = new TH2F("phi_top",";;#frac{dN}{d#phi}",1,-TMath::Pi(),TMath::Pi(),1,0.,phireco->GetMaximum()*1.1);
	TH2F* hemptyphi_ratio = new TH2F("phi_r",";#phi;corrected/gen",1,-TMath::Pi(),TMath::Pi(),1,0.996,1.004);
	//top->SetLogy();
	top->Draw();
	bottom->Draw();

	top->cd();
	xjjroot::sethempty(hemptypt_top,0,0.15);
	hemptypt_top->Draw();
	xjjroot::setthgrstyle(ptreco,kRed,20,1.,kRed,1,1,-1,-1,-1);
	xjjroot::setthgrstyle(ptgenfactor,kBlack,20,1.,kBlack,1,1,-1,-1,-1);
	xjjroot::setthgrstyle(ptuncorr,kBlue,20,1.,kBlue,1,1,-1,-1,-1);
	ptreco->Draw("pe same");
	ptgenfactor->Draw("pe same");
	ptuncorr->Draw("pe same");
	lpt->AddEntry(ptgenfactor,"gen","lp");
	lpt->AddEntry(ptuncorr,"reco uncorrected","lp");
	lpt->AddEntry(ptreco,"reco corrected","lp");
	xjjroot::setlegndraw(lpt);
	top->Update();

	bottom->cd();
	xjjroot::sethempty(hemptypt_ratio,0,-0.05);
	hemptypt_ratio->GetXaxis()->SetTitleSize(hemptypt_ratio->GetXaxis()->GetTitleSize() / 0.4);
    hemptypt_ratio->GetYaxis()->SetTitleSize(hemptypt_ratio->GetYaxis()->GetTitleSize() / 0.4);
    hemptypt_ratio->GetXaxis()->SetLabelSize(hemptypt_ratio->GetXaxis()->GetLabelSize() / 0.4);
    hemptypt_ratio->GetYaxis()->SetLabelSize(hemptypt_ratio->GetYaxis()->GetLabelSize() / 0.4);
    hemptypt_ratio->GetXaxis()->SetTitleOffset(hemptypt_ratio->GetXaxis()->GetTitleOffset() * 0.9);
    hemptypt_ratio->GetYaxis()->SetTitleOffset(hemptypt_ratio->GetYaxis()->GetTitleOffset() * 0.5);
    hemptypt_ratio->GetXaxis()->SetLabelOffset(hemptypt_ratio->GetXaxis()->GetLabelOffset() * 0.9);
    hemptypt_ratio->GetYaxis()->SetLabelOffset(hemptypt_ratio->GetYaxis()->GetLabelOffset() * 0.5);
	hemptypt_ratio->Draw();
	xjjroot::setthgrstyle(ptratio,kBlack,20,1,kBlack,1,1,-1,-1,-1);
	ptratio->Draw("pe same");
	bottom->Update();
	c->SaveAs("ptratio.pdf");

	top->cd();
	xjjroot::sethempty(hemptytheta_top,0,0.15);
	hemptytheta_top->Draw();
	xjjroot::setthgrstyle(thetareco,kRed,20,1.,kRed,1,1,-1,-1,-1);
	xjjroot::setthgrstyle(thetagenfactor,kBlack,20,1.,kBlack,1,1,-1,-1,-1);
	xjjroot::setthgrstyle(thetauncorr,kBlue,20,1.,kBlue,1,1,-1,-1,-1);
	thetareco->Draw("pe same");
	thetagenfactor->Draw("pe same");
	thetauncorr->Draw("pe same");
	ltheta->AddEntry(thetagenfactor,"gen","lp");
	ltheta->AddEntry(thetauncorr,"reco uncorrected","lp");
	ltheta->AddEntry(thetareco,"reco corrected","lp");
	xjjroot::setlegndraw(ltheta);
	top->Update();

	bottom->cd();
	xjjroot::sethempty(hemptytheta_ratio,0,-0.05);
	hemptytheta_ratio->GetXaxis()->SetTitleSize(hemptytheta_ratio->GetXaxis()->GetTitleSize() / 0.4);
    hemptytheta_ratio->GetYaxis()->SetTitleSize(hemptytheta_ratio->GetYaxis()->GetTitleSize() / 0.4);
    hemptytheta_ratio->GetXaxis()->SetLabelSize(hemptytheta_ratio->GetXaxis()->GetLabelSize() / 0.4);
    hemptytheta_ratio->GetYaxis()->SetLabelSize(hemptytheta_ratio->GetYaxis()->GetLabelSize() / 0.4);
    hemptytheta_ratio->GetXaxis()->SetTitleOffset(hemptytheta_ratio->GetXaxis()->GetTitleOffset() * 0.9);
    hemptytheta_ratio->GetYaxis()->SetTitleOffset(hemptytheta_ratio->GetYaxis()->GetTitleOffset() * 0.5);
    hemptytheta_ratio->GetXaxis()->SetLabelOffset(hemptytheta_ratio->GetXaxis()->GetLabelOffset() * 0.9);
    hemptytheta_ratio->GetYaxis()->SetLabelOffset(hemptytheta_ratio->GetYaxis()->GetLabelOffset() * 0.5);
	hemptytheta_ratio->Draw();
	xjjroot::setthgrstyle(thetaratio,kBlack,20,1,kBlack,1,1,-1,-1,-1);
	thetaratio->Draw("pe same");
	bottom->Update();
	c->SaveAs("thetaratio.pdf");

	top->cd();
	xjjroot::sethempty(hemptyphi_top,0,0.15);
	hemptyphi_top->Draw();
	xjjroot::setthgrstyle(phireco,kRed,20,1.,kRed,1,1,-1,-1,-1);
	xjjroot::setthgrstyle(phigenfactor,kBlack,20,1.,kBlack,1,1,-1,-1,-1);
	xjjroot::setthgrstyle(phiuncorr,kBlue,20,1.,kBlue,1,1,-1,-1,-1);
	phireco->Draw("pe same");
	phigenfactor->Draw("pe same");
	phiuncorr->Draw("pe same");
	lphi->AddEntry(phigenfactor,"gen","lp");
	lphi->AddEntry(phiuncorr,"reco uncorrected","lp");
	lphi->AddEntry(phireco,"reco corrected","lp");
	xjjroot::setlegndraw(lphi);
	top->Update();

	bottom->cd();
	xjjroot::sethempty(hemptyphi_ratio,0,-0.05);
	hemptyphi_ratio->GetXaxis()->SetTitleSize(hemptyphi_ratio->GetXaxis()->GetTitleSize() / 0.4);
    hemptyphi_ratio->GetYaxis()->SetTitleSize(hemptyphi_ratio->GetYaxis()->GetTitleSize() / 0.4);
    hemptyphi_ratio->GetXaxis()->SetLabelSize(hemptyphi_ratio->GetXaxis()->GetLabelSize() / 0.4);
    hemptyphi_ratio->GetYaxis()->SetLabelSize(hemptyphi_ratio->GetYaxis()->GetLabelSize() / 0.4);
    hemptyphi_ratio->GetXaxis()->SetTitleOffset(hemptyphi_ratio->GetXaxis()->GetTitleOffset() * 0.9);
    hemptyphi_ratio->GetYaxis()->SetTitleOffset(hemptyphi_ratio->GetYaxis()->GetTitleOffset() * 0.5);
    hemptyphi_ratio->GetXaxis()->SetLabelOffset(hemptyphi_ratio->GetXaxis()->GetLabelOffset() * 0.9);
    hemptyphi_ratio->GetYaxis()->SetLabelOffset(hemptyphi_ratio->GetYaxis()->GetLabelOffset() * 0.5);
	hemptyphi_ratio->Draw();
	xjjroot::setthgrstyle(phiratio,kBlack,20,1,kBlack,1,1,-1,-1,-1);
	phiratio->Draw("pe same");
	bottom->Update();
	c->SaveAs("phiratio.pdf");

	f->Close();
	delete eff;
	return 0;
}
