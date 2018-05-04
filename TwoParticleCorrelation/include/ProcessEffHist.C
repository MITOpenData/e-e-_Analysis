#include "TH3.h"
#include "TH2.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TVector.h"
#include "TStyle.h"
#include "TFile.h"

int ProcessEffHist()
{
	gStyle->SetOptStat(0);

	std::string fullPath = std::getenv("STUDYMULTDIR");
	fullPath.append("/DataProcessing/tables/efficiency_hist.root");
	TFile* f = TFile::Open(fullPath.c_str(),"update");
	TH3F* eff = (TH3F*)f->Get("eff3D");

	TH1F* pteff = (TH1F*)eff->Project3D("xe");
	TH1F* thetaeff = (TH1F*)eff->Project3D("ye");
	TH1F* phieff = (TH1F*)eff->Project3D("ze");
	TH2F* ptthetaeff = (TH2F*)eff->Project3D("xye");
	TH2F* ptphieff = (TH2F*)eff->Project3D("xze");
	TH2F* thetaphieff = (TH2F*)eff->Project3D("yze");

	TH1F* thetaclean = (TH1F*)f->Get("efftheta_clean");
	TH1F* ptclean = (TH1F*)f->Get("effpt_clean");
	TH1F* phiclean = (TH1F*)f->Get("effphi_clean");

	int nptbins = pteff->GetNbinsX();
	int nthetabins = thetaeff->GetNbinsX();
	int nphibins = phieff->GetNbinsX();

	pteff->Scale(1./(nthetabins*nphibins));
	thetaeff->Scale(1./(nptbins*nphibins));
	phieff->Scale(1./(nptbins*nthetabins));
	ptthetaeff->Scale(1./nphibins);
	ptphieff->Scale(1./nthetabins);
	thetaphieff->Scale(1./nptbins);

	f->cd();
	pteff->Write("pteff",TObject::kOverwrite);
	thetaeff->Write("thetaeff",TObject::kOverwrite);
	phieff->Write("phieff",TObject::kOverwrite);
	ptthetaeff->Write("ptthetaeff",TObject::kOverwrite);
	ptphieff->Write("ptphieff",TObject::kOverwrite);
	thetaphieff->Write("thetaphieff",TObject::kOverwrite);

	TCanvas* c = new TCanvas("c","c",800,800);
	pteff->GetXaxis()->SetTitle("p_{t}");
	pteff->GetYaxis()->SetTitle("efficiency");
	pteff->Draw("pe");
	c->SaveAs("pteff.pdf");
	thetaeff->GetXaxis()->SetTitle("#theta");
	thetaeff->GetYaxis()->SetTitle("efficiency");
	thetaeff->Draw("pe");
	c->SaveAs("thetaeff.pdf");
	phieff->GetXaxis()->SetTitle("#phi");
	phieff->GetYaxis()->SetTitle("efficiency");
	phieff->Draw("pe");
	c->SaveAs("phieff.pdf");
	ptthetaeff->SetTitle("(#theta, pt) efficiency");
	ptthetaeff->Draw("COLZ");
	c->SaveAs("ptthetaeff.pdf");
	ptphieff->SetTitle("(#phi, pt) efficiency");
	ptphieff->Draw("COLZ");
	c->SaveAs("ptphieff.pdf");
	thetaphieff->SetTitle("(#phi, #theta) efficiency");
	thetaphieff->Draw("COLZ");
	c->SaveAs("thetaphieff.pdf");
	thetaclean->SetTitle("highPurity, pwflag<=2, pt>1");
	thetaclean->GetXaxis()->SetTitle("#theta");
	thetaclean->GetYaxis()->SetTitle("efficiency");
	thetaclean->Draw("pe");
	c->SaveAs("thetaclean.pdf");
	ptclean->SetTitle("highPurity, pwflag<=2, 0.35<#theta<2.8");
	ptclean->GetXaxis()->SetTitle("p_{t}");
	ptclean->GetYaxis()->SetTitle("efficiency");
	ptclean->Draw("pe");
	c->SaveAs("ptclean.pdf");
	phiclean->SetTitle("highPurity, pwflag<=2, pt>1, 0.35<#theta<2.8");
	phiclean->GetXaxis()->SetTitle("#phi");
	phiclean->GetYaxis()->SetTitle("efficiency");
	phiclean->Draw("pe");
	c->SaveAs("phiclean.pdf");
	f->Close();
	return 0;
}
