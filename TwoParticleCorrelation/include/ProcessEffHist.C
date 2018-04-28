#include "TH3.h"
#include "TH2.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TVector.h"

int ProcessEffHist()
{
	gStyle->SetOptStat(0);

	TFile* f = TFile::Open("efficiency_hist.root","update");
	TH3F* eff = (TH3F*)f->Get("eff");

	TH1F* pteff = (TH1F*)eff->Project3D("xe");
	TH1F* thetaeff = (TH1F*)eff->Project3D("ye");
	TH1F* phieff = (TH1F*)eff->Project3D("ze");
	TH2F* ptthetaeff = (TH2F*)eff->Project3D("xye");
	TH2F* ptphieff = (TH2F*)eff->Project3D("xze");
	TH2F* thetaphieff = (TH2F*)eff->Project3D("yze");

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
	pteff->SetTitle("pt efficiency");
	pteff->Draw("pe");
	c->SaveAs("pteff.png");
	thetaeff->SetTitle("#theta efficiency");
	thetaeff->Draw("pe");
	c->SaveAs("thetaeff.png");
	phieff->SetTitle("#phi efficiency");
	phieff->Draw("pe");
	c->SaveAs("phieff.png");
	ptthetaeff->SetTitle("(#theta, pt) efficiency");
	ptthetaeff->Draw("COLZ");
	c->SaveAs("ptthetaeff.png");
	ptphieff->SetTitle("(#phi, pt) efficiency");
	ptphieff->Draw("COLZ");
	c->SaveAs("ptphieff.png");
	thetaphieff->SetTitle("(#phi, #theta) efficiency");
	thetaphieff->Draw("COLZ");
	c->SaveAs("thetaphieff.png");
	return 0;
}