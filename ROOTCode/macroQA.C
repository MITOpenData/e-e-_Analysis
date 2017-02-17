#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;


#define PI 3.1415926
enum SIMPLEPID {PHOTON, ELECTRON, PION, MUON, KAON, PROTON};
TString suffix[6]={"photon","electron","pion","muon","kaon","proton"};
int colors[6]={1,1,1,1,1,1};

void macroQA(int maxevt=1000000,int mult=0,int nbin=20){


	TFile *f = new TFile("../Inputs/output-2.root");
	TTree *t1 = (TTree*)f->Get("t");
	Int_t nParticle;
	Float_t pt[10000];
	Float_t eta[10000];
	Float_t pid[10000];
	Float_t phi[10000];
	Float_t mass[10000];


	TH1F*hpt[6];
	TH1F*heta[6];
	for (int i=0;i<6;i++){
		hpt[i] = new TH1F(Form("h%s",suffix[i].Data()),Form("pt %s",suffix[i].Data()),100,-0,30); 
		heta[i] = new TH1F(Form("h%s",suffix[i].Data()),Form("#eta %s",suffix[i].Data()),100,-5,5); 
		hpt[i]->SetName(Form("hpt%s",suffix[i].Data()));
		heta[i]->SetName(Form("heta%s",suffix[i].Data()));
	}  

	t1->SetBranchAddress("nParticle",&nParticle);
	t1->SetBranchAddress("pt",pt);
	t1->SetBranchAddress("eta",eta);
	t1->SetBranchAddress("pid",pid);
	t1->SetBranchAddress("phi",phi);
	t1->SetBranchAddress("mass",mass);

	// all entries and fill the histograms
	Int_t nevent = (Int_t)t1->GetEntries();

	int nevent_process = nevent;
	if( maxevt>0 && maxevt<nevent ) nevent_process = maxevt;  

	for (Int_t i=0;i<nevent_process;i++) {
		if (i%1000==0) cout <<i<<"/"<<nevent_process<<endl;
		t1->GetEntry(i);

		int nparticles = nParticle;
		if (nparticles<30) continue;

		for ( int j=0;j<nparticles;j++ ) {

			int pid1 = pid[j];
			float eta1 = eta[j];
			float phi1 = phi[j];
			float pt1 = pt[j];
			float mass1 = mass[j];
			//if(pt1<0.||pt1>4.) continue;
			for (int j=0;j<6;j++){
				if(pid1==j) {
					hpt[j]->Fill(pt1);      
					heta[j]->Fill(eta1);               
				}      
			}
		}//end of particle loop
	}// end of loop over events  
	TFile* foutput=new TFile("foutputQA.root","recreate");
	foutput->cd();
	for (int j=0;j<6;j++) {
	  hpt[j]->Scale(1./hpt[j]->GetEntries());
	  heta[j]->Scale(1./heta[j]->GetEntries());
	  hpt[j]->Write();
      heta[j]->Write();
    }
}
void plots(){

    TFile* finput=new TFile("foutputQA.root");
	finput->cd();
	TH1F*hpt[6];
	TH1F*heta[6];
	
	for (int j=0;j<6;j++) {
	    hpt[j]=(TH1F*)finput->Get(Form("hpt%s",suffix[j].Data()));
	    heta[j]=(TH1F*)finput->Get(Form("heta%s",suffix[j].Data()));
	}
	TCanvas*canvasPt;
    TH2F* hemptyPt=new TH2F("hemptyPt","",100,0.1,40,10,0.00000001,10.);  
	hemptyPt->GetXaxis()->SetTitle("p_{T}");
	hemptyPt->GetXaxis()->CenterTitle();
	hemptyPt->GetYaxis()->CenterTitle();
	hemptyPt->GetXaxis()->SetTitleOffset(1.);
	hemptyPt->GetYaxis()->SetTitle("dN/dp_{T}(GeV^{-1})");
	hemptyPt->GetYaxis()->SetTitleOffset(1.32);
	hemptyPt->GetXaxis()->SetTitleSize(0.045);
	hemptyPt->GetYaxis()->SetTitleSize(0.045);
	hemptyPt->GetXaxis()->SetTitleFont(42);
	hemptyPt->GetYaxis()->SetTitleFont(42);
	hemptyPt->GetXaxis()->SetLabelFont(42);
	hemptyPt->GetYaxis()->SetLabelFont(42);
	hemptyPt->GetXaxis()->SetLabelSize(0.04);
	hemptyPt->GetYaxis()->SetLabelSize(0.04);  
	hemptyPt->SetMaximum(10.);
	hemptyPt->SetMinimum(0.00001);

	TCanvas*canvasEta;
    TH2F* hemptyEta=new TH2F("hemptyEta","",100,-10.,10.,10,0.00000001,10.);  
	hemptyEta->GetXaxis()->SetTitle("p_{T}");
	hemptyEta->GetXaxis()->CenterTitle();
	hemptyEta->GetYaxis()->CenterTitle();
	hemptyEta->GetXaxis()->SetTitleOffset(1.);
	hemptyEta->GetYaxis()->SetTitle("dN/d#eta");
	hemptyEta->GetYaxis()->SetTitleOffset(1.32);
	hemptyEta->GetXaxis()->SetTitleSize(0.045);
	hemptyEta->GetYaxis()->SetTitleSize(0.045);
	hemptyEta->GetXaxis()->SetTitleFont(42);
	hemptyEta->GetYaxis()->SetTitleFont(42);
	hemptyEta->GetXaxis()->SetLabelFont(42);
	hemptyEta->GetYaxis()->SetLabelFont(42);
	hemptyEta->GetXaxis()->SetLabelSize(0.04);
	hemptyEta->GetYaxis()->SetLabelSize(0.04);  
	hemptyEta->SetMaximum(10.);
	hemptyEta->SetMinimum(0.00001);


	for (int j=0;j<6;j++) {
	
	  canvasPt=new TCanvas("canvasPt","canvasPt",500,500);
	  canvasPt->cd();
	  canvasPt->SetLogy();
	  hemptyPt->GetXaxis()->SetTitle(Form("p_{T} %s",suffix[j].Data()));
	  hemptyPt->Draw();
	  hpt[j]->SetLineColor(colors[j]);
	  hpt[j]->Draw("same");
	  canvasPt->SaveAs(Form("Plots/canvasPt_%s.pdf",suffix[j].Data()));
	  
	  canvasEta=new TCanvas("canvasEta","canvasEta",500,500);
	  canvasEta->cd();
	  canvasEta->SetLogy();
	  hemptyEta->GetXaxis()->SetTitle(Form("#eta %s",suffix[j].Data()));
	  hemptyEta->Draw();
	  heta[j]->SetLineColor(colors[j]);
	  heta[j]->Draw("same");
	  canvasEta->SaveAs(Form("Plots/canvasEta_%s.pdf",suffix[j].Data()));
	}
}
