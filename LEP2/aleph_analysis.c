//
//  aleph_analysis.c
//  
//
//  Created by Anthony Badea and Bibek K Pandit on 5/31/17.
//
//

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <iostream>
#include <stdlib.h>
#include <TMath.h>

using namespace std;


#define PI 3.1415926
#define plank (4.135667516 * pow(10,-18))
enum SIMPLEPWFLAG {CHARGED_TRACK, CHARGED_LEPTONS1, CHARGED_LEPTONS2, V0, PHOTON, NEUTRAL_HADRON};
// LEP2mcggbby200e


void analysis(int isBelle=1, int maxevt=0,int mult=50,int nbin=40,bool verbose=0){
    
    TString filename;
    //if(isBelle) filename="/Users/anthony/Documents/StudyMult/LEP2/ROOTfiles/cleaned_ALEPH_DATA-all.aleph.root";
    if(isBelle) filename="/data/flowex/Datasamples/LEP2_MAIN/ROOTfiles/cleaned_ALEPH_DATA_all.aleph.root";
    
    TFile *f = new TFile(filename.Data());
    TTree *t1 = (TTree*)f->Get("t");
    Int_t nParticle;
    Float_t pt[50000];
    Float_t eta[50000];
    Float_t pid[50000];
    Float_t phi[50000];
    Float_t mass[50000];
    Float_t pwflag[50000];
    Float_t Energy[50000];
    
    t1->SetBranchAddress("nParticle",&nParticle);
    t1->SetBranchAddress("pt",pt);
    t1->SetBranchAddress("eta",eta);
    t1->SetBranchAddress("pid",pid);
    t1->SetBranchAddress("phi",phi);
    t1->SetBranchAddress("mass",mass);
    t1->SetBranchAddress("pwflag",pwflag);
    t1->SetBranchAddress("Energy",Energy);
    
    TCanvas *c1 = new TCanvas("eta_photon", "eta_photon", 600,  600);
    TH1F *h1 = new TH1F("eta_photon", "eta_photon", 100, -PI,PI);
    TCanvas *c2 = new TCanvas("eta_charged", "eta_charged", 600,  600);
    TH1F *h2 = new TH1F("eta_charged", "eta_charged", 100, -PI,PI);
    TCanvas *c3 = new TCanvas("eta_neutral", "eta_neutral", 600,  600);
    TH1F *h3 = new TH1F("eta_neutral", "eta_neutral", 100, -PI,PI);
    
    TCanvas *c4 = new TCanvas("phi_photon", "phi_photon", 600,  600);
    TH1F *h4 = new TH1F("phi_photon", "phi_photon", 100, -PI,PI);
    TCanvas *c5 = new TCanvas("phi_charged", "phi_charged", 600,  600);
    TH1F *h5 = new TH1F("phi_charged", "phi_charged", 100, -PI,PI);
    TCanvas *c6 = new TCanvas("phi_neutral", "phi_neutral", 600,  600);
    TH1F *h6 = new TH1F("phi_neutral", "phi_neutral", 100, -PI,PI);
    
    TCanvas *c7 = new TCanvas("2d_photon", "2d_photon", 900,  800);
    TH2F *h7 = new TH2F("", "Multiplicity of Photons in #eta-#phi", 100, -PI,PI,100, -PI, PI);
    h7->GetXaxis()->SetTitle("#eta");
    h7->GetYaxis()->SetTitle("#phi");
    h7->SetStats(0);
    TCanvas *c8 = new TCanvas("2d_charged", "2d_charged", 900,  800);
    TH2F *h8 = new TH2F("", "Multiplicity of Charged Particles in #eta-#phi", 100, -PI,PI,100, -PI, PI);
    h8->GetXaxis()->SetTitle("#eta");
    h8->GetYaxis()->SetTitle("#phi");
    h8->SetStats(0);
    TCanvas *c9 = new TCanvas("2d_neutral", "2d_neutral", 900,  800);
    TH2F *h9 = new TH2F("", "Multiplicity of Neutral Particles in #eta-#phi", 100, -PI,PI,100, -PI,PI);
    h9->GetXaxis()->SetTitle("#eta");
    h9->GetYaxis()->SetTitle("#phi");
    h9->SetStats(0);
    
    char *h_photon_energy = new char[100];
    char *h_charged_energy = new char[100];
    char *h_neutral_energy = new char[100];
    char *h_all_energy = new char[100];
    char *energy_range = new char[100];
    char *name = new char[100];
    
    TCanvas *c10 = new TCanvas("all_mult", "all_mult", 600,  600);
    TCanvas *c11 = new TCanvas("photon_mult", "photon_mult", 600,  600);
    TCanvas *c12 = new TCanvas("charged_mult", "charged_mult", 600,  600);
    TCanvas *c13 = new TCanvas("neutral_mult", "neutral_mult", 600,  600);
    
    TCanvas *c14 = new TCanvas("total_mult all energies","total_mult all energies",600,600);
    TCanvas *c15 = new TCanvas("average multiplicity","average multiplicity",600,600);
    
    // all entries and fill the histograms
    Int_t nevent = (Int_t)t1->GetEntries();

    if (maxevt > 0 && maxevt < nevent) nevent_process = maxevt;
    int nevent_process = nevent;
    
    for (int j = 130; j<210; j = j + 10)
    {
      sprintf(energy_range,"%d<=Energy && Energy <= %d",j,j+10);

      c10->cd();
      t1->Draw("nParticle>>h_all",energy_range);
      TH1F *h_all = (TH1F*)gDirectory->Get("h_all");
      h_all->SetName("");
      h_all->SetTitle(Form("Multiplicity for All Particles Between Energies %d and %d GeV",j,j+10));
      h_all->GetXaxis()->SetTitle("Multiplicity");
      h_all->GetYaxis()->SetTitle("Number of Events");
      h_all->GetYaxis()->SetTitleOffset(1.5);
      h_all->Sumw2();
      h_all->SetMarkerStyle(20);
      sprintf(h_all_energy, "h_all_%d.pdf", j);
      c10->SaveAs(Form("Images/%s",h_all_energy));
    
      c11->cd();
      t1->Draw("Sum$(pwflag==4)>>h_photon",energy_range);
      TH1F *h_photon = (TH1F*)gDirectory->Get("h_photon");
      h_photon->SetName("");
      h_photon->SetTitle(Form("Multiplicity for Photons Between Energies %d and %d GeV",j,j+10));
      h_photon->GetXaxis()->SetTitle("Multiplicity");
      h_photon->GetYaxis()->SetTitle("Number of Events");
      h_photon->GetYaxis()->SetTitleOffset(1.5);
      h_photon->Sumw2();
      h_photon->SetMarkerStyle(20);
      sprintf(h_photon_energy, "h_photon_%d.pdf", j);
      c11->SaveAs(Form("Images/%s",h_photon_energy));
        
      c12->cd();
      t1->Draw("Sum$(pwflag==0 || pwflag==1 || pwflag == 2)>>h_charged",energy_range);
      TH1F *h_charged = (TH1F*)gDirectory->Get("h_charged");
      h_charged->SetName("");
      h_charged->SetTitle(Form("Multiplicity for Charged Particles Between Energies %d and %d GeV",j,j+10));
      h_charged->GetXaxis()->SetTitle("Multiplicity");
      h_charged->GetYaxis()->SetTitle("Number of Events");
      h_charged->GetYaxis()->SetTitleOffset(1.5);
      h_charged->Sumw2();
      h_charged->SetMarkerStyle(20);
      sprintf(h_charged_energy, "h_charged_%d.pdf", j);
      c12->SaveAs(Form("Images/%s",h_charged_energy));
        
      c13->cd();
      t1->Draw("Sum$(pwflag==5)>>h_neutral",energy_range);
      TH1F *h_neutral = (TH1F*)gDirectory->Get("h_neutral");
      h_neutral->SetName("");
      h_neutral->SetTitle(Form("Multiplicity for Neutral Particles Between Energies %d and %d GeV", j,j+10));
      h_neutral->GetXaxis()->SetTitle("Multiplicity");
      h_neutral->GetYaxis()->SetTitle("Number of Events");
      h_neutral->GetYaxis()->SetTitleOffset(1.5);
      h_neutral->Sumw2();
      h_neutral->SetMarkerStyle(20);
      sprintf(h_neutral_energy, "h_neutral_%d.pdf", j);
      c13->SaveAs(Form("Images/%s",h_neutral_energy));
        
    }
    
    for (Int_t i = 0; i<nevent_process; i++)
    {
        t1->GetEntry(i);
        int nparticles = nParticle;
        for (Int_t j =0;j<nparticles;j++)
        {
            if(pwflag[j]==PHOTON)
            {
                cout<<"hi"<<endl;
                h1->Fill(eta[j]);
                h4->Fill(phi[j]);
                h7->Fill(eta[j],phi[j]);
            }
            
            if(pwflag[j] == CHARGED_TRACK || pwflag[j] == CHARGED_LEPTONS1 || pwflag[j] == CHARGED_LEPTONS2)
            {
                h2->Fill(eta[j]);
                h5->Fill(phi[j]);
                h8->Fill(eta[j],phi[j]);
            }
            
            if(pwflag[j] == NEUTRAL_HADRON) //pwflag[j] == V0 ||
            {
                h3->Fill(eta[j]);
                h6->Fill(phi[j]);
                h9->Fill(eta[j],phi[j]);
            }
        }
        
    }
    
    c1->cd();
    h1->Draw();
    c1->SaveAs(Form("Images/eta_photon.pdf"));
    c2->cd();
    h2->Draw();
    c2->SaveAs(Form("Images/eta_charged.pdf"));
    c3->cd();
    h3->Draw();
    c3->SaveAs(Form("Images/eta_neutral.pdf"));
    c4->cd();
    h4->Draw();
    c4->SaveAs(Form("Images/phi_photon.pdf"));
    c5->cd();
    h5->Draw();
    c5->SaveAs(Form("Images/phi_charged.pdf"));
    c6->cd();
    h6->Draw();
    c6->SaveAs(Form("Images/phi_neutral.pdf"));
    c7->cd();
    h7->Draw("colz");
    c7->SaveAs(Form("Images/2d_photon.pdf"));
    c8->cd();
    h8->Draw("colz");
    c8->SaveAs(Form("Images/2d_charged.pdf"));
    c9->cd();
    h9->Draw("colz");
    c9->SaveAs(Form("Images/2d_neutral.pdf"));
    
    c14->cd();
    t1->Draw("nParticle>>h_14");
    TH1F *h_14 = (TH1F*)gDirectory->Get("h_14");
    h_14->SetName("");
    h_14->SetTitle("Number of Events vs Multiplicity for All Energies");
    h_14->GetXaxis()->SetTitle("Multiplicity");
    h_14->GetYaxis()->SetTitle("Number of Events");
    h_14->GetYaxis()->SetTitleOffset(1.55);
    c14->SaveAs(Form("Images/mult_all_energies.pdf"));
    
    c15->cd();
    t1->Draw("nParticle:Energy>>h_15(8,130,210)","","prof");
    TH1F *h_15 = (TH1F*)gDirectory->Get("h_15");
    h_15->SetName("");
    h_15->SetTitle("Average Multiplicity vs Energy");
    h_15->GetXaxis()->SetTitle("Energy (GeV)");
    h_15->GetYaxis()->SetTitle("<Multiplicity>");
    h_15->GetYaxis()->SetTitleOffset(1.2);
    c15->SaveAs(Form("Images/mean_mult_all_energies.pdf"));
}

