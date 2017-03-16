\
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>

using namespace std;


#define PI 3.1415926
enum SIMPLEPID {PHOTON, ELECTRON, PION, MUON, KAON, PROTON};

void macroQA(int maxevt=0,int mult=0,int nbin=20){

  TFile *f = new TFile("output-2.root");
  TTree *t1 = (TTree*)f->Get("t");
  Int_t nParticle;
  Float_t pt[10000];
  Float_t eta[10000];
  Float_t pid[10000];
  Float_t phi[10000];
  Float_t mass[10000];
  Float_t theta[10000];
  
  t1->SetBranchAddress("nParticle",&nParticle);
  t1->SetBranchAddress("pt",pt);
  t1->SetBranchAddress("eta",eta);
  t1->SetBranchAddress("pid",pid);
  t1->SetBranchAddress("phi",phi);
  t1->SetBranchAddress("mass",mass);
  t1->SetBranchAddress("theta",theta);

  // all entries and fill the histograms
  Int_t nevent = (Int_t)t1->GetEntries();
  
  int nevent_process = nevent;
  if( maxevt>0 && maxevt<nevent ) nevent_process = maxevt;  
  
  //

  //Pion ETA
  TCanvas *c1 = new TCanvas("Pion Eta","Pion Eta",600,600);
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  TH1F *pion_eta = new TH1F ("pion_eta","pion eta distribution",100,-20.,20.);
  pion_eta->SetStats(0);
  pion_eta->Draw();
  pion_eta->GetXaxis()->SetTitle("eta");
  pion_eta->GetYaxis()->SetTitle("Number of Events");
  pion_eta->GetYaxis()->SetTitleOffset(1.4);


  // Proton ETA
  TCanvas *c8 = new TCanvas("Proton Eta","Proton Eta",600,600);
  c8->cd();
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  TH1F *proton_eta = new TH1F ("proton_eta","proton eta distribution",100,-20.,20.);
  proton_eta->SetStats(0);
  proton_eta->Draw();
  proton_eta->GetXaxis()->SetTitle("eta");
  proton_eta->GetYaxis()->SetTitle("Number of Events");
  proton_eta->GetYaxis()->SetTitleOffset(1.4);

  // Kaon ETA
  TCanvas *c9 = new TCanvas("Kaon Eta","Kaon Eta",600,600);
  c9->cd();
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  TH1F *kaon_eta = new TH1F ("kaon_eta","kaon eta distribution",100,-20.,20.);
  kaon_eta->SetStats(0);
  kaon_eta->Draw();
  kaon_eta->GetXaxis()->SetTitle("eta");
  kaon_eta->GetYaxis()->SetTitle("Number of Events");
  kaon_eta->GetYaxis()->SetTitleOffset(1.4);
    
  // PT
  TCanvas *c2 = new TCanvas("PT","PT",600,600);
  c2->cd();
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  TH1F *pion_pt = new TH1F ("pion_pt","pion",100,0.,20.);
  TH1F *proton_pt = new TH1F ("proton_pt","proton",100,0.,20.);
  TH1F *kaon_pt = new TH1F ("proton_pt","pion",100,0.,20.);
  pion_pt->SetStats(0);
  pion_pt->Draw();
  pion_pt->GetXaxis()->SetTitle("Transverse Momentum");
  pion_pt->GetYaxis()->SetTitle("Number of Events");
  pion_pt->GetYaxis()->SetTitleOffset(1.4);
  proton_pt->SetLineColor(kRed);
  proton_pt->Draw("same");
  kaon_pt->SetLineColor(kGreen);
  kaon_pt->Draw("same");
    
  // Pion Multiplicity
  TCanvas *c3 = new TCanvas("Pion Multiplicity","Pion Multiplicity",600,600);
  c3->cd();
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  TH1F *pion_mult = new TH1F ("pion_mult", "pion_mult", 100, 0, 100);
  pion_mult->GetXaxis()->SetTitle("Pion Multiplicity");
  pion_mult->GetYaxis()->SetTitle("Number of Events");
  pion_mult->GetYaxis()->SetTitleOffset(1.4);
  pion_mult->SetTitle("BELLE e^{+}e^{-}, 0.1 < p_{T}< 4 GeV/c");
  pion_mult->SetStats(0);
  pion_mult->Draw();
    
  // Proton Multiplicity
  TCanvas *c4 = new TCanvas("Proton Multiplicity","Proton Multiplicity",600,600);
  c4->cd();
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  TH1F *proton_mult = new TH1F ("proton_mult", "proton_mult", 20, 0, 20);
  proton_mult->GetXaxis()->SetTitle("Proton Multiplicity");
  proton_mult->GetYaxis()->SetTitleOffset(1.4);
  proton_mult->GetYaxis()->SetTitle("Number of Events");
  proton_mult->SetTitle("BELLE e^{+}e^{-}, 0.1 < p_{T}< 4 GeV/c");
  proton_mult->SetStats(0);
  proton_mult->Draw();
    
  // Kaon Multiplicity
  TCanvas *c5 = new TCanvas("Kaon Multiplicity","Kaon Multiplicity",600,600);
  c5->cd();
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  TH1F *kaon_mult = new TH1F ("kaon_mult", "kaon_mult", 20, 0, 20);
  kaon_mult->GetXaxis()->SetTitle("Kaon Multiplicity");
  kaon_mult->GetYaxis()->SetTitle("Number of Events");
  kaon_mult->GetYaxis()->SetTitleOffset(1.4);
  kaon_mult->SetTitle("BELLE e^{+}e^{-}, 0.1 < p_{T}< 4 GeV/c");
  kaon_mult->SetStats(0);
  kaon_mult->Draw();
    
  // Total Multiplicity
  TCanvas *c6 = new TCanvas("Total Multiplicity"," Total Multiplicity",600,600);
  c6->cd();
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  TH1F *total_mult = new TH1F ("total_mult", "total_mult", 100, 0, 100);
  total_mult->GetXaxis()->SetTitle("Total Multiplicity");
  total_mult->GetYaxis()->SetTitle("Number of Events");
  total_mult->GetYaxis()->SetTitleOffset(1.4);
  total_mult->SetTitle("BELLE e^{+}e^{-}, 0.1 < p_{T}< 4 GeV/c");
  total_mult->SetStats(0);
  total_mult->Draw();
  
  // THETA
  TCanvas *c7 = new TCanvas("Theta","Theta",600,600);
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  TH1F *pion_theta = new TH1F ("pion_theta","#theta distribution of particles",100,-PI,0);
  TH1F *proton_theta = new TH1F ("proton_theta","proton",100,-PI,0);
  TH1F *kaon_theta = new TH1F ("proton_theta","pion",100,-PI,0);
  pion_theta->SetStats(0);
  pion_theta->Draw();
  pion_theta->GetXaxis()->SetTitle("#theta");
  pion_theta->GetYaxis()->SetTitle("Number of Events");
  pion_theta->GetYaxis()->SetTitleOffset(1.4);
  proton_theta->SetLineColor(kRed);
  proton_theta->Draw("same");
  kaon_theta->SetLineColor(kGreen);
  kaon_theta->Draw("same");
  
 
  for (Int_t i=0;i<nevent_process;i++) {
    if (i%1000==0) cout <<i<<"/"<<nevent_process<<endl;
    t1->GetEntry(i);
    
    int nparticles = nParticle;
    if (nparticles<mult) continue;
      
      int count_pion = 0;
      int count_proton = 0;
      int count_kaon = 0;
      int count_total = 0;
    for ( int j=0;j<nparticles;j++ )
    {
        if(pid[j] == PION){pion_eta->Fill(eta[j]);pion_pt->Fill(pt[j]);count_pion++;count_total++; pion_theta->Fill(theta[j]);}
        if(pid[j] == PROTON){proton_eta->Fill(eta[j]);proton_pt->Fill(pt[j]);count_proton++;count_total++; proton_theta->Fill(theta[j]);}
        if(pid[j] == KAON){kaon_eta->Fill(eta[j]);kaon_pt->Fill(pt[j]);count_kaon++;count_total++; kaon_theta->Fill(theta[j]);}
    }
      c3->cd();
      pion_mult -> Fill(count_pion);
      c4->cd();
      proton_mult -> Fill(count_proton);
      c5->cd();
      kaon_mult -> Fill(count_kaon);
      c6->cd();
      total_mult-> Fill(count_total);
  }
    // end of loop over events

    pion_eta->Scale(1./pion_eta->GetEntries());
    proton_eta->Scale(1./proton_eta->GetEntries());
    kaon_eta->Scale(1./kaon_eta->GetEntries());
    
    pion_pt->Scale(1./pion_pt->GetEntries());
    proton_pt->Scale(1./proton_pt->GetEntries());
    kaon_pt->Scale(1./kaon_pt->GetEntries());
    
    pion_theta->Scale(1./pion_theta->GetEntries());
    proton_theta->Scale(1./proton_theta->GetEntries());
    kaon_theta->Scale(1./kaon_theta->GetEntries());
    /*
    c1->cd();
    
    TLegend *leg = new TLegend(0.6526846,0.6474695,0.8808725,0.8010471,NULL,"brNDC");
    leg->SetBorderSize(1);
    leg->SetTextSize(0.055);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(1001);
    leg->SetHeader("Particles","C"); // option "C" allows to center the header
    //leg->SetTextFont(132);
    leg->AddEntry(pion_eta,"Pion","l");
    //leg->SetTextFont(132);
    leg->AddEntry(proton_eta,"Proton","l");
    //leg->SetTextFont(132);
    leg->AddEntry(kaon_eta,"Kaon","l");
    //leg->SetTextFont(132);
    //leg->SetTextSize(0.05);
    leg->Draw();   
    c1->Modified();
    c1->Update();
    */
    c7->cd();
    TLegend* leg2 = new TLegend(0.8,0.7,0.7,0.8);
    leg2->SetHeader("Particles","C"); // option "C" allows to center the header
    leg2->SetTextFont(132);
    leg2->AddEntry(pion_theta,"Pion","l");
    leg2->SetTextFont(132);
    leg2->AddEntry(proton_theta,"Proton","l");
    leg2->SetTextFont(132);
    leg2->AddEntry(kaon_theta,"Kaon","l");
    leg2->SetTextFont(132);
    leg2->Draw();
    
    c7->Modified();
    c7->Update();

    c2->cd();
    TLegend* leg1 = new TLegend(0.8,0.7,0.7,0.8);
    leg1->SetHeader("Particles","C"); // option "C" allows to center the header
    leg1->AddEntry(pion_pt,"Pion","l");
    leg1->AddEntry(proton_pt,"Proton","l");
    leg1->AddEntry(kaon_pt,"Kaon","l");
    leg1->Draw();
    
    c2->Modified();
    c2->Update();
 
    c1->Draw();
    c8->Draw();
    c9->Draw();
    c1->SaveAs("pion_eta.pdf");
    c8->SaveAs("proton_eta.pdf");
    c9->SaveAs("kaon_eta.pdf");
    //c2->SaveAs("pt.pdf");
    //c3->SaveAs("pion_mult.pdf");
    //c4->SaveAs("proton_mult.pdf");
    //c5->SaveAs("kaon_mult.pdf");
    //c6->SaveAs("total_mult.pdf");
    //c7->SaveAs("theta.pdf");
    

}




