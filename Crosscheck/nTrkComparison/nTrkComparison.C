#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TMath.h"
#include <string>
#include <vector>
#include <iostream>
#include "TProfile.h"

void nTrkComparison(int nTrkLow = 0, int nTrkHigh = 999){

  TH1::SetDefaultSumw2();

  std::string alephCut = "TMath::Abs(eta)<1.8 && pt>0.2 && pwflag==0 && charge!=0";
  std::string cmsGENCut = "TMath::Abs(eta)<2.4 && pt>0.4 && chg!=0";
  
  TFile * alephData = TFile::Open("/data/cmcginn/StudyMultSamples/ALEPH/LEP1/20171016/alephDataPaths_LEP1_1995.root","read");
  TFile * alephMC = TFile::Open("/data/cmcginn/StudyMultSamples/ALEPH/MC/20171016/alephMCRecoAfterCutPaths_1997_20171016.root","read");
  TFile * cmsMC = TFile::Open("/data/abaty/EpEmStudies/CMS_5TeVpp_MC_MBPythia8_Forests/HiForestAOD_1.root","read");

  TTree * Taleph = (TTree*) alephData->Get("t");
  TTree * TalephMCReco = (TTree*) alephMC->Get("t");
  TTree * TalephMCGen = (TTree*)  alephMC->Get("tgen");
  TTree * TcmsMCGen = (TTree*) cmsMC->Get("HiGenParticleAna/hi");

  TH1D * Ntrk_AlephOffline = new TH1D("Ntrk_AlephOffline","Ntrk_AlephOffline",300,0,300);
  Taleph->Draw(Form("Sum$(%s)>>Ntrk_AlephOffline",alephCut.c_str()),Form("Sum$(%s)>%d && Sum$(%s)<%d",alephCut.c_str(),nTrkLow,alephCut.c_str(),nTrkHigh));
  std::cout << "<n_Trk^Offline> for bin from " << nTrkLow << " to " << nTrkHigh << ":    " << Ntrk_AlephOffline->GetMean() << std::endl; 

  TProfile * nTrk_AlephMCGen = new TProfile("Ntrk_AlephMCGen","Ntrk_AlephMCGen",300,0,300);   
  TProfile * nTrk_AlephMCGen_Cut = new TProfile("Ntrk_AlephMCGen_Cut","Ntrk_AlephMCGen_Cut",300,0,300);   

  float pt[10000], ptG[10000];
  float eta[10000], etaG[10000];
  float pwflag[10000], pwflagG[10000];
  float charge[10000], chargeG[10000];
  int nParticle, nParticleG;
  int process;

  TalephMCReco->SetBranchAddress("nParticle",&nParticle);
  TalephMCReco->SetBranchAddress("pt",&pt);
  TalephMCReco->SetBranchAddress("eta",&eta);
  TalephMCReco->SetBranchAddress("pwflag",&pwflag);
  TalephMCReco->SetBranchAddress("charge",&charge);
  TalephMCGen->SetBranchAddress("nParticle",&nParticleG);
  TalephMCGen->SetBranchAddress("pt",&ptG);
  TalephMCGen->SetBranchAddress("eta",&etaG);
  TalephMCGen->SetBranchAddress("pwflag",&pwflagG);
  TalephMCGen->SetBranchAddress("charge",&chargeG);
  TalephMCGen->SetBranchAddress("process",&process);


  for(int i = 0; i<TalephMCReco->GetEntries(); i++){
  //for(int i = 0; i<100000; i++){
    if(i%1000==0) std::cout << i << "/" << TalephMCReco->GetEntries() << std::endl;
    TalephMCReco->GetEntry(i);
    TalephMCGen->GetEntry(i);
    if(process!=5) continue;
    int sumReco = 0;
    for(int j = 0; j<nParticle; j++){
      if(TMath::Abs(eta[j])>1.8 || pt[j]<0.2 || charge[j]==0 || pwflag[j]!=0) continue;
      sumReco++;
    }
    int sumGen = 0;
    int sumGenCut = 0;
    for(int j = 0; j<nParticleG; j++){
      sumGen++;
      if(TMath::Abs(etaG[j])>1.8 || ptG[j]<0.2 || chargeG[j]==0 || pwflagG[j]!=0) continue;
      sumGenCut++;
    }
    if(sumReco>0){
      nTrk_AlephMCGen->Fill(sumReco,sumReco/(double)sumGen);
      nTrk_AlephMCGen_Cut->Fill(sumReco,sumReco/(double)sumGenCut);
    }
  }
  nTrk_AlephMCGen->Print("All");
  nTrk_AlephMCGen_Cut->Print("All");
  nTrk_AlephMCGen->Draw();

 
  TProfile * nTrk_CmsMCGen = new TProfile("Ntrk_CmsMCGen","Ntrk_CmsMCGen",300,0,300);   
  std::vector<float> * ptV = 0;
  std::vector<float> * etaV = 0;
  std::vector<float> * chgV = 0;
 
  TcmsMCGen->SetBranchAddress("pt",&ptV);
  TcmsMCGen->SetBranchAddress("eta",&etaV);
  TcmsMCGen->SetBranchAddress("chg",&chgV);
  
  for(int i = 0; i<TcmsMCGen->GetEntries(); i++){
    if(i%1000==0) std::cout << i << "/" << TcmsMCGen->GetEntries() << std::endl;
    TcmsMCGen->GetEntry(i);
    int sumGen = 0;
    int sumGenCut = 0;
    for(unsigned int j = 0; j<ptV->size(); j++){
      sumGen++;
      if(TMath::Abs(etaV->at(j))>2.4 || ptV->at(j)<0.4 || chgV->at(j)==0) continue;
      sumGenCut++;
    }
    if(sumGen>0)  nTrk_CmsMCGen->Fill(sumGen,sumGenCut/(double)sumGen);
  }
  nTrk_CmsMCGen->Print("All");
 
   
  TH1D * conversion = new TH1D("conversion","conversion",300,0,300);
  TH1D * conversion_weighted = new TH1D("conversion_weighted","conversion_weighted",300,0,300);
  for(int i = 1; i<conversion->GetSize()-1; i++){
    std::cout << nTrk_CmsMCGen->GetBinContent(i)*nTrk_AlephMCGen->GetBinContent(i) << std::endl;
    if(nTrk_AlephMCGen->GetBinContent(i)>0){
      conversion->SetBinContent(i,nTrk_CmsMCGen->GetBinContent((int)(1.0/nTrk_AlephMCGen->GetBinContent(i)*i))*1.0/nTrk_AlephMCGen->GetBinContent(i));
      conversion->SetBinError(i,TMath::Power(TMath::Power(nTrk_CmsMCGen->GetBinError((int)(i/nTrk_AlephMCGen->GetBinContent(i))),2)+TMath::Power(nTrk_AlephMCGen->GetBinError(i)/nTrk_AlephMCGen->GetBinContent(i),2),0.5));
      conversion_weighted->SetBinContent(i,i*nTrk_CmsMCGen->GetBinContent((int)(1.0/nTrk_AlephMCGen->GetBinContent(i)*i))*1.0/nTrk_AlephMCGen->GetBinContent(i));
      conversion_weighted->SetBinError(i,i*TMath::Power(TMath::Power(nTrk_CmsMCGen->GetBinError((int)(i/nTrk_AlephMCGen->GetBinContent(i))),2)+TMath::Power(nTrk_AlephMCGen->GetBinError(i)/nTrk_AlephMCGen->GetBinContent(i),2),0.5));
    }
  }
  conversion->Print("All");

  TCanvas * c1 = new TCanvas("c1","c1",800,600);
  Ntrk_AlephOffline->GetXaxis()->SetRangeUser(0,100); 
  Ntrk_AlephOffline->GetXaxis()->SetTitle("nTrk"); 
  Ntrk_AlephOffline->SetMarkerStyle(8); 
  Ntrk_AlephOffline->Draw("p"); 
  c1->SaveAs("img/nTrkOffline_ALEPH.png");
  c1->SaveAs("img/nTrkOffline_ALEPH.pdf");
  c1->SaveAs("img/nTrkOffline_ALEPH.C");
   
  nTrk_AlephMCGen_Cut->GetXaxis()->SetRangeUser(0,100);
  nTrk_AlephMCGen_Cut->GetYaxis()->SetRangeUser(0,1);
  nTrk_AlephMCGen_Cut->GetYaxis()->SetTitle("efficiency");
  nTrk_AlephMCGen_Cut->GetXaxis()->SetTitle("nTrk"); 
  nTrk_AlephMCGen_Cut->SetMarkerStyle(8);
  nTrk_AlephMCGen_Cut->Draw("p"); 
  c1->SaveAs("img/nTrk_AlephMCGen_Cut.png");
  c1->SaveAs("img/nTrk_AlephMCGen_Cut.pdf");
  c1->SaveAs("img/nTrk_AlephMCGen_Cut.C");
  
  nTrk_AlephMCGen->GetXaxis()->SetRangeUser(0,100);
  nTrk_AlephMCGen->GetYaxis()->SetRangeUser(0,1);
  nTrk_AlephMCGen->GetYaxis()->SetTitle("efficiency+acceptance");
  nTrk_AlephMCGen->GetXaxis()->SetTitle("nTrk"); 
  nTrk_AlephMCGen->SetMarkerStyle(8);
  nTrk_AlephMCGen->Draw("p"); 
  c1->SaveAs("img/nTrk_AlephMCGen.png");
  c1->SaveAs("img/nTrk_AlephMCGen.pdf");
  c1->SaveAs("img/nTrk_AlephMCGen.C");
  
  nTrk_CmsMCGen->GetXaxis()->SetRangeUser(0,100);
  nTrk_CmsMCGen->GetYaxis()->SetRangeUser(0,1);
  nTrk_CmsMCGen->GetYaxis()->SetTitle("CMS acceptance");
  nTrk_CmsMCGen->GetXaxis()->SetTitle("nGen"); 
  nTrk_CmsMCGen->SetMarkerStyle(8);
  nTrk_CmsMCGen->Draw("p"); 
  c1->SaveAs("img/nTrk_CMSMCGen.png");
  c1->SaveAs("img/nTrk_CMSMCGen.pdf");
  c1->SaveAs("img/nTrk_CMSMCGen.C");

  conversion->GetXaxis()->SetRangeUser(0,100);
  conversion->GetYaxis()->SetRangeUser(0.5,3);
  conversion->GetYaxis()->SetTitle("efficiency+acceptance+CMS acceptance");
  conversion->GetXaxis()->SetTitle("nTrk"); 
  conversion->SetMarkerStyle(8);
  conversion->Draw("p"); 
  c1->SaveAs("img/conversion.png");
  c1->SaveAs("img/conversion.pdf");
  c1->SaveAs("img/conversion.C");
  
  conversion_weighted->GetXaxis()->SetRangeUser(0,100);
  conversion_weighted->GetYaxis()->SetRangeUser(0,50);
  conversion_weighted->GetYaxis()->SetTitle("nTrk_{CMS}^{Corr}");
  conversion_weighted->GetXaxis()->SetTitle("nTrk_{ALEPH}"); 
  conversion_weighted->SetMarkerStyle(8);
  conversion_weighted->Draw("p"); 
  c1->SaveAs("img/conversion_weighted.png");
  c1->SaveAs("img/conversion_weighted.pdf");
  c1->SaveAs("img/conversion_weighted.C");
  
}
