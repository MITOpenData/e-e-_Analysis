//c and c++ dependencies
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

//root dependencies
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TPad.h"
#include "TStyle.h"

//local headers
#include "../include/xjjrootuti.h"
#include "../include/getLogBins.h"
#include "../include/getLinBins.h"


////////// Current Plots /////////
// multiplicity, pt, eta, phi, rapidity, eta-phi scatter plot (all in beam and thrust axis)
//////////////////////////////////
int eeplots
(
 const TString inFileName, // input file
 const TString datalabel // Text in upper-left corner
)
{
    // setting ROOT defaults
    xjjroot::setgstyle();
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    // declare log binning for probability of pt/multiplicity
    const int nBinsP = 80;
    Double_t binsPLog[nBinsP+1];
    const Double_t logLow = .0000005;
    const Double_t logHi = .3;
    getLogBins(logLow, logHi, nBinsP, binsPLog);
    
    // declare linear binning for probability of eta
    Double_t binsYLin[nBinsP+1];
    const Double_t linLow = 0;
    const Double_t linHi = 0.15;
    getLinBins(linLow,linHi,nBinsP,binsYLin);
    
    // declare binning for x-axis
    int nBinsMult = 60;
    Double_t multHi = 0;
    Double_t multLow = 0;
    
    int nBinsPt = 60;
    Double_t ptHi = 0;
    Double_t ptLow = 0;
    
    int nBinsEta = 60;
    Double_t etaHi = 0;
    Double_t etaLow = 0;
    
    int nBinsY = 60;
    Double_t yHi = 8;
    Double_t yLow = 0;

    int nBinsPhi = 60;
    Double_t phiHi = TMath::Pi();
    Double_t phiLow = -TMath::Pi()/2.;
    
    // assign binning ranges for different data sets
    if(datalabel == "PYTHIA8"){nBinsMult = 100; multHi = 100; etaHi = 3; ptHi = 60;}
    if(datalabel == "LEP1"){nBinsMult = 80; multHi = 80; etaHi = 2.5; ptHi = 60;}
    if(datalabel == "LEP2"){nBinsMult = 90; multHi = 90; etaHi = 1.75; ptHi = 60;}
    etaLow = -etaHi;
    
    // declare bins
    Double_t binsMult[nBinsMult+1];
    Double_t binsPt[nBinsPt+1];
    Double_t binsEta[nBinsEta+1];
    Double_t binsY[nBinsY+1];
    Double_t binsPhi[nBinsPhi+1];
    
    // generate the binning
    getLinBins(multLow, multHi, nBinsMult, binsMult);
    getLinBins(ptLow, ptHi, nBinsPt, binsPt);
    getLinBins(etaLow, etaHi, nBinsEta, binsEta);
    getLinBins(yLow, yHi, nBinsY, binsY);
    getLinBins(phiLow, phiHi, nBinsPhi, binsPhi);
    
    // Declare Histograms
    TH2F *hempty_mult = new TH2F("",";Multiplicity;Probability",nBinsMult,binsMult,nBinsP,binsPLog);
    TH1F* a_mult = new TH1F("a_mult","",nBinsMult,binsMult);
    TH1F* c_mult = new TH1F("c_mult","",nBinsMult,binsMult);
    TH1F* n_mult = new TH1F("n_mult","",nBinsMult,binsMult);
    //TH1F* lmult = new TH1F("lmult","",nBinsMult,binsMult); // leptons

    TH2F *hempty_pt = new TH2F("",";pt;Probability",nBinsPt,binsPt,nBinsP,binsPLog);
    TH1F* a_pt = new TH1F("a_pt","",nBinsPt,binsPt);
    TH1F* c_pt = new TH1F("c_pt","",nBinsPt,binsPt);
    TH1F* n_pt = new TH1F("n_pt","",nBinsPt,binsPt);
    
    TH2F *hempty_pt_wrtThr = new TH2F("",";pt_wrtThr;Probability",nBinsPt,binsPt,nBinsP,binsPLog);
    TH1F* a_pt_wrtThr = new TH1F("a_pt_wrtThr","",nBinsPt,binsPt);
    TH1F* c_pt_wrtThr = new TH1F("c_pt_wrtThr","",nBinsPt,binsPt);
    TH1F* n_pt_wrtThr = new TH1F("n_pt_wrtThr","",nBinsPt,binsPt);
    
    TH2F *hempty_eta = new TH2F("",";eta;Probability",nBinsEta,binsEta,nBinsP,binsYLin);
    TH1F* a_eta = new TH1F("a_eta","",nBinsEta,binsEta);
    TH1F* c_eta = new TH1F("c_eta","",nBinsEta,binsEta);
    TH1F* n_eta = new TH1F("n_eta","",nBinsEta,binsEta);
    
    TH2F *hempty_eta_wrtThr = new TH2F("",";eta_wrtThr;Probability",nBinsEta,binsEta,nBinsP,binsYLin);
    TH1F* a_eta_wrtThr = new TH1F("a_eta_wrtThr","",nBinsEta,binsEta);
    TH1F* c_eta_wrtThr = new TH1F("c_eta_wrtThr","",nBinsEta,binsEta);
    TH1F* n_eta_wrtThr = new TH1F("n_eta_wrtThr","",nBinsEta,binsEta);
    
    TH2F *hempty_phi = new TH2F("",";phi;Probability",nBinsPhi,binsPhi,nBinsP,binsYLin);
    TH1F* a_phi = new TH1F("a_phi","",nBinsPhi,binsPhi);
    TH1F* c_phi = new TH1F("c_phi","",nBinsPhi,binsPhi);
    TH1F* n_phi = new TH1F("n_phi","",nBinsPhi,binsPhi);
    
    TH2F *hempty_phi_wrtThr = new TH2F("",";phi_wrtThr;Probability",nBinsPhi,binsPhi,nBinsP,binsYLin);
    TH1F* a_phi_wrtThr = new TH1F("a_phi_wrtThr","",nBinsPhi,binsPhi);
    TH1F* c_phi_wrtThr = new TH1F("c_phi_wrtThr","",nBinsPhi,binsPhi);
    TH1F* n_phi_wrtThr = new TH1F("n_phi_wrtThr","",nBinsPhi,binsPhi);
    
    TH2F *hempty_y = new TH2F("",";y;Probability",nBinsY,binsY,nBinsP,binsYLin);
    TH1F* a_y = new TH1F("a_y","",nBinsY,binsY);
    TH1F* c_y = new TH1F("c_y","",nBinsY,binsY);
    TH1F* n_y = new TH1F("n_y","",nBinsY,binsY);
    
    TH2F *hempty_y_wrtThr = new TH2F("",";y_wrtThr;Probability",nBinsY,binsY,nBinsP,binsYLin);
    TH1F* a_y_wrtThr = new TH1F("a_y_wrtThr","",nBinsY,binsY);
    TH1F* c_y_wrtThr = new TH1F("c_y_wrtThr","",nBinsY,binsY);
    TH1F* n_y_wrtThr = new TH1F("n_y_wrtThr","",nBinsY,binsY);
    
    TH2F *a_eta_phi = new TH2F("a_eta_phi","#eta-#phi of all particles;#eta;#phi;",nBinsEta,binsEta,nBinsPhi,binsPhi);
    a_eta_phi->SetTitle("#eta-#phi of all particles");
    TH2F *c_eta_phi = new TH2F("c_eta_phi","#eta-#phi of charged hadrons;#eta;#phi;",nBinsEta,binsEta,nBinsPhi,binsPhi);
    c_eta_phi->SetTitle("#eta-#phi of charged hadrons");
    TH2F *n_eta_phi = new TH2F("n_eta_phi","#eta-#phi of neutral hadrons and photons;#eta;#phi;",nBinsEta,binsEta,nBinsPhi,binsPhi);
    n_eta_phi->SetTitle("#eta-#phi of neutral hadrons and photons");
    
    TH2F *a_eta_phi_wrtThr = new TH2F("a_eta_phi_wrtThr","#eta-#phi of all particles w/respect to thrust axis;#eta;#phi;",nBinsEta,binsEta,nBinsPhi,binsPhi);
    a_eta_phi_wrtThr->SetTitle("#eta-#phi of all particles");
    TH2F *c_eta_phi_wrtThr = new TH2F("c_eta_phi_wrtThr","#eta-#phi of charged hadrons w/respect to thrust axis;#eta;#phi;",nBinsEta,binsEta,nBinsPhi,binsPhi);
    c_eta_phi_wrtThr->SetTitle("#eta-#phi of charged hadrons");
    TH2F *n_eta_phi_wrtThr = new TH2F("n_eta_phi_wrtThr","#eta-#phi of neutral hadrons and photons w/respect to thrust axis;#eta;#phi;",nBinsEta,binsEta,nBinsPhi,binsPhi);
    n_eta_phi_wrtThr->SetTitle("#eta-#phi of neutral hadrons and photons");
    
    // fill weights
    static const double weight = 0.0000000000001;
    
    // Load the data
    TFile *f = new TFile(inFileName,"READ");
    TTree *t1 = (TTree*)f->Get("t");
    TTree *ak4ESchemeJetTree = (TTree*)f->Get("ak4ESchemeJetTree");
    t1->AddFriend(ak4ESchemeJetTree);
    //TTree *ak4JetTree = (TTree*)f->Get("ak4JetTree");
    //TTree *ak8JetTree = (TTree*)f->Get("ak8JetTree");
    
    //// main Tree ////
    static const int nMaxPart = 10000;
    
    Int_t nParticle;
    Int_t nChargedHadrons;
    Float_t Energy;
    Float_t mass[nMaxPart];
    
    Float_t pt[nMaxPart];       Float_t pt_wrtThr[nMaxPart];
    Float_t eta[nMaxPart];      Float_t eta_wrtThr[nMaxPart];
    //Float_t theta[nMaxPart];    Float_t theta_wrtThr[nMaxPart];
    Float_t phi[nMaxPart];      Float_t phi_wrtThr[nMaxPart];
    
    Float_t px[nMaxPart];
    Float_t py[nMaxPart];
    Float_t pz[nMaxPart];
    Float_t rap[nMaxPart];
    Int_t pwflag[nMaxPart];
    
    t1->SetBranchAddress("nParticle",&nParticle);
    t1->SetBranchAddress("nChargedHadrons",&nChargedHadrons);
    t1->SetBranchAddress("Energy",&Energy);
    t1->SetBranchAddress("pt",pt);      t1->SetBranchAddress("pt_wrtThr",pt_wrtThr);
    t1->SetBranchAddress("eta",eta);        t1->SetBranchAddress("eta_wrtThr",eta_wrtThr);
    //t1->SetBranchAddress("theta",theta);      t1->SetBranchAddress("theta",theta);
    t1->SetBranchAddress("phi",phi);        t1->SetBranchAddress("phi_wrtThr",phi_wrtThr);
    t1->SetBranchAddress("pwflag",pwflag);
    t1->SetBranchAddress("px",px);
    t1->SetBranchAddress("py",py);
    t1->SetBranchAddress("pz",pz);
    t1->SetBranchAddress("rap",rap);
    t1->SetBranchAddress("mass",mass);
    
    // all entries and fill the histograms
    Int_t nevent = (Int_t)t1->GetEntries();
    
    // basic variable declarations
    //Int_t nLepton = 0;
    Float_t min_ZPole_Energy = 91.0;
    Float_t max_ZPole_Energy = 91.5;
    TLorentzVector v;
    // fill the histograms
    for (Int_t i = 0; i<nevent; i++)
    {
        t1->GetEntry(i);
        
        // for LEP2 Energy (91,91.5) to focus on Z pole
        if(datalabel == "LEP2"){if(Energy<min_ZPole_Energy || Energy>max_ZPole_Energy){continue;}}
 
        // fill the multiplicity
        a_mult->Fill(nParticle);
        c_mult->Fill(nChargedHadrons);
        n_mult->Fill(nParticle-nChargedHadrons);
        
        for (Int_t j=0;j<nParticle;j++)
        {
            v.SetPtEtaPhiM(pt_wrtThr[j],eta_wrtThr[j],phi_wrtThr[j],mass[j]);
            
            // fill all particles
            a_pt->Fill(pt[j],weight);    a_pt_wrtThr->Fill(pt_wrtThr[j],weight);
            a_eta->Fill(eta[j],weight);  a_eta_wrtThr->Fill(eta_wrtThr[j],weight);
            a_phi->Fill(phi[j],weight);  a_phi_wrtThr->Fill(phi_wrtThr[j],weight);
            a_y->Fill(rap[j],weight);    a_y_wrtThr->Fill(v.Rapidity(),weight);
            a_eta_phi->Fill(eta[j],phi[j],weight);        a_eta_phi_wrtThr->Fill(eta_wrtThr[j],phi_wrtThr[j],weight);
            
            // fill charged hadrons
            if(pwflag[j]==0)
            {
                c_pt->Fill(pt[j],weight);    c_pt_wrtThr->Fill(pt_wrtThr[j],weight);
                c_eta->Fill(eta[j],weight);  c_eta_wrtThr->Fill(eta_wrtThr[j],weight);
                c_phi->Fill(phi[j],weight);  c_phi_wrtThr->Fill(phi_wrtThr[j],weight);
                c_y->Fill(rap[j],weight);    c_y_wrtThr->Fill(v.Rapidity(),weight);
                c_eta_phi->Fill(eta[j],phi[j],weight);        c_eta_phi_wrtThr->Fill(eta_wrtThr[j],phi_wrtThr[j],weight);
            }
            
            // fill the neutral hadrons and photons
            if(pwflag[j]!=0)
            {
                n_pt->Fill(pt[j],weight);    n_pt_wrtThr->Fill(pt_wrtThr[j],weight);
                n_eta->Fill(eta[j],weight);  n_eta_wrtThr->Fill(eta_wrtThr[j],weight);
                n_phi->Fill(phi[j],weight);  n_phi_wrtThr->Fill(phi_wrtThr[j],weight);
                n_y->Fill(rap[j],weight);    n_y_wrtThr->Fill(v.Rapidity(),weight);
                n_eta_phi->Fill(eta[j],phi[j],weight);        n_eta_phi_wrtThr->Fill(eta_wrtThr[j],phi_wrtThr[j],weight);
            }
            // fill the leptons
            //if(pwflag[j] == 1 || pwflag[j] == 2){nLepton++;}
        }
        //lmult->Fill(nLepton);
    }
    
    // Normalize all to unity
    a_mult->Scale(1./a_mult->Integral());
    c_mult->Scale(1./c_mult->Integral());
    n_mult->Scale(1./n_mult->Integral());
    //lmult->Scale(1./lmult->Integral());
    
    a_pt->Scale(1./a_pt->Integral());     a_pt_wrtThr->Scale(1./a_pt_wrtThr->Integral());
    c_pt->Scale(1./c_pt->Integral());     c_pt_wrtThr->Scale(1./c_pt_wrtThr->Integral());
    n_pt->Scale(1./n_pt->Integral());     n_pt_wrtThr->Scale(1./n_pt_wrtThr->Integral());
    
    a_eta->Scale(1./a_eta->Integral());       a_eta_wrtThr->Scale(1./a_eta_wrtThr->Integral());
    c_eta->Scale(1./c_eta->Integral());       c_eta_wrtThr->Scale(1./c_eta_wrtThr->Integral());
    n_eta->Scale(1./n_eta->Integral());       n_eta_wrtThr->Scale(1./n_eta_wrtThr->Integral());
    
    a_phi->Scale(1./a_phi->Integral());       a_phi_wrtThr->Scale(1./a_phi_wrtThr->Integral());
    c_phi->Scale(1./c_phi->Integral());       c_phi_wrtThr->Scale(1./c_phi_wrtThr->Integral());
    n_phi->Scale(1./n_phi->Integral());       n_phi_wrtThr->Scale(1./n_phi_wrtThr->Integral());
    
    a_y->Scale(1./a_y->Integral());       a_y_wrtThr->Scale(1./a_y_wrtThr->Integral());
    c_y->Scale(1./c_y->Integral());       c_y_wrtThr->Scale(1./c_y_wrtThr->Integral());
    n_y->Scale(1./n_y->Integral());       n_y_wrtThr->Scale(1./n_y_wrtThr->Integral());
    
    a_eta_phi->Scale(1./a_eta_phi->Integral());     a_eta_phi_wrtThr->Scale(1./a_eta_phi_wrtThr->Integral());
    c_eta_phi->Scale(1./c_eta_phi->Integral());     c_eta_phi_wrtThr->Scale(1./c_eta_phi_wrtThr->Integral());
    n_eta_phi->Scale(1./n_eta_phi->Integral());     n_eta_phi_wrtThr->Scale(1./n_eta_phi_wrtThr->Integral());
    
    // Plot all, neutral, charged multiplicity
    TCanvas *allmult = new TCanvas ("allmult","allmult",600,600);
    gPad->SetLogy();
    xjjroot::sethempty(hempty_mult,0,0.3);
    hempty_mult->Draw();
    xjjroot::setthgrstyle(a_mult, kBlack, 20, 1.2, kBlack, 1, 1, -1, -1, -1);
    a_mult->Draw("pe same");
    xjjroot::setthgrstyle(c_mult, kRed, 21, 1.2, kRed, 1, 1, -1, -1, -1);
    c_mult->Draw("pe same");
    xjjroot::setthgrstyle(n_mult, kBlue, 22, 1.2, kBlue, 1, 1, -1, -1, -1);
    n_mult->Draw("pe same");
    //xjjroot::setthgrstyle(lmult, kGreen, 22, 1.2, kGreen, 1, 1, -1, -1, -1);
    //lmult->Draw("pe same");
    TLegend* multleg = new TLegend(0.42,0.80,0.85,0.895);
    xjjroot::setleg(multleg);
    multleg->AddEntry(a_mult,"All","p");
    multleg->AddEntry(c_mult,"Charged Hadrons","p");
    multleg->AddEntry(n_mult,"Neutral Hadrons + Photons","p");
    //leg->AddEntry(lmult,"Leptons","p");
    multleg->Draw();
    xjjroot::drawtex(0.2,0.876,datalabel);
    allmult->SaveAs(Form("pdfDir/%s_mult.pdf",datalabel.Data()));
    
    // Plot all, neutral, charged pt
    TCanvas *allpt = new TCanvas("allpt","allpt",600,600);
    gPad->SetLogy();
    xjjroot::sethempty(hempty_pt,0,0.3);
    hempty_pt->Draw();
    xjjroot::setthgrstyle(a_pt, kBlack, 20, 1.2, kRed, 1, 1, -1, -1, -1);
    a_pt->Draw("pe same");
    xjjroot::setthgrstyle(c_pt, kRed, 21, 1.2, kRed, 1, 1, -1, -1, -1);
    c_pt->Draw("pe same");
    xjjroot::setthgrstyle(n_pt, kBlue, 22, 1.2, kBlue, 1, 1, -1, -1, -1);
    n_pt->Draw("pe same");
    TLegend* ptleg = new TLegend(0.42,0.7,0.85,0.88);
    xjjroot::setleg(ptleg);
    ptleg->AddEntry(a_pt,"All","p");
    ptleg->AddEntry(c_pt,"Charged Hadrons","p");
    ptleg->AddEntry(n_pt,"Neutral Hadrons + Photons","p");
    ptleg->Draw();
    xjjroot::drawtex(0.25,0.876,datalabel);
    allpt->SaveAs(Form("pdfDir/%s_pt.pdf",datalabel.Data()));
    
    // Plot all, neutral, charged eta
    TCanvas *alleta = new TCanvas("alleta","alleta",600,600);
    xjjroot::sethempty(hempty_eta,0,0.3);
    hempty_eta->Draw();
    xjjroot::setthgrstyle(a_eta, kBlack, 20, 1.2, kRed, 1, 1, -1, -1, -1);
    a_eta->Draw("pe same");
    xjjroot::setthgrstyle(c_eta, kRed, 21, 1.2, kRed, 1, 1, -1, -1, -1);
    c_eta->Draw("pe same");
    xjjroot::setthgrstyle(n_eta, kBlue, 22, 1.2, kBlue, 1, 1, -1, -1, -1);
    n_eta->Draw("pe same");
    TLegend* etaleg = new TLegend(0.42,0.7,0.85,0.88);
    xjjroot::setleg(etaleg);
    etaleg->AddEntry(a_eta,"All","p");
    etaleg->AddEntry(c_eta,"Charged Hadrons","p");
    etaleg->AddEntry(n_eta,"Neutral Hadrons + Photons","p");
    etaleg->Draw();
    xjjroot::drawtex(0.2,0.876,datalabel);
    alleta->SaveAs(Form("pdfDir/%s_eta.pdf",datalabel.Data()));
    
    // Plot all, neutral, charged rapidity
    TCanvas *ally = new TCanvas("ally","ally",600,600);
    xjjroot::sethempty(hempty_y,0,0.3);
    hempty_y->Draw();
    xjjroot::setthgrstyle(a_y, kBlack, 20, 1.2, kRed, 1, 1, -1, -1, -1);
    a_y->Draw("pe same");
    xjjroot::setthgrstyle(c_y, kRed, 21, 1.2, kRed, 1, 1, -1, -1, -1);
    c_y->Draw("pe same");
    xjjroot::setthgrstyle(n_y, kBlue, 22, 1.2, kBlue, 1, 1, -1, -1, -1);
    n_y->Draw("pe same");
    TLegend* yleg = new TLegend(0.42,0.7,0.85,0.88);
    xjjroot::setleg(yleg);
    yleg->AddEntry(a_y,"All","p");
    yleg->AddEntry(c_y,"Charged Hadrons","p");
    yleg->AddEntry(n_y,"Neutral Hadrons + Photons","p");
    yleg->Draw();
    xjjroot::drawtex(0.2,0.876,datalabel);
    ally->SaveAs(Form("pdfDir/%s_y.pdf",datalabel.Data()));
    
    TCanvas *alleta_phi = new TCanvas("alleta_phi","alleta_phi",600,600);
    xjjroot::sethempty(a_eta_phi,0,0);
    a_eta_phi->Draw("colz");
    xjjroot::drawtex(0.2,0.876,datalabel);
    alleta_phi->SetRightMargin(0.13);
    alleta_phi->SaveAs(Form("pdfDir/%s_a_eta_phi.pdf",datalabel.Data()));
    
    TCanvas *chareta_phi = new TCanvas("chareta_phi","chareta_phi",600,600);
    xjjroot::sethempty(c_eta_phi,0,0);
    c_eta_phi->Draw("colz same");
    xjjroot::drawtex(0.2,0.876,datalabel);
    chareta_phi->SetRightMargin(0.13);
    chareta_phi->SaveAs(Form("pdfDir/%s_c_eta_phi.pdf",datalabel.Data()));
    
    TCanvas *neuteta_phi = new TCanvas("neuteta_phi","neuteta_phi",600,600);
    xjjroot::sethempty(n_eta_phi,0,0);
    n_eta_phi->Draw("colz");
    xjjroot::drawtex(0.2,0.876,datalabel);
    neuteta_phi->SetRightMargin(0.13);
    neuteta_phi->SaveAs(Form("pdfDir/%s_n_eta_phi.pdf",datalabel.Data()));
    
    TFile* outFile_p = new TFile(Form("inputs/qualityCheck/outFile_%s.root",datalabel.Data()), "RECREATE");
    //write all, first arg name ("" == declaration name), second arg overwrites buffer saves in file
    
    a_mult->Write("", TObject::kOverwrite);
    c_mult->Write("", TObject::kOverwrite);
    n_mult->Write("", TObject::kOverwrite);
    
    a_pt->Write("", TObject::kOverwrite);   a_pt_wrtThr->Write("", TObject::kOverwrite);
    c_pt->Write("", TObject::kOverwrite);   c_pt_wrtThr->Write("", TObject::kOverwrite);
    n_pt->Write("", TObject::kOverwrite);   n_pt_wrtThr->Write("", TObject::kOverwrite);
    
    a_eta->Write("", TObject::kOverwrite);  a_eta_wrtThr->Write("", TObject::kOverwrite);
    c_eta->Write("", TObject::kOverwrite);  c_eta_wrtThr->Write("", TObject::kOverwrite);
    n_eta->Write("", TObject::kOverwrite);  n_eta_wrtThr->Write("", TObject::kOverwrite);
    
    a_phi->Write("", TObject::kOverwrite);  a_phi_wrtThr->Write("", TObject::kOverwrite);
    c_phi->Write("", TObject::kOverwrite);  c_phi_wrtThr->Write("", TObject::kOverwrite);
    n_phi->Write("", TObject::kOverwrite);  n_phi_wrtThr->Write("", TObject::kOverwrite);
    
    a_y->Write("", TObject::kOverwrite);    a_y_wrtThr->Write("", TObject::kOverwrite);
    c_y->Write("", TObject::kOverwrite);    c_y_wrtThr->Write("", TObject::kOverwrite);
    n_y->Write("", TObject::kOverwrite);    n_y_wrtThr->Write("", TObject::kOverwrite);
    
    a_eta_phi->Write("", TObject::kOverwrite);  a_eta_phi_wrtThr->Write("", TObject::kOverwrite);
    c_eta_phi->Write("", TObject::kOverwrite);  c_eta_phi_wrtThr->Write("", TObject::kOverwrite);
    n_eta_phi->Write("", TObject::kOverwrite);  n_eta_phi_wrtThr->Write("", TObject::kOverwrite);
    
    outFile_p->Close();
    delete outFile_p;
    
    delete n_eta_phi_wrtThr;
    delete c_eta_phi_wrtThr;
    delete a_eta_phi_wrtThr;
    
    delete n_eta_phi;
    delete c_eta_phi;
    delete a_eta_phi;
    
    //delete yleg_wrtThr;
    delete n_y_wrtThr;
    delete c_y_wrtThr;
    delete a_y_wrtThr;
    delete hempty_y_wrtThr;
    //delete ally_wrtThr;
    
    delete yleg;
    delete n_y;
    delete c_y;
    delete a_y;
    delete hempty_y;
    delete ally;
    
    //delete phileg_wrtThr;
    delete n_phi_wrtThr;
    delete c_phi_wrtThr;
    delete a_phi_wrtThr;
    delete hempty_phi_wrtThr;
    //delete allphi_wrtThr;
    
    //delete phileg;
    delete n_phi;
    delete c_phi;
    delete a_phi;
    delete hempty_phi;
    //delete allphi;
    
    //delete etaleg_wrtThr;
    delete n_eta_wrtThr;
    delete c_eta_wrtThr;
    delete a_eta_wrtThr;
    delete hempty_eta_wrtThr;
    //delete alleta_wrtThr;
    
    delete etaleg;
    delete n_eta;
    delete c_eta;
    delete a_eta;
    delete hempty_eta;
    delete alleta;
    
    //delete ptleg_wrtThr;
    delete n_pt_wrtThr;
    delete c_pt_wrtThr;
    delete a_pt_wrtThr;
    delete hempty_pt_wrtThr;
    //delete allpt_wrtThr;
    
    delete ptleg;
    delete n_pt;
    delete c_pt;
    delete a_pt;
    delete hempty_pt;
    delete allpt;
    
    delete multleg;
    delete n_mult;
    delete c_mult;
    delete a_mult;
    delete hempty_mult;
    delete allmult;
    
    f->Close();
    delete f;
    
    return 0;
}

int main(int argc, char* argv[])
{
    if(argc != 3)
    {
        std::cout << "Usage: ./eeplots.exe <inFileName> <datalabel>" << std::endl;
        return 1;
    }
    
    int retVal = 0;
    if(argc == 3) retVal += eeplots(argv[1], argv[2]);
    return retVal;
}


