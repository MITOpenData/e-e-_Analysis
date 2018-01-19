//
// Plots variables from the LEP trees.
//

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
    const Double_t linHi = 0.1;
    getLinBins(linLow,linHi,nBinsP,binsYLin);
    
    // declare binning for x-axis
    int nBinsMult = 0;
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
    TH1F* amult = new TH1F("amult","",nBinsMult,binsMult);
    TH1F* cmult = new TH1F("cmult","",nBinsMult,binsMult);
    TH1F* nmult = new TH1F("nmult","",nBinsMult,binsMult);
    //TH1F* lmult = new TH1F("lmult","",nBinsMult,binsMult); // leptons
    std::cout<< __FILE__<<","<< __LINE__ <<std::endl;
    TH2F *hempty_pt = new TH2F("",";pt;Probability",nBinsPt,binsPt,nBinsP,binsPLog);
    TH1F* apt = new TH1F("apt","",nBinsPt,binsPt);
    TH1F* cpt = new TH1F("cpt","",nBinsPt,binsPt);
    TH1F* npt = new TH1F("npt","",nBinsPt,binsPt);
    std::cout<< __FILE__<<","<< __LINE__ <<std::endl;
    TH2F *hempty_eta = new TH2F("",";eta;Probability",nBinsEta,binsEta,nBinsP,binsYLin);
    TH1F* aeta = new TH1F("aeta","",nBinsEta,binsEta);
    TH1F* ceta = new TH1F("ceta","",nBinsEta,binsEta);
    TH1F* neta = new TH1F("neta","",nBinsEta,binsEta);
    std::cout<< __FILE__<<","<< __LINE__ <<std::endl;
    TH2F *hempty_y = new TH2F("",";y;Probability",nBinsY,binsY,nBinsP,binsYLin);
    TH1F* ay = new TH1F("ay","",nBinsY,binsY);
    TH1F* cy = new TH1F("cy","",nBinsY,binsY);
    TH1F* ny = new TH1F("ny","",nBinsY,binsY);
    
    std::cout<< __FILE__<<","<< __LINE__ <<std::endl;
    
    // don't plot in this file just fill
    //TH2F *hempty_aetaphi = new TH2F("","All Particles;#eta;#phi",nBinsEta,binsEta,nBinsPhi,binsPhi);
    TH2F *aetaphi = new TH2F("aetaphi","#eta-#phi of all particles;#eta;#phi;",nBinsEta,binsEta,nBinsPhi,binsPhi);
    aetaphi->SetTitle("#eta-#phi of all particles");
    //TH2F *hempty_cetaphi = new TH2F("","Charged Hadrons;#eta;#phi",nBinsEta,binsEta,nBinsPhi,binsPhi);
    TH2F *cetaphi = new TH2F("cetaphi","#eta-#phi of charged hadrons;#eta;#phi;",nBinsEta,binsEta,nBinsPhi,binsPhi);
    cetaphi->SetTitle("#eta-#phi of charged hadrons");
    //TH2F *hempty_netaphi = new TH2F("","Neutral Hadrons and Photons;#eta;#phi",nBinsEta,binsEta,nBinsPhi,binsPhi);
    TH2F *netaphi = new TH2F("netaphi","#eta-#phi of neutral hadrons and photons;#eta;#phi;",nBinsEta,binsEta,nBinsPhi,binsPhi);
    netaphi->SetTitle("#eta-#phi of neutral hadrons and photons");
    
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
    Float_t pt[nMaxPart];
    Float_t eta[nMaxPart];
    Float_t theta[nMaxPart];
    Float_t phi[nMaxPart];
    Int_t pwflag[nMaxPart];
    Float_t px[nMaxPart];
    Float_t py[nMaxPart];
    Float_t pz[nMaxPart];
    Float_t rap[nMaxPart];
    
    
    t1->SetBranchAddress("nParticle",&nParticle);
    t1->SetBranchAddress("nChargedHadrons",&nChargedHadrons);
    t1->SetBranchAddress("Energy",&Energy);
    t1->SetBranchAddress("pt",pt);
    t1->SetBranchAddress("eta",eta);
    t1->SetBranchAddress("theta",theta);
    t1->SetBranchAddress("phi",phi);
    t1->SetBranchAddress("pwflag",pwflag);
    t1->SetBranchAddress("px",px);
    t1->SetBranchAddress("py",py);
    t1->SetBranchAddress("pz",pz);
    t1->SetBranchAddress("rap",rap);
    
    // all entries and fill the histograms
    Int_t nevent = (Int_t)t1->GetEntries();
    
    // basic variable declarations
    //Int_t nLepton = 0;
    Float_t min_ZPole_Energy = 91.0;
    Float_t max_ZPole_Energy = 91.5;
    
    // fill the histograms
    for (Int_t i = 0; i<nevent; i++)
    {
        t1->GetEntry(i);
        
        // for LEP2 Energy (91,91.5) to focus on Z pole
        if(datalabel == "LEP2"){if(Energy<min_ZPole_Energy || Energy>max_ZPole_Energy){continue;}}
 
        // fill the multiplicity
        amult->Fill(nParticle);
        cmult->Fill(nChargedHadrons);
        nmult->Fill(nParticle-nChargedHadrons);
        
        for (Int_t j=0;j<nParticle;j++)
        {
            // fill all particles
            aeta->Fill(eta[j],weight); apt->Fill(pt[j],weight); ay->Fill(rap[j],weight); aetaphi->Fill(eta[j],phi[j],weight);
            // fill charged hadrons
            if(pwflag[j]==0){ ceta->Fill(eta[j],weight); cpt->Fill(pt[j],weight); cy->Fill(rap[j],weight); cetaphi->Fill(eta[j],phi[j],weight);}
            // fill the neutral hadrons and photons
            if(pwflag[j]!=0){ neta->Fill(eta[j],weight); npt->Fill(pt[j],weight); ny->Fill(rap[j],weight); netaphi->Fill(eta[j],phi[j],weight);}
            // fill the leptons
            //if(pwflag[j] == 1 || pwflag[j] == 2){nLepton++;}
        }
        //lmult->Fill(nLepton);
    }
    
    // Normalize all to unity
    amult->Scale(1./amult->Integral());
    cmult->Scale(1./cmult->Integral());
    nmult->Scale(1./nmult->Integral());
    //lmult->Scale(1./lmult->Integral());
    
    apt->Scale(1./apt->Integral());
    cpt->Scale(1./cpt->Integral());
    npt->Scale(1./npt->Integral());
    
    aeta->Scale(1./aeta->Integral());
    ceta->Scale(1./ceta->Integral());
    neta->Scale(1./neta->Integral());
    
    ay->Scale(1./ay->Integral());
    cy->Scale(1./cy->Integral());
    ny->Scale(1./ny->Integral());
    
    aetaphi->Scale(1./aetaphi->Integral());
    cetaphi->Scale(1./cetaphi->Integral());
    netaphi->Scale(1./netaphi->Integral());
    
    // Plot all, neutral, charged multiplicity
    TCanvas *allmult = new TCanvas ("allmult","allmult",600,600);
    gPad->SetLogy();
    xjjroot::sethempty(hempty_mult,0,0.3);
    hempty_mult->Draw();
    xjjroot::setthgrstyle(amult, kBlack, 20, 1.2, kBlack, 1, 1, -1, -1, -1);
    amult->Draw("pe same");
    xjjroot::setthgrstyle(cmult, kRed, 21, 1.2, kRed, 1, 1, -1, -1, -1);
    cmult->Draw("pe same");
    xjjroot::setthgrstyle(nmult, kBlue, 22, 1.2, kBlue, 1, 1, -1, -1, -1);
    nmult->Draw("pe same");
    //xjjroot::setthgrstyle(lmult, kGreen, 22, 1.2, kGreen, 1, 1, -1, -1, -1);
    //lmult->Draw("pe same");
    TLegend* multleg = new TLegend(0.42,0.80,0.85,0.895);
    xjjroot::setleg(multleg);
    multleg->AddEntry(amult,"All","p");
    multleg->AddEntry(cmult,"Charged Hadrons","p");
    multleg->AddEntry(nmult,"Neutral Hadrons + Photons","p");
    //leg->AddEntry(lmult,"Leptons","p");
    multleg->Draw();
    xjjroot::drawtex(0.2,0.876,datalabel);
    allmult->SaveAs(Form("pdfDir/%s_mult.pdf",datalabel.Data()));
    
    // Plot all, neutral, charged pt
    TCanvas *allpt = new TCanvas("allpt","allpt",600,600);
    gPad->SetLogy();
    xjjroot::sethempty(hempty_pt,0,0.3);
    hempty_pt->Draw();
    xjjroot::setthgrstyle(apt, kBlack, 20, 1.2, kRed, 1, 1, -1, -1, -1);
    apt->Draw("pe same");
    xjjroot::setthgrstyle(cpt, kRed, 21, 1.2, kRed, 1, 1, -1, -1, -1);
    cpt->Draw("pe same");
    xjjroot::setthgrstyle(npt, kBlue, 22, 1.2, kBlue, 1, 1, -1, -1, -1);
    npt->Draw("pe same");
    TLegend* ptleg = new TLegend(0.42,0.7,0.85,0.88);
    xjjroot::setleg(ptleg);
    ptleg->AddEntry(apt,"All","p");
    ptleg->AddEntry(cpt,"Charged Hadrons","p");
    ptleg->AddEntry(npt,"Neutral Hadrons + Photons","p");
    ptleg->Draw();
    xjjroot::drawtex(0.25,0.876,datalabel);
    allpt->SaveAs(Form("pdfDir/%s_pt.pdf",datalabel.Data()));
    
    // Plot all, neutral, charged eta
    TCanvas *alleta = new TCanvas("alleta","alleta",600,600);
    xjjroot::sethempty(hempty_eta,0,0.3);
    hempty_eta->Draw();
    xjjroot::setthgrstyle(aeta, kBlack, 20, 1.2, kRed, 1, 1, -1, -1, -1);
    aeta->Draw("pe same");
    xjjroot::setthgrstyle(ceta, kRed, 21, 1.2, kRed, 1, 1, -1, -1, -1);
    ceta->Draw("pe same");
    xjjroot::setthgrstyle(neta, kBlue, 22, 1.2, kBlue, 1, 1, -1, -1, -1);
    neta->Draw("pe same");
    TLegend* etaleg = new TLegend(0.42,0.7,0.85,0.88);
    xjjroot::setleg(etaleg);
    etaleg->AddEntry(aeta,"All","p");
    etaleg->AddEntry(ceta,"Charged Hadrons","p");
    etaleg->AddEntry(neta,"Neutral Hadrons + Photons","p");
    etaleg->Draw();
    xjjroot::drawtex(0.2,0.876,datalabel);
    alleta->SaveAs(Form("pdfDir/%s_eta.pdf",datalabel.Data()));
    
    // Plot all, neutral, charged rapidity
    TCanvas *ally = new TCanvas("ally","ally",600,600);
    xjjroot::sethempty(hempty_y,0,0.3);
    hempty_y->Draw();
    xjjroot::setthgrstyle(ay, kBlack, 20, 1.2, kRed, 1, 1, -1, -1, -1);
    ay->Draw("pe same");
    xjjroot::setthgrstyle(cy, kRed, 21, 1.2, kRed, 1, 1, -1, -1, -1);
    cy->Draw("pe same");
    xjjroot::setthgrstyle(ny, kBlue, 22, 1.2, kBlue, 1, 1, -1, -1, -1);
    ny->Draw("pe same");
    TLegend* yleg = new TLegend(0.42,0.7,0.85,0.88);
    xjjroot::setleg(yleg);
    yleg->AddEntry(ay,"All","p");
    yleg->AddEntry(cy,"Charged Hadrons","p");
    yleg->AddEntry(ny,"Neutral Hadrons + Photons","p");
    yleg->Draw();
    xjjroot::drawtex(0.2,0.876,datalabel);
    ally->SaveAs(Form("pdfDir/%s_y.pdf",datalabel.Data()));
    
    TCanvas *alletaphi = new TCanvas("alletaphi","alletaphi",600,600);
    xjjroot::sethempty(aetaphi,0,0);
    aetaphi->Draw("colz");
    xjjroot::drawtex(0.2,0.876,datalabel);
    alletaphi->SetRightMargin(0.13);
    alletaphi->SaveAs(Form("pdfDir/%s_aetaphi.pdf",datalabel.Data()));
    
    TCanvas *charetaphi = new TCanvas("charetaphi","charetaphi",600,600);
    xjjroot::sethempty(cetaphi,0,0);
    cetaphi->Draw("colz same");
    xjjroot::drawtex(0.2,0.876,datalabel);
    charetaphi->SetRightMargin(0.13);
    charetaphi->SaveAs(Form("pdfDir/%s_cetaphi.pdf",datalabel.Data()));
    
    TCanvas *neutetaphi = new TCanvas("neutetaphi","neutetaphi",600,600);
    xjjroot::sethempty(netaphi,0,0);
    netaphi->Draw("colz");
    xjjroot::drawtex(0.2,0.876,datalabel);
    neutetaphi->SetRightMargin(0.13);
    neutetaphi->SaveAs(Form("pdfDir/%s_netaphi.pdf",datalabel.Data()));
    
    TFile* outFile_p = new TFile(Form("inputs/qualityCheck/outFile_%s.root",datalabel.Data()), "RECREATE");
    //write all, first arg name ("" == declaration name), second arg overwrites buffer saves in file
    
    amult->Write("", TObject::kOverwrite);
    cmult->Write("", TObject::kOverwrite);
    nmult->Write("", TObject::kOverwrite);
    apt->Write("", TObject::kOverwrite);
    cpt->Write("", TObject::kOverwrite);
    npt->Write("", TObject::kOverwrite);
    aeta->Write("", TObject::kOverwrite);
    ceta->Write("", TObject::kOverwrite);
    neta->Write("", TObject::kOverwrite);
    ay->Write("", TObject::kOverwrite);
    cy->Write("", TObject::kOverwrite);
    ny->Write("", TObject::kOverwrite);
    aetaphi->Write("", TObject::kOverwrite);
    cetaphi->Write("", TObject::kOverwrite);
    netaphi->Write("", TObject::kOverwrite);
    
    outFile_p->Close();
    delete outFile_p;
    
    delete netaphi;
    delete cetaphi;
    delete aetaphi;
    
    delete yleg;
    delete ny;
    delete cy;
    delete ay;
    delete hempty_y;
    delete ally;
    
    delete etaleg;
    delete neta;
    delete ceta;
    delete aeta;
    delete hempty_eta;
    delete alleta;
    
    delete ptleg;
    delete npt;
    delete cpt;
    delete apt;
    delete hempty_pt;
    delete allpt;
    
    delete multleg;
    delete nmult;
    delete cmult;
    delete amult;
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


