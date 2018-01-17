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
    const int nBinsY = 80;
    Double_t binsYLog[nBinsY+1];
    const Double_t logLow = .0000005;
    const Double_t logHi = .3;
    getLogBins(logLow, logHi, nBinsY, binsYLog);
    
    // declare linear binning for probability of eta
    Double_t binsYLin[nBinsY+1];
    const Double_t linLow = 0;
    const Double_t linHi = 0.1;
    getLinBins(linLow,linHi,nBinsY,binsYLin);
    
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
    
    // assign binning ranges for different data sets
    if(datalabel == "PYTHIA8"){nBinsMult = 100; multHi = 100; etaHi = 3; ptHi = 60;}
    if(datalabel == "LEP1"){nBinsMult = 80; multHi = 80; etaHi = 2.5; ptHi = 60;}
    if(datalabel == "LEP2"){nBinsMult = 90; multHi = 90; etaHi = 1.75; ptHi = 60;}
    
    // declare bins
    Double_t binsMult[nBinsMult+1];
    Double_t binsPt[nBinsPt+1];
    Double_t binsEta[nBinsEta+1];
    etaLow = -etaHi;
    
    // generate the binning
    getLinBins(multLow, multHi, nBinsMult, binsMult);
    getLinBins(ptLow, ptHi, nBinsPt, binsPt);
    getLinBins(etaLow, etaHi, nBinsEta, binsEta);
    
    // Declare Histograms
    TH2F *hemptymult = new TH2F("",";Multiplicity;Probability",nBinsMult,binsMult,nBinsY,binsYLog);
    TH1F* amult = new TH1F("amult","",nBinsMult,binsMult);
    TH1F* cmult = new TH1F("cmult","",nBinsMult,binsMult);
    TH1F* nmult = new TH1F("nmult","",nBinsMult,binsMult);
    TH1F* lmult = new TH1F("lmult","",nBinsMult,binsMult); // leptons
    
    TH2F *hemptypt = new TH2F("",";pt;Probability",nBinsPt,binsPt,nBinsY,binsYLog);
    TH1F* apt = new TH1F("apt","",nBinsPt,binsPt);
    TH1F* cpt = new TH1F("cpt","",nBinsPt,binsPt);
    TH1F* npt = new TH1F("npt","",nBinsPt,binsPt);
    
    TH2F *hemptyeta = new TH2F("",";eta;Probability",nBinsEta,binsEta,nBinsEta,binsYLin);
    TH1F* aeta = new TH1F("aeta","",nBinsEta,binsEta);
    TH1F* ceta = new TH1F("ceta","",nBinsEta,binsEta);
    TH1F* neta = new TH1F("neta","",nBinsEta,binsEta);
    
    // fill weights
    static const double weight = 0.000000001;
    
    // Load the data
    TFile *f = new TFile(inFileName,"READ");
    TTree *t1 = (TTree*)f->Get("t");
    TTree *ak4JetTree = (TTree*)f->Get("ak4JetTree");
    t1->AddFriend(ak4JetTree);
    //TTree *ak4JetTree = (TTree*)f->Get("ak4JetTree");
    //TTree *ak8JetTree = (TTree*)f->Get("ak8JetTree");
    
    //// main Tree ////
    Int_t nParticle;
    Int_t nChargedHadrons;
    Int_t nParticleChg;
    Float_t pt[50000];
    Float_t eta[50000];
    Float_t theta[50000];
    Float_t phi[50000];
    Int_t pwflag[50000];
    
    t1->SetBranchAddress("nParticle",&nParticle);
    t1->SetBranchAddress("nChargedHadrons",&nChargedHadrons);
    if(datalabel == "PYTHIA8")t1->SetBranchAddress("nParticleChg",&nParticleChg);
    t1->SetBranchAddress("pt",pt);
    t1->SetBranchAddress("eta",eta);
    t1->SetBranchAddress("theta",theta);
    t1->SetBranchAddress("phi",phi);
    t1->SetBranchAddress("pwflag",pwflag);
    
    // all entries and fill the histograms
    Int_t nevent = (Int_t)t1->GetEntries();
    Int_t nChg = 0;
    Int_t nLepton = 0;
    // fill the histograms
    for (Int_t i = 0; i<nevent; i++)
    {
        t1->GetEntry(i);
        
        if(datalabel == "PYTHIA8")nChg = nParticleChg;
        else nChg = nChargedHadrons;
        
        // fill the multiplicity
        amult->Fill(nParticle);
        cmult->Fill(nChg);
        nmult->Fill(nParticle-nChg);
        
        for (Int_t j=0;j<nParticle;j++)
        {
            // fill all particles
            aeta->Fill(eta[j],weight); apt->Fill(pt[j],weight);
            // fill charged hadrons
            if(pwflag[j]==0){ ceta->Fill(eta[j],weight); cpt->Fill(pt[j],weight);}
            // fill the neutral hadrons and photons
            if(pwflag[j]!=0){ neta->Fill(eta[j],weight); npt->Fill(pt[j],weight);}
            // fill the leptons
            if(pwflag[j] == 1 || pwflag[j] == 2){nLepton++;}
        }
        lmult->Fill(nLepton);
    }
    
    // Normalize all to unity
    amult->Scale(1./amult->Integral());
    cmult->Scale(1./cmult->Integral());
    nmult->Scale(1./nmult->Integral());
    lmult->Scale(1./lmult->Integral());
    
    apt->Scale(1./apt->Integral());
    cpt->Scale(1./cpt->Integral());
    npt->Scale(1./npt->Integral());
    
    aeta->Scale(1./aeta->Integral());
    ceta->Scale(1./ceta->Integral());
    neta->Scale(1./neta->Integral());
    
    // Plot all, neutral, charged multiplicity
    TCanvas *allmult = new TCanvas ("allmult","allmult",600,600);
    gPad->SetLogy();
    xjjroot::sethempty(hemptymult,0,0.3);
    hemptymult->Draw();
    xjjroot::setthgrstyle(amult, kBlack, 20, 1.2, kBlack, 1, 1, -1, -1, -1);
    amult->Draw("pe same");
    xjjroot::setthgrstyle(cmult, kRed, 21, 1.2, kRed, 1, 1, -1, -1, -1);
    cmult->Draw("pe same");
    xjjroot::setthgrstyle(nmult, kBlue, 22, 1.2, kBlue, 1, 1, -1, -1, -1);
    nmult->Draw("pe same");
    xjjroot::setthgrstyle(lmult, kGreen, 22, 1.2, kGreen, 1, 1, -1, -1, -1);
    lmult->Draw("pe same");
    TLegend* leg = new TLegend(0.42,0.80,0.85,0.895);
    xjjroot::setleg(leg);
    leg->AddEntry(amult,"All","p");
    leg->AddEntry(cmult,"Charged Hadrons","p");
    leg->AddEntry(nmult,"Neutral Hadrons + Photons","p");
    leg->AddEntry(lmult,"Leptons","p");
    leg->Draw();
    xjjroot::drawtex(0.2,0.876,datalabel);
    allmult->SaveAs(Form("pdfDir/%s_mult.pdf",datalabel.Data()));
    
    // Plot all, neutral, charged pt
    TCanvas *allpt = new TCanvas("allpt","allpt",600,600);
    gPad->SetLogy();
    xjjroot::sethempty(hemptypt,0,0.3);
    hemptypt->Draw();
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
    xjjroot::sethempty(hemptyeta,0,0.3);
    hemptyeta->Draw();
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
    
    outFile_p->Close();
    delete outFile_p;
    
    delete etaleg;
    delete neta;
    delete ceta;
    delete aeta;
    delete hemptyeta;
    delete alleta;
    
    delete ptleg;
    delete npt;
    delete cpt;
    delete apt;
    delete hemptypt;
    delete allpt;
    
    delete leg;
    delete nmult;
    delete cmult;
    delete amult;
    delete hemptymult;
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


