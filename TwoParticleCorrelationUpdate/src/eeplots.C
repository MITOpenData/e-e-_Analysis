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
    // setting style
    xjjroot::setgstyle();
    
    // declare binning for probability
    const int nBinsY = 80;
    Double_t binsY[nBinsY+1];
    const Double_t logLow = .00005;
    const Double_t logHi = .3;
    getLogBins(logLow, logHi, nBinsY, binsY);
    
    // declare binning for axis
    const int nBinsMult = 60;
    Double_t binsMult[nBinsMult+1];
    const Double_t multLow = 0;
    const Double_t multHi = 60;
    getLinBins(multLow, multHi, nBinsMult, binsMult);
    
    const int nBinsPt = 60;
    Double_t binsPt[nBinsPt+1];
    const Double_t ptLow = 0;
    const Double_t ptHi = 60;
    getLinBins(ptLow, ptHi, nBinsPt, binsPt);
    
    
    //
    // open the file
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
    
    /*
    //// ak4JetTree ////
    Int_t nref4;
    Float_t jtpt4[100000];
    Float_t jteta4[100000];
    Float_t jtphi4[100000];
    
    t1->SetBranchAddress("nref",&nref4);
    t1->SetBranchAddress("jtpt",jtpt4);
    t1->SetBranchAddress("jteta",jteta4);
    t1->SetBranchAddress("jtphi",jtphi4);

    //// ak8JetTree ////
    Int_t nref8;
    Float_t jtpt8[100000];
    Float_t jteta8[100000];
    Float_t jtphi8[100000];
    
    ak8JetTree->SetBranchAddress("nref",&nref8);
    ak8JetTree->SetBranchAddress("jtpt",jtpt8);
    ak8JetTree->SetBranchAddress("jteta",jteta8);
    ak8JetTree->SetBranchAddress("jtphi",jtphi8);
     */
    
    // Plot all, neutral, charged multiplicity
    TCanvas *allmult = new TCanvas ("allmult","allmult",600,600);
    gPad->SetLogy();
	//TH2F *hempty = new TH2F("",";Multiplicity;Probability",1,0,60,1,0,0.15);
    TH2F *hempty = new TH2F("",";Multiplicity;Probability",nBinsMult,binsMult,nBinsY,binsY);
    xjjroot::sethempty(hempty,0,0.3);
    hempty->Draw();
    t1->Draw("nParticle>>amult","","goff");
    TH1F* amult = (TH1F*)gDirectory->Get("amult");
    amult->Scale(1./amult->Integral());
    xjjroot::setthgrstyle(amult, kBlack, 20, 1.2, kBlack, 1, 1, -1, -1, -1);
    amult->Draw("pe same");
    if(datalabel == "PYTHIA8")t1->Draw("nParticleChg>>cmult","","goff");
    else t1->Draw("nChargedHadrons>>cmult","","goff");
    TH1F* cmult = (TH1F*)gDirectory->Get("cmult");
    cmult->Scale(1./cmult->Integral());
    xjjroot::setthgrstyle(cmult, kRed, 21, 1.2, kRed, 1, 1, -1, -1, -1);
    cmult->Draw("pe same");
    if(datalabel == "PYTHIA8")t1->Draw("(nParticle - nParticleChg)>>nmult","","goff");
    else t1->Draw("(nParticle - nChargedHadrons)>>nmult","","goff");
    TH1F* nmult = (TH1F*)gDirectory->Get("nmult");
    nmult->Scale(1./nmult->Integral());
    xjjroot::setthgrstyle(nmult, kBlue, 22, 1.2, kBlue, 1, 1, -1, -1, -1);
    nmult->Draw("pe same");
    TLegend* leg = new TLegend(0.42,0.76,0.85,0.89);
    xjjroot::setleg(leg);
    leg->AddEntry(amult,"All","p");
    leg->AddEntry(cmult,"Charged Hadrons","p");
    leg->AddEntry(nmult,"Neutral Hadrons + Photons","p");
    leg->Draw();
    xjjroot::drawtex(0.2,0.876,datalabel);
    allmult->SaveAs(Form("pdfDir/%s_mult.pdf",datalabel.Data()));
    
    // Plot all, neutral, charged pt
    TCanvas *allmom = new TCanvas("allmom","allmom",600,600);
    gPad->SetLogy();
    //TH2F *hemptyp = new TH2F("",";pt;Probability",1,0,60,1,0,0.3);
    TH2F *hemptyp = new TH2F("",";pt;Probability",nBinsPt,binsPt,nBinsY,binsY);
    xjjroot::sethempty(hemptyp,0,0.3);
    hemptyp->Draw();
    t1->Draw(("pt>>amom(" + std::to_string(nBinsPt) + ", " + std::to_string(binsPt[0]) + ", " + std::to_string(binsPt[nBinsPt]) + ")").c_str(),"(1)*.000001","goff");
    TH1F* amom = (TH1F*)gDirectory->Get("amom");
    amom->Scale(1./amom->Integral());
    xjjroot::setthgrstyle(amom, kBlack, 20, 1.2, kRed, 1, 1, -1, -1, -1);
    amom->Draw("pe same");
    t1->Draw(("pt>>cmom(" + std::to_string(nBinsPt) + ", " + std::to_string(binsPt[0]) + ", " + std::to_string(binsPt[nBinsPt]) + ")").c_str(),"(pwflag==0)*.000001","goff");
    TH1F* cmom = (TH1F*)gDirectory->Get("cmom");
    cmom->Scale(1./cmom->Integral());
    xjjroot::setthgrstyle(cmom, kRed, 21, 1.2, kRed, 1, 1, -1, -1, -1);
    cmom->Draw("pe same");
    t1->Draw(("pt>>nmom(" + std::to_string(nBinsPt) + ", " + std::to_string(binsPt[0]) + ", " + std::to_string(binsPt[nBinsPt]) + ")").c_str(),"(pwflag!=0)*.000001","goff");
    TH1F* nmom = (TH1F*)gDirectory->Get("nmom");
    nmom->Scale(1./nmom->Integral());
    xjjroot::setthgrstyle(nmom, kBlue, 22, 1.2, kBlue, 1, 1, -1, -1, -1);
    nmom->Draw("pe same");
    TLegend* pleg = new TLegend(0.42,0.7,0.85,0.88);
    xjjroot::setleg(pleg);
    pleg->AddEntry(amom,"All","p");
    pleg->AddEntry(cmom,"Charged Hadrons","p");
    pleg->AddEntry(nmom,"Neutral Hadrons + Photons","p");
    pleg->Draw();
    xjjroot::drawtex(0.2,0.876,datalabel);
    allmom->SaveAs(Form("pdfDir/%s_mom.pdf",datalabel.Data()));
    
    // Plot all, neutral, charged eta
    TCanvas *alleta = new TCanvas("alleta","alleta",600,600);
    TH2F *hemptyeta = 0;
    if(datalabel == "LEP2") hemptyeta = new TH2F("",";eta;Probability",1,-2.5,2.5,1,0,0.1);
    else hemptyeta = new TH2F("",";eta;Probability",1,-1.7,1.7,1,0,0.1);
    xjjroot::sethempty(hemptyeta,0,0.3);
    hemptyeta->Draw();
    t1->Draw("eta>>aeta","(1)*.000001","goff");
    TH1F* aeta = (TH1F*)gDirectory->Get("aeta");
    aeta->Scale(1./aeta->Integral());
    xjjroot::setthgrstyle(aeta, kBlack, 20, 1.2, kRed, 1, 1, -1, -1, -1);
    aeta->Draw("pe same");
    t1->Draw("eta>>ceta","(pwflag==0)*.000001","goff");
    TH1F* ceta = (TH1F*)gDirectory->Get("ceta");
    ceta->Scale(1./ceta->Integral());
    xjjroot::setthgrstyle(ceta, kRed, 21, 1.2, kRed, 1, 1, -1, -1, -1);
    ceta->Draw("pe same");
    t1->Draw("eta>>neta","(pwflag!=0)*.000001","goff");
    TH1F* neta = (TH1F*)gDirectory->Get("neta");
    neta->Scale(1./neta->Integral());
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
    amom->Write("", TObject::kOverwrite);
    cmom->Write("", TObject::kOverwrite);
    nmom->Write("", TObject::kOverwrite);
    aeta->Write("", TObject::kOverwrite);
    ceta->Write("", TObject::kOverwrite);
    neta->Write("", TObject::kOverwrite);
    
    /*
     
    // Plot jet eta
    TCanvas *jeteta = new TCanvas("jeteta","jeteta",600,600);
    TH2F *hemptye = new TH2F("",";Jet #eta;Probability",1,-3,3,1,0,0.05);
    xjjroot::sethempty(hemptye,0,0.5);
    hemptye->Draw();
    t1->Draw("jteta>>jeta","","goff");
    TH1F* jeta = (TH1F*)gDirectory->Get("jeta");
    jeta->Scale(1./jeta->GetEntries());
    xjjroot::setthgrstyle(jeta, kBlack, 20, 1.2, kBlack, 1, 1, -1, -1, -1);
    jeta->Draw("pe same");
    xjjroot::drawtex(0.2,0.876,datalabel);

    // Plot jet pt
    TCanvas *jetpt = new TCanvas("jetpt","jetpt",600,600);
    TH2F *hemptyj = new TH2F("",";Jet pt;Probability",1,0,80,1,0,1);
    xjjroot::sethempty(hemptyj,0,0.3);
    hemptyj->Draw();
    t1->Draw("jtpt>>jpt","","goff");
    TH1F* jpt = (TH1F*)gDirectory->Get("jpt");
    jpt->Scale(1./jpt->GetEntries());
    xjjroot::setthgrstyle(jpt, kBlack, 20, 1.2, kBlack, 1, 1, -1, -1, -1);
    jpt->Draw("pe same");
    xjjroot::drawtex(0.2,0.876,datalabel);
     
     */
    
    outFile_p->Close();
    delete outFile_p;
    
    delete etaleg;
    delete neta;
    delete ceta;
    delete aeta;
    delete hemptyeta;
    delete alleta;
    
    delete pleg;
    delete nmom;
    delete cmom;
    delete amom;
    delete hemptyp;
    delete allmom;
    
    delete leg;
    delete nmult;
    delete cmult;
    delete amult;
    delete hempty;
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

