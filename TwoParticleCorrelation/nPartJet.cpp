//
//  npartJet.cpp
//  
//
//  Created by Anthony Badea on 10/4/17.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include "TLorentzVector.h"
#include "TCanvas.h"
void nPartJet(TString filename = "/home/abadea/Documents/20171002/alephDataPaths_LEP1.root")
{
    TFile *f = new TFile(filename);
    TTree *t1 = (TTree*)f->Get("t");
    t1->AddFriend("ak4JetTree");
    TTree *ak4JetTree = (TTree*)f->Get("ak4JetTree");
    TTree *ak8JetTree = (TTree*)f->Get("ak8JetTree");
    
    //// main Tree ////
    Int_t nParticle;
    Float_t pt[50000];
    Float_t eta[50000];
    Float_t theta[50000];
    Float_t phi[50000];
    Int_t pwflag[50000];

    t1->SetBranchAddress("nParticle",&nParticle);
    t1->SetBranchAddress("pt",pt);
    t1->SetBranchAddress("eta",eta);
    t1->SetBranchAddress("theta",theta);
    t1->SetBranchAddress("phi",phi);
    t1->SetBranchAddress("pwflag",pwflag);
    
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
    
    
    // Get number of events
    Int_t nevent = (Int_t)t1->GetEntries();
    Int_t nevent_process = nevent;

    
    // 4 cuts of nref:nParticle 
    // A = (10,15):(55,65)
    // B = (25,30):(55,65)
    // C = (20,25):(20,30)
    // D = (5,10):(20,30)
    
    static const  int nRegions = 4;
    Int_t nref_l[nRegions]= {10,25,20,5};
    Int_t nref_h[nRegions]= {15,30,25,10};
    Int_t nPart_l[nRegions] = {55,55,20,20};
    Int_t nPart_h[nRegions] = {65,65,30,30};
    
    static const  int nEvents = 3;
    Int_t eventNo[nRegions][nEvents] = {{4919,7415,5262},{1984,5376,1396},{1819,796,7882},{11,12,36}};
    Int_t runNo[nRegions][nEvents] = {{12182,12466,12551},{12326,12679,11114},{12051,12053,12109},{11447,11590,11449}};
    TCanvas * c0 = new TCanvas("c", "c", 600, 600);
    for(Int_t i = 0; i<nRegions;i++)
    {
        for(Int_t j=0;j<nEvents;j++)
        {
            
            // Note that all we really need is the EventNo but we got the event number using the other parameters so the other parameters serve as a check to make sure we are picking up the correct event
            
            // for all jets
            c0->Clear();
            TString param = Form("EventNo==%d&&RunNo==%d&&nref>%d&&nref<%d&&nParticle>%d&&nParticle<%d",eventNo[i][j],runNo[i][j],nref_l[i],nref_h[i],nPart_l[i],nPart_h[i]);
	    cout<<"hi"<<endl;
            t1->Draw("jteta:jtphi",param,"COLZ");
            c0->SaveAs(Form("jet_region%d_%d_%d.pdf",i,eventNo[i][j],runNo[i][j]));
            
            // all hadronic particles
            c0->Clear();
            TString param1 = Form("EventNo==%d&&RunNo==%d&&pwflag==0&&nParticle>%d&&nParticle<%d",eventNo[i][j],runNo[i][j],nPart_l[i],nPart_h[i]);
            t1->Draw("eta:phi",param1,"COLZ");
            c0->SaveAs(Form("hadronic_region%d_%d_%d.pdf",i,eventNo[i][j],runNo[i][j]));
            
            // all non-hadronic particles
            c0->Clear();
            TString param2 = Form("EventNo==%d&&RunNo==%d&&pwflag!=0&&nParticle>%d&&nParticle<%d",eventNo[i][j],runNo[i][j],nPart_l[i],nPart_h[i]);
            t1->Draw("eta:phi",param2,"COLZ");
            c0->SaveAs(Form("nonhadronic_region%d_%d_%d.pdf",i,eventNo[i][j],runNo[i][j]));
            
            // all particles
            c0->Clear();
            TString param3 = Form("EventNo==%d&&RunNo==%d&&nParticle>%d&&nParticle<%d",eventNo[i][j],runNo[i][j],nPart_l[i],nPart_h[i]);
            t1->Draw("eta:phi",param3,"COLZ");
            c0->SaveAs(Form("allparticle_region%d_%d_%d.pdf",i,eventNo[i][j],runNo[i][j]));
            
            // jet pt
            c0->Clear();
            TString param4 = Form("EventNo==%d&&RunNo==%d&&nref>%d&&nref<%d&&nParticle>%d&&nParticle<%d",eventNo[i][j],runNo[i][j],nref_l[i],nref_h[i],nPart_l[i],nPart_h[i]);
            t1->Draw("jtpt",param4);
            c0->SaveAs(Form("jtpt_region%d_%d_%d.pdf",i,eventNo[i][j],runNo[i][j]));
            
            // particles pt
            c0->Clear();
            t1->Draw("pt",param4);
            c0->SaveAs(Form("pt_region%d_%d_%d.pdf",i,eventNo[i][j],runNo[i][j]));
        }
    }
}
