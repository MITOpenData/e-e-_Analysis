//c and c++ dependencies
#include <sstream>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <fstream>

//root dependencies
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TFormula.h>
#include <TNtuple.h>
#include "TDatime.h"

//local headers
#include "include/fourier.h"
#include "include/TPCNtupleData.h"
#include "include/Selection.h"
#include "include/smartJetName.h"

/********************************************************************************************************************/
// Two particle correlation analysis
//
// ridge_check_parallel.c
//
// Yen-Jie Lee, Gian Michele Innocenti, Anthony Badea
//
//
/********************************************************************************************************************/

using namespace std;

/********************************************************************************************************************/
// Main Analysis Routine
/********************************************************************************************************************/

int ridge_check_parallel
    (
     const std::string inFileName, // input file
     std::string outFileName    // output file
    )
{
    
    //////////////////// CREATE OUTFILE NAME ////////////////
    // create a outFileName if the one supplied is nothing
    // load the date
    TDatime* date = new TDatime();
    if(outFileName.size() == 0) outFileName = inFileName;
    if(outFileName.find(".root") != std::string::npos){outFileName.replace(outFileName.find(".root"), 5, "");}
    outFileName = outFileName + "_PLOTS_" + std::to_string(date->GetDate());
    
    // ROOT Global setting
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    // set up plots
    Selection s = Selection();
    TFile * output = TFile::Open(Form("%s_%d.root",outFileName.c_str(),(Int_t)s.doThrust),"recreate");
    TH2F * signal2PC[s.nMultBins];
    TH2F * bkgrnd2PC[s.nMultBins];
    TH2F * ratio2PC[s.nMultBins];
    TH1F * longRangeYield[s.nMultBins];
    float nSignalEvts[s.nMultBins] = {0};
    float nBkgrndEvts[s.nMultBins] = {0};
    
    TH1D * h_eff = new TH1D("h_eff","Selection Efficiency",100,0.0,1.0);;
    TH1D * h_phi = new TH1D("phi","phi",100,-TMath::Pi(),TMath::Pi());
    TH1D * h_eta = new TH1D("eta","eta",100,-5,5);
    TH1D * h_theta = new TH1D("theta","theta",100,0,TMath::Pi());
    TH1D * h_pt = new TH1D("pt","pt",100,0,10);
    TH1D * h_Ttheta = new TH1D("T_theta","T_theta",100,0,TMath::Pi());
    TH1D * h_Tphi = new TH1D("T_phi","T_phi",100,-TMath::Pi(),TMath::Pi());
    TH1D * h_Aj = new TH1D("h_Aj","h_Aj",50,0,0.5);
    TH1D * nEvtSigHist = new TH1D("nEvtSigHisto","nEvtSigHisto",10,0,10);
    TH1D * nEvtBkgHist = new TH1D("nEvtBkgHisto","nEvtBkgHisto",10,0,10);
    
    Float_t etaPlotRange = s.getEtaPlotRange();
    for(int i = 0; i<s.nMultBins; i++)
    {
        signal2PC[i] = new TH2F(Form("signal2PC_%d_%d",s.multBinsLow[i],s.multBinsHigh[i]),";#Delta#eta;#Delta#Phi",s.dEtaBins,-etaPlotRange,etaPlotRange,s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
        signal2PC[i]->Sumw2();
        bkgrnd2PC[i] = new TH2F(Form("bkgrnd2PC_%d_%d",s.multBinsLow[i],s.multBinsHigh[i]),";#Delta#eta;#Delta#Phi",s.dEtaBins,-etaPlotRange,etaPlotRange,s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
        bkgrnd2PC[i]->Sumw2();
        ratio2PC[i] = new TH2F(Form("ratio2PC_%d_%d",s.multBinsLow[i],s.multBinsHigh[i]),";#Delta#eta;#Delta#Phi",s.dEtaBins,-etaPlotRange,etaPlotRange,s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
        ratio2PC[i]->Sumw2();
        longRangeYield[i] = new TH1F(Form("longRangeYield_%d_%d",s.multBinsLow[i],s.multBinsHigh[i]),";#Delta#phi;Y(#Delta#Phi)",s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
        longRangeYield[i]->Sumw2();
    }
    TH1F * multiplicity = new TH1F("multiplicity",";nTrk;nEvents",200,0,200);
    
    std::string jtTreeName = "";
    if (s.jttree == 0) jtTreeName = "ak4ESchemeJetTree";
    if (s.jttree == 1) jtTreeName = "ak4WTAmodpSchemeJetTree";
    if (s.jttree == 2) jtTreeName = "ak8ESchemeJetTree";
    if (s.jttree == 3) jtTreeName = "ak8WTAmodpSchemeJetTree";
    
    if(inFileName.find(".root") != std::string::npos){
      TFile* temp_p = new TFile(inFileName.c_str(), "READ");
      jtTreeName = smartJetName(jtTreeName, temp_p);
      temp_p->Close();
      delete temp_p;
    }
    
    // files and variables for input
    // get the correct tree for the analysis
    std::string sideTree = "t";
    if (s.side_tree == 1) sideTree = "BoostedWTAR8Evt";
    
    TChain * t = new TChain("t");       t->Add(inFileName.c_str());
    TChain * side_t = new TChain(sideTree.c_str());       side_t->Add(inFileName.c_str());
    TChain * jt = new TChain(jtTreeName.c_str());       jt->Add(inFileName.c_str());

    TPCNtupleData data(s.doBelle, s.doThrust, s.side_tree);      setupTPCTree(t,side_t,jt,data);       data.setTPCTreeStatus(t,side_t);
    
    TChain * t_mix = new TChain("t");       t_mix->Add(inFileName.c_str());
    TChain * side_t_mix = new TChain(sideTree.c_str());       side_t_mix->Add(inFileName.c_str());
    TChain * jt_mix = new TChain(jtTreeName.c_str());       jt_mix->Add(inFileName.c_str());
    
    TPCNtupleData mix(s.doBelle, s.doThrust, s.side_tree);       setupTPCTree(t_mix,side_t_mix,jt_mix,mix);        mix.setTPCTreeStatus(t_mix,side_t_mix);
    
    // analysis
    Int_t nevent = (Int_t)t->GetEntries();
    if(s.doOneEvent) nevent = s.numEvents;

    /****************************************/
    // Main Event Loop
    /****************************************/
    for (Int_t i=0;i<nevent;i++) {
        t->GetEntry(i);
        side_t->GetEntry(i);
        jt->GetEntry(i);
        data.update();
        if (i%10000==0) std::cout <<i<<"/"<<nevent<<std::endl;
        
        // nTrk calculation
        Int_t nTrk = 0;
        
        //
        if(!s.donTrkThrust) nTrk = s.ridge_eventSelection(data.passesWW, data.nParticle, data.missP, data.pt, data.eta, data.nTPC, data.pwflag, data.nref, data.jtpt, data.jteta,s.donTrkThrust);
        if(s.donTrkThrust) nTrk = s.ridge_eventSelection(data.passesWW, data.nParticle, data.missP, data.pt_wrtThr, data.eta_wrtThr, data.nTPC, data.pwflag, data.nref, data.jtpt, data.jteta,s.donTrkThrust);
        //if(!s.doThrust) nTrk = s.ridge_eventSelection(data.passesWW, data.nParticle, data.missP, data.pt, data.eta, data.nTPC, data.pwflag, data.nref, data.jtpt, data.jteta);
        //if(s.doThrust) nTrk = s.ridge_eventSelection(data.passesWW, data.nParticle, data.missP, data.pt_wrtThr, data.eta_wrtThr, data.nTPC, data.pwflag, data.nref, data.jtpt, data.jteta);
        
        if( nTrk < 0) continue;
        Int_t histNum = s.histNum(nTrk);
	
        h_eff->Fill(nTrk/data.nParticle);
        h_Aj->Fill(s.fillAj);
        multiplicity->Fill(nTrk);
        nSignalEvts[histNum] += 1;
        
        h_Ttheta->Fill(data.TTheta);
        h_Tphi->Fill(data.TPhi);

        Float_t fillNumerator = 1.0;
        if (s.doPP) fillNumerator = data.pthatWeight;
        /****************************************/
        // S calculation using multiplicity cut //
        /****************************************/
        for ( Int_t j=0;j<data.nParticle;j++ )
        {
            if(!s.ridge_trackSelection(data.getPt(j),data.getEta(j),data.nTPC[j],data.pwflag[j],s.doThrust)) continue;
            
            h_phi->Fill(data.getPhi(j));
            h_eta->Fill(data.getEta(j));
            h_theta->Fill(data.getTheta(j));
            h_pt->Fill(data.getPt(j));
            
            Float_t angle1;
            if (s.doTheta) angle1 = data.getTheta(j); else angle1 = data.getEta(j);
            Float_t phi1 = data.getPhi(j);
            
            // Signal loop, calculate S correlation function
            for ( Int_t k=j+1;k<data.nParticle;k++ )
            {
                if(!s.ridge_trackSelection(data.getPt(k),data.getEta(k),data.nTPC[k],data.pwflag[k],s.doThrust)) continue;
                
                Float_t angle2;
                if (s.doTheta) angle2 = data.getTheta(k); else angle2 = data.getEta(k);
                Float_t phi2 = data.getPhi(k);
                
                signal2PC[histNum]->Fill(angle1-angle2,dphi(phi1,phi2),fillNumerator/(s.getDifferential())/nTrk);
                signal2PC[histNum]->Fill(angle1-angle2,dphi(phi2,phi1),fillNumerator/(s.getDifferential())/nTrk);
                signal2PC[histNum]->Fill(angle2-angle1,dphi(phi1,phi2),fillNumerator/(s.getDifferential())/nTrk);
                signal2PC[histNum]->Fill(angle2-angle1,dphi(phi2,phi1),fillNumerator/(s.getDifferential())/nTrk);
            }
        }
        
        /****************************************/
        // B calculation using multiplicity cut //
        /****************************************/
        for (Int_t nMix = 0; nMix<s.bkgrd_runs; nMix++)
        {
            Int_t selected=i+1;
            if (selected >= t->GetEntries()) selected = 0;
            t_mix->GetEntry(selected);
            side_t_mix->GetEntry(selected);
            jt_mix->GetEntry(selected);
        
            Int_t nTrk_mix;
            if(!s.donTrkThrust) nTrk_mix = s.ridge_eventSelection(mix.passesWW, mix.nParticle, mix.missP, mix.pt, mix.eta, mix.nTPC, mix.pwflag, mix.nref, mix.jtpt, mix.jteta,s.donTrkThrust);
            if(s.donTrkThrust) nTrk_mix = s.ridge_eventSelection(mix.passesWW, mix.nParticle, mix.missP, mix.pt_wrtThr, mix.eta_wrtThr, mix.nTPC, mix.pwflag, mix.nref, mix.jtpt, mix.jteta,s.donTrkThrust);
            Int_t histNum_mix = s.histNum(nTrk_mix);
            
            // Select a mixed event
            while (histNum_mix != histNum && s.isMixedEvent(nTrk, nTrk_mix, data.jteta[0], mix.jteta[0]))
            {
                selected++;
                if (selected >= t->GetEntries()) selected = 0;
                t_mix->GetEntry(selected);
                side_t_mix->GetEntry(selected);
                jt_mix->GetEntry(selected);
                
                if(!s.donTrkThrust) nTrk_mix = s.ridge_eventSelection(mix.passesWW, mix.nParticle, mix.missP, mix.pt, mix.eta, mix.nTPC, mix.pwflag, mix.nref, mix.jtpt, mix.jteta,s.donTrkThrust);
                if(s.donTrkThrust) nTrk_mix = s.ridge_eventSelection(mix.passesWW, mix.nParticle, mix.missP, mix.pt_wrtThr, mix.eta_wrtThr, mix.nTPC, mix.pwflag, mix.nref, mix.jtpt, mix.jteta,s.donTrkThrust);
                histNum_mix = s.histNum(nTrk_mix);
            }
            
            // Use the Thrust axis from the signal event instead of mixed event
            mix.TTheta=data.TTheta;
            mix.TPhi=data.TPhi;
            mix.update();
            
            if (i==selected)
            {
                std::cout <<"Error in Mixing!"<<std::endl;
                continue;
            }
            
            nBkgrndEvts[histNum] += 1;
            for ( Int_t j=0;j<data.nParticle;j++ )
            {
                if(!s.ridge_trackSelection(data.getPt(j),data.getEta(j),data.nTPC[j],data.pwflag[j],s.doThrust)) continue;
                Float_t angle;
                if (s.doTheta) angle = data.getTheta(j); else angle = data.getEta(j);
                Float_t phi = data.getPhi(j);
                
                // Background loop, calculate B correlation function from mixed event
                for ( Int_t k=0;k<mix.nParticle;k++ )
                {
                    if(!s.ridge_trackSelection(mix.getPt(k),mix.getEta(k),mix.nTPC[k],mix.pwflag[k],s.doThrust)) continue;
                    Float_t angle_mix;
                    if (s.doTheta) angle_mix = mix.getTheta(k); else angle_mix = mix.getEta(k);
                    Float_t phi_mix = mix.getPhi(k);
                    
                    bkgrnd2PC[histNum]->Fill(angle-angle_mix,dphi(phi,phi_mix),fillNumerator/(s.getDifferential())/nTrk);
                    bkgrnd2PC[histNum]->Fill(angle-angle_mix,dphi(phi_mix,phi),fillNumerator/(s.getDifferential())/nTrk);
                    bkgrnd2PC[histNum]->Fill(angle_mix-angle,dphi(phi,phi_mix),fillNumerator/(s.getDifferential())/nTrk);
                    bkgrnd2PC[histNum]->Fill(angle_mix-angle,dphi(phi_mix,phi),fillNumerator/(s.getDifferential())/nTrk);
                } //end of mixed event loop
            } // end of working event loop
            //std::cout<<i<<" "<<selected<<" "<<nTrk<<" "<<nTrk_mix<<" "<<std::endl;
            
        } // end of nMix loop
    } // end of loop over events
    
    // write the histograms
    output->cd();
    for(int i = 0; i<s.nMultBins; i++)
    {
        std::cout << "nSignalEvts "<< s.multBinsLow[i]<<" "<<s.multBinsHigh[i]<<" : "<<nSignalEvts[i] <<" " << std::endl;
        std::cout << "nBkgrndEvts "<< s.multBinsLow[i]<<" "<<s.multBinsHigh[i]<<" : "<<nBkgrndEvts[i] <<" " << std::endl;
        if(!s.doParallel) signal2PC[i]->Scale(1.0/(float)nSignalEvts[i]);
        if(!s.doParallel) bkgrnd2PC[i]->Scale(1.0/(float)nBkgrndEvts[i]);
        nEvtSigHist->Fill(i,nSignalEvts[i]);
 +      nEvtBkgHist->Fill(i,nBkgrndEvts[i]);
        calculateRatio(signal2PC[i],bkgrnd2PC[i],ratio2PC[i]);
        getLongRangeYield(s,ratio2PC[i],longRangeYield[i]);

        signal2PC[i]->Write();
        bkgrnd2PC[i]->Write();
        ratio2PC[i]->Write();
        longRangeYield[i]->Write();
        nEvtSigHist->Write();
        nEvtBkgHist->Write();
    }
    
    h_eff->Write();
    h_Aj->Write();
    h_Tphi->Write();
    h_Ttheta->Write();
    h_pt->Write();
    h_theta->Write();
    h_eta->Write();
    h_phi->Write();
    multiplicity->Write();
    
    // cleanup
    delete h_Aj;
    delete h_Tphi;
    delete h_Ttheta;
    delete h_pt;
    delete h_theta;
    delete h_eta;
    delete h_phi;
    delete h_eff;
    delete multiplicity;
    for(int i = 0; i<s.nMultBins; i++)
    {
        delete ratio2PC[i];
        delete bkgrnd2PC[i];
        delete signal2PC[i];
    }
    delete jt_mix;
    delete t_mix;
    delete jt;
    delete t;
    
    output->Close();
    delete output;
    
    return 0;
}

int main(int argc, char* argv[])
{
    if(argc < 3 || argc > 19)
    {
        std::cout << "Usage: ./ridge_check_parallel.exe <inFileName> <outFileName>" <<std::endl;
        return 1;
    }
    
    int retVal = 0;
    if(argc == 3) retVal += ridge_check_parallel(argv[1], argv[2]);
    return retVal;
}