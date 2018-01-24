//c and c++ dependencies
#include <sstream>
#include <string>
#include <iostream>
#include <stdlib.h>

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

//local headers
#include "../include/fourier.h"
#include "../include/TPCNtupleData.h"
#include "../include/Selection.h"

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
     const std::string outFileName    // output file
    )
{
    // ROOT Global setting
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    // set up plots
    Selection s = Selection();
    TFile * output = TFile::Open("RidgeCheck_Output.root","recreate");
    TH2F * signal2PC[s.nMultBins];
    TH2F * bkgrnd2PC[s.nMultBins];
    TH2F * ratio2PC[s.nMultBins];
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
    
    for(int i = 0; i<s.nMultBins; i++)
    {
        signal2PC[i] = new TH2F(Form("signal2PC_%d_%d",s.multBinsLow[i],s.multBinsHigh[i]),";#Delta#eta;#Delta#Phi",s.dEtaBins,-2*s.etaPlotRange,2*s.etaPlotRange,s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
        bkgrnd2PC[i] = new TH2F(Form("bkgrnd2PC_%d_%d",s.multBinsLow[i],s.multBinsHigh[i]),";#Delta#eta;#Delta#Phi",s.dEtaBins,-2*s.etaPlotRange,2*s.etaPlotRange,s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
        ratio2PC[i] = new TH2F(Form("ratio2PC_%d_%d",s.multBinsLow[i],s.multBinsHigh[i]),";#Delta#eta;#Delta#Phi",s.dEtaBins,-2*s.etaPlotRange,2*s.etaPlotRange,s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
    }
    TH1F * multiplicity = new TH1F("multiplicity",";nTrk;nEvents",200,0,200);
    
    // files and variables for input
    TChain * t = new TChain("t");       t->Add(inFileName.c_str());
    TChain * jt = new TChain("ak4ESchemeJetTree");       jt->Add(inFileName.c_str());
    TPCNtupleData data(s.doBelle, s.doThrust);      setupTPCTree(t,jt,data);       data.setTPCTreeStatus(t);
    
    TChain * t_mix = new TChain("t");       t_mix->Add(inFileName.c_str());
    TChain * jt_mix = new TChain("ak4ESchemeJetTree");       jt_mix->Add(inFileName.c_str());
    TPCNtupleData mix(s.doBelle, s.doThrust);       setupTPCTree(t_mix,jt_mix,mix);        mix.setTPCTreeStatus(t_mix);
    
    // analysis
    Int_t nevent = (Int_t)t->GetEntries();
    
    /****************************************/
    // Main Event Loop
    /****************************************/
    for (Int_t i=0;i<nevent;i++) {
        t->GetEntry(i);
        jt->GetEntry(i);
        data.update();
        if (i%10000==0) std::cout <<i<<"/"<<nevent<<std::endl;
        
        // event cut
        Int_t N = s.ridge_eventSelection(data.passesWW, data.nParticle, data.missP, data.pt, data.eta, data.nTPC, data.pwflag, data.nref, data.jtpt, data.jteta);
        if( N < 0) continue;
        Int_t histNum = s.histNum(N);
        
        h_eff->Fill(N/data.nParticle);
        h_Aj->Fill(s.fillAj);
        multiplicity->Fill(N);
        nSignalEvts[histNum] += 1;
        
        h_Ttheta->Fill(data.TTheta);
        h_Tphi->Fill(data.TPhi);
        /****************************************/
        // S calculation using multiplicity cut //
        /****************************************/
        for ( Int_t j=0;j<data.nParticle;j++ )
        {
            if(!s.ridge_trackSelection(data.getPt(j),data.getEta(j),data.nTPC[j],data.pwflag[j])) continue;
            
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
                if(!s.ridge_trackSelection(data.getPt(k),data.getEta(k),data.nTPC[k],data.pwflag[k])) continue;
                
                Float_t angle2;
                if (s.doTheta) angle2 = data.getTheta(k); else angle2 = data.getEta(k);
                Float_t phi2 = data.getPhi(k);
                
                signal2PC[histNum]->Fill(angle1-angle2,dphi(phi1,phi2),1./(s.differential)/N);
                signal2PC[histNum]->Fill(angle1-angle2,dphi(phi2,phi1),1./(s.differential)/N);
                signal2PC[histNum]->Fill(angle2-angle1,dphi(phi1,phi2),1./(s.differential)/N);
                signal2PC[histNum]->Fill(angle2-angle1,dphi(phi2,phi1),1./(s.differential)/N);
            }
        }
        
        /****************************************/
        // B calculation using multiplicity cut //
        /****************************************/
        for (Int_t nMix = 0; nMix<s.bkgrd_runs; nMix++)
        {
            Int_t selected=i+1;
            t_mix->GetEntry(selected);
            jt_mix->GetEntry(selected);
            
            Int_t N_mix = s.ridge_eventSelection(mix.passesWW, mix.nParticle, mix.missP, mix.pt, mix.eta, mix.nTPC, mix.pwflag, mix.nref, mix.jtpt, mix.jteta);
            
            // Select a mixed event
            Int_t flag=0;
            while (N_mix < 0 && !s.mixedEvent(N, N_mix, data.jteta[0], mix.jteta[0]))
            {
                selected++;
                if (selected > nevent) break;
                t_mix->GetEntry(selected);
                jt_mix->GetEntry(selected);
                N_mix = s.ridge_eventSelection(mix.passesWW, mix.nParticle, mix.missP, mix.pt, mix.eta, mix.nTPC, mix.pwflag, mix.nref, mix.jtpt, mix.jteta);
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
                if(!s.ridge_trackSelection(data.getPt(j),data.getEta(j),data.nTPC[j],data.pwflag[j])) continue;
                
                Float_t angle;
                if (s.doTheta) angle = data.getTheta(j); else angle = data.getEta(j);
                Float_t phi = data.getPhi(j);
                
                // Background loop, calculate B correlation function from mixed event
                for ( Int_t k=0;k<mix.nParticle;k++ )
                {
                    if(!s.ridge_trackSelection(mix.getPt(k),mix.getEta(k),mix.nTPC[k],mix.pwflag[k])) continue;
                    
                    Float_t angle_mix;
                    if (s.doTheta) angle_mix = data.getTheta(j); else angle_mix = mix.getEta(k);
                    Float_t phi_mix = mix.getPhi(j);
                    
                    bkgrnd2PC[histNum]->Fill(angle-angle_mix,dphi(phi,phi_mix),1./(s.differential)/N);
                    bkgrnd2PC[histNum]->Fill(angle-angle_mix,dphi(phi_mix,phi),1./(s.differential)/N);
                    bkgrnd2PC[histNum]->Fill(angle_mix-angle,dphi(phi,phi_mix),1./(s.differential)/N);
                    bkgrnd2PC[histNum]->Fill(angle_mix-angle,dphi(phi_mix,phi),1./(s.differential)/N);
                } //end of mixed event loop
            } // end of working event loop
        } // end of nMix loop
    } // end of loop over events
    
    // write the histograms
    output->cd();
    for(int i = 0; i<s.nMultBins; i++)
    {
        std::cout << "nSignalEvts "<< s.multBinsLow[i]<<" "<<s.multBinsHigh[i]<<" : "<<nSignalEvts[i] <<" " << std::endl;
        std::cout << "nBkgrndEvts "<< s.multBinsLow[i]<<" "<<s.multBinsHigh[i]<<" : "<<nBkgrndEvts[i] <<" " << std::endl;
        signal2PC[i]->Scale(1./ (float)nSignalEvts[i]);
        bkgrnd2PC[i]->Scale(1./(float)nBkgrndEvts[i]);
        calculateRatio(signal2PC[i],bkgrnd2PC[i],ratio2PC[i]);
        signal2PC[i]->Write();
        bkgrnd2PC[i]->Write();
        ratio2PC[i]->Write();
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
