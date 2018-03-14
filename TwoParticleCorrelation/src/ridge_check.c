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
#include "../../DataProcessing/include/smartJetName.h"

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

int ridge_check
    (
     const std::string inFileName, // input file
     std::string outFileName    // output file
    )
{
    // ROOT Global setting
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    // set up plots
    Selection s = Selection();
    if(inFileName.find("JobNum") != std::string::npos){s.doNTPC = false; s.dod0 = false; s.doz0 = false; s.doWW = false;} // for MC

    if(outFileName.find(".root") == std::string::npos) outFileName = outFileName + ".root";
    TFile * output = TFile::Open(outFileName.c_str(),"recreate");

    /// Initialize the histograms
    std::cout<<"Initializing histograms..."<<std::endl;
    static const Int_t nEnergyBins = s.nEnergyBins;
    static const Int_t nMultBins = s.nMultBins;
    static const Int_t nptBins = s.nptBins;
    static const Int_t netaBins = s.netaBins;

    
    TH2F * signal2PC[nEnergyBins][nMultBins][nptBins][netaBins];
    TH2F * bkgrnd2PC[nEnergyBins][nMultBins][nptBins][netaBins];
    TH2F * ratio2PC[nEnergyBins][nMultBins][nptBins][netaBins];
    TH1F * longRangeYield[nEnergyBins][nMultBins][nptBins][netaBins];
    float nSignalEvts[nMultBins] = {0};
    float nBkgrndEvts[nMultBins] = {0};
    
    TH1D * h_eff = new TH1D("h_eff","Selection Efficiency",100,0.0,1.0);;
    TH1D * h_phi = new TH1D("phi","phi",100,-TMath::Pi(),TMath::Pi());
    TH1D * h_eta = new TH1D("eta","eta",100,-5,5);
    TH1D * h_theta = new TH1D("theta","theta",100,0,TMath::Pi());
    TH1D * h_pt = new TH1D("pt","pt",100,0,10);
    TH1D * h_Ttheta = new TH1D("T_theta","T_theta",100,0,TMath::Pi());
    TH1D * h_Tphi = new TH1D("T_phi","T_phi",100,-TMath::Pi(),TMath::Pi());
    TH1D * h_Aj = new TH1D("h_Aj","h_Aj",50,0,0.5);
    TH1D * nEvtSigHist = new TH1D("nEvtSigHisto","nEvtSigHisto",10,0,10);  // this only depends on the event selection criteria so it is independent of the pt and eta binning
    TH1D * nEvtBkgHist = new TH1D("nEvtBkgHisto","nEvtBkgHisto",10,0,10);
    
    Float_t etaPlotRange = s.getEtaPlotRange();
    for(int e = 0; e<nEnergyBins; e++)
    {
        for(int i = 0; i<nMultBins; i++)
        {
            for(int j = 0; j<nptBins; j++)
            {
                for(int k = 0; k<netaBins; k++)
                {
                    if(s.doThrust)
                    {
                        signal2PC[e][i][j][k] = new TH2F(Form("signal2PC_%d_%d_%d_%d_%d",e,s.multBinsLow[i],s.multBinsHigh[i], j,k),";#Delta#eta;#Delta#Phi",s.dEtaBins,-etaPlotRange,etaPlotRange,s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
                        signal2PC[e][i][j][k]->Sumw2();
                        bkgrnd2PC[e][i][j][k] = new TH2F(Form("bkgrnd2PC_%d_%d_%d_%d_%d",e,s.multBinsLow[i],s.multBinsHigh[i], j,k),";#Delta#eta;#Delta#Phi",s.dEtaBins,-etaPlotRange,etaPlotRange,s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
                        bkgrnd2PC[e][i][j][k]->Sumw2();
                        //nEvtSigHist[i][j][k] = new TH1D(Form("nEvtSigHisto_%d_%d_%d_%d_%d",e,s.multBinsLow[i],s.multBinsHigh[i], j,k),"nEvtSigHisto",1,0,1);
                        //nEvtBkgHist[i][j][k] = new TH1D(Form("nEvtBkgHisto_%d_%d_%d_%d_%d",e,s.multBinsLow[i],s.multBinsHigh[i], j,k),"nEvtBkgHisto",1,0,1);
                        if(!s.doParallel)
                        {
                            ratio2PC[e][i][j][k] = new TH2F(Form("ratio2PC_%d_%d_%d_%d_%d",e,s.multBinsLow[i],s.multBinsHigh[i], j,k),";#Delta#eta;#Delta#Phi",s.dEtaBins,-etaPlotRange,etaPlotRange,s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
                            ratio2PC[e][i][j][k]->Sumw2();
                            longRangeYield[e][i][j][k] = new TH1F(Form("longRangeYield_%d_%d_%d_%d_%d",e,s.multBinsLow[i],s.multBinsHigh[i], j,k),";#Delta#phi;Y(#Delta#Phi)",s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
                            longRangeYield[e][i][j][k]->Sumw2();
                        } 
                    }
                    else if (s.doWTA)
                    {
                        signal2PC[e][i][j][k] = new TH2F(Form("signal2PC_%d_%d_%d_%d_%d",e,s.multBinsLow[i],s.multBinsHigh[i], j,k),";#Delta#eta;#Delta#Phi",s.dEtaBins,-etaPlotRange,etaPlotRange,s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
                        signal2PC[e][i][j][k]->Sumw2();
                        bkgrnd2PC[e][i][j][k] = new TH2F(Form("bkgrnd2PC_%d_%d_%d_%d_%d",e,s.multBinsLow[i],s.multBinsHigh[i], j,k),";#Delta#eta;#Delta#Phi",s.dEtaBins,-etaPlotRange,etaPlotRange,s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
                        bkgrnd2PC[e][i][j][k]->Sumw2();
                        //nEvtSigHist[i][j][k] = new TH1D(Form("nEvtSigHisto_%d_%d_%d_%d_%d",e,s.multBinsLow[i],s.multBinsHigh[i], s.ptBinsLow_wrtWTA[j],s.ptBinsHigh_wrtWTA[j],s.etaBinsLow_wrtWTA[k]),"nEvtSigHisto",1,0,1);
                        //nEvtBkgHist[i][j][k] = new TH1D(Form("nEvtBkgHisto_%d_%d_%d_%d_%d",e,s.multBinsLow[i],s.multBinsHigh[i], s.ptBinsLow_wrtWTA[j],s.ptBinsHigh_wrtWTA[j],s.etaBinsLow_wrtWTA[k]),"nEvtBkgHisto",1,0,1);
                        if(!s.doParallel)
                        {
                            ratio2PC[e][i][j][k] = new TH2F(Form("ratio2PC_%d_%d_%d_%d_%d",e,s.multBinsLow[i],s.multBinsHigh[i], j,k),";#Delta#eta;#Delta#Phi",s.dEtaBins,-etaPlotRange,etaPlotRange,s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
                            ratio2PC[e][i][j][k]->Sumw2();
                            longRangeYield[e][i][j][k] = new TH1F(Form("longRangeYield_%d_%d_%d_%d_%d",e,s.multBinsLow[i],s.multBinsHigh[i], j,k),";#Delta#phi;Y(#Delta#Phi)",s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
                            longRangeYield[e][i][j][k]->Sumw2();
                        }  
                    }
                    else /* Beam */
                    {
                        signal2PC[e][i][j][k] = new TH2F(Form("signal2PC_%d_%d_%d_%d_%d",e,s.multBinsLow[i],s.multBinsHigh[i], j,k),";#Delta#eta;#Delta#Phi",s.dEtaBins,-etaPlotRange,etaPlotRange,s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
                        signal2PC[e][i][j][k]->Sumw2();
                        bkgrnd2PC[e][i][j][k] = new TH2F(Form("bkgrnd2PC_%d_%d_%d_%d_%d",e,s.multBinsLow[i],s.multBinsHigh[i], j,k),";#Delta#eta;#Delta#Phi",s.dEtaBins,-etaPlotRange,etaPlotRange,s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
                        bkgrnd2PC[e][i][j][k]->Sumw2();
                        //nEvtSigHist[i][j][k] = new TH1D(Form("nEvtSigHisto_%d_%d_%d_%d_%d",e,s.multBinsLow[i],s.multBinsHigh[i], s.ptBinsLow_wrtWTA[j],s.ptBinsHigh_wrtWTA[j],s.etaBinsLow_wrtWTA[k]),"nEvtSigHisto",1,0,1);
                        //nEvtBkgHist[i][j][k] = new TH1D(Form("nEvtBkgHisto_%d_%d_%d_%d_%d",e,s.multBinsLow[i],s.multBinsHigh[i], s.ptBinsLow_wrtWTA[j],s.ptBinsHigh_wrtWTA[j],s.etaBinsLow_wrtWTA[k]),"nEvtBkgHisto",1,0,1);
                        if(!s.doParallel)
                        {
                            ratio2PC[e][i][j][k] = new TH2F(Form("ratio2PC_%d_%d_%d_%d_%d",e,s.multBinsLow[i],s.multBinsHigh[i], j,k),";#Delta#eta;#Delta#Phi",s.dEtaBins,-etaPlotRange,etaPlotRange,s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
                            ratio2PC[e][i][j][k]->Sumw2();
                            longRangeYield[e][i][j][k] = new TH1F(Form("longRangeYield_%d_%d_%d_%d_%d",e,s.multBinsLow[i],s.multBinsHigh[i], j,k),";#Delta#phi;Y(#Delta#Phi)",s.dPhiBins,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
                            longRangeYield[e][i][j][k]->Sumw2();
                        }  
                    }   
                }
            }   
        }
    }

    TH1F * multiplicity = new TH1F("multiplicity",";nTrk;nEvents",200,0,200);
    
    /// Initialize the trees for use
    // std::cout<<"Initializing trees for use..."<<std::endl;
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
    if (s.doWTA) sideTree = "BoostedWTAR8Evt";
    
    TChain * t = new TChain("t");       t->Add(inFileName.c_str());
    TChain * side_t = new TChain(sideTree.c_str());       side_t->Add(inFileName.c_str());
    TChain * jt = new TChain(jtTreeName.c_str());       jt->Add(inFileName.c_str());

    TPCNtupleData data(s.doBelle, s.doThrust, s.doWTA);      setupTPCTree(t,side_t,jt,data);       data.setTPCTreeStatus(t,side_t);
    
    TChain * t_mix = new TChain("t");       t_mix->Add(inFileName.c_str());
    TChain * side_t_mix = new TChain(sideTree.c_str());       side_t_mix->Add(inFileName.c_str());
    TChain * jt_mix = new TChain(jtTreeName.c_str());       jt_mix->Add(inFileName.c_str());
    
    TPCNtupleData mix(s.doBelle, s.doThrust, s.doWTA);       setupTPCTree(t_mix,side_t_mix,jt_mix,mix);        mix.setTPCTreeStatus(t_mix,side_t_mix);
    
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
        
        // find the event energy histogram(s)
        std::vector<Int_t> histE;
        if(inFileName.find("JobNum") != std::string::npos) histE.push_back(0); // the PYTHIA8 is e+e-→Z on Z pole ~ 91.5 GeV
        else histE = s.histEnergy(data.Energy);
        if(histE.size() == 0) { std::cout<<"event energy does not fit in any of the specified ranges...skipping event"<<std::endl; continue;}

        // nTrk calculation
        Int_t nTrk = 0;
        if(!s.donTrkThrust) nTrk = s.ridge_eventSelection(data.passesWW, data.missP, data.nParticle, data.nref, data.jtpt, data.jteta, data.STheta, data.mass, data.nTPC, data.theta, data.pmag, data.d0, data.z0, data.pwflag);
        if(s.donTrkThrust) nTrk = s.ridge_eventSelection(data.passesWW, data.missP, data.nParticle, data.nref, data.jtpt, data.jteta, data.STheta, data.mass, data.nTPC, data.theta_wrtThr, data.pmag, data.d0, data.z0, data.pwflag);
        if( nTrk < 0) continue;
        //std::cout<<data.RunNo<<","<<data.EventNo<<std::endl;
        
        // find the event nTrk histogram(s)
        std::vector<Int_t> histNtrk = s.histNtrk(nTrk);
	    if(histNtrk.size() == 0) { std::cout<<"calculated nTrk does not fit in any of the specified ranges...skipping event"<<std::endl; continue;}

        h_eff->Fill(nTrk/data.nParticle);
        h_Aj->Fill(s.fillAj);
        multiplicity->Fill(nTrk);
        
        // increment the number of signal events bin for each nTrk histogram going to be filled
        for(unsigned int nI = 0; nI < histNtrk.size(); ++nI) nSignalEvts[histNtrk.at(nI)] += 1;
        
        h_Ttheta->Fill(data.TTheta);
        h_Tphi->Fill(data.TPhi);

        Float_t fillNumerator = 1.0;
        if (s.doPP) fillNumerator = data.pthatWeight;
        /****************************************/
        // S calculation using multiplicity cut //
        /****************************************/
        for ( Int_t j=0;j<data.nParticle;j++ )
        {
            if(!s.ridge_trackSelection(data.nTPC[j],data.getTheta(j),data.getPmag(j), data.d0[j], data.z0[j], data.pwflag[j])) continue;
            h_phi->Fill(data.getPhi(j));
            h_eta->Fill(data.getEta(j));
            h_theta->Fill(data.getTheta(j));
            h_pt->Fill(data.getPt(j));

            // Decide which pt and eta range to fill
            std::vector<Int_t> histPt1 = s.histPt(data.getPt(j));
            std::vector<Int_t> histEta1 = s.histEta(data.getEta(j));
            if(histPt1.size() == 0 || histEta1.size() == 0) continue; 
            //std::cout<<"1st particle "<<histEta1.size()<<std::endl;
            Float_t angle1;
            if (s.getThetaAngle) angle1 = data.getTheta(j); else angle1 = data.getEta(j);
            Float_t phi1 = data.getPhi(j);
            //std::cout<<"eta for first particle "<<angle1<<std::endl;
            // Signal loop, calculate S correlation function
            for ( Int_t k=j+1;k<data.nParticle;k++ )
            {
                if(!s.ridge_trackSelection(data.nTPC[k],data.getTheta(k),data.getPmag(k), data.d0[k], data.z0[k], data.pwflag[k])) continue;
                // Check if the second particle is in the same range of pt and eta 

                std::vector<Int_t> histPt2 = s.histPt(data.getPt(k));
                std::vector<Int_t> histEta2 = s.histEta(data.getEta(k));
                //std::cout<<histEta2.size()<<std::endl;
                //initial check that histPt2(Eta2) are not empty
                if(histPt2.size() == 0 || histEta2.size() == 0) continue; 
                // Loop over elements of histPt1(Eta1) and check if the element is present in histPt2(Eta2). if yes then leave the element, if no then remove it. 
                // This ensures we fill the correct histograms since both particles should be within the same pt range 
                for(unsigned int pI = 0; pI <histPt1.size(); pI++) if(std::find(histPt2.begin(), histPt2.end(), histPt1.at(pI)) == histPt2.end()) histPt1.erase(histPt1.begin() + pI);
                for(unsigned int eI = 0; eI <histEta1.size(); eI++) if(std::find(histEta2.begin(), histEta2.end(), histEta1.at(eI)) == histEta2.end()) histEta1.erase(histEta1.begin() + eI);
                //check after the loops that histPt1(Eta1) are not empty
                if(histPt1.size() == 0 || histEta1.size() == 0) continue; 
                
                Float_t angle2;
                if (s.getThetaAngle) angle2 = data.getTheta(k); else angle2 = data.getEta(k);
                Float_t phi2 = data.getPhi(k);
                
                for(unsigned int eI = 0; eI< histE.size(); eI++)
                {
                    for(unsigned int nI = 0; nI< histNtrk.size(); nI++)
                    {
                        for(unsigned int pI = 0; pI< histPt1.size(); pI++)
                        {
                            for(unsigned int etI = 0; etI< histEta1.size(); etI++)
                            {
                                //std::cout<<"histEta "<<histEta1.at(etI)<<" eta1 = "<<angle1<<" eta2 = "<<angle2<<std::endl;
                                signal2PC[histE.at(eI)][histNtrk.at(nI)][histPt1.at(pI)][histEta1.at(etI)]->Fill(angle1-angle2,dphi(phi1,phi2),fillNumerator/(s.getDifferential())/nTrk);
                                signal2PC[histE.at(eI)][histNtrk.at(nI)][histPt1.at(pI)][histEta1.at(etI)]->Fill(angle1-angle2,dphi(phi2,phi1),fillNumerator/(s.getDifferential())/nTrk);
                                signal2PC[histE.at(eI)][histNtrk.at(nI)][histPt1.at(pI)][histEta1.at(etI)]->Fill(angle2-angle1,dphi(phi1,phi2),fillNumerator/(s.getDifferential())/nTrk);
                                signal2PC[histE.at(eI)][histNtrk.at(nI)][histPt1.at(pI)][histEta1.at(etI)]->Fill(angle2-angle1,dphi(phi2,phi1),fillNumerator/(s.getDifferential())/nTrk);
                            }
                        }
                    }
                }
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
            if(!s.donTrkThrust) nTrk_mix = s.ridge_eventSelection(mix.passesWW, mix.missP, mix.nParticle, mix.nref, mix.jtpt, mix.jteta, mix.STheta, mix.mass, mix.nTPC, mix.theta, mix.pmag, mix.d0, mix.z0, mix.pwflag);
            if(s.donTrkThrust) nTrk_mix = s.ridge_eventSelection(mix.passesWW, mix.missP, mix.nParticle, mix.nref, mix.jtpt, mix.jteta, mix.STheta, mix.mass, mix.nTPC, mix.theta, mix.pmag, mix.d0, mix.z0, mix.pwflag);
            // find the event nTrk histogram(s)
            std::vector<Int_t> histNtrk_mix = s.histNtrk(nTrk_mix);
            // loop over the background event nTrk histos and determine if the signal event has the same histos. If yes then leave the entry. If no then remove the entry from histNtrk_mix. 
            for(unsigned int nI = 0; nI <histNtrk_mix.size(); nI++) if(std::find(histNtrk.begin(), histNtrk.end(), histNtrk_mix.at(nI)) == histNtrk.end()) histNtrk_mix.erase(histNtrk_mix.begin() + nI);

            // We continue to search until we find an event that passes the mixed event selection criteria and has a non-null histNtrk_mix intersection with histNtrk.
            // s.isMixedEvent returns false if the mixed event fails the criteria in which case we want to enter the while loop.
            // If the size of histNtrk_mix is 0 then the mixed event is a different nTrk range than the signal event and we do not want mix in which case we want to enter the while loop.
            while (histNtrk_mix.size() == 0 || !s.isMixedEvent(nTrk, nTrk_mix, data.jteta[0], mix.jteta[0], data.getTTheta(), mix.getTTheta(), data.getTPhi(), mix.getTPhi()))
            {
                selected++;
                if (selected >= t->GetEntries()) selected = 0;
                t_mix->GetEntry(selected);
                side_t_mix->GetEntry(selected);
                jt_mix->GetEntry(selected);
                
                if(!s.donTrkThrust) nTrk_mix = s.ridge_eventSelection(mix.passesWW, mix.missP, mix.nParticle, mix.nref, mix.jtpt, mix.jteta, mix.STheta, mix.mass, mix.nTPC, mix.theta, mix.pmag, mix.d0, mix.z0, mix.pwflag);
                if(s.donTrkThrust) nTrk_mix = s.ridge_eventSelection(mix.passesWW, mix.missP, mix.nParticle, mix.nref, mix.jtpt, mix.jteta, mix.STheta, mix.mass, mix.nTPC, mix.theta, mix.pmag, mix.d0, mix.z0, mix.pwflag);
                // find the event nTrk histogram(s)
                histNtrk_mix = s.histNtrk(nTrk_mix);
                // loop over the background event nTrk histos and determine if the signal event has the same histos. If yes then leave the entry. If no then remove the entry from histNtrk_mix. 
                for(unsigned int nI = 0; nI <histNtrk_mix.size(); nI++) if(std::find(histNtrk.begin(), histNtrk.end(), histNtrk_mix.at(nI)) == histNtrk.end()) histNtrk_mix.erase(histNtrk_mix.begin() + nI);
            }
            
            // Use the Thrust axis from the signal event instead of mixed event
            //mix.TTheta=data.TTheta;
            //mix.TPhi=data.TPhi;
            mix.update();
            
            if (i==selected)
            {
                std::cout <<"Error in Mixing!"<<std::endl;
                continue;
            }
            
            // increment the number of background events bin for each nTrk histogram going to be filled (note that this should give the same result as the signal events increment)
            for(unsigned int nI = 0; nI < histNtrk.size(); ++nI) nBkgrndEvts[histNtrk.at(nI)] += 1;

            for ( Int_t j=0;j<data.nParticle;j++ )
            {
                // decide if valid track
                if(!s.ridge_trackSelection(data.nTPC[j],data.getTheta(j),data.getPmag(j), data.d0[j], data.z0[j], data.pwflag[j])) continue;

                // Decide which pt and eta range to fill
                std::vector<Int_t> histPt_bkg1 = s.histPt(data.getPt(j));
                std::vector<Int_t> histEta_bkg1 = s.histEta(data.getEta(j));
                if(histPt_bkg1.size() == 0 || histEta_bkg1.size() == 0) continue; 

                // load the angles to be used
                Float_t angle;
                if (s.getThetaAngle) angle = data.getTheta(j); else angle = data.getEta(j);
                Float_t phi = data.getPhi(j);

                // Background loop, calculate B correlation function from mixed event
                for ( Int_t k=0;k<mix.nParticle;k++ )
                {
                    // decide if valid track
                    if(!s.ridge_trackSelection(mix.nTPC[k],mix.getTheta(k),mix.getPmag(k), mix.d0[k], mix.z0[k], mix.pwflag[k])) continue;
                    
                    // Check if the second particle is in the same range of pt and eta    
                    std::vector<Int_t> histPt_bkg2 = s.histPt(mix.getPt(k));
                    std::vector<Int_t> histEta_bkg2 = s.histEta(mix.getEta(k));
                    //initial check that histPt2(Eta2) are not empty
                    if(histPt_bkg2.size() == 0 || histEta_bkg2.size() == 0) continue; 
                    // Loop over elements of histPt1(Eta1) and check if the element is present in histPt2(Eta2). if yes then leave the element, if no then remove it. 
                    // This ensures we fill the correct histograms since both particles should be within the same pt range 
                    for(unsigned int pI = 0; pI <histPt_bkg1.size(); pI++) if(std::find(histPt_bkg2.begin(), histPt_bkg2.end(), histPt_bkg1.at(pI)) == histPt_bkg2.end()) histPt_bkg1.erase(histPt_bkg1.begin() + pI);
                    for(unsigned int eI = 0; eI <histEta_bkg1.size(); eI++) if(std::find(histEta_bkg2.begin(), histEta_bkg2.end(), histEta_bkg1.at(eI)) == histEta_bkg2.end()) histEta_bkg1.erase(histEta_bkg1.begin() + eI);
                    //check after the loops that histPt1(Eta1) are not empty
                    if(histPt_bkg1.size() == 0 || histEta_bkg1.size() == 0) continue; 

                    // load the angles to be used 
                    Float_t angle_mix;
                    if (s.getThetaAngle) angle_mix = mix.getTheta(k); else angle_mix = mix.getEta(k);
                    Float_t phi_mix = mix.getPhi(k);

                    for(unsigned int eI = 0; eI< histE.size(); eI++)
                    {
                        for(unsigned int nI = 0; nI< histNtrk_mix.size(); nI++) // histNtrk_mix contains the intersection of histNtrk_mix and histNtrk
                        {
                            for(unsigned int pI = 0; pI< histPt_bkg1.size(); pI++)
                            {
                                for(unsigned int etI = 0; etI< histEta_bkg1.size(); etI++)
                                {
                                    bkgrnd2PC[histE.at(eI)][histNtrk.at(nI)][histPt_bkg1.at(pI)][histEta_bkg1.at(etI)]->Fill(angle-angle_mix,dphi(phi,phi_mix),fillNumerator/(s.getDifferential())/nTrk);
                                    bkgrnd2PC[histE.at(eI)][histNtrk.at(nI)][histPt_bkg1.at(pI)][histEta_bkg1.at(etI)]->Fill(angle-angle_mix,dphi(phi_mix,phi),fillNumerator/(s.getDifferential())/nTrk);
                                    bkgrnd2PC[histE.at(eI)][histNtrk.at(nI)][histPt_bkg1.at(pI)][histEta_bkg1.at(etI)]->Fill(angle_mix-angle,dphi(phi,phi_mix),fillNumerator/(s.getDifferential())/nTrk);
                                    bkgrnd2PC[histE.at(eI)][histNtrk.at(nI)][histPt_bkg1.at(pI)][histEta_bkg1.at(etI)]->Fill(angle_mix-angle,dphi(phi_mix,phi),fillNumerator/(s.getDifferential())/nTrk);
                                }
                            }
                        }
                    }

                } //end of mixed event loop
            } // end of working event loop
            //std::cout<<i<<" "<<selected<<" "<<nTrk<<" "<<nTrk_mix<<" "<<std::endl;
        } // end of nMix loop
    } // end of loop over events
    
    output->cd();
    // write the histograms
    for(int e = 0; e<nEnergyBins;e++)
    {
        for(int i = 0; i<nMultBins; i++)
        {
            for(int j = 0; j<nptBins; j++)
            {
                for(int k = 0; k<netaBins; k++)
                {
                    signal2PC[e][i][j][k]->Write();
                    bkgrnd2PC[e][i][j][k]->Write();
                    if(!s.doParallel)
                    {
                        signal2PC[e][i][j][k]->Scale(1.0/(float)nSignalEvts[i]);
                        bkgrnd2PC[e][i][j][k]->Scale(1.0/(float)nBkgrndEvts[i]);
                        calculateRatio(signal2PC[e][i][j][k],bkgrnd2PC[e][i][j][k],ratio2PC[e][i][j][k]);
                        getLongRangeYield(s,ratio2PC[e][i][j][k],longRangeYield[e][i][j][k]);
                    }
                }
            }
            nEvtSigHist->Fill(i,nSignalEvts[i]);
            nEvtBkgHist->Fill(i,nBkgrndEvts[i]);
        }
    }

    nEvtSigHist->Write();
    nEvtBkgHist->Write();
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

    for(int e = 0; e<nEnergyBins; e++)
    {
        for(int i = 0; i<nMultBins; i++)
        {
            for(int j = 0; j<nptBins; j++)
            {
                for(int k = 0; k<netaBins; k++)
                {
                    delete signal2PC[e][i][j][k];
                    delete bkgrnd2PC[e][i][j][k];
                    if(!s.doParallel)
                    {
                       delete ratio2PC[e][i][j][k];
                       delete longRangeYield[e][i][j][k]; 
                    } 
                }
            }
        }
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
        std::cout << "Usage: ./ridge_check.exe <inFileName> <outFileName>" <<std::endl;
        return 1;
    }
    
    int retVal = 0;
    if(argc == 3) retVal += ridge_check(argv[1], argv[2]);
    return retVal;
}
