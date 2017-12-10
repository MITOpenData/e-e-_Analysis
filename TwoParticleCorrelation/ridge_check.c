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
#include <iostream>
#include <stdlib.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TFormula.h>
#include <sstream>
#include <string>
#include <TNtuple.h>

#include "fourier.h"
#include "TPCNtupleData.h"

/**************************************************************************************/
// Two particle correlation analysis
//
// ridge_check.c
//
// Yen-Jie Lee, Gian Michele Innocenti, Anthony Badea, Bibek Pandit, Tzu-An Sheng
//
// LOG
//     2017/08/15 Added Thurst Axis based correlation function to the TPCNtupleData 
//                Class (based on Austin's code). Added filename to the class.
//     2017/08/16 Cleanup and change the plotting style. Add Fourier Fit.
//     2017/08/17 Add Thetabased analysis.
//
/**************************************************************************************/

using namespace std;

/**************************************************************************************/
// Main Analysis Routine
/**************************************************************************************/

void analysis(TString filename = "/data/flowex/CMSsample/TPCNtuple_MinBias_TuneCUETP8M1_5p02TeV-pythia8-HINppWinter16DR-NoPU_75X_mcRun2_asymptotic_ppAt5TeV_forest_v2_track.root",
              TString savename = "TPCNtuple_MinBias_TuneCUETP8M1_5p02TeV-pythia8-HINppWinter16DR-NoPU_75X_mcRun2_asymptotic_ppAt5TeV_forest_v2_track", // used for saving the files
              Int_t isBelle      = 0,		// BELLE analysis = 1, CMS/ALEPH analysis = 0
              Int_t isThrust     = 0, 		// Leading Jet Axis analysis = 2, Thurst Axis analysis = 1, Beam Axis analysis = 0
              Int_t isTheta      = 0, 		// Use Theta angle = 1, Use Eta = 0
              Int_t isGen       = 0,         // Use isGen = 1 if gen level so no cut on charged particles, 0 if not gen level
	      Int_t maxevt       = 1000000,	// Max number of events to be processed, 0 = all events
	      Int_t mult_low     = 0,		// Lower cut on the event multiplicity
	      Int_t mult_high    = 100,		// Upper cut on the event multiplicity
	      Int_t nbin         = 20,		// Number of bins in the correlation function
	      bool verbose     = 0,		// Verbose mode
	      Int_t num_runs     = 5,		// 
              double ptMin     = 0,           // min pT of the particles used for correlation function
              double ptMax     = 4.0,             // max pT of the particles used for correlation function
              double detaRange = 3.2,             // deta window of the correlation function
              Float_t ptMinForN = 0.4,          // pT min for the N_{trk}^{offline} calculation (used for event classification) 
              Float_t ptMaxForN = 100,          // pT max for the N_{trk}^{offline} calculation (used for event classification)
              Float_t etaCutForN = 2.4          // eta window for the N_{trk}^{offline} calculation (used for event classification)
	     ) {
    // ROOT Global setting
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    TChain *t1 = new TChain("t");
    t1->Add(filename);

    TChain *t2 = new TChain("ak4JetTree");
    t2->Add(filename);

    /* TString friendname =  "qqmc-e07-00_flavor.root"; */
    /* TFile *ff = new TFile(friendname); */
    /* TTree *tf = (TTree*)ff->Get("t"); */
    /* t1->AddFriend(tf); */

    TPCNtupleData data(isBelle, isThrust);
    setupTPCTree(t1,t2,data);
    data.setTPCTreeStatus(t1);
    
    // File for event mixing, use the same file for the moment
    /* TFile *f_mix = new TFile(filename.Data()); */
    /* TTree *t1_mix = (TTree*)f_mix->Get("t"); */
    TChain *t1_mix = new TChain("t");
    t1_mix->Add(filename);
             
    // Not necessary
    TChain *t2_mix = new TChain("ak4JetTree");
    t2_mix->Add(filename);
    
    TPCNtupleData mix(isBelle, isThrust);
    setupTPCTree(t1_mix,t2_mix,mix);
    mix.setTPCTreeStatus(t1_mix);
    
    //TNtuple *nt = new TNtuple("nt","","pEta:pTheta:pPhi:theta:phi:TTheta:TPhi");
    
    // Define 2D histograms
    Float_t dthetaRange = PI;
    Float_t normalization = detaRange*2/nbin*2*3.14159/nbin;
    TH2F *h_2D, *h_2Dmix, *h_ratio;

    if (isTheta){
      h_2D    = CFTH2F ( "h_2D",    "#theta-#phi of all particles;#Delta#theta;#Delta#phi;#frac{1}{N_{trig}}#frac{d^{2}N^{pair}}{d#Delta#thetad#Delta#phi} ",nbin, -detaRange, detaRange,nbin, -PI/2., PI*1.5);
      h_2Dmix = CFTH2F ( "h_2Dmix", "#theta-#phi of all particles;#Delta#theta;#Delta#phi;#frac{1}{N_{trig}}#frac{d^{2}N^{pair}}{d#Delta#thetad#Delta#phi} ",nbin, -detaRange, detaRange,nbin, -PI/2., PI*1.5);
      h_ratio = CFTH2F ( "h_ratio", "#theta-#phi of all particles;#Delta#theta;#Delta#phi;#frac{1}{N_{trig}}#frac{d^{2}N^{pair}}{d#Delta#thetad#Delta#phi} ",nbin, -detaRange, detaRange,nbin, -PI/2., PI*1.5);
    } else {
      h_2D    = CFTH2F ( "h_2D",    "#eta-#phi of all particles;#Delta#eta;#Delta#phi;#frac{1}{N_{trig}}#frac{d^{2}N^{pair}}{d#Delta#etad#Delta#phi} ",nbin, -detaRange, detaRange,nbin, -PI/2., PI*1.5);
      h_2Dmix = CFTH2F ( "h_2Dmix", "#eta-#phi of all particles;#Delta#eta;#Delta#phi;#frac{1}{N_{trig}}#frac{d^{2}N^{pair}}{d#Delta#etad#Delta#phi} ",nbin, -detaRange, detaRange,nbin, -PI/2., PI*1.5);
      h_ratio = CFTH2F ( "h_ratio", "#eta-#phi of all particles;#Delta#eta;#Delta#phi;#frac{1}{N_{trig}}#frac{d^{2}N^{pair}}{d#Delta#etad#Delta#phi} ",nbin, -detaRange, detaRange,nbin, -PI/2., PI*1.5);
    }
    TH1D *h_mult_dist = new TH1D("h_multi","multiplicity distribution of charged particles;N_{trk}^{offline};Count",60,0,60);
    TH1D *h_phi = new TH1D("h_phi","Phi distribution",60,-PI,PI);
    
    // all entries and fill the histograms
    Int_t nevent = (Int_t)t1->GetEntries();
    
    Int_t nevent_process = nevent;
    //if( maxevt>0 && maxevt<nevent ) nevent_process = maxevt;
    
    Float_t averageN=0;
    Int_t nEventProcessed=0;
    Int_t nEventInMultBin=0;
    
    /* Int_t flavor; */
    /* t1->SetBranchAddress("flavor", &flavor); */

    // check number of fills
    Int_t fill_count=0;
    Int_t fill_count_mix=0;
    /****************************************/
    // Main Event Loop
    /****************************************/
    for (Int_t i=0;i<nevent_process;i++) {
       t1->GetEntry(i);
       t2->GetEntry(i);
       data.update();
       if (i%10000==0) cout <<i<<"/"<<nevent_process<<endl;
       if (verbose) cout<<"nparticles="<<data.nParticle<<endl;
       
       if (isThrust==1&&fabs(data.TTheta-3.14159/2.)>0.8) continue;
       if (isThrust==2&&acos(cos(data.jtphi[0]-data.jtphi[1]))<3.14159*2/3.) continue;
       
       /* if (flavor == 0) {break;} */
       /* if (flavor == 1) {continue;} */

       Int_t selected=i+1;  //questo diventa il numero di evento estratto nel file +1
       Int_t flag=0;	  //definisco una flag a zero
       
       // Yen-Jie: we would like to use a definition similar to CMS publication
       // all particles with pT > 0.4 GeV/c for the calculation of event multiplicity N        

       Int_t N=0;	// N_{trk}^{offline}
       Int_t N_TP=0;	// Number of particles used in the two particle correlation function calculation, can be different from N for event classification

       // calculate the number of tracks in the passing selection
       for ( Int_t j=0;j<data.nParticle;j++ ) {
           Float_t pt1 = data.getPt(j);
           // if isGen = 0 then we want to cut on charged particles so !0->1 then isCharged matter
           if (!isGen&&!data.isChargedHadron(j)) continue;
	   if (fabs(data.getEta(j))>=etaCutForN) continue;
	   if (pt1>ptMinForN&&pt1<ptMaxForN) N++;
           if (pt1>ptMin&&pt1<ptMax) N_TP++;
	   h_phi->Fill(data.getPhi(j));
	   //nt->Fill(data.eta[j],data.theta[j],data.phi[j],data.getTheta(j),data.getPhi(j),data.TTheta,data.TPhi);
       }
       h_mult_dist->Fill(N);
       averageN+=N;
       nEventProcessed++;
       
       if (N<mult_low||N>=mult_high||N_TP<1) continue;
       nEventInMultBin++;
       
       /****************************************/
       // S calculation using multiplicity cut //
       /****************************************/
       for ( Int_t j=0;j<data.nParticle;j++ ) {
           Float_t angle1;
	   if (isTheta) angle1 = data.getTheta(j); else angle1 = data.getEta(j);
           Float_t phi1 = data.getPhi(j);
           Float_t pt1 = data.getPt(j);
	   if (!isGen&&!data.isChargedHadron(j)) continue;
	   if (pt1<=ptMin||pt1>=ptMax) continue;
           
	   // Signal loop, calculate S correlation function
           for ( Int_t k=j+1;k<data.nParticle;k++ ) {
               Float_t angle2;
	       if (isTheta) angle2 = data.getTheta(k); else angle2 = data.getEta(k);
               Float_t phi2 = data.getPhi(k);
               Float_t pt2 = data.getPt(k);
               if (!isGen&&!data.isChargedHadron(k)) continue;
	       if (pt2<=ptMin||pt2>=ptMax) continue;
               if (N!=0) {
        	   h_2D->Fill(angle1-angle2,dphi(phi1,phi2),1./N_TP);
        	   h_2D->Fill(angle1-angle2,dphi(phi2,phi1),1./N_TP);
        	   h_2D->Fill(angle2-angle1,dphi(phi1,phi2),1./N_TP);
        	   h_2D->Fill(angle2-angle1,dphi(phi2,phi1),1./N_TP);
                   fill_count+=4;
               }
           }
       }
           
       /****************************************/
       // B calculation using multiplicity cut //
       /****************************************/
       for (Int_t nMix = 0; nMix<num_runs; nMix++) {
          t1_mix->GetEntry ( selected );
          t2_mix->GetEntry ( selected );

          // Select a matched event
	  // Currently we are matching the total multiplicity in th event
	  
	  flag=0;
          while ((fabs(mix.nParticle-data.nParticle)>4&&data.nParticle<1000&&fabs(mix.jteta[0]-data.jteta[0])>0.2)||i==selected){
              selected++;
              if (selected>nevent_process&&flag==2) break;
              if (selected>nevent_process) flag++;
              selected = selected % nevent_process;
              t1_mix->GetEntry ( selected );
              t2_mix->GetEntry ( selected );
          }

          // Use the Thrust axis from the signal event instead of mixed event
          mix.TTheta=data.TTheta;
	  mix.TPhi=data.TPhi;
	  mix.update();
	  
          Int_t N2_TP=0;
          
          // calculate the number of tracks in the mixed event passing selection
          for ( Int_t j=0;j<data.nParticle;j++ ) {
             Float_t pt1 = mix.getPt(j);
             if (!isGen&&!mix.isChargedHadron(j)) continue;
	     if (pt1>ptMin&&pt1<ptMax) N2_TP++;
          }
	  if (N2_TP==0) continue;
       
          if (i==selected) {
             cout <<"Error in Mixing!"<<endl;
             continue;
          }
	  
          for ( Int_t j=0;j<data.nParticle;j++ ) {
              Float_t angle1;
	      if (isTheta) angle1 = data.getTheta(j); else angle1 = data.getEta(j);
              Float_t phi1 = data.getPhi(j);
              Float_t pt1 = data.getPt(j);
              if (!isGen&&!data.isChargedHadron(j)) continue;
	      if(pt1<=ptMin||pt1>=ptMax) continue;
           
              // Background loop, calculate B correlation function from mixed event
              for ( Int_t k=0;k<mix.nParticle;k++ ) {
        	 Float_t anglemix;
		 if (isTheta) anglemix = mix.getTheta(k); else anglemix = mix.getEta(k);
        	 Float_t phimix = mix.getPhi(k);
        	 Float_t ptmix = mix.getPt(k);
        	 if (!!isGen&&!mix.isChargedHadron(k)) continue;
		 if(ptmix<=ptMin||ptmix>=ptMax) continue;

		 if (N!=0) {
        	   h_2Dmix->Fill(angle1-anglemix,dphi(phi1,phimix),1./N_TP);
        	   h_2Dmix->Fill(angle1-anglemix,dphi(phimix,phi1),1./N_TP);
        	   h_2Dmix->Fill(anglemix-angle1,dphi(phi1,phimix),1./N_TP);
        	   h_2Dmix->Fill(anglemix-angle1,dphi(phimix,phi1),1./N_TP);
             fill_count_mix+=4;
        	 }
               } //end of mixed event loop
           } // end of working event loop
       } // end of nMix loop      
    } // end of loop over events
    
    cout<<"nEventInMultBin "<<nEventInMultBin<<endl;
    
    h_2Dmix->Scale(1./nEventInMultBin);
    h_2D->Scale(1./nEventInMultBin);
    
    averageN=averageN/nEventProcessed;
    cout <<"Average N = "<<averageN<<endl;
    
    // calculate the  correlation function
    
    calculateRatio(h_2D,h_2Dmix,h_ratio);
        
    // Perform 1D projection
    double etaranges[8]={2,10,2.2,10,2.4,10,2.6,10};
    Int_t minbin,maxbin;
    Int_t thetaranges[4] = {-3,-2,-1,0};
    
    TH1D*h_deltaphi[7];
    
    for (Int_t i=0;i<7;i =i+2){
        // if (i==0) h_deltaphi[i]->SetFillColor(kRed);
        minbin =  h_ratio->GetXaxis()->FindBin(etaranges[i]);
        maxbin =  h_ratio->GetXaxis()->FindBin(etaranges[i+1]);
        
        h_deltaphi[i]  = (TH1D*) h_ratio->ProjectionY(Form("h_deltaphi%d",i),minbin,maxbin);
        //h_deltaphi[i]->Sumw2();
        h_deltaphi[i]->SetName(Form("h_deltaphi_%d",i));
        h_deltaphi[i]->GetXaxis()->SetTitle("#Delta#phi");
        if (isTheta)  h_deltaphi[i]->SetTitle(Form("#Delta#phi, #Delta#theta (%f, %f), Multipliplicity (%d, %d)",etaranges[i],etaranges[i+1], mult_low, mult_high));
        else          h_deltaphi[i]->SetTitle(Form("#Delta#phi, #Delta#eta (%f, %f), Multipliplicity (%d, %d)",etaranges[i],etaranges[i+1], mult_low, mult_high));
        h_deltaphi[i]->GetYaxis()->SetTitle("Y(#Delta#phi)");
        
	h_deltaphi[i]->Scale(1./(maxbin-minbin+1));
    }

    /**************************************************************************************/
    // Save and plot the results
    /**************************************************************************************/
    
    
    TCanvas *c0 = new TCanvas("c0","Multiplicity",600,600);
    h_mult_dist->Draw("e");     

    TCanvas *c = CFViewer("c","S",600,600);
    h_2D->Draw("surf1 fb");    
    
    TCanvas *c1 = CFViewer("c1","B",600,600);
    h_2Dmix->Draw("surf1 fb");    
    
    gStyle->SetOptStat(0);

    TCanvas *c2 = CFViewer("c2","Ratio",600,600);
    h_ratio->Draw("surf1 fb");  
    //c2->SaveAs(Form("ratio_%d_%.1f_%.1f_%d_%d.pdf", isThrust,ptMin,ptMax,mult_low, mult_high));
    
    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextAlign(13);  //align at top

    if (filename.Index("charm") >= 0) {
      latex.DrawLatexNDC(.07, .9, "e^{+}e^{-}#rightarrow c#bar{c} MC");
    } else if (filename.Index("uds") >= 0) {
      latex.DrawLatexNDC(.07, .9, "e^{+}e^{-}#rightarrow q#bar{q} (q=u,d,s) MC");
    } else if (filename.Index("bb") >= 0) {
      latex.DrawLatexNDC(.07, .9, "#Upsilon {4S}#rightarrow B#bar{B} MC");
    }

    char nevt[50];
    sprintf(nevt, "%.2g M events", nevent_process/1000000.);
    latex.DrawLatexNDC(.73, .9, nevt);

    TCanvas *c3 = new TCanvas("c3","dphi",600,600);
    c3->Divide(2,2);
    
    // Fourier decomposition
    TF1 *f = new TF1("f","[0]*(1+2*([1]*cos(1*x)+[2]*cos(2*x)+[3]*cos(3*x)+[4]*cos(4*x)+[5]*cos(5*x)+[6]*cos(6*x)))");
    for (Int_t i=0;i<4;i++) {
      c3->cd(i+1);
      h_deltaphi[i*2]->Draw();
      if (i==0){
         h_deltaphi[0]->Fit("f");
         h_deltaphi[0]->Fit("f");
	 h_deltaphi[0]->SetStats(0);
      }
    }
    c3->SaveAs(Form("proj_%s_%d_%d_%d.pdf",savename.Data(),mult_low, mult_high,isThrust));
             
    // Save the results
    TFile *background = new TFile(Form("ridge_%s_%d_%d_%d.root",savename.Data(),mult_low, mult_high,isThrust), "recreate");
    h_deltaphi[0]->Write();
    h_deltaphi[2]->Write();
    h_deltaphi[4]->Write();
    h_deltaphi[6]->Write();
    h_ratio->Write();
    h_2D->Write();
    h_2Dmix->Write();
    h_mult_dist->Write();
    h_phi->Write();
    //nt->Write();
    f->Write();
    
    TNtuple *ntResult = new TNtuple("ntResult","result","mL:mH:ZYAM:ZYAMError");
    cout <<"Minimum = "<<f->GetMinimum()<<" at dphi = "<<f->GetMinimumX()<<endl;
    cout <<"ZYAM =  "<< f->Integral(0,f->GetMinimumX())-f->GetMinimumX()*f->GetMinimum() <<" +- "<<f->IntegralError(0,f->GetMinimumX())<<endl;
    ntResult->Fill(mult_low,mult_high,f->Integral(0,f->GetMinimumX())-f->GetMinimumX()*f->GetMinimum(),f->IntegralError(0,f->GetMinimumX()));
    ntResult->Write();
    background->Close();
             
             cout<<"FILL COUNT = "<<fill_count<<endl;
             cout<<"FILL COUNT = "<<fill_count_mix<<endl;
    
}

/**************************************************************************************/
// For interactive session in ROOT
/**************************************************************************************/
void ridge_check()
{
   analysis();
}
