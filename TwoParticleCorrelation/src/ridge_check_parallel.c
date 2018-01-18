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

/********************************************************************************************************************/
// Two particle correlation analysis
//
// ridge_check_parallel.c
//
// Yen-Jie Lee, Gian Michele Innocenti, Anthony Badea, Austin Baty, Chris McGinn, Michael Peters, Tzu-An Sheng
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
     const std::string outFileName,    // output file
     Int_t isThrust     = 0,         // Leading Jet Axis ridge_check = 2, Thurst Axis ridge_check = 1, Beam Axis ridge_check = 0
     Int_t isBelle      = 0,        // BELLE ridge_check = 1, CMS/ALEPH ridge_check = 0
     Int_t isTheta      = 0,         // Use Theta angle = 1, Use Eta = 0
     Int_t isGen       = 0,         // Use isGen = 1 if gen level so no cut on charged particles, 0 if not gen level
     Int_t maxevt       = 1000000,    // Max number of events to be processed, 0 = all events
     Int_t mult_low     = 0,        // Lower cut on the event multiplicity
     Int_t mult_high    = 100,        // Upper cut on the event multiplicity
     Int_t nbin         = 20,        // Number of bins in the correlation function
     bool verbose     = 0,        // Verbose mode
     Int_t num_runs     = 5,        //
     double ptMin     = 0.4,           // min pT of the particles used for correlation function
     double ptMax     = 100.0,             // max pT of the particles used for correlation function
     double detaRange = 3.2,             // deta window of the correlation function
     Float_t ptMinForN = 0.4,          // pT min for the N_{trk}^{offline} calculation (used for event classification)
     Float_t ptMaxForN = 100,          // pT max for the N_{trk}^{offline} calculation (used for event classification)
     Float_t etaCutForN = 1.8          // eta window for the N_{trk}^{offline} calculation (used for event classification)
    )
{
    // ROOT Global setting
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    TChain *t1 = new TChain("t");
    t1->Add(inFileName.c_str());
    
    TChain *t2 = new TChain("ak4ESchemeJetTree");
    t2->Add(inFileName.c_str());
    
    /* TString friendname =  "qqmc-e07-00_flavor.root"; */
    /* TFile *ff = new TFile(friendname); */
    /* TTree *tf = (TTree*)ff->Get("t"); */
    /* t1->AddFriend(tf); */
    
    TPCNtupleData data(isBelle, isThrust);
    setupTPCTree(t1,t2,data);
    data.setTPCTreeStatus(t1);
    
    // File for event mixing, use the same file for the moment
    /* TFile *f_mix = new TFile(inFileName.Data()); */
    /* TTree *t1_mix = (TTree*)f_mix->Get("t"); */
    TChain *t1_mix = new TChain("t");
    t1_mix->Add(inFileName.c_str());
    
    // Not necessary
    TChain *t2_mix = new TChain("ak4ESchemeJetTree");
    t2_mix->Add(inFileName.c_str());
    
    TPCNtupleData mix(isBelle, isThrust);
    setupTPCTree(t1_mix,t2_mix,mix);
    mix.setTPCTreeStatus(t1_mix);
    
    //TNtuple *nt = new TNtuple("nt","","pEta:pTheta:pPhi:theta:phi:TTheta:TPhi");
    
    // Define 2D histograms
    //Float_t dthetaRange = PI;
    //Float_t normalization = detaRange*2/nbin*2*3.14159/nbin;
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
    //cout<<"maxevt "<<maxevt<<endl;
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
        if (i%10000==0) std::cout <<i<<"/"<<nevent_process<<std::endl;
        if (verbose) std::cout<<"nparticles="<<data.nParticle<<std::endl;
        
        //if (isThrust==1&&fabs(data.TTheta-3.14159/2.)>0.8) continue;
        //if (isThrust==2&&acos(cos(data.jtphi[0]-data.jtphi[1]))<3.14159*2/3.) continue;
        
        /* if (flavor == 0) {break;} */
        /* if (flavor == 1) {continue;} */
        
        Int_t selected=i+1;  //questo diventa il numero di evento estratto nel file +1
        Int_t flag=0;      //definisco una flag a zero
        
        // Yen-Jie: we would like to use a definition similar to CMS publication
        // all particles with pT > 0.4 GeV/c for the calculation of event multiplicity N
        
        Int_t N=0;    // N_{trk}^{offline}
        Int_t N_TP=0;    // Number of particles used in the two particle correlation function calculation, can be different from N for event classification
        
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
                std::cout <<"Error in Mixing!"<<std::endl;
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
    
    std::cout<<"nEventInMultBin "<<nEventInMultBin<<std::endl;
    
    h_2Dmix->Scale(1./nEventInMultBin);
    h_2D->Scale(1./nEventInMultBin);
    
    averageN=averageN/nEventProcessed;
    std::cout <<"Average N = "<<averageN<<std::endl;
    
    // calculate the  correlation function
    
    calculateRatio(h_2D,h_2Dmix,h_ratio);
    
    // Save the results
    TFile *background = new TFile(Form("../pdfDir/%s.root",outFileName.c_str()), "recreate");
    h_ratio->Write();
    h_2D->Write();
    h_2Dmix->Write();
    h_phi->Write();
    h_mult_dist->Write();
    
    std::cout<<"FILL COUNT = "<<fill_count<<std::endl;
    std::cout<<"FILL COUNT = "<<fill_count_mix<<std::endl;
    
    background->Close();
    delete background;
    delete h_mult_dist;
    delete h_ratio;
    delete h_2Dmix;
    delete h_2D;
    delete h_phi;
    delete t2_mix;
    delete t1_mix;
    delete t2;
    delete t1;
    return 0;
}

int main(int argc, char* argv[])
{
    if(argc < 3 || argc > 19)
    {
        std::cout << "Usage: ./ridge_check_parallel.exe <inFileName> <outFileName> <isThrust-optional> <isBelle-optional> <isTheta-optional> <isGen-optional> <maxevt-optional> <mult_low-optional> <mult_high-optional> <nbin-optional> <verbose-optional> <num_runs-optional> <ptMin-optional> <ptMax-optional> <detaRange-optional> <ptMinForN-optional> <ptMaxForN-optional> <etaCutForN-optional>" <<std::endl;
        return 1;
    }
    
    int retVal = 0;
    if(argc == 3) retVal += ridge_check_parallel(argv[1], argv[2]);
    else if(argc == 4) retVal += ridge_check_parallel(argv[1], argv[2], std::stof(argv[3]));
    else if(argc == 5) retVal += ridge_check_parallel(argv[1], argv[2], std::stof(argv[3]), std::stof(argv[4]));
    else if(argc == 6) retVal += ridge_check_parallel(argv[1], argv[2], std::stof(argv[3]), std::stof(argv[4]), std::stof(argv[5]));
    else if(argc == 7) retVal += ridge_check_parallel(argv[1], argv[2], std::stof(argv[3]), std::stof(argv[4]), std::stof(argv[5]), std::stof(argv[6]));
    else if(argc == 8) retVal += ridge_check_parallel(argv[1], argv[2], std::stof(argv[3]), std::stof(argv[4]), std::stof(argv[5]), std::stof(argv[6]), std::stof(argv[7]));
    else if(argc == 9) retVal += ridge_check_parallel(argv[1], argv[2], std::stof(argv[3]), std::stof(argv[4]), std::stof(argv[5]), std::stof(argv[6]), std::stof(argv[7]), std::stoi(argv[8]));
    else if(argc == 10) retVal += ridge_check_parallel(argv[1], argv[2], std::stof(argv[3]), std::stof(argv[4]), std::stof(argv[5]), std::stof(argv[6]), std::stof(argv[7]), std::stoi(argv[8]), std::stoi(argv[9]));
    else if(argc == 11) retVal += ridge_check_parallel(argv[1], argv[2], std::stof(argv[3]), std::stof(argv[4]), std::stof(argv[5]), std::stof(argv[6]), std::stof(argv[7]), std::stoi(argv[8]), std::stoi(argv[9]), std::stoi(argv[10]));
    else if(argc == 12) retVal += ridge_check_parallel(argv[1], argv[2], std::stof(argv[3]), std::stof(argv[4]), std::stof(argv[5]), std::stof(argv[6]), std::stof(argv[7]), std::stoi(argv[8]), std::stoi(argv[9]), std::stoi(argv[10]),std::stoi(argv[11]));
    else if(argc == 13) retVal += ridge_check_parallel(argv[1], argv[2], std::stof(argv[3]), std::stof(argv[4]), std::stof(argv[5]), std::stof(argv[6]), std::stof(argv[7]), std::stoi(argv[8]), std::stoi(argv[9]), std::stoi(argv[10]),std::stoi(argv[11]),std::stoi(argv[12]));
    else if(argc == 14) retVal += ridge_check_parallel(argv[1], argv[2], std::stof(argv[3]), std::stof(argv[4]), std::stof(argv[5]), std::stof(argv[6]), std::stof(argv[7]), std::stoi(argv[8]), std::stoi(argv[9]), std::stoi(argv[10]),std::stoi(argv[11]),std::stoi(argv[12]),std::stoi(argv[13]));
    else if(argc == 15) retVal += ridge_check_parallel(argv[1], argv[2], std::stof(argv[3]), std::stof(argv[4]), std::stof(argv[5]), std::stof(argv[6]), std::stof(argv[7]), std::stoi(argv[8]), std::stoi(argv[9]), std::stoi(argv[10]),std::stoi(argv[11]),std::stoi(argv[12]),std::stoi(argv[13]),std::stoi(argv[14]));
    else if(argc == 16) retVal += ridge_check_parallel(argv[1], argv[2], std::stof(argv[3]), std::stof(argv[4]), std::stof(argv[5]), std::stof(argv[6]), std::stof(argv[7]), std::stoi(argv[8]), std::stoi(argv[9]), std::stoi(argv[10]),std::stoi(argv[11]),std::stoi(argv[12]),std::stoi(argv[13]),std::stoi(argv[14]),std::stoi(argv[15]));
    else if(argc == 17) retVal += ridge_check_parallel(argv[1], argv[2], std::stof(argv[3]), std::stof(argv[4]), std::stof(argv[5]), std::stof(argv[6]), std::stof(argv[7]), std::stoi(argv[8]), std::stoi(argv[9]), std::stoi(argv[10]),std::stoi(argv[11]),std::stoi(argv[12]),std::stoi(argv[13]),std::stoi(argv[14]),std::stoi(argv[15]),std::stoi(argv[16]));
    else if(argc == 18) retVal += ridge_check_parallel(argv[1], argv[2], std::stof(argv[3]), std::stof(argv[4]), std::stof(argv[5]), std::stof(argv[6]), std::stof(argv[7]), std::stoi(argv[8]), std::stoi(argv[9]), std::stoi(argv[10]),std::stoi(argv[11]),std::stoi(argv[12]),std::stoi(argv[13]),std::stoi(argv[14]),std::stoi(argv[15]),std::stoi(argv[16]),std::stoi(argv[17]));
    else if(argc == 19) retVal += ridge_check_parallel(argv[1], argv[2], std::stof(argv[3]), std::stof(argv[4]), std::stof(argv[5]), std::stof(argv[6]), std::stof(argv[7]), std::stoi(argv[8]), std::stoi(argv[9]), std::stoi(argv[10]),std::stoi(argv[11]),std::stoi(argv[12]),std::stoi(argv[13]),std::stoi(argv[14]),std::stoi(argv[15]),std::stoi(argv[16]),std::stoi(argv[17]),std::stoi(argv[18]));
    
    return retVal;
}
