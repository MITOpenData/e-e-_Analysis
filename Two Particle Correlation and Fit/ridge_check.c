#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>
#include <stdlib.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TFormula.h>
#include "fourier.h"
using namespace std;

#define PI 3.14159265358979
//enum SIMPLEPID {PHOTON, ELECTRON, PION, MUON, KAON, PROTON};
enum SIMPLEPWLAG {CHARGED_TRACK, CHARGED_LEPTONS1, CHARGED_LEPTONS2, V0, PHOTON, NEUTRAL_HADRON};

// DataFormat
class TPCNtupleData{
    public:
    Int_t nParticle;
    Float_t pt[50000];
    Float_t eta[50000];
    Float_t theta[50000];
    Float_t pid[50000];
    Float_t phi[50000];
    Float_t mass[50000];
    Float_t pwflag[50000];
    
    Float_t px[100000];
    Float_t py[100000];
    Float_t pz[100000];
    Float_t TTheta;
    Float_t TPhi;
};

// Calculate Delta Phi
double dphi(double phi1,double phi2)
{
    double a=phi1-phi2;
    while (a<-PI) a+=2*PI;
    while (a>PI) a-=2*PI;
    
    if (a<-PI/2) a=2*PI+a;
    return a;
}

// Set the branch addresses
void setupTPCTree(TTree *t1, TPCNtupleData &data)
{
    t1->SetBranchAddress("nParticle",&data.nParticle);
    t1->SetBranchAddress("pt",data.pt);
    t1->SetBranchAddress("eta",data.eta);
    t1->SetBranchAddress("theta",data.theta);
    t1->SetBranchAddress("pid",data.pid);
    t1->SetBranchAddress("phi",data.phi);
    t1->SetBranchAddress("mass",data.mass);
    t1->SetBranchAddress("pwflag",data.pwflag);
    
    t1->SetBranchAddress("px",data.px);
    t1->SetBranchAddress("py",data.py);
    t1->SetBranchAddress("pz",data.pz);
    t1->SetBranchAddress("TTheta", &data.TTheta);
    t1->SetBranchAddress("TPhi", &data.TPhi);
}

// Main Analysis 
void analysis(int isBelle    = 0,		// 
              int isThrust   = 0, 		//
	      int maxevt     = 1000000,		// Max number of events To be Processed
	      int mult_low   = 0,		//
	      int mult_high  = 100,		//
	      int nbin       = 26,		//
	      bool verbose   = 0,		//
	      int num_runs   = 1,		//
              double ptMin   = 0.4,             // min pT of the particles
              double ptMax   = 4                // max pT of the particles
	     )
{
    // Input File
    TString filename;
    //filename = "cleaned_ALEPH_Data-all.aleph.root";
    //filename = "/Users/anthony/desktop/LEP2MCGGCCY1997E183_mctrue_aftercut-001.aleph.root";
    filename = "/data/flowex/CMSsample/cleaned_MinBias_TuneCUETP8M1_5p02TeV-pythia8-HINppWinter16DR-NoPU_75X_mcRun2_asymptotic_ppAt5TeV_forest_v2_track.root";
        
    TFile *f = new TFile(filename.Data());
    TTree *t1 = (TTree*)f->Get("t");

    TPCNtupleData data;
    setupTPCTree(t1,data);
    
    // File for event mixing, use the same file for the moment
    TFile *f_mix = new TFile(filename.Data());
    TTree *t1_mix = (TTree*)f_mix->Get("t");
    
    TPCNtupleData mix;
    setupTPCTree(t1_mix,mix);
    
    // two histograms    
    double detaRange = 4;//1.8
    double dthetaRange = PI;
    double normalization = detaRange*2/nbin*2*3.14159/nbin;
    TH2F *h_2D = new TH2F ( "h_2D", "eta-phi of all particles ",nbin, -detaRange, detaRange,nbin, -3.1416/2., 3.1416*1.5);
    TH2F *h_2Dmix = new TH2F ( "h_2Dmix", "eta-phi of all particles ",nbin, -detaRange, detaRange,nbin, -3.1416/2., 3.1416*1.5);
    TH2F *h_ratio = new TH2F ( "h_ratio", "eta-phi of all particles ", nbin, -detaRange, detaRange,nbin, -3.1416/2.,3.1416*1.5);
    TH1D *mult_dist = new TH1D("multi","multiplicity distribution of charged particles",20,0,20);
    h_2D->Sumw2();
    h_2Dmix->Sumw2();
    h_ratio->Sumw2();
    
    // all entries and fill the histograms
    Int_t nevent = (Int_t)t1->GetEntries();
    
    int nevent_process = nevent;
    if( maxevt>0 && maxevt<nevent ) nevent_process = maxevt;
    
    double averageN=0;
    double nEventProcessed=0;
    double nEventInMultBin=0;
    
    
    // S calculation using multiplicity cut
    for (Int_t i=0;i<nevent_process;i++) {
        
        t1->GetEntry(i);
        if (i%10000==0) cout <<i<<"/"<<nevent_process<<endl;
        
        if (verbose) cout<<"nparticles="<<data.nParticle<<endl;
        
        int selected=i+1;  //questo diventa il numero di evento estratto nel file +1
        int flag=0;        //definisco una flag a zero
        
        // Yen-Jie: cut on maximum number of particles to avoid infinite loop, for the moment it is 100
        // if (nparticles>mult_high) continue;
        
        double N=0;
        for (int p = 0;p<num_runs;p++)
        {
            // calculate the number of tracks in the passing selection
            for ( int j=0;j<data.nParticle;j++ ) {
                float pt1 = data.pt[j];
                int pid1 = data.pid[j];
                int pwflag1 = data.pwflag[j];
                //if (pid1!=PION&&pid1!=PROTON&&pid1!=KAON&&!(!isBelle&&pwflag1==0)) continue;
                if (pwflag1!=CHARGED_TRACK){continue;}
                if(pt1<ptMin||pt1>ptMax) {continue;}
                N++;
            }
//            cout<<"N"<<N<<endl;
            mult_dist->Fill(N);
            averageN+=N;
            nEventProcessed++;
            
            //if (nparticles<mult_low) continue;
            if (N<mult_low) continue;
            if (N>mult_high) continue;
            nEventInMultBin++;
            
            // find a mixed event
            while ((fabs(mix.nParticle-data.nParticle)>5&&data.nParticle<mult_high)||i==selected){
                selected++;
                if (selected>nevent_process&&flag==1) break;
                if (selected>nevent_process) flag=1;
                selected = selected % nevent_process;
                t1_mix->GetEntry ( selected );
            }
            
            double N2=0;
            // calculate the number of tracks in the mixed event passing selection
            for ( int j=0;j<data.nParticle;j++ ) {
                float pt1 = mix.pt[j];
                int pid1 = mix.pid[j];
                int pwflag1 = mix.pwflag[j];
                if (pwflag1!=CHARGED_TRACK) continue;
                //if (pid1!=PION&&pid1!=PROTON&&pid1!=KAON&&!(!isBelle&&pwflag1==0)) continue;
                if(pt1<ptMin||pt1>ptMax) continue;
                N2++;
            }
            
            TVector3 thrust(1, 1, 1);
            thrust.SetTheta(data.TTheta);
            thrust.SetPhi(data.TPhi);
            
            
            TVector3 rotVec;
            rotVec.SetX(-thrust.Y());
            rotVec.SetY(thrust.X());
            
            TVector3 thrust_mix(1, 1, 1);
            thrust_mix.SetTheta(mix.TTheta);
            thrust_mix.SetPhi(mix.TPhi);
            
            TVector3 rotVec_mix;
            rotVec_mix.SetX(-thrust_mix.Y());
            rotVec_mix.SetY(thrust_mix.X());
            
            for ( int j=0;j<data.nParticle;j++ ) {
                int pid1 = data.pid[j];
                float eta1 = data.eta[j];
                float phi1 = data.phi[j];
                float pt1 = data.pt[j];
                float mass1 = data.mass[j];
                float pwflag1 = data.pwflag[j];
                float theta1 = data.theta[j];
                if (isThrust == 1)
                {
                    TVector3 v1;
                    float px1 = data.px[j];
                    float py1 = data.py[j];
                    float pz1 = data.pz[j];
                    v1.SetXYZ(px1,py1,pz1);
                    v1.Rotate(-data.TTheta, rotVec);
                    phi1 = v1.Phi();
                    eta1 = v1.PseudoRapidity();
                    pt1 = v1.Pt();
                    theta1 = v1.Theta();
                    //
                }
                //cout<<phi1<<endl;
                if(pwflag1!=CHARGED_TRACK)continue;
                //if (pid1!=PION&&pid1!=PROTON&&pid1!=KAON&&!(!isBelle&&pwflag1==0)) continue;
                if(pt1<ptMin||pt1>ptMax) continue;
                if(p==0)
                {
                    // Signal loop, calculate S correlation function
                    for ( int k=j+1;k<data.nParticle;k++ ) {
                        int pid2 = data.pid[k];
                        float eta2 = data.eta[k];
                        float phi2 = data.phi[k];
                        float pt2 = data.pt[k];
                        float mass2 = data.mass[k];
                        float pwflag2 = data.pwflag[k];
                        float theta2 = data.theta[k];
                        if (isThrust == 1)
                        {
                            TVector3 v2;
                            float px2 = data.px[k];
                            float py2 = data.py[k];
                            float pz2 = data.pz[k];
                            v2.SetXYZ(px2,py2,pz2);
                            v2.Rotate(-data.TTheta, rotVec);
                            phi2 = v2.Phi();
                            eta2 = v2.PseudoRapidity();
                            pt2 = v2.Pt();
                            theta2 = v2.Theta();
                            //
                        }
                        if(pwflag2!=CHARGED_TRACK)continue;
                        if(pt2<ptMin||pt2>ptMax) continue;
                        if (N!=0)
                        {
                            //cout<<"HELLO"<<endl;
                            h_2D->Fill(eta1-eta2,dphi(phi1,phi2),1./N);
                            h_2D->Fill(eta1-eta2,dphi(phi2,phi1),1./N);
                            h_2D->Fill(eta2-eta1,dphi(phi1,phi2),1./N);
                            h_2D->Fill(eta2-eta1,dphi(phi2,phi1),1./N);
                        }
                    }
                }//end of second loop
            }
        }// end of first loop
    }// end of loop over events
    
    // BACKGROUND LOOP using no multiplicity cut
    for (Int_t i=0;i<nevent_process;i++) {
        
        if (i%10000==0) cout <<i<<"/"<<nevent_process<<endl;
        t1->GetEntry(i);
        
        if (verbose) cout<<"nparticles="<<data.nParticle<<endl;
        
        int selected=i+1;  //questo diventa il numero di evento estratto nel file +1
        int flag=0;    //definisco una flag a zero
        
        // Yen-Jie: cut on maximum number of particles to avoid infinite loop, for the moment it is 100
        if (data.nParticle>1000) continue;
        
        
        double N=0;
        for (int p = 0;p<num_runs;p++)
        {
            // calculate the number of tracks in the passing selection
            for ( int j=0;j<data.nParticle;j++ ) {
                float pt1 = data.pt[j];
                int pid1 = data.pid[j];
                int pwflag1 = data.pwflag[j];
                //if (pid1!=PION&&pid1!=PROTON&&pid1!=KAON&&!(!isBelle&&pwflag1==0)) continue;
                if (pwflag1!=CHARGED_TRACK)continue;
                if(pt1<ptMin||pt1>ptMax) continue;
                N++;
            }
            
            averageN+=N;
            nEventProcessed++;
            
            //if (nparticles<mult_low) continue;
            if (N<0) continue;
            
            nEventInMultBin++;
            
            // find a mixed event
            //	cout <<N<<endl;
            while ((fabs(mix.nParticle-data.nParticle)>40&&data.nParticle<1000)||i==selected){
                //cout <<nparticles<<" "<<selected<<endl;
                selected++;
                if (selected>nevent_process&&flag==1) break;
                if (selected>nevent_process) flag=1;
                selected = selected % nevent_process;
                t1_mix->GetEntry ( selected );
            }
            
            double N2=0;
            // calculate the number of tracks in the mixed event passing selection
            for ( int j=0;j<data.nParticle;j++ ) {
                float pt1 = mix.pt[j];
                int pid1 = mix.pid[j];
                int pwflag1 = mix.pwflag[j];
                if (pwflag1!=CHARGED_TRACK) continue;
                //if (pid1!=PION&&pid1!=PROTON&&pid1!=KAON&&!(!isBelle&&pwflag1==0)) continue;
                if(pt1<ptMin||pt1>ptMax) continue;
                N2++;
            }
            
            TVector3 thrust(1, 1, 1);
            thrust.SetTheta(data.TTheta);
            thrust.SetPhi(data.TPhi);
            
            
            TVector3 rotVec;
            rotVec.SetX(-thrust.Y());
            rotVec.SetY(thrust.X());
            
            TVector3 thrust_mix(1, 1, 1);
            thrust_mix.SetTheta(mix.TTheta);
            thrust_mix.SetPhi(mix.TPhi);
            
            TVector3 rotVec_mix;
            rotVec_mix.SetX(-thrust_mix.Y());
            rotVec_mix.SetY(thrust_mix.X());
            
            for ( int j=0;j<data.nParticle;j++ ) {
                int pid1 = data.pid[j];
                float eta1 = data.eta[j];
                float phi1 = data.phi[j];
                float pt1 = data.pt[j];
                float mass1 = data.mass[j];
                float pwflag1 = data.pwflag[j];
                float theta1 = data.theta[j];
                if (isThrust == 1)
                {
                    TVector3 v1;
                    float px1 = data.px[j];
                    float py1 = data.py[j];
                    float pz1 = data.pz[j];
                    v1.SetXYZ(px1,py1,pz1);
                    v1.Rotate(-data.TTheta, rotVec);
                    phi1 = v1.Phi();
                    eta1 = v1.PseudoRapidity();
                    pt1 = v1.Pt();
                    theta1 = v1.Theta();
                }
                //cout<<phi1<<endl;
                if(pwflag1!=CHARGED_TRACK)continue;
                //if (pid1!=PION&&pid1!=PROTON&&pid1!=KAON&&!(!isBelle&&pwflag1==0)) continue;
                if(pt1<ptMin||pt1>ptMax) continue;
                
                // Background loop, calculate B correlation function from mixed event
                for ( int k=0;k<mix.nParticle;k++ ) {
                    int pidmix = mix.pid[k];
                    float etamix = mix.eta[k];
                    float phimix = mix.phi[k];
                    float ptmix = mix.pt[k];
                    float massmix = mix.mass[k];
                    float pwflagmix = mix.pwflag[k];
                    float thetamix = mix.theta[k];
                    if (isThrust == 1)
                    {
                        TVector3 vmix;
                        float pxmix = mix.px[j];
                        float pymix = mix.py[j];
                        float pzmix = mix.pz[j];
                        vmix.SetXYZ(pxmix,pymix,pzmix);
                        vmix.Rotate(-mix.TTheta, rotVec_mix);
                        phimix = vmix.Phi();
                        etamix = vmix.PseudoRapidity();
                        ptmix = vmix.Pt();
                        thetamix = vmix.Theta();
                        
                    }
                    if(pwflagmix!= CHARGED_TRACK)continue;
                    //if (pidmix!=PION&&pidmix!=PROTON&&pidmix!=KAON&&!(!isBelle&&pwflagmix==0)) continue;
                    if(ptmix<ptMin||ptmix>ptMax) continue;
                    if(N!=0)
                    {
                        h_2Dmix->Fill(eta1-etamix,dphi(phi1,phimix),1./N);
                        h_2Dmix->Fill(eta1-etamix,dphi(phimix,phi1),1./N);
                        h_2Dmix->Fill(etamix-eta1,dphi(phi1,phimix),1./N);
                        h_2Dmix->Fill(etamix-eta1,dphi(phimix,phi1),1./N);
                        
                    }
                }//end of second loop
                
            }
        }// end of first loop
    }
    
    
    cout<<"nEventInMultBin"<<nEventInMultBin<<endl;
    
    h_2Dmix->Scale(1./nEventInMultBin);
    h_2D->Scale(1./nEventInMultBin);
    
    averageN=averageN/nEventProcessed;
    cout <<"Average N = "<<averageN<<endl;
    double ratio;
    double errrel_ratio;
    double errrel_num;
    double errrel_den;
    
    // calculate the  correlation function
    
    double b00_x=h_2Dmix->GetXaxis()->FindBin(0.);
    double b00_y=h_2Dmix->GetYaxis()->FindBin(0.);
    double B00=h_2Dmix->GetBinContent(b00_x,b00_y);
    double errrel_B00=h_2Dmix->GetBinError(b00_x,b00_y)/B00;
    
    cout<<"value of B(0,0)="<<B00<<endl;
    
    
    cout<<"x axis "<<h_2Dmix->GetXaxis()->GetBinCenter(b00_x)<<endl;
    cout<<"y axis "<<h_2Dmix->GetYaxis()->GetBinCenter(b00_y)<<endl;
    
    TH2F *h_ratio_theta = new TH2F ( "h_ratio_theta", "theta-phi of all particles ", nbin, -PI*1.1, PI/6.,nbin, -3.1416/2.,3.1416*1.5);
    
    h_ratio_theta->Sumw2();

    double ratio_phi;
    double ratio_eta;
    
    //float deta_count = 2*detaRange/(float)nbin;
    //float dphi_count = (3.1416*1.5 - -3.1416/2.)/(float)nbin;
    
    for (int x=0;x<=h_2D->GetNbinsX();x++){
        
        for (int y=0;y<=h_2D->GetNbinsY();y++){
            
            
            if(h_2Dmix->GetBinContent(x,y)>0){
                
                ratio=B00*(h_2D->GetBinContent(x,y)/h_2Dmix->GetBinContent(x,y));
                errrel_num=h_2D->GetBinError(x,y)/h_2D->GetBinContent(x,y);
                errrel_den=h_2Dmix->GetBinError(x,y)/h_2Dmix->GetBinContent(x,y);
                errrel_ratio=TMath::Sqrt(errrel_num*errrel_num+errrel_den*errrel_den+errrel_B00*errrel_B00);
                
                ratio_eta = h_2D->GetXaxis()->GetBinCenter(x);
                ratio_phi = h_2D->GetYaxis()->GetBinCenter(y);
                h_ratio->Fill(ratio_eta,ratio_phi,ratio/normalization);
            }
        }
        
    }
    
    
    
    for(int binIterX = 0; binIterX < h_ratio->GetNbinsX(); ++binIterX){
        for(int binIterY = 0; binIterY < h_ratio->GetNbinsX(); ++binIterY){
            double eta = h_ratio->GetXaxis()->GetBinCenter(binIterX+1);
            double phi = h_ratio->GetYaxis()->GetBinCenter(binIterY+1);
            double val = h_ratio->GetBinContent(binIterX,binIterY);
            double err = h_ratio->GetBinError(binIterX,binIterY);
            double theta = -2*TMath::ATan(exp(-eta));
            h_ratio_theta->Fill(theta,phi,val);
        }
        
    }
    cout<<"bin"<<h_ratio_theta->GetBinContent(0,0)<<endl;
    double etaranges[8]={1.5,4,1.5,3, 2,3, 2,4};
    int minbin,maxbin;
    int thetaranges[4] = {-3,-2,-1,0};
    
    TH1D*h_deltaphi[7];
    TH1D*h_deltaphi_theta[3];
    
    //h_deltaphi  = (TH1D*) h_ratio->ProjectionY("h_deltaphi",0,-1);
    //h_deltaphi->Sumw2();
    
    for (int i=0;i<7;i =i+2){
        // if (i==0) h_deltaphi[i]->SetFillColor(kRed);
        minbin =  h_ratio->GetXaxis()->FindBin(etaranges[i]);
        maxbin =  h_ratio->GetXaxis()->FindBin(etaranges[i+1]);
        cout<<"1"<<" "<<minbin<<" "<<maxbin<<endl;
        
        h_deltaphi[i]  = (TH1D*) h_ratio->ProjectionY(Form("h_deltaphi%d",i),minbin,maxbin);
        h_deltaphi[i]->Sumw2();
        
        h_deltaphi[i]->SetName(Form("h_deltaphi_%d",i));
        
        h_deltaphi[i]->GetXaxis()->SetTitle("#Delta#phi");
        h_deltaphi[i]->SetTitle(Form("#Delta#phi, #Delta#eta (%f, %f), Multipliplicity (%d, %d)",etaranges[i],etaranges[i+1], mult_low, mult_high));
        h_deltaphi[i]->GetYaxis()->SetTitle("Y(#Delta#phi)");
        cout<<1./(maxbin-minbin+1)<<endl;
        h_deltaphi[i]->Scale(1./(maxbin-minbin+1));
    }
    
    
    for (int i=0;i<3;i++){
        minbin =  h_ratio_theta->GetXaxis()->FindBin(thetaranges[i]);
        maxbin =  h_ratio_theta->GetXaxis()->FindBin(thetaranges[i+1]);
        cout<<minbin<<" "<<maxbin<<endl;
        h_deltaphi_theta[i]  = (TH1D*) h_ratio_theta->ProjectionY(Form("h_deltaphi_thetamin%d_max%d",thetaranges[i],thetaranges[i+1]),minbin,maxbin);
        
        h_deltaphi_theta[i]->Sumw2();
        
        h_deltaphi_theta[i]->GetXaxis()->SetTitle("#Delta#phi For #theta");
        h_deltaphi_theta[i]->SetTitle(Form("deltaphi theta (%d, %d)",thetaranges[i],thetaranges[i+1]));
        
        h_deltaphi_theta[i]->Scale(1./(maxbin-minbin+1));
    }
    
    
    h_ratio->GetXaxis()->SetTitle("#Delta#eta");
    h_ratio->GetYaxis()->SetTitle("#Delta#phi");
    
    TFile *background = new TFile(Form("correlation_%d_%d_%d.root", isThrust,mult_low, mult_high), "recreate");
    h_deltaphi[0]->Write();
    h_ratio->Write();
    h_2D->Write();
    h_2Dmix->Write();
    mult_dist->Write();
    background->Close();
    
    TCanvas *c = new TCanvas("c","S",600,600);
    c->SetTheta(60.839);
    c->SetPhi(38.0172);
    h_2D->Draw("surf1");    
    
    TCanvas *c1 = new TCanvas("c1","B",600,600);
    c1->SetTheta(60.839);
    c1->SetPhi(38.0172);
    h_2Dmix->Draw("surf1");    
    
    TCanvas *c2 = new TCanvas("c2","Ratio",600,600);
    c2->SetTheta(60.839);
    c2->SetPhi(38.0172);
    h_ratio->Draw("surf1");    
}
