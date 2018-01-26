#include "TFile.h"
#include "TTree.h"
#include "TAttLine.h"
#include "THelix.h"
#include "TView.h"
#include <string>
#include "TAttPad.h"
#include "TMath.h"
#include "TVector3.h"
#include "TView3D.h"
#include "TCanvas.h"

#include "include/smartJetName.h"

//nasty global variables for nextEvent() function
std::string currentFile;
int currentEvtIndx;
float currentPtCut;

float dR(float eta1, float phi1, float eta2, float phi2){
  return TMath::Power(TMath::Power(TMath::ACos(TMath::Cos(phi1-phi2)),2)+TMath::Power(eta1-eta2,2),0.5);
}

void formatHelix(THelix * h, float pz, float pt, int color = 0, bool doWTA = false){
  if(color==0) h->SetLineColor(kBlue);
  if(color==1) h->SetLineColor(kWhite);
  if(color==2) h->SetLineColor(kRed);
  if(color==3) h->SetLineColor(6);
  if(color==4) h->SetLineColor(7);
  if(color==5) h->SetLineColor(kGreen);
  h->SetLineWidth(1);
  
  float rangeBound = 1;
  if(!doWTA){
    if(pt<2.5 && TMath::Abs(pz)<0.5) rangeBound = TMath::Abs(pz);
  }
  h->SetRange(0,rangeBound);
  if(pz<0) h->SetRange(-rangeBound,0);
}


void EventDisplay(std::string inputFile = "/data/cmcginn/StudyMultSamples/ALEPH/LEP1/20180119/LEP1Data1995_recons_aftercut-MERGED.root" ,int eventIndx = 574, float ptCut = 0, bool doWTA = false){

  currentFile=inputFile;
  currentEvtIndx=eventIndx;
  currentPtCut=ptCut;

  TFile * f = TFile::Open(inputFile.c_str(),"read");
  TTree * t = (TTree*)f->Get("t");
  TTree * WTA_t=0;
  if(doWTA){
    WTA_t = (TTree*)f->Get("BoostedWTAR8Evt");
  }
  TTree * jt = (TTree*)f->Get(smartJetName("ak8ESchemeJetTree", f).c_str()); 
  TTree * wta = (TTree*)f->Get(smartJetName("ak8WTAmodpSchemeJetTree", f).c_str());

  int nParticle;
  float pt[500];
  float eta[500];
  float phi[500];
  float px[500];
  float py[500];
  float pz[500];
  float charge[500];
  int pwflag[500], pwflagMix[500];


  float wtapt[500];
  float wtaeta[500];
  float wtaphi[500];

  int nref;
  float jtpt[50];
  float jteta[50];
  float jtphi[50];
  
  int wtanref;
  float wtajtpt[50];
  float wtajteta[50];
  float wtajtphi[50];

  float TTheta;
  float TPhi;
  t->SetBranchAddress("TTheta_charged",&TTheta);
  t->SetBranchAddress("TPhi_charged",&TPhi);    
  t->SetBranchAddress("nParticle",&nParticle);
  t->SetBranchAddress("pt",&pt);
  t->SetBranchAddress("phi",&phi);
  t->SetBranchAddress("eta",&eta);
  t->SetBranchAddress("px",&px);              
  t->SetBranchAddress("py",&py);              
  t->SetBranchAddress("pz",&pz);     
  if(doWTA){
    WTA_t->SetBranchAddress("pt",&wtapt);
    WTA_t->SetBranchAddress("phi",&wtaphi);
    WTA_t->SetBranchAddress("eta",&wtaeta);
  }        
  t->SetBranchAddress("pwflag",&pwflag);     
  t->SetBranchAddress("charge",&charge); 
  jt->SetBranchAddress("nref",&nref);
  jt->SetBranchAddress("jtpt",&jtpt);
  jt->SetBranchAddress("jteta",&jteta);
  jt->SetBranchAddress("jtphi",&jtphi);
  wta->SetBranchAddress("nref",&wtanref);
  wta->SetBranchAddress("jtpt",&wtajtpt);
  wta->SetBranchAddress("jteta",&wtajteta);
  wta->SetBranchAddress("jtphi",&wtajtphi);
  t->GetEntry(eventIndx);
  if(doWTA) WTA_t->GetEntry(eventIndx);
  jt->GetEntry(eventIndx);
  wta->GetEntry(eventIndx);

  TVector3 thrust = TVector3(0,0,0);
  thrust.SetMagThetaPhi(10,TTheta,TPhi);
  TVector3 wta1 = TVector3(0,0,0);
  if(wtanref>0) wta1.SetMagThetaPhi(10,2*TMath::ATan(TMath::Exp(-wtajteta[0])),wtajtphi[0]);
  TVector3 wta2 = TVector3(0,0,0);
  if(wtanref>1) wta2.SetMagThetaPhi(10,2*TMath::ATan(TMath::Exp(-wtajteta[1])),wtajtphi[1]);


  TCanvas* c1 = new TCanvas("AzimuthalView","AzimuthalView",600,600);
  c1->SetFillColor(kBlack);  
  TView *view = TView::CreateView(1);
  view->SetRange(-1,-1,-1,1,1,1);
  view->TopView();

  TCanvas* c2 = new TCanvas("TopView","TopView",600,600);
  c2->SetFillColor(kBlack);  
  TView *view2 = TView::CreateView(1);
  view2->SetRange(-1,-1,-1,1,1,1);
  view2->SideView();
  
  TCanvas* c3 = new TCanvas("AzimuthalThrustView","AzimuthalThrustView",600,600);
  c3->SetFillColor(kBlack);  
  TView *view3 = TView::CreateView(1);
  view3->SetRange(-1,-1,-1,1,1,1);
  view3->RotateView(TPhi*180/TMath::Pi(),TTheta*180/TMath::Pi());
  view3->ZoomView(c3,10);

  THelix * helix[1000];
  if(!doWTA){
    helix[999] = new THelix(0,0,0,thrust.Px(),thrust.Py(),thrust.Pz(),0.000001);
    helix[999]->SetRange(-1,1);
    helix[999]->SetLineColor(8);
    helix[999]->SetLineWidth(3);
    c1->cd();
    helix[999]->Draw();
    c2->cd();
    helix[999]->Draw();
    c3->cd();
    helix[999]->Draw();
    
    helix[998] = new THelix(0,0,0,wta1.Px(),wta1.Py(),wta1.Pz(),0.000001);
    if(wta1.Pz()<0) helix[998]->SetRange(-1,0);
    if(wta1.Pz()>=0) helix[998]->SetRange(0,1);
    helix[998]->SetLineColor(38);
    helix[998]->SetLineWidth(3);
    c1->cd();
    helix[998]->Draw();
    c2->cd();
    helix[998]->Draw();
    c3->cd();
    helix[998]->Draw();
    
    helix[997] = new THelix(0,0,0,wta2.Px(),wta2.Py(),wta2.Pz(),0.000001);
    if(wta2.Pz()<0) helix[997]->SetRange(-1,0);
    if(wta2.Pz()>=0) helix[997]->SetRange(0,1);
    helix[997]->SetLineColor(38);
    helix[997]->SetLineWidth(3);
    c1->cd();
    helix[997]->Draw();
    c2->cd();
    helix[997]->Draw();
    c3->cd();
    helix[997]->Draw();
  }

  int nHelix = 0;
  for(int i = 0; i<nParticle; i++){
    if(pwflag[i]!=0) continue;
    if(pt[i]<ptCut) continue;
    //jet tracks
    int trackColor = 0;
    if(jtpt[0]*TMath::CosH(jteta[0])>5 && dR(jteta[0],jtphi[0],eta[i],phi[i])<0.8) trackColor = 1;
    if(jtpt[1]*TMath::CosH(jteta[1])>5 && dR(jteta[1],jtphi[1],eta[i],phi[i])<0.8) trackColor = 2;
    if(jtpt[2]*TMath::CosH(jteta[2])>5 && dR(jteta[2],jtphi[2],eta[i],phi[i])<0.8) trackColor = 3;
    if(jtpt[3]*TMath::CosH(jteta[3])>5 && dR(jteta[3],jtphi[3],eta[i],phi[i])<0.8) trackColor = 4;
    if(doWTA && wtapt[i]<0.01) trackColor = 5;
    //std::cout <<  jteta[0] << " " << jtphi[0] <<" " <<  eta[i] << " " << phi[i] << std::endl;
    //std::cout <<  dR(jteta[0],jtphi[0],eta[i],phi[i]) << std::endl;

    //std::cout << i<<" "  << px[i]<<" "  << py[i]<<" " << pz[i]<<" " << charge[i] << " " << trackColor <<  std::endl; 
    if(doWTA){
      TVector3 trk3(0,0,0);
      trk3.SetPtEtaPhi(wtapt[i],wtaeta[i],wtaphi[i]);
      px[i] = trk3.x();
      py[i] = trk3.y();
      pz[i] = trk3.z();
    }

    helix[nHelix] = new THelix(0,0,0,px[i],py[i],pz[i],charge[i]*1.5);
    formatHelix(helix[nHelix],pz[i],pt[i],trackColor,doWTA);
    c1->cd();
    helix[nHelix]->Draw("same");
    c2->cd();
    helix[nHelix]->Draw("same");
    c3->cd();
    helix[nHelix]->Draw("same");
    nHelix++;
  }

}

int nextEvent(){
  EventDisplay(currentFile,currentEvtIndx+1,currentPtCut);
  return currentEvtIndx;
}
