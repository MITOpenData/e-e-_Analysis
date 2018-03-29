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

// Dataformat
#include "include/smartJetName.h"
#include "include/particleData.h"
#include "include/jetData.h"
#include "include/boostedEvtData.h"
#include "include/eventData.h"

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
  if(!doWTA) && pt<2.5 && TMath::Abs(pz)<0.5) rangeBound = TMath::Abs(pz);
  h->SetRange(0,rangeBound);
  if(pz<0) h->SetRange(-rangeBound,0);
}

void EventDisplay(std::string inputFile = "/data/cmcginn/StudyMultSamples/ALEPH/LEP1/20180125/LEP1Data1995_recons_aftercut-MERGED.root" ,
                  int eventIndx = 574, 
		  float ptCut = 0, 
		  bool doWTA = false){
  currentFile=inputFile;
  currentEvtIndx=eventIndx;
  currentPtCut=ptCut;

  TFile * f 	= TFile::Open(inputFile.c_str(),"read");
  TTree * t 	= (TTree*)f->Get("t");
  TTree * WTA_t = (TTree*)f->Get("BoostedWTAR8Evt");
  TTree * jt 	= (TTree*)f->Get(smartJetName("ak8ESchemeJetTree", f).c_str()); 
  TTree * wta 	= (TTree*)f->Get(smartJetName("ak8WTAmodpSchemeJetTree", f).c_str());

  // Define dataformat
  particleData particle;
  jetData jet;
  jetData WTAjet;
  boostedEvtData boosted;
  eventData event;

  // Setup branches
  std::vector<std::string> list;  
  particle.SetStatusAndAddressRead(t,list);
  event.SetStatusAndAddressRead(t,list);
  boosted.SetStatusAndAddressRead(WTA_t,list);
  jet.SetStatusAndAddressRead(jt,list);
  WTAjet.SetStatusAndAddressRead(jt,list);

  // Take the event of interest
  t->GetEntry(eventIndx);
  if(doWTA) WTA_t->GetEntry(eventIndx);
  jt->GetEntry(eventIndx);
  wta->GetEntry(eventIndx);

  // Draw Thrust and WTA axis
  TVector3 thrust = TVector3(0,0,0);
  thrust.SetMagThetaPhi(10,event.TTheta,event.TPhi);
  TVector3 wta1 = TVector3(0,0,0);
  if(WTAjet.nref>0) wta1.SetMagThetaPhi(10,2*TMath::ATan(TMath::Exp(-WTAjet.jteta[0])),WTAjet.jtphi[0]);
  TVector3 wta2 = TVector3(0,0,0);
  if(WTAjet.nref>1) wta2.SetMagThetaPhi(10,2*TMath::ATan(TMath::Exp(-WTAjet.jteta[1])),WTAjet.jtphi[1]);


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
  view3->RotateView(event.TPhi*180/TMath::Pi(),event.TTheta*180/TMath::Pi());
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

  // Draw charged particles
  int nHelix = 0;
  for(int i = 0; i<particle.nParticle; i++){
    if(particle.pwflag[i]!=0) continue;
    if(particle.pt[i]<ptCut) continue;

    //jet tracks
    int trackColor = 0;
    if(jet.jtpt[0]*TMath::CosH(jet.jteta[0])>5 && dR(jet.jteta[0],jet.jtphi[0],particle.eta[i],particle.phi[i])<0.8) trackColor = 1;
    if(jet.jtpt[1]*TMath::CosH(jet.jteta[1])>5 && dR(jet.jteta[1],jet.jtphi[1],particle.eta[i],particle.phi[i])<0.8) trackColor = 2;
    if(jet.jtpt[2]*TMath::CosH(jet.jteta[2])>5 && dR(jet.jteta[2],jet.jtphi[2],particle.eta[i],particle.phi[i])<0.8) trackColor = 3;
    if(jet.jtpt[3]*TMath::CosH(jet.jteta[3])>5 && dR(jet.jteta[3],jet.jtphi[3],particle.eta[i],particle.phi[i])<0.8) trackColor = 4;
    if(doWTA && boosted.pt[i]<0.01) trackColor = 5;
    //std::cout <<  jteta[0] << " " << jtphi[0] <<" " <<  eta[i] << " " << phi[i] << std::endl;
    //std::cout <<  dR(jteta[0],jtphi[0],eta[i],phi[i]) << std::endl;

    //std::cout << i<<" "  << px<<" "  << py<<" " << pz<<" " << particle.charge[i] << " " << trackColor <<  std::endl; 

    Float_t px,py,pz,pt;
    if(doWTA){
      TVector3 trk3(0,0,0);
      trk3.SetPtEtaPhi(boosted.pt[i],boosted.eta[i],boosted.phi[i]);
      px = trk3.x();
      py = trk3.y();
      pz = trk3.z();
      pt = boosted.pt[i];
    } else {
      px = particle.px[i];
      py = particle.py[i];
      pz = particle.pz[i];
      pt = particle.pt[i];
    }

    helix[nHelix] = new THelix(0,0,0,px,py,pz,particle.charge[i]*1.5);
    formatHelix(helix[nHelix],pz,pt,trackColor,doWTA);
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
