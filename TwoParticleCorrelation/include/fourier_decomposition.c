
#include "fourier.h"
using namespace std;

#define PI 3.1415926
//enum SIMPLEPID {PHOTON, ELECTRON, PION, MUON, KAON, PROTON};
//enum SIMPLEPWLAG {CHARGED_TRACK, CHARGED_LEPTONS1, CHARGED_LEPTONS2, V0, PHOTON, NEUTRAL_HADRON};


void subtract()
{
  TH1D* hist1;
  TH1D* hist2;
  
  analysis(0,0, 0, 30, 40, 50, 0, 1, hist1);
  
  analysis ( 0, 0,  0, 0, 10, 50, 0, 1, hist2);
  
  
  cout<<"hi"<<" "<<"1"<<endl;
  //hist1->FillRandom("gaus", 10000);
  //hist2->FillRandom("gaus", 10000);
  
  cout<<"hi"<<" "<<"2"<<endl;
  hist1->Add(hist2, -1);
  
  TCanvas *c1 = new TCanvas("c1","",600,600);
  c1->cd();
  cout<<"hi"<<" "<<"3"<<endl;
  hist1->Draw();
  
  c1->SaveAs("fit.pdf");
  cout<<"hi"<<" "<<"4"<<endl;
}