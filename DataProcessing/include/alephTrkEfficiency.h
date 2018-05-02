//Based on efficiency.h by Michael Peters -> updated to be a class
#ifndef ALEPHTRKEFFICIENCY_H
#define ALEPHTRKEFFICIENCY_H

#include <iostream>

#include <TFile.h>
#include <TH3F.h>

class alephTrkEfficiency{
 public:
  TFile *_effInf=NULL;
  TH3F *_heff=NULL;
  
  alephTrkEfficiency();
  ~alephTrkEfficiency();

  Float_t efficiency(Float_t theta, Float_t phi, Float_t pt);
  Int_t counter;
};

alephTrkEfficiency::alephTrkEfficiency()
{
  _effInf = new TFile("DataProcessing/tables/efficiency_hist.root","READ");
  _heff = (TH3F*)_effInf->Get("eff");
  counter = 0;
  return;
}

alephTrkEfficiency::~alephTrkEfficiency(){if(_effInf != NULL) delete _effInf;}

Float_t alephTrkEfficiency::efficiency(Float_t theta, Float_t phi, Float_t pt)
{
  Float_t e = _heff->GetBinContent(_heff->FindBin(pt,theta,phi));
  if(e < 0.00000000001){
    if(counter < 10){
      std::cout << "!!!Error on efficiency correction! Zero efficiency!!! theta=" << theta << " phi=" << phi << " pt=" << pt << std::endl << std::endl;
      if(counter == 9) std::cout << " !!!Greater than ten calls of this error... TERMINATING OUTPUT, PLEASE FIX" << std::endl;
      ++counter;
    }
    
    e = 1;
  }  
  return e;
}

#endif
