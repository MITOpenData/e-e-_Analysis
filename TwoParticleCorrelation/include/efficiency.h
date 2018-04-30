#include <TH3F.h>
#include <TFile.h>

TFile *_effInf;
TH3F *_heff;

Float_t efficiency(Float_t theta, Float_t phi, Float_t pt)
{
        static int i=0;
	if (i==0){
	   _effInf = new TFile("../include/efficiency_hist.root","read");
   	   _heff = (TH3F*)_effInf->Get("eff");
	}
	Float_t e = _heff->GetBinContent(_heff->FindBin(pt,theta,phi));
        if (e==0) cout <<"!!!Error on efficiency correction! Zero efficiency!!!"<<endl<<endl;
	  
	return e;
}