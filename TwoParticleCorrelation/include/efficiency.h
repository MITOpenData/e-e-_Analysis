#include "TH3.h"

Float_t efficiency(Float_t theta, Float_t phi, Float_t pt)
{
	TFile* f = TFile::Open("efficiency_hist.root","read");
	TH3F* eff = (TH3F*)f->Get("eff");
	Float_t e = eff->GetBinContent(eff->FindBin(pt,theta,phi));
	f->Close();
	return e;
}