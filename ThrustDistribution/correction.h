#ifndef correction_h
#define correction_h

#include "TH1.h"
#include "TMath.h"

template <class hist>
// Template works for any of the TH1 types.
hist* correct_hist(hist* h)
{
	// Determine number of bins above 0.6
	Int_t nbins = h->GetNbinsX();
	Int_t nbins_crop=0;
	Double_t first_edge;
	for(int i=1;i<=nbins;i++) 
	{
		if(h->GetBinLowEdge(i)>0.6) nbins_crop++;
		if(h->GetBinLowEdge(i-1)<0.6) first_edge = h->GetBinLowEdge(i);
	}
	// Create and fill new "cropped" histogram that starts at 0.6
	hist* h_crop = new hist("thrust_cropped","thrust_cropped",nbins_crop,first_edge,1);
	for(int i=1;i<=nbins_crop;i++) 
	{
		h_crop->SetBinContent(i,h->GetBinContent(i+(nbins-nbins_crop)));
		h_crop->SetBinError(i,h->GetBinError(i+(nbins-nbins_crop)));
	}
	// Digitized correction data below
	Double_t bins[22] = {.6,.65,.7,.75,.8,.82,.84,.86,.88,.9,.92,.94,.95,.96,.965,.97,.975,.98,.985,.99,.995,1};
	Double_t corrs[21] = {1.054,0.910,1.013,0.944,1.055,1.034,1.100,1.105,1.070,1.079,1.104,1.104,1.083,1.096,1.048,1.055,1.007,0.926,0.803,0.664,0.425};
	Double_t errs[21] = {.192,.053,.040,.027,.041,.034,.034,.030,.023,.023,.018,.022,.019,.023,.021,.019,.019,.016,.016,.016,.022};
	// Create and fill correction data histogram
	hist* corr = new hist("corr","corr",21,bins);
	for(int i=1;i<=21;i++)
	{
		corr->SetBinContent(i,corrs[i-1]);
		corr->SetBinError(i,errs[i-1]);
	}
	// Apply corrections.
	// Differences in bin size necessitate the use of interpolation.
	// Uncertainty is calculated by propagating the uncertainty in the product of thrust and
	// correction factor (where uncertainty in correction factor is taken to be the
	// uncertainty of the nearest data point).
	for(int i=1;i<=nbins_crop;i++)
	{
		
		Double_t bcorr = corr->Interpolate(h_crop->GetBinCenter(i));
		h_crop->SetBinContent(i,h_crop->GetBinContent(i)*bcorr);
		Int_t corrbin = corr->FindBin(h_crop->GetBinCenter(i));
		h_crop->SetBinError(i,h_crop->GetBinContent(i)*TMath::Sqrt(pow(h_crop->GetBinError(i)/h_crop->GetBinContent(i),2)+pow(corr->GetBinError(corrbin)/corr->GetBinContent(corrbin),2)));
	}
	delete(corr);
	delete(h);
	return h_crop;
}

Double_t correct_entry(Double_t t)
{
	Double_t bins[22] = {.6,.65,.7,.75,.8,.82,.84,.86,.88,.9,.92,.94,.95,.96,.965,.97,.975,.98,.985,.99,.995,1};
	Double_t corrs[21] = {1.054,0.910,1.013,0.944,1.055,1.034,1.100,1.105,1.070,1.079,1.104,1.104,1.083,1.096,1.048,1.055,1.007,0.926,0.803,0.664,0.425};
	Double_t errs[21] = {.192,.053,.040,.027,.041,.034,.034,.030,.023,.023,.018,.022,.019,.023,.021,.019,.019,.016,.016,.016,.022};
	TH1D* corr = new TH1D("corr","corr",21,bins);
	for(int i=1;i<=21;i++)
	{
		corr->SetBinContent(i,corrs[i-1]);
		corr->SetBinError(i,errs[i-1]);
	}
	Double_t res = corr->Interpolate(t);
	delete(corr);
	return res;
}
#endif
