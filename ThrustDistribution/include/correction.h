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
	/*
	Double_t bins[22] = {.6,.65,.7,.75,.8,.82,.84,.86,.88,.9,.92,.94,.95,.96,.965,.97,.975,.98,.985,.99,.995,1};
	Double_t corrs[21] = {1.054,0.910,1.013,0.944,1.055,1.034,1.100,1.105,1.070,1.079,1.104,1.104,1.083,1.096,1.048,1.055,1.007,0.926,0.803,0.664,0.425};
	Double_t errs[21] = {.192,.053,.040,.027,.041,.034,.034,.030,.023,.023,.018,.022,.019,.023,.021,.019,.019,.016,.016,.016,.022};
	*/
	
	// Create and fill correction data histogram
	//hist* corr = new hist("corr","corr",21,bins);
	 TH1F *corr = new TH1F("hGen1","Thrust",100,0.5,1);
        corr->SetBinContent(16,5);
        corr->SetBinContent(17,1.5);
        corr->SetBinContent(18,3.5);
        corr->SetBinContent(19,2);
        corr->SetBinContent(20,1.5);
        corr->SetBinContent(21,1.222222);
        corr->SetBinContent(22,1);
        corr->SetBinContent(23,1.666667);
        corr->SetBinContent(24,1.380952);
        corr->SetBinContent(25,0.9574468);
        corr->SetBinContent(26,1.342105);
        corr->SetBinContent(27,1.153846);
        corr->SetBinContent(28,1.191489);
        corr->SetBinContent(29,1.674157);
        corr->SetBinContent(30,1.311258);
        corr->SetBinContent(31,1.188119);
        corr->SetBinContent(32,1.391837);
        corr->SetBinContent(33,1.263158);
        corr->SetBinContent(34,1.141791);
        corr->SetBinContent(35,1.200436);
        corr->SetBinContent(36,1.112339);
        corr->SetBinContent(37,1.205527);
        corr->SetBinContent(38,1.092219);
        corr->SetBinContent(39,1.027778);
        corr->SetBinContent(40,1.081608);
        corr->SetBinContent(41,1.155131);
        corr->SetBinContent(42,1.061492);
        corr->SetBinContent(43,1.087584);
        corr->SetBinContent(44,1.082902);
        corr->SetBinContent(45,1.01672);
        corr->SetBinContent(46,1.082192);
        corr->SetBinContent(47,1.124528);
        corr->SetBinContent(48,1.107627);
        corr->SetBinContent(49,1.096583);
        corr->SetBinContent(50,1.08645);
        corr->SetBinContent(51,1.146988);
        corr->SetBinContent(52,1.161093);
        corr->SetBinContent(53,1.124805);
        corr->SetBinContent(54,1.115137);
        corr->SetBinContent(55,1.131768);
        corr->SetBinContent(56,1.163498);
        corr->SetBinContent(57,1.088826);
        corr->SetBinContent(58,1.186255);
        corr->SetBinContent(59,1.105679);
        corr->SetBinContent(60,1.122695);
        corr->SetBinContent(61,1.141532);
        corr->SetBinContent(62,1.0751);
        corr->SetBinContent(63,1.124241);
        corr->SetBinContent(64,1.111548);
        corr->SetBinContent(65,1.12716);
        corr->SetBinContent(66,1.105945);
        corr->SetBinContent(67,1.202689);
        corr->SetBinContent(68,1.149346);
        corr->SetBinContent(69,1.121861);
        corr->SetBinContent(70,1.137408);
        corr->SetBinContent(71,1.161931);
        corr->SetBinContent(72,1.178187);
        corr->SetBinContent(73,1.134632);
        corr->SetBinContent(74,1.123188);
        corr->SetBinContent(75,1.121468);
        corr->SetBinContent(76,1.163248);
        corr->SetBinContent(77,1.131413);
        corr->SetBinContent(78,1.138463);
        corr->SetBinContent(79,1.135695);
        corr->SetBinContent(80,1.141768);
        corr->SetBinContent(81,1.128625);
        corr->SetBinContent(82,1.123521);
        corr->SetBinContent(83,1.134837);
        corr->SetBinContent(84,1.142383);
        corr->SetBinContent(85,1.113584);
        corr->SetBinContent(86,1.146445);
        corr->SetBinContent(87,1.115501);
        corr->SetBinContent(88,1.130316);
        corr->SetBinContent(89,1.120865);
        corr->SetBinContent(90,1.122441);
        corr->SetBinContent(91,1.116792);
        corr->SetBinContent(92,1.103873);
        corr->SetBinContent(93,1.114389);
        corr->SetBinContent(94,1.112109);
        corr->SetBinContent(95,1.105644);
        corr->SetBinContent(96,1.036745);
        corr->SetBinContent(97,0.8128948);
        corr->SetBinContent(98,0.5197755);
        corr->SetBinContent(99,0.2701077);
        corr->SetBinContent(100,0.1466109);
        corr->SetEntries(102.6972);
	corr->Rebin(2);
	corr->Scale(1./2);
	/*
	for(int i=1;i<=21;i++)
	{
		corr->SetBinContent(i,corrs[i-1]);
		corr->SetBinError(i,errs[i-1]);
	}*/
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
