#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TStyle.h>
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "getLogBins.h"
#include "thrust_distribution.c"

void thrust_vs_hep(TString thrust_file = "~/Downloads/StudyMult/TwoParticleCorrelation/alephDataPaths_LEP1.root",
TString hep_file = "~/Downloads/HEPData-ins636645-v1-Table54.root")
{
	TFile *hdata = new TFile(hep_file);
	TH1F *hep;
	hdata->cd("Table 54");
	hep = (TH1F*)gDirectory->Get("Hist1D_y1");
	thrust_distribution(thrust_file);
	TH1F *thrust = (TH1F*)gPad->GetPrimitive("h_thrust");
	hep->GetYaxis()->SetRange(0,22);	
	hep->SetLineColor(kRed);
	hep->Draw("HIST same");
	TLegend *leg = new TLegend(.15,.15,.4,.4);
	leg->AddEntry(thrust,"Badea","l");
	leg->AddEntry(hep,"ALEPH","l");
	leg->Draw();
}


