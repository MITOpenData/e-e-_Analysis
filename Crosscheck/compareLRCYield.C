#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLegend.h"
#include "TAttMarker.h"
#include <iostream>
#include <string>

void formatTH1D(TH1D * h, float offsetx, float offsety){
   h->SetTitleOffset(offsetx,"X");
   h->SetTitleOffset(offsety,"Y");
   h->GetXaxis()->CenterTitle();
   h->GetYaxis()->CenterTitle();
   std::cout << h->GetMaximum()*1.2 << " " << h->GetMinimum()*0.8 << std::endl;
   h->GetYaxis()->SetRangeUser(h->GetMaximum()*1.2,h->GetMinimum()*0.8);
   h->SetStats(0);
}

void compareLRCYield(std::string first, std::string second, std::string plot, std::string legText, std::string saveText, std::string axis){
  TH1::SetDefaultSumw2();

  TFile * f1 = TFile::Open(first.c_str(),"read");
  TH1D * h1 = (TH1D*) f1->Get(plot.c_str());
  TFile * f2 = TFile::Open(second.c_str(),"read");
  TH1D * h2 = (TH1D*) f2->Get(plot.c_str());  

  TCanvas * c1 = new TCanvas("c1","c1",800,600);
  c1->SetLeftMargin(0.2);
  formatTH1D(h1,1.0,1.5); 
  h1->SetMarkerStyle(8);
  h1->SetMarkerColor(kBlack);  
  h1->SetLineColor(kBlack);  
  h1->Draw("p");

  h2->SetMarkerStyle(25);
  h2->SetMarkerColor(kRed);
  h2->SetLineColor(kRed);
  h2->Draw("same p");

  TLegend * l = new TLegend(0.25,0.55,0.55,0.85);
  l->AddEntry((TObject*)0, legText.c_str(),"");
  l->AddEntry((TObject*)0, axis.c_str(),"");
  l->AddEntry(h1,"Pythia 8","p");
  l->AddEntry(h2,"Pythia 8 w/ Ropewalk","p");
  l->Draw("same");

  c1->SaveAs(Form("img/LRCYieldComparisonPlots/%s.png",saveText.c_str()));
  c1->SaveAs(Form("img/LRCYieldComparisonPlots/%s.pdf",saveText.c_str()));
  c1->SaveAs(Form("img/LRCYieldComparisonPlots/%s.C",saveText.c_str()));

}

void runCompare(){
  compareLRCYield("img/Dec14/beamline_pythia8/Analyzer_Output.root","img/Dec14/beamline_pythia8_rope/Analyzer_Output_Rope.root","longRangeYield_0_20","Multiplicity: [0,20)","pythia_beam_0_20","Beam Axis");
  compareLRCYield("img/Dec14/beamline_pythia8/Analyzer_Output.root","img/Dec14/beamline_pythia8_rope/Analyzer_Output_Rope.root","longRangeYield_20_30","Multiplicity: [20,30)","pythia_beam_20_30","Beam Axis");
  compareLRCYield("img/Dec14/beamline_pythia8/Analyzer_Output_HM.root","img/Dec14/beamline_pythia8_rope/Analyzer_Output_RopeHM.root","longRangeYield_30_999","Multiplicity: [30,999)","pythia_beam_30_999","Beam Axis");


  compareLRCYield("img/Dec14/thrust_pythia8/Analyzer_Output.root","img/Dec14/thrust_pythia8_rope/Analyzer_Output_Rope.root","longRangeYield_0_20","Multiplicity: [0,20)","pythia_thrust_0_20","Thrust Axis");
  compareLRCYield("img/Dec14/thrust_pythia8/Analyzer_Output.root","img/Dec14/thrust_pythia8_rope/Analyzer_Output_Rope.root","longRangeYield_20_30","Multiplicity: [20,30)","pythia_thrust_20_30","Thrust Axis");
  compareLRCYield("img/Dec14/thrust_pythia8/Analyzer_Output_HM.root","img/Dec14/thrust_pythia8_rope/Analyzer_Output_RopeHM.root","longRangeYield_30_999","Multiplicity: [30,999)","pythia_thrust_30_999","Thrust Axis");
}
