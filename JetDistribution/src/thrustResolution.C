// Computing delta eta and phi between gen and reco for thrust resolution 
// Original author Anthony Badea

// cpp dependencies
#include <iostream>
#include <string>

// ROOT dependencies
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TDatime.h"
#include "TMath.h"
#include "TVector2.h"

// Non-Local StudyMult dependencies
#include "DataProcessing/include/checkMakeDir.h"
#include "DataProcessing/include/histDefUtility.h"
#include "DataProcessing/include/particleData.h"
#include "DataProcessing/include/eventData.h"

// Local StudyMult (JetDistribution) dependencies

int thrustResolution(const std::string inFileName)
{
  if(!checkFile(inFileName)){
    std::cout << "Given inFileName \'" << inFileName << "\' is invalid. return \
1" << std::endl;
    return 1;
  }
  else if(inFileName.find(".root") == std::string::npos){
    std::cout << "Given inFileName \'" << inFileName << "\' is not a root file.\
 return 1" << std::endl;
    return 1;
  }

  checkMakeDir("output");

  TDatime* date = new TDatime();
  std::string outFileName = inFileName;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  outFileName.replace(outFileName.find(".root"), 5, "");
  outFileName = "output/" + outFileName + "_ThrustResolution_" + std::to_string(date->GetDate()) + ".root";
  delete date;

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  
  std::cout<<"Initializing histograms and loading data..."<<std::endl;
  /************************************************/
  // Initialize Histograms                         /
  /************************************************/
  
  // Subranges of NtrkOffline
  static const int numRanges = 5;
  int nTrkMin[numRanges] = {4,10,20,30,35};
  int nTrkMax[numRanges] = {10,20,30,999,999};
  int colors[numRanges] = {0,1,2,3,4}; // gold,blue,red,orange,green,peach
    
  TH1F* dTTheta_h = new TH1F("dTTheta_h",";#Delta TTheta;Counts",600,-TMath::Pi(),TMath::Pi());
  TH1F* dTEta_h = new TH1F("dTEta_h",";#Delta TEta;Counts",1200,-10,10);
  TH1F* dTPhi_h = new TH1F("dTPhi_h",";#Delta TPhi;Counts",360,-TMath::Pi()/2.0,TMath::Pi()/2.0);
  centerTitles({dTTheta_h, dTEta_h,dTPhi_h});
    
  // nTrk cuts //
  TH1F * dTTheta_Cut[numRanges];
  TH1F * dTEta_Cut[numRanges];
  TH1F * dTPhi_Cut[numRanges];
    
  for(unsigned int i = 0; i < numRanges; ++i)
  {
      dTTheta_Cut[i] = new TH1D(Form("dTTheta%d",nTrkMin[i]), Form("nTrkOffline_%d;N_{Trk}^{Offline};Entries",nTrkMin[i]), 600,-TMath::Pi(),TMath::Pi());
      dTEta_Cut[i] = new TH1F(Form("dTEta_%d",nTrkMin[i]), Form("Thrust_%d;Thrust;Entries",1200,-10,10);
      dTPhi_Cut[i] = new TH1F(Form("dTPhi_%d",nTrkMin[i]), Form("nJets_%d;N_{Jets}^{Offline};Entries",nTrkMin[i]),360,-TMath::Pi()/2.0,TMath::Pi()/2.0);
                              
      centerTitles({dTTheta_Cut[i], dTEta_Cut[i],dTPhi_Cut[i]});
  }

  /************************************************/
  // Initialize Tree and Variables                 /
  /************************************************/
                              
  TFile* inFile_p = new TFile(inFileName.c_str(),"READ");
  TTree* genTree_p = (TTree*)inFile_p->Get("tgen");
  TTree* recoTree_p = (TTree*)inFile_p->Get("t");
  
  std::vector<std::string> particleList;
  std::vector<std::string> eventList;

  particleList.push_back("nChargedHadronsHP");
  eventList.push_back("TTheta");
  eventList.push_back("TPhi");

  particleData gen_particle; eventData gen_event; 
  particleData reco_particle; eventData reco_event;

  gen_particle.SetStatusAndAddressRead(genTree_p,particleList);   gen_event.SetStatusAndAddressRead(genTree_p,eventList);
  reco_particle.SetStatusAndAddressRead(recoTree_p,particleList);   reco_event.SetStatusAndAddressRead(recoTree_p,eventList);

  const Int_t nEntries = genTree_p->GetEntries();

  std::cout << "Processing " << nEntries << " events..." << std::endl;

  Float nChargedHadronsHP; // use gen particles to do this
                              
  Float_t TTheta_diff;
  Float_t TEta_diff;
  Float_t TPhi_diff;

  TVector2 vec_recoPhi;                                                                                                                                                     
  TVector2 vec_recoPhi_rotatePi;
  TVector2 vec_genPhi;
  TVector2 vec_genPhi_rotatePi;

  for(Int_t entry = 0; entry < nEntries; ++entry)
  {
    if(entry%10000 == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;
    genTree_p->GetEntry(entry);
    recoTree_p->GetEntry(entry);
    nChargedHadronsHP = gen_particle.nChargedHadronsHP;
      
    TTheta_diff = gen_event.TTheta - reco_event.TTheta;
    TEta_diff = -1*std::log(TMath::Tan(gen_event.TTheta/2.0)) - -1*std::log(TMath::Tan(reco_event.TTheta/2.0));
    
    vec_recoPhi.SetMagPhi(1.0,reco_event.TPhi);
    vec_recoPhi_rotatePi.SetMagPhi(1.0,reco_event.TPhi+TMath::Pi());
    vec_genPhi.SetMagPhi(1.0,gen_event.TPhi);
    vec_genPhi_rotatePi.SetMagPhi(1.0,gen_event.TPhi+TMath::Pi());

    float ang1 = vec_recoPhi.DeltaPhi(vec_genPhi);
    float ang2 = vec_recoPhi.DeltaPhi(vec_genPhi_rotatePi);
    float ang3 = vec_recoPhi_rotatePi.DeltaPhi(vec_genPhi);
    float ang4 = vec_recoPhi_rotatePi.DeltaPhi(vec_genPhi_rotatePi);

    TPhi_diff = ang1;
    if(fabs(ang2)<fabs(TPhi_diff)) TPhi_diff = ang2;
    if(fabs(ang3)<fabs(TPhi_diff)) TPhi_diff = ang3;
    if(fabs(ang4)<fabs(TPhi_diff)) TPhi_diff = ang4;
    
    // no cut
    dTTheta_h->Fill(TTheta_diff);
    dTEta_h->Fill(TEta_diff);
    dTPhi_h->Fill(TPhi_diff);
      
    // 4-10
    if(nChargedHadronsHP >= nTrkMin[0] && nChargedHadronsHP < nTrkMax[0])
    {
        dTTheta_Cut[0]->Fill(TTheta_diff);
        dTEta_Cut[0]->Fill(TEta_diff);
        dTPhi_Cut[0]->Fill(TPhi_diff);
    }
    // 10-20
    else if(nChargedHadronsHP >= nTrkMin[1] && nChargedHadronsHP < nTrkMax[1])
    {
        dTTheta_Cut[1]->Fill(TTheta_diff);
        dTEta_Cut[1]->Fill(TEta_diff);
        dTPhi_Cut[1]->Fill(TPhi_diff);
    }
    // 20-30
    else if(nChargedHadronsHP >= nTrkMin[2] && nChargedHadronsHP < nTrkMax[2])
    {
        dTTheta_Cut[2]->Fill(TTheta_diff);
        dTEta_Cut[2]->Fill(TEta_diff);
        dTPhi_Cut[2]->Fill(TPhi_diff);
    }
    // 30-999
    else if(nChargedHadronsHP >= nTrkMin[3] && nChargedHadronsHP < nTrkMax[3])
    {
        dTTheta_Cut[3]->Fill(TTheta_diff);
        dTEta_Cut[3]->Fill(TEta_diff);
        dTPhi_Cut[3]->Fill(TPhi_diff);
        // 35-999
        if (nChargedHadronsHP >= nTrkMin[4] && nChargedHadronsHP < nTrkMax[4])
        {
            dTTheta_Cut[4]->Fill(TTheta_diff);
            dTEta_Cut[4]->Fill(TEta_diff);
            dTPhi_Cut[4]->Fill(TPhi_diff);
        }
    }
  }
  
  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();
  dTTheta_h->Write("",TObject::kOverwrite);
  dTEta_h->Write("",TObject::kOverwrite);
  dTPhi_h->Write("",TObject::kOverwrite);
                              
  for(unsigned int i = 0; i < numRanges; ++i)
  {
      dTTheta_Cut[i]->Write("",TObject::kOverwrite);
      dTEta_Cut[i]->Write("",TObject::kOverwrite);
      dTPhi_Cut[i]->Write("",TObject::kOverwrite);
      
      delete dTTheta_Cut[i];
      delete dTPhi_Cut[i];
      delete dTEta_Cut[i];
  }
                                                          
  delete dTTheta_h;
  delete dTPhi_h;
  delete dTEta_h;

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/thrustResolution.exe <inName>" << std::endl\
      ;
    return 1;
  }

  int retVal = 0;
  retVal += thrustResolution(argv[1]);
  return retVal;
}
