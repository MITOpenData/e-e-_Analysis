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
  TH1F* dTTheta_h = new TH1F("dTTheta_h",";#Delta TTheta;Counts",600,-TMath::Pi(),TMath::Pi());
  TH1F* dTEta_h = new TH1F("dTEta_h",";#Delta TEta;Counts",1200,-10,10);
  TH1F* dTPhi_h = new TH1F("dTPhi_h",";#Delta TPhi;Counts",360,-TMath::Pi()/2.0,TMath::Pi()/2.0);

  centerTitles({dTTheta_h, dTEta_h,dTPhi_h});

  TFile* inFile_p = new TFile(inFileName.c_str(),"READ");
  TTree* genTree_p = (TTree*)inFile_p->Get("tgen");
  TTree* recoTree_p = (TTree*)inFile_p->Get("t");
  
  std::vector<std::string> particleList;
  std::vector<std::string> eventList;

  particleList.push_back("nParticle");
  eventList.push_back("TTheta");
  eventList.push_back("TPhi");

  particleData gen_particle; eventData gen_event; 
  particleData reco_particle; eventData reco_event;

  gen_particle.SetStatusAndAddressRead(genTree_p,particleList);   gen_event.SetStatusAndAddressRead(genTree_p,eventList);
  reco_particle.SetStatusAndAddressRead(recoTree_p,particleList);   reco_event.SetStatusAndAddressRead(recoTree_p,eventList);

  const Int_t nEntries = genTree_p->GetEntries();

  std::cout << "Processing " << nEntries << " events..." << std::endl;

  Float_t TTheta_diff;
  Float_t TEta_diff;
  Float_t TPhi_diff;

  TVector2 vec_recoPhi;                                                                                                                                                     
  TVector2 vec_recoPhi_rotatePi;
  TVector2 vec_genPhi;
  TVector2 vec_genPhi_rotatePi;

  for(Int_t entry = 0; entry < nEntries; ++entry)
  {
    if(entry%10000 == 0) std::cout << " Entry " << entry << "/" << nEntries << \
			   std::endl;
    genTree_p->GetEntry(entry);
    recoTree_p->GetEntry(entry);
    
    // std::cout<<gen_event.TTheta << " "<<reco_event.TTheta<<std::endl;
    //  std::cout<<gen_event.TPhi<< " "<<reco_event.TPhi<<std::endl;

    TTheta_diff = gen_event.TTheta - reco_event.TTheta;
    TEta_diff = -1*std::log(TMath::Tan(gen_event.TTheta/2.0)) - -1*std::log(TMath::Tan(reco_event.TTheta/2.0));
    
    vec_recoPhi.SetMagPhi(1.0,reco_event.TPhi);
    vec_recoPhi_rotatePi.SetMagPhi(1.0,reco_event.TPhi+TMath::Pi());
    vec_genPhi.SetMagPhi(1.0,gen_event.TPhi);
    vec_genPhi_rotatePi.SetMagPhi(1.0,gen_event.TPhi+TMath::Pi());
    //    Float_t tempDTPhi = vec_recoPhi.DeltaPhi(vec_genPhi);
    float ang1 = vec_recoPhi.DeltaPhi(vec_genPhi);
    float ang2 = vec_recoPhi.DeltaPhi(vec_genPhi_rotatePi);
    float ang3 = vec_recoPhi_rotatePi.DeltaPhi(vec_genPhi);
    float ang4 = vec_recoPhi_rotatePi.DeltaPhi(vec_genPhi_rotatePi);
    //std::cout<<ang1<<" "<<ang2<<" "<<ang3<<" "<<ang4<<std::endl;
    TPhi_diff = ang1;
    if(fabs(ang2)<fabs(TPhi_diff)) TPhi_diff = ang2;
    if(fabs(ang3)<fabs(TPhi_diff)) TPhi_diff = ang3;
    if(fabs(ang4)<fabs(TPhi_diff)) TPhi_diff = ang4;
    
    dTTheta_h->Fill(TTheta_diff);
    dTEta_h->Fill(TEta_diff);
    dTPhi_h->Fill(TPhi_diff);
  }
  
  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();
  dTTheta_h->Write("",TObject::kOverwrite);
  dTEta_h->Write("",TObject::kOverwrite);
  dTPhi_h->Write("",TObject::kOverwrite);

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
