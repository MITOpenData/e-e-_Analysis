//c and cpp dependencies
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

//root dependencies
#include <TFile.h>
#include <TTree.h>
#include "TLorentzVector.h"
#include "TNamed.h"

//fastjet dependencies
#include "fastjet/ClusterSequence.hh"

//local dependencies
#include "include/checkMakeDir.h"
#include "include/particleData.h"
#include "include/jetData.h"

int scan(std::string inFileName="LEP2MCGGCCY1997E183_recons_aftercut-001.aleph")
{
  if(!checkFile(inFileName)){
    std::cout << "Given inFileName \'" << inFileName << "\' is invalid, return 1." << std::endl;
    return 1;
  }

  // "/data/flowex/Datasamples/LEP2_MAIN/ROOTfiles/cleaned_ALEPH_Data-all.aleph"
  FILE *fp=fopen(Form("%s",inFileName.c_str()),"r");
  float _px,_py,_pz,_m,_charge,_pwflag;
  
  TFile *hf = new TFile(Form("%s.root",inFileName.c_str()), "RECREATE");
  TTree *tout = new TTree("t","");
  TTree *jout = new TTree("jetTree","");

  particleData pData;
  jetData jData;

  tout->Branch("nParticle", &pData.nParticle,"nParticle/I");
  tout->Branch("EventNo", &pData.EventNo,"EventNo/I");
  tout->Branch("RunNo", &pData.RunNo,"RunNo/I");
  tout->Branch("Energy", &pData.Energy,"Energy/F");
  tout->Branch("px", pData.px,"px[nParticle]/F");
  tout->Branch("py", pData.py,"py[nParticle]/F");
  tout->Branch("pz", pData.pz,"pz[nParticle]/F");
  tout->Branch("pt", pData.pt,"pt[nParticle]/F");
  tout->Branch("pmag", pData.pmag,"pmag[nParticle]/F");//Added later on
  tout->Branch("mass", pData.mass,"mass[nParticle]/F");
  tout->Branch("eta", pData.eta,"eta[nParticle]/F");
  tout->Branch("theta", pData.theta,"theta[nParticle]/F");
  tout->Branch("phi", pData.phi,"phi[nParticle]/F");
  tout->Branch("charge", pData.charge,"charge[nParticle]/I");
  tout->Branch("pwflag", pData.pwflag,"pwflag[nParticle]/I");
  tout->Branch("pid", pData.pid,"pid[nParticle]/I");

  jout->Branch("nref", &jData.nref,"nref/I");
  jout->Branch("jtpt", jData.jtpt,"jtpt[nref]/F");
  jout->Branch("jteta", jData.jteta,"jteta[nref]/F");
  jout->Branch("jtphi", jData.jtphi,"jtphi[nref]/F");

  int counterEntries=0;
  int counterParticles=0;
  TLorentzVector v;

  std::vector<fastjet::PseudoJet> particles;
  const double rParam = 0.4;
  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rParam);

  while(fscanf(fp,"%f %f %f %f %f %f",&_px,&_py,&_pz,&_m,&_charge,&_pwflag)!=EOF){
    if(_px==-999. && _py==-999. && _pz==-999.){ 
      pData.nParticle=counterParticles;

      if(counterEntries>0) tout->Fill(); 


      //Processing particles->jets
      jData.nref = 0;
      if(particles.size() > 0){
	fastjet::ClusterSequence cs(particles, jetDef);
	std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());

	for(unsigned int j = 0; j < jets.size(); ++j){
	  jData.jtpt[jData.nref] = jets.at(j).pt();
	  jData.jtphi[jData.nref] = jets.at(j).phi_std();
	  jData.jteta[jData.nref] = jets.at(j).eta();
	  ++jData.nref;
	}
      }
      if(counterEntries>0) jout->Fill();

      //clear particles for next iteration
      particles.clear();

      pData.RunNo=_m; 
      pData.EventNo=_charge; 
      pData.Energy=_pwflag; 
      counterParticles=0;   

      continue;
    }  


    pData.px[counterParticles]=_px;
    pData.py[counterParticles]=_py;
    pData.pz[counterParticles]=_pz;
    pData.mass[counterParticles]=_m;
    v.SetXYZM(_px,_py,_pz,_m);
    particles.push_back(fastjet::PseudoJet(_px,_py,_pz,v.E()));
    pData.pt[counterParticles]=v.Pt();
    pData.pmag[counterParticles]=v.Rho(); //Added later on
    pData.eta[counterParticles]=v.PseudoRapidity();
    pData.theta[counterParticles]=v.Theta();
    pData.phi[counterParticles]=v.Phi();
    
    pData.charge[counterParticles]=_charge;
    pData.pwflag[counterParticles]=_pwflag;
    pData.pid[counterParticles]=0;
    ++counterParticles;	
    ++counterEntries;	
  }

  hf->cd();

  tout->Write("", TObject::kOverwrite);
  delete tout;
  jout->Write("", TObject::kOverwrite);
  delete jout;

  hf->Close();
  delete hf;

  return 0;
}

int main(int argc, char *argv[])
{
  if(argc != 1 && argc != 2){
    std::cout << "Usage: ./bin/scan.exe <OPT-inFileName>" << std::endl;
    return 1;
  }

  std::cout << "Begin processing..." << std::endl;
  int retVal = 0;
  if(argc == 1) retVal += scan();
  else if(argc == 2) retVal += scan(argv[1]);
  return retVal;
}
