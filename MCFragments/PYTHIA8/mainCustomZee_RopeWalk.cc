// main01.cc is a part of the PYTHIA event generator.
// Copyright (C) 2017 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk.
// It studies the charged multiplicity distribution at the LHC.

#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"

//fastjet dependencies
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

//local dependencies
#include "include/jetData.h"
#include "include/particleData.h"
#include "include/thrustTools.h"

void processJets(std::vector<fastjet::PseudoJet> p, fastjet::JetDefinition jDef, fastjet::JetDefinition jDefReclust, jetData *d, const double ptCut = 0.1)
{
  d->nref = 0;
  if(p.size() > 0){
    fastjet::ClusterSequence cs(p, jDef);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
    for(unsigned int i = 0; i < jets.size(); ++i){
      if(jets.at(i).pt() < ptCut) break; //Arbitrarily low cut on jets, removes spike at phi zero when things become ill defined         
      d->jtpt[d->nref] = jets.at(i).pt();
      d->jtphi[d->nref] = jets.at(i).phi_std();
      d->jteta[d->nref] = jets.at(i).eta();


      std::vector<fastjet::PseudoJet> jetConst = jets.at(i).constituents();
      d->jtN[d->nref] = jetConst.size();
      std::vector< std::vector<fastjet::PseudoJet> > subJets;

      for(int j = 0; j < 6; ++j){
        d->jtNPW[d->nref][j] = 0;
        d->jtptFracPW[d->nref][j] = 0;

	std::vector<fastjet::PseudoJet> tempSubJets;
        subJets.push_back(tempSubJets);
      }

      for(unsigned int k = 0; k < jetConst.size(); ++k){
        if(jetConst.at(k).user_index() < 0 || jetConst.at(k).user_index() > 5) continue;
        subJets.at(jetConst.at(k).user_index()).push_back(jetConst.at(k));
        d->jtNPW[d->nref][jetConst.at(k).user_index()]++;
      }

      for(int j = 0; j < 6; ++j){
	fastjet::ClusterSequence csSub(subJets.at(j), jDefReclust);

	std::vector<fastjet::PseudoJet> constTot = fastjet::sorted_by_pt(csSub.inclusive_jets());
        if(constTot.size() > 1){
	  std::cout << "WARNING - RECLUSTER OF CONSTITUENTS YIELDS GREATER THAN 1 JET" << std::endl;
	  std::cout << "Top jet: " << jets.at(i).pt() << ", " << jets.at(i).phi() << ", " << jets.at(i).eta() << std::endl;
	  std::cout << "Const: " << std::endl;
          for(unsigned int k = 0; k < subJets.at(j).size(); ++k){
	    std::cout << " " << k << "/" << subJets.at(j).size() << ": " << subJets.at(j).at(k).pt() << ", " << subJets.at(j).at(k).phi() << ", " << subJets.at(j).at(k).eta() << std::endl;
          }

	  std::cout << "Reclust: " << std::endl;
          for(unsigned int k = 0; k < constTot.size(); ++k){
	    std::cout << " " << k << "/" << constTot.size() << ": " << constTot.at(k).pt() << ", " << constTot.at(k).phi() << ", " << constTot.at(k).eta() << std::endl;
          }

        }
        else if(constTot.size() == 1) d->jtptFracPW[d->nref][j] = constTot.at(0).pt()/jets.at(i).pt();
        else d->jtptFracPW[d->nref][j] = 0.;
      }

      ++d->nref;
    }
    jets.clear();
    std::vector<fastjet::PseudoJet>().swap(jets);
  }
  return;
}


int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./mainCustom.exe <inFileName>" << std::endl;
    return 1;
  }
  
  const double jtPtCut = .01;
  const int nJtAlgo = 4;
  const double rParam[nJtAlgo] = {0.4, 0.4, 0.8, 0.8};
  const double recombScheme[nJtAlgo] = {fastjet::E_scheme, fastjet::WTA_modp_scheme, fastjet::E_scheme, fastjet::WTA_modp_scheme};
  fastjet::JetDefinition jDef[nJtAlgo];
  fastjet::JetDefinition jDefReclust[nJtAlgo];
  std::string jetTreeName[nJtAlgo];
  for(int i = 0; i < nJtAlgo; ++i){
    std::string recombSchemeStr = "EScheme";
    if(recombScheme[i] == fastjet::WTA_modp_scheme) recombSchemeStr = "WTAmodpScheme";

    jetTreeName[i] = "ak" + std::to_string(int(rParam[i]*10)) + recombSchemeStr + "JetTree";
  }

  for(int i = 0; i < nJtAlgo; ++i){
    jDef[i] = fastjet::JetDefinition(fastjet::ee_genkt_algorithm, rParam[i], -1, fastjet::RecombinationScheme(recombScheme[i]));
    jDefReclust[i] = fastjet::JetDefinition(fastjet::ee_genkt_algorithm, 5, -1, fastjet::RecombinationScheme(recombScheme[i]));
  }
  
  
  TFile* inFile_p = new TFile(argv[1], "RECREATE");
  TTree* genTree_p = new TTree("t", "t");
  TTree* jetTree_p[nJtAlgo];
  for(int i = 0; i < nJtAlgo; ++i){jetTree_p[i] = new TTree(jetTreeName[i].c_str(), jetTreeName[i].c_str());}

  jetData jData[nJtAlgo];

  Float_t pthat_;
  Float_t pthatWeight_;

  const int nMaxPart = 100000;
  Int_t nParticle_;
  Float_t px_[nMaxPart];
  Float_t py_[nMaxPart];
  Float_t pz_[nMaxPart];
  Float_t mass_[nMaxPart];
  Float_t theta_[nMaxPart];
  Float_t pt_[nMaxPart];
  Float_t phi_[nMaxPart];
  Float_t eta_[nMaxPart];
  Int_t pid_[nMaxPart];
  Int_t pwflag_[nMaxPart];

  Float_t TTheta_;
  Float_t TPhi_;

  genTree_p->Branch("pthat", &pthat_, "pthat/F");
  genTree_p->Branch("pthatWeight", &pthatWeight_, "pthatWeight/F");
  genTree_p->Branch("nParticle", &nParticle_, "nParticle/I");
  genTree_p->Branch("px", px_, "px[nParticle]/F");
  genTree_p->Branch("py", py_, "py[nParticle]/F");
  genTree_p->Branch("pz", pz_, "pz[nParticle]/F");
  genTree_p->Branch("mass", mass_, "mass[nParticle]/F");
  genTree_p->Branch("theta", theta_, "theta[nParticle]/F");
  genTree_p->Branch("pt", pt_, "pt[nParticle]/F");
  genTree_p->Branch("phi", phi_, "phi[nParticle]/F");
  genTree_p->Branch("eta", eta_, "eta[nParticle]/F");
  genTree_p->Branch("pid", pid_, "pid[nParticle]/I");
  genTree_p->Branch("pwflag", pwflag_, "pwflag[nParticle]/I");
  genTree_p->Branch("TTheta", &TTheta_, "TTheta/F");
  genTree_p->Branch("TPhi", &TPhi_, "TPhi/F");

  for(int i = 0; i < nJtAlgo; ++i){
    jetTree_p[i]->Branch("nref", &jData[i].nref,"nref/I");
    jetTree_p[i]->Branch("jtpt", jData[i].jtpt,"jtpt[nref]/F");
    jetTree_p[i]->Branch("jteta", jData[i].jteta,"jteta[nref]/F");
    jetTree_p[i]->Branch("jtphi", jData[i].jtphi,"jtphi[nref]/F");
    jetTree_p[i]->Branch("jtN", jData[i].jtN, "jtN[nref]/I");
    jetTree_p[i]->Branch("jtNPW", jData[i].jtNPW, "jtNPW[nref][6]/I");
    jetTree_p[i]->Branch("jtptFracPW", jData[i].jtptFracPW, "jtptFracPW[nref][6]/F");
  }



  Pythia8::Pythia pythia;
  double mZ = pythia.particleData.m0(23);
  pythia.settings.parm("Beams:eCM", mZ);
  pythia.readString("Beams:idA = 11");
  pythia.readString("Beams:idB = -11");
  pythia.readString("PDF:lepton = off");
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  pythia.readString("PhaseSpace:pTHatMin = 0.0");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");

  pythia.readString("Ropewalk:RopeHadronization = on");
  pythia.readString("Ropewalk:doShoving = on");
  pythia.readString("Ropewalk:doFlavour = off");
  pythia.readString("Ropewalk:rCutOff = 10.0");
  pythia.readString("Ropewalk:limitMom = on");
  pythia.readString("Ropewalk:pTcut = 2.0");
  pythia.readString("Ropewalk:r0 = 0.41");
  pythia.readString("Ropewalk:m0 = 0.2");
  pythia.readString("Ropewalk:gAmplitude = 10.0");
  pythia.readString("Ropewalk:gExponent = 1.0");
  pythia.readString("Ropewalk:deltat = 0.1");
  pythia.readString("Ropewalk:tShove = 1.");
  pythia.readString("Ropewalk:deltay = 0.1");
  pythia.readString("Ropewalk:tInit = 1.5");
  // Enabling setting of vertex information.
  pythia.readString("PartonVertex:setVertex = on");
  pythia.readString("PartonVertex:protonRadius = 0.7");
  pythia.readString("PartonVertex:emissionWidth = 0.1");

  pythia.init();  

  const int maxEvent = 100000;
  while(genTree_p->GetEntries() < maxEvent){
    if(!pythia.next()) continue;

    pthat_ = pythia.info.pTHat();
    pthatWeight_ = pythia.info.weight();
    nParticle_ = 0;
    particleData pData;
    std::vector<fastjet::PseudoJet> particles;

    for(int i = 0; i < pythia.event.size(); ++i){
      if(!pythia.event[i].isFinal()) continue;
      
      TLorentzVector temp(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
      
      if(TMath::Abs(pythia.event[i].id()) == 12) continue;
      if(TMath::Abs(pythia.event[i].id()) == 14) continue;
      if(TMath::Abs(pythia.event[i].id()) == 16) continue;
      if(TMath::Abs(temp.Eta()) > 10.) continue;
      if(temp.Pt() < .1) continue;

      px_[nParticle_] = pythia.event[i].px();
      py_[nParticle_] = pythia.event[i].py();
      pz_[nParticle_] = pythia.event[i].pz();
      mass_[nParticle_] = pythia.event[i].m();
      theta_[nParticle_] = pythia.event[i].theta();
      pt_[nParticle_] = temp.Pt();
      phi_[nParticle_] = temp.Phi();
      eta_[nParticle_] = temp.Eta();
      pid_[nParticle_] = pythia.event[i].id();
      pwflag_[nParticle_] = pythia.event[i].charge();

      pData.px[nParticle_] = pythia.event[i].px();
      pData.py[nParticle_] = pythia.event[i].py();
      pData.pz[nParticle_] = pythia.event[i].pz();
      
      fastjet::PseudoJet particle(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
      int tempFlag = TMath::Abs(pythia.event[i].charge());
      particle.set_user_index(tempFlag);
      particles.push_back(particle);

      ++nParticle_;
    }

    pData.nParticle = nParticle_;

    if(nParticle_ > 0){
      TVector3 thrust = getThrust(pData.nParticle, pData.px, pData.py, pData.pz, THRUST::OPTIMAL);
      TTheta_ = thrust.Theta();
      TPhi_ = thrust.Phi();

      for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	processJets(particles, jDef[jIter], jDefReclust[jIter], &(jData[jIter]), jtPtCut);
      }

      genTree_p->Fill();

      for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	jetTree_p[jIter]->Fill();
      }
    }
  }

  inFile_p->cd();
  genTree_p->Write("", TObject::kOverwrite);
  delete genTree_p;

  for(int jIter = 0; jIter < nJtAlgo; ++jIter){
    jetTree_p[jIter]->Write("", TObject::kOverwrite);
    delete jetTree_p[jIter];
  }

  inFile_p->Close();
  delete inFile_p;

  return 0;
}
