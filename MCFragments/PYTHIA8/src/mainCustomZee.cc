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
#include "include/doGlobalDebug.h"

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
  if(argc != 2 && argc != 3 && argc != 4 && argc != 5){
    std::cout << "Usage: ./mainCustom.exe <outFileName> <maxEvt> <nMinPartChgCut> <doRopeWalk>" << std::endl;
    return 1;
  }

  int tempPartChg_ = 0;
  if(argc >= 4) tempPartChg_ = std::stoi(argv[3]);
  const int nMinPartChgCut_ = tempPartChg_;

  int tempMaxEvent = 1000;
  if(argc >= 3) tempMaxEvent = std::stoi(argv[2]);
  const int maxEvent = tempMaxEvent;
  
  bool tempRopeWalk = false;
  if(argc >= 5) tempRopeWalk = std::stoi(argv[4]);
  const bool doRopeWalk = tempRopeWalk;

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
  
  std::string outFileName = argv[1];
  if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), 5, "");
  outFileName = outFileName + "_nEvt" + std::to_string(maxEvent) + "_nMinChgPart" + std::to_string(nMinPartChgCut_) + "_RopeWalk" + std::to_string(doRopeWalk) + ".root";
  
  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TTree* genTree_p = new TTree("t", "t");
  TTree* jetTree_p[nJtAlgo];
  for(int i = 0; i < nJtAlgo; ++i){jetTree_p[i] = new TTree(jetTreeName[i].c_str(), jetTreeName[i].c_str());}

  particleData pData;
  jetData jData[nJtAlgo];

  Float_t pthat_;
  Float_t pthatWeight_;
  Int_t scatterMom_;
  std::vector<int> scatterMomDaughter_;

  Int_t nParticleChg_;
  Int_t nParticleChg_Pt0p4_Eta1p8_;

  Float_t Thrust_;
  Float_t TTheta_;
  Float_t TPhi_;

  Float_t Thrust_charged_;
  Float_t TTheta_charged_;
  Float_t TPhi_charged_;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  genTree_p->Branch("pthat", &pthat_, "pthat/F");
  genTree_p->Branch("pthatWeight", &pthatWeight_, "pthatWeight/F");
  genTree_p->Branch("scatterMom", &scatterMom_, "scatterMom/I");
  genTree_p->Branch("scatterMomDaughter", &scatterMomDaughter_);
  genTree_p->Branch("nParticleChg", &nParticleChg_, "nParticleChg/I");
  genTree_p->Branch("nParticleChg_Pt0p4_Eta1p8", &nParticleChg_Pt0p4_Eta1p8_, "nParticleChg_Pt0p4_Eta1p8/I");
  pData.SetBranchWrite(genTree_p);
  genTree_p->Branch("Thrust", &Thrust_, "Thrust/F");
  genTree_p->Branch("TTheta", &TTheta_, "TTheta/F");
  genTree_p->Branch("TPhi", &TPhi_, "TPhi/F");
  genTree_p->Branch("Thrust_charged", &Thrust_charged_, "Thrust_charged/F");
  genTree_p->Branch("TTheta_charged", &TTheta_charged_, "TTheta_charged/F");
  genTree_p->Branch("TPhi_charged", &TPhi_charged_, "TPhi_charged/F");

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

  if(doRopeWalk){
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

    pythia.readString("PartonVertex:setVertex = on");
    pythia.readString("PartonVertex:protonRadius = 0.7");
    pythia.readString("PartonVertex:emissionWidth = 0.1");
  }

  pythia.init();

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  int currTreeFills = -1;
  while(genTree_p->GetEntries() < maxEvent){
    if(!pythia.next()) continue;
    if(pythia.event.size() < nMinPartChgCut_) continue;

    if(genTree_p->GetEntries()%(maxEvent/100) == 0 && genTree_p->GetEntries() != currTreeFills){
      std::cout << "GenTree fills: " << genTree_p->GetEntries() << std::endl;
      currTreeFills = genTree_p->GetEntries();
    }

    pthat_ = pythia.info.pTHat();
    pthatWeight_ = pythia.info.weight();
    pData.nParticle = 0;
    nParticleChg_ = 0;
    nParticleChg_Pt0p4_Eta1p8_ = 0;

    particleData pDataCh;
    std::vector<fastjet::PseudoJet> particles;

    scatterMom_ = -999;
    Int_t scatterMomPos = -1;
    scatterMomDaughter_.clear();

    for(int i = 0; i < pythia.event.size(); ++i){
      if(pythia.event[i].mother1() == 3 && pythia.event[i].mother2() == 4){
	scatterMom_ = pythia.event[i].id();
	scatterMomPos = i;
      }
      if(pythia.event[i].mother1() == scatterMomPos) scatterMomDaughter_.push_back(pythia.event[i].id());

      if(!pythia.event[i].isFinal()) continue;
      
      TLorentzVector temp(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
      
      if(TMath::Abs(pythia.event[i].id()) == 12) continue;
      if(TMath::Abs(pythia.event[i].id()) == 14) continue;
      if(TMath::Abs(pythia.event[i].id()) == 16) continue;
      if(TMath::Abs(temp.Eta()) > 10.) continue;
      if(temp.Pt() < .1) continue;

      pData.px[pData.nParticle] = pythia.event[i].px();
      pData.py[pData.nParticle] = pythia.event[i].py();
      pData.pz[pData.nParticle] = pythia.event[i].pz();
      pData.mass[pData.nParticle] = pythia.event[i].m();
      pData.theta[pData.nParticle] = pythia.event[i].theta();
      pData.pt[pData.nParticle] = temp.Pt();
      pData.pmag[pData.nParticle] = temp.P();
      pData.phi[pData.nParticle] = temp.Phi();
      pData.eta[pData.nParticle] = temp.Eta();
      pData.rap[pData.nParticle] = temp.Rapidity();
      pData.pid[pData.nParticle] = pythia.event[i].id();
      pData.charge[pData.nParticle] = pythia.event[i].charge();

      //done based on here: https://github.com/ginnocen/StudyMult/blob/master/DataProcessing/src/createMC.c#L13
      pData.pwflag[pData.nParticle] = -999;
      if(TMath::Abs(pythia.event[i].id()) == 11) pData.pwflag[pData.nParticle] = 1;
      else if(TMath::Abs(pythia.event[i].id()) == 13) pData.pwflag[pData.nParticle] = 2;
      else if(TMath::Abs(pythia.event[pythia.event[i].mother1()].id()) == 310 && pythia.event[i].mother2() == 0) pData.pwflag[pData.nParticle] = 3;
      else if(TMath::Abs(pythia.event[pythia.event[i].mother1()].id()) == 3122 && pythia.event[i].mother2() == 0) pData.pwflag[pData.nParticle] = 3;
      else if(TMath::Abs(pythia.event[i].id()) == 22) pData.pwflag[pData.nParticle] = 4;
      else if(pythia.event[i].charge() == 0) pData.pwflag[pData.nParticle] = 5;
      else if(pythia.event[i].charge() != 0) pData.pwflag[pData.nParticle] = 0;

      fastjet::PseudoJet particle(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
      int tempFlag = TMath::Abs(pythia.event[i].charge());
      particle.set_user_index(tempFlag);
      particles.push_back(particle);

      if(pData.pwflag[pData.nParticle] == 0){
	pDataCh.px[nParticleChg_] = pythia.event[i].px();
        pDataCh.py[nParticleChg_] = pythia.event[i].py();
        pDataCh.pz[nParticleChg_] = pythia.event[i].pz();

	++nParticleChg_;
	if(temp.Pt() > 0.40 && TMath::Abs(temp.Eta()) < 1.8) ++nParticleChg_Pt0p4_Eta1p8_;
      }

      ++pData.nParticle;
    }

    if(nParticleChg_Pt0p4_Eta1p8_ < nMinPartChgCut_) continue;

    pDataCh.nParticle = nParticleChg_;

    if(pData.nParticle > 0){
      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      TVector3 thrust = getThrust(pData.nParticle, pData.px, pData.py, pData.pz, THRUST::OPTIMAL);
      TVector3 thrustCh = getThrust(pDataCh.nParticle, pDataCh.px, pDataCh.py, pDataCh.pz, THRUST::OPTIMAL);

      Thrust_ = thrust.Mag();
      TTheta_ = thrust.Theta();
      TPhi_ = thrust.Phi();

      Thrust_charged_ = thrustCh.Mag();
      TTheta_charged_ = thrustCh.Theta();
      TPhi_charged_ = thrustCh.Phi();
      
      for(int pI = 0; pI < pData.nParticle; ++pI){
	pData.pt_wrtThr[pI] = ptFromThrust(thrust, TVector3(pData.px[pI], pData.py[pI], pData.pz[pI]));
	pData.theta_wrtThr[pI] = thetaFromThrust(thrust, TVector3(pData.px[pI], pData.py[pI], pData.pz[pI]));
	pData.eta_wrtThr[pI] = etaFromThrust(thrust, TVector3(pData.px[pI], pData.py[pI], pData.pz[pI]));
	pData.phi_wrtThr[pI] = phiFromThrust(thrust, TVector3(pData.px[pI], pData.py[pI], pData.pz[pI]));

	pData.pt_wrtChThr[pI] = ptFromThrust(thrustCh, TVector3(pData.px[pI], pData.py[pI], pData.pz[pI]));
	pData.theta_wrtChThr[pI] = thetaFromThrust(thrustCh, TVector3(pData.px[pI], pData.py[pI], pData.pz[pI]));
	pData.eta_wrtChThr[pI] = etaFromThrust(thrustCh, TVector3(pData.px[pI], pData.py[pI], pData.pz[pI]));
	pData.phi_wrtChThr[pI] = phiFromThrust(thrustCh, TVector3(pData.px[pI], pData.py[pI], pData.pz[pI]));
      }

      for(int jIter = 0; jIter < nJtAlgo; ++jIter){processJets(particles, jDef[jIter], jDefReclust[jIter], &(jData[jIter]), jtPtCut);}

      genTree_p->Fill();
      for(int jIter = 0; jIter < nJtAlgo; ++jIter){jetTree_p[jIter]->Fill();}
    }
  }

  outFile_p->cd();
  genTree_p->Write("", TObject::kOverwrite);
  delete genTree_p;

  for(int jIter = 0; jIter < nJtAlgo; ++jIter){
    jetTree_p[jIter]->Write("", TObject::kOverwrite);
    delete jetTree_p[jIter];
  }

  outFile_p->Close();
  delete outFile_p;

  std::cout << "Job Complete" << std::endl;

  return 0;
}
