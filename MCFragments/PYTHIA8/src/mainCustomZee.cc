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
#include "include/eventData.h"
#include "include/boostedEvtData.h"
#include "include/thrustTools.h"
#include "include/boostTools.h"
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
      d->jtm[d->nref] = jets.at(i).m();
      d->jteta[d->nref] = jets.at(i).eta();
      d->fourJet[d->nref] = TLorentzVector(jets.at(i).px(), jets.at(i).py(), jets.at(i).pz(), jets.at(i).E());

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
  TTree* boostedGenTree_p = new TTree("BoostedWTAR8Evt", "BoostedWTAR8Evt");
  TTree* jetTree_p[nJtAlgo];
  for(int i = 0; i < nJtAlgo; ++i){jetTree_p[i] = new TTree(jetTreeName[i].c_str(), jetTreeName[i].c_str());}

  particleData pData;
  jetData jData[nJtAlgo];
  eventData eData;
  boostedEvtData bData;

  Float_t pthat_;
  Float_t pthatWeight_;
  Int_t scatterMom_;
  std::vector<int> scatterMomDaughter_;

  genTree_p->Branch("pthat", &pthat_, "pthat/F");
  genTree_p->Branch("pthatWeight", &pthatWeight_, "pthatWeight/F");
  genTree_p->Branch("scatterMom", &scatterMom_, "scatterMom/I");
  genTree_p->Branch("scatterMomDaughter", &scatterMomDaughter_);
  pData.SetBranchWrite(genTree_p);
  eData.SetBranchWrite(genTree_p);
  bData.SetBranchWrite(boostedGenTree_p);

  for(int i = 0; i < nJtAlgo; ++i){jData[i].SetBranchWrite(jetTree_p[i]);}

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
    eData.nChargedHadrons = 0;
    eData.nChargedHadrons_GT0p4 = 0;
    eData.nChargedHadrons_GT0p4Thrust = 0;

    TVector3 netP(0, 0, 0);
    TVector3 netP_charged(0, 0, 0);
    particleData pDataCh;
    pDataCh.nParticle = 0;
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

      netP -= TVector3(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz());

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
	netP_charged -= TVector3(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz());

	pDataCh.px[eData.nChargedHadrons] = pythia.event[i].px();
        pDataCh.py[eData.nChargedHadrons] = pythia.event[i].py();
        pDataCh.pz[eData.nChargedHadrons] = pythia.event[i].pz();

	++(pDataCh.nParticle);
	++(eData.nChargedHadrons);
	if(temp.Pt() > 0.40) ++eData.nChargedHadrons_GT0p4;
	if(temp.Pt() > 0.40) ++eData.nChargedHadrons_GT0p4Thrust;
      }

      ++(pData.nParticle);
    }

    if(eData.nChargedHadrons_GT0p4 < nMinPartChgCut_) continue;

    if(pData.nParticle > 0){
      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      TVector3 thrust = getThrust(pData.nParticle, pData.px, pData.py, pData.pz, THRUST::OPTIMAL);
      TVector3 thrustCh = getThrust(pDataCh.nParticle, pDataCh.px, pDataCh.py, pDataCh.pz, THRUST::OPTIMAL);

      eData.Thrust = thrust.Mag();
      eData.TTheta = thrust.Theta();
      eData.TPhi = thrust.Phi();
      eData.Thrust_charged = thrustCh.Mag();
      eData.TTheta_charged = thrustCh.Theta();
      eData.TPhi_charged = thrustCh.Phi();
      eData.missP = netP.Mag();
      eData.missPt = netP.Perp();
      eData.missTheta = netP.Theta();
      eData.missPhi = netP.Phi();
      eData.missChargedP = netP_charged.Mag();
      eData.missChargedPt = netP_charged.Perp();
      eData.missChargedTheta = netP_charged.Theta();
      eData.missChargedPhi = netP_charged.Phi();


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
      for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	jetTree_p[jIter]->Fill();
	if(rParam[jIter]==0.8 && recombScheme[jIter]==fastjet::WTA_modp_scheme){
	  if(jData[jIter].nref<2){setBoostedVariables(false, &pData, &bData);}
	  else{
	    TVector3 wtaBoost = findBack2BackBoost(jData[jIter].fourJet[0],jData[jIter].fourJet[1]);
	    setBoostedVariables(true, &pData, &bData, jData[jIter].fourJet[0], wtaBoost);
	  }
	  boostedGenTree_p->Fill();
	}
      }
    }
  }

  outFile_p->cd();
  genTree_p->Write("", TObject::kOverwrite);
  delete genTree_p;

  for(int jIter = 0; jIter < nJtAlgo; ++jIter){
    jetTree_p[jIter]->Write("", TObject::kOverwrite);
    delete jetTree_p[jIter];
  }

  boostedGenTree_p->Write("", TObject::kOverwrite);
  delete boostedGenTree_p;

  outFile_p->Close();
  delete outFile_p;

  std::cout << "Job Complete" << std::endl;

  return 0;
}
