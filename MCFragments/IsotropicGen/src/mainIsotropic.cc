// main01.cc is a part of the PYTHIA event generator.
// Copyright (C) 2017 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk.
// It studies the charged multiplicity distribution at the LHC.

#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TNamed.h"

#include <time.h>
#include <sys/time.h>

//fastjet dependencies
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

//mixing
#include "src/mixFile.cc"

//local dependencies
#include "include/jetData.h"
#include "include/particleData.h"
#include "include/eventData.h"
#include "include/boostedEvtData.h"
#include "include/thrustTools.h"
#include "include/boostTools.h"
#include "include/dataTools.h"
#include "include/doGlobalDebug.h"
#include "include/processJets.h"

inline bool passesTheta(float t){return TMath::Abs(TMath::Cos(t))<0.8;}

int main(int argc, char* argv[])
{
  double startTime = get_wall_time();

  if(argc != 4 && argc != 5 && argc != 6 && argc != 7 && argc != 8){
    std::cout << "Usage: ./mainCustom.exe <outFileName> <isSysPP> <jobNum> <maxEvt> <poissonMean> <doRegularThrustWTAAxis> <doAcceptanceCut>" << std::endl;
    return 1;
  }

  bool isSysPP = bool(std::stoi(argv[2]));
  std::string isSysPPStr = "pp";
  if(!isSysPP) isSysPPStr = "ee";

  int tempMaxEvent = 1000;
  if(argc >= 5) tempMaxEvent = std::stoi(argv[4]);
  const int maxEvent = tempMaxEvent;

  int tempPoissonMean_ = 15;
  if(argc >= 6) tempPoissonMean_ = std::stoi(argv[5]);
  const int poissonMean_ = tempPoissonMean_;

  int tempDoRegularAxis = 0;
  if(argc>=7) tempDoRegularAxis = std::stoi(argv[6]);
  const bool doRegularAxis = (bool)tempDoRegularAxis;

  int tempDoAcceptanceCut = 0;
  if(argc>=8) tempDoAcceptanceCut = std::stoi(argv[7]);
  const bool doAcceptanceCut = (bool)tempDoAcceptanceCut;

  const double jtPtCut = .01;
  const int nJtAlgo = 8;
  const double rParam[nJtAlgo] = {0.4, 0.4, 0.8, 0.8, -1., -1., -1., -1.};
  const int nFinalClust[nJtAlgo] = {-1, -1, -1, -1, 2, 2, 3, 3};

  const double recombScheme[nJtAlgo] = {fastjet::E_scheme, fastjet::WTA_modp_scheme, fastjet::E_scheme, fastjet::WTA_modp_scheme, fastjet::E_scheme, fastjet::WTA_modp_scheme, fastjet::E_scheme, fastjet::WTA_modp_scheme};
  fastjet::JetDefinition jDef[nJtAlgo];
  fastjet::JetDefinition jDefReclust[nJtAlgo];
  std::string jetTreeName[nJtAlgo];

  for(int i = 0; i < nJtAlgo; ++i){
    std::string recombSchemeStr = "EScheme";
    if(recombScheme[i] == fastjet::WTA_modp_scheme) recombSchemeStr = "WTAmodpScheme";

    std::string rOrJetStr = "akR" + std::to_string(int(rParam[i]*10));
    if(rParam[i] < 0) rOrJetStr = "ktN" + std::to_string(nFinalClust[i]);

    jetTreeName[i] = rOrJetStr + recombSchemeStr + "JetTree";
  }

  for(int i = 0; i < nJtAlgo; ++i){
    if(isSysPP){
      if(rParam[i] < 0) continue;

      jDef[i] = fastjet::JetDefinition(fastjet::antikt_algorithm, rParam[i], fastjet::RecombinationScheme(recombScheme[i]));
      jDefReclust[i] = fastjet::JetDefinition(fastjet::antikt_algorithm, 5, fastjet::RecombinationScheme(recombScheme[i]));
    }
    else{
      if(rParam[i] > 0){
	jDef[i] = fastjet::JetDefinition(fastjet::ee_genkt_algorithm, rParam[i], -1, fastjet::RecombinationScheme(recombScheme[i]));
	jDefReclust[i] = fastjet::JetDefinition(fastjet::ee_genkt_algorithm, 5, -1, fastjet::RecombinationScheme(recombScheme[i]))
;
      }
      else{
	jDef[i] = fastjet::JetDefinition(fastjet::ee_kt_algorithm, fastjet::RecombinationScheme(recombScheme[i]));
	jDefReclust[i] = fastjet::JetDefinition(fastjet::ee_kt_algorithm, fastjet::RecombinationScheme(recombScheme[i]));
      }
    }
  }
  
  std::string outFileName = argv[1];
  if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), 5, "");
  outFileName = outFileName + "_" + isSysPPStr + "_JobNum" + argv[3] + "_nEvt" + std::to_string(maxEvent) + "_poissonMean" + std::to_string(poissonMean_) + ".root";
  
  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TTree* genTree_p = new TTree("t", "t");
  TTree* boostedGenTree_p = new TTree("BoostedWTAR8Evt", "BoostedWTAR8Evt");
  TTree* jetTree_p[nJtAlgo];
  for(int i = 0; i < nJtAlgo; ++i){
    if(rParam[i] < 0 && isSysPP) continue;
    jetTree_p[i] = new TTree(jetTreeName[i].c_str(), jetTreeName[i].c_str());
  }

  particleData pData;
  jetData jData[nJtAlgo];
  eventData eData;
  boostedEvtData bData;

  Int_t seed_;
  Float_t pthat_;
  Float_t pthatWeight_;
  Int_t scatterMom_;
  std::vector<int> scatterMomDaughter_;

  genTree_p->Branch("seed", &seed_, "seed/I");
  genTree_p->Branch("pthat", &pthat_, "pthat/F");
  genTree_p->Branch("pthatWeight", &pthatWeight_, "pthatWeight/F");
  genTree_p->Branch("scatterMom", &scatterMom_, "scatterMom/I");
  genTree_p->Branch("scatterMomDaughter", &scatterMomDaughter_);
  pData.SetBranchWrite(genTree_p);
  eData.SetBranchWrite(genTree_p);
  bData.SetBranchWrite(boostedGenTree_p);

  for(int i = 0; i < nJtAlgo; ++i){
    if(rParam[i] < 0 && isSysPP) continue;
    jData[i].SetBranchWrite(jetTree_p[i]);
  }

  const Int_t seedStart = 55217;//Chosen by drawing cards from a deck, will be used to do versioning validation from now on
  const Int_t totalSeed = seedStart + std::stoi(argv[3]);

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  TRandom3 * r = new TRandom3(totalSeed);

  int currTreeFills = -1;
  while(genTree_p->GetEntries() < maxEvent){
   if(currTreeFills%1000 == 0) std::cout << currTreeFills << std::endl;
   //isotropic generator here
   //size of event from a Poissonian of mean 15
   int eventSize = r->Poisson(poissonMean_); 
   std::vector< TLorentzVector > particleList;

   //generate particles
   for(int i = 0; i<eventSize; i++){
     float phi = r->Rndm()*TMath::Pi()*2;
     float theta = TMath::ACos(r->Rndm()*2-1);
     if(doAcceptanceCut && !passesTheta(theta)){
       eventSize--;
       i--;
       continue;
     }
     float p = r->Exp(6.0);
     TVector3 temp = TVector3(0,0,0);
     temp.SetMagThetaPhi(p,theta,phi);
     TLorentzVector tempLV = TLorentzVector(temp,TMath::Sqrt(p*p+0.14*0.14));
     particleList.push_back(tempLV);
   }

    seed_ = totalSeed;
    pthat_ = -1;
    pthatWeight_ = 1;
    pData.nParticle = 0;
    bData.nParticle = 0;
    eData.nChargedHadrons = 0;
    eData.nChargedHadronsHP = eventSize;
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

    for(int i = 0; i < eventSize; ++i){
      TLorentzVector temp = particleList[i];

      netP -= TVector3(temp.X(), temp.Y(), temp.Z());

      pData.px[pData.nParticle] = temp.X();
      pData.py[pData.nParticle] = temp.Y();
      pData.pz[pData.nParticle] = temp.Z();
      pData.mass[pData.nParticle] = 0.14;
      pData.theta[pData.nParticle] = temp.Theta();
      pData.pt[pData.nParticle] = temp.Pt();
      pData.pmag[pData.nParticle] = temp.P();
      pData.phi[pData.nParticle] = temp.Phi();
      pData.eta[pData.nParticle] = temp.Eta();
      pData.rap[pData.nParticle] = temp.Rapidity();
      pData.pid[pData.nParticle] = -1;
      pData.charge[pData.nParticle] = 1;

      //done based on here: https://github.com/ginnocen/StudyMult/blob/master/DataProcessing/src/createMC.c#L13
      pData.pwflag[pData.nParticle] = 0;

      fastjet::PseudoJet particle(temp.X(), temp.Y(), temp.Z(), temp.T());
      int tempFlag = TMath::Abs(1);
      particle.set_user_index(tempFlag);
      particles.push_back(particle);

      if(pData.pwflag[pData.nParticle] == 0){
	netP_charged -= TVector3(temp.X(), temp.Y(), temp.Z());

	pDataCh.px[eData.nChargedHadrons] = temp.X();
        pDataCh.py[eData.nChargedHadrons] = temp.Y();
        pDataCh.pz[eData.nChargedHadrons] = temp.Z();

	++(pDataCh.nParticle);
	++(eData.nChargedHadrons);
	if(temp.Pt() > 0.40 && TMath::Abs(temp.Eta()) < 2.4) ++eData.nChargedHadrons_GT0p4;
	if(temp.Pt() > 0.40 && TMath::Abs(temp.Eta()) < 2.4) ++eData.nChargedHadrons_GT0p4Thrust;
      }

      ++(pData.nParticle);
      ++(bData.nParticle);
    }

    if(pData.nParticle > 0){
      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      TVector3 thrust = TVector3(0,0,0);
      thrust.SetMagThetaPhi(1,TMath::ACos(r->Rndm()*2-1),r->Rndm()*TMath::Pi()*2);
      TVector3 thrustCh = thrust;
      if(tempDoRegularAxis){
        thrust = getThrust(pData.nParticle, pData.px, pData.py, pData.pz, THRUST::OPTIMAL);
        thrustCh = getThrust(pDataCh.nParticle, pDataCh.px, pDataCh.py, pDataCh.pz, THRUST::OPTIMAL);
      }

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
      eData.passesAll = 1;


      for(int pI = 0; pI < pData.nParticle; ++pI){
	pData.pt_wrtThr[pI] = ptFromThrust(thrust, TVector3(pData.px[pI], pData.py[pI], pData.pz[pI]));
	pData.theta_wrtThr[pI] = thetaFromThrust(thrust, TVector3(pData.px[pI], pData.py[pI], pData.pz[pI]));
	pData.eta_wrtThr[pI] = etaFromThrust(thrust, TVector3(pData.px[pI], pData.py[pI], pData.pz[pI]));
	pData.phi_wrtThr[pI] = phiFromThrust(thrust, TVector3(pData.px[pI], pData.py[pI], pData.pz[pI]));
	pData.pt_wrtThrPerp[pI] = ptFromThrust(thrust, TVector3(pData.px[pI], pData.py[pI], pData.pz[pI]), true);
	pData.theta_wrtThrPerp[pI] = thetaFromThrust(thrust, TVector3(pData.px[pI], pData.py[pI], pData.pz[pI]), true);
	pData.eta_wrtThrPerp[pI] = etaFromThrust(thrust, TVector3(pData.px[pI], pData.py[pI], pData.pz[pI]), true);
	pData.phi_wrtThrPerp[pI] = phiFromThrust(thrust, TVector3(pData.px[pI], pData.py[pI], pData.pz[pI]), true);


	pData.pt_wrtChThr[pI] = ptFromThrust(thrustCh, TVector3(pData.px[pI], pData.py[pI], pData.pz[pI]));
	pData.theta_wrtChThr[pI] = thetaFromThrust(thrustCh, TVector3(pData.px[pI], pData.py[pI], pData.pz[pI]));
	pData.eta_wrtChThr[pI] = etaFromThrust(thrustCh, TVector3(pData.px[pI], pData.py[pI], pData.pz[pI]));
	pData.phi_wrtChThr[pI] = phiFromThrust(thrustCh, TVector3(pData.px[pI], pData.py[pI], pData.pz[pI]));
	pData.pt_wrtChThrPerp[pI] = ptFromThrust(thrustCh, TVector3(pData.px[pI], pData.py[pI], pData.pz[pI]),true);
	pData.theta_wrtChThrPerp[pI] = thetaFromThrust(thrustCh, TVector3(pData.px[pI], pData.py[pI], pData.pz[pI]),true);
	pData.eta_wrtChThrPerp[pI] = etaFromThrust(thrustCh, TVector3(pData.px[pI], pData.py[pI], pData.pz[pI]),true);
	pData.phi_wrtChThrPerp[pI] = phiFromThrust(thrustCh, TVector3(pData.px[pI], pData.py[pI], pData.pz[pI]),true);
      }

      for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	if(rParam[jIter] < 0 && isSysPP) continue;
	processJets(particles, jDef[jIter], jDefReclust[jIter], &(jData[jIter]), jtPtCut, rParam[jIter], nFinalClust[jIter]);
      }

      genTree_p->Fill();
      for(int jIter = 0; jIter < nJtAlgo; ++jIter){
	if(rParam[jIter] < 0 && isSysPP) continue;

	jetTree_p[jIter]->Fill();
	if(rParam[jIter]==0.8 && recombScheme[jIter]==fastjet::WTA_modp_scheme){
	  if(jData[jIter].nref<2){setBoostedVariables(false, &pData, &bData);}
	  else{
            TVector3 jet1 = TVector3(0,0,0);
            jet1.SetMagThetaPhi(r->Rndm()*20+15,TMath::ACos(r->Rndm()*2-1),r->Rndm()*TMath::Pi()*2);  
            TLorentzVector jet1LV = TLorentzVector(jet1,jet1.Mag());          
            TVector3 jet2 = TVector3(0,0,0);
            while(true){
              jet2.SetMagThetaPhi(r->Rndm()*20+15,TMath::ACos(r->Rndm()*2-1),r->Rndm()*TMath::Pi()*2);  
              if(jet1.Angle(jet2)>TMath::Pi()*3.0/4.0) break;
            }
            TLorentzVector jet2LV = TLorentzVector(jet2,jet2.Mag());          

	    TVector3 wtaBoost = findBack2BackBoost(jet1LV, jet2LV);
            if(tempDoRegularAxis){
              wtaBoost = findBack2BackBoost(jData[jIter].fourJet[0],jData[jIter].fourJet[1]);
	      setBoostedVariables(true, &pData, &bData, jData[jIter].fourJet[0], wtaBoost);
            } else {
	      setBoostedVariables(true, &pData, &bData, jet1LV, wtaBoost);
            }
	  }
	  boostedGenTree_p->Fill();
	}
      }
    }
    currTreeFills++;
  }

  outFile_p->cd();
  genTree_p->Write("", TObject::kOverwrite);
  delete genTree_p;

  for(int jIter = 0; jIter < nJtAlgo; ++jIter){
    if(rParam[jIter] < 0 && isSysPP) continue;

    jetTree_p[jIter]->Write("", TObject::kOverwrite);
    delete jetTree_p[jIter];
  }

  boostedGenTree_p->Write("", TObject::kOverwrite);
  delete boostedGenTree_p;

  outFile_p->mkdir("infoDir");
  outFile_p->cd("infoDir");

  TNamed jobNum(("jobNum" + std::string(argv[3]) + "_Seed").c_str(), (std::to_string(totalSeed)).c_str());
  jobNum.Write("", TObject::kOverwrite);

  outFile_p->Close();
  delete outFile_p;

  std::cout << "Job Complete" << std::endl;

  double endTime = get_wall_time();

  std::cout << "Start to end: " << startTime << ", " << endTime << std::endl;

  std::cout << "Mixing File..." << std::endl;
  makeMixFile(outFileName, "", 1);
  std::cout << "Job done" << std::endl;

  return 0;
}
