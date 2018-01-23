#ifndef PROCESSJETS_SCAN_H
#define PROCESSJETS_SCAN_H

#include <vector>
#include <string>

#include "TLorentzVector.h"

#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include "include/jetData.h"

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

#endif
