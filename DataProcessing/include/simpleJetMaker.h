//Chris McGinn; 2017.08.17
#ifndef SIMPLEJETMAKER_H
#define SIMPLEJETMAKER_H

//cpp dependencies
#include <vector>
#include <string>

//fastjet dependencies
#include "fastjet/ClusterSequence.hh"

class simpleJetMaker{
 public:
  simpleJetMaker(fastjet::JetDefinition inDef);
  ~simpleJetMaker(){};

  fastjet::JetDefinition jetDef;
  std::vector<fastjet::PseudoJet> getSimpleJets(std::vector<fastjet::PseudoJet> particles);
};

simpleJetMaker::simpleJetMaker(fastjet::JetDefinition inDef)
{
  jetDef = inDef;
  return;
}

std::vector<fastjet::PseudoJet> simpleJetMaker::getSimpleJets(std::vector<fastjet::PseudoJet> particles)
{
  fastjet::ClusterSequence cs(particles, jetDef);
  std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
  return jets;
}

#endif
