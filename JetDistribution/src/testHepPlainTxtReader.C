
//cpp dependencies
#include <iostream>
#include <string>

//ROOT dependencies
#include "TH1D.h"

//Local StudyMult (JetDistribution) dependencies
#include "JetDistribution/include/hepPlainTxtReader.h"

int testHepPlainTxtReader(const std::string inFileName)
{
  hepPlainTxtReader testReader(inFileName, "/HepData/6834/d15-x1-y1");
  testReader.Print();
  testReader.Reset();
  testReader.Print();
  testReader.Init(inFileName, "BadName");
  testReader.Reset();

  testReader.Init(inFileName, "/HepData/6834/d15-x1-y1");
  TH1D* hist_p = NULL;
  testReader.GetHistogram(&hist_p, "histTest");
  hist_p->Print("ALL");
  testReader.GetHistogram(&hist_p, "histTest");

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/testHepPlainTextReader.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += testHepPlainTxtReader(argv[1]);
  return retVal;
}
