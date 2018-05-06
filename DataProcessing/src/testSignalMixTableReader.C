#include <iostream>

#include "DataProcessing/include/signalMixTableReader.h"

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/testSignalMixTableReader.exe <inFileName>" << std::endl;
    return 1;
  }
  
  signalMixTableReader test(argv[1]);
  std::cout << test.getSigOverMixFactor(19, .4) << std::endl;
  std::cout << test.getSigOverMixFactor(59, 5) << std::endl;

  return 0;
}
