#include "DataProcessing/src/mixFile.cc"

int main(int argc, char *argv[])
{
  if(argc != 3 && argc != 4 && argc !=5){
    std::cout << "Usage: ./bin/scan.exe <inFileName> <isMC> <OPT-outFileName> <OPT-NEvts2Mix>" << std::endl;
    return 1;
  }

  std::cout << "Begin processing..." << std::endl;
  int retVal = 0;
  if(argc == 3) retVal += makeMixFile(argv[1], std::stoi(argv[2]));
  else if(argc == 4) retVal += makeMixFile(argv[1], std::stoi(argv[2]), argv[3]);
  else if(argc == 5) retVal += makeMixFile(argv[1], std::stoi(argv[2]), argv[3], std::atoi(argv[4]));
  return retVal;
}
