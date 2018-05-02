#include "DataProcessing/src/mixFile.cc"

int main(int argc, char *argv[])
{
  if(argc != 2 && argc != 3 && argc !=4){
    std::cout << "Usage: ./bin/scan.exe <inFileName> <OPT-outFileName> <OPT-NEvts2Mix>" << std::endl;
    return 1;
  }

  std::cout << "Begin processing..." << std::endl;
  int retVal = 0;
  if(argc == 2) retVal += makeMixFile(argv[1]);
  else if(argc == 3) retVal += makeMixFile(argv[1], argv[2]);
  else if(argc == 4) retVal += makeMixFile(argv[1], argv[2], std::atoi(argv[3]));
  return retVal;
}
