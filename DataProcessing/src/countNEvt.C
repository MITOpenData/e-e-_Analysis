#include <iostream>
#include <string>
#include <fstream>

int countNEvt(const std::string inFileName)
{
  std::ifstream file(inFileName.c_str());
  std::string line;

  int counter = 0;
  while(std::getline(file,line)){
    if(line.size() == 0) continue;
    while(line.substr(0,1).find(" ") != std::string::npos){line.replace(0,1,"");}
    if(line.size() == 0) continue;

    if(line.substr(0,5).find("-999.") != std::string::npos) ++counter;
  }

  file.close();

  std::cout << counter << std::endl;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage ./countNEvt.exe <inFileName>" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += countNEvt(argv[1]);
  return retVal;
}
