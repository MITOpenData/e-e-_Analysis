//
//  handleMultipleFiles.h
//  
//
//
//

#ifndef handleMultipleFiles_h
#define handleMultipleFiles_h

#include <string>
#include <vector>
#include <fstream>

// used to return a vector of strings if the input is a text file containing file name
// or just the single root file if that was what was passed in
std::vector<std::string> handleMultipleFiles(const std::string inFileName)
{
    std::vector<std::string> fileList;
    if(inFileName.find(".root") != std::string::npos) fileList.push_back(inFileName);
    else if(inFileName.find(".txt") != std::string::npos){
        std::ifstream file(inFileName.c_str());
        std::string tempStr;
        
        while(std::getline(file, tempStr)){
            while(tempStr.substr(0,1).find(" ") != std::string::npos) tempStr.replace(0,1,"");
                if(tempStr.size() == 0) continue;
            if(tempStr.find(".root") == std::string::npos) continue;
            
            fileList.push_back(tempStr);
        }
        
        file.close();
    }
    return fileList;
}

/*
 INSERT THIS INTO YOUR CODE TO CHECK IF FILELIST IS EMPTY
 
 if(fileList.size() == 0){
 std::cout << "Given input \'" << inFileName << "\' doesn't produce valid input. return 1" << std::endl;
 return 1;
 }
 */
#endif /* handleMultipleFiles_h */
