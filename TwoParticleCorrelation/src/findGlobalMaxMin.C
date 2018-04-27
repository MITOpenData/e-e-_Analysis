#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TDatime.h"
#include "TObjArray.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"

#include "../../DataProcessing/include/particleData.h"
#include "../../DataProcessing/include/eventData.h"
#include "../../DataProcessing/include/jetData.h"
#include "../../DataProcessing/include/returnRootFileContentsList.h"
#include "../../DataProcessing/include/removeVectorDuplicates.h"
#include "../../DataProcessing/include/histDefUtility.h"
#include "include/handleMultipleFiles.h"
#include "include/Selection.h"

int findGlobalMaxMin(const std::string inFileName)
{
	///// check if there are multiple files /////
	std::vector<std::string> fileList = handleMultipleFiles(inFileName);
	if(fileList.size() == 0)
	{
		std::cout << "Given input \'" << inFileName << "\' doesn't produce valid input. return 1" << std::endl;
		return 1;
	}

	//std::cout<<"Initializing Mins and Maxs..."<<std::endl;
  /////// LOADING THE INITIAL BRANCHES/RANGES FROM FILE 1 ///////

  	TFile* inFile_p = new TFile(fileList.at(0).c_str(), "READ");
    TTree* inTree_p = (TTree*)inFile_p->Get("t");

    // get list of trees in the file, we assume that each file has the same list of trees because they are made
    // from the same processing. if not the case then the user needs to take only the intersection
    std::vector<std::string> listofTrees = returnRootFileContentsList(inFile_p,"TTree");
    removeVectorDuplicates(&listofTrees);
    // initialize branch names for each tree
    std::vector< std::vector<std::string> > listOfCompBranches;
    // initialize mins and maxes for each tree
    std::vector< std::vector<double> > branchMins;
  	std::vector< std::vector<double> > branchMaxs;

  	for (unsigned int tI = 0; tI < listofTrees.size(); tI++)
  	{
  		// pick a tree
  		inTree_p = (TTree*)inFile_p->Get(listofTrees.at(tI).c_str());

  		// load the branch names, mins/maxes
  		std::vector<std::string> tIBranches;
  		std::vector<double> tIMin;
  		std::vector<double> tIMax;
  		listOfCompBranches.push_back(tIBranches);
  		branchMins.push_back(tIMin);
  		branchMaxs.push_back(tIMax);
  		
  		TObjArray* list1_p = (TObjArray*)inTree_p->GetListOfBranches();
  		for(Int_t i = 0; i < list1_p->GetEntries(); ++i)
  		{
	    	listOfCompBranches.at(tI).push_back(list1_p->At(i)->GetName());
	    	// load in the minimum and maximums from the branches
	    	branchMins.at(tI).push_back(inTree_p->GetMinimum(list1_p->At(i)->GetName()));
	    	branchMaxs.at(tI).push_back(inTree_p->GetMaximum(list1_p->At(i)->GetName()));
	    	//std::cout<<listOfCompBranches.at(tI).at(i)<<", "<<branchMins.at(tI).at(i)<<", "<<branchMaxs.at(tI).at(i)<<std::endl;
  		}
  	}
  	// clean up
  	delete inTree_p;
  	inFile_p->Close();
  	delete inFile_p;

  	////////////////////// LOOP OVER OTHER FILES //////////////////////
	for (unsigned int fI = 1; fI < fileList.size(); ++fI)
	{
		inFile_p = new TFile(fileList.at(fI).c_str(), "READ");
		//std::cout<<fileList.at(fI).c_str()<<std::endl;
		for (unsigned int tI = 0; tI < listofTrees.size(); tI++)
  		{
  			// pick a tree
  			//std::cout<<listofTrees.at(tI).c_str()<<std::endl;
  			inTree_p = (TTree*)inFile_p->Get(listofTrees.at(tI).c_str());
  			TObjArray* list2_p = (TObjArray*)inTree_p->GetListOfBranches();
  			// load the branch names, mins/maxes
  			unsigned int pos = 0;		
			while(pos < listOfCompBranches.at(tI).size())
			{
				std::string tempStr = listOfCompBranches.at(tI).at(pos);
				bool isInArr = false;
				for(Int_t i = 0; i < list2_p->GetEntries(); ++i)
				{
					// find the tempStr branch in the next file
					std::string tempStr2 = list2_p->At(i)->GetName();
					if(tempStr.size() == tempStr2.size() && tempStr.find(tempStr2) != std::string::npos)
					{
						isInArr = true;
						break;
					}
				}
				// if the branch is not in both branches then delete it from the list of things to be plotted
				if(!isInArr)
				{
					listOfCompBranches.at(tI).erase(listOfCompBranches.at(tI).begin()+pos);
					branchMins.at(tI).erase(branchMins.at(tI).begin()+pos);
					branchMaxs.at(tI).erase(branchMaxs.at(tI).begin()+pos);
				}
				else
				{
					// if not then update the min and max for that variable
					Double_t tempMin = inTree_p->GetMinimum(listOfCompBranches.at(tI).at(pos).c_str());
					Double_t tempMax = inTree_p->GetMaximum(listOfCompBranches.at(tI).at(pos).c_str());

					if(branchMins.at(tI).at(pos) > tempMin) branchMins.at(tI).at(pos) = tempMin;
					if(branchMaxs.at(tI).at(pos) < tempMax) branchMaxs.at(tI).at(pos) = tempMax;

					//std::cout<<listOfCompBranches.at(tI).at(pos)<<", "<<branchMins.at(tI).at(pos)<<", "<<branchMaxs.at(tI).at(pos)<<std::endl;
					++pos;
				}
			}
			delete inTree_p;
  		}
  		inFile_p->Close();
		delete inFile_p;
	}
	for (unsigned int tI = 0; tI < listofTrees.size(); tI++)
  	{
  		std::cout<<listofTrees.at(tI)<<std::endl;
  		for (unsigned int lI = 0; lI < listOfCompBranches.at(tI).size(); lI++)
  		{
  			std::cout<<listOfCompBranches.at(tI).at(lI)<<", "<<branchMins.at(tI).at(lI)<<", "<<branchMaxs.at(tI).at(lI)<<std::endl;
  		}
  	}
	return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage ./findGlobalMaxMin.exe <inFileName1>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 2) retVal += findGlobalMaxMin(argv[1]);
  return retVal;
}
