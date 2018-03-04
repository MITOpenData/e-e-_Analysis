// Created by Anthony Badea 03/04/18


//c and c++ dependencies
#include <sstream>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm>
 
//root dependencies 
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TList.h>

//local headers
#include "../include/xjjrootuti.h"
#include "../include/returnRootFileContentsList.h"
#include "../../DataProcessing/include/removeVectorDuplicates.h"
#include "../include/getLinBins.h"
#include "../include/getLogBins.h"
#include "../include/handleMultipleFiles.h"

#include "TMath.h"

/*
Merges TH1F histograms with different limits by taking the max and min limits and setting them as the global limits.
*/

int histoHadd_TH1F(const std::string inFileName, const std::string outFileName)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    std::vector<std::string> fileList = handleMultipleFiles(inFileName);
    if(fileList.size() == 0)
    {
      std::cout << "Given input \'" << inFileName << "\' doesn't produce valid input. return 1" << std::endl;
      return 1;
    }

    std::vector< TFile* > inFileList;
    inFileList.push_back(new TFile(fileList.at(0).c_str(), "READ"));
    std::vector<std::string> inTH1FNames = returnRootFileContentsList(inFileList.at(0), "TH1F");
    removeVectorDuplicates(&inTH1FNames);
    for (unsigned int fI = 1; fI < fileList.size() ; ++fI)
    {
        inFileList.push_back(new TFile(fileList.at(fI).c_str(), "READ"));
        std::vector<std::string> temp = returnRootFileContentsList(inFileList.at(fI),"TH1F");
        removeVectorDuplicates(&temp);
        for(unsigned int tI = 0;tI<inTH1FNames.size();tI++)
        {
            if(std::find(temp.begin(),temp.end(),inTH1FNames.at(tI).c_str()) == temp.end())
            {
                std::cout<<Form("Histogram missing %s from file 2.",inTH1FNames.at(tI).c_str())<<std::endl;
                inTH1FNames.erase(inTH1FNames.begin() + tI);
            }
        }

    }

    const Int_t nTH1F = (Int_t)inTH1FNames.size();
    const Int_t numFiles = fileList.size();

    TH1F* inTH1F_p[nTH1F][numFiles];
    TH1F * mergedHist[nTH1F];

    TList * list;
    TFile * outFile_p = new TFile(Form("%s.root",outFileName.c_str()),"RECREATE");
    for(Int_t iter = 0; iter < nTH1F; iter++)
    {
        list = new TList;

        // Find the xMin,xMax of all histograms
        inTH1F_p[iter][0] = (TH1F*)inFileList.at(0)->Get(inTH1FNames.at(iter).c_str());
        Double_t xMin = inTH1F_p[iter][0]->GetXaxis()->GetBinLowEdge(1); // since first bin is for underflow
        Double_t xMax = inTH1F_p[iter][0]->GetXaxis()->GetBinLowEdge(inTH1F_p[iter][0]->GetXaxis()->GetNbins()+1); // since last bin is overflow

        for(unsigned int fI = 1; fI < fileList.size(); ++fI)
        {
            inTH1F_p[iter][fI] = (TH1F*)inFileList.at(fI)->Get(inTH1FNames.at(iter).c_str());
            // check if xMin,xMax, or yMax needs to be adjusted 
            if( xMin > inTH1F_p[iter][fI]->GetBinLowEdge(1)) xMin = inTH1F_p[iter][fI]->GetBinLowEdge(1);
            if ( xMax < inTH1F_p[iter][fI]->GetBinLowEdge(inTH1F_p[iter][fI]->GetNbinsX()+1)) xMax = inTH1F_p[iter][fI]->GetBinLowEdge(inTH1F_p[iter][fI]->GetNbinsX()+1);
        }
        // Force all histograms to have the same xMin,xMax
        for(unsigned int fI = 0; fI < fileList.size(); ++fI)
        {
            inTH1F_p[iter][fI]->GetXaxis()->SetLimits(xMin,xMax);
            list->Add(inTH1F_p[iter][fI]);
        }
        mergedHist[iter] = (TH1F*)inTH1F_p[iter][0]->Clone(Form("%s",inTH1FNames.at(iter).c_str()));
        mergedHist[iter]->Reset();
        mergedHist[iter]->Merge(list);
        mergedHist[iter]->Write("", TObject::kOverwrite);

        delete list;
    }

    outFile_p->Close();
    delete outFile_p;

    return 0;
}











