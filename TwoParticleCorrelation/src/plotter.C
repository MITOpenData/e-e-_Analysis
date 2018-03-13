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

//local headers
#include "../include/xjjrootuti.h"
#include "../include/returnRootFileContentsList.h"
#include "../../DataProcessing/include/removeVectorDuplicates.h"
#include "../include/getLinBins.h"
#include "../include/getLogBins.h"
#include "../include/handleMultipleFiles.h"


std::string getDataName(std::string fileName)
{
    std::string dataName = "";
    if(fileName.find("LEP1") != std::string::npos && fileName.find("MC") == std::string::npos) return "LEP1";
    if(fileName.find("LEP1MC") != std::string::npos) return "LEP1 MC";
    if(fileName.find("LEP2") != std::string::npos && fileName.find("MC") == std::string::npos && fileName.find("ECUT") == std::string::npos) return "LEP2";
    if(fileName.find("LEP2MC") != std::string::npos) return "LEP2 MC";
    if(fileName.find("LEP2") != std::string::npos && fileName.find("ECUT") != std::string::npos) return "LEP2 (<100GeV)";
    if(fileName.find("PYTHIA8") != std::string::npos) return "PYTHIA8";
    if(fileName.find("pythia8") != std::string::npos) return "PYTHIA8";

    std::cout<<"No valid data name"<<std::endl;
    return "";
}

int plotAllTH1F(const std::string inFileName, const std::string dataName)
{
    //xjjroot::setgstyle();
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    // List of plots to draw in log scale
    static const std::string arr[] = {"jtm","jtN","pt","px","py","pz","pmag", "missP","missPt","missChargedP","missChargedPt","mass","ebx","eby","nTrk","aplanarity"};
    std::vector<std::string> logPlots;
    logPlots.assign(arr, arr + sizeof(arr) / sizeof(arr[0]));
    // colors for plotting LEP1,LEP2,PYTHIA8
    Int_t colors[3] = {1,632,600}; // kBlack,kRed,kBlue
    Int_t markerStyles[3] = {20,21,22}; // circle, square, triangle

      ///// check if there are multiple files /////
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
    
    /// anthony you are here and are checking that all files have the same list of histograms to plot
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
                tI--; // since the list shortens 
            }

        }

    }
    const Int_t nTH1F = (Int_t)inTH1FNames.size();
    static const int numFiles = fileList.size();

    TH1F* inTH1F_p[nTH1F][numFiles];
    TCanvas* inTH1F_c[nTH1F];

    const Int_t nBinsX = 50; // should be a divisor of 100 = number of bins of inputted histograms
    const Int_t nBinsY = 20;

    for(Int_t iter = 0; iter < nTH1F; iter++)
    {
        std::string histName = inTH1FNames.at(iter);
        int start_position_to_erase = histName.find("_file");
        histName.erase(start_position_to_erase,histName.length());
        if(histName.find("nParticle") != std::string::npos) histName = "nTrk";
        //std::cout<<histName<<std::endl;
        // the default # bins of the input files is 100. Reduce by factor of 5 to 20 bins 
        inTH1F_p[iter][0] = (TH1F*)inFileList.at(0)->Get(inTH1FNames.at(iter).c_str());
        // rebin to the nBinsX
        inTH1F_p[iter][0]->Rebin((inTH1F_p[iter][0]->GetSize()-2)/nBinsX);
        inTH1F_p[iter][0]->SetName(Form("%s",histName.c_str()));
        inTH1F_p[iter][0]->Scale(1.0/inTH1F_p[iter][0]->Integral()); // scale to probability distribution

        // load the histogram min's and max's
        Double_t xMin = inTH1F_p[iter][0]->GetXaxis()->GetBinLowEdge(1); // since first bin is for underflow
        Double_t xMax = inTH1F_p[iter][0]->GetXaxis()->GetBinLowEdge(inTH1F_p[iter][0]->GetXaxis()->GetNbins()+1); // since last bin is overflow
        Double_t yMax = inTH1F_p[iter][0]->GetMaximum();

        if (histName.find("eta") != std::string::npos && xMin < -10.0) xMin = -10;
        if (histName.find("eta") != std::string::npos && xMax > 10.0) xMax = 10;

        for(unsigned int fI = 1; fI < fileList.size(); ++fI)
        {
            std::cout<<inTH1FNames.at(iter).c_str()<<std::endl;
            inTH1F_p[iter][fI] = (TH1F*)inFileList.at(fI)->Get(inTH1FNames.at(iter).c_str());
            inTH1F_p[iter][fI]->Rebin((inTH1F_p[iter][fI]->GetSize()-2)/nBinsX);
            inTH1F_p[iter][fI]->SetName(Form("%s",histName.c_str()));
            inTH1F_p[iter][fI]->Scale(1.0/inTH1F_p[iter][fI]->Integral()); // scale to probability distribution

            // check if xMin,xMax, or yMax needs to be adjusted 
            if( xMin > inTH1F_p[iter][fI]->GetBinLowEdge(1) && inTH1F_p[iter][fI]->GetBinLowEdge(1) > -999.0 /* check for -999 */) 
            {
                if (histName.find("eta") != std::string::npos && inTH1F_p[iter][fI]->GetBinLowEdge(1) < -10.0) xMin = xMin;
                else xMin = inTH1F_p[iter][fI]->GetBinLowEdge(1);
            }
            if ( xMax < inTH1F_p[iter][fI]->GetBinLowEdge(inTH1F_p[iter][fI]->GetNbinsX()+1)) 
            {
                if (histName.find("eta") != std::string::npos && inTH1F_p[iter][fI]->GetBinLowEdge(inTH1F_p[iter][fI]->GetNbinsX()+1) > 10.0) xMax = xMax;
                else xMax = inTH1F_p[iter][fI]->GetBinLowEdge(inTH1F_p[iter][fI]->GetNbinsX()+1);
            }
            if ( yMax < inTH1F_p[iter][fI]->GetMaximum()) yMax = inTH1F_p[iter][fI]->GetMaximum();
        }

        // std::cout<<"xMin = "<<xMin<<" xMax = "<<xMax<<" yMax = "<<yMax<<std::endl;
        inTH1F_c[iter] = new TCanvas(histName.c_str(), histName.c_str(), 1000, 1000);
        gPad->SetLeftMargin(gPad->GetLeftMargin()*1.5);
        gPad->SetRightMargin(gPad->GetLeftMargin());
        gPad->SetBottomMargin(gPad->GetLeftMargin());
        gPad->SetTopMargin(gPad->GetLeftMargin());
        std::string searchName = histName;
        if(histName.find("pt") != std::string::npos ) searchName = "pt";
        if(histName.find("jtm") != std::string::npos ) searchName = "jtm";
        if(histName.find("jtN") != std::string::npos ) searchName = "jtN";

        TH2F *hempty;
        // Declare binning for hempty histogram
        if(histName.find("phi") != std::string::npos) {xMin = -TMath::Pi(); xMax = TMath::Pi();}
        Double_t binsX[nBinsX+1];
        if(searchName == "pt") {xMin = 0; xMax = 100;}
        if(histName.find("eta") != std::string::npos) std::cout<<"ETA"<<yMax<<std::endl;
        const Double_t xLow = xMin;
        const Double_t xHi = xMax;
        getLinBins(xLow,xHi,nBinsX,binsX);
        // rebin all of the histograms to the same binning as hempty
        if (std::find(logPlots.begin(), logPlots.end(), searchName) != logPlots.end())
        {
            gPad->SetLogy();
            // save the log version
            Double_t binsYLog[nBinsY+1];
            const Double_t logLow = .000000005;
            const Double_t logHi = 5.0*yMax;
            getLogBins(logLow, logHi, nBinsY, binsYLog);
            hempty = new TH2F("hempty",Form(";%s;Probability",histName.c_str()),nBinsX, binsX, nBinsY, binsYLog);
        }
        else 
        {   
            Double_t binsY[nBinsY+1];
            const Double_t yLow = 0;
            const Double_t yHi = 1.2*yMax;
            getLinBins(yLow,yHi,nBinsY,binsY);
            hempty = new TH2F("hempty",Form(";%s;Probability",histName.c_str()),nBinsX, binsX, nBinsY, binsY);
        }
    	xjjroot::sethempty(hempty,0,0);
    	hempty->Draw();
        TLegend *l;
        if(histName.find("phi") != std::string::npos || histName.find("Phi") != std::string::npos) l = new TLegend(0.56,0.7-0.15,0.96,0.83-0.15);
        else l = new TLegend(0.56,0.7,0.88,0.83);
        xjjroot::setleg(l);

        // Manually rebin the old histograms to the new binning with new min and max because the default root rebinning is sketchy and sucks
        TH1F* inTH1F_rebin[numFiles];
        for (unsigned int fI = 0; fI < fileList.size(); fI++)
        {
            inTH1F_rebin[fI] = new TH1F(Form("%s_rebin",histName.c_str()),"",nBinsX,binsX);
            std::cout<<inTH1F_rebin[fI]->GetSize()<<std::endl;
            //inTH1F_rebin[fI] = (TH1F*)inTH1F_p[iter][fI]->Rebin(nBinsX,Form("%s_rebin",histName.c_str()),binsX);
            //inTH1F_p[iter][fI] = dynamic_cast<TH1F*>(inTH1F_p[iter][fI]->Rebin(nBinsX,Form("%s_rebin",histName.c_str()),binsX));
            for (unsigned int nI = 1; nI<=inTH1F_p[fI]->GetNbinsX(); nI++)
            {
                Float_t binCenter = inTH1F_rebin[fI]->GetBinCenter(nI);
                Int_t y = (Int_t) inTH1F_p[iter][fI]->GetBinContent(inTH1F_p[iter][fI]->FindBin(binCenter));
                inTH1F_rebin[fI]->SetBinContent(nI,)
                inTH1F_rebin[fI]->Fill(inTH1F_p[iter][fI]->GetBinCenter(nI),inTH1F_p[iter][fI]->GetBinContent(nI));
            }
            xjjroot::setthgrstyle(inTH1F_rebin[fI], colors[fI], 21, 1.2, colors[fI], 1, 1, -1, -1, -1);
            inTH1F_rebin[fI]->SetMarkerStyle(markerStyles[fI]);
            inTH1F_rebin[fI]->Draw("p same");
            l->AddEntry(inTH1F_rebin[fI],getDataName(fileList.at(fI)).c_str(),"p");
        }
        l->Draw("SAME");
        // save the normal plot
        std::string saveNamePDF = "../../AnalysisNote/ALEPH/images/DQC/" + dataName + "/" + histName;
        saveNamePDF = saveNamePDF + ".pdf";
        inTH1F_c[iter]->SaveAs(saveNamePDF.c_str());
        delete l;
        delete hempty;
        delete inTH1F_c[iter];
        for(unsigned int fI = 1; fI < fileList.size(); ++fI) delete inTH1F_p[iter][fI];
    }

    return 0;
}