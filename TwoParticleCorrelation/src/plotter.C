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
#include "../include/getLinBins.h"
#include "../include/getLogBins.h"


int plotAllTH1F(const std::string inFileName, const std::string dataName)
{
    xjjroot::setgstyle();
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    // List of plots to draw in log scale
    static const std::string arr[] = {"jtm","jtN","pt","px","py","pz","pmag", "missP","missPt","missChargedP","missChargedPt","mass","ebx","eby"};
    std::vector<std::string> logPlots;
    logPlots.assign(arr, arr + sizeof(arr) / sizeof(arr[0]));


	TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
 
    std::vector<std::string> inTH1FNames = returnRootFileContentsList(inFile_p, "TH1F");
    
    const Int_t nTH1F = (Int_t)inTH1FNames.size();

    TH1F* inTH1F_p[nTH1F];
    TCanvas* inTH1F_c[nTH1F];

    const Int_t nBinsX = 10;
    const Int_t nBinsY = 10;

    for(Int_t iter = 0; iter < nTH1F; iter++){
        inTH1F_p[iter] = (TH1F*)inFile_p->Get(inTH1FNames.at(iter).c_str());
        
        std::string histName = inTH1FNames.at(iter);
        int start_position_to_erase = histName.find("_file");
        histName.erase(start_position_to_erase,histName.length());
        inTH1F_p[iter]->SetName(Form("%s %s", dataName.c_str(),histName.c_str()));

        inTH1F_c[iter] = new TCanvas(histName.c_str(), histName.c_str(), 1000, 1000);
        gPad->SetLeftMargin(gPad->GetLeftMargin()*1.25);
        gPad->SetRightMargin(gPad->GetLeftMargin());
        gPad->SetBottomMargin(gPad->GetLeftMargin());
        gPad->SetTopMargin(.1);

        std::string searchName = histName;
        if(histName.find("pt") != std::string::npos ) searchName = "pt";
        if(histName.find("jtm") != std::string::npos ) searchName = "jtm";
        if(histName.find("jtN") != std::string::npos ) searchName = "jtN";

        TH2F *hempty;

        // Declare binning for hempty histogram
        Double_t binsX[nBinsX+1];
        const Double_t xLow = inTH1F_p[iter]->GetXaxis()->GetBinLowEdge(1);
        const Double_t xHi = inTH1F_p[iter]->GetXaxis()->GetBinLowEdge(inTH1F_p[iter]->GetXaxis()->GetNbins()+1);
        getLinBins(xLow,xHi,nBinsX,binsX);

        if (std::find(logPlots.begin(), logPlots.end(), searchName) != logPlots.end())
        {
            gPad->SetLogy();
            // save the log version
            Double_t binsYLog[nBinsY+1];
            const Double_t logLow = .000000005;
            const Double_t logHi = 5.0*inTH1F_p[iter]->GetMaximum();
            getLogBins(logLow, logHi, nBinsY, binsYLog);
            hempty = new TH2F(histName.c_str(),Form(";%s;Counts",histName.c_str()),nBinsX, binsX, nBinsY, binsYLog);
        }
        else 
        {   
            Double_t binsY[nBinsY+1];
            const Double_t yLow = 0;
            const Double_t yHi = 1.2*inTH1F_p[iter]->GetMaximum();
            getLinBins(yLow,yHi,nBinsY,binsY);
            hempty = new TH2F(histName.c_str(),Form(";%s;Counts",histName.c_str()),nBinsX, binsX, nBinsY, binsY);
        }

    	xjjroot::sethempty(hempty,0,0);
    	hempty->Draw();
    	xjjroot::setthgrstyle(inTH1F_p[iter], kBlack, 21, 1.2, kBlack, 1, 1, -1, -1, -1);
   		inTH1F_p[iter]->Draw("pe sames");
        
        // save the normal plot
        std::string saveNamePDF = "../../AnalysisNote/ALEPH/images/DQC/" + dataName + "/" + histName;
        saveNamePDF = saveNamePDF + ".pdf";
        inTH1F_c[iter]->SaveAs(saveNamePDF.c_str());

        /*
        hempty = new TH2F(histName.c_str(),Form(";%s;Counts",histName.c_str()),nBinsX, binsX, nBinsY, binsYLog);
        xjjroot::sethempty(hempty,0,0);
        hempty->Draw();
        xjjroot::setthgrstyle(inTH1F_p[iter], kBlack, 21, 1.2, kBlack, 1, 1, -1, -1, -1);

        inTH1F_p[iter]->Draw("pe sames");
        gPad->SetLogy();
        gPad->Modified();

        saveNamePDF = "../../AnalysisNote/ALEPH/images/DQC/" + dataName + "/" + histName + "_log";
        saveNamePDF = saveNamePDF + ".pdf";
        inTH1F_c[iter]->SaveAs(saveNamePDF.c_str());
        */
        delete hempty;
        delete inTH1F_c[iter];
    }

    return 0;
}