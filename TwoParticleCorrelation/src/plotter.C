#include "../include/xjjrootuti.h"
#include "../include/returnRootFileContentsList.h"
#include "../include/getLinBins.h"
#include "../include/getLogBins.h"


int plotAllTH1F(const std::string inFileName)
{
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

        inTH1F_c[iter] = new TCanvas(histName.c_str(), histName.c_str(), 1000, 1000);
        gPad->SetLeftMargin(gPad->GetLeftMargin()*1.25);
        gPad->SetRightMargin(gPad->GetLeftMargin());
        gPad->SetBottomMargin(gPad->GetLeftMargin());
        gPad->SetTopMargin(.01);
        gStyle->SetOptStat(1);

        // Declare binning for hempty histogram
        Double_t binsX[nBinsX+1];
    	const Double_t xLow = inTH1F_p[iter]->GetXaxis()->GetBinLowEdge(1);
    	const Double_t xHi = inTH1F_p[iter]->GetXaxis()->GetBinLowEdge(inTH1F_p[iter]->GetXaxis()->GetNbins()+1);
    	getLinBins(xLow,xHi,nBinsX,binsX);

    	Double_t binsY[nBinsY+1];
    	const Double_t yLow = 0;
    	const Double_t yHi = 1.2*inTH1F_p[iter]->GetMaximum();
    	getLinBins(yLow,yHi,nBinsY,binsY);

        TH2F *hempty = new TH2F(histName.c_str(),Form(";%s;Counts",histName.c_str()),nBinsX, binsX, nBinsY, binsY);
    	xjjroot::sethempty(hempty,0,0.3);
    	hempty->DrawCopy();
    	xjjroot::setthgrstyle(inTH1F_p[iter], kBlack, 21, 1.2, kBlack, 1, 1, -1, -1, -1);
   		inTH1F_p[iter]->DrawCopy("pe same");

        gPad->Modified();
        
        TDatime* date = new TDatime();
        
        std::string saveNamePDF = "../pdfDir/" + histName + std::to_string(date->GetDate());
        saveNamePDF = saveNamePDF + ".pdf";
        
        inTH1F_c[iter]->SaveAs(saveNamePDF.c_str());
        
        delete date;
        delete hempty;
        delete inTH1F_c[iter];
    }

    return 0;
}