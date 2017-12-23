#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TDatime.h"
#include "TObjArray.h"
#include "TH1F.h"

#include "include/particleData.h"

int doComparison(const std::string inFileName1, const std::string inFileName2, std::string outFileName = "")
{
  TDatime* date = new TDatime();
  if(outFileName.size() == 0){
    outFileName = inFileName1;
    while(outFileName.find("/") != std::string::npos) outFileName.replace(0, outFileName.find("/")+1, "");
    while(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), 5, "");

    std::string appOutFileName = inFileName2;
    while(appOutFileName.find("/") != std::string::npos) appOutFileName.replace(0, appOutFileName.find("/")+1, "");
    outFileName = outFileName + "_" + appOutFileName;
  }
  while(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), 5, "");    
  outFileName = outFileName + "_DOCOMP_" + std::to_string(date->GetDate()) + ".root";
  delete date;

  std::vector<std::string> listOfCompBranches;
  std::vector<double> branchMins;
  std::vector<double> branchMaxs;

  particleData pData1;
  particleData pData2;

  std::map<ULong64_t, Int_t> f1RunEvtToEntry;

  TFile* inFile1_p = new TFile(inFileName1.c_str(), "READ");
  TTree* inTree1_p = (TTree*)inFile1_p->Get("t");
  inTree1_p->SetBranchStatus("*", 0);
  inTree1_p->SetBranchStatus("RunNo", 1);
  inTree1_p->SetBranchStatus("EventNo", 1);

  inTree1_p->SetBranchAddress("RunNo", &pData1.RunNo);
  inTree1_p->SetBranchAddress("EventNo", &pData1.EventNo);

  for(Int_t entry = 0; entry < inTree1_p->GetEntries(); ++entry){
    inTree1_p->GetEntry(entry);

    ULong64_t key = pData1.RunNo*10000000 + pData1.EventNo;
    if(f1RunEvtToEntry.find(key) != f1RunEvtToEntry.end()){
      std::cout << "Uhoh key duplication \'" << key << "\'..." << std::endl;
    }
    f1RunEvtToEntry[key] = entry;
  }

  TObjArray* list1_p = (TObjArray*)inTree1_p->GetListOfBranches();
  for(Int_t i = 0; i < list1_p->GetEntries(); ++i){
    listOfCompBranches.push_back(list1_p->At(i)->GetName());
    branchMins.push_back(inTree1_p->GetMinimum(list1_p->At(i)->GetName()));
    branchMaxs.push_back(inTree1_p->GetMaximum(list1_p->At(i)->GetName()));
  }
  inFile1_p->Close();
  delete inFile1_p;

  TFile* inFile2_p = new TFile(inFileName2.c_str(), "READ");
  TTree* inTree2_p = (TTree*)inFile2_p->Get("t");
  TObjArray* list2_p = (TObjArray*)inTree2_p->GetListOfBranches();

  unsigned int pos = 0;
  while(pos < listOfCompBranches.size()){
    std::string tempStr = listOfCompBranches.at(pos);

    bool isInArr = false;
    for(Int_t i = 0; i < list2_p->GetEntries(); ++i){
      std::string tempStr2 = list2_p->At(i)->GetName();

      if(tempStr.size() == tempStr2.size() && tempStr.find(tempStr2) != std::string::npos){
	isInArr = true;
	break;
      }
    }

    if(!isInArr){
      listOfCompBranches.erase(listOfCompBranches.begin()+pos);
      branchMins.erase(branchMins.begin()+pos);
      branchMaxs.erase(branchMaxs.begin()+pos);      
    }
    else{
      Double_t tempMin = inTree2_p->GetMinimum(listOfCompBranches.at(pos).c_str());
      Double_t tempMax = inTree2_p->GetMaximum(listOfCompBranches.at(pos).c_str());

      if(branchMins.at(pos) > tempMin) branchMins.at(pos) = tempMin;
      if(branchMaxs.at(pos) < tempMax) branchMaxs.at(pos) = tempMax;
      ++pos;
    }
  }

  inFile2_p->Close();
  delete inFile2_p;

  const Int_t nVarToComp = listOfCompBranches.size();

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1F* hist1_p[nVarToComp];
  TH1F* hist2_p[nVarToComp];
  TH1F* hist_Rat1Over2_p[nVarToComp];

  TH1F* hist_Delta1From2_EvtByEvt_p[nVarToComp];

  for(Int_t i = 0; i < nVarToComp; ++i){
    double tempInterval = branchMaxs.at(i) - branchMins.at(i);
    branchMaxs.at(i) += tempInterval/10.;
    branchMins.at(i) -= tempInterval/10.;

    if(tempInterval == 0){
      branchMaxs.at(i) += 1;
      branchMins.at(i) -= 1;
    }

    std::cout << branchMins.at(i) << ", " << branchMaxs.at(i) << std::endl;

    hist1_p[i] = new TH1F((listOfCompBranches.at(i) + "_file1_h").c_str(), (";" + listOfCompBranches.at(i) + ";Counts").c_str(), 100, branchMins.at(i), branchMaxs.at(i));
    hist2_p[i] = new TH1F((listOfCompBranches.at(i) + "_file2_h").c_str(), (";" + listOfCompBranches.at(i) + ";Counts").c_str(), 100, branchMins.at(i), branchMaxs.at(i));

    hist_Rat1Over2_p[i] = new TH1F((listOfCompBranches.at(i) + "_file_Rat1Over2_h").c_str(), (";" + listOfCompBranches.at(i) + ";File1/File2").c_str(), 100, branchMins.at(i), branchMaxs.at(i));

    hist_Delta1From2_EvtByEvt_p[i] = new TH1F((listOfCompBranches.at(i) + "_Delta1From2_EvtByEvt_h").c_str(), (";#Delta" + listOfCompBranches.at(i) + "_{1,2};Counts").c_str(), 100, -1, 1);
  }

  inFile1_p = new TFile(inFileName1.c_str(), "READ");
  inTree1_p = (TTree*)inFile1_p->Get("t");
  pData1.SetStatusAndAddressRead(inTree1_p, listOfCompBranches);

  inFile2_p = new TFile(inFileName2.c_str(), "READ");
  inTree2_p = (TTree*)inFile2_p->Get("t");
  pData2.SetStatusAndAddressRead(inTree2_p, listOfCompBranches);

  const Int_t nEntries = inTree2_p->GetEntries();

  std::cout << "Doing full processing..." << std::endl;

  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%1000 == 0) std::cout << " Entry: " << entry << "/" << nEntries << std::endl;
    inTree2_p->GetEntry(entry);

    ULong64_t key = pData2.RunNo*10000000 + pData2.EventNo;
    if(f1RunEvtToEntry.find(key) == f1RunEvtToEntry.end()){
      std::cout << "Uhoh missing key \'" << key << "\'..." << std::endl;
      continue;
    }
    
    inTree1_p->GetEntry(f1RunEvtToEntry[key]);

    for(unsigned int lI = 0; lI < listOfCompBranches.size(); ++lI){
      std::string tempS = listOfCompBranches.at(lI);
      if(tempS.find("nParticle") != std::string::npos && tempS.size() == std::string("nParticle").size()){
	hist1_p[lI]->Fill(pData1.nParticle);
	hist2_p[lI]->Fill(pData2.nParticle);
	hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.nParticle - pData2.nParticle);
      }
      else if(tempS.find("EventNo") != std::string::npos && tempS.size() == std::string("EventNo").size()){
	hist1_p[lI]->Fill(pData1.EventNo);
	hist2_p[lI]->Fill(pData2.EventNo);
	hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.EventNo - pData2.EventNo);
      }
      else if(tempS.find("RunNo") != std::string::npos && tempS.size() == std::string("RunNo").size()){
	hist1_p[lI]->Fill(pData1.RunNo);
	hist2_p[lI]->Fill(pData2.RunNo);
	hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.RunNo - pData2.RunNo);
      }
      else if(tempS.find("year") != std::string::npos && tempS.size() == std::string("year").size()){
	hist1_p[lI]->Fill(pData1.year);
	hist2_p[lI]->Fill(pData2.year);
	hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.year - pData2.year);
      }
      else if(tempS.find("process") != std::string::npos && tempS.size() == std::string("process").size()){
	hist1_p[lI]->Fill(pData1.process);
	hist2_p[lI]->Fill(pData2.process);
	hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.process - pData2.process);
      }
      else if(tempS.find("Energy") != std::string::npos && tempS.size() == std::string("Energy").size()){
	hist1_p[lI]->Fill(pData1.Energy);
	hist2_p[lI]->Fill(pData2.Energy);
	hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.Energy - pData2.Energy);
      }
      else if(tempS.find("bFlag") != std::string::npos && tempS.size() == std::string("bFlag").size()){
	hist1_p[lI]->Fill(pData1.bFlag);
	hist2_p[lI]->Fill(pData2.bFlag);
	hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.bFlag - pData2.bFlag);
      }
      else if(tempS.find("bx") != std::string::npos && tempS.size() == std::string("bx").size()){
	hist1_p[lI]->Fill(pData1.bx);
	hist2_p[lI]->Fill(pData2.bx);
	hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.bx - pData2.bx);
      }
      else if(tempS.find("by") != std::string::npos && tempS.size() == std::string("by").size()){
	hist1_p[lI]->Fill(pData1.by);
	hist2_p[lI]->Fill(pData2.by);
	hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.by - pData2.by);
      }
      else if(tempS.find("ebx") != std::string::npos && tempS.size() == std::string("ebx").size()){
	hist1_p[lI]->Fill(pData1.ebx);
	hist2_p[lI]->Fill(pData2.ebx);
	hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.ebx - pData2.ebx);
      }
      else if(tempS.find("eby") != std::string::npos && tempS.size() == std::string("eby").size()){
	hist1_p[lI]->Fill(pData1.eby);
	hist2_p[lI]->Fill(pData2.eby);
	hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.eby - pData2.eby);
      }
      else if(tempS.find("px") != std::string::npos && tempS.size() == std::string("px").size()){
	for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist1_p[lI]->Fill(pData1.px[pI]);}
	for(Int_t pI = 0; pI < pData2.nParticle; ++pI){hist2_p[lI]->Fill(pData2.px[pI]);}
	if(pData1.nParticle == pData2.nParticle){
	  for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.px[pI] - pData2.px[pI]);}
	}
      }
      else if(tempS.find("py") != std::string::npos && tempS.size() == std::string("py").size()){
	for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist1_p[lI]->Fill(pData1.py[pI]);}
	for(Int_t pI = 0; pI < pData2.nParticle; ++pI){hist2_p[lI]->Fill(pData2.py[pI]);}
	if(pData1.nParticle == pData2.nParticle){
	  for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.py[pI] - pData2.py[pI]);}
	}
      }
      else if(tempS.find("pz") != std::string::npos && tempS.size() == std::string("pz").size()){
	for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist1_p[lI]->Fill(pData1.pz[pI]);}
	for(Int_t pI = 0; pI < pData2.nParticle; ++pI){hist2_p[lI]->Fill(pData2.pz[pI]);}
	if(pData1.nParticle == pData2.nParticle){
	  for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.pz[pI] - pData2.pz[pI]);}
	}
      }
      else if(tempS.find("pt") != std::string::npos && tempS.size() == std::string("pt").size()){
	for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist1_p[lI]->Fill(pData1.pt[pI]);}
	for(Int_t pI = 0; pI < pData2.nParticle; ++pI){hist2_p[lI]->Fill(pData2.pt[pI]);}
	if(pData1.nParticle == pData2.nParticle){
	  for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.pt[pI] - pData2.pt[pI]);}
	}
      }
      else if(tempS.find("pmag") != std::string::npos && tempS.size() == std::string("pmag").size()){
	for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist1_p[lI]->Fill(pData1.pmag[pI]);}
	for(Int_t pI = 0; pI < pData2.nParticle; ++pI){hist2_p[lI]->Fill(pData2.pmag[pI]);}
	if(pData1.nParticle == pData2.nParticle){
	  for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.pmag[pI] - pData2.pmag[pI]);}
	}
      }
      else if(tempS.find("eta") != std::string::npos && tempS.size() == std::string("eta").size()){
	for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist1_p[lI]->Fill(pData1.eta[pI]);}
	for(Int_t pI = 0; pI < pData2.nParticle; ++pI){hist2_p[lI]->Fill(pData2.eta[pI]);}
	if(pData1.nParticle == pData2.nParticle){
	  for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.eta[pI] - pData2.eta[pI]);}
	}
      }
      else if(tempS.find("theta") != std::string::npos && tempS.size() == std::string("theta").size()){
	for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist1_p[lI]->Fill(pData1.theta[pI]);}
	for(Int_t pI = 0; pI < pData2.nParticle; ++pI){hist2_p[lI]->Fill(pData2.theta[pI]);}
	if(pData1.nParticle == pData2.nParticle){
	  for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.theta[pI] - pData2.theta[pI]);}
	}
      }
      else if(tempS.find("phi") != std::string::npos && tempS.size() == std::string("phi").size()){
	for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist1_p[lI]->Fill(pData1.phi[pI]);}
	for(Int_t pI = 0; pI < pData2.nParticle; ++pI){hist2_p[lI]->Fill(pData2.phi[pI]);}
	if(pData1.nParticle == pData2.nParticle){
	  for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.phi[pI] - pData2.phi[pI]);}
	}
      }
      else if(tempS.find("mass") != std::string::npos && tempS.size() == std::string("mass").size()){
	for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist1_p[lI]->Fill(pData1.mass[pI]);}
	for(Int_t pI = 0; pI < pData2.nParticle; ++pI){hist2_p[lI]->Fill(pData2.mass[pI]);}
	if(pData1.nParticle == pData2.nParticle){
	  for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.mass[pI] - pData2.mass[pI]);}
	}
      }
      else if(tempS.find("charge") != std::string::npos && tempS.size() == std::string("charge").size()){
	for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist1_p[lI]->Fill(pData1.charge[pI]);}
	for(Int_t pI = 0; pI < pData2.nParticle; ++pI){hist2_p[lI]->Fill(pData2.charge[pI]);}
	if(pData1.nParticle == pData2.nParticle){
	  for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.charge[pI] - pData2.charge[pI]);}
	}
      }
      else if(tempS.find("pwflag") != std::string::npos && tempS.size() == std::string("pwflag").size()){
	for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist1_p[lI]->Fill(pData1.pwflag[pI]);}
	for(Int_t pI = 0; pI < pData2.nParticle; ++pI){hist2_p[lI]->Fill(pData2.pwflag[pI]);}
	if(pData1.nParticle == pData2.nParticle){
	  for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.pwflag[pI] - pData2.pwflag[pI]);}
	}
      }
      else if(tempS.find("pid") != std::string::npos && tempS.size() == std::string("pid").size()){
	for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist1_p[lI]->Fill(pData1.pid[pI]);}
	for(Int_t pI = 0; pI < pData2.nParticle; ++pI){hist2_p[lI]->Fill(pData2.pid[pI]);}
	if(pData1.nParticle == pData2.nParticle){
	  for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.pid[pI] - pData2.pid[pI]);}
	}
      }
      else if(tempS.find("d0") != std::string::npos && tempS.size() == std::string("d0").size()){
	for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist1_p[lI]->Fill(pData1.d0[pI]);}
	for(Int_t pI = 0; pI < pData2.nParticle; ++pI){hist2_p[lI]->Fill(pData2.d0[pI]);}
	if(pData1.nParticle == pData2.nParticle){
	  for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.d0[pI] - pData2.d0[pI]);}
	}
      }
      else if(tempS.find("z0") != std::string::npos && tempS.size() == std::string("z0").size()){
	for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist1_p[lI]->Fill(pData1.z0[pI]);}
	for(Int_t pI = 0; pI < pData2.nParticle; ++pI){hist2_p[lI]->Fill(pData2.z0[pI]);}
	if(pData1.nParticle == pData2.nParticle){
	  for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.z0[pI] - pData2.z0[pI]);}
	}
      }
      else if(tempS.find("ntpc") != std::string::npos && tempS.size() == std::string("ntpc").size()){
	for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist1_p[lI]->Fill(pData1.ntpc[pI]);}
	for(Int_t pI = 0; pI < pData2.nParticle; ++pI){hist2_p[lI]->Fill(pData2.ntpc[pI]);}
	if(pData1.nParticle == pData2.nParticle){
	  for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.ntpc[pI] - pData2.ntpc[pI]);}
	}
      }
      else if(tempS.find("pt_wrtThr") != std::string::npos && tempS.size() == std::string("pt_wrtThr").size()){
	for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist1_p[lI]->Fill(pData1.pt_wrtThr[pI]);}
	for(Int_t pI = 0; pI < pData2.nParticle; ++pI){hist2_p[lI]->Fill(pData2.pt_wrtThr[pI]);}
	if(pData1.nParticle == pData2.nParticle){
	  for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.pt_wrtThr[pI] - pData2.pt_wrtThr[pI]);}
	}
      }
      else if(tempS.find("eta_wrtThr") != std::string::npos && tempS.size() == std::string("eta_wrtThr").size()){
	for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist1_p[lI]->Fill(pData1.eta_wrtThr[pI]);}
	for(Int_t pI = 0; pI < pData2.nParticle; ++pI){hist2_p[lI]->Fill(pData2.eta_wrtThr[pI]);}
	if(pData1.nParticle == pData2.nParticle){
	  for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.eta_wrtThr[pI] - pData2.eta_wrtThr[pI]);}
	}
      }
      else if(tempS.find("theta_wrtThr") != std::string::npos && tempS.size() == std::string("theta_wrtThr").size()){
	for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist1_p[lI]->Fill(pData1.theta_wrtThr[pI]);}
	for(Int_t pI = 0; pI < pData2.nParticle; ++pI){hist2_p[lI]->Fill(pData2.theta_wrtThr[pI]);}
	if(pData1.nParticle == pData2.nParticle){
	  for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.theta_wrtThr[pI] - pData2.theta_wrtThr[pI]);}
	}
      }
      else if(tempS.find("phi_wrtThr") != std::string::npos && tempS.size() == std::string("phi_wrtThr").size()){
	for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist1_p[lI]->Fill(pData1.phi_wrtThr[pI]);}
	for(Int_t pI = 0; pI < pData2.nParticle; ++pI){hist2_p[lI]->Fill(pData2.phi_wrtThr[pI]);}
	if(pData1.nParticle == pData2.nParticle){
	  for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.phi_wrtThr[pI] - pData2.phi_wrtThr[pI]);}
	}
      }
      else if(tempS.find("pt_wrtChThr") != std::string::npos && tempS.size() == std::string("pt_wrtChThr").size()){
	for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist1_p[lI]->Fill(pData1.pt_wrtChThr[pI]);}
	for(Int_t pI = 0; pI < pData2.nParticle; ++pI){hist2_p[lI]->Fill(pData2.pt_wrtChThr[pI]);}
	if(pData1.nParticle == pData2.nParticle){
	  for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.pt_wrtChThr[pI] - pData2.pt_wrtChThr[pI]);}
	}
      }
      else if(tempS.find("eta_wrtChThr") != std::string::npos && tempS.size() == std::string("eta_wrtChThr").size()){
	for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist1_p[lI]->Fill(pData1.eta_wrtChThr[pI]);}
	for(Int_t pI = 0; pI < pData2.nParticle; ++pI){hist2_p[lI]->Fill(pData2.eta_wrtChThr[pI]);}
	if(pData1.nParticle == pData2.nParticle){
	  for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.eta_wrtChThr[pI] - pData2.eta_wrtChThr[pI]);}
	}
      }
      else if(tempS.find("theta_wrtChThr") != std::string::npos && tempS.size() == std::string("theta_wrtChThr").size()){
	for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist1_p[lI]->Fill(pData1.theta_wrtChThr[pI]);}
	for(Int_t pI = 0; pI < pData2.nParticle; ++pI){hist2_p[lI]->Fill(pData2.theta_wrtChThr[pI]);}
	if(pData1.nParticle == pData2.nParticle){
	  for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.theta_wrtChThr[pI] - pData2.theta_wrtChThr[pI]);}
	}
      }
      else if(tempS.find("phi_wrtChThr") != std::string::npos && tempS.size() == std::string("phi_wrtChThr").size()){
	for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist1_p[lI]->Fill(pData1.phi_wrtChThr[pI]);}
	for(Int_t pI = 0; pI < pData2.nParticle; ++pI){hist2_p[lI]->Fill(pData2.phi_wrtChThr[pI]);}
	if(pData1.nParticle == pData2.nParticle){
	  for(Int_t pI = 0; pI < pData1.nParticle; ++pI){hist_Delta1From2_EvtByEvt_p[lI]->Fill(pData1.phi_wrtChThr[pI] - pData2.phi_wrtChThr[pI]);}
	}
      }
    }
  }

  inFile2_p->Close();
  delete inFile2_p;

  inFile1_p->Close();
  delete inFile1_p;

  outFile_p->cd();

  for(Int_t i = 0; i < nVarToComp; ++i){
    hist1_p[i]->Write("", TObject::kOverwrite);
    hist2_p[i]->Write("", TObject::kOverwrite);

    hist1_p[i]->Sumw2();
    hist2_p[i]->Sumw2();
    hist_Rat1Over2_p[i]->Sumw2();

    hist_Rat1Over2_p[i]->Divide(hist1_p[i], hist2_p[i]);
    hist_Rat1Over2_p[i]->Write("", TObject::kOverwrite);

    hist_Delta1From2_EvtByEvt_p[i]->Write("", TObject::kOverwrite);

    delete hist1_p[i];
    delete hist2_p[i];
    delete hist_Rat1Over2_p[i];

    delete hist_Delta1From2_EvtByEvt_p[i];
  }

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3 && argc != 4){
    std::cout << "Usage ./doComparison.exe <inFileName1> <inFileName2> <outFileName-optional>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 3) retVal += doComparison(argv[1], argv[2]);
  else if(argc == 4) retVal += doComparison(argv[1], argv[2], argv[3]);
  return retVal;
}
