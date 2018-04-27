//  Created by Anthony Badea on 4/18/18

#include <TCanvas.h>
#include <TVector3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>
#include <TMath.h>

#include "include/Selection.h"
#include "include/TPCNtupleData.h"
#include "include/smartJetName.h"
#include "../../DataProcessing/include/particleData.h"
#include "../../DataProcessing/include/eventData.h"
#include "../../DataProcessing/include/jetData.h"

// first loop over the tracks and calculate total energy 
// then find the jet energy 
int jetEnergyRatio( const std::string inFileName, // Input file
                    std::string outFileName,        // Output file
		    int jttree = 0,
		    float rCut = 0, // include radius greater than rCut into the jetEnergy sum
                    float delta = 13    // half angle of cone formed around thrust axis to allow in particles (measured in angle)
                    ) 		
{
	std::cout<<"Initializing trees for use..."<<std::endl;
    std::string jtTreeName = "";
    if (jttree == 0) jtTreeName = "akR4ESchemeJetTree";
    if (jttree == 1) jtTreeName = "akR4WTAmodpSchemeJetTree";
    if (jttree == 2) jtTreeName = "akR8ESchemeJetTree";
    if (jttree == 3) jtTreeName = "akR8WTAmodpSchemeJetTree";
    if (jttree == 4) jtTreeName = "ktN2WTAmodpSchemeJetTree";

	if(inFileName.find(".root") != std::string::npos)
	{
      TFile* temp_p = new TFile(inFileName.c_str(), "READ");
      jtTreeName = smartJetName(jtTreeName, temp_p);
      temp_p->Close();
      delete temp_p;
    }

    TFile * output = TFile::Open(outFileName.c_str(),"recreate");

    TFile * f = TFile::Open(inFileName.c_str(),"read");
    TTree * t = (TTree*)f->Get("t"); 
  	TTree * jt = (TTree*)f->Get(smartJetName("ak4ESchemeJetTree", f).c_str()); 

  	int nParticle;		t->SetBranchAddress("nParticle",&nParticle); 
    float pmag[1000];	t->SetBranchAddress("pmag",&pmag);  
    float mass[1000];	t->SetBranchAddress("mass",&mass); 
    float theta_wrtThr[1000];   t->SetBranchAddress("theta_wrtThr",&theta_wrtThr); 
    float phi_wrtThr[1000];     t->SetBranchAddress("phi_wrtThr",&phi_wrtThr); 
    float eta_wrtThr[1000];     t->SetBranchAddress("eta_wrtThr",&eta_wrtThr);
    float TTheta;       t->SetBranchAddress("TTheta",&TTheta); 
    float TPhi;         t->SetBranchAddress("TPhi",&TPhi); 

    int nref;			jt->SetBranchAddress("nref",&nref); 
    float jtpt[100];	jt->SetBranchAddress("jtpt",&jtpt); 
    float jteta[100];	jt->SetBranchAddress("jteta",&jteta); 
    float jtphi[100];	jt->SetBranchAddress("jtphi",&jtphi); 

    TH1D * h_jtEnergyRatio = new TH1D("h_jtEnergyRatio","Dijet Energy Ratio",100,0.0,1.0);
    TH1D * h_radius = new TH1D("h_radius","Radius of tracks wrt Thrust Axis", 100, 0.0, 10.0);
    for(unsigned int nE = 0; nE < t->GetEntries(); nE++)
    {
    	t->GetEntry(nE);
        jt->GetEntry(nE);

        // look at all particles
    	Float_t E = 0;
        // look at a cone of half angle delta around thrust axis
        TVector3 vThrust(0,0,1); vThrust.SetTheta(TTheta); vThrust.SetPhi(TPhi);
        TVector3 vTrack(0,0,1);
        Float_t jetE = 0;
	// ANTHONY YOU ARE RIGHT HERE AND LOOKING AT WHY THE ANGLES ARE SO CLOSE TO 2PI AND PI AND NOT VERY SMALL???
	Float_t radius = 0.0;
	for(unsigned int nP = 0; nP < nParticle; nP++)
    	{
    	    E += sqrt(pmag[nP]*pmag[nP] +mass[nP]*mass[nP]);
            //vTrack.SetTheta(theta_wrtThr[nP]); vTrack.SetPhi(phi_wrtThr[nP]); 
	    //std::cout<<"Angle "<<fabs( vTrack.Angle(vThrust) )<<std::endl;
            //if( fabs( vTrack.Angle(vThrust) ) <= (delta*TMath::Pi()/180.0)) jetE += sqrt(pmag[nP]*pmag[nP] +mass[nP]*mass[nP]);
	    radius = sqrt(eta_wrtThr[nP]*eta_wrtThr[nP] + theta_wrtThr[nP]*theta_wrtThr[nP]);
	    h_radius->Fill(radius);
	    if ( radius > rCut ) jetE += sqrt(pmag[nP]*pmag[nP] +mass[nP]*mass[nP]); 
	}

	// std::cout<<"jetE = "<<jetE<<" total energy = "<<E<<std::endl;
    	h_jtEnergyRatio->Fill( jetE/E);
    }

    output->cd();
    h_jtEnergyRatio->Write();
    h_radius->Write();
    delete  h_radius;
    delete h_jtEnergyRatio;
    delete jt;
    delete t;
    f->Close();
    delete f;
    output->Close();
    delete output;

    return 0;
}

