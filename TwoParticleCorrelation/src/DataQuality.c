//c and c++ dependencies
#include <sstream>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <fstream>

//root dependencies
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TFormula.h>
#include <TNtuple.h>
#include "TDatime.h"

//local headers
#include "include/fourier.h"
#include "include/Selection.h"
#include "DataProcessing/include/trackSelection.h"
#include "include/smartJetName.h"
#include "include/ProgressBar.h"
#include "DataQuality.h"

/********************************************************************************************************************/
// Two particle correlation analysis
//
// ridge_check_parallel.c
//
// Yen-Jie Lee, Gian Michele Innocenti, Anthony Badea
//
//
/********************************************************************************************************************/


using namespace std;

/********************************************************************************************************************/
// Data Quality macro
/********************************************************************************************************************/

int DataQuality( const std::string inFileName, 		// Input file
                       std::string outFileName,    	// Output file
      		       std::string inMixFileName="", 	// Input mix file
		       bool overwrite = false,		// if we overwrite the doWTA, doThrust and doPerp from selection.h
		       bool owThrust = false,		// overwrite the value of doThrust from selection.h
		       bool owWTA = false,		// overwrite the value of doWTA from selection.h
		       bool owPerp = false,		// overwrite the value of doPerp from selection.h
		       bool owDoGen = false,            // overwrite the value of doGen from selection.h 
		       int verbose = 1,                 // Verbose level
		       bool owAjCut=false,      // dijet selection flag
		       double _AjCut=-999,     // dijet  cut selection
		       bool ow3jetEvtCut=false,      // three jet selection flag
		       double _thirdJetCut=-999,     // three jet cut selection
		       bool owBarrel = false,            // overwrite the values related to barrel selection and barrel multiplicity 
		       int _anatyperegion=0,  // option==0 -> no eta selection on tracks, option==1 -> reject the tracks at small eta, option==2 -> reject the tracks at large rapidity
		       double _etabarrelcut=-1,  // eta cut values that defines the three regions defined above
		       int _typeEnergyBarrelSel=0,      // type of selection on the total energy inside the barrel vs total. 0=no selection, 1=rejecting events with large energy in barrel, 2= rejecting events with small energy in barrel 
		       double _etabarrelcutforEselection=2.0, //define the eta region for the Ebarrel selection
		       double _maxrelenergyinsidebarrel=0., // define the cut on the Ebarrel, _maxrelenergyinsidebarrel=Ebarrel/Etotal 
		       int _typemultiplicity=0 // 0=total charged track multiplicity, 1=charged track multiplicity in barrel
               ) 
{

   cout<<"overwrite="<<overwrite<<endl;
   cout<<"owThrust="<<owThrust<<endl;
   cout<<"owWTA="<<owWTA<<endl;
   cout<<"owPerp="<<owPerp<<endl;
   cout<<"owDoGen="<<owDoGen<<endl;
   cout<<"owAjCut="<<owAjCut<<endl;
   cout<<"AjCut="<<_AjCut<<endl;
   cout<<"ow3jetEvtCut="<<ow3jetEvtCut<<endl;
   cout<<"thirdJetCut="<<_thirdJetCut<<endl;
   cout<<"owBarrel="<<owBarrel<<endl;
   cout<<"anatyperegion="<<_anatyperegion<<endl;
   cout<<"etabarrelcut="<<_etabarrelcut<<endl;
   cout<<"typeEnergyBarrelSel="<<_typeEnergyBarrelSel<<endl;
   cout<<"etabarrelcutforEselection="<<_etabarrelcutforEselection<<endl;
   cout<<"maxrelenergyinsidebarrel="<<_maxrelenergyinsidebarrel<<endl;
   cout<<"typemultiplicity="<<_typemultiplicity<<endl;
   

    // ROOT Global setting
    TH1::SetDefaultSumw2();    TH2::SetDefaultSumw2();

    // Setup event selection (Selection) and TrackSelector (TrackSelection)
    Selection s = Selection();
    TrackSelection trackSelector;
        
    // Print general information
    if (verbose) cout <<"Input File : "<<inFileName<<endl;
    if (verbose) cout <<"Output File: "<<outFileName<<endl;
    if (verbose) cout <<"Mix File: "<<inMixFileName<<endl;

    // Check if the user overwrite the setting in Selection.h    
    if (overwrite) {
       if (verbose) cout <<"Overwrite the default value from Selection.h"<<endl;
       if (owThrust) s.doThrust = true; else s.doThrust = false;
       if (owWTA)    s.doWTA = true; else s.doWTA = false;
       if (owPerp)   s.doPerp = true; else s.doPerp = false;
       if (owDoGen)  s.doGen = true; else s.doGen = false;
       if (owAjCut) {s.doAjCut = true; s.AjCut=_AjCut;}
       else s.doAjCut = false;
       if (ow3jetEvtCut) {s.do3jetEvtCut = true; s.thirdJetCut=_thirdJetCut;}
       else s.do3jetEvtCut = false;
       if (owBarrel){
         s.anatyperegion=_anatyperegion;
         s.etabarrelcut=_etabarrelcut;
         s.typeEnergyBarrelSel=_typeEnergyBarrelSel;
         s.etabarrelcutforEselection=_etabarrelcutforEselection;
         s.maxrelenergyinsidebarrel=_maxrelenergyinsidebarrel;
         s.typemultiplicity=_typemultiplicity;
         if(owThrust==false) {std::cout<<"the barrel selection and multiplicity is defined now in the thrust axis only!!!"<<std::endl; return 0;} 
       }
    }
    
    
    if (verbose) {
      if (s.doThrust && s.doWTA) { cout <<"Can't have both doWTA and doThrust on in Selection.h!"<<endl; return 0; }
      if (s.doThrust) 	cout <<"Thrust axis analysis"<<endl;
      else if (s.doWTA) 	cout <<"WTA axis analysis"<<endl;
      else 		cout <<"Beam axis analysis"<<endl;
      if ((s.doThrust||s.doWTA) && s.doPerp) cout <<"Reference axis rotated by 90 degree"<<endl;
      if (s.do3jetEvtCut) 	cout <<"Applying three jet rejection with value="<<s.thirdJetCut<<endl;
      if (s.doGen) cout <<"This is a generator level analysis, i.e., no track or event selection applied!"<<endl;    
      if (owBarrel){
        cout <<"you are modifying the parameters of the barrel selection"<<endl;    
        cout<<"_anatyperegion="<<s.anatyperegion<<endl;
        cout<<"_etabarrelcut="<<s.etabarrelcut<<endl;
        cout<<"_typeEnergyBarrelSel="<<s.typeEnergyBarrelSel<<endl;
        cout<<"_etabarrelcutforEselection="<<s.etabarrelcutforEselection<<endl;
        cout<<"_maxrelenergyinsidebarrel="<<s.maxrelenergyinsidebarrel<<endl;
        }
    }
    
    /********************************************************************************************************************/
    // Define the output file
    /********************************************************************************************************************/
    if(createhists("savehist")) return -1;
    
    TFile * output = TFile::Open(outFileName.c_str(),"recreate");
        
    /// Initialize the trees for use
    std::cout<<"Initializing trees for use..."<<std::endl;
    std::string jtTreeName = "";
    if (s.jttree == 0) jtTreeName = "akR4ESchemeJetTree";
    if (s.jttree == 1) jtTreeName = "akR4WTAmodpSchemeJetTree";
    if (s.jttree == 2) jtTreeName = "akR8ESchemeJetTree";
    if (s.jttree == 3) jtTreeName = "akR8WTAmodpSchemeJetTree";
    if (s.jttree == 4) jtTreeName = "ktN2WTAmodpSchemeJetTree";
    
    if(inFileName.find(".root") != std::string::npos){
      TFile* temp_p = new TFile(inFileName.c_str(), "READ");
      jtTreeName = smartJetName(jtTreeName, temp_p);
      temp_p->Close();
      delete temp_p;
    }
    
    // files and variables for input
    // get the correct tree for the analysis
    std::string boostTree = "BoostedWTAR8Evt";
    
    TChain * t = new TChain("t");       			t->Add(inFileName.c_str());
    TChain * boost_t = new TChain(boostTree.c_str());       	boost_t->Add(inFileName.c_str());
    TChain * jt = new TChain(jtTreeName.c_str());       	jt->Add(inFileName.c_str());

    TPCNtupleData data(s.doBelle, s.doThrust, s.doWTA);      	data.setupTPCTree(t,boost_t,jt);   
    
    TChain * t_mix = new TChain("t"); 

    // analysis
    Int_t nevent = (Int_t)t->GetEntries();
    if(s.doOneEvent) nevent = s.numEvents;

    // Setup Progress bar
    ProgressBar Bar(cout, nevent);
    Bar.SetStyle(1);

    unsigned int entryDiv = (nevent > 200) ? nevent / 200 : 1;
 
    /********************************************************************************************************************/
    // Main Event Loop
    /********************************************************************************************************************/
    nevent=20000;
    for (Int_t i=0;i<nevent;i++) {
        t->GetEntry(i); 
        boost_t->GetEntry(i);
        jt->GetEntry(i);
	    Bar.Update(i);
        Bar.PrintWithMod(entryDiv);


        for ( Int_t j=0;j<data.particle.nParticle;j++ )
        {
	       fillHisto(data,j,0);
        }
        
        // nTrk calculation
        s.domissPCut=false;
        s.doAjCut=false;
        s.doE=false;
        s.typeEnergyBarrelSel=0;
        s.typemultiplicity=0;
        Int_t nTrk= s.ridge_eventSelection(&data.event, &data.jet, &data.particle);
        if( nTrk < 0) continue;
        
        for ( Int_t j=0;j<data.particle.nParticle;j++ )
        {
	       fillHisto(data,j,1);
        }
        
        s.typeEnergyBarrelSel=1;
        s.etabarrelcutforEselection=2.0;
        s.maxrelenergyinsidebarrel=0.4;
        s.anatyperegion=1;
        s.etabarrelcut=2.0;

        nTrk= s.ridge_eventSelection(&data.event, &data.jet, &data.particle);
        if( nTrk < 0) continue;
        
        for ( Int_t j=0;j<data.particle.nParticle;j++ )
        {
	       fillHisto(data,j,2);
        }
        
        s.typemultiplicity=1;
        Int_t nTrkBarrel= s.ridge_eventSelection(&data.event, &data.jet, &data.particle);
        if( nTrkBarrel < 0) continue;
    

        for ( Int_t j=0;j<data.particle.nParticle;j++ )
        {  
            if (!trackSelector.highPurity(&data.particle,j)&&s.doGen==0) continue;
	         fillHisto(data,j,3);
            if (s.anatyperegion==1 && fabs(data.getEta(j))>s.etabarrelcut) continue;
	         fillHisto(data,j,4);
        }      
    }
    
    output->cd();
    if(writehists("savehist")) return -1;
    delete jt;
    delete t;
    output->Close();
    delete output;
    return 0;
}

