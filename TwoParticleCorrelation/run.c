#include "ridge_check.c"

void run(     TString filename = "aaa",
              Int_t isBelle      = 0,		// BELLE analysis = 1, CMS/ALEPH analysis = 0
              Int_t isThrust     = 0, 		// Thurst Axis analysis = 1, Beam Axis analysis = 0
              Int_t isTheta      = 0, 		// Use Theta angle = 1, Use Eta = 0
	      int maxevt       = 1000000,	// Max number of events to be processed, 0 = all events
	      int mult_low     = 38,		// Lower cut on the event multiplicity
	      int mult_high    = 100,		// Upper cut on the event multiplicity
	      int nbin         = 16,		// Number of bins in the correlation function
	      bool verbose     = 0,		// Verbose mode
	      int num_runs     = 10,		// 
              double ptMin     = 0.,           // min pT of the particles used for correlation function
              double ptMax     = 4,             // max pT of the particles used for correlation function
              double detaRange = 3.2,             // deta window of the correlation function
              Float_t ptMinForN = 0.4,          // pT min for the N_{trk}^{offline} calculation (used for event classification) 
              Float_t ptMaxForN = 100,          // pT max for the N_{trk}^{offline} calculation (used for event classification)
              Float_t etaCutForN = 2.4          // eta window for the N_{trk}^{offline} calculation (used for event classification)
	     )
{
    // Input File
    /* TString filename; */
    
    //BELLE data
    //filename = "/data/flowex/Datasamples/Belle_Data/TPCNtuple-HadronBJ-e07-00.root";
    
    // BELLE MC
    //filename = "/data/flowex/MCsamples/Belle_MC/TPCNtuple-bbmc-e07-00.root";
    //filename = "/data/flowex/MCsamples/Belle_MC/TPCNtuple-qqmc-e07-00.root";
    //filename = "/data/flowex/MCsamples/Belle_MC/TPCNtuple-qq+bbmc-e07-00.root";
    
    //ALEPH data
    //filename = "cleaned_ALEPH_Data-all.aleph.root";

       // With new Thrust code from Austin
    filename = "cleaned_ALEPH_Data2-v3_Aug11_2017.root";


   analysis(filename, isBelle, isThrust, isTheta,maxevt,mult_low,mult_high,nbin,verbose,num_runs,ptMin,ptMax,detaRange,ptMinForN,ptMaxForN,etaCutForN); 
}

void scan()
{
 for (int i=1;i<50;i++) {
   run("",0,0,0,i);
 }
}
