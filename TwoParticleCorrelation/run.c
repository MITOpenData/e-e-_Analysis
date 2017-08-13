#include "ridge_check.c"

void run(int isBelle      = 0,		// 
              int isThrust     = 0, 		//
	      int maxevt       = 1000000,	// Max number of events to be processed, 0 = all events
	      int mult_low     = 38,		// Lower cut on the event multiplicity
	      int mult_high    = 100,		// Upper cut on the event multiplicity
	      int nbin         = 16,		// Number of bins in the correlation function
	      bool verbose     = 0,		// Verbose mode
	      int num_runs     = 10,		// 
              double ptMin     = 0.1,           // min pT of the particles used for correlation function
              double ptMax     = 4,             // max pT of the particles used for correlation function
              double detaRange = 3.5              //  deta window
	     )
{
   analysis(isBelle, isThrust,maxevt,mult_low,mult_high,nbin,verbose,num_runs,ptMin,ptMax,detaRange); 
}
