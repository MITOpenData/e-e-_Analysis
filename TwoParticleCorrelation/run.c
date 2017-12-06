#include "ridge_check.c"

void scan(TString filename,TString savename,Int_t isGen,Int_t isThrust,Int_t mult_low[], Int_t mult_high[])
{
    Int_t isBelle      = 0;        // BELLE analysis = 1, CMS/ALEPH analysis = 0
    Int_t isTheta      = 0;         // Use Theta angle = 1, Use Eta = 0
    Int_t isGen       = 0;         // Use isGen = 1 if gen level so no cut on charged particles, 0 if not gen level
    Int_t maxevt       = 1000000;    // Max number of events to be processed, 0 = all events
    Int_t nbin         = 20;        // Number of bins in the correlation function
    bool verbose     = 0;        // Verbose mode
    Int_t num_runs     = 5;        //
    double ptMin     = 0.4;           // min pT of the particles used for correlation function
    double ptMax     = 100;             // max pT of the particles used for correlation function
    double detaRange = 3.2;             // deta window of the correlation function
    Float_t ptMinForN = 0.4;          // pT min for the N_{trk}^{offline} calculation (used for event classification)
    Float_t ptMaxForN = 100;          // pT max for the N_{trk}^{offline} calculation (used for event classification)
    Float_t etaCutForN = 1.8;          // eta window for the N_{trk}^{offline} calculation (used for event classification)
    
    int n = sizeof(mult_low)/sizeof(mult_low[0]);
    for (int i=0; i<n; i++)
    {
        analysis(filename,savename, isBelle, isThrust, isTheta, isGen, maxevt, mult_low[i], mult_high[i], nbin, verbose, num_runs, ptMin, ptMax, detaRange, ptMinForN, ptMaxForN,etaCutForN);
    }
}

void run()
{
    // Minimum Bias //
    Int_t mult_low[4] = {20,30,50};
    Int_t mult_high[4] = {9999,9999,9999};
    TString pthat1_Zee = "/data/cmcginn/GeneratorsHEPMC/PYTHIA8/0p0912/pthat1_Zee/outFile_PYTHIA8_0p0912_pthat1_Zee_MERGED.root";
    TString pthat1_Zee_save = "pthat1_Zee";
    scan(pthat1_Zee,pthat1_Zee_save,1,1,mult_low,mult_high);
    TString pthat1_Zee_RopeWalk = "/data/cmcginn/GeneratorsHEPMC/PYTHIA8/0p0912/pthat1_Zee_RopeWalk/outFile_PYTHIA8_0p0912_pthat1_Zee_RopeWalk_MERGED.root";
    TString pthat1_Zee_RopeWalk_save = "pthat1_Zee_RopeWalk";
    scan(pthat1_Zee_RopeWalk,pthat1_Zee_RopeWalk_save,1,1,mult_low,mult_high);
    
    // nParticle min 60 //
    Int_t mult_low_minNPart60[1] = {60};
    Int_t mult_high_minNPart60[1] = {9999};
    TString pthat1_Zee_minNPart60 = "/data/cmcginn/GeneratorsHEPMC/PYTHIA8/0p0912/pthat1_Zee_minNPart60/outFile_PYTHIA8_0p0912_pthat1_Zee_minNPart60_MERGED.root";
    TString pthat1_Zee_minNPart60_save = "pthat1_Zee_minNPart60";
    scan(pthat1_Zee_minNPart60,pthat1_Zee_minNPart60_save,1,1,mult_low_minNPart60,mult_high_minNPart60);
    TString pthat1_Zee_RopeWalk_minNPart60 = "/data/cmcginn/GeneratorsHEPMC/PYTHIA8/0p0912/pthat1_Zee_RopeWalk_minNPart60/outFile_PYTHIA8_0p0912_pthat1_Zee_RopeWalk_minNPart60_MERGED.root";
    TString pthat1_Zee_RopeWalk_minNPart60_save = "pthat1_Zee_RopeWalk_minNPart60";
    scan(pthat1_Zee_RopeWalk_minNPart60,pthat1_Zee_RopeWalk_minNPart60_save,1,1,mult_low_minNPart60,mult_high_minNPart60);
    
    /// LEP DATA//
    Int_t mult_low_LEP[3] = {20,30,50};
    Int_t mult_high_LEP[3] = {9999,9999,9999};
    TString LEP1 = "/home/abadea/Documents/20171022/alephDataPaths_LEP1.root";
    TString LEP1_save = "LEP1";
    scan(LEP1,LEP1_save,0,1,mult_low_LEP,mult_high_LEP);
    TString LEP2 = "/home/abadea/Documents/20171022/alephDataPaths_LEP2_1995to2000.root";
    TString LEP2_save = "LEP2";
    scan(LEP2,LEP2_save,0,1,mult_low_LEP,mult_high_LEP);
    
}
