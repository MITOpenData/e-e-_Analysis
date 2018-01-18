
//local headers
#include "ridge_check.c"

void scan(TString filename,TString savename,Int_t isGen,Int_t isThrust,Int_t mult_low[], Int_t mult_high[], int n)
{
    Int_t isBelle      = 0;        // BELLE analysis = 1, CMS/ALEPH analysis = 0
    Int_t isTheta      = 0;         // Use Theta angle = 1, Use Eta = 0
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
    
    for (int i=0; i<n; i++)
    {
        analysis(filename,savename, isThrust, isBelle, isTheta, isGen, maxevt, mult_low[i], mult_high[i], nbin, verbose, num_runs, ptMin, ptMax, detaRange, ptMinForN, ptMaxForN,etaCutForN);
    }
}

void run(int MB = 0, int min60 = 0, int LEP = 0, int pythia8_12_8_17 = 0)
{
    if(MB)
    {
    // Minimum Bias //
    int n_mb = 4;
    Int_t mult_low[4] = {20,30,50};
    Int_t mult_high[4] = {9999,9999,9999};
    TString pthat1_Zee = "/data/cmcginn/GeneratorsHEPMC/PYTHIA8/0p0912/pthat1_Zee/outFile_PYTHIA8_0p0912_pthat1_Zee_MERGED.root";
    TString pthat1_Zee_save = "pthat1_Zee";
    scan(pthat1_Zee,pthat1_Zee_save,1,1,mult_low,mult_high,n_mb);
    TString pthat1_Zee_RopeWalk = "/data/cmcginn/GeneratorsHEPMC/PYTHIA8/0p0912/pthat1_Zee_RopeWalk/outFile_PYTHIA8_0p0912_pthat1_Zee_RopeWalk_MERGED.root";
    TString pthat1_Zee_RopeWalk_save = "pthat1_Zee_RopeWalk";
    scan(pthat1_Zee_RopeWalk,pthat1_Zee_RopeWalk_save,1,1,mult_low,mult_high,n_mb);
    }
    
    
    if(min60)
    {
    // nParticle min 60 //
    int n_60 = 1;
    Int_t mult_low_minNPart60[1] = {60};
    Int_t mult_high_minNPart60[1] = {9999};
    TString pthat1_Zee_minNPart60 = "/data/cmcginn/GeneratorsHEPMC/PYTHIA8/0p0912/pthat1_Zee_minNPart60/outFile_PYTHIA8_0p0912_pthat1_Zee_minNPart60_MERGED.root";
    TString pthat1_Zee_minNPart60_save = "pthat1_Zee_minNPart60";
    scan(pthat1_Zee_minNPart60,pthat1_Zee_minNPart60_save,1,1,mult_low_minNPart60,mult_high_minNPart60,n_60);
    TString pthat1_Zee_RopeWalk_minNPart60 = "/data/cmcginn/GeneratorsHEPMC/PYTHIA8/0p0912/pthat1_Zee_RopeWalk_minNPart60/outFile_PYTHIA8_0p0912_pthat1_Zee_RopeWalk_minNPart60_MERGED.root";
    TString pthat1_Zee_RopeWalk_minNPart60_save = "pthat1_Zee_RopeWalk_minNPart60";
    scan(pthat1_Zee_RopeWalk_minNPart60,pthat1_Zee_RopeWalk_minNPart60_save,1,1,mult_low_minNPart60,mult_high_minNPart60,n_60);
    }
    
    
    if(LEP)
    {
    /// LEP DATA//
    int n_LEP = 3;
    Int_t mult_low_LEP[3] = {20,30,50};
    Int_t mult_high_LEP[3] = {9999,9999,9999};
    TString LEP1 = "/home/abadea/Documents/20171022/alephDataPaths_LEP1.root";
    TString LEP1_save = "LEP1";
    scan(LEP1,LEP1_save,0,1,mult_low_LEP,mult_high_LEP,n_LEP);
    TString LEP2 = "/home/abadea/Documents/20171022/alephDataPaths_LEP2_1995to2000.root";
    TString LEP2_save = "LEP2";
    scan(LEP2,LEP2_save,0,1,mult_low_LEP,mult_high_LEP,n_LEP);
    }
    
    
    //
    if(pythia8_12_8_17)
    {
        int n_12_8_17 = 4;
        Int_t mult_low_12_8_17[4] = {20,30,40,50};
        Int_t mult_high_12_8_17[4] = {9999,9999,9999,9999};
        TString pthat1_Zee = "/data/cmcginn/GeneratorsHEPMC/PYTHIA8/0p0912/pthat1_Zee/20171208/outFile_MERGED_nEvt2000000_nMinChgPart0_RopeWalk0.root";
        TString pthat1_Zee_save = "pthat1_Zee_nMinChgPart0_RopeWalk0";
        scan(pthat1_Zee,pthat1_Zee_save,1,1,mult_low_12_8_17,mult_high_12_8_17,n_12_8_17);
        TString pthat1_Zee_RopeWalk = "/data/cmcginn/GeneratorsHEPMC/PYTHIA8/0p0912/pthat1_Zee_RopeWalk/20171208/outFile_MERGED_nEvt2000000_nMinChgPart0_RopeWalk1.root";
        TString pthat1_Zee_RopeWalk_save = "pthat1_Zee_nMinChgPart0_RopeWalk1";
        scan(pthat1_Zee_RopeWalk,pthat1_Zee_RopeWalk_save,1,1,mult_low_12_8_17,mult_high_12_8_17,n_12_8_17);
        
        int n_12_8_17_nMinChgPart30 = 3;
        Int_t mult_low_12_8_17_nMinChgPart30[3] = {30,40,50};
        Int_t mult_high_12_8_17_nMinChgPart30[3] = {9999,9999,9999};
        
        TString pthat1_Zee_30 = "/data/cmcginn/GeneratorsHEPMC/PYTHIA8/0p0912/pthat1_Zee_minChgNPart30/20171208/outFile_MERGED_nEvt20000_nMinChgPart30_RopeWalk0.root";
        TString pthat1_Zee_save_30 = "pthat1_Zee_nMinChgPart30_RopeWalk0";
        scan(pthat1_Zee_30,pthat1_Zee_save_30,1,1,mult_low_12_8_17_nMinChgPart30,mult_high_12_8_17_nMinChgPart30,n_12_8_17_nMinChgPart30);
        
        TString pthat1_Zee_30_RopeWalk = "/data/cmcginn/GeneratorsHEPMC/PYTHIA8/0p0912/pthat1_Zee_RopeWalk_minChgNPart30/20171208/outFile_MERGED_nEvt20000_nMinChgPart30_RopeWalk1.root";
        TString pthat1_Zee_save_30_RopeWalk = "pthat1_Zee_RopeWalk_minChgNPart30_RopeWalk1";
        scan(pthat1_Zee_30_RopeWalk,pthat1_Zee_save_30_RopeWalk,1,1,mult_low_12_8_17_nMinChgPart30,mult_high_12_8_17_nMinChgPart30,n_12_8_17_nMinChgPart30);
    }
}
