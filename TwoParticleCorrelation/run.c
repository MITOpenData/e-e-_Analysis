#include "ridge_check.c"

void scan(TString filename,TString savename,Int_t isGen,Int_t isThrust,Int_t mult_low[], Int_t mult_high[])
{
    for (int i=0;i<4;i++)
    {
        // check max event above 10000000
        analysis(filename,savename, 0, isThrust, 0, isGen, 1000000,mult_low[i],mult_high[i],20,0,5,0,4.0,3.2,0.4,100,2.4);
    }
}

void run()
{
    // Minimum Bias //
    Int_t mult_low[4] = {20,30,50};
    Int_t mult_high[4] = {9999,9999,9999};
    TString pthat1_Zee = "/data/cmcginn/GeneratorsHEPMC/PYTHIA8/0p0912/pthat1_Zee/outFile_MERGED.root";
    TString pthat1_Zee_save = "pthat1_Zee";
    scan(pthat1_Zee,pthat1_Zee_save,1,1,mult_low,mult_high);
    TString pthat1_Zee_RopeWalk = "/data/cmcginn/GeneratorsHEPMC/PYTHIA8/0p0912/pthat1_Zee_RopeWalk/outFile_MERGED.root";
    TString pthat1_Zee_RopeWalk_save = "pthat1_Zee_RopeWalk";
    scan(pthat1_Zee_RopeWalk,pthat1_Zee_RopeWalk_save,1,1,mult_low,mult_high);
    
    // nParticle min 60 //
    Int_t mult_low_minNPart60[1] = {60};
    Int_t mult_high_minNPart60[1] = {9999};
    TString pthat1_Zee_minNPart60 = "/data/cmcginn/GeneratorsHEPMC/PYTHIA8/0p0912/pthat1_Zee_minNPart60/outFile_MERGED.root";
    TString pthat1_Zee_minNPart60_save = "pthat1_Zee_minNPart60";
    scan(pthat1_Zee_minNPart60,pthat1_Zee_minNPart60_save,1,1,mult_low_minNPart60,mult_high_minNPart60);
    TString pthat1_Zee_RopeWalk_minNPart60 = "/data/cmcginn/GeneratorsHEPMC/PYTHIA8/0p0912/pthat1_Zee_RopeWalk_minNPart60/outFile_MERGED.root";
    TString pthat1_Zee_RopeWalk_minNPart60_save = "pthat1_Zee_RopeWalk_minNPart60";
    scan(pthat1_Zee_RopeWalk_minNPart60,pthat1_Zee_RopeWalk_minNPart60_save,1,1,mult_low_minNPart60,mult_high_minNPart60);
    
    /// LEP DATA//
    /*
    TString LEP1 = "/home/abadea/Documents/20171022/alephDataPaths_LEP1.root";
    TString LEP1_save = "LEP1";
    scan(LEP1,LEP1_save,0);
    TString LEP2 = "/home/abadea/Documents/20171022/alephDataPaths_LEP2_1995to2000.root";
    TString LEP2_save = "LEP2";
    scan(LEP2,LEP2_save,0);
    */
}
