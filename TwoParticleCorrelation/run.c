#include "ridge_check.c"

void scan(TString filename = "/data/flowex/CMSsample/TPCNtuple_MinBias_TuneCUETP8M1_5p02TeV-pythia8-HINppWinter16DR-NoPU_75X_mcRun2_asymptotic_ppAt5TeV_forest_v2_track.root",)
{
    Int_t mult_low[3] = {0,20,30};
    Int_t mult_hig[3] = {20,30,9999};
    
    for (int i=0;i<3;i++)
    {
        // check max event above 10000000
        analysis(filename, 0, 0, 0, 1000000,mult_low[i],mult_high[i],20,0,5,0,4.0,3.2,0.4,100,2.4);
    }
}

void run()
{
    TString pthat1_Zee = "/data/cmcginn/GeneratorsHEPMC/PYTHIA8/0p0912/pthat1_Zee/outFile_MERGED.root";
    scan(pthat1_Zee);
    TString pthat1_Zee_RopeWalk = "/data/cmcginn/GeneratorsHEPMC/PYTHIA8/0p0912/pthat1_Zee_RopeWalk/outFile_MERGED.root";
    scan(pthat1_Zee_RopeWalk);
}
