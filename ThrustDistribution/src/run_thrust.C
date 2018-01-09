#include "thrust_distribution.c"

void scan(TString filename,std::string dataname,Float_t min_E,Float_t max_E,Int_t mult_low[], Int_t mult_high[])
{
    for (int i=0;i<4;i++)
    {
        thrust_distribution(filename,dataname,0.3,0.0,3.5,min_E,max_E,mult_low[i],mult_high[i],0,0);
    }
}

void run_thrust()
{
    Int_t mult_low[4] = {20,30,50,60};
    Int_t mult_high[4] = {9999,9999,9999,9999};

    /// LEP DATA//

    TString LEP1 = "~/Downloads/StudyMult-backup/TwoParticleCorrelation/alephDataPaths_LEP1.root";
    std::string LEP1_save = "LEP1";
    Float_t LEP1_E[2] = {91.,91.5};
    scan(LEP1,LEP1_save,LEP1_E[0],LEP1_E[1],mult_low,mult_high);

    std::string LEP2 = "~/Downloads/StudyMult-backup/alephDataPaths_LEP2_1995to2000.root";
    std::string LEP2_save = "LEP2";
    Float_t LEP2_E[2] = {188.,190.};
    scan(LEP2,LEP2_save,LEP2_E[0],LEP2_E[1],mult_low,mult_high);
    
}
