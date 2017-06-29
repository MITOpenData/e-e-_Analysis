#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include "TLorentzVector.h"


/** 
The process: 
 
 STEP 1
 Go into loopMC.sh and edit the file so that the top is as follows:
 folder="LEP2_MC/GGUD"
 filelist="LEP2MCGGUDY1997E183_mctrue_aftercut-001.list"
 output="cleaned_LEP2MCGGUDY1997E183_mctrue_aftercut-001.aleph"
 suffix="LEP2MCGGUDY1997E183_mctrue_aftercut-001"   

This means that the aleph monte-carlo file has to be in the specified folder and that there should exist a filelist in LEP2 with the name of aleph file in it. Then it outputs the cleaned file. 
 
 STEP 2
 Replace the infile into scanMC with the correct file from LEP2 and run scanMC. This will output a root file into ROOTFiles
 
 */

using namespace std;

class particleData
{
public:
    int nParticle_gen;
    int EventNo_gen;
    int RunNo_gen;
    float Energy_gen;
    
    Float_t px_gen[100000];
    Float_t py_gen[100000];
    Float_t pz_gen[100000];
    Float_t pt_gen[100000];
    Float_t eta_gen[100000];
    Float_t phi_gen[100000];
    Float_t mass_gen[100000];
    Float_t charge_gen[10000];
    Float_t pwflag_gen[10000];
    Float_t pid_gen[10000];
};

void scanMC (TString infile="/data/flowex/MCsamples/LEP2_MC_MAIN/LEP2_MC/GGUD/cleaned_LEP2MCGGUDY1997E183_mctrue_beforecut-001.aleph"){
    
    FILE *fp=fopen(Form("%s",infile.Data()),"r");
    
    float _px,_py,_pz,_m,_charge,_pwflag,_pid;
    
    TFile *hf = new TFile(Form("ROOTfiles/%s.root",infile.Data()), "RECREATE" );
    TTree *tout = new TTree("tgen","");
    TLorentzVector v;
    cout<<"step1"<<endl;
    particleData pData;
    tout->Branch("nParticle_gen", &pData.nParticle_gen,"nParticle_gen/I");
    tout->Branch("EventNo_gen", &pData.EventNo_gen,"EventNo_gen/I");
    tout->Branch("RunNo_gen", &pData.RunNo_gen,"RunNo_gen/I");
    tout->Branch("Energy_gen", &pData.Energy_gen,"Energy_gen/F");
    tout->Branch("px_gen", pData.px_gen,"px_gen[nParticle_gen]/F");
    tout->Branch("py_gen", pData.py_gen,"py_gen[nParticle_gen]/F");
    tout->Branch("pz_gen", pData.pz_gen,"pz_gen[nParticle_gen]/F");
    tout->Branch("pt_gen", pData.pt_gen,"pt_gen[nParticle_gen]/F");
    tout->Branch("mass_gen", pData.mass_gen,"mass_gen[nParticle_gen]/F");
    tout->Branch("eta_gen", pData.eta_gen,"eta_gen[nParticle_gen]/F");
    tout->Branch("phi_gen", pData.phi_gen,"phi_gen[nParticle_gen]/F");
    tout->Branch("charge_gen", pData.charge_gen,"charge_gen[nParticle_gen]/F");
    tout->Branch("pwflag_gen", pData.pwflag_gen,"pwflag_gen[nParticle_gen]/F");
    tout->Branch("pid_gen", pData.pid_gen,"pid_gen[nParticle_gen]/F");
    
    int counterEntries=0;
    int counterParticles=0;
    
    cout<<"step2"<<endl;
    while(fscanf(fp,"%f %f %f %f %f %f %f",&_px,&_py,&_pz,&_m,&_charge,&_pwflag,&_pid)!=EOF) {
        if (_px==-999.&&_py==-999.&&_pz==-999.) {
            cout<<"counterEntries="<<counterEntries<<endl;
            pData.nParticle_gen=counterParticles;
            if(counterEntries>0) tout->Fill();
            pData.RunNo_gen=_charge;
            pData.EventNo_gen=_pwflag;
            pData.Energy_gen=_pid;
            counterParticles=0;
            cout<<"------------------------------------------------------------------------"<<endl;
            continue;
        }
        cout<<_px<<"    "<<_py<<"    "<<_pz<<"    "<<_charge<<"   "<<_pwflag<<"   "<<_pid<<endl;
        pData.px_gen[counterParticles]=_px;
        pData.py_gen[counterParticles]=_py;
        pData.pz_gen[counterParticles]=_pz;
        pData.mass_gen[counterParticles]=_m;
        v.SetXYZM(_px,_py,_pz,_m);
        pData.pt_gen[counterParticles]=v.Pt();
        pData.eta_gen[counterParticles]=v.PseudoRapidity();
        pData.phi_gen[counterParticles]=v.Phi();
        pData.charge_gen[counterParticles]=_charge;
        pData.pwflag_gen[counterParticles]=_pwflag;
        pData.pid_gen[counterParticles]=_pid;
        counterParticles++;
        counterEntries++;
    }
    hf->Write();
}


void checkMC(){
    
    // FILE LOCATION NEEDS TO BE FIXED TO svmithi02 location
    // NOT ALREADY FIXED BECAUSE CURRENT FILE NOT IN SPECIFIED LOCATION
    TFile *f = new TFile("/data/flowex/Datasamples/LEP2_MAIN/ROOTfiles/cleaned_LEP2MCGGBBY2000E207_mctrue_beforecut-001.aleph.root");
    TTree *t1 = (TTree*)f->Get("tgen");
    Int_t nParticle_gen,EventNo_gen,RunNo_gen;
    Float_t Energy_gen;
    Float_t px_gen[100000];
    Float_t py_gen[100000];
    Float_t pz_gen[100000];
    Float_t pt_gen[100000];
    Float_t eta_gen[100000];
    Float_t pid_gen[100000];
    Float_t phi_gen[100000];
    Float_t mass_gen[100000];
    
    t1->SetBranchAddress("nParticle_gen",&nParticle_gen);
    t1->SetBranchAddress("EventNo_gen",&EventNo_gen);
    t1->SetBranchAddress("RunNo_gen",&RunNo_gen);
    t1->SetBranchAddress("Energy_gen",&Energy_gen);
    t1->SetBranchAddress("px_gen",px_gen);
    t1->SetBranchAddress("py_gen",py_gen);
    t1->SetBranchAddress("pz_gen",pz_gen);
    t1->SetBranchAddress("pt_gen",pt_gen);
    t1->SetBranchAddress("eta_gen",eta_gen);
    t1->SetBranchAddress("pid_gen",pid_gen);
    t1->SetBranchAddress("phi_gen",phi_gen);
    t1->SetBranchAddress("mass_gen",mass_gen);
    
    Int_t nevent = (Int_t)t1->GetEntries();
    
    for (Int_t i=0;i<nevent;i++) {
        t1->GetEntry(i);
        
        int nparticles = nParticle_gen;
        cout<<"------------------------------------------------------------------------"<<endl;
        cout<<"nparticles="<<nparticles<<endl;
        
        for ( int j=0;j<nparticles;j++ ) {
            
            int pid1 = pid_gen[j];
            // if (pid1!=PION&&pid1!=PROTON&&pid1!=KAON) continue;
            float eta1 = eta_gen[j];
            float phi1 = phi_gen[j];
            float px1 = px_gen[j];
            float py1 = py_gen[j];
            float pz1 = pz_gen[j];
            float pt1 = pt_gen[j];
            float mass1 = mass_gen[j];
            cout<<"EventNo_gen="<<EventNo_gen<<" px_gen="<<px1<<"py_gen="<<py1<<"pz_gen="<<pz1<<"pt_gen="<<pt1<<"mass_gen="<<mass1<<endl;
        }
    }
}



int main(int argc, char *argv[])
{
    if((argc != 2))
    {
        std::cout << "Wrong number of inputs" << std::endl;
        return 1;
    }
    
    if(argc == 2)
        std::cout << "Running for the sake of god" << std::endl;
    scanMC(argv[1]);
    return 0;
}
