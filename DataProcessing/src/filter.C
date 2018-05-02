#include <TChain.h>
#include <TFile.h>
#include <iostream>

#include "DataProcessing/include/eventData.h"

//=====================================
// filter.C
//
// To filter out interesting events
//  
// ====================================

void filter(char* fname = "/data/cmcginn/StudyMultSamples/ALEPH/MC/20180323/alephMCRecoAfterCutPaths_1994.root",char *outfile="output.root")
{

string trees[33] = {
"t",
"akR4ESchemeJetTree",
"akR4WTAmodpSchemeJetTree",
"akR8ESchemeJetTree",
"akR8WTAmodpSchemeJetTree",
"BoostedWTAR8Evt",
"ktN2ESchemeJetTree",
"ktN2WTAmodpSchemeJetTree",
"BoostedWTAktN2Evt",
"ktN3ESchemeJetTree",
"ktN3WTAmodpSchemeJetTree",
"tgen",
"tgenBefore",
"akR4ESchemeGenJetTree",
"akR4ESchemeGenJetBeforeTree",
"akR4WTAmodpSchemeGenJetTree",
"akR4WTAmodpSchemeGenJetBeforeTree",
"akR8ESchemeGenJetTree",
"akR8ESchemeGenJetBeforeTree",
"akR8WTAmodpSchemeGenJetTree",
"akR8WTAmodpSchemeGenJetBeforeTree",
"genBoostedWTAR8Evt",
"genBoostedWTAR8BeforeEvt",
"ktN2ESchemeGenJetTree",
"ktN2ESchemeGenJetBeforeTree",
"ktN2WTAmodpSchemeGenJetTree",
"ktN2WTAmodpSchemeGenJetBeforeTree",
"genBoostedWTAktN2Evt",
"genBoostedWTAktN2BeforeEvt",
"ktN3ESchemeGenJetTree",
"ktN3ESchemeGenJetBeforeTree",
"ktN3WTAmodpSchemeGenJetTree",
"ktN3WTAmodpSchemeGenJetBeforeTree"
};


   TTree* ch[33];
   TTree *originalTrees[33];
   int exists[33];
   int N = 33;
   
   TFile *inf = new TFile(fname);
   TFile *outf = new TFile(outfile, "RECREATE");
   for(int i = 0; i < N; ++i){
      exists[i]=0;
      originalTrees[i]=(TTree*) inf->Get(trees[i].data());
      cout<<"Tree loading : "<<string(trees[i]).data()<<endl;
      if (originalTrees[i]==0) continue;
      exists[i]=1;
      ch[i] = originalTrees[i]->CloneTree(0);
      cout<<"Entries : "<<((TTree*) inf->Get(trees[i].data()))->GetEntries()<<endl;
      ch[i]->SetMaxTreeSize(4000000000);
   }
   
   TTree *t = (TTree*)inf->Get("t");
   eventData data;
   std::vector<std::string> inList;
   inList.push_back("nChargedHadrons");
   data.SetStatusAndAddressRead(t,inList);
  
   for (int i=0;i<t->GetEntries();i++) {
     t->GetEntry(i);
     if (i%1000 == 0) cout <<i<<" / "<<t->GetEntries()<<endl;
     if (data.nChargedHadrons>=30) {
        for (int j=0;j<N;j++) {
	  if (exists[j]==0)continue;
	  originalTrees[j]->GetEntry(i);
	  ch[j]->Fill();
	}  
     }
   } 
   
   for (int j=0;j<N;j++){
	  if (exists[j]==0)continue;
      ch[j]->Write();
   }
   outf->ls();
   outf->Close();
/*
   TFile* file = new TFile(outfile, "RECREATE");
   
   for(int i = 0; i < N; ++i){
      file->cd();
      cout <<string(trees[i]).data()<<endl;
      //ch[i]->Merge(file,0,"keep");
      
   }
   cout <<"Good"<<endl;
   //file->Write();
   file->Close();
*/
}


