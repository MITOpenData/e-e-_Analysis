void MCvalidation(){

TString filetrue="/data/flowex/MCsamples/LEP2_MC_MAIN/ROOTfiles_MCcleaned_LEP2MCGGBBY2000E207_mctrue_aftercut-001.aleph.root";
TString filereco="/data/flowex/MCsamples/LEP2_MC_MAIN/ROOTfiles_MCcleaned_LEP2MCGGBBY2000E207_recons_aftercut-001.aleph.root";
TString fileoutput="test.root";

TFile*ftrue=new TFile(filetrue.Data());
TTree*ttrue=(TTree*)ftrue->Get("tgen");
TFile*freco=new TFile(filereco.Data());
TTree*treco=(TTree*)freco->Get("t");
ttrue->AddFriend(treco);

ttrue->Draw("nParticle_gen:nParticle>>h2d");
//ttrue->Draw("Sum$(pt):Sum$(pt_gen)");

TFile*foutput=new TFile(fileoutput.Data(),"recreate");
foutput->cd();
ttrue->Write();
treco->Write();
foutput->Close();

}
