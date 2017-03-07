void MCvalidation(){

TString filetrue="cleaned_LEP2MCGGBBY2000E207_mctrue_aftercut-001.aleph.root";
TString filereco="cleaned_LEP2MCGGBBY2000E207_recons_aftercut-001.aleph.root";
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
