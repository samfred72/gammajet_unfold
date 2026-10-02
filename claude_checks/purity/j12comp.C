void j12comp(){
  const char* fs[3]={"gammajet_pythia_Jet12_long","gammajet_pythia_Jet12","gammajet_pythia_Photon10"};
  for(auto fn:fs){
    TFile*f=TFile::Open(Form("/home/samson72/sphnx/gammajet_unfold/trees/%s.root",fn)); TTree*t=(TTree*)f->Get("towerntup");
    double n=t->GetEntries();
    double nclu=t->GetEntries("cluster_pt>15&&cluster_pt<20&&abs(cluster_eta)<1.1");
    double ntru=t->GetEntries("truth_cluster_pt>15&&truth_cluster_pt<20");
    double ntruiso=t->GetEntries("truth_cluster_pt>15&&truth_cluster_pt<20&&truth_cluster_iso4<2");
    // truth balance: truth photon vs leading truth jet (R=0.4), isolated truth photons 15-20
    t->Draw("truth_jet_pt[2]/truth_cluster_pt>>hb(40,0,2)","truth_cluster_pt>15&&truth_cluster_pt<20&&truth_cluster_iso4<2&&truth_jet_pt[2]>5","goff");
    TH1*hb=(TH1*)gDirectory->Get("hb");
    // how many reco 15-20 clusters are within dR<0.1 of truth photon
    double nmatch=t->GetEntries("cluster_pt>15&&cluster_pt<20&&abs(cluster_eta)<1.1&&truth_cluster_pt>0&&sqrt(pow(cluster_eta-truth_cluster_eta,2)+pow(TVector2::Phi_mpi_pi(cluster_phi-truth_cluster_phi),2))<0.1");
    printf("%-28s entries %9.0f | reco clu 15-20: %8.0f, truth-matched %.3f | truth gamma 15-20: %8.0f (iso<2: %.3f) | per entry %.4f | <jet/gamma>_truth %.3f (N %.0f)\n",
      fn,n,nclu,nmatch/nclu,ntru,ntruiso/std::max(1.,ntru),ntru/n,hb->GetMean(),hb->GetEntries());
    t->Draw("truth_jet_pt[2]>>hj(50,0,50)","","goff"); TH1*hj=(TH1*)gDirectory->Get("hj");
    printf("     truth jet (R=0.4) pT: mean %.1f, frac>12 GeV %.3f;  events with no truth jet (>=5): %.3f\n",hj->GetMean(),hj->Integral(hj->FindBin(12.001),51)/std::max(1.,hj->Integral(0,51)), t->GetEntries("truth_jet_pt[2]<5")/n);
  }
}
