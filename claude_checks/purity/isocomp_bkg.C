R__LOAD_LIBRARY(libgammajet_unfold.so);
#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
// iso4 (showershape[11], topo iso used by pho_object) for tight-BDT (0.8-1.0) clusters, |eta|<1.1, |vz|<60
void isocomp_bkg(){
  struct S { const char* f; const char* sel; };
  const char* base = "abs(vz)<60 && abs(cluster_eta)<1.1 && cluster_bdt_scores[9]>0.2 && cluster_bdt_scores[9]<0.6 && cluster_showershape[11]>-999";
  const char* match = "truth_cluster_pt>0 && sqrt(pow(cluster_eta-truth_cluster_eta,2)+pow(TVector2::Phi_mpi_pi(cluster_phi-truth_cluster_phi),2))<0.1";
  double bins[3][2]={{15,20},{20,25},{25,35}};
  TFile*fd=TFile::Open("/home/samson72/sphnx/gammajet_unfold/trees/gammajet_Data.root"); TTree*td=(TTree*)fd->Get("towerntup");
  TFile*f10=TFile::Open("/home/samson72/sphnx/gammajet_unfold/trees/gammajet_pythia_Photon10.root"); TTree*t10=(TTree*)f10->Get("towerntup");
  TFile*f20=TFile::Open("/home/samson72/sphnx/gammajet_unfold/trees/gammajet_pythia_Photon20.root"); TTree*t20=(TTree*)f20->Get("towerntup");
  TFile*fj=TFile::Open("/home/samson72/sphnx/gammajet_unfold/trees/gammajet_pythia_Jet12_long.root"); TTree*tj=(TTree*)fj->Get("towerntup");
  printf("NON-TIGHT (BDT 0.2-0.6, regions C+D) clusters: iso4 [GeV] -- fraction <2 (A-like) | 2-4 (gap) | >4 (B-like) | median\n");
  for(auto&b:bins){
    TTree* tsig = b[0]<20 ? t10 : t20; // Photon10: truth 12-24, Photon20: truth > 24 (treeuser.h windows)
    TString pt=Form("cluster_pt>%g&&cluster_pt<%g",b[0],b[1]);
    TString tw = b[0]<20 ? "truth_cluster_pt>12&&truth_cluster_pt<24" : (b[0]<25 ? "truth_cluster_pt>12&&truth_cluster_pt<24" : "truth_cluster_pt>24");
    if (b[0]>=20 && b[0]<25) tsig=t10;
    struct Q{const char*name; TTree*t; TString sel;};
    Q qs[4]={{"Data (all non-tight)",td,TString(base)+"&&"+pt},
             {"gamma+jet MC, truth-tagged",tsig,TString(base)+"&&"+pt+"&&"+match+"&&"+tw},
             {"Jet12 MC, truth-tagged",tj,TString(base)+"&&"+pt+"&&"+match},
             {"Jet12 MC, untagged (bkg)",tj,TString(base)+"&&"+pt+"&&!("+match+")"}};
    printf("--- %g-%g GeV ---\n",b[0],b[1]);
    for(auto&q:qs){
      q.t->Draw("cluster_showershape[11]>>hi(240,-4,20)",q.sel,"goff"); TH1*h=(TH1*)gDirectory->Get("hi");
      double n=h->Integral(0,241); if(n<=0){printf("%-30s empty\n",q.name);continue;}
      double lo=h->Integral(0,h->FindBin(2-1e-4))/n, gap=h->Integral(h->FindBin(2+1e-4),h->FindBin(4-1e-4))/n, hi=h->Integral(h->FindBin(4+1e-4),241)/n;
      double qv[1],pr[1]={0.5}; h->GetQuantiles(1,qv,pr);
      printf("%-30s N=%8.0f  <2: %.3f  2-4: %.3f  >4: %.3f  median %.2f\n",q.name,n,lo,gap,hi,qv[0]);
      delete h;
    }
  }
}
