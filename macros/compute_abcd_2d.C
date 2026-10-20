#include "../src/ana.h"
#include "TVectorD.h"

R__LOAD_LIBRARY(libgammajet_unfold.so);

// BDT score vs isolation of the leading Data photon cluster, for drawing/draw_abcd_2d.C (which adds
// the nominal ABCD regions). Calculation only.
// Selection: unfolder's maxpho, |vz| < ana::vzcut, |eta| < ana::photonEtaCut, ptLow <= pT < ptHigh
// (default 15-20 GeV); no jet pairing. Clusters without topo iso (<= -999) are skipped.
// Output: hists/abcd_2d_data.root (TH2D habcd, plus the pT range as TVectorD ptRange)
void compute_abcd_2d(double ptLow = ana::ptBinsUsed[0], double ptHigh = ana::ptBinsUsed[1],
                  bool showProgress = false) {

  TFile * f = TFile::Open(ana::path("trees/gammajet_Data.root"), "read");
  TTree * t = (TTree*)f->Get("towerntup");
  Float_t vz, cluster_pt, cluster_eta, cluster_bdt_scores[11], cluster_showershape[12];
  t->SetBranchStatus("*", 0);
  for (const char * b : {"vz", "cluster_pt", "cluster_eta", "cluster_bdt_scores", "cluster_showershape"})
    t->SetBranchStatus(b, 1);
  t->SetBranchAddress("vz", &vz);
  t->SetBranchAddress("cluster_pt", &cluster_pt);
  t->SetBranchAddress("cluster_eta", &cluster_eta);
  t->SetBranchAddress("cluster_bdt_scores", cluster_bdt_scores);
  t->SetBranchAddress("cluster_showershape", cluster_showershape);

  const double isoMin = -2, isoMax = 15;
  TH2D * h = new TH2D("habcd", ";#it{E}_{T}^{iso} (#it{R}=0.4) [GeV];Photon ID BDT score",
                      136, isoMin, isoMax, 100, 0, 1);  // 0.125 GeV x 0.01
  Long64_t n = t->GetEntries();
  for (Long64_t e = 0; e < n; e++) {
    t->GetEntry(e);
    if (showProgress && e % 100000 == 0) cout << "entry " << e << "/" << n << "\r" << flush;
    if (fabs(vz) > ana::vzcut) continue;
    if (fabs(cluster_eta) >= ana::photonEtaCut) continue;
    if (cluster_pt < ptLow || cluster_pt >= ptHigh) continue;
    float iso = cluster_showershape[11], bdt = cluster_bdt_scores[9];
    if (iso <= -999 || !std::isfinite(iso) || !std::isfinite(bdt)) continue;
    h->Fill(std::min(std::max((double)iso, isoMin + 1e-6), isoMax - 1e-6), bdt);
  }
  cout << "filled " << (Long64_t)h->GetEntries() << " clusters" << endl;

  const char * out = ana::path("hists/abcd_2d_data.root");
  TFile * fout = TFile::Open(out, "RECREATE");
  h->Write();
  TVectorD ptRange(2);
  ptRange[0] = ptLow;
  ptRange[1] = ptHigh;
  ptRange.Write("ptRange");
  fout->Close();
  cout << "Wrote " << out << endl;
}
