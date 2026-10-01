#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"

R__LOAD_LIBRARY(libgammajet_unfold.so);

// BDT score vs isolation of the leading photon cluster in Data, with the nominal ABCD
// regions drawn on top - a clean version of the ABCD-method figure for talks.
//
// Selection: the leading cluster as unfolder.cc builds maxpho (iso = iso_topo_04 =
// cluster_showershape[11], bdt = cluster_bdt_scores[9]), |vz| < ana::vzcut,
// |eta| < ana::photonEtaCut, and ptLow <= pT < ptHigh (default: the first reported pT
// bin, 15-20 GeV). Clusters
// without a computed topo isolation (iso <= -999, the PhotonClusterBuilder sentinel) are
// skipped, as ana::findabcdBin does. No jet pairing is required here (that needs the
// jet calibration chain), so this is the photon selection only; the region boundaries
// are read straight from ana.h (nominal bin 0), not hard-coded.
//
// Output: pdfs/abcd_2d_data.pdf
void draw_abcd_2d(double ptLow = ana::ptBinsUsed[0], double ptHigh = ana::ptBinsUsed[1],
                  bool showProgress = false) {
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);

  TFile * f = TFile::Open("/home/samson72/sphnx/gammajet_unfold/trees/gammajet_Data.root", "read");
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

  TCanvas * c = new TCanvas("c", "", 900, 700);
  c->SetLeftMargin(0.12);
  c->SetRightMargin(0.14);
  c->SetBottomMargin(0.13);
  c->SetTopMargin(0.05);
  c->SetLogz();
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.045);
  h->GetYaxis()->SetLabelSize(0.045);
  h->GetXaxis()->SetTitleOffset(1.1);
  h->GetYaxis()->SetTitleOffset(1.1);
  h->GetZaxis()->SetTitle("Clusters");
  h->GetZaxis()->SetTitleSize(0.045);
  h->GetZaxis()->SetLabelSize(0.04);
  h->Draw("colz");

  // nominal ABCD boundaries (ana.h, bin 0)
  const double isoA = ana::isoBins[0], isoB = ana::isoBinsHigh[0];
  const double bdtA = ana::bdtGoodLow[0], bdtAhi = ana::bdtGoodHigh[0];
  const double bdtC = ana::bdtBadLow[0], bdtChi = ana::bdtBadHigh[0];

  auto region = [](double x1, double y1, double x2, double y2, const char * name) {
    TBox * b = new TBox(x1, y1, x2, y2);
    b->SetFillStyle(0);
    b->SetLineColor(kRed + 1);
    b->SetLineWidth(4);
    b->Draw();
    TLatex * l = new TLatex(x1 + 0.35, y2 - 0.035, name);
    l->SetTextFont(62);
    l->SetTextSize(0.06);
    l->SetTextAlign(13);
    l->SetTextColor(kRed + 1);
    l->Draw();
  };
  region(isoMin, bdtA, isoA, bdtAhi, "A");
  region(isoB,   bdtA, isoMax, bdtAhi, "B");
  region(isoMin, bdtC, isoA, bdtChi, "C");
  region(isoB,   bdtC, isoMax, bdtChi, "D");

  // white panel for the label block, inside region D (upper right of it)
  TBox * panel = new TBox(7.4, bdtChi - 0.191, isoMax - 0.25, bdtChi - 0.015);
  panel->SetFillColor(kWhite);
  panel->SetFillStyle(1001);  // TBox defaults to hollow
  panel->SetLineWidth(0);
  panel->Draw();
  TBox * frame = (TBox*)panel->Clone();  // a filled TBox draws no outline: add one
  frame->SetFillStyle(0);
  frame->SetLineColor(kBlack);
  frame->SetLineWidth(2);
  frame->Draw();
  drawer d;
  d.drawAll({"p+p Run24 Data"},
            {Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", ptLow, ptHigh),
             Form("|#eta^{#gamma}| < %.1f, leading cluster", ana::photonEtaCut)},
            .58, .588, 18, 700);
  c->RedrawAxis();

  const char * out = "/home/samson72/sphnx/gammajet_unfold/pdfs/abcd_2d_data.pdf";
  c->SaveAs(out);
  cout << "Wrote " << out << endl;
}
