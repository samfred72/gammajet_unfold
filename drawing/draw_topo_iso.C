#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/treeuser.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/pho_object.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/reweight_utility.h"
// The original gammajet project builds its OWN, differently-laid-out drawer/ana classes
// into /home/samson72/root/lib/libgammajet.so, sitting on the same library search path
// as this project's libgammajet_unfold.so. Without forcing which one loads first, ROOT's
// implicit symbol autoload can bind drawer/ana calls to the wrong (mismatched-layout)
// library and segfault - explicit load removes the ambiguity.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Diagnostic: topo-isolation energy (cluster_showershape[11], i.e. iso_topo_04 -
// src/treeuser.h:132) for TIGHT photon clusters only. Tightness reuses the pipeline's own
// definition, pho_object::get_showershape() == 2 (src/pho_object.cc:6-25) - not
// re-derived here, so it can't drift from the production selection.
//
// pythia_Photon10 only, |v_z| < 30 cm. That vz cut is looser than the pipeline's own
// ana::vzcut = 60 (src/ana.h:51) - this plot uses the tighter 30 cm value on explicit
// request, not the standard analysis cut. No other selection (no pT/eta/pairing/purity
// cuts) - this is a raw, unselected look at the isolation shower-shape variable itself.
void draw_topo_iso() {
  gStyle->SetOptStat(0);
  const double vzcut = 30;
  const char * pdfPath = "/home/samson72/sphnx/gammajet_unfold/pdfs/topo_iso_tight_Photon10.pdf";

  treeuser tu("Photon10", "pythia");
  // MC-only vz/cluster-pT reweighting (reweight/make_vz_pt_reweight.C) - this is a
  // pythia MC sample, so every Fill below carries this event's weight.
  Reweighter rw;

  // Range covers the full observed spread for tight clusters in this sample (checked
  // to run roughly -32 to +31 GeV) - a narrower range clipped a large fraction of
  // entries into the under/overflow bins and silently truncated the plotted shape.
  TH1D * h = new TH1D("hTopoIsoTight",
      ";Topo-isolation energy, iso_{topo}^{R=0.4} [GeV];Tight clusters", 100, -35, 35);
  h->Sumw2();

  Long64_t nentries = tu.t->GetEntriesFast();
  for (Long64_t e = 0; e < nentries; e++) {
    tu.t->GetEntry(e);

    if (fabs(tu.vz) > vzcut) continue;

    int showershape = pho_object::get_showershape(tu.cluster_showershape, tu.cluster_pt);
    if (showershape != 2) continue; // tight only

    h->Fill(tu.cluster_showershape[11], rw.GetWeight(tu.vz, tu.cluster_pt));
  }

  // Sanitize NaN/Inf bin content before drawing - NaN poisons ROOT's axis auto-ranging
  // and silently produces a blank canvas.
  for (int b = 0; b <= h->GetNbinsX() + 1; b++) {
    double content = h->GetBinContent(b);
    if (std::isnan(content) || std::isinf(content)) h->SetBinContent(b, 0);
  }

  cout << "Tight clusters filled: " << h->GetEntries() << endl;

  TCanvas * c = new TCanvas("c", "", 700, 700);
  TPad * p1 = new TPad("p1", "", 0, 0, 1, 1);
  p1->Draw();
  p1->cd();
  gPad->SetLogy();
  p1->SetLeftMargin(.15);
  p1->SetBottomMargin(.13);
  gPad->SetTicks(1, 1);

  h->SetLineColor(kBlack);
  h->SetLineWidth(2);
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetRangeUser(0.5, h->GetMaximum() * 1.4);
  h->Draw("hist");

  drawer d("pythia", "nominal");
  d.drawAll({"Pythia8 #gamma+jet MC"},
      {"Tight clusters", "|v_{z}| < 30 cm"}, .55, .85, 14, gPad->GetWh() * 0.8);

  c->SaveAs(pdfPath);
  cout << "Saved " << pdfPath << endl;
}
