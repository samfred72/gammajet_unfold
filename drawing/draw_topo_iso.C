#include "../src/ana.h"
#include "../src/drawer.h"
#include "../src/treeuser.h"
#include "../src/pho_object.h"
#include "../src/reweight_utility.h"
// Load explicitly: the sibling gammajet project's libgammajet.so has same-named classes.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Topo isolation (iso_topo_04) of tight photon clusters (pho_object::get_showershape() == 2),
// pythia Photon10, |v_z| < 30 cm (tighter than ana::vzcut, on request), no other selection.
void draw_topo_iso() {
  gStyle->SetOptStat(0);
  const double vzcut = 30;
  const char * pdfPath = ana::path("pdfs/topo_iso_tight_Photon10.pdf");

  treeuser tu("Photon10", "pythia");
  // MC vz/cluster-pT weight.
  Reweighter rw;

  // Covers the full spread (about -32 to +31 GeV) so nothing lands in under/overflow.
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

  // Sanitize NaN/Inf before drawing.
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
