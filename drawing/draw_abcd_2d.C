#include "../src/ana.h"
#include "../src/drawer.h"
#include "TVectorD.h"

R__LOAD_LIBRARY(libgammajet_unfold.so);

// BDT score vs isolation of the leading Data photon cluster with the nominal ABCD regions, from
// macros/compute_abcd_2d.C's hists/abcd_2d_data.root. Drawing only.
// Output: pdfs/abcd_2d_data.pdf
void draw_abcd_2d() {
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);
  const char * in = ana::path("hists/abcd_2d_data.root");
  TFile * fin = TFile::Open(in);
  if (!fin || fin->IsZombie()) {
    cout << "ERROR: cannot open " << in << " - run macros/compute_abcd_2d.C first." << endl;
    return;
  }
  TH2D * h = (TH2D*)fin->Get("habcd");
  TVectorD * ptRange = (TVectorD*)fin->Get("ptRange");
  const double ptLow = (*ptRange)[0], ptHigh = (*ptRange)[1];
  const double isoMin = h->GetXaxis()->GetXmin(), isoMax = h->GetXaxis()->GetXmax();

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
    b->SetLineColor(kBlack);
    b->SetLineWidth(4);
    b->Draw();
    TLatex * l = new TLatex(x1 + 0.35, y2 - 0.035, name);
    l->SetTextFont(62);
    l->SetTextSize(0.06);
    l->SetTextAlign(13);
    l->SetTextColor(kBlack);
    l->Draw();
  };
  region(isoMin, bdtA, isoA, bdtAhi, "A");
  region(isoB,   bdtA, isoMax, bdtAhi, "B");
  region(isoMin, bdtC, isoA, bdtChi, "C");
  region(isoB,   bdtC, isoMax, bdtChi, "D");

  // White panel for the label block inside region D, labelDrop below its top edge in BDT score
  // (the text moves by the same distance in NDC).
  const double labelDrop = 0.1;
  const double labelDropNDC = labelDrop * (1 - c->GetBottomMargin() - c->GetTopMargin());
  TBox * panel = new TBox(7.15, bdtChi - 0.232 - labelDrop, isoMax - 0.25, bdtChi - 0.015 - labelDrop);
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
  d.drawAll({ana::dataLabel()},
            {Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", ptLow, ptHigh),
             Form("|#eta^{#gamma}| < %.1f, leading cluster", ana::photonEtaCut)},
            .567, .574 - labelDropNDC, 18, 700);
  c->RedrawAxis();

  const char * out = ana::path("pdfs/abcd_2d_data.pdf");
  c->SaveAs(out);
  cout << "Wrote " << out << endl;
}
