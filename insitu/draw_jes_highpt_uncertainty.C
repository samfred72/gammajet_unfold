#include "../src/ana.h"
#include "../src/insitu_utility.h"
#include <string>
#include <vector>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TBox.h"
#include "TH1F.h"
using namespace std;

R__LOAD_LIBRARY(libgammajet_unfold.so);

// pT-dependent data-to-MC JES uncertainty for R=0.4 jets: the in-situ gamma+jet p_a uncertainty
// (draw_jes_summary.C) as the anchor, extended in pT with the two terms of the JES leakage note
// (Nagle, Liechty, Huang, "Data-to-MC Single-Hadron Response from the T-1044 Test Beam and its
// Propagation to the sPHENIX Jet Energy Scale", v3.0 draft, Sec. 8.4 / Fig. 23):
//   A: in-situ anchor. The note assumes a flat 2%; here it is this repo's total (stat+syst) p_a
//      uncertainty divided by p_a, flat in pT.
//   B: T-1044 single-hadron data/MC ratio, +-1% tilt envelope propagated through Pythia8 jets.
//   C: +-50% mis-modelling of the sPHENIX (aluminium inner HCal) single-pion leakage.
// B and C are re-anchored in the note to vanish in the 15-25 GeV in-situ window, so they add on top
// of A. They are only defined for anti-kt R=0.4, |eta^jet|<0.7, over the note's 10.5-91 GeV points;
// nothing is extrapolated outside that range.
// Combination: B and C are symmetric envelopes and enter both totals with their full magnitude.
// A is asymmetric and sign-split as in draw_systematics.C. Data jet pT is divided by p_a, so a lower
// p_a means a higher jet scale: the JES "up" total uses p_a's low-side error, "down" the high-side.
// Run draw_jes_summary.C first.

const char * insitu_output_dir = ana::path("insitu/output");
const char * insitu_pdf_dir    = ana::path("insitu/pdfs");

// B and C (percent, symmetric half-widths) at the note's jet pT points, read off the vector
// graphics of Fig. 23 (identical to the Fig. 18 band for B). With A = 2% these reproduce the
// note's total to 0.001%, including 3.05% at 91 GeV.
const int nNotePts = 13;
const double notePt[nNotePts]   = {10.5, 13.5, 17.0, 21.5, 26.5, 32.0, 38.5, 46.0, 54.0, 62.5, 72.0, 82.0, 91.0};
const double noteB[nNotePts]    = {0.839, 0.457, 0.177, 0.115, 0.453, 0.706, 0.960, 1.241, 1.467, 1.650, 1.817, 1.979, 2.090};
const double noteC[nNotePts]    = {0.274, 0.157, 0.068, 0.044, 0.157, 0.263, 0.378, 0.500, 0.610, 0.712, 0.810, 0.898, 0.971};

void draw_jes_highpt_uncertainty() {
  const int ir = 2; // R=0.4, the only radius the note propagates
  const char * summaryPath = Form("%s/draw_jes_summary.root", insitu_output_dir);
  TFile * fin = TFile::Open(summaryPath);
  if (!fin || fin->IsZombie()) {
    cout << "ERROR: cannot open " << summaryPath << " - run draw_jes_summary.C first." << endl;
    return;
  }
  TTree * t = (TTree*)fin->Get("jes_summary");
  int tIr; float pa, totalLow, totalHigh;
  t->SetBranchAddress("ir", &tIr);
  t->SetBranchAddress("pa", &pa);
  t->SetBranchAddress("totalLow", &totalLow);
  t->SetBranchAddress("totalHigh", &totalHigh);
  bool found = false;
  for (Long64_t i = 0; i < t->GetEntries(); i++) {
    t->GetEntry(i);
    if (tIr == ir) { found = true; break; }
  }
  fin->Close();
  if (!found) {
    cout << "ERROR: no R=" << ana::JetRs[ir] << " entry in " << summaryPath << endl;
    return;
  }

  const double anchorUp   = 100.0 * totalLow  / pa;
  const double anchorDown = 100.0 * totalHigh / pa;
  cout << Form("In-situ anchor R=%.1f: p_a = %.4f -%.4f/+%.4f -> JES +%.3f%% / -%.3f%%",
               ana::JetRs[ir], pa, totalLow, totalHigh, anchorUp, anchorDown) << endl;

  TGraphAsymmErrors * gTotal = new TGraphAsymmErrors(nNotePts);
  gTotal->SetName("gJESUnc_total");
  gTotal->SetTitle("Total data-to-MC JES uncertainty [%] vs jet p_{T}, R=0.4");
  TGraph * gAup = new TGraph(nNotePts), * gAdn = new TGraph(nNotePts);
  TGraph * gBup = new TGraph(nNotePts), * gBdn = new TGraph(nNotePts);
  TGraph * gCup = new TGraph(nNotePts), * gCdn = new TGraph(nNotePts);
  TGraph * gTup = new TGraph(nNotePts), * gTdn = new TGraph(nNotePts);
  cout << Form("%8s %8s %8s %8s %8s %8s", "pT", "A+", "A-", "B", "C", "total +/-") << endl;
  for (int i = 0; i < nNotePts; i++) {
    double up = sqrt(anchorUp*anchorUp     + noteB[i]*noteB[i] + noteC[i]*noteC[i]);
    double dn = sqrt(anchorDown*anchorDown + noteB[i]*noteB[i] + noteC[i]*noteC[i]);
    gTotal->SetPoint(i, notePt[i], 0);
    gTotal->SetPointError(i, 0, 0, dn, up);
    gAup->SetPoint(i, notePt[i],  anchorUp); gAdn->SetPoint(i, notePt[i], -anchorDown);
    gBup->SetPoint(i, notePt[i],  noteB[i]); gBdn->SetPoint(i, notePt[i], -noteB[i]);
    gCup->SetPoint(i, notePt[i],  noteC[i]); gCdn->SetPoint(i, notePt[i], -noteC[i]);
    gTup->SetPoint(i, notePt[i],  up);       gTdn->SetPoint(i, notePt[i], -dn);
    cout << Form("%8.1f %8.3f %8.3f %8.3f %8.3f %8.3f/%.3f", notePt[i], anchorUp, anchorDown, noteB[i], noteC[i], up, dn) << endl;
  }

  TCanvas * c = new TCanvas("c", "", 800, 600);
  c->SetLeftMargin(0.13);
  c->SetRightMargin(0.04);
  c->SetBottomMargin(0.13);
  c->SetTopMargin(0.05);
  c->SetTicks(1,1);
  TH1F * frame = c->DrawFrame(5, -5.5, 100, 6.5);
  frame->SetTitle(";Jet #it{p}_{T} [GeV];Data-to-MC JES uncertainty [%]");
  frame->GetXaxis()->SetTitleSize(0.05);
  frame->GetYaxis()->SetTitleSize(0.05);
  frame->GetXaxis()->SetTitleOffset(1.1);
  frame->GetYaxis()->SetTitleOffset(1.1);

  gTotal->SetFillStyle(1001);
  gTotal->SetFillColor(kAzure-9);
  gTotal->SetLineColor(kAzure-9);
  gTotal->Draw("3 same");

  // In-situ window where B and C are anchored; hatched so the band stays visible (batch-mode
  // output drops alpha fills).
  TBox * window = new TBox(15, -5.5, 25, 6.5);
  window->SetFillStyle(3004);
  window->SetFillColor(kGray+1);
  window->Draw("same");

  TLine * zero = new TLine(5, 0, 100, 0);
  zero->SetLineStyle(2);
  zero->SetLineColor(kGray+2);
  zero->Draw();

  auto style = [](TGraph * g, int color, int lstyle, int width) {
    g->SetLineColor(color); g->SetLineStyle(lstyle); g->SetLineWidth(width);
  };
  style(gTup, kBlack,     1, 3); style(gTdn, kBlack,     1, 3);
  style(gAup, kBlue+1,    2, 2); style(gAdn, kBlue+1,    2, 2);
  style(gBup, kGreen+2,   1, 2); style(gBdn, kGreen+2,   1, 2);
  style(gCup, kRed+1,     7, 2); style(gCdn, kRed+1,     7, 2);
  for (TGraph * g : {gAup, gAdn, gBup, gBdn, gCup, gCdn, gTup, gTdn}) g->Draw("L same");

  TLegend * leg = new TLegend(0.50, 0.72, 0.94, 0.935);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(43);
  leg->SetTextSize(15);
  leg->AddEntry(gTup, "Total (quadrature)", "l");
  leg->AddEntry(gAup, Form("A: in-situ #gamma+jet, ^{+%.2f}_{#minus%.2f}%%", anchorUp, anchorDown), "l");
  leg->AddEntry(gBup, "B: T-1044 single hadron (#pm1% tilt)", "l");
  leg->AddEntry(gCup, "C: #pm50% sPHENIX leakage", "l");
  leg->AddEntry(window, "In-situ anchor window", "f");
  leg->Draw();

  insitu_utility::drawSPhenixLabel({"p+p Run24 Data"},
      {Form("Jet R=%.1f, |#eta^{jet}| < 0.7", ana::JetRs[ir]), "B, C: JES leakage note Fig. 23"},
      .17, .87, 16, c->GetWh());

  c->RedrawAxis();
  string pdfPath = Form("%s/draw_jes_highpt_uncertainty.pdf", insitu_pdf_dir);
  c->SaveAs(pdfPath.c_str());

  string outPath = Form("%s/draw_jes_highpt_uncertainty.root", insitu_output_dir);
  TFile * fout = TFile::Open(outPath.c_str(), "RECREATE");
  gTotal->Write();
  fout->Close();
  cout << "Wrote " << pdfPath << " and " << outPath << endl;
}
