#include "../src/ana.h"
#include "../src/insitu_utility.h"
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TGraphAsymmErrors.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
using namespace std;

R__LOAD_LIBRARY(libgammajet_unfold.so);

// In-situ JES p_a vs jet radius: statistical error as points, systematic as boxes, from
// insitu/compute_jes_summary.C's output/jes_summary.root. Drawing only.
// Output: insitu/pdfs/draw_jes_summary.pdf

void draw_jes_summary() {
  gStyle->SetOptStat(0);
  const char * insitu_output_dir = ana::path("insitu/output");
  const char * insitu_pdf_dir    = ana::path("insitu/pdfs");

  string infile = Form("%s/jes_summary.root", insitu_output_dir);
  TFile * fin = TFile::Open(infile.c_str());
  if (!fin || fin->IsZombie()) {
    cout << "ERROR: cannot open " << infile << " - run insitu/compute_jes_summary.C first." << endl;
    return;
  }
  TTree * t = (TTree*)fin->Get("jes_summary");
  float R, pa, systLow, systHigh, totalLow, totalHigh;
  t->SetBranchAddress("R", &R);
  t->SetBranchAddress("pa", &pa);
  t->SetBranchAddress("systLow", &systLow);
  t->SetBranchAddress("systHigh", &systHigh);
  t->SetBranchAddress("totalLow", &totalLow);
  t->SetBranchAddress("totalHigh", &totalHigh);
  vector<float> vR, vPa, vSystLow, vSystHigh, vTotalLow, vTotalHigh;
  for (Long64_t e = 0; e < t->GetEntries(); e++) {
    t->GetEntry(e);
    vR.push_back(R); vPa.push_back(pa);
    vSystLow.push_back(systLow); vSystHigh.push_back(systHigh);
    vTotalLow.push_back(totalLow); vTotalHigh.push_back(totalHigh);
  }
  TGraphAsymmErrors * gStat = (TGraphAsymmErrors*)fin->Get("gJES_stat_vs_radius");
  int nOk = (int)vR.size();

  TCanvas * c = new TCanvas("c", "", 800, 600);
  gPad->SetTicks(1,1);

  double ymin = 1e9, ymax = -1e9;
  for (int i = 0; i < nOk; i++) {
    ymin = std::min(ymin, (double)(vPa[i] - vTotalLow[i]));
    ymax = std::max(ymax, (double)(vPa[i] + vTotalHigh[i]));
  }
  double pad = (ymax-ymin)*0.25;
  TH1F * frame = c->DrawFrame(vR.front()-0.05, ymin-pad, vR.back()+0.05, ymax+pad);
  frame->SetTitle(";Jet R;Data-to-MC JES Correction");
  frame->GetXaxis()->SetNdivisions(8);

  // Solid fill: batch-mode output drops alpha fills. SetFillStyle(1001) is required for a TBox.
  const float boxHalfWidth = 0.04;
  vector<TBox*> boxes;
  for (int i = 0; i < nOk; i++) {
    TBox * b = new TBox(vR[i]-boxHalfWidth, vPa[i]-vSystLow[i], vR[i]+boxHalfWidth, vPa[i]+vSystHigh[i]);
    b->SetFillStyle(1001);
    b->SetFillColor(kAzure-9);
    b->SetLineColor(kAzure+2);
    b->Draw("same");
    boxes.push_back(b);
  }

  // Statistical points on top of the boxes.
  gStat->SetMarkerStyle(20);
  gStat->SetMarkerSize(1.2);
  gStat->SetLineWidth(2);
  gStat->SetLineColor(kBlack);
  gStat->SetMarkerColor(kBlack);
  gStat->Draw("PZ SAME");

  TLegend * leg = new TLegend(0.45, 0.72, 0.88, 0.85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(gStat, "Statistical uncertainty", "lep");
  leg->AddEntry(boxes[0], "Systematic uncertainty", "f");
  leg->Draw();

  insitu_utility::drawSPhenixLabel({ana::dataLabel()}, {ana::mcLabel(), "Purity-corrected"}, .18, .85, 16, c->GetWh());

  c->RedrawAxis();
  const char * pdfPath = Form("%s/draw_jes_summary.pdf", insitu_pdf_dir);
  c->SaveAs(pdfPath);
  cout << "Wrote " << pdfPath << endl;
}
