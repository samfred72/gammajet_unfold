#include "../src/ana.h"
#include "../src/insitu_utility.h"
#include <string>
#include "TFile.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TBox.h"
#include "TH1F.h"
using namespace std;

R__LOAD_LIBRARY(libgammajet_unfold.so);

// pT-dependent data-to-MC JES uncertainty for R=0.4 jets (in-situ anchor A, T-1044 single hadron B,
// sPHENIX leakage C, and their quadrature total), from insitu/compute_jes_highpt_uncertainty.C's
// output/jes_highpt_uncertainty.root. Drawing only.
// Output: insitu/pdfs/draw_jes_highpt_uncertainty.pdf

void draw_jes_highpt_uncertainty() {
  const int ir = 2; // R=0.4, the only radius the curves exist for
  const char * insitu_output_dir = ana::path("insitu/output");
  const char * insitu_pdf_dir    = ana::path("insitu/pdfs");
  string inPath = Form("%s/jes_highpt_uncertainty.root", insitu_output_dir);
  TFile * fin = TFile::Open(inPath.c_str());
  if (!fin || fin->IsZombie()) {
    cout << "ERROR: cannot open " << inPath << " - run insitu/compute_jes_highpt_uncertainty.C first." << endl;
    return;
  }
  TGraphAsymmErrors * gTotal = (TGraphAsymmErrors*)fin->Get("gJESUnc_total");
  TGraph * gAup = (TGraph*)fin->Get("gA_up"), * gAdn = (TGraph*)fin->Get("gA_down");
  TGraph * gBup = (TGraph*)fin->Get("gB_up"), * gBdn = (TGraph*)fin->Get("gB_down");
  TGraph * gCup = (TGraph*)fin->Get("gC_up"), * gCdn = (TGraph*)fin->Get("gC_down");
  TGraph * gTup = (TGraph*)fin->Get("gTotal_up"), * gTdn = (TGraph*)fin->Get("gTotal_down");
  TVectorD * anchor = (TVectorD*)fin->Get("anchor_percent");
  const double anchorUp = (*anchor)[0], anchorDown = (*anchor)[1];

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

  insitu_utility::drawSPhenixLabel({ana::dataLabel()},
      {Form("Jet R=%.1f, |#eta^{jet}| < 0.7", ana::JetRs[ir]), "B, C: JES leakage note Fig. 23"},
      .17, .87, 16, c->GetWh());

  c->RedrawAxis();
  string pdfPath = Form("%s/draw_jes_highpt_uncertainty.pdf", insitu_pdf_dir);
  c->SaveAs(pdfPath.c_str());
  cout << "Wrote " << pdfPath << endl;
}
