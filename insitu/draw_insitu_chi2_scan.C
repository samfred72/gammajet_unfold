#include "../src/ana.h"
#include "../src/insitu_utility.h"
#include <string>
#include <cmath>
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TStyle.h"

R__LOAD_LIBRARY(libgammajet_unfold.so);

// chi2 vs p_a for the purity-corrected mean in-situ fit (the one behind ana::jesNominal), read
// from grid_insitu.C's gchisq_puritycorrected. Best fit and errors as grid_insitu.C (chi2_min + 1
// crossings on the full graph); the band is the Delta chi2 < 1 interval.
const char * scan_output_dir = ana::path("insitu/output");
const char * scan_pdf_dir    = ana::path("insitu/pdfs");

void draw_insitu_chi2_scan(string systag = "nominal", int ir = 2) {
  gStyle->SetOptStat(0);

  const char * filename = Form("%s/grid_insitu_%s.root", scan_output_dir, systag.c_str());
  TFile * f = TFile::Open(filename);
  if (!f || f->IsZombie()) {
    cout << "Missing " << filename << " - run grid_insitu.C(\"" << systag << "\") first." << endl;
    return;
  }
  TGraph * g = (TGraph*)f->Get(Form("%s/gchisq_puritycorrected", ana::rnames[ir]));
  if (!g) {
    cout << "No " << ana::rnames[ir] << "/gchisq_puritycorrected in " << filename << endl;
    return;
  }

  // Minimum over the full-resolution scan, then the chi2_min+1 crossings.
  int ibest = 0;
  double xbest, minchisq;
  g->GetPoint(0, xbest, minchisq);
  for (int i = 1; i < g->GetN(); i++) {
    double x, y;
    g->GetPoint(i, x, y);
    if (std::isfinite(y) && y < minchisq) { minchisq = y; xbest = x; ibest = i; }
  }
  float errLow, errHigh;
  insitu_utility::findError(g, ibest, minchisq, errLow, errHigh);

  // Display window +-halfWindow around the minimum, so the chi2_min+1 line is legible.
  const double halfWindow = 0.025;
  double xlow  = std::max((double)insitu_utility::scanLow,  xbest - halfWindow);
  double xhigh = std::min((double)insitu_utility::scanHigh, xbest + halfWindow);
  double ymax = 0;
  for (int i = 0; i < g->GetN(); i++) {
    double x, y;
    g->GetPoint(i, x, y);
    if (x < xlow || x > xhigh) continue;
    if (std::isfinite(y) && y > ymax) ymax = y;
  }
  if (ymax <= 0) ymax = 1;
  ymax *= 1.7;

  TCanvas * c = new TCanvas("c", "", 900, 700);
  c->SetLeftMargin(.14);
  c->SetRightMargin(.04);
  c->SetTopMargin(.05);
  c->SetBottomMargin(.13);
  c->SetTicks(1, 1);

  TH1F * frame = c->DrawFrame(xlow, 0, xhigh, ymax);
  frame->SetTitle(";Data-to-MC JES correction p_{a};#chi^{2} (mean x_{J#gamma}, Data/MC vs. 1)");
  frame->GetXaxis()->SetTitleSize(.05);
  frame->GetYaxis()->SetTitleSize(.05);
  frame->GetXaxis()->SetLabelSize(.042);
  frame->GetYaxis()->SetLabelSize(.042);
  frame->GetYaxis()->SetTitleOffset(1.25);

  // Band first, below 55% of the frame height (clear of the labels).
  TBox * band = new TBox(xbest - errLow, 0, xbest + errHigh, 0.55*ymax);
  // Opaque fill (alpha fills do not render in PDF); RedrawAxis() restores the ticks.
  band->SetFillColor(kRed-10);
  band->SetFillStyle(1001);
  band->SetLineWidth(0);
  band->Draw();

  // Display-only thinning (as draw_grid_chi2.C).
  const int drawStride = 5;
  TGraph * gDisplay = new TGraph();
  for (int i = 0; i < g->GetN(); i += drawStride) {
    double x, y;
    g->GetPoint(i, x, y);
    if (x < xlow || x > xhigh) continue;
    if (std::isfinite(y)) gDisplay->SetPoint(gDisplay->GetN(), x, y);
  }
  gDisplay->SetMarkerColor(kBlack);
  gDisplay->SetMarkerStyle(kFullCircle);
  gDisplay->SetMarkerSize(.8);
  gDisplay->Draw("P SAME");

  TLine * vline = new TLine(xbest, 0, xbest, ymax);
  vline->SetLineColor(kGreen+2);
  vline->SetLineStyle(2);
  vline->SetLineWidth(2);
  vline->Draw();

  TLine * hline = new TLine(xlow, minchisq+1, xhigh, minchisq+1);
  hline->SetLineColor(kRed+1);
  hline->SetLineStyle(3);
  hline->SetLineWidth(2);
  hline->Draw();

  c->RedrawAxis();

  TLegend * leg = new TLegend(.53, .64, .90, .92);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(.034);
  leg->AddEntry(gDisplay, "Purity-corrected scan", "p");
  leg->AddEntry(vline, Form("Minimum: p_{a} = %.4f^{+%.4f}_{-%.4f}", xbest, errHigh, errLow), "l");
  leg->AddEntry(hline, "#chi^{2}_{min} + 1", "l");
  leg->AddEntry(band, "#Delta#chi^{2} < 1", "f");
  leg->Draw();

  insitu_utility::drawSPhenixLabel({"p+p Run24 Data"},
      {"Pythia8 #gamma+jet MC reference", Form("Jet R=%.1f", ana::JetRs[ir])},
      .18, .90, 20, c->GetWh());
  // Label block sits left of the minimum line; the legend sits right of it (see TLegend).

  string pdfPath = Form("%s/insitu_chi2_scan_%s_%s.pdf", scan_pdf_dir, systag.c_str(), ana::rnames[ir]);
  c->SaveAs(pdfPath.c_str());
  cout << "Best-fit p_a = " << xbest << " +" << errHigh << "/-" << errLow
       << " (chi2_min=" << minchisq << ")" << endl;
}
