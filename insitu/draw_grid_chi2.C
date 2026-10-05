#include "../src/ana.h"
#include "../src/insitu_utility.h"
#include <string>
#include <vector>
#include <cfloat>
#include <cmath>
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Overlays the chi2-vs-p_a curves written by the grid_insitu*.C methods for one systag, one page
// per jet radius. Reads their output only (run them first). The Jet12 inputs (pads 3 and 6) are
// optional; a missing file or radius subdirectory skips that pad.
const char * insitu_output_dir = ana::path("insitu/output");
const char * insitu_pdf_dir    = ana::path("insitu/pdfs");

// Index and chi2 of a TGraph's minimum point.
void findMin(TGraph * g, int & ibest, double & minchisq) {
  ibest = 0;
  double x0, y0;
  g->GetPoint(0, x0, y0);
  minchisq = y0;
  for (int i = 1; i < g->GetN(); i++) {
    double x, y;
    g->GetPoint(i, x, y);
    if (y < minchisq) { minchisq = y; ibest = i; }
  }
}

// One pad: 1-2 chi2 curves on the fixed scan window, each with a dashed line at its best p_a and
// a dotted line at its minchisq+1, color-matched (the curves need not share a chi2 scale).
void drawChi2Pad(TVirtualPad * p, const char * padTitle, const char * yaxisTitle,
    const string & systag, int ir, vector<TGraph*> graphs, vector<string> labels, vector<int> colors) {
  p->SetLeftMargin(.15);
  p->SetBottomMargin(.13);
  p->SetTopMargin(.09);
  p->SetRightMargin(.04);

  double ymax = 0;
  for (TGraph * g : graphs) {
    for (int i = 0; i < g->GetN(); i++) {
      double x, y;
      g->GetPoint(i, x, y);
      if (std::isfinite(y) && y > ymax) ymax = y;
    }
  }
  if (ymax <= 0) ymax = 1;
  ymax *= 1.15;

  TH1F * frame = p->DrawFrame(insitu_utility::scanLow, 0, insitu_utility::scanHigh, ymax);
  frame->SetTitle(Form(";p_{a} (jet_{pt,corrected} = jet_{pt}/p_{a});%s", yaxisTitle));
  frame->GetXaxis()->SetTitleSize(.045);
  frame->GetYaxis()->SetTitleSize(.045);
  frame->GetXaxis()->SetLabelSize(.035);
  frame->GetYaxis()->SetLabelSize(.035);

  TLegend * leg = new TLegend(.24, .48, .91, .70);
  leg->SetBorderSize(0);
  leg->SetTextSize(.034);

  for (size_t k = 0; k < graphs.size(); k++) {
    TGraph * g = graphs[k];
    // Markers, not a line: these are scanned points. The display is thinned so the markers are
    // visible; findMin/findError use the full graph.
    const int drawStride = 10;
    TGraph * gDisplay = new TGraph();
    for (int i = 0; i < g->GetN(); i += drawStride) {
      double xi, yi;
      g->GetPoint(i, xi, yi);
      gDisplay->SetPoint(gDisplay->GetN(), xi, yi);
    }
    gDisplay->SetMarkerColor(colors[k]);
    gDisplay->SetMarkerStyle(kFullCircle);
    gDisplay->SetMarkerSize(.6);
    gDisplay->Draw("P SAME");
    g->SetMarkerColor(colors[k]);
    g->SetMarkerStyle(kFullCircle);
    g->SetMarkerSize(.6);

    int ibest;
    double minchisq;
    findMin(g, ibest, minchisq);
    float errLow, errHigh;
    insitu_utility::findError(g, ibest, minchisq, errLow, errHigh);
    double xbest, ybest;
    g->GetPoint(ibest, xbest, ybest);

    TLine * vline = new TLine(xbest, 0, xbest, ymax);
    vline->SetLineColor(colors[k]);
    vline->SetLineStyle(2);
    vline->Draw();

    TLine * hline = new TLine(insitu_utility::scanLow, minchisq+1, insitu_utility::scanHigh, minchisq+1);
    hline->SetLineColor(colors[k]);
    hline->SetLineStyle(3);
    hline->Draw();

    leg->AddEntry(g, Form("%s: p_{a}=%.4f^{+%.4f}_{-%.4f}", labels[k].c_str(), xbest, errHigh, errLow), "p");
  }
  leg->Draw();

  TLatex title;
  title.SetNDC();
  title.SetTextSize(.045);
  title.SetTextFont(62);
  title.DrawLatex(.15, .93, padTitle);

  insitu_utility::drawSPhenixLabel({"p+p Run24 Data"},
      {Form("Jet R=%.1f", ana::JetRs[ir]), Form("systag=%s", systag.c_str())},
      .40, .85, 14, p->GetWh()*0.55);
}

void draw_grid_chi2(string systag = "nominal") {
  gStyle->SetOptStat(0);

  TFile * fMean        = TFile::Open(Form("%s/grid_insitu_%s.root", insitu_output_dir, systag.c_str()));
  TFile * fShape        = TFile::Open(Form("%s/grid_insitu_shapechi2_%s.root", insitu_output_dir, systag.c_str()));
  TFile * fUnfold       = TFile::Open(Form("%s/grid_insitu_unfolded_%s.root", insitu_output_dir, systag.c_str()));
  TFile * fUnfoldShape  = TFile::Open(Form("%s/grid_insitu_unfolded_shapechi2_%s.root", insitu_output_dir, systag.c_str()));

  if (!fMean || fMean->IsZombie() || !fShape || fShape->IsZombie() ||
      !fUnfold || fUnfold->IsZombie() || !fUnfoldShape || fUnfoldShape->IsZombie()) {
    cout << "Missing one of the four grid_insitu*.C output files for systag=" << systag
         << " - run grid_insitu.C(\"" << systag << "\") with method \"mean\" and \"shape\", "
         << "and grid_insitu_unfolded.C(\"" << systag << "\") with both methods first."
         << endl;
    return;
  }

  // Optional: the Jet12_long-referenced scan.
  const char * jet12filename = Form("%s/grid_insitu_jet12_%s.root", insitu_output_dir, systag.c_str());
  TFile * fJet12 = TFile::Open(jet12filename);
  if (!fJet12 || fJet12->IsZombie()) {
    cout << "WARNING: missing " << jet12filename << " - run grid_insitu_jet12.C(\"" << systag
         << "\") first if you want the Region-A-vs-Jet12 pad. Skipping pad 3 on every page." << endl;
    fJet12 = nullptr;
  }

  // Its shape counterpart, also optional.
  const char * jet12ShapeFilename = Form("%s/grid_insitu_jet12_shapechi2_%s.root", insitu_output_dir, systag.c_str());
  TFile * fJet12Shape = TFile::Open(jet12ShapeFilename);
  if (!fJet12Shape || fJet12Shape->IsZombie()) {
    cout << "WARNING: missing " << jet12ShapeFilename << " - run grid_insitu_jet12.C(\"" << systag << "\",\"shape\") (\"" << systag
         << "\") first if you want the Region-A-vs-Jet12 shape-chi2 pad. Skipping pad 6 on every page." << endl;
    fJet12Shape = nullptr;
  }

  string pdfPathStr = Form("%s/grid_chi2_comparison_%s.pdf", insitu_pdf_dir, systag.c_str());
  TCanvas * c = new TCanvas("c", "", 1500, 900);
  c->SaveAs(Form("%s[", pdfPathStr.c_str()));

  for (int ir = 0; ir < ana::nJetR; ir++) {
    // Skip radii missing from a partial rerun.
    const char * rname = ana::rnames[ir];
    TDirectory * dMean       = fMean->GetDirectory(rname);
    TDirectory * dShape      = fShape->GetDirectory(rname);
    TDirectory * dUnfold     = fUnfold->GetDirectory(rname);
    TDirectory * dUnfoldShape = fUnfoldShape->GetDirectory(rname);
    if (!dMean || !dShape || !dUnfold || !dUnfoldShape) {
      cout << "Missing " << rname << " subdirectory in one of the four mandatory grid_insitu*.C "
           << "files for systag=" << systag << " - skipping this radius's page." << endl;
      continue;
    }

    TGraph * gJet12 = fJet12 ? (TGraph*)fJet12->Get(Form("%s/gchisq_regionA_jet12ref", rname)) : nullptr;
    TGraph * gJet12Shape = fJet12Shape ? (TGraph*)fJet12Shape->Get(Form("%s/gchisq_regionA_jet12ref", rname)) : nullptr;

    TGraph * gMeanA        = (TGraph*)dMean->Get("gchisq_regionA");
    TGraph * gMeanCorr     = (TGraph*)dMean->Get("gchisq_puritycorrected");
    TGraph * gShapeA       = (TGraph*)dShape->Get("gchisq_regionA");
    TGraph * gShapeCorr    = (TGraph*)dShape->Get("gchisq_puritycorrected");
    TGraph * gUnfold       = (TGraph*)dUnfold->Get("gchisq_unfolded");
    TGraph * gUnfoldShape  = (TGraph*)dUnfoldShape->Get("gchisq_unfolded");

    c->Clear();
    c->Divide(3, 2);

    // Columns: region A / purity-corrected, unfolded vs truth, Jet12_long; mean chi2 on top, shape below.
    drawChi2Pad(c->cd(1), "Region A / Purity-corrected", "#chi^{2} (mean x_{J})", systag, ir,
        {gMeanA, gMeanCorr}, {"Region A", "Purity-corrected"}, {kBlue+1, kRed+1});
    drawChi2Pad(c->cd(2), "Unfolded vs. truth", "#chi^{2} (mean x_{J})", systag, ir,
        {gUnfold}, {"Purity-corrected + unfolded"}, {kGreen+2});
    if (gJet12) {
      drawChi2Pad(c->cd(3), "Region A vs. Jet12_long (dijet MC)", "#chi^{2} (mean x_{J})", systag, ir,
          {gJet12}, {"Region A (Jet12_long ref.)"}, {kMagenta+1});
    }
    drawChi2Pad(c->cd(4), "Region A / Purity-corrected", "Shape #chi^{2}", systag, ir,
        {gShapeA, gShapeCorr}, {"Region A", "Purity-corrected"}, {kBlue+1, kRed+1});
    drawChi2Pad(c->cd(5), "Unfolded vs. truth", "Shape #chi^{2}", systag, ir,
        {gUnfoldShape}, {"Purity-corrected + unfolded"}, {kOrange+7});
    if (gJet12Shape) {
      drawChi2Pad(c->cd(6), "Region A vs. Jet12_long (dijet MC)", "Shape #chi^{2}", systag, ir,
          {gJet12Shape}, {"Region A (Jet12_long ref.)"}, {kMagenta+1});
    }

    c->SaveAs(pdfPathStr.c_str());
  }

  c->SaveAs(Form("%s]", pdfPathStr.c_str()));
  cout << "Wrote " << pdfPathStr << endl;
}
