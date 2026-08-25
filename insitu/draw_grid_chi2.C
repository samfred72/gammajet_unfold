#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/insitu_utility.h"
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

// Compares the chi2(-analogue)-vs-p_a scan curves from the five grid_insitu*.C JES-scan
// methods side by side, for one systag - one page per jet radius, all seven in one
// multi-page PDF (this macro no longer takes an `ir` argument; it loops
// ana::nJetR internally, matching every grid_insitu*.C macro it reads). None of
// grid_insitu.C, grid_insitu_shapechi2.C, grid_insitu_unfolded.C,
// grid_insitu_unfolded_shapechi2.C, or grid_insitu_jet12.C actually draws its own
// gchisq_regionA/gchisq_puritycorrected/gchisq_unfolded/gchisq_regionA_jet12ref graph as
// a PDF page - each only Write()s it into its own output .root, one ana::rnames[ir]
// subdirectory per radius (see that macro's "// Save" block at the end) - so this reads
// those already-written graphs back out and draws them together on one page per radius,
// instead of duplicating each macro's own scan loop. Run the grid_insitu*.C macros for
// `systag` first (see run_grid.sh); this macro only reads their output and does no
// scanning itself.
//
// Pads 3 and 6 read grid_insitu_jet12.C's and grid_insitu_jet12_shapechi2.C's output
// (Data Region A vs. Jet12_long Region A, no truth-level jet-pT cut - the legacy
// "Jet12_full" cross-check those two macros used to also support has been retired).
// Since those two macros have to be re-run through the production pipeline
// (runall_unfold_allsys.sh) and then themselves before their output exists for a given
// systag, both are treated as optional here like the other three: if
// grid_insitu_jet12{,_shapechi2}_<systag>.root is missing entirely, pad 3/6 is skipped
// on every page; if a given radius's subdirectory is missing from an otherwise-present
// file (a partial rerun), it's skipped just for that radius's page.
// insitu/ is split into output/ (the grid_insitu*.C macros' .root output, which this
// macro only reads) and pdfs/ (this macro's own .pdf output) - it reads no input
// ntuples of its own.
const char * insitu_output_dir = "/home/samson72/sphnx/gammajet_unfold/insitu/output";
const char * insitu_pdf_dir    = "/home/samson72/sphnx/gammajet_unfold/insitu/pdfs";

// Index and chi2 of a TGraph's minimum point - grid_insitu*.C track this inline during
// their own scan loops; reading the graph back from file after the fact needs its own
// small re-derivation.
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

// One pad: draws 1-2 chi2-vs-p_a graphs together, x-axis fixed to
// [insitu_utility::scanLow, scanHigh] (the shared scan window - see src/insitu_utility.h)
// regardless of what other pads on the same page need. Each curve gets a dashed
// vertical line at its own best-fit p_a and a dotted horizontal line at its own
// minchisq+1 (the 68% CL threshold - same convention as insitu_utility::findError,
// used by every grid_insitu*.C's own FINAL RESULT printout), both color-matched to the
// curve - drawn separately per curve rather than shared, since the two curves on the
// Region-A/Purity-corrected pads don't sit on the same absolute chi2 scale.
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
    // Markers, not a connecting line: these are the na actually-scanned (pa,chi2)
    // points, not an interpolation/fit - a line here (even though TGraph::Draw("L")
    // only connects real points, drawing nothing in between) reads as a smooth fitted
    // curve, which is misleading given how discontinuous the shape-chi2 curves turn out
    // to be (see debug_shapechi2_spike.C).
    //
    // g itself has na=insitu_utility::scanN (1000) points spaced 1e-4 apart in p_a -
    // far denser than this pad's ~450px width (2+ points per pixel), so even drawing
    // real markers for every point renders as an apparently continuous curve rather
    // than visibly discrete points, especially for scans with little point-to-point
    // statistical noise to break that continuity (e.g. the Jet12-referenced pad, which
    // otherwise looked just like a smooth fitted function despite being markers-only).
    // Thin the *display* only, onto its own TGraph at a real, visible marker size -
    // findMin/findError below still run over the full-resolution g, so the reported
    // best-fit p_a/chi2/errors are unaffected by this display-only subsampling.
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

  // Each input is now one file for the whole systag (all seven jet radii, one
  // ana::rnames[ir] subdirectory each - see grid_insitu.C) - opened once, instead of a
  // different file per radius.
  TFile * fMean        = TFile::Open(Form("%s/grid_insitu_%s.root", insitu_output_dir, systag.c_str()));
  TFile * fShape        = TFile::Open(Form("%s/grid_insitu_shapechi2_%s.root", insitu_output_dir, systag.c_str()));
  TFile * fUnfold       = TFile::Open(Form("%s/grid_insitu_unfolded_%s.root", insitu_output_dir, systag.c_str()));
  TFile * fUnfoldShape  = TFile::Open(Form("%s/grid_insitu_unfolded_shapechi2_%s.root", insitu_output_dir, systag.c_str()));

  if (!fMean || fMean->IsZombie() || !fShape || fShape->IsZombie() ||
      !fUnfold || fUnfold->IsZombie() || !fUnfoldShape || fUnfoldShape->IsZombie()) {
    cout << "Missing one of the four grid_insitu*.C output files for systag=" << systag
         << " - run grid_insitu.C(\"" << systag << "\"), grid_insitu_shapechi2.C(...), "
         << "grid_insitu_unfolded.C(\"" << systag << "\"), and grid_insitu_unfolded_shapechi2.C(...) first."
         << endl;
    return;
  }

  // Optional fifth method: grid_insitu_jet12.C's Region-A(Data)-vs-Region-A(Jet12_long
  // MC) scan (see comment at the top of this file for why this stays optional rather
  // than joining the four-file mandatory check above).
  const char * jet12filename = Form("%s/grid_insitu_jet12_%s.root", insitu_output_dir, systag.c_str());
  TFile * fJet12 = TFile::Open(jet12filename);
  if (!fJet12 || fJet12->IsZombie()) {
    cout << "WARNING: missing " << jet12filename << " - run grid_insitu_jet12.C(\"" << systag
         << "\") first if you want the Region-A-vs-Jet12 pad. Skipping pad 3 on every page." << endl;
    fJet12 = nullptr;
  }

  // Shape-chi2 counterpart, from grid_insitu_jet12_shapechi2.C - same optional
  // treatment as fJet12 above.
  const char * jet12ShapeFilename = Form("%s/grid_insitu_jet12_shapechi2_%s.root", insitu_output_dir, systag.c_str());
  TFile * fJet12Shape = TFile::Open(jet12ShapeFilename);
  if (!fJet12Shape || fJet12Shape->IsZombie()) {
    cout << "WARNING: missing " << jet12ShapeFilename << " - run grid_insitu_jet12_shapechi2.C(\"" << systag
         << "\") first if you want the Region-A-vs-Jet12 shape-chi2 pad. Skipping pad 6 on every page." << endl;
    fJet12Shape = nullptr;
  }

  string pdfPathStr = Form("%s/grid_chi2_comparison_%s.pdf", insitu_pdf_dir, systag.c_str());
  TCanvas * c = new TCanvas("c", "", 1500, 900);
  c->SaveAs(Form("%s[", pdfPathStr.c_str()));

  for (int ir = 0; ir < ana::nJetR; ir++) {
    // A given systag's files might have been produced for only some radii so far (a
    // partial rerun) - check the radius subdirectory exists in each already-open file,
    // not just that the file itself isn't a zombie.
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

    // Column-major layout: col 1 = Region A/Purity-corrected (mean chi2 on top, shape
    // chi2 below), col 2 = Unfolded vs. truth (mean chi2 on top, shape chi2 below),
    // col 3 = Region A vs. Jet12_long (mean chi2 on top, shape chi2 below).
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
