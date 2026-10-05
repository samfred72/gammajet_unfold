#include "../src/ana.h"
#include "../src/insitu_utility.h"
#include <string>
#include <vector>
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TBox.h"
#include "TLine.h"
#include "TLegend.h"
#include "TStyle.h"
using namespace std;

R__LOAD_LIBRARY(libgammajet_unfold.so);

// In-situ JES scale factor p_a under every systematic variation, per jet radius.
//
// Reads grid_insitu.C's purity-corrected best fit (results tree: pa_puritycorrected,
// errLow/errHigh_puritycorrected) from insitu/output/grid_insitu_<systag>.root for every
// ana::systags entry - the same numbers draw_jes_summary.C writes into ana.h's
// jesBySystag table, which unfolder.cc uses to correct Data in each variation (each
// variation is propagated with its own p_a - see ana.h). This plot shows those values
// directly:
//   page 1  one pad per radius: p_a for each systag (point = fit, bar = its own
//           statistical Delta chi2 = 1 error), the nominal value +- its statistical error
//           (grey band), and the lower edge of the scan range (red dashed). A value within
//           edgeTol of the scan edge is drawn open red: the fit's minimum is at or beyond
//           the range, so the number is a bound, not a measurement.
//   page 2  fractional shift (p_a_systag - p_a_nominal)/p_a_nominal vs jet radius, one
//           series per variation (high/low pairs share a colour: filled up-triangle =
//           high/first, open down-triangle = low/second), with the nominal statistical
//           uncertainty as the grey band. Points where the variation or the nominal fit is
//           at the scan edge are ringed (their shift is not meaningful). jes_high/jes_low
//           are omitted here (see below).
// Also prints the table. Run after insitu/run_grid.sh (reads only; runs no scan).

const char * insitu_output_dir = ana::path("insitu/output");
const char * insitu_pdf_dir    = ana::path("insitu/pdfs");
const double edgeTol = 5e-4; // |p_a - scanLow| below this counts as "at the scan edge"

namespace {
  struct Pa { bool ok = false; float pa = 0, lo = 0, hi = 0; };
  Pa readPa(const string & systag, int ir) {
    Pa r;
    string fname = string(insitu_output_dir) + "/grid_insitu_" + systag + ".root";
    TFile * f = TFile::Open(fname.c_str(), "READ");
    if (!f || f->IsZombie()) { cout << "WARNING: could not open " << fname << endl; return r; }
    TTree * t = (TTree*)f->Get(Form("%s/results", ana::rnames[ir]));
    if (!t) { f->Close(); return r; }
    float pa, lo, hi;
    t->SetBranchAddress("pa_puritycorrected", &pa);
    t->SetBranchAddress("errLow_puritycorrected", &lo);
    t->SetBranchAddress("errHigh_puritycorrected", &hi);
    t->GetEntry(0);
    f->Close();
    r.ok = true; r.pa = pa; r.lo = lo; r.hi = hi;
    return r;
  }
  bool atEdge(float pa) { return fabs(pa - insitu_utility::scanLow) < edgeTol; }
}

void draw_jes_variations() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  const int nR = ana::nJetR;
  const vector<string> & tags = ana::systags;
  const int nT = tags.size();

  // p[tag][ir]
  map<string, vector<Pa>> p;
  for (const string & t : tags) { p[t] = vector<Pa>(nR); for (int ir = 0; ir < nR; ir++) p[t][ir] = readPa(t, ir); }

  // ---- table ----
  printf("\nIn-situ p_a (purity-corrected) per systag and radius; '*' = at the scan edge (%.3f)\n%-14s", insitu_utility::scanLow, "systag");
  for (int ir = 0; ir < nR; ir++) printf("   R=%.1f  ", ana::JetRs[ir]);
  printf("\n");
  for (const string & t : tags) {
    printf("%-14s", t.c_str());
    for (int ir = 0; ir < nR; ir++) {
      const Pa & a = p[t][ir];
      if (!a.ok) { printf("   ---    "); continue; }
      printf("  %.4f%s ", a.pa, atEdge(a.pa) ? "*" : " ");
    }
    printf("\n");
  }

  // colours: pairs share a colour; singles get their own
  map<string,int> col, mk;
  const int palette[] = {kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+1, kCyan+2, kViolet+1, kPink+6, kAzure+7, kSpring-6, kGray+2};
  int ic = 0;
  for (auto & pr : ana::asymmetricSystagPairs) { col[pr.first] = col[pr.second] = palette[ic++ % 11]; mk[pr.first] = 22; mk[pr.second] = 26; }
  for (const string & t : tags) if (t != "nominal" && !col.count(t)) { col[t] = palette[ic++ % 11]; mk[t] = 21; }
  col["nominal"] = kBlack; mk["nominal"] = 20;

  string pdf = string(insitu_pdf_dir) + "/draw_jes_variations.pdf";
  TCanvas * c = new TCanvas("c", "", 1600, 900);
  c->SaveAs((pdf + "[").c_str());

  // ---- page 1: p_a per systag, one pad per radius ----
  c->Divide(4, 2, 0.002, 0.002);
  for (int ir = 0; ir < nR; ir++) {
    c->cd(ir+1);
    gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.06); gPad->SetBottomMargin(0.3);
    gPad->SetTicks(1, 1);
    const Pa & nom = p["nominal"][ir];
    double ymin = insitu_utility::scanLow - 0.01, ymax = insitu_utility::scanLow + 0.07;
    for (const string & t : tags) if (p[t][ir].ok) ymax = std::max(ymax, (double)(p[t][ir].pa + p[t][ir].hi + 0.01));
    TH1D * fr = new TH1D(Form("fr%d", ir), ";;in-situ p_{a}", nT, 0, nT);
    for (int i = 0; i < nT; i++) fr->GetXaxis()->SetBinLabel(i+1, tags[i].c_str());
    fr->GetXaxis()->LabelsOption("v"); fr->GetXaxis()->SetLabelSize(0.065);
    fr->GetYaxis()->SetLabelSize(0.05); fr->GetYaxis()->SetTitleSize(0.06); fr->GetYaxis()->SetTitleOffset(1.25);
    fr->GetYaxis()->SetNdivisions(505);
    fr->SetMinimum(ymin); fr->SetMaximum(ymax);
    fr->Draw("axis");
    if (nom.ok) {
      TBox * band = new TBox(0, nom.pa - nom.lo, nT, nom.pa + nom.hi);
      band->SetFillColor(kGray); band->SetFillStyle(1001); band->Draw();
      TLine * ln = new TLine(0, nom.pa, nT, nom.pa); ln->SetLineColor(kGray+2); ln->Draw();
    }
    TLine * edge = new TLine(0, insitu_utility::scanLow, nT, insitu_utility::scanLow);
    edge->SetLineColor(kRed+1); edge->SetLineStyle(2); edge->SetLineWidth(2); edge->Draw();
    fr->Draw("axis same");
    for (int i = 0; i < nT; i++) {
      const Pa & a = p[tags[i]][ir];
      if (!a.ok) continue;
      TGraphAsymmErrors * g = new TGraphAsymmErrors(1);
      g->SetPoint(0, i + 0.5, a.pa); g->SetPointError(0, 0, 0, a.lo, a.hi);
      bool e = atEdge(a.pa);
      g->SetMarkerStyle(e ? 24 : 20); g->SetMarkerColor(e ? kRed+1 : kBlack); g->SetLineColor(e ? kRed+1 : kBlack); g->SetMarkerSize(1.0);
      g->Draw("p same");
    }
    // label size scaled to this pad's own height, not the whole canvas'
    insitu_utility::drawSPhenixLabel({}, {Form("Jet R=%.1f", ana::JetRs[ir])}, .2, .86, 16, gPad->GetWh()*gPad->GetHNDC());
  }
  c->cd(8);
  {
    TLegend * l = new TLegend(0.03, 0.3, 0.98, 0.72); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.055);
    TGraph * gf = new TGraph(); gf->SetMarkerStyle(20);
    TGraph * ge = new TGraph(); ge->SetMarkerStyle(24); ge->SetMarkerColor(kRed+1); ge->SetLineColor(kRed+1);
    TBox * bb = new TBox(); bb->SetFillColor(kGray); bb->SetFillStyle(1001);
    TLine * le = new TLine(); le->SetLineColor(kRed+1); le->SetLineStyle(2); le->SetLineWidth(2);
    l->AddEntry(gf, "p_{a} under the variation (stat.)", "p");
    l->AddEntry(ge, "at the scan edge (bound)", "p");
    l->AddEntry(bb, "nominal #pm stat.", "f");
    l->AddEntry(le, Form("scan lower edge (%.2f)", insitu_utility::scanLow), "l");
    l->Draw();
    insitu_utility::drawSPhenixLabel({"p+p Run24 Data vs Pythia8 #gamma+jet"}, {"purity-corrected mean-x_{J} fit"}, .05, .9, 16, gPad->GetWh()*gPad->GetHNDC());
  }
  c->SaveAs(pdf.c_str());

  // ---- page 2: fractional shift vs R ----
  c->Clear();
  c->cd();
  gPad->SetLeftMargin(0.1); gPad->SetRightMargin(0.26); gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.11); gPad->SetTicks(1, 1);
  TH1D * fr2 = new TH1D("frShift", ";jet radius R;(p_{a}^{var} - p_{a}^{nom}) / p_{a}^{nom}", 1, 0.1, 0.9);
  double smin = -0.01, smax = 0.01;
  for (const string & t : tags) for (int ir = 0; ir < nR; ir++) if (p[t][ir].ok && p["nominal"][ir].ok) {
    double s = (p[t][ir].pa - p["nominal"][ir].pa)/p["nominal"][ir].pa; smin = std::min(smin, s); smax = std::max(smax, s);
  }
  fr2->SetMinimum(smin - 0.2*(smax-smin)); fr2->SetMaximum(smax + 0.35*(smax-smin));
  fr2->GetYaxis()->SetTitleOffset(1.1); fr2->GetXaxis()->SetTitleSize(0.045); fr2->GetYaxis()->SetTitleSize(0.045); fr2->GetYaxis()->CenterTitle();
  fr2->Draw("axis");
  for (int ir = 0; ir < nR; ir++) {
    const Pa & nom = p["nominal"][ir];
    if (!nom.ok) continue;
    TBox * b = new TBox(ana::JetRs[ir] - 0.045, -nom.lo/nom.pa, ana::JetRs[ir] + 0.045, nom.hi/nom.pa);
    b->SetFillColor(kGray); b->SetFillStyle(1001); b->Draw();
  }
  TLine * zero = new TLine(0.1, 0, 0.9, 0); zero->SetLineStyle(2); zero->Draw();
  fr2->Draw("axis same");
  TLegend * l2 = new TLegend(0.75, 0.12, 0.99, 0.95); l2->SetBorderSize(0); l2->SetFillStyle(0); l2->SetTextSize(0.028);
  int k = 0;
  // jes_high/jes_low are left out: they vary the JES correction applied in the unfolding,
  // not any input of the in-situ scan (which uses the uncorrected jet pT), so their scan
  // results differ from nominal only by scan noise - and unfolder.cc does not use them.
  TGraph * gEdge = new TGraph(); // ring around every point where var or nominal is at the scan edge
  for (const string & t : tags) {
    if (t == "nominal" || t == "jes_high" || t == "jes_low") continue;
    TGraph * g = new TGraph();
    for (int ir = 0; ir < nR; ir++) if (p[t][ir].ok && p["nominal"][ir].ok) {
      double x = ana::JetRs[ir] + (k - nT/2.0)*0.004, y = (p[t][ir].pa - p["nominal"][ir].pa)/p["nominal"][ir].pa;
      g->SetPoint(g->GetN(), x, y);
      if (atEdge(p[t][ir].pa) || atEdge(p["nominal"][ir].pa)) gEdge->SetPoint(gEdge->GetN(), x, y);
    }
    g->SetMarkerStyle(mk[t]); g->SetMarkerColor(col[t]); g->SetLineColor(col[t]); g->SetMarkerSize(1.1);
    g->Draw("p same");
    l2->AddEntry(g, t.c_str(), "p");
    k++;
  }
  gEdge->SetMarkerStyle(24); gEdge->SetMarkerSize(2.0); gEdge->SetMarkerColor(kBlack);
  if (gEdge->GetN()) { gEdge->Draw("p same"); l2->AddEntry(gEdge, "at scan edge (var. or nom.)", "p"); }
  TBox * bl = new TBox(); bl->SetFillColor(kGray); bl->SetFillStyle(1001);
  l2->AddEntry(bl, "nominal stat. unc.", "f");
  l2->Draw();
  insitu_utility::drawSPhenixLabel({"p+p Run24 Data vs Pythia8 #gamma+jet"},
      {"in-situ p_{a} shift under each variation", "each is propagated with its own p_{a}"}, .13, .9, 18, gPad->GetWh());
  c->SaveAs(pdf.c_str());
  c->SaveAs((pdf + "]").c_str());
  cout << "Wrote " << pdf << endl;
}
