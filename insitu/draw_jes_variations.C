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
#include "TArrow.h"
#include "TLatex.h"
#include "TStyle.h"
using namespace std;

R__LOAD_LIBRARY(libgammajet_unfold.so);

// In-situ p_a under every variation, per jet radius, from grid_insitu.C's purity-corrected fit
// (the values in ana.h's jesBySystag). Run after insitu/run_grid.sh.
//   page 1: one pad per radius, p_a per systag with the nominal +- stat band; values within
//           edgeTol of the scan edge are open red (a bound, not a measurement).
//   page 2: fractional shift vs R per variation (high filled up-triangle, low open down-triangle);
//           edge points are ringed.

const char * insitu_output_dir = ana::path("insitu/output");
const char * insitu_pdf_dir    = ana::path("insitu/pdfs");
const double edgeTol = 5e-4; // |p_a - scanLow| below this counts as "at the scan edge"
const double page1Ymin = 0.88, page1Ymax = 0.955; // page 1 p_a axis range
// The scan-edge line is drawn only if scanLow is inside page 1's range (a TLine is not clipped).
const bool edgeInRange = insitu_utility::scanLow > page1Ymin && insitu_utility::scanLow < page1Ymax;

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
  // Total (stat (+) syst, sign-split as draw_jes_summary.C) per radius, from draw_jes_summary.root.
  struct Total { bool ok = false; float lo = 0, hi = 0; };
  vector<Total> readTotals() {
    vector<Total> tot(ana::nJetR);
    string fname = string(insitu_output_dir) + "/draw_jes_summary.root";
    TFile * f = TFile::Open(fname.c_str(), "READ");
    TTree * t = (f && !f->IsZombie()) ? (TTree*)f->Get("jes_summary") : nullptr;
    if (!t) { cout << "WARNING: no total uncertainty (run draw_jes_summary.C first): " << fname << endl; if (f) f->Close(); return tot; }
    int ir; float lo, hi;
    t->SetBranchAddress("ir", &ir); t->SetBranchAddress("totalLow", &lo); t->SetBranchAddress("totalHigh", &hi);
    for (Long64_t e = 0; e < t->GetEntries(); e++) { t->GetEntry(e); if (ir >= 0 && ir < ana::nJetR) tot[ir] = {true, lo, hi}; }
    f->Close();
    return tot;
  }
  const int totalLineStyle = 7, totalLineWidth = 2;
}

void draw_jes_variations() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  const int nR = ana::nJetR;
  // Herwig (grid_insitu.C("herwig")) is not shown: the current Herwig production has wrong settings.
  const vector<string> & tags = ana::systags;
  const int nT = tags.size();

  // p[tag][ir]
  map<string, vector<Pa>> p;
  const vector<Total> tot = readTotals();
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
  bool anyOutside = false;
  for (int ir = 0; ir < nR; ir++) {
    c->cd(ir+1);
    gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.06); gPad->SetBottomMargin(0.3);
    gPad->SetTicks(1, 1);
    const Pa & nom = p["nominal"][ir];
    const double ymin = page1Ymin, ymax = page1Ymax;
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
      // nominal +- total; a TLine is not clipped, so skip an edge outside the range
      if (tot[ir].ok) for (double y : {nom.pa - tot[ir].lo, nom.pa + tot[ir].hi}) {
        if (y < ymin || y > ymax) { cout << "WARNING: R=" << ana::JetRs[ir] << " total edge " << y << " outside page 1's range" << endl; continue; }
        TLine * lt = new TLine(0, y, nT, y); lt->SetLineStyle(totalLineStyle); lt->SetLineWidth(totalLineWidth); lt->Draw();
      }
    }
    if (edgeInRange) {
      TLine * edge = new TLine(0, insitu_utility::scanLow, nT, insitu_utility::scanLow);
      edge->SetLineColor(kRed+1); edge->SetLineStyle(2); edge->SetLineWidth(2); edge->Draw();
    }
    fr->Draw("axis same");
    for (int i = 0; i < nT; i++) {
      const Pa & a = p[tags[i]][ir];
      if (!a.ok) continue;
      // Outside the zoomed range: an arrow at that edge, labelled with the value.
      if (a.pa < ymin || a.pa > ymax) {
        anyOutside = true;
        const bool up = a.pa > ymax;
        const double tip = up ? ymax - 0.001 : ymin + 0.001, tail = up ? ymax - 0.012 : ymin + 0.012;
        TArrow * ar = new TArrow(i + 0.5, tail, i + 0.5, tip, 0.012, "|>");
        ar->SetLineColor(kBlack); ar->SetFillColor(kBlack); ar->SetLineWidth(2); ar->Draw();
        TLatex * tv = new TLatex(i + 0.5, up ? tail - 0.001 : tail + 0.001, Form("%.3f", a.pa));
        tv->SetTextAngle(90); tv->SetTextAlign(up ? 32 : 12); tv->SetTextSize(0.045); tv->Draw();
        continue;
      }
      TGraphAsymmErrors * g = new TGraphAsymmErrors(1);
      g->SetPoint(0, i + 0.5, a.pa); g->SetPointError(0, 0, 0, a.lo, a.hi);
      bool e = atEdge(a.pa);
      g->SetMarkerStyle(e ? 24 : 20); g->SetMarkerColor(e ? kRed+1 : kBlack); g->SetLineColor(e ? kRed+1 : kBlack); g->SetMarkerSize(1.0);
      g->Draw("p same");
    }
    // label size scaled to this pad's own height, not the whole canvas'
    // top-left up to R=0.5; bottom from R=0.6, where the points sit high in the zoomed range
    insitu_utility::drawSPhenixLabel({}, {Form("Jet R=%.1f", ana::JetRs[ir])}, .2, ana::JetRs[ir] < 0.55 ? .88 : .39, 16, gPad->GetWh()*gPad->GetHNDC());
  }
  c->cd(8);
  {
    TLegend * l = new TLegend(0.03, 0.3, 0.98, 0.72); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.055);
    TGraph * gf = new TGraph(); gf->SetMarkerStyle(20);
    TGraph * ge = new TGraph(); ge->SetMarkerStyle(24); ge->SetMarkerColor(kRed+1); ge->SetLineColor(kRed+1);
    TBox * bb = new TBox(); bb->SetFillColor(kGray); bb->SetFillStyle(1001);
    TLine * le = new TLine(); le->SetLineColor(kRed+1); le->SetLineStyle(2); le->SetLineWidth(2);
    l->AddEntry(gf, "p_{a} under the variation (stat.)", "p");
    if (edgeInRange) l->AddEntry(ge, "at the scan edge (bound)", "p");
    l->AddEntry(bb, "nominal #pm stat.", "f");
    TLine * lt = new TLine(); lt->SetLineStyle(totalLineStyle); lt->SetLineWidth(totalLineWidth);
    l->AddEntry(lt, "nominal #pm total (stat. #oplus syst.)", "l");
    if (edgeInRange) l->AddEntry(le, Form("scan lower edge (%.2f)", insitu_utility::scanLow), "l");
    TArrow * la = new TArrow(); la->SetLineWidth(2);
    if (anyOutside) l->AddEntry(la, Form("outside %.2f-%.2f (value printed)", page1Ymin, page1Ymax), "l");
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
  for (int ir = 0; ir < nR; ir++) if (tot[ir].ok && p["nominal"][ir].ok) {
    smin = std::min(smin, -(double)tot[ir].lo/p["nominal"][ir].pa); smax = std::max(smax, (double)tot[ir].hi/p["nominal"][ir].pa);
  }
  fr2->SetMinimum(smin - 0.2*(smax-smin)); fr2->SetMaximum(smax + 0.35*(smax-smin));
  fr2->GetYaxis()->SetTitleOffset(1.1); fr2->GetXaxis()->SetTitleSize(0.045); fr2->GetYaxis()->SetTitleSize(0.045); fr2->GetYaxis()->CenterTitle();
  fr2->Draw("axis");
  for (int ir = 0; ir < nR; ir++) {
    const Pa & nom = p["nominal"][ir];
    if (!nom.ok) continue;
    TBox * b = new TBox(ana::JetRs[ir] - 0.045, -nom.lo/nom.pa, ana::JetRs[ir] + 0.045, nom.hi/nom.pa);
    b->SetFillColor(kGray); b->SetFillStyle(1001); b->Draw();
    if (tot[ir].ok) for (double y : {-tot[ir].lo/nom.pa, tot[ir].hi/nom.pa}) {
      TLine * lt = new TLine(ana::JetRs[ir] - 0.045, y, ana::JetRs[ir] + 0.045, y);
      lt->SetLineStyle(totalLineStyle); lt->SetLineWidth(totalLineWidth); lt->Draw();
    }
  }
  TLine * zero = new TLine(0.1, 0, 0.9, 0); zero->SetLineStyle(2); zero->Draw();
  fr2->Draw("axis same");
  TLegend * l2 = new TLegend(0.75, 0.12, 0.99, 0.95); l2->SetBorderSize(0); l2->SetFillStyle(0); l2->SetTextSize(0.028);
  int k = 0;
  // jes_high/low vary the correction applied in unfolding, not the scan input: omitted.
  TGraph * gEdge = new TGraph(); // ring points at the scan edge
  for (const string & t : tags) {
    if (t == "nominal" || t == "jes_high" || t == "jes_low") continue;
    TGraph * g = new TGraph();
    for (int ir = 0; ir < nR; ir++) if (p[t][ir].ok && p["nominal"][ir].ok) {
      double x = ana::JetRs[ir] + (k - nT/2.0)*0.004, y = (p[t][ir].pa - p["nominal"][ir].pa)/p["nominal"][ir].pa;
      g->SetPoint(g->GetN(), x, y);
      if (atEdge(p[t][ir].pa) || atEdge(p["nominal"][ir].pa)) gEdge->SetPoint(gEdge->GetN(), x, y);
    }
    if (g->GetN() == 0) continue; // no scan output for this systag yet
    g->SetMarkerStyle(mk[t]); g->SetMarkerColor(col[t]); g->SetLineColor(col[t]); g->SetMarkerSize(1.1);
    g->Draw("p same");
    l2->AddEntry(g, t.c_str(), "p");
    k++;
  }
  gEdge->SetMarkerStyle(24); gEdge->SetMarkerSize(2.0); gEdge->SetMarkerColor(kBlack);
  if (gEdge->GetN()) { gEdge->Draw("p same"); l2->AddEntry(gEdge, "at scan edge (var. or nom.)", "p"); }
  TBox * bl = new TBox(); bl->SetFillColor(kGray); bl->SetFillStyle(1001);
  l2->AddEntry(bl, "nominal stat. unc.", "f");
  TLine * lt2 = new TLine(); lt2->SetLineStyle(totalLineStyle); lt2->SetLineWidth(totalLineWidth);
  l2->AddEntry(lt2, "total unc. (stat. #oplus syst.)", "l");
  l2->Draw();
  insitu_utility::drawSPhenixLabel({"p+p Run24 Data vs Pythia8 #gamma+jet"},
      {"in-situ p_{a} shift under each variation", "each is propagated with its own p_{a}"}, .13, .9, 18, gPad->GetWh());
  c->SaveAs(pdf.c_str());
  c->SaveAs((pdf + "]").c_str());
  cout << "Wrote " << pdf << endl;
}
