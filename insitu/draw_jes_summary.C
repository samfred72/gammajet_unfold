#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/insitu_utility.h"
#include <string>
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TGraphAsymmErrors.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
using namespace std;

R__LOAD_LIBRARY(libgammajet_unfold.so);

// Final in-situ JES scale factor per jet radius, with full systematic uncertainty -
// reads grid_insitu.C's purity-corrected best-fit p_a ("results" tree,
// pa_puritycorrected/errLow_puritycorrected/errHigh_puritycorrected) for every
// (systag, jet radius) in ana::systags x ana::JetRs (src/ana.h) - the same systag set
// unfolder.cc/unfold_allsys.C already produce an insitu_tree for, and the same
// asymmetric/symmetric sign-split combination rule drawing/draw_systematics.C uses
// (gammajet_unfold/CLAUDE.md's ground rule): each ana::asymmetricSystagPairs member
// contributes its own signed (p_a_systag-p_a_nominal)/p_a_nominal to the "up" total if
// positive or "down" if negative; every other non-nominal systags entry is symmetric,
// contributing its full magnitude to both. Adding a systag to ana::systags (and, if it's
// a two-point pair, to ana::asymmetricSystagPairs) is the only thing needed for it to
// show up here automatically - see src/ana.h's comment on those two members.
//
// Scoped to grid_insitu.C's purity-corrected mean(x_J) fit specifically (not the shape-
// chi2 or unfolded methods) - it's the most statistically robust of the four in-situ
// methods (see the shape-chi2 spike/sawtooth investigation this session), and it's the
// one unfolder.cc's jesCorrection is actually derived from.
//
// Plot style matches gammajet/drawing/newdraw_all.C's final xj_vs_R.pdf page (systematic
// band as a semi-transparent TBox per radius, statistical uncertainty as a black
// TGraphAsymmErrors point+bar) rather than this project's usual band-graph convention,
// per explicit request to match that macro's look.
//
// Run insitu/run_grid.sh for every radius/systag first - this macro only reads existing
// grid_insitu_<systag>_<rname>.root files, it does not run any scan itself.

// insitu/ is split into output/ (grid_insitu.C's .root output, which this macro reads,
// plus this macro's own .root output) and pdfs/ (this macro's own .pdf output) - it
// reads no input ntuples of its own.
const char * insitu_output_dir = "/home/samson72/sphnx/gammajet_unfold/insitu/output";
const char * insitu_pdf_dir    = "/home/samson72/sphnx/gammajet_unfold/insitu/pdfs";

// Rewrites src/ana.h's jesNominal/jesTotalErrLow/jesTotalErrHigh array literals in
// place with this scan's freshly measured values - this is the "generate, don't
// hand-copy" replacement for what used to be a manual transcription of this macro's
// own console table into ana.h after every in-situ re-scan. Only touches the numeric
// literal inside each "= {...};" - the variable name, alignment padding before "=",
// and everything else in the file is left byte-for-byte untouched. Does NOT rebuild
// libgammajet_unfold.so itself (see CLAUDE.md's Build & Run section on why driver-style
// side effects like a full recompile aren't triggered automatically here) - the caller
// (run_grid.sh) or the user still needs to run src/make.sh afterward for anything
// linking ana::jesNominal/jesTotalErrLow/jesTotalErrHigh to see the new numbers.
void updateAnaHeader(const float pa[ana::nJetR], const float totalLow[ana::nJetR], const float totalHigh[ana::nJetR]) {
  const char * anaHeaderPath = "/home/samson72/sphnx/gammajet_unfold/src/ana.h";
  ifstream fin(anaHeaderPath);
  if (!fin) {
    cout << "WARNING: could not open " << anaHeaderPath << " to update the JES constants - left unchanged." << endl;
    return;
  }
  stringstream sbuf;
  sbuf << fin.rdbuf();
  string content = sbuf.str();
  fin.close();

  auto formatArray = [](const float * v) {
    string s = "{";
    for (int i = 0; i < ana::nJetR; i++) {
      s += Form("%.4f", v[i]);
      if (i+1 < ana::nJetR) s += ", ";
    }
    s += "}";
    return s;
  };

  bool allFound = true;
  auto replaceArrayLiteral = [&](const string & varName, const string & newLiteral) {
    string declTag = "double " + varName + "[nJetR]";
    size_t declPos = content.find(declTag);
    if (declPos == string::npos) {
      cout << "WARNING: could not find `" << declTag << "` in ana.h - " << varName << " not updated." << endl;
      allFound = false;
      return;
    }
    size_t eq = content.find('=', declPos);
    size_t semi = content.find(';', eq);
    if (eq == string::npos || semi == string::npos) {
      cout << "WARNING: malformed declaration for " << varName << " in ana.h - not updated." << endl;
      allFound = false;
      return;
    }
    content = content.substr(0, eq+1) + " " + newLiteral + content.substr(semi);
  };

  replaceArrayLiteral("jesNominal",      formatArray(pa));
  replaceArrayLiteral("jesTotalErrLow",  formatArray(totalLow));
  replaceArrayLiteral("jesTotalErrHigh", formatArray(totalHigh));

  if (!allFound) {
    cout << "ana.h left unchanged due to the warning(s) above." << endl;
    return;
  }

  ofstream fout(anaHeaderPath);
  fout << content;
  fout.close();
  cout << "Updated " << anaHeaderPath << "'s jesNominal/jesTotalErrLow/jesTotalErrHigh "
       << "with this scan's results. Run src/make.sh to rebuild before trusting anything "
       << "downstream." << endl;
}

struct PaResult { bool ok = false; float pa = 0, errLow = 0, errHigh = 0; };

PaResult readPa(const string & systag, int ir) {
  PaResult r;
  // grid_insitu.C now writes one file per systag (all seven jet radii inside, one
  // ana::rnames[ir] subdirectory each) instead of one file per (systag,radius) - see
  // that macro's header comment.
  const char * filename = Form("%s/grid_insitu_%s.root", insitu_output_dir, systag.c_str());
  TFile * f = TFile::Open(filename, "READ");
  if (!f || f->IsZombie()) {
    cout << "WARNING: could not open " << filename << " - run grid_insitu.C(\"" << systag
         << "\") first." << endl;
    return r;
  }
  TTree * t = (TTree*)f->Get(Form("%s/results", ana::rnames[ir]));
  if (!t) { f->Close(); return r; }
  float pa, errLow, errHigh;
  t->SetBranchAddress("pa_puritycorrected", &pa);
  t->SetBranchAddress("errLow_puritycorrected", &errLow);
  t->SetBranchAddress("errHigh_puritycorrected", &errHigh);
  t->GetEntry(0);
  f->Close();
  r.ok = true; r.pa = pa; r.errLow = errLow; r.errHigh = errHigh;
  return r;
}

void draw_jes_summary() {
  gStyle->SetOptStat(0);

  // Flatten ana::asymmetricSystagPairs into a lookup set of individual systag names -
  // derived, not duplicated, so a new pair added there is picked up here automatically.
  set<string> asymmetricSystags;
  for (auto & pr : ana::asymmetricSystagPairs) {
    asymmetricSystags.insert(pr.first);
    asymmetricSystags.insert(pr.second);
  }

  cout << Form("%-6s %10s %12s %12s %12s\n", "radius", "p_a", "stat -/+", "syst -/+", "total -/+");

  // Captured as a string immediately: the many Form() calls inside readPa() (called up
  // to 70 times below) rotate through and overwrite Form()'s static buffer pool long
  // before outfilename is used again at the very end of this function.
  string outfilename = Form("%s/draw_jes_summary.root", insitu_output_dir);
  TFile * fout = TFile::Open(outfilename.c_str(), "RECREATE");
  TTree * wt = new TTree("jes_summary", "in-situ JES scale factor per jet radius, with full systematic uncertainty");
  int wir; float wR, wpa, wstatLow, wstatHigh, wsystLow, wsystHigh, wtotalLow, wtotalHigh;
  wt->Branch("ir", &wir);
  wt->Branch("R", &wR);
  wt->Branch("pa", &wpa);
  wt->Branch("statLow", &wstatLow);
  wt->Branch("statHigh", &wstatHigh);
  wt->Branch("systLow", &wsystLow);
  wt->Branch("systHigh", &wsystHigh);
  wt->Branch("totalLow", &wtotalLow);
  wt->Branch("totalHigh", &wtotalHigh);

  // Per-radius results, collected here and used for both the tree and the plot below -
  // avoids re-deriving systLow/systHigh from the graphs after the fact.
  vector<float> vR, vPa, vStatLow, vStatHigh, vSystLow, vSystHigh, vTotalLow, vTotalHigh;

  // Per-systag coverage count across the radius loop below - used to gate the ana.h
  // auto-update at the end: a partial sweep (e.g. run_grid.sh --systag nominal, or a
  // radius/systag combination that just hasn't been scanned yet) must never silently
  // bake an underestimated systematic into ana.h, so every non-nominal systag needs a
  // valid result at every radius, not just "some".
  map<string,int> systagCoverage;
  for (const string & systag : ana::systags) if (systag != "nominal") systagCoverage[systag] = 0;

  for (int ir = 0; ir < ana::nJetR; ir++) {
    PaResult nom = readPa("nominal", ir);
    if (!nom.ok || nom.pa <= 0) {
      cout << "Skipping R=" << ana::JetRs[ir] << " (ir=" << ir << "): no nominal result." << endl;
      continue;
    }

    double sumsqUp = 0, sumsqDown = 0;
    for (const string & systag : ana::systags) {
      if (systag == "nominal") continue;
      PaResult var = readPa(systag, ir);
      if (!var.ok) continue; // already warned in readPa()
      systagCoverage[systag]++;
      double fracDiff = (var.pa - nom.pa)/nom.pa;
      if (asymmetricSystags.count(systag)) {
        if (fracDiff > 0) sumsqUp   += fracDiff*fracDiff;
        else              sumsqDown += fracDiff*fracDiff;
      } else {
        sumsqUp   += fracDiff*fracDiff;
        sumsqDown += fracDiff*fracDiff;
      }
    }
    double systLow  = nom.pa * sqrt(sumsqDown);
    double systHigh = nom.pa * sqrt(sumsqUp);
    double totalLow  = sqrt(nom.errLow*nom.errLow   + systLow*systLow);
    double totalHigh = sqrt(nom.errHigh*nom.errHigh + systHigh*systHigh);

    cout << Form("R=%.1f  %10.4f  -%.4f/+%.4f  -%.4f/+%.4f  -%.4f/+%.4f\n",
        ana::JetRs[ir], nom.pa, nom.errLow, nom.errHigh, systLow, systHigh, totalLow, totalHigh);

    vR.push_back(ana::JetRs[ir]);
    vPa.push_back(nom.pa);
    vStatLow.push_back(nom.errLow);
    vStatHigh.push_back(nom.errHigh);
    vSystLow.push_back(systLow);
    vSystHigh.push_back(systHigh);
    vTotalLow.push_back(totalLow);
    vTotalHigh.push_back(totalHigh);

    wir = ir; wR = ana::JetRs[ir]; wpa = nom.pa;
    wstatLow = nom.errLow; wstatHigh = nom.errHigh;
    wsystLow = systLow; wsystHigh = systHigh;
    wtotalLow = totalLow; wtotalHigh = totalHigh;
    wt->Fill();
  }

  int nOk = (int)vR.size();
  if (nOk == 0) {
    cout << "No radii had a nominal grid_insitu.C result - nothing to plot. Run "
         << "insitu/run_grid.sh first." << endl;
    fout->Close();
    return;
  }

  // Auto-update ana.h's jesNominal/jesTotalErrLow/jesTotalErrHigh iff this was a
  // complete sweep: every radius has a nominal result AND every ana::systags entry
  // (the full systematic set the totalLow/totalHigh above are quadrature-summed over)
  // has a result at every radius. Anything short of that is a partial/quick-check run
  // (e.g. run_grid_nominal.sh, or --systag) whose systLow/systHigh would understate the
  // real systematic - report what's missing instead of overwriting a real physics
  // constant with an incomplete one.
  bool complete = (nOk == ana::nJetR);
  vector<string> incompleteSystags;
  for (auto & kv : systagCoverage) if (kv.second != ana::nJetR) incompleteSystags.push_back(kv.first);
  complete = complete && incompleteSystags.empty();

  if (complete) {
    float paArr[ana::nJetR], totalLowArr[ana::nJetR], totalHighArr[ana::nJetR];
    for (int i = 0; i < ana::nJetR; i++) { paArr[i] = vPa[i]; totalLowArr[i] = vTotalLow[i]; totalHighArr[i] = vTotalHigh[i]; }
    updateAnaHeader(paArr, totalLowArr, totalHighArr);
  } else {
    cout << "Not updating src/ana.h: incomplete sweep (";
    if (nOk != ana::nJetR) cout << "only " << nOk << "/" << ana::nJetR << " radii have a nominal result";
    if (nOk != ana::nJetR && !incompleteSystags.empty()) cout << "; ";
    if (!incompleteSystags.empty()) {
      cout << "missing/partial systags: ";
      for (size_t i = 0; i < incompleteSystags.size(); i++) cout << (i?", ":"") << incompleteSystags[i];
    }
    cout << "). Run the full insitu/run_grid.sh sweep (all systags) first." << endl;
  }

  // -----------------------------
  // Summary plot: JES scale factor vs jet radius - same box(syst)+point(stat) style as
  // gammajet/drawing/newdraw_all.C's xj_vs_R.pdf - directly answers whether the in-situ
  // correction is consistent across radii (it was previously only ever measured at
  // R=0.4 and applied uniformly to every radius - see this session's weakness-finding
  // discussion).
  // -----------------------------
  TGraphAsymmErrors * gStat = new TGraphAsymmErrors(nOk);
  gStat->SetName("gJES_stat_vs_radius");
  for (int i = 0; i < nOk; i++) {
    gStat->SetPoint(i, vR[i], vPa[i]);
    gStat->SetPointError(i, 0, 0, vStatLow[i], vStatHigh[i]);
  }

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

  // Systematic uncertainty boxes - one per radius, centered on its point. Solid fill,
  // not SetFillColorAlpha - this ROOT build's batch-mode PNG/PDF output silently drops
  // alpha-blended fills entirely (verified with a standalone test), so
  // newdraw_all.C's semi-transparent style isn't reproducible here. SetFillStyle(1001)
  // is likewise required explicitly - a fresh TBox's default fill style renders as
  // fully hollow regardless of SetFillColor until this is set.
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

  // Statistical uncertainty points+bars, drawn on top of the systematic boxes.
  gStat->SetMarkerStyle(20);
  gStat->SetMarkerSize(1.2);
  gStat->SetLineWidth(2);
  gStat->SetLineColor(kBlack);
  gStat->SetMarkerColor(kBlack);
  gStat->Draw("PZ SAME");

  TLegend * leg = new TLegend(0.45, 0.72, 0.88, 0.85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(gStat, "Stat. unc.", "lep");
  leg->AddEntry(boxes[0], "Syst. unc.", "f");
  leg->Draw();

  insitu_utility::drawSPhenixLabel({"p+p Run24 Data"}, {"Pythia8 #gamma+jet MC", "Purity-corrected"}, .18, .85, 16, c->GetWh());

  c->RedrawAxis();
  const char * pdfPath = Form("%s/draw_jes_summary.pdf", insitu_pdf_dir);
  c->SaveAs(pdfPath);
  cout << "Wrote " << pdfPath << endl;

  // readPa()'s many TFile::Open()/Close() calls during the loop above leave gDirectory
  // pointing at gROOT ("Rint"), not fout, once the last of those files closes - cd()
  // back explicitly or Write() silently no-ops with a "not associated with a file" error
  // instead of writing into fout.
  fout->cd();
  gStat->Write();
  wt->Write();
  fout->Close();
  cout << "Wrote " << outfilename << endl;
}
