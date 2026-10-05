#include "../src/ana.h"
#include "../src/insitu_utility.h"
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

// Final in-situ JES p_a per jet radius, from grid_insitu.C's purity-corrected mean fit for every
// (systag, radius). Systematics combine as in draw_systematics.C (asymmetricSystagPairs sign-split
// per source, others symmetric). Writes into src/ana.h the nominal p_a, its statistical error
// (jesStatErrLow/High, used by jes_high/low) and the per-systag table jesBySystag. The quadrature
// systematic is drawn but not fed back: each variation already carries its JES effect.
// Run insitu/run_grid.sh first.

const char * insitu_output_dir = ana::path("insitu/output");
const char * insitu_pdf_dir    = ana::path("insitu/pdfs");

// Rewrites the jesNominal/jesStatErrLow/jesStatErrHigh/jesBySystag literals in src/ana.h in
// place. Does not rebuild; run src/make.sh afterwards.
void updateAnaHeader(const float pa[ana::nJetR], const float statLow[ana::nJetR], const float statHigh[ana::nJetR],
    const map<string, vector<float>> & paBySystag) {
  const char * anaHeaderPath = ana::path("src/ana.h");
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

  replaceArrayLiteral("jesNominal",     formatArray(pa));
  replaceArrayLiteral("jesStatErrLow",  formatArray(statLow));
  replaceArrayLiteral("jesStatErrHigh", formatArray(statHigh));

  // Regenerate the whole jesBySystag block (between its markers) in ana::systags order.
  {
    const string beginTag = "    // BEGIN jesBySystag\n", endTag = "    // END jesBySystag\n";
    size_t b = content.find(beginTag), e = content.find(endTag);
    if (b == string::npos || e == string::npos || e < b) {
      cout << "WARNING: could not find the BEGIN/END jesBySystag markers in ana.h - table not updated." << endl;
      allFound = false;
    } else {
      string block = beginTag;
      block += Form("    static constexpr int nJesSystags = %d;\n", (int)ana::systags.size());
      block += "    static constexpr const char * jesSystagNames[nJesSystags] = {";
      for (size_t i = 0; i < ana::systags.size(); i++) block += (i ? ", \"" : "\"") + ana::systags[i] + "\"";
      block += "};\n    static constexpr double jesBySystag[nJesSystags][nJetR] = {\n";
      for (size_t i = 0; i < ana::systags.size(); i++) {
        block += "      " + formatArray(paBySystag.at(ana::systags[i]).data()) + ", // " + ana::systags[i] + "\n";
      }
      block += "    };\n";
      content = content.substr(0, b) + block + content.substr(e);
    }
  }

  if (!allFound) {
    cout << "ana.h left unchanged due to the warning(s) above." << endl;
    return;
  }

  ofstream fout(anaHeaderPath);
  fout << content;
  fout.close();
  cout << "Updated " << anaHeaderPath << "'s jesNominal/jesStatErrLow/jesStatErrHigh/jesBySystag "
       << "with this scan's results. Run src/make.sh to rebuild before trusting anything "
       << "downstream." << endl;
}

struct PaResult { bool ok = false; float pa = 0, errLow = 0, errHigh = 0; };

PaResult readPa(const string & systag, int ir) {
  PaResult r;
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

  set<string> asymmetricSystags;
  for (auto & pr : ana::asymmetricSystagPairs) {
    asymmetricSystags.insert(pr.first);
    asymmetricSystags.insert(pr.second);
  }

  cout << Form("%-6s %10s %12s %12s %12s\n", "radius", "p_a", "stat -/+", "syst -/+", "total -/+");

  // Copy now: readPa()'s Form() calls overwrite Form()'s buffer pool.
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

  vector<float> vR, vPa, vStatLow, vStatHigh, vSystLow, vSystHigh, vTotalLow, vTotalHigh;

  // Per-systag radius coverage: ana.h is updated only from a complete sweep.
  map<string,int> systagCoverage;
  for (const string & systag : ana::systags) if (systag != "nominal") systagCoverage[systag] = 0;
  // p_a per (systag, radius) for jesBySystag, nominal included.
  map<string, vector<float>> paBySystag;
  for (const string & systag : ana::systags) paBySystag[systag] = vector<float>(ana::nJetR, 0);

  for (int ir = 0; ir < ana::nJetR; ir++) {
    PaResult nom = readPa("nominal", ir);
    if (!nom.ok || nom.pa <= 0) {
      cout << "Skipping R=" << ana::JetRs[ir] << " (ir=" << ir << "): no nominal result." << endl;
      continue;
    }

    paBySystag["nominal"][ir] = nom.pa;
    double sumsqUp = 0, sumsqDown = 0;
    for (const string & systag : ana::systags) {
      if (systag == "nominal") continue;
      PaResult var = readPa(systag, ir);
      if (!var.ok) continue; // already warned in readPa()
      systagCoverage[systag]++;
      paBySystag[systag][ir] = var.pa;
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

  // Update ana.h only if every ana::systags entry has a result at every radius.
  bool complete = (nOk == ana::nJetR);
  vector<string> incompleteSystags;
  for (auto & kv : systagCoverage) if (kv.second != ana::nJetR) incompleteSystags.push_back(kv.first);
  complete = complete && incompleteSystags.empty();

  if (complete) {
    float paArr[ana::nJetR], statLowArr[ana::nJetR], statHighArr[ana::nJetR];
    for (int i = 0; i < ana::nJetR; i++) { paArr[i] = vPa[i]; statLowArr[i] = vStatLow[i]; statHighArr[i] = vStatHigh[i]; }
    updateAnaHeader(paArr, statLowArr, statHighArr, paBySystag);
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
  // Summary plot: JES scale vs jet radius (box = syst, point = stat)
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

  insitu_utility::drawSPhenixLabel({"p+p Run24 Data"}, {"Pythia8 #gamma+jet MC", "Purity-corrected"}, .18, .85, 16, c->GetWh());

  c->RedrawAxis();
  const char * pdfPath = Form("%s/draw_jes_summary.pdf", insitu_pdf_dir);
  c->SaveAs(pdfPath);
  cout << "Wrote " << pdfPath << endl;

  // The TFile::Open/Close calls in readPa() leave gDirectory at gROOT; cd() back before Write().
  fout->cd();
  gStat->Write();
  wt->Write();
  fout->Close();
  cout << "Wrote " << outfilename << endl;
}
