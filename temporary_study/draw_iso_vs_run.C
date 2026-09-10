// Data-only QA plot: mean photon isolation energy vs. RunNumber, old iso
// (cluster_showershape[9]) vs. new iso (cluster_showershape[11]) - same two indices
// compare_old_new.C's case 1/4 ("old") and case 2/5 ("new") read, from the single default
// (unsuffixed) cluster_showershape array - see that file's header comment on why iso is
// always read from there regardless of which cluster/BDT branch set is otherwise in play.
//
// Selection is deliberately just a basic photon-candidate quality cut (vz cut,
// cluster_pt>0, cluster_pt in ana::ptBins range, |cluster_eta|<ana::etacut) - NOT the full
// photon+jet pairing selection (no jet pt cut, no dphi cut, no ABCD region cut) used
// elsewhere in this study. Isolation is a photon-level shower-shape variable independent
// of any jet leg, and conditioning on the ABCD region in particular would be circular
// (region A is defined by a cut on iso itself), so this macro looks at the same
// photon-candidate population old/new iso would each see, not the reduced photon+jet
// sample. Events where a given iso variant reports the "topocluster iso not computed"
// sentinel (<=-999, ana::findabcdBin's isiso==-1 case) are excluded from that variant's
// mean only - old/new validity is not required to agree event-by-event.
//
// One point per distinct RunNumber actually present in the file (1565 runs found,
// spanning 47289-53864, with per-run statistics ranging from 1 to ~5200 clusters) -
// plotted at its true run number (not a sequential run index), so gaps in the x-axis
// reflect genuine gaps in which runs are in this file.

#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"
R__LOAD_LIBRARY(libgammajet_unfold.so)

#include <map>
#include <cmath>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

// Per-run mean/error accumulator (unbinned) - same convention as compare_old_new.C's
// Accum, unweighted (Data, weight always 1).
struct Accum {
  double sum = 0, sum2 = 0;
  Long64_t n = 0;
  void fill(double x) { sum += x; sum2 += x * x; n++; }
  double mean() const { return n > 0 ? sum / n : 0; }
  double meanErr() const {
    if (n < 2) return 0;
    double var = sum2 / n - mean() * mean();
    if (var < 0) var = 0;
    return std::sqrt(var / n);
  }
};

void draw_iso_vs_run(Long64_t maxEntries = -1) {
  gStyle->SetOptStat(0);

  int oldErrLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  drawer d("pythia", "nominal");
  gErrorIgnoreLevel = oldErrLevel;

  std::string fname = "/home/samson72/sphnx/gammajet/trees/gammajet_Data.root";
  TFile * f = TFile::Open(fname.c_str(), "read");
  if (!f || f->IsZombie()) { std::cout << "Could not open " << fname << std::endl; return; }
  TTree * t = (TTree*)f->Get("towerntup");
  Long64_t nentries = t->GetEntries();
  if (maxEntries >= 0 && maxEntries < nentries) nentries = maxEntries;
  std::cout << "Processing " << fname << " (" << nentries << " entries, Data)" << std::endl;

  t->SetBranchStatus("*", 0);
  for (const char * bn : {"vz", "cluster_pt", "cluster_eta", "cluster_showershape", "RunNumber"})
    t->SetBranchStatus(bn, 1);

  Float_t vz, cluster_pt, cluster_eta, cluster_showershape[12];
  Int_t runNumber;
  t->SetBranchAddress("vz", &vz);
  t->SetBranchAddress("cluster_pt", &cluster_pt);
  t->SetBranchAddress("cluster_eta", &cluster_eta);
  t->SetBranchAddress("cluster_showershape", cluster_showershape);
  t->SetBranchAddress("RunNumber", &runNumber);

  std::map<int, Accum> oldByRun, newByRun;

  for (Long64_t e = 0; e < nentries; e++) {
    t->GetEntry(e);
    if (e % 2000000 == 0)
      std::cout << "  entry " << e << "/" << nentries
                << " (" << (float)e / nentries * 100. << "%)" << std::endl;

    if (fabs(vz) > ana::vzcut) continue;
    if (cluster_pt <= 0 || !std::isfinite(cluster_pt)) continue;
    if (!(cluster_pt >= ana::ptBins[0] && cluster_pt < ana::ptBins[ana::nPtBins])) continue;
    if (fabs(cluster_eta) > ana::etacut) continue;

    float isoOld = cluster_showershape[9];
    float isoNew = cluster_showershape[11];
    if (isoOld > -999 && std::isfinite(isoOld)) oldByRun[runNumber].fill(isoOld);
    if (isoNew > -999 && std::isfinite(isoNew)) newByRun[runNumber].fill(isoNew);
  }
  f->Close();

  std::cout << "\nDistinct runs: old iso = " << oldByRun.size() << ", new iso = " << newByRun.size() << std::endl;

  TGraphErrors * gOld = new TGraphErrors();
  gOld->SetName("gIsoOld");
  int i = 0;
  for (auto & kv : oldByRun) {
    gOld->SetPoint(i, kv.first, kv.second.mean());
    gOld->SetPointError(i, 0, kv.second.meanErr());
    i++;
  }
  TGraphErrors * gNew = new TGraphErrors();
  gNew->SetName("gIsoNew");
  i = 0;
  for (auto & kv : newByRun) {
    gNew->SetPoint(i, kv.first, kv.second.mean());
    gNew->SetPointError(i, 0, kv.second.meanErr());
    i++;
  }

  gOld->SetLineColor(kBlack);
  gOld->SetMarkerColor(kBlack);
  gOld->SetMarkerStyle(20);
  gOld->SetMarkerSize(0.6);
  gNew->SetLineColor(kRed + 1);
  gNew->SetMarkerColor(kRed + 1);
  gNew->SetMarkerStyle(21);
  gNew->SetMarkerSize(0.6);

  TCanvas * c = new TCanvas("c_iso_vs_run", "", 900, 600);
  c->SetLeftMargin(.12);
  c->SetBottomMargin(.13);
  std::string pdfPath = "/home/samson72/sphnx/gammajet_unfold/temporary_study/draw_iso_vs_run.pdf";
  c->SaveAs((pdfPath + "[").c_str());

  // ---------- page 1: old + new overlaid ----------
  c->Clear();
  TMultiGraph * mg = new TMultiGraph();
  mg->Add(gOld, "P");
  mg->Add(gNew, "P");
  mg->SetTitle(";Run Number;#LTIsolation Energy#GT [GeV]");
  mg->Draw("A P");

  TLegend * leg = new TLegend(0.6, 0.75, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(gOld, "Old isolation", "p");
  leg->AddEntry(gNew, "New isolation", "p");
  leg->Draw();

  d.drawAll({"p+p Run24 Data"}, {"Photon candidates, no jet requirement"}, .15, .87, 14, gPad->GetWh() * 0.8);
  c->SaveAs(pdfPath.c_str());

  // ---------- page 2: old alone ----------
  c->Clear();
  gOld->SetTitle(";Run Number;#LTIsolation Energy#GT [GeV]");
  gOld->Draw("A P");
  TLegend * legOld = new TLegend(0.6, 0.8, 0.88, 0.88);
  legOld->SetBorderSize(0);
  legOld->SetFillStyle(0);
  legOld->AddEntry(gOld, "Old isolation", "p");
  legOld->Draw();
  d.drawAll({"p+p Run24 Data"}, {"Photon candidates, no jet requirement"}, .15, .87, 14, gPad->GetWh() * 0.8);
  c->SaveAs(pdfPath.c_str());

  // ---------- page 3: new alone ----------
  c->Clear();
  gNew->SetTitle(";Run Number;#LTIsolation Energy#GT [GeV]");
  gNew->Draw("A P");
  TLegend * legNew = new TLegend(0.6, 0.8, 0.88, 0.88);
  legNew->SetBorderSize(0);
  legNew->SetFillStyle(0);
  legNew->AddEntry(gNew, "New isolation", "p");
  legNew->Draw();
  d.drawAll({"p+p Run24 Data"}, {"Photon candidates, no jet requirement"}, .15, .87, 14, gPad->GetWh() * 0.8);
  c->SaveAs(pdfPath.c_str());

  c->SaveAs((pdfPath + "]").c_str());
  std::cout << "Wrote " << pdfPath << std::endl;
}
