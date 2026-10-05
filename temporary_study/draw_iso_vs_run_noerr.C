// Same as draw_iso_vs_run.C (mean photon isolation energy vs. RunNumber, Data only, old
// iso cluster_showershape[9] vs. new iso cluster_showershape[11]) but without error bars -
// plain TGraph instead of TGraphErrors. See draw_iso_vs_run.C's header comment for the
// full selection/binning rationale (basic photon-candidate quality cut only, no ABCD
// region cut, one point per distinct RunNumber at its true run number).

#include "../src/ana.h"
#include "../src/drawer.h"
R__LOAD_LIBRARY(libgammajet_unfold.so)

#include <map>
#include <cmath>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

void draw_iso_vs_run_noerr(Long64_t maxEntries = -1) {
  gStyle->SetOptStat(0);

  int oldErrLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  drawer d("pythia", "nominal");
  gErrorIgnoreLevel = oldErrLevel;

  std::string fname = ana::path("trees/gammajet_Data.root");
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

  // sum only - mean() is all this version needs (no error bars drawn).
  std::map<int, std::pair<double, Long64_t>> oldByRun, newByRun;

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
    if (isoOld > -999 && std::isfinite(isoOld)) {
      auto & acc = oldByRun[runNumber];
      acc.first += isoOld; acc.second++;
    }
    if (isoNew > -999 && std::isfinite(isoNew)) {
      auto & acc = newByRun[runNumber];
      acc.first += isoNew; acc.second++;
    }
  }
  f->Close();

  std::cout << "\nDistinct runs: old iso = " << oldByRun.size() << ", new iso = " << newByRun.size() << std::endl;

  TGraph * gOld = new TGraph();
  gOld->SetName("gIsoOld");
  int i = 0;
  for (auto & kv : oldByRun) {
    gOld->SetPoint(i, kv.first, kv.second.first / kv.second.second);
    i++;
  }
  TGraph * gNew = new TGraph();
  gNew->SetName("gIsoNew");
  i = 0;
  for (auto & kv : newByRun) {
    gNew->SetPoint(i, kv.first, kv.second.first / kv.second.second);
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

  TCanvas * c = new TCanvas("c_iso_vs_run_noerr", "", 900, 600);
  c->SetLeftMargin(.12);
  c->SetBottomMargin(.13);
  std::string pdfPath = ana::path("temporary_study/draw_iso_vs_run_noerr.pdf");
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
