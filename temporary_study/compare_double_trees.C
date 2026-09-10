// Compares the reco-level xJgamma shape between the plain gammajet_pythia_Photon{5,10,20}
// MC trees and the gammajet_pythia_Photon{5,10,20}_double trees, both in
// /home/samson72/sphnx/gammajet/trees/, one page per reported photon-pT bin
// (ana::ptBinsUsed, the 3 "main" bins: 15-20, 20-25, 25-35 GeV).
//
// Both sets of trees have the standard (non-temporary_study) branch layout - confirmed by
// inspection, the _double trees do NOT carry the _old/_nosat branches that compare_old_new.C
// needs, so this only reproduces draw_xj_data_mc.C's single "Everything new" nominal reco
// case (cluster_pt/eta/phi, cluster_bdt_scores[9], cluster_showershape[11],
// jet_pt_smear_truth[ir] - MC-only, both sides are MC here so there's no Data-side jet_pt_calib
// branch needed). Same event selection (vz cut, truth-pt trigger-stitching window, xJ floor,
// |eta| cuts, dphi cut, ABCD region-A cut), same representative jet radius (ir=2, R=0.4), and
// same cross-section/vz/cluster-pT reweighting as draw_xj_data_mc.C/compare_old_new.C - see
// those files' header comments for the full rationale.
//
// Both histograms are unit-normalized (independently) before drawing - this is a shape
// comparison, not an absolute-rate one (the two productions may differ in raw statistics).

#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/reweight_utility.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"
R__LOAD_LIBRARY(libgammajet_unfold.so)

#include <string>
#include <map>
#include <utility>
#include <cmath>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"

const int ir = 2; // R=0.4, representative radius - same convention as compare_old_new.C

// Measured double-interaction rate in real data - used to build a third "weighted" xJ
// shape/mean, (1 - rate)*single + rate*double, representing the realistic mixture rather
// than either pure MC production on its own.
const double doubleInteractionRate = 0.185;

// truth-level cluster-pT trigger-stitching window, pythia sim, from
// treeuser.cc's threshmap[-1]/threshmap_high[-1] (Photon5/10/20 rows).
std::map<std::string, std::pair<float, float>> truthPtWindow = {
  {"Photon5",  {0.f,  12.f}},
  {"Photon10", {12.f, 24.f}},
  {"Photon20", {24.f, 100.f}}
};

// Per-sample cross-section/luminosity normalization - literal copy of drawer.h's
// scalemap[/*isphoton=*/1][sample] for sim=="pythia" (src/drawer.h's `drawer` ctor).
std::map<std::string, double> sampleWeight = {
  {"Photon5",  146359.3},
  {"Photon10", 6944.675},
  {"Photon20", 130.4461}
};

// Exact unbinned mean/error accumulator - same as draw_xj_data_mc.C's Accum.
struct Accum {
  double sumw = 0, sumw2 = 0, sumwx = 0, sumwx2 = 0;
  Long64_t n = 0;
  void fill(double x, double w) {
    sumw += w; sumw2 += w * w; sumwx += w * x; sumwx2 += w * x * x; n++;
  }
  double mean() const { return sumw > 0 ? sumwx / sumw : 0; }
};

// Mirrors unfolder::check_pair (src/unfolder.cc) at fixed ir=2, floorScale=1, testPt=-1.
bool checkPair(float jetPt, float jetEta, float phoPt, float phoEta, float dphi) {
  int ptbin = ana::findPtBin(phoPt);
  if (ptbin == -1) return false;
  float val = jetPt / phoPt;
  float lowval = ana::jet_calib_pt_cut[ir] / ana::ptBins[ptbin];
  float lowbin = ana::unfoldXjBins[ana::findUnfoldXjBin(lowval) + 1];
  if (val < lowbin) return false;
  if (fabs(phoEta) > ana::etacut) return false;
  if (fabs(jetEta) > ana::etacut - ana::JetRs[ir]) return false;
  if (dphi < ana::oppcut) return false;
  return true;
}

float deltaPhi(float p1, float p2) {
  float dphi = fabs(p1 - p2);
  if (dphi > M_PI) dphi = 2 * M_PI - dphi;
  return dphi;
}

// Fills hxjPt[ip] (ana::nPtBinsUsed histograms, ana::unfoldXjBins binning) and accPt[ip]
// (exact unbinned mean) with the "Everything new" nominal-reco xJ, reading trigger's tree
// from fname.
void processFile(const std::string & fname, const std::string & trigger,
    std::vector<TH1D*> & hxjPt, std::vector<Accum> & accPt, Reweighter & rw,
    Long64_t maxEntries = -1) {

  TFile * f = TFile::Open(fname.c_str(), "read");
  if (!f || f->IsZombie()) { std::cout << "Could not open " << fname << std::endl; return; }
  TTree * t = (TTree*)f->Get("towerntup");
  Long64_t nentries = t->GetEntries();
  if (maxEntries >= 0 && maxEntries < nentries) nentries = maxEntries;
  std::cout << "Processing " << fname << " (" << nentries << " entries, MC)" << std::endl;

  t->SetBranchStatus("*", 0);
  std::vector<std::string> branches = {"vz", "cluster_pt", "cluster_eta", "cluster_phi",
       "cluster_showershape", "cluster_bdt_scores", "jet_eta", "jet_phi",
       "truth_cluster_pt", "jet_pt_smear_truth"};
  for (const auto & bn : branches) t->SetBranchStatus(bn.c_str(), 1);

  Float_t vz, cluster_pt, cluster_eta, cluster_phi;
  Float_t cluster_showershape[12], cluster_bdt_scores[11];
  Float_t truth_cluster_pt = 0;
  Float_t jet_pt_smear_truth[ana::nJetR] = {0};
  Float_t jet_eta[ana::nJetR], jet_phi[ana::nJetR];

  t->SetBranchAddress("vz", &vz);
  t->SetBranchAddress("cluster_pt", &cluster_pt);
  t->SetBranchAddress("cluster_eta", &cluster_eta);
  t->SetBranchAddress("cluster_phi", &cluster_phi);
  t->SetBranchAddress("cluster_showershape", cluster_showershape);
  t->SetBranchAddress("cluster_bdt_scores", cluster_bdt_scores);
  t->SetBranchAddress("jet_eta", jet_eta);
  t->SetBranchAddress("jet_phi", jet_phi);
  t->SetBranchAddress("truth_cluster_pt", &truth_cluster_pt);
  t->SetBranchAddress("jet_pt_smear_truth", jet_pt_smear_truth);

  auto windowIt = truthPtWindow.find(trigger);
  bool haveWindow = windowIt != truthPtWindow.end();
  double xsecWeight = sampleWeight.at(trigger);

  for (Long64_t e = 0; e < nentries; e++) {
    t->GetEntry(e);
    if (e % 2000000 == 0)
      std::cout << "  entry " << e << "/" << nentries
                << " (" << (float)e / nentries * 100. << "%)" << std::endl;

    if (fabs(vz) > ana::vzcut) continue;
    if (haveWindow && !(truth_cluster_pt > windowIt->second.first && truth_cluster_pt < windowIt->second.second))
      continue;

    float mcWeight = rw.GetWeight(vz, cluster_pt) * xsecWeight;

    float jetEta = jet_eta[ir];
    float jetPhi = jet_phi[ir];
    // "Everything new": cluster_pt/eta/phi (with-saturation), cluster_bdt_scores[9]
    // (new bdt), cluster_showershape[11] (new iso), jet_pt_smear_truth[ir] (new jet pt) -
    // same nominal reco definition as draw_xj_data_mc.C's MC case.
    float phoPt = cluster_pt, phoEta = cluster_eta, phoPhi = cluster_phi;
    float jetPt = jet_pt_smear_truth[ir];
    float iso = cluster_showershape[11];
    float bdt = cluster_bdt_scores[9];

    if (!std::isfinite(phoPt) || phoPt <= 0 || !std::isfinite(jetPt)) continue;
    if (!(phoPt >= ana::ptBins[0] && phoPt < ana::ptBins[ana::nPtBins] && jetPt > ana::jet_calib_pt_cut[ir]))
      continue;
    float dphi = deltaPhi(phoPhi, jetPhi);
    if (!checkPair(jetPt, jetEta, phoPt, phoEta, dphi)) continue;

    int iabcd = ana::findabcdBin(iso, bdt, 0);
    if (iabcd != 0) continue;

    float xj = jetPt / phoPt;
    if (!std::isfinite(xj)) continue;

    int ptbin = ana::findPtBin(phoPt);
    int iused = ptbin - ana::firstUsedPtBin;
    if (iused >= 0 && iused < ana::nPtBinsUsed) {
      hxjPt[iused]->Fill(xj, mcWeight);
      accPt[iused].fill(xj, mcWeight);
    }
  }
  f->Close();
}

void compare_double_trees(Long64_t maxEntriesPerFile = -1) {
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();

  Reweighter rw;

  // drawer::drawAll (used below for the sPHENIX label) only needs its member function, not
  // the per-sample unfolding-output files its constructor opens - see compare_old_new.C's
  // identical comment for why the resulting TFile::Open error spam is silenced here.
  int oldErrLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  drawer d("pythia", "nominal");
  gErrorIgnoreLevel = oldErrLevel;

  int nb = ana::nUnfoldXjBins;
  const double * edges = ana::unfoldXjBins;
  std::string jetFeature = Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV", ana::JetRs[ir], ana::jet_calib_pt_cut[ir]);

  std::vector<TH1D*> hxjOrig(ana::nPtBinsUsed), hxjDouble(ana::nPtBinsUsed);
  std::vector<Accum> accOrig(ana::nPtBinsUsed), accDouble(ana::nPtBinsUsed);
  for (int ip = 0; ip < ana::nPtBinsUsed; ip++) {
    hxjOrig[ip] = new TH1D(Form("hxjOrig_pt%d", ip), ";x_{J#gamma};normalized counts", nb, edges);
    hxjOrig[ip]->SetStatOverflows(TH1::kConsider);
    hxjOrig[ip]->SetLineColor(kMagenta + 1);
    hxjOrig[ip]->SetLineWidth(2);

    hxjDouble[ip] = new TH1D(Form("hxjDouble_pt%d", ip), ";x_{J#gamma};normalized counts", nb, edges);
    hxjDouble[ip]->SetStatOverflows(TH1::kConsider);
    hxjDouble[ip]->SetLineColor(kGreen + 2);
    hxjDouble[ip]->SetLineWidth(2);
  }

  for (const std::string & trigger : {"Photon5", "Photon10", "Photon20"}) {
    processFile("/home/samson72/sphnx/gammajet/trees/gammajet_pythia_" + trigger + ".root",
        trigger, hxjOrig, accOrig, rw, maxEntriesPerFile);
    processFile("/home/samson72/sphnx/gammajet/trees/gammajet_pythia_" + trigger + "_double.root",
        trigger, hxjDouble, accDouble, rw, maxEntriesPerFile);
  }

  std::cout << "\n" << Form("%-22s %10s %10s %10s", "pT bin", "<xJ> orig", "<xJ> double", "double/orig") << std::endl;
  for (int ip = 0; ip < ana::nPtBinsUsed; ip++) {
    double mOrig = accOrig[ip].mean(), mDouble = accDouble[ip].mean();
    std::cout << Form("%-22s %10.4f %10.4f %10.4f",
        Form("%.0f-%.0f GeV", ana::ptBinsUsed[ip], ana::ptBinsUsed[ip + 1]), mOrig, mDouble,
        (mOrig > 0 && mDouble > 0) ? mDouble / mOrig : 0) << std::endl;
  }

  std::string pdfPath = "/home/samson72/sphnx/gammajet_unfold/temporary_study/compare_double_trees.pdf";
  TCanvas * c = new TCanvas("c_compare_double", "", 700, 600);
  c->SetLeftMargin(.13);
  c->SetBottomMargin(.13);
  c->SaveAs((pdfPath + "[").c_str());

  for (int ip = 0; ip < ana::nPtBinsUsed; ip++) {
    c->Clear();
    // ana::unfoldXjBins is variable-width - Integral() omits "width" (total event/weight
    // count, correct normalization denominator), Scale's "width" option makes bin content
    // a density so wider tail bins aren't shown too tall relative to narrower ones.
    double integOrig = hxjOrig[ip]->Integral(0, hxjOrig[ip]->GetNbinsX() + 1);
    double integDouble = hxjDouble[ip]->Integral(0, hxjDouble[ip]->GetNbinsX() + 1);
    if (integOrig > 0) hxjOrig[ip]->Scale(1.0 / integOrig, "width");
    if (integDouble > 0) hxjDouble[ip]->Scale(1.0 / integDouble, "width");

    // Weighted mixture: (1 - rate)*single + rate*double, both terms already unit-density
    // normalized above, so the combination integrates to 1 as well - see file header
    // comment on doubleInteractionRate.
    TH1D * hxjWeighted = (TH1D*)hxjOrig[ip]->Clone(Form("hxjWeighted_pt%d", ip));
    hxjWeighted->Scale(1.0 - doubleInteractionRate);
    hxjWeighted->Add(hxjDouble[ip], doubleInteractionRate);
    hxjWeighted->SetLineColor(kRed + 1);
    hxjWeighted->SetLineWidth(2);
    double meanWeighted = (1.0 - doubleInteractionRate) * accOrig[ip].mean()
                         + doubleInteractionRate * accDouble[ip].mean();

    double ymax = std::max({hxjOrig[ip]->GetMaximum(), hxjDouble[ip]->GetMaximum(), hxjWeighted->GetMaximum()});
    hxjOrig[ip]->SetTitle(";x_{J#gamma};(1/N) dN/dx_{J#gamma}");
    hxjOrig[ip]->GetYaxis()->SetRangeUser(0, ymax * 1.5);
    hxjOrig[ip]->Draw("hist");
    hxjDouble[ip]->Draw("hist same");
    hxjWeighted->Draw("hist same");

    TLegend * leg = new TLegend(0.47, 0.55, 0.88, 0.83);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    // <x_J> from the exact unbinned Accum, NOT TH1::GetMean() on the width-scaled
    // histogram above (that would implicitly re-weight bins by 1/width).
    leg->AddEntry(hxjOrig[ip], Form("Single interaction  #LTx_{J}#GT=%.3f", accOrig[ip].mean()), "l");
    leg->AddEntry(hxjDouble[ip], Form("Double interaction  #LTx_{J}#GT=%.3f", accDouble[ip].mean()), "l");
    leg->AddEntry(hxjWeighted, Form("Weighted (%.1f%% double)  #LTx_{J}#GT=%.3f", doubleInteractionRate * 100, meanWeighted), "l");
    leg->Draw();

    double meanRatio = (accOrig[ip].mean() > 0) ? accDouble[ip].mean() / accOrig[ip].mean() : 0;
    double weightedRatio = (accOrig[ip].mean() > 0) ? meanWeighted / accOrig[ip].mean() : 0;
    TLatex latexRatio;
    latexRatio.SetNDC();
    latexRatio.SetTextSize(0.028);
    latexRatio.DrawLatex(0.52, 0.515, Form("Double/Single #LTx_{J}#GT ratio = %.3f", meanRatio));
    latexRatio.DrawLatex(0.52, 0.48, Form("Weighted/Single #LTx_{J}#GT ratio = %.3f", weightedRatio));

    d.drawAll({"Pythia8 #gamma+jet"},
        {Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", ana::ptBinsUsed[ip], ana::ptBinsUsed[ip + 1]), jetFeature},
        .17, .87, 14, gPad->GetWh() * 0.8);
    c->SaveAs(pdfPath.c_str());
  }

  c->SaveAs((pdfPath + "]").c_str());
  std::cout << "Wrote " << pdfPath << std::endl;
}
