// Reco-level xJ comparison across six photon/jet reconstruction variants, using the
// extra "_old"/"_nosat" branches present only in temporary_study's ttrees (they are NOT
// in the standard gammajet trees src/treeuser.cc reads). Reimplements the same event
// selection unfolder.cc::check_pair/fill_matrix apply (vz cut, truth-pt trigger-stitching
// window, xJ floor, |eta| cuts, dphi cut, ABCD region-A cut) at a single representative
// jet radius (ir=2, R=0.4 - same "representative radius" convention as ana::getPurity's
// default ir and unfolder.cc's hpurity_num/den comment), rather than instantiating the
// full `unfolder` class (whose constructor opens the standard gammajet trees + insitu
// output files, not these one-off temporary_study trees).
//
// Branch mapping for the five cases (locked down with the user - see conversation, the
// old-BDT/nosat-kinematics pairing in cases 1-3 is a known, accepted inconsistency:
// cluster_bdt_score_old was trained on with-saturation inputs, not the _nosat ones):
//   1. AllOld:       cluster_{pt,eta,phi}_nosat, cluster_bdt_score_old,
//                     cluster_showershape[9] (iso), jet_pt_old
//   2. NewBDT+ISO:   same cluster_nosat kinematics + jet_pt_old, but
//                     cluster_bdt_scores[9] (bdt) and cluster_showershape[11] (iso)
//   3. NewJES+JER:   same as AllOld, but the "new" jet pt (see below) instead of jet_pt_old
//   4. WithPixelSat: same as AllOld, but cluster_{pt,eta,phi} (default, with-saturation)
//                     instead of the _nosat branches
//   5. AllNew:       cluster_{pt,eta,phi} (default, with-saturation), cluster_bdt_scores[9],
//                     cluster_showershape[11] (iso), the "new" jet pt - i.e. every axis
//                     flipped to new at once; this is the actual unfolder.cc nominal reco
//                     definition (see its fill_matrix()'s maxpho/maxjet construction)
//
// iso is always read from the default (unsuffixed) cluster_showershape[12] array,
// independent of which cluster/BDT branch set is otherwise in play (also locked down
// with the user).
//
// MC vs. Data ("new" jet pt, and how the two samples are combined): MC's "new" jet pt
// (cases 3/5) is jet_pt_smear_truth[ir] (no Data equivalent - there's no truth to smear
// against); Data's is jet_pt_calib[ir] directly, taken as-is from the ttree with no
// further correction applied (deliberately NOT unfolder.cc's own jet_pt_calib[ir]/
// ana::jesNominal[ir] - this study compares what's actually in the branches, not a
// rederived/recorrected quantity). jet_pt_old (cases 1/2/4) is assumed to be a genuine
// old-calibration branch present for both MC and Data (unlike jet_pt_smear_truth, it
// isn't inherently MC-only).
// Similarly, MC-only inputs (truth_cluster_pt's trigger-stitching window, the per-sample
// cross-section weight combining Photon5/10/20 - see truthPtWindow/sampleWeight below) are
// skipped entirely for Data (weight 1, single file, no truth branch read at all).
//
// Data support: once ~/sphnx/gammajet_unfold/trees/gammajet_Data.root exists (same tree format/branch
// names as the MC files here, "towerntup"), this macro auto-detects it (gSystem::
// AccessPathName) and additionally writes the same set of plots for Data
// (compare_old_new_data.pdf) plus a Data/MC ratio set (compare_old_new_ratio.pdf) - see
// compare_old_new()'s tail. Until then it silently skips both and behaves exactly as the
// MC-only version did.
//
// Data has NO cluster_{pt,eta,phi}_nosat branches at all (confirmed: only cluster_pt_old,
// unlike MC's separate cluster_pt_old/cluster_pt_nosat pair) - there was never a
// no-saturation reprocessing done for Data. cluster_pt_old == cluster_pt bit-for-bit in
// Data (checked all 1.3M nonzero entries), so photon kinematics genuinely cannot vary
// across cases for Data; only bdt/iso/jet-pt choice can. Locked down with the user:
// for Data only, cases 1/3 (both pair with the old BDT/iso) read cluster_pt_old/eta_old/
// phi_old instead of the nonexistent _nosat branches; case 2 (new BDT/iso) reads plain
// cluster_pt/eta/phi. Case 3 still uses jet_pt_calib (the "new" jet pt) for Data, same as
// before - only the photon-kinematics substitution changed, not the jet-pt one. MC is
// completely unaffected (all three still read cluster_pt_nosat, exactly as before).
//
// Case 6 (kPurityCorr, "Everything new, purity-corrected"): locked down with the user -
// case 5's ("Everything new") own kinematics/bdt/iso are reused unchanged; MC "does
// nothing different" (case 6 == case 5 exactly, literal duplicate, no correction applied -
// there's no ABCD background to subtract in MC truth-matched samples). Data instead gets
// the two-purity background subtraction insitu/grid_insitu.C applies: case 5's selection
// but also keeping its Region-C (findabcdBin==2: good iso, bad bdt) sibling events, then
// combining Region A and Region C per used-pT-bin via unfold_utility::purityCorrectCoeffs/
// purityCorrect (same functions grid_insitu.C's purityCorrectByPtBin calls), with purity
// P_A/P_C from ana::getPurity/getPurityC (+ErrorLow/ErrorHigh) at systag="nominal", this
// file's ir=2 - i.e. the same purity_nominal.root inputs, not re-derived here. Purity is
// only measured per ana::ptBinsUsed bin (see insitu_utility.h's comment on why the
// buffer/overflow pT bins are excluded from the in-situ study) so case 6's Data
// "Inclusive" row is the sum of just the 3 reported pT bins, NOT the full
// ptBins[0]-ptBins[nPtBins] range cases 1-5's inclusive row covers - see
// writeComparisonOutput's report() footnote.
//
// Plotted histograms use ana.h's own binning (ana::unfoldXjBins, 0-2.0) per the user's
// request; the reported means (in the summary table and the legend) are still accumulated
// unbinned event-by-event (not from TH1::GetMean() on a range-clipped histogram), so a
// real high-xJ tail beyond 2.0 still shows up in the mean even though it's not resolved
// bin-by-bin in the plot (it lands in the overflow bin, counted via SetStatOverflows so
// it's not silently dropped either).
//
// Photon5/Photon10/Photon20 are separate pT-hat-sliced MC productions with wildly
// different cross sections/statistics (see truthPtWindow below for how double-counting
// across the slices is avoided) - combining them with equal per-event weight would let
// whichever sample has the most raw entries dominate. Instead each sample's events are
// scaled by the same per-sample cross-section/luminosity normalization drawer.cc's
// combineMC() applies via scalemap[isphoton=1][...] (src/drawer.h) before summing -
// see sampleWeight below. This is on top of, not instead of, the per-event vz/cluster-pT
// Data/MC reweighting (Reweighter, same as unfolder.cc's mcWeight; always 1 for Data).

#include "../src/ana.h"
#include "../src/reweight_utility.h"
#include "../src/drawer.h"
#include "../src/unfold_utility.h"
R__LOAD_LIBRARY(libgammajet_unfold.so)

#include <vector>
#include <string>
#include <map>
#include <utility>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"

const int nCases = 6;
// Only cases 0-4 are directly built per event from the branches (nDirectCases below);
// case 5 (kPurityCorr) is derived from case 4's ("Everything new") own selection - see
// file header comment.
const int nDirectCases = 5;
enum CaseIdx { kAllOld = 0, kNewBdtIso = 1, kNewJesJer = 2, kWithSat = 3, kAllNew = 4, kPurityCorr = 5 };
const char * caseNames[nCases]  = {"AllOld", "NewBDT+ISO", "NewJES+JER", "WithPixelSat", "AllNew", "PurityCorrected"};
const char * caseLabels[nCases] = {"1. Everything old", "2. New BDT + new iso only",
                                    "3. New JES+JER only", "4. With pixel-saturation only",
                                    "5. Everything new", "6. Everything new, purity-corr. (MC=case 5)"};
const int caseColors[nCases] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 2, kOrange + 7};

const int ir = 2; // R=0.4, representative radius (see file header comment)

// truth-level cluster-pT trigger-stitching window, pythia sim, from
// treeuser.cc's threshmap[-1]/threshmap_high[-1] (Photon5/10/20 rows). MC only.
std::map<std::string, std::pair<float, float>> truthPtWindow = {
  {"Photon5",  {0.f,  12.f}},
  {"Photon10", {12.f, 24.f}},
  {"Photon20", {24.f, 100.f}}
};

// Per-sample cross-section/luminosity normalization - literal copy of drawer.h's
// scalemap[/*isphoton=*/1][sample] for sim=="pythia" (src/drawer.h's `drawer` ctor).
// MC only - Data is a single unified stream, not a pT-hat-sliced production.
std::map<std::string, double> sampleWeight = {
  {"Photon5",  146359.3},
  {"Photon10", 6944.675},
  {"Photon20", 130.4461}
};

struct Accum {
  double sumw = 0, sumw2 = 0, sumwx = 0, sumwx2 = 0;
  Long64_t n = 0;
  void fill(double x, double w) {
    sumw += w; sumw2 += w * w; sumwx += w * x; sumwx2 += w * x * x; n++;
  }
  double mean() const { return sumw > 0 ? sumwx / sumw : 0; }
  double meanErr() const {
    if (sumw <= 0 || n < 2) return 0;
    double var = sumwx2 / sumw - mean() * mean();
    if (var < 0) var = 0;
    // effective N for weighted stat error, standard sumw^2/sumw2
    double neff = sumw2 > 0 ? sumw * sumw / sumw2 : n;
    return std::sqrt(var / neff);
  }
};

// One case's worth of histograms/accumulators for one sample (MC or Data) - see
// makeResultSet() below.
struct ResultSet {
  TH1D * hxjIncl[nCases];
  std::vector<std::vector<TH1D*>> hxjPt; // [case][reported pT bin]
  Accum accIncl[nCases];
  std::vector<std::vector<Accum>> accPt; // [case][reported pT bin]
};

ResultSet makeResultSet(const char * tag, int nb, const double * edges) {
  ResultSet r;
  r.hxjPt.assign(nCases, std::vector<TH1D*>(ana::nPtBinsUsed));
  r.accPt.assign(nCases, std::vector<Accum>(ana::nPtBinsUsed));
  for (int ic = 0; ic < nCases; ic++) {
    r.hxjIncl[ic] = new TH1D(Form("hxjIncl_%s_%s", tag, caseNames[ic]), ";x_{J#gamma};counts", nb, edges);
    r.hxjIncl[ic]->SetStatOverflows(TH1::kConsider);
    r.hxjIncl[ic]->SetLineColor(caseColors[ic]);
    r.hxjIncl[ic]->SetLineWidth(2);
    for (int ip = 0; ip < ana::nPtBinsUsed; ip++) {
      r.hxjPt[ic][ip] = new TH1D(Form("hxjPt_%s_%s_%d", tag, caseNames[ic], ip), ";x_{J#gamma};counts", nb, edges);
      r.hxjPt[ic][ip]->SetStatOverflows(TH1::kConsider);
      r.hxjPt[ic][ip]->SetLineColor(caseColors[ic]);
      r.hxjPt[ic][ip]->SetLineWidth(2);
    }
  }
  return r;
}

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

// trigger: one of "Photon5"/"Photon10"/"Photon20" when isMC, or "Data" when not (only
// used to pick the input filename/log text for Data - there's no per-trigger split).
void processFile(const std::string & trigger, bool isMC, TH1D * hxjIncl[nCases],
    std::vector<std::vector<TH1D*>> & hxjPt, Accum accIncl[nCases],
    std::vector<std::vector<Accum>> & accPt, Reweighter & rw, Long64_t maxEntries = -1) {

  std::string fname = isMC ?
      ana::path("trees/gammajet_pythia_") + trigger + ".root" :
      ana::path("trees/gammajet_Data.root");
  TFile * f = TFile::Open(fname.c_str(), "read");
  if (!f || f->IsZombie()) { std::cout << "Could not open " << fname << std::endl; return; }
  TTree * t = (TTree*)f->Get("towerntup");
  Long64_t nentries = t->GetEntries();
  if (maxEntries >= 0 && maxEntries < nentries) nentries = maxEntries;
  std::cout << "Processing " << fname << " (" << nentries << " entries, "
            << (isMC ? "MC" : "Data") << ")" << std::endl;

  t->SetBranchStatus("*", 0);
  std::vector<std::string> branches = {"vz", "cluster_pt", "cluster_eta", "cluster_phi",
       "cluster_showershape", "cluster_bdt_scores", "cluster_bdt_score_old",
       "jet_pt_old", "jet_eta", "jet_phi"};
  if (isMC) {
    branches.push_back("cluster_pt_nosat");
    branches.push_back("cluster_eta_nosat");
    branches.push_back("cluster_phi_nosat");
    branches.push_back("truth_cluster_pt");
    branches.push_back("jet_pt_smear_truth");
  } else {
    // no _nosat branches exist for Data - see file header comment.
    branches.push_back("cluster_pt_old");
    branches.push_back("cluster_eta_old");
    branches.push_back("cluster_phi_old");
    branches.push_back("jet_pt_calib");
  }
  for (const auto & bn : branches) t->SetBranchStatus(bn.c_str(), 1);

  Float_t vz, cluster_pt, cluster_eta, cluster_phi;
  Float_t cluster_pt_nosat = 0, cluster_eta_nosat = 0, cluster_phi_nosat = 0;
  Float_t cluster_pt_old = 0, cluster_eta_old = 0, cluster_phi_old = 0;
  Float_t cluster_showershape[12], cluster_bdt_scores[11], cluster_bdt_score_old;
  Float_t truth_cluster_pt = 0, jet_pt_old;
  Float_t jet_pt_smear_truth[ana::nJetR] = {0}, jet_pt_calib[ana::nJetR] = {0};
  Float_t jet_eta[ana::nJetR], jet_phi[ana::nJetR];

  t->SetBranchAddress("vz", &vz);
  t->SetBranchAddress("cluster_pt", &cluster_pt);
  t->SetBranchAddress("cluster_eta", &cluster_eta);
  t->SetBranchAddress("cluster_phi", &cluster_phi);
  t->SetBranchAddress("cluster_showershape", cluster_showershape);
  t->SetBranchAddress("cluster_bdt_scores", cluster_bdt_scores);
  t->SetBranchAddress("cluster_bdt_score_old", &cluster_bdt_score_old);
  t->SetBranchAddress("jet_pt_old", &jet_pt_old);
  t->SetBranchAddress("jet_eta", jet_eta);
  t->SetBranchAddress("jet_phi", jet_phi);
  if (isMC) {
    t->SetBranchAddress("cluster_pt_nosat", &cluster_pt_nosat);
    t->SetBranchAddress("cluster_eta_nosat", &cluster_eta_nosat);
    t->SetBranchAddress("cluster_phi_nosat", &cluster_phi_nosat);
    t->SetBranchAddress("truth_cluster_pt", &truth_cluster_pt);
    t->SetBranchAddress("jet_pt_smear_truth", jet_pt_smear_truth);
  } else {
    t->SetBranchAddress("cluster_pt_old", &cluster_pt_old);
    t->SetBranchAddress("cluster_eta_old", &cluster_eta_old);
    t->SetBranchAddress("cluster_phi_old", &cluster_phi_old);
    t->SetBranchAddress("jet_pt_calib", jet_pt_calib);
  }

  auto windowIt = truthPtWindow.find(trigger);
  bool haveWindow = isMC && windowIt != truthPtWindow.end();
  double xsecWeight = isMC ? sampleWeight.at(trigger) : 1.0;

  // Region A/Region C bookkeeping for case 6's Data purity correction (see file header
  // comment) - unused for MC, so only allocated for the (single) Data call to avoid
  // TH1D name collisions across processFile's repeated MC-trigger calls. Binning is
  // cloned from an existing per-pT-bin histogram purely to reuse ana::unfoldXjBins.
  Accum accA_incl, accC_incl;
  std::vector<Accum> accA_pt(ana::nPtBinsUsed), accC_pt(ana::nPtBinsUsed);
  TH1D * hA_incl = nullptr, * hC_incl = nullptr;
  std::vector<TH1D*> hA_pt(ana::nPtBinsUsed, nullptr), hC_pt(ana::nPtBinsUsed, nullptr);
  if (!isMC) {
    hA_incl = (TH1D*)hxjIncl[kAllOld]->Clone("hPurityRegionA_incl");
    hA_incl->Reset();
    hC_incl = (TH1D*)hxjIncl[kAllOld]->Clone("hPurityRegionC_incl");
    hC_incl->Reset();
    for (int ip = 0; ip < ana::nPtBinsUsed; ip++) {
      hA_pt[ip] = (TH1D*)hxjPt[kAllOld][ip]->Clone(Form("hPurityRegionA_pt%d", ip));
      hA_pt[ip]->Reset();
      hC_pt[ip] = (TH1D*)hxjPt[kAllOld][ip]->Clone(Form("hPurityRegionC_pt%d", ip));
      hC_pt[ip]->Reset();
    }
  }

  for (Long64_t e = 0; e < nentries; e++) {
    t->GetEntry(e);
    if (e % 2000000 == 0)
      std::cout << "  entry " << e << "/" << nentries
                << " (" << (float)e / nentries * 100. << "%)" << std::endl;

    if (fabs(vz) > ana::vzcut) continue;
    if (haveWindow && !(truth_cluster_pt > windowIt->second.first && truth_cluster_pt < windowIt->second.second))
      continue;

    float mcWeight = (isMC ? rw.GetWeight(vz, cluster_pt) : 1.0f) * xsecWeight;

    float jetEta = jet_eta[ir];
    float jetPhi = jet_phi[ir];
    // "new" jet pt for cases 3/5 - see the file header comment for the isMC/Data split.
    float newJetPt = isMC ? jet_pt_smear_truth[ir] : jet_pt_calib[ir];

    // "old kinematics" for cases 1/3 (paired with the old BDT/iso) and case 2 (paired with
    // the new BDT/iso): MC uses cluster_pt_nosat for all three, unchanged from before. Data
    // has no _nosat branch (see file header comment), so cases 1/3 fall back to
    // cluster_pt_old and case 2 falls back to plain cluster_pt - both locked down with the
    // user.
    float case13Pt  = isMC ? cluster_pt_nosat  : cluster_pt_old;
    float case13Eta = isMC ? cluster_eta_nosat : cluster_eta_old;
    float case13Phi = isMC ? cluster_phi_nosat : cluster_phi_old;
    float case2Pt   = isMC ? cluster_pt_nosat  : cluster_pt;
    float case2Eta  = isMC ? cluster_eta_nosat : cluster_eta;
    float case2Phi  = isMC ? cluster_phi_nosat : cluster_phi;

    // {kAllOld, kNewBdtIso, kNewJesJer, kWithSat, kAllNew}
    float phoPtArr[nDirectCases]  = {case13Pt, case2Pt, case13Pt, cluster_pt, cluster_pt};
    float phoEtaArr[nDirectCases] = {case13Eta, case2Eta, case13Eta, cluster_eta, cluster_eta};
    float phoPhiArr[nDirectCases] = {case13Phi, case2Phi, case13Phi, cluster_phi, cluster_phi};
    float isoArr[nDirectCases]    = {cluster_showershape[9], cluster_showershape[11],
                                cluster_showershape[9], cluster_showershape[9], cluster_showershape[11]};
    float bdtArr[nDirectCases]    = {cluster_bdt_score_old, cluster_bdt_scores[9],
                                cluster_bdt_score_old, cluster_bdt_score_old, cluster_bdt_scores[9]};
    float jetPtArr[nDirectCases]  = {jet_pt_old, jet_pt_old, newJetPt, jet_pt_old, newJetPt};

    for (int ic = 0; ic < nDirectCases; ic++) {
      float phoPt = phoPtArr[ic], phoEta = phoEtaArr[ic], phoPhi = phoPhiArr[ic];
      float jetPt = jetPtArr[ic];
      if (!std::isfinite(phoPt) || phoPt <= 0 || !std::isfinite(jetPt)) continue;

      if (!(phoPt >= ana::ptBins[0] && phoPt < ana::ptBins[ana::nPtBins] && jetPt > ana::jet_calib_pt_cut[ir]))
        continue;
      float dphi = deltaPhi(phoPhi, jetPhi);
      if (!checkPair(jetPt, jetEta, phoPt, phoEta, dphi)) continue;

      int iabcd = ana::findabcdBin(isoArr[ic], bdtArr[ic], 0);

      // Region-C bookkeeping for case 6's Data purity correction (see file header
      // comment) - reuses case 5's ("Everything new") kinematics/cuts just computed
      // above. Combined into case 6 after the full event loop, once NA/NC (and hence
      // the purity coefficients) are known for every used pT bin.
      if (!isMC && ic == kAllNew && (iabcd == 0 || iabcd == 2)) {
        float xjPurity = jetPt / phoPt;
        if (std::isfinite(xjPurity)) {
          Accum & acc = (iabcd == 0) ? accA_incl : accC_incl;
          TH1D * h = (iabcd == 0) ? hA_incl : hC_incl;
          acc.fill(xjPurity, mcWeight);
          h->Fill(xjPurity, mcWeight);
          int ptbinPurity = ana::findPtBin(phoPt);
          int iusedPurity = ptbinPurity - ana::firstUsedPtBin;
          if (iusedPurity >= 0 && iusedPurity < ana::nPtBinsUsed) {
            (iabcd == 0 ? accA_pt[iusedPurity] : accC_pt[iusedPurity]).fill(xjPurity, mcWeight);
            (iabcd == 0 ? hA_pt[iusedPurity] : hC_pt[iusedPurity])->Fill(xjPurity, mcWeight);
          }
        }
      }

      if (iabcd != 0) continue;

      float xj = jetPt / phoPt;
      if (!std::isfinite(xj)) continue;

      accIncl[ic].fill(xj, mcWeight);
      hxjIncl[ic]->Fill(xj, mcWeight);

      int ptbin = ana::findPtBin(phoPt);
      // only the 3 "reported" bins (ana::firstUsedPtBin..+nPtBinsUsed) get their own plot
      int iused = ptbin - ana::firstUsedPtBin;
      if (iused >= 0 && iused < ana::nPtBinsUsed) {
        accPt[ic][iused].fill(xj, mcWeight);
        hxjPt[ic][iused]->Fill(xj, mcWeight);
      }

      // MC case 6 = literal duplicate of case 5 ("do nothing different", locked down
      // with the user) - only reached here inside case 5/kAllNew's own iabcd==0
      // (Region A) selection, i.e. exactly case 5's own passing events.
      if (isMC && ic == kAllNew) {
        accIncl[kPurityCorr].fill(xj, mcWeight);
        hxjIncl[kPurityCorr]->Fill(xj, mcWeight);
        if (iused >= 0 && iused < ana::nPtBinsUsed) {
          accPt[kPurityCorr][iused].fill(xj, mcWeight);
          hxjPt[kPurityCorr][iused]->Fill(xj, mcWeight);
        }
      }
    }
  }

  if (!isMC) {
    // Combine case 5's Region A/Region C bookkeeping into case 6's purity-corrected
    // result, per used pT bin - see file header comment. Mirrors
    // insitu_utility::computeCorrectedMeans/purityCorrectByPtBin exactly (same
    // unfold_utility::purityCorrectCoeffs/purityCorrect calls, same ana::getPurity/
    // getPurityC "nominal"/ir=2 inputs), just fed from this file's own event loop
    // instead of a pre-built insitutree.
    for (int ip = 0; ip < ana::nPtBinsUsed; ip++) {
      double ptlo = ana::ptBinsUsed[ip], pthi = ana::ptBinsUsed[ip + 1];
      float pA        = ana::getPurity(ptlo, pthi, "nominal", ir);
      float pAErrLow  = ana::getPurityErrorLow(ptlo, pthi, "nominal", ir);
      float pAErrHigh = ana::getPurityErrorHigh(ptlo, pthi, "nominal", ir);
      float pC        = ana::getPurityC(ptlo, pthi, "nominal", ir);
      float pCErrLow  = ana::getPurityCErrorLow(ptlo, pthi, "nominal", ir);
      float pCErrHigh = ana::getPurityCErrorHigh(ptlo, pthi, "nominal", ir);

      double NA = accA_pt[ip].sumw, NC = accC_pt[ip].sumw;
      accPt[kPurityCorr][ip].n = accA_pt[ip].n + accC_pt[ip].n; // raw Region-A+C event count feeding the correction, not directly comparable to the other cases' region-A-only N_pass
      float coeffA, coeffC;
      unfold_utility::purityCorrectCoeffs(pA, pC, NA, NC, coeffA, coeffC);
      double Ncorr = coeffA * NA - coeffC * NC;
      if (Ncorr > 0) {
        accPt[kPurityCorr][ip].sumw   = Ncorr;
        accPt[kPurityCorr][ip].sumw2  = Ncorr; // effective-N == Ncorr, matching insitu_utility::computeCorrectedMeans' err = sqrt(var/Ncorr)
        accPt[kPurityCorr][ip].sumwx  = coeffA * accA_pt[ip].sumwx  - coeffC * accC_pt[ip].sumwx;
        accPt[kPurityCorr][ip].sumwx2 = coeffA * accA_pt[ip].sumwx2 - coeffC * accC_pt[ip].sumwx2;
      }

      TH1D * hcorr = unfold_utility::purityCorrect(hA_pt[ip], hC_pt[ip], pA, pAErrLow, pAErrHigh,
          pC, pCErrLow, pCErrHigh, Form("hPurityCorrTmp_pt%d", ip));
      if (hcorr) {
        hxjPt[kPurityCorr][ip]->Add(hcorr);
        delete hcorr;
      }
      delete hA_pt[ip];
      delete hC_pt[ip];
    }
    // "Inclusive" row for case 6: sum of the 3 used-pT-bin corrected results (see file
    // header comment) - NOT ana::ptBins[0]-ptBins[nPtBins] like cases 1-5's inclusive
    // row, since purity is only measured per ana::ptBinsUsed bin.
    for (int ip = 0; ip < ana::nPtBinsUsed; ip++) {
      accIncl[kPurityCorr].sumw   += accPt[kPurityCorr][ip].sumw;
      accIncl[kPurityCorr].sumw2  += accPt[kPurityCorr][ip].sumw2;
      accIncl[kPurityCorr].sumwx  += accPt[kPurityCorr][ip].sumwx;
      accIncl[kPurityCorr].sumwx2 += accPt[kPurityCorr][ip].sumwx2;
      accIncl[kPurityCorr].n      += accPt[kPurityCorr][ip].n;
    }
    delete hA_incl;
    delete hC_incl;
  }

  f->Close();
}

// Writes the summary table + the mean-vs-pT/shape-comparison PDF for one ResultSet
// (either MC or Data) - the exact same plot set either way, just re-labeled.
void writeComparisonOutput(drawer & d, const std::string & jetFeature, const std::string & sampleLabel,
    ResultSet & rs, const std::string & pdfPath, const std::string & summaryPath) {

  std::ofstream out(summaryPath);
  auto report = [&](std::ostream & os) {
    os << "=== " << sampleLabel << ": Inclusive (ana::ptBins[0]-ptBins[nPtBins] = 13-100 GeV) ===\n";
    os << Form("%-28s %12s %14s %12s\n", "case", "N_pass", "sum(weight)", "<xJ> +/- err");
    for (int ic = 0; ic < nCases; ic++) {
      os << Form("%-28s %12lld %14.1f %6.4f +/- %6.4f\n", caseLabels[ic], rs.accIncl[ic].n,
                  rs.accIncl[ic].sumw, rs.accIncl[ic].mean(), rs.accIncl[ic].meanErr());
    }
    os << "\n=== " << sampleLabel << ": Per reported pT bin (ana::ptBinsUsed) ===\n";
    for (int ip = 0; ip < ana::nPtBinsUsed; ip++) {
      double lo = ana::ptBinsUsed[ip], hi = ana::ptBinsUsed[ip + 1];
      os << Form("-- %.0f < p_{T}^{gamma} < %.0f GeV --\n", lo, hi);
      os << Form("%-28s %12s %14s %12s\n", "case", "N_pass", "sum(weight)", "<xJ> +/- err");
      for (int ic = 0; ic < nCases; ic++) {
        os << Form("%-28s %12lld %14.1f %6.4f +/- %6.4f\n", caseLabels[ic], rs.accPt[ic][ip].n,
                    rs.accPt[ic][ip].sumw, rs.accPt[ic][ip].mean(), rs.accPt[ic][ip].meanErr());
      }
    }
    os << "\nNote: case 6's N_pass/sum(weight) are the purity-corrected effective yield "
          "(Region A minus the Region-C-derived background estimate), not a raw event "
          "count - for MC it's identical to case 5 by construction; for Data, its "
          "Inclusive row above sums only the 3 reported pT bins, not the full 13-100 GeV "
          "range (purity is only measured per ana::ptBinsUsed bin).\n";
  };
  report(std::cout);
  report(out);
  out.close();
  std::cout << "\nWrote " << sampleLabel << " summary to " << summaryPath << std::endl;

  TCanvas * c = new TCanvas(Form("c_%s", sampleLabel.c_str()), "", 700, 600);
  c->SetLeftMargin(.13);
  c->SetBottomMargin(.13);
  c->SaveAs((pdfPath + "[").c_str());

  auto drawOverlay = [&](std::vector<TH1D*> hs, std::vector<Accum*> accs, const char * title,
      const char * ptFeature) {
    c->Clear();
    double ymax = 0;
    for (auto h : hs) {
      // ana::unfoldXjBins is variable-width (0.1 wide up to 1.3, then 0.2, then 0.3 near
      // the tail) - Integral() here deliberately omits "width" (it's the total
      // event/weight count, the correct normalization denominator; under/overflow
      // included so the shown shape stays normalized to the true total), but Scale's
      // "width" option is required so bin content is a density, not a bare bin
      // fraction - otherwise the wider tail bins would be shown too tall relative to the
      // narrower ones for the same underlying number of events.
      double integ = h->Integral(0, h->GetNbinsX() + 1);
      if (integ > 0) h->Scale(1.0 / integ, "width");
      if (h->GetMaximum() > ymax) ymax = h->GetMaximum();
    }
    TLegend * leg = new TLegend(0.47, 0.60, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.028);
    for (int ic = 0; ic < nCases; ic++) {
      hs[ic]->SetTitle(title);
      hs[ic]->GetYaxis()->SetRangeUser(0, ymax * 1.35);
      hs[ic]->GetYaxis()->SetTitle("(1/N) dN/dx_{J#gamma}");
      hs[ic]->Draw(ic == 0 ? "hist" : "hist same");
      leg->AddEntry(hs[ic], Form("%s  #LTx_{J}#GT=%.3f#pm%.3f", caseLabels[ic],
                    accs[ic]->mean(), accs[ic]->meanErr()), "l");
    }
    leg->Draw();
    d.drawAll({sampleLabel}, {ptFeature, jetFeature}, .17, .87, 14, gPad->GetWh() * 0.8);
    c->SaveAs(pdfPath.c_str());
  };

  // ---------- mean xJ vs pT, all six cases together ----------
  c->Clear();
  TMultiGraph * mg = new TMultiGraph();
  TLegend * legMean = new TLegend(0.5, 0.65, 0.88, 0.88);
  legMean->SetBorderSize(0);
  legMean->SetFillStyle(0);
  // stagger points sharing a pT bin (GeV) so their error bars don't overlap
  double xOffsets[nCases] = {-0.75, -0.45, -0.15, 0.15, 0.45, 0.75};
  for (int ic = 0; ic < nCases; ic++) {
    TGraphErrors * g = new TGraphErrors(ana::nPtBinsUsed);
    for (int ip = 0; ip < ana::nPtBinsUsed; ip++) {
      double center = 0.5 * (ana::ptBinsUsed[ip] + ana::ptBinsUsed[ip + 1]);
      g->SetPoint(ip, center + xOffsets[ic], rs.accPt[ic][ip].mean());
      g->SetPointError(ip, 0, rs.accPt[ic][ip].meanErr());
    }
    g->SetLineColor(caseColors[ic]);
    g->SetMarkerColor(caseColors[ic]);
    g->SetMarkerStyle(20 + ic);
    g->SetLineWidth(2);
    mg->Add(g, "P");
    legMean->AddEntry(g, caseLabels[ic], "p");
  }
  mg->SetTitle(";p_{T}^{#gamma} [GeV];#LTx_{J#gamma}#GT");
  // range from the actual points (mean +/- err) rather than a fixed MC-tuned window -
  // Data's lower-stat 20-25 GeV bin can sit well below MC's typical means and got clipped
  // by the old hardcoded 0.78-0.9 range.
  double gmin = 1e9, gmax = -1e9;
  for (int ic = 0; ic < nCases; ic++) {
    for (int ip = 0; ip < ana::nPtBinsUsed; ip++) {
      if (rs.accPt[ic][ip].sumw <= 0) continue;
      double m = rs.accPt[ic][ip].mean(), e = rs.accPt[ic][ip].meanErr();
      gmin = std::min(gmin, m - e);
      gmax = std::max(gmax, m + e);
    }
  }
  if (gmin > gmax) { gmin = 0.78; gmax = 0.9; } // fallback: no points at all
  // asymmetric padding - legMean (y=0.65-0.88) and drawAll's label both sit at the top of
  // the pad, so give that side much more headroom than the bottom to keep points/error
  // bars from running into the text.
  double range = gmax - gmin;
  double padBottom = std::max(0.05 * range, 0.005);
  double padTop = std::max(0.45 * range, 0.02);
  mg->SetMinimum(gmin - padBottom);
  mg->SetMaximum(gmax + padTop);
  mg->Draw("A P");
  mg->GetXaxis()->SetLimits(ana::ptBinsUsed[0] - 2, ana::ptBinsUsed[ana::nPtBinsUsed] + 2);
  gPad->Modified();
  legMean->Draw();
  d.drawAll({sampleLabel}, {jetFeature}, .17, .87, 14, gPad->GetWh() * 0.8);
  c->SaveAs(pdfPath.c_str());

  for (int ip = 0; ip < ana::nPtBinsUsed; ip++) {
    drawOverlay({rs.hxjPt[0][ip], rs.hxjPt[1][ip], rs.hxjPt[2][ip], rs.hxjPt[3][ip], rs.hxjPt[4][ip], rs.hxjPt[5][ip]},
                {&rs.accPt[0][ip], &rs.accPt[1][ip], &rs.accPt[2][ip], &rs.accPt[3][ip], &rs.accPt[4][ip], &rs.accPt[5][ip]},
                ";x_{J#gamma};(1/N) dN/dx_{J#gamma}",
                Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", ana::ptBinsUsed[ip], ana::ptBinsUsed[ip + 1]));
  }

  c->SaveAs((pdfPath + "]").c_str());
  std::cout << "Wrote " << sampleLabel << " plots to " << pdfPath << std::endl;
}

// Data/MC ratio: mean-ratio-vs-pT only (one page), all six cases overlaid - per user
// request, dropped the earlier per-pT-bin normalized-shape ratio pages (so this now
// compares only the ratio of the means, not the full xJ shapes) and reports Data/MC
// rather than MC/Data.
void writeRatioOutput(drawer & d, const std::string & jetFeature, ResultSet & mc, ResultSet & data,
    const std::string & pdfPath) {

  TCanvas * c = new TCanvas("c_ratio", "", 700, 600);
  c->SetLeftMargin(.13);
  c->SetBottomMargin(.13);
  c->SaveAs((pdfPath + "[").c_str());

  // ---------- mean ratio vs pT (Data/MC) ----------
  c->Clear();
  TMultiGraph * mg = new TMultiGraph();
  TLegend * leg = new TLegend(0.5, 0.65, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  double xOffsets[nCases] = {-0.75, -0.45, -0.15, 0.15, 0.45, 0.75};
  std::vector<std::vector<double>> ratios(nCases, std::vector<double>(ana::nPtBinsUsed, 0)),
      ratioErrs(nCases, std::vector<double>(ana::nPtBinsUsed, 0));
  for (int ic = 0; ic < nCases; ic++) {
    TGraphErrors * g = new TGraphErrors(ana::nPtBinsUsed);
    for (int ip = 0; ip < ana::nPtBinsUsed; ip++) {
      double center = 0.5 * (ana::ptBinsUsed[ip] + ana::ptBinsUsed[ip + 1]);
      double mMC = mc.accPt[ic][ip].mean(), eMC = mc.accPt[ic][ip].meanErr();
      double mData = data.accPt[ic][ip].mean(), eData = data.accPt[ic][ip].meanErr();
      double ratio = (mMC > 0 && mData > 0) ? mData / mMC : 0;
      double ratioErr = (ratio > 0) ?
          ratio * std::sqrt(std::pow(eMC / mMC, 2) + std::pow(eData / mData, 2)) : 0;
      ratios[ic][ip] = ratio;
      ratioErrs[ic][ip] = ratioErr;
      g->SetPoint(ip, center + xOffsets[ic], ratio);
      g->SetPointError(ip, 0, ratioErr);
    }
    g->SetLineColor(caseColors[ic]);
    g->SetMarkerColor(caseColors[ic]);
    g->SetMarkerStyle(20 + ic);
    g->SetLineWidth(2);
    mg->Add(g, "P");
    leg->AddEntry(g, caseLabels[ic], "p");
  }
  mg->SetTitle(";p_{T}^{#gamma} [GeV];Data / MC  #LTx_{J#gamma}#GT");
  // range from the actual points, same rationale as writeComparisonOutput's mean-vs-pT
  // plot - a fixed 0.9-1.1 window could clip case 6's purity-corrected Data points.
  double gmin = 1e9, gmax = -1e9;
  for (int ic = 0; ic < nCases; ic++) {
    for (int ip = 0; ip < ana::nPtBinsUsed; ip++) {
      if (ratios[ic][ip] <= 0) continue;
      gmin = std::min(gmin, ratios[ic][ip] - ratioErrs[ic][ip]);
      gmax = std::max(gmax, ratios[ic][ip] + ratioErrs[ic][ip]);
    }
  }
  if (gmin > gmax) { gmin = 0.9; gmax = 1.1; } // fallback: no points at all
  // asymmetric padding - leg (y=0.65-0.88) and drawAll's label both sit at the top of the
  // pad, same rationale as writeComparisonOutput's mean-vs-pT plot above.
  double range = gmax - gmin;
  double padBottom = std::max(0.05 * range, 0.005);
  double padTop = std::max(0.45 * range, 0.02);
  mg->SetMinimum(gmin - padBottom);
  mg->SetMaximum(gmax + padTop);
  mg->Draw("A P");
  mg->GetXaxis()->SetLimits(ana::ptBinsUsed[0] - 2, ana::ptBinsUsed[ana::nPtBinsUsed] + 2);
  gPad->Modified();
  TLine * refLine1 = new TLine(ana::ptBinsUsed[0] - 2, 1.0, ana::ptBinsUsed[ana::nPtBinsUsed] + 2, 1.0);
  refLine1->SetLineStyle(2);
  refLine1->SetLineColor(kGray + 2);
  refLine1->Draw();
  leg->Draw();
  d.drawAll({"Data / MC"}, {jetFeature}, .17, .87, 14, gPad->GetWh() * 0.8);
  c->SaveAs(pdfPath.c_str());

  c->SaveAs((pdfPath + "]").c_str());
  std::cout << "Wrote Data/MC ratio plots to " << pdfPath << std::endl;
}

void compare_old_new(Long64_t maxEntriesPerFile = -1) {
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();

  Reweighter rw;

  // drawer::drawAll (used below for the sPHENIX label) only needs its member function,
  // not the per-sample unfolding-output files its constructor opens - some of those
  // (Jet8/20/30/50/60/80_pythia_nominal_unfolding.root) don't exist on this box, so
  // silence the resulting (harmless) TFile::Open error spam.
  int oldErrLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  drawer d("pythia", "nominal");
  gErrorIgnoreLevel = oldErrLevel;

  int nb = ana::nUnfoldXjBins;
  const double * edges = ana::unfoldXjBins;
  // jet-radius/pT-cut feature line for drawer::drawAll - fixed at the representative
  // ir=2/R=0.4 radius (see file header comment), same on every page/sample.
  std::string jetFeature = Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV", ana::JetRs[ir], ana::jet_calib_pt_cut[ir]);

  ResultSet mcSet = makeResultSet("MC", nb, edges);
  for (const std::string & trigger : {"Photon5", "Photon10", "Photon20"}) {
    processFile(trigger, /*isMC=*/true, mcSet.hxjIncl, mcSet.hxjPt, mcSet.accIncl, mcSet.accPt, rw, maxEntriesPerFile);
  }
  writeComparisonOutput(d, jetFeature, "Pythia8 #gamma+jet MC", mcSet,
      ana::path("temporary_study/compare_old_new.pdf"),
      ana::path("temporary_study/compare_old_new_summary.txt"));

  std::string dataFname = ana::path("trees/gammajet_Data.root");
  bool haveData = !gSystem->AccessPathName(dataFname.c_str()); // AccessPathName returns 0 (false) iff the file exists
  if (haveData) {
    ResultSet dataSet = makeResultSet("Data", nb, edges);
    processFile("Data", /*isMC=*/false, dataSet.hxjIncl, dataSet.hxjPt, dataSet.accIncl, dataSet.accPt, rw, maxEntriesPerFile);
    writeComparisonOutput(d, jetFeature, "p+p Run24 Data", dataSet,
        ana::path("temporary_study/compare_old_new_data.pdf"),
        ana::path("temporary_study/compare_old_new_data_summary.txt"));
    writeRatioOutput(d, jetFeature, mcSet, dataSet,
        ana::path("temporary_study/compare_old_new_ratio.pdf"));
  } else {
    std::cout << "\nData tree not found at " << dataFname
              << " - skipping Data processing and MC/Data ratio plots."
              << " Re-run this macro once it's available." << std::endl;
  }
}
