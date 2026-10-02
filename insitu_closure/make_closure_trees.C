#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/treeuser.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/pho_object.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/jet_object.h"
#include <cmath>
#include <map>
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"
#include "TRandom3.h"
using namespace std;

// ana::findPtBin/findabcdBin/etc. are implemented in ana.cc, compiled into
// libgammajet_unfold.so - load it explicitly so cling resolves those symbols against
// the real compiled definitions instead of misbinding against the sibling gammajet
// project's own ana/drawer classes on the same library path (see insitu/grid_insitu.C).
R__LOAD_LIBRARY(libgammajet_unfold.so);

// In-situ JES closure test, stage 1: build a pair of insitu-style ("pho_pt","jet_pt",
// "abcd","weight","ir") trees for ONE MC sample (trigger, e.g. "Photon20"), split
// event-by-event into two random halves - same `rand.Integer(2) % 2` coin flip
// unfolder.cc's fill_matrix() uses to split MC into response-matrix-training vs.
// closure-test halves (src/unfolder.cc). One half plays "Data" (jet_pt scaled by
// `injectedScale`) and the other plays the fixed MC reference ("Sim", left at nominal
// reco scale). Each event is weighted by this sample's own cross-section weight
// (photon_scale below, matching insitu/grid_insitu.C's photon_scale map /
// drawer.h's scalemap[isphoton=1][sim="pythia"]) so fit_closure.C can combine this
// sample's output with Photon5/10/20's other two the same cross-section-weighted way
// grid_insitu.C's referenceMeans/buildMCXjByPtBin combine them for the real MC
// reference - the actual "stitching" step, done at read time in fit_closure.C rather
// than by this macro trying to process all three samples itself.
//
// One sample per macro call (rather than looping Photon5+10+20 internally) so
// insitu_closure/run_closure.sh can launch the three samples as separate, parallel ROOT
// processes - looping the ~9.8M-entry Photon20 and ~8.7M-entry Photon10 raw ntuples
// sequentially in one process is the dominant cost of this closure test, and they're
// otherwise fully independent (this macro touches no shared state besides its own
// output files, one pair per sample). `injectedScale` MUST then be the same value
// passed to every sample's invocation - see run_closure.sh, which draws it once and
// passes it to all three - since the whole point is one single, known JES
// miscalibration to recover, not three different ones.
//
// Deliberately does NOT reuse the `unfolder` class: its constructor immediately opens
// per-systag insitu_tree output files in RECREATE mode as a side effect (see
// unfolder.h), which would collide with/clobber the real production insitu/inputs/
// files. `treeuser` (unfolder's base class) has no such side effect - it only opens the
// raw ntuple and sets branch addresses - so it's instantiated directly here instead.

const char * closure_input_dir = "/home/samson72/sphnx/gammajet_unfold/insitu_closure/inputs";

// Cross-section weights for stitching the Photon5/10/20 pythia MC samples - same
// numbers as insitu/grid_insitu.C's photon_scale map / drawer.h's
// scalemap[isphoton=1][sample] for sim="pythia". Only pythia has a Photon5 sample at
// all (see treeuser.h's threshmap - herwig substitutes Photon20-only thresholds in its
// [-1] entry), so this closure test is pythia-only, matching the primary in-situ study.
// A `trigger` not in this map (e.g. a future non-photon cross-check sample) falls back
// to weight 1.0 with a printed warning, rather than failing outright.
map<string,double> photon_scale = {{"Photon5",146359.3},{"Photon10",6944.675},{"Photon20",130.4461}};

// Same pairing/low-xJ-floor cut as unfolder::check_pair (src/unfolder.cc), minus the
// unused `isreco` parameter - reimplemented here rather than pulled from an `unfolder`
// instance (see the file-header comment above for why) since the formula itself has no
// class-state dependency, only ana:: statics. testPt/floorScale follow check_pair's own
// convention: testPt is the "as measured" jet pt to test the floor against (here, the
// tree's own value - already scaled for the closure "data" half), floorScale loosens
// the floor around jet.pt's nominal lowbin edge so marginal events aren't lost before
// fit_closure.C's own pa scan gets to re-evaluate them at each trial pa.
bool closure_check_pair(jet_object & jet, int ir, pho_object & pho, float testPt, float floorScale) {
  float dphi = jet.deltaPhi(pho);
  int ptbin = ana::findPtBin(pho.pt);
  if (ptbin == -1) return false;
  float val = testPt/pho.pt;
  float lowval = ana::jet_calib_pt_cut[ir]/ana::ptBins[ptbin];
  float lowbin = ana::unfoldXjBins[ana::findUnfoldXjBin(lowval)+1];
  if (val < lowbin*floorScale) return false;
  if (fabs(pho.eta) > ana::photonEtaCut) return false;
  if (fabs(jet.eta) > ana::etacut - ana::JetRs[ir]) return false;
  if (dphi < ana::oppcut) return false;
  return true;
}

// trigger/sim: which raw MC ntuple to split (must be MC, e.g. "Photon20"/"pythia").
// injectedScale: the "data"-half jet-energy-mis-scale to apply - if < 0 (default), one
// is drawn uniformly from [scaleLow,scaleHigh) here for a convenient single-sample
// standalone run; for the real (parallel, multi-sample) closure test, pass an explicit
// value shared across all three samples' invocations - see run_closure.sh, which draws
// it once and passes it to every sample so they all get the SAME injected
// miscalibration instead of three different, uncombinable ones. seed=0 lets TRandom3
// auto-seed (a different event split each run); pass a fixed nonzero seed for a
// reproducible closure run (run_closure.sh uses distinct fixed seeds per sample so its
// three parallel processes don't share a random stream).
void make_closure_trees(string trigger = "Photon20", string sim = "pythia",
    float injectedScale = -1, float scaleLow = 0.9, float scaleHigh = 1.1, ULong_t seed = 0) {

  treeuser tu(trigger, sim);
  if (!tu.isMC) {
    cout << "make_closure_trees: trigger \"" << trigger << "\" is Data, not MC - the "
         << "closure test needs a truth-known MC sample to split. Aborting." << endl;
    return;
  }

  TRandom3 rand(seed);
  if (injectedScale < 0) {
    injectedScale = rand.Uniform(scaleLow, scaleHigh);
    cout << "No injectedScale given - drew one locally (fine for a single-sample run, "
         << "but NOT for combining with other samples - see run_closure.sh)." << endl;
  }
  cout << "Injected closure JES scale (\"data\"-half jet_pt *= " << injectedScale << ")" << endl;

  double sampleWeight = 1.0;
  auto it = photon_scale.find(trigger);
  if (it != photon_scale.end()) sampleWeight = it->second;
  else cout << "WARNING: no cross-section weight known for \"" << trigger << "\" - using 1.0." << endl;

  // Loose enough that both the widest possible injected scale (scaleHigh) and
  // fit_closure.C's own [0.80,1.20] scan window can still re-evaluate the low-xJ floor
  // correctly instead of losing marginal events to this tree-filling stage first - same
  // role insitu_utility::scanLow plays for unfolder.cc's ispairedInsitu, just widened
  // here since this closure test's injected/scan range is wider than the production
  // Data/MC JES gap that scanLow/scanHigh are tuned to.
  const float floorScale = 0.75;

  TFile * fdata = TFile::Open(Form("%s/ClosureData_%s_insitu.root", closure_input_dir, trigger.c_str()), "RECREATE");
  TTree * tdata = new TTree("insitutree", "closure-test \"data\" half (jet pt scaled by injected JES factor)");
  Float_t d_pho_pt, d_jet_pt, d_weight; Int_t d_abcd, d_ir;
  tdata->Branch("pho_pt", &d_pho_pt);
  tdata->Branch("jet_pt", &d_jet_pt);
  tdata->Branch("abcd", &d_abcd);
  tdata->Branch("weight", &d_weight);
  tdata->Branch("ir", &d_ir);

  TFile * fsim = TFile::Open(Form("%s/ClosureSim_%s_insitu.root", closure_input_dir, trigger.c_str()), "RECREATE");
  TTree * tsim = new TTree("insitutree", "closure-test \"sim\" half (nominal reco, fixed reference)");
  Float_t s_pho_pt, s_jet_pt, s_weight; Int_t s_abcd, s_ir;
  tsim->Branch("pho_pt", &s_pho_pt);
  tsim->Branch("jet_pt", &s_jet_pt);
  tsim->Branch("abcd", &s_abcd);
  tsim->Branch("weight", &s_weight);
  tsim->Branch("ir", &s_ir);

  Long64_t nentries = tu.t->GetEntriesFast();
  for (Long64_t e = 0; e < nentries; e++) {
    tu.t->GetEntry(e);
    if (e % 200000 == 0) cout << trigger << " entry " << e << "/" << nentries << " ("
        << (float)e/nentries*100. << "%)\t\r" << flush;

    if (fabs(tu.vz) > ana::vzcut) continue;

    // Per-sample pT-hat window (threshmap/threshmap_high, src/treeuser.h) - keeps the
    // three samples' truth-photon-pT ranges non-overlapping, so summing their
    // cross-section-weighted contributions in fit_closure.C is a real stitch, not
    // double-counting.
    vector<bool> keepMC = tu.check_keep_MC(tu.truth_cluster_pt, tu.cluster_pt, tu.truth_jet_pt, tu.jet_pt_smear_truth, trigger);
    if (!keepMC.at(keepMC.size()-1)) continue;

    // Same per-event coin flip unfolder.cc's fill_matrix() uses to split MC into its
    // two response-matrix-training/closure-test halves.
    bool isDataHalf = rand.Integer(2) % 2;

    // Nominal reco chain only (no JER/emscale/EMR systag variants) - but the nominal EM
    // resolution smearing (ana::emResolutionSigma) IS applied here even at "nominal", matching
    // every MC systag in unfolder.cc's fill_matrix() (systagEmrVariantArr defaults to nominal): it's how MC's
    // cluster-pt resolution is made to match Data's, not itself a systematic variation.
    float recoClusterPt = tu.cluster_pt + rand.Gaus(0, ana::emResolutionSigma(tu.truth_cluster_pt)*tu.truth_cluster_pt);
    pho_object maxpho(
        recoClusterPt, tu.cluster_e, tu.cluster_eta, tu.cluster_phi,
        tu.cluster_showershape[10], tu.cluster_showershape[11], tu.cluster_time,
        tu.cluster_bdt_scores[9],
        pho_object::get_showershape(tu.cluster_showershape, recoClusterPt));

    if (maxpho.pt < ana::ptBins[0] || maxpho.pt >= ana::ptBins[ana::nPtBins]) continue;

    for (int ir = 0; ir < ana::nJetR; ir++) {
      if (!keepMC[ir]) continue;

      float recoJetPt = tu.jet_pt_smear_truth[ir];
      jet_object maxjet(recoJetPt, tu.jet_e[ir], tu.jet_eta[ir], tu.jet_phi[ir], tu.jet_emfrac[ir], 0, 0, tu.jet_time[ir]);

      // The injected miscalibration - applied only to the "data" half, and only to what
      // gets written to the tree as the "measured" jet pt (the underlying reco chain
      // above, shared by both halves, is untouched) - mirrors unfolder.cc's Data-only
      // jesCorrectionArr convention (fill_matrix()'s recoJetPt branch), just with a
      // single known injected factor standing in for the real, unknown Data/MC gap.
      float measuredJetPt = isDataHalf ? recoJetPt*injectedScale : recoJetPt;

      if (!closure_check_pair(maxjet, ir, maxpho, measuredJetPt, floorScale)) continue;

      int iabcd = ana::findabcdBin(maxpho.iso4, maxpho.bdt, 0);
      if (iabcd == -1) continue;

      if (isDataHalf) {
        d_pho_pt = maxpho.pt; d_jet_pt = measuredJetPt; d_abcd = iabcd; d_ir = ir; d_weight = sampleWeight;
        tdata->Fill();
      } else {
        s_pho_pt = maxpho.pt; s_jet_pt = measuredJetPt; s_abcd = iabcd; s_ir = ir; s_weight = sampleWeight;
        tsim->Fill();
      }
    }
  }
  cout << endl;

  Long64_t nData = tdata->GetEntries(), nSim = tsim->GetEntries();

  fdata->cd();
  tdata->Write();
  TParameter<double> pScale("injectedScale", injectedScale);
  pScale.Write();
  fdata->Close();

  fsim->cd();
  tsim->Write();
  fsim->Close();

  cout << "Wrote " << closure_input_dir << "/ClosureData_" << trigger << "_insitu.root (" << nData << " rows)" << endl;
  cout << "Wrote " << closure_input_dir << "/ClosureSim_" << trigger << "_insitu.root (" << nSim << " rows)" << endl;
}
