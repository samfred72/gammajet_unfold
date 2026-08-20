// Diagnostic: for each truth photon in nikhil_only_targets.csv (pt/eta/phi of truth
// photons that nikhil's list keeps but sam.txt / print_truth_photon_pt25_35.C does not),
// walk the exact same cut sequence as print_truth_photon_pt25_35.C and report the first
// cut that removes it.
//
// Cut chain mirrors print_truth_photon_pt25_35.C (macros/print_truth_photon_pt25_35.C)
// byte-for-byte: vz cut -> check_keep_MC trigger window -> jet_calib_pt_cut precondition
// -> check_pair (xJ floor, photon |eta|, jet |eta|, back-to-back dphi) -> ABCD isolation
// region -> target pT bin.

#include "/home/samson72/sphnx/gammajet_unfold/src/treeuser.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/pho_object.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/jet_object.h"
#include <fstream>
#include <sstream>
R__LOAD_LIBRARY(libgammajet_unfold.so);

struct Target {
  int idx;
  float pt, eta, phi;
  bool found = false;
};

void diagnose_nikhil_only(string trigger = "Photon20", string sim = "pythia") {
  const int ir = 2; // R=0.4, same as print_truth_photon_pt25_35.C
  const int target_ptbin = 3;

  // Load targets
  vector<Target> targets;
  ifstream fin("/tmp/claude-1000/-home-samson72-sphnx/96a98a79-a3d7-4e7f-95dc-e5cdaaefefee/scratchpad/nikhil_only_targets.csv");
  string line;
  while (getline(fin, line)) {
    stringstream ss(line);
    string tok;
    Target t;
    getline(ss, tok, ','); t.idx = stoi(tok);
    getline(ss, tok, ','); t.pt = stof(tok);
    getline(ss, tok, ','); t.eta = stof(tok);
    getline(ss, tok, ','); t.phi = stof(tok);
    targets.push_back(t);
  }
  cout << "Loaded " << targets.size() << " target photons" << endl;

  treeuser tu(trigger, sim);
  Long64_t nentries = tu.t->GetEntriesFast();

  const float TOL = 0.001;

  for (Long64_t e = 0; e < nentries; e++) {
    tu.t->GetEntry(e);

    // Quick pre-filter: is this event's truth cluster near ANY target?
    int match = -1;
    for (size_t k = 0; k < targets.size(); k++) {
      if (targets[k].found) continue;
      if (fabs(tu.truth_cluster_pt - targets[k].pt) < TOL &&
          fabs(tu.truth_cluster_eta - targets[k].eta) < TOL &&
          fabs(tu.truth_cluster_phi - targets[k].phi) < TOL) {
        match = k;
        break;
      }
    }
    if (match < 0) continue;

    targets[match].found = true;
    cout << "\n=== nik#" << targets[match].idx
         << "  pt=" << tu.truth_cluster_pt
         << " eta=" << tu.truth_cluster_eta
         << " phi=" << tu.truth_cluster_phi
         << "  (event " << e << ") ===" << endl;

    // --- Step 1: vz cut ---
    if (fabs(tu.vz) > ana::vzcut) {
      cout << "  CUT: |vz| = " << fabs(tu.vz) << " > vzcut(" << ana::vzcut << ")" << endl;
      continue;
    }

    // --- Step 2: check_keep_MC trigger truth-pt window ---
    vector<bool> keepMC = tu.check_keep_MC(tu.truth_cluster_pt, tu.cluster_pt, tu.truth_jet_pt, tu.jet_pt_smear, trigger);
    if (!keepMC.at(keepMC.size() - 1)) {
      cout << "  CUT: truth_cluster_pt=" << tu.truth_cluster_pt
           << " fails trigger window (Photon20/pythia: 24 < pt < 100)" << endl;
      continue;
    }
    if (!keepMC[ir]) {
      cout << "  CUT: keepMC[" << ir << "] false (jet trigger window - should be trivial=1 for photon trigger)" << endl;
      continue;
    }

    // --- Build truth photon/jet objects exactly as in print_truth_photon_pt25_35.C ---
    pho_object maxpho_truth(
        tu.truth_cluster_pt, tu.truth_cluster_e, tu.truth_cluster_eta, tu.truth_cluster_phi,
        tu.truth_cluster_iso3, tu.truth_cluster_iso4, 0, 0.99, 2);
    jet_object maxjet_truth(
        tu.truth_jet_pt[ir], tu.truth_jet_e[ir], tu.truth_jet_eta[ir], tu.truth_jet_phi[ir],
        0, 0, 0, 0);

    cout << "  truth jet[R=" << ana::JetRs[ir] << "]: pt=" << maxjet_truth.pt
         << " eta=" << maxjet_truth.eta << " phi=" << maxjet_truth.phi << endl;
    cout << "  truth photon iso4=" << maxpho_truth.iso4 << endl;

    // --- Step 3: pairing precondition ---
    bool precond = (maxpho_truth.pt >= ana::ptBins[0] && maxpho_truth.pt < ana::ptBins[ana::nPtBins]
                     && maxjet_truth.pt > ana::jet_calib_pt_cut[ir]);
    bool ispaired_truth = false;
    if (!precond) {
      if (!(maxpho_truth.pt >= ana::ptBins[0] && maxpho_truth.pt < ana::ptBins[ana::nPtBins])) {
        cout << "  CUT: maxpho_truth.pt=" << maxpho_truth.pt << " outside ptBins[0.." << ana::nPtBins << "] range" << endl;
      }
      if (!(maxjet_truth.pt > ana::jet_calib_pt_cut[ir])) {
        cout << "  CUT: truth jet pt=" << maxjet_truth.pt
             << " <= jet_calib_pt_cut[" << ir << "]=" << ana::jet_calib_pt_cut[ir]
             << "  (no matching truth jet of adequate pT was found/paired for this photon)" << endl;
      }
    } else {
      // --- Step 4: check_pair ---
      float dphi_val = maxjet_truth.deltaPhi(maxpho_truth);
      int ptbin = ana::findPtBin(maxpho_truth.pt);
      if (ptbin == -1) {
        cout << "  CUT: check_pair - findPtBin(" << maxpho_truth.pt << ") == -1" << endl;
      } else {
        float val = maxjet_truth.pt / maxpho_truth.pt;
        float lowval = ana::jet_calib_pt_cut[ir] / ana::ptBins[ptbin];
        float lowbin = ana::unfoldXjBins[ana::findUnfoldXjBin(lowval) + 1];
        bool pass = true;
        if (val < lowbin) {
          cout << "  CUT: check_pair - xJ = jetpt/phopt = " << val
               << " < lowbin(" << lowbin << ")  [xJ unfolding floor]" << endl;
          pass = false;
        }
        if (fabs(maxpho_truth.eta) > ana::etacut) {
          cout << "  CUT: check_pair - |photon eta| = " << fabs(maxpho_truth.eta)
               << " > etacut(" << ana::etacut << ")" << endl;
          pass = false;
        }
        if (fabs(maxjet_truth.eta) > ana::etacut - ana::JetRs[ir]) {
          cout << "  CUT: check_pair - |jet eta| = " << fabs(maxjet_truth.eta)
               << " > etacut-R(" << (ana::etacut - ana::JetRs[ir]) << ")" << endl;
          pass = false;
        }
        if (dphi_val < ana::oppcut) {
          cout << "  CUT: check_pair - dphi(jet,photon) = " << dphi_val
               << " < oppcut(" << ana::oppcut << ")  [not back-to-back]" << endl;
          pass = false;
        }
        ispaired_truth = pass;
      }
    }

    // --- Step 5: ABCD isolation region ---
    int iabcd_truth = ana::findabcdBin(maxpho_truth.iso4, maxpho_truth.bdt, 0);
    bool ispaired_final = ispaired_truth && (iabcd_truth == 0);
    if (ispaired_truth && iabcd_truth != 0) {
      cout << "  CUT: ABCD region - iso4=" << maxpho_truth.iso4
           << " -> findabcdBin=" << iabcd_truth << " != 0 (fails signal-region isolation cut, isoBins[0]=2)" << endl;
    }

    if (!ispaired_final) continue;

    // --- Step 6: target pt bin ---
    if (ana::findPtBin(maxpho_truth.pt) != target_ptbin) {
      cout << "  CUT: findPtBin(" << maxpho_truth.pt << ") != target_ptbin(" << target_ptbin << ")" << endl;
      continue;
    }

    cout << "  PASSES ALL CUTS - should have appeared in sam.txt (unexpected!)" << endl;
  }

  cout << "\n\n=== SUMMARY ===" << endl;
  int nnotfound = 0;
  for (auto & t : targets) {
    if (!t.found) {
      cout << "nik#" << t.idx << " (pt=" << t.pt << " eta=" << t.eta << " phi=" << t.phi
           << ") -- NO MATCHING TREE ENTRY FOUND (pt/eta/phi didn't match any truth_cluster in the tree within tol)" << endl;
      nnotfound++;
    }
  }
  cout << "Matched to tree entries: " << (targets.size() - nnotfound) << "/" << targets.size() << endl;
}
