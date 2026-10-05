// Print truth photons passing the full selection in the [25,35) GeV truth bin. Uses treeuser
// directly because unfolder's constructor recreates the in-situ output file.
// check_pair is a copy of unfolder::check_pair.

#include "../src/treeuser.h"
#include "../src/pho_object.h"
#include "../src/jet_object.h"
R__LOAD_LIBRARY(libgammajet_unfold.so);

bool check_pair(jet_object jet, int ir, pho_object pho) {
  float dphi = jet.deltaPhi(pho);
  int ptbin = ana::findPtBin(pho.pt);
  if (ptbin == -1) return false;

  float val = jet.pt / pho.pt;
  float lowval = ana::jet_calib_pt_cut[ir] / ana::ptBins[ptbin];
  float lowbin = ana::unfoldXjBins[ana::findUnfoldXjBin(lowval) + 1];

  if (val < lowbin) return false;
  if (fabs(pho.eta) > ana::photonEtaCut) return false;
  if (fabs(jet.eta) > ana::etacut - ana::JetRs[ir]) return false;
  if (dphi < ana::oppcut) return false;
  return true;
}

void print_truth_photon_pt25_35(string trigger = "Photon20", string sim = "pythia") {
  const int ir = 2; // nominal R=0.4
  const int target_ptbin = 3; // ana::ptBins = {13,15,20,25,35,100} -> index 3 == [25,35)

  treeuser tu(trigger, sim);
  if (!tu.isMC) {
    cout << "trigger=" << trigger << " is Data - no truth branches, nothing to print." << endl;
    return;
  }

  Long64_t nentries = tu.t->GetEntriesFast();
  int nfound = 0;

  for (Long64_t e = 0; e < nentries; e++) {
    tu.t->GetEntry(e);

    if (fabs(tu.vz) > ana::vzcut) continue;

    vector<bool> keepMC = tu.check_keep_MC(tu.truth_cluster_pt, tu.cluster_pt, tu.truth_jet_pt, tu.jet_pt_smear_reco, trigger);
    if (!keepMC.at(keepMC.size() - 1)) continue;
    if (!keepMC[ir]) continue;

    pho_object maxpho_truth(
        tu.truth_cluster_pt,
        tu.truth_cluster_e,
        tu.truth_cluster_eta,
        tu.truth_cluster_phi,
        tu.truth_cluster_iso3,
        tu.truth_cluster_iso4,
        0,    // no time object for truth
        0.99, // truth photon is guaranteed a photon
        2     // truth photon is guaranteed a photon
    );
    jet_object maxjet_truth(
        tu.truth_jet_pt[ir],
        tu.truth_jet_e[ir],
        tu.truth_jet_eta[ir],
        tu.truth_jet_phi[ir],
        0, 0, 0, 0
    );

    bool ispaired_truth = false;
    if (maxpho_truth.pt >= ana::ptBins[0] && maxpho_truth.pt < ana::ptBins[ana::nPtBins]
        && maxjet_truth.pt > ana::jet_calib_pt_cut[ir]) {
      ispaired_truth = check_pair(maxjet_truth, ir, maxpho_truth);
    }
    int iabcd_truth = ana::findabcdBin(maxpho_truth.iso4, maxpho_truth.bdt, 0);
    ispaired_truth = ispaired_truth && (iabcd_truth == 0);

    if (!ispaired_truth) continue;
    if (ana::findPtBin(maxpho_truth.pt) != target_ptbin) continue;

    if (nfound < 1000) {
    cout << nfound << ": "
         << "  pt=" << maxpho_truth.pt
         << "  eta=" << maxpho_truth.eta
         << "  phi=" << maxpho_truth.phi
         << "  | jet pt=" << maxjet_truth.pt
         << "  e=" << maxjet_truth.e
         << "  eta=" << maxjet_truth.eta
         << "  phi=" << maxjet_truth.phi
         << "  | z=" << tu.vz
         << endl;
    }
    else break;
    nfound++;
  }

  cout << "Found " << nfound << " truth photons passing all cuts in the ["
       << ana::ptBins[target_ptbin] << "," << ana::ptBins[target_ptbin + 1]
       << ") GeV pT bin, jet R=" << ana::JetRs[ir] << endl;
}
