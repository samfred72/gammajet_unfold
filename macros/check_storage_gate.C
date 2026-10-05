// Check whether sam's tree-storage-gate rejections (nikhil's "CLUSTER leg failed") are
// actually explained by the reco cluster pt failing, while the truth cluster pt clears
// the same threshold via Sam's OR-based accept logic (reco OR truth pt > cut).
#include "../src/treeuser.h"
R__LOAD_LIBRARY(libgammajet_unfold.so);

void check_storage_gate(string trigger = "Photon20", string sim = "pythia") {
  struct Target { float pt, eta, phi; };
  vector<Target> targets = {
    {25.1907f,  0.994574f, -0.0400101f}, // sam#8
    {28.8779f,  0.149523f, -1.10169f},   // sam#48
    {27.4255f, -0.816893f, -2.90489f},   // sam#57
    {25.4376f,  0.90752f,   2.30416f},   // sam#59
    {33.0944f, -0.167785f,  0.758571f},  // sam#79
    {25.7458f,  0.345114f, -2.57063f},   // sam#72 (zvtx case)
  };

  treeuser tu(trigger, sim);
  Long64_t nentries = tu.t->GetEntriesFast();
  const float TOL = 0.001;
  int nleft = targets.size();

  for (Long64_t e = 0; e < nentries && nleft > 0; e++) {
    tu.t->GetEntry(e);
    for (auto & t : targets) {
      if (t.pt == 0) continue;
      if (fabs(tu.truth_cluster_pt - t.pt) < TOL &&
          fabs(tu.truth_cluster_eta - t.eta) < TOL &&
          fabs(tu.truth_cluster_phi - t.phi) < TOL) {
        cout << "truth_cluster_pt=" << tu.truth_cluster_pt
             << "  reco cluster_pt=" << tu.cluster_pt
             << "  vz=" << tu.vz
             << "  jet_pt_calib[2]=" << tu.jet_pt_calib[2]
             << endl;
        t.pt = 0; // mark done
        nleft--;
        break;
      }
    }
  }
}
