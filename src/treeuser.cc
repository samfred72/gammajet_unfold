#include "treeuser.h"

void treeuser::treesetup() {
  // Set branch addresses and branch pointers
  cout << "Setting up tree " << t->GetEntries() << endl;
  if (!t) return;

  t->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
  t->SetBranchAddress("vz", &vz, &b_vz);
  if (!isMC) {
    t->SetBranchAddress("ScaledTriggerBit", ScaledTriggerBit, &b_ScaledTriggerBit);
    t->SetBranchAddress("LiveTriggerBit", LiveTriggerBit, &b_LiveTriggerBit);
    t->SetBranchAddress("Scaledowns", Scaledowns, &b_Scaledowns);
  }
  t->SetBranchAddress("mbd_time", &mbd_time, &b_mbd_time);

  t->SetBranchAddress("cluster_pt" , &cluster_pt , &b_cluster_pt );
  t->SetBranchAddress("cluster_e"  , &cluster_e  , &b_cluster_e  );
  t->SetBranchAddress("cluster_eta", &cluster_eta, &b_cluster_eta);
  t->SetBranchAddress("cluster_phi", &cluster_phi, &b_cluster_phi);
  t->SetBranchAddress("cluster_showershape", cluster_showershape, &b_cluster_showershape);
  t->SetBranchAddress("cluster_bdt_scores", cluster_bdt_scores, &b_cluster_bdt_scores);
  t->SetBranchAddress("cluster_time", &cluster_time, &b_cluster_time);

  t->SetBranchAddress("jet_pt" , jet_pt , &b_jet_pt );
  t->SetBranchAddress("jet_pt_calib" , jet_pt_calib , &b_jet_pt_calib );
  t->SetBranchAddress("jet_pt_recalib" , jet_pt_recalib , &b_jet_pt_recalib );
  t->SetBranchAddress("jet_e"  , jet_e  , &b_jet_e  );
  t->SetBranchAddress("jet_eta", jet_eta, &b_jet_eta);
  t->SetBranchAddress("jet_phi", jet_phi, &b_jet_phi);
  t->SetBranchAddress("jet_emfrac", jet_emfrac, &b_jet_emfrac);
  t->SetBranchAddress("jet_time", jet_time, &b_jet_time);
  
  t->SetBranchAddress("thirdjet_pt", thirdjet_pt, &b_thirdjet_pt);

  if (isMC) {   
    t->SetBranchAddress("truth_cluster_pt" , &truth_cluster_pt , &b_truth_cluster_pt );
    t->SetBranchAddress("truth_cluster_e"  , &truth_cluster_e  , &b_truth_cluster_e  );
    t->SetBranchAddress("truth_cluster_eta", &truth_cluster_eta, &b_truth_cluster_eta);
    t->SetBranchAddress("truth_cluster_phi", &truth_cluster_phi, &b_truth_cluster_phi);
    t->SetBranchAddress("truth_cluster_iso3", &truth_cluster_iso3, &b_truth_cluster_iso3);
    t->SetBranchAddress("truth_cluster_iso4", &truth_cluster_iso4, &b_truth_cluster_iso4);
    
    t->SetBranchAddress("jet_pt_smear_reco" , jet_pt_smear_reco , &b_jet_pt_smear_reco );
    t->SetBranchAddress("jet_pt_smear_high_reco" , jet_pt_smear_high_reco , &b_jet_pt_smear_high_reco );
    t->SetBranchAddress("jet_pt_smear_low_reco"  , jet_pt_smear_low_reco  , &b_jet_pt_smear_low_reco  );
    t->SetBranchAddress("jet_pt_smear_truth" , jet_pt_smear_truth , &b_jet_pt_smear_truth );
    t->SetBranchAddress("jet_pt_smear_high_truth" , jet_pt_smear_high_truth , &b_jet_pt_smear_high_truth );
    t->SetBranchAddress("jet_pt_smear_low_truth"  , jet_pt_smear_low_truth  , &b_jet_pt_smear_low_truth  );
    
    t->SetBranchAddress("truth_jet_pt" , truth_jet_pt , &b_truth_jet_pt );
    t->SetBranchAddress("truth_jet_e"  , truth_jet_e  , &b_truth_jet_e  );
    t->SetBranchAddress("truth_jet_eta", truth_jet_eta, &b_truth_jet_eta);
    t->SetBranchAddress("truth_jet_phi", truth_jet_phi, &b_truth_jet_phi);
    
    t->SetBranchAddress("hadron_p", hadron_p, &b_hadron_p);
  }
}


// Turns off deserialization (TTree::SetBranchStatus(...,0)) for branches treesetup()
// reads into memory - or, for hasthirdjet/thirdjet_eta/thirdjet_phi/thirdjet_dr, never
// calls SetBranchAddress for at all - that src/unfolder.cc and src/unfolder.h
// never read (verified by grepping every treeuser member name against both files; every
// hit traced back to an actual read, not just the member's declaration - RunNumber,
// ScaledTriggerBit, LiveTriggerBit, Scaledowns, mbd_time, the raw uncalibrated jet_pt
// array (only jet_pt_calib is read), jet_pt_recalib, jet_pt_smear_reco/high_reco/
// low_reco, hasthirdjet, thirdjet_eta/phi/dr, and hadron_p all came up empty). Skipping these
// branches' I/O/decompression entirely is a real win since the underlying towerntup
// files are large.
//
// Opt-in (called explicitly by unfolder's constructor), NOT folded into treesetup()
// itself - several other treeuser consumers outside the main unfolding pipeline read
// some of these same branches for one-off checks (e.g. macros/print_truth_photon_pt25_35.C
// and macros/diagnose_nikhil_only.C both read jet_pt_smear_reco), and each of those
// builds its own independent treeuser and TTree, so disabling branches on one instance
// never affects another. If unfolder.cc/.h's cuts ever grow to need one of these
// branches, re-grep before removing it from the list below.
void treeuser::disableBranchesUnusedByUnfolder() {
  if (!t) return;
  t->SetBranchStatus("RunNumber", 0);
  t->SetBranchStatus("mbd_time", 0);
  t->SetBranchStatus("jet_pt", 0);
  t->SetBranchStatus("jet_pt_recalib", 0);
  // Third-jet branches differ between tree versions (thirdjet_dr before Oct 2026,
  // thirdjet_eta/phi after) - only disable the ones this tree actually has.
  for (const char *b : {"hasthirdjet", "thirdjet_eta", "thirdjet_phi", "thirdjet_dr"})
    if (t->GetBranch(b)) t->SetBranchStatus(b, 0);
  if (!isMC) {
    t->SetBranchStatus("ScaledTriggerBit", 0);
    t->SetBranchStatus("LiveTriggerBit", 0);
    t->SetBranchStatus("Scaledowns", 0);
  }
  if (isMC) {
    t->SetBranchStatus("jet_pt_smear_reco", 0);
    t->SetBranchStatus("jet_pt_smear_high_reco", 0);
    t->SetBranchStatus("jet_pt_smear_low_reco", 0);
    t->SetBranchStatus("hadron_p", 0);
  }
}

vector<bool> treeuser::check_keep_MC(float pt_pho, float pt_reco_pho, float pt_jet[], float pt_reco[], string trigger) {
  bool isphoton = (trigger == "Photon5" || trigger == "Photon10" || trigger == "Photon20");
  vector<bool> keep(ana::nJetR + 1); 
  // check photon
  keep[keep.size()-1] = (isphoton ? (pt_pho > threshmap[-1][trigger] && pt_pho < threshmap_high[-1][trigger]) : 1);

  // Check the jets
  //std::cout << "Form keep: " << !isphoton << " " << trigger << " " << pt_reco_pho << std::endl;
  for (int i = 0; i < ana::nJetR; i++) {
    keep[i] = (isphoton ? 1 : (pt_jet[i] > threshmap[i][trigger] && pt_jet[i] < threshmap_high[i][trigger]));// && pt_reco[i] < reco_threshmap_high[i][trigger]));
    if (!isphoton && trigger == "Jet5"  && pt_reco_pho > 12) keep[i] = 0;
    //if (!isphoton && trigger == "Jet8" && pt_reco_pho > 17) keep[i] = 0;
  }

  return keep;
}


