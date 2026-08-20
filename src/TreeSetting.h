#ifndef TreeSetting_h
#define TreeSetting_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "vector"

void treesetup(TTree *tree, bool isMC = 0)
{
  cluster_pt = 0;
  cluster_e = 0;
  cluster_eta = 0;
  cluster_phi = 0;
  cluster_time = 0;
  cluster_bdt_scores = 0;
  cluster_showershape = 0;

  truth_cluster_pt = 0;
  truth_cluster_e = 0;
  truth_cluster_eta = 0;
  truth_cluster_phi = 0;

  jet_pt = 0;
  jet_pt_calib = 0;
  jet_pt_smear = 0;
  jet_e = 0;
  jet_eta = 0;
  jet_phi = 0;
  jet_emfrac = 0;
  jet_ihfrac = 0;
  jet_ohfrac = 0;
  jet_time = 0;  

  truth_jet_pt = 0;
  truth_jet_e = 0;
  truth_jet_eta = 0;
  truth_jet_phi = 0;

  // Set branch addresses and branch pointers
  if (!tree) return;

  tree->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
  tree->SetBranchAddress("vz", &vz, &b_vz);
  tree->SetBranchAddress("ScaledTriggerBit", ScaledTriggerBit, &b_ScaledTriggerBit);
  tree->SetBranchAddress("LiveTriggerBit", LiveTriggerBit, &b_LiveTriggerBit);
  tree->SetBranchAddress("Scaledowns", Scaledowns, &b_Scaledowns);
  tree->SetBranchAddress("mbd_nhits_south", &mbd_nhits_south, &b_mbd_nhits_south);
  tree->SetBranchAddress("mbd_nhits_north", &mbd_nhits_north, &b_mbd_nhits_north);
  tree->SetBranchAddress("mbd_time_south", &mbd_time_south, &b_mbd_time_south);
  tree->SetBranchAddress("mbd_time_north", &mbd_time_north, &b_mbd_time_north);

  tree->SetBranchAddress("nClusters", &nClusters, &b_nClusters);
  tree->SetBranchAddress("cluster_pt" , &cluster_pt , &b_cluster_pt );
  tree->SetBranchAddress("cluster_e"  , &cluster_e  , &b_cluster_e  );
  tree->SetBranchAddress("cluster_eta", &cluster_eta, &b_cluster_eta);
  tree->SetBranchAddress("cluster_phi", &cluster_phi, &b_cluster_phi);
  tree->SetBranchAddress("cluster_showershape", &cluster_showershape, &b_cluster_showershape);
  tree->SetBranchAddress("cluster_time", &cluster_time, &b_cluster_time);
  tree->SetBranchAddress("cluster_bdt_scores", &cluster_bdt_scores, &b_cluster_bdt_scores);

  tree->SetBranchAddress("nJets", nJets, &b_nJets);
  tree->SetBranchAddress("jet_pt" , &jet_pt , &b_jet_pt );
  tree->SetBranchAddress("jet_pt_calib" , &jet_pt_calib , &b_jet_pt_calib );
  tree->SetBranchAddress("jet_e"  , &jet_e  , &b_jet_e  );
  tree->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
  tree->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
  tree->SetBranchAddress("jet_emfrac", &jet_emfrac, &b_jet_emfrac);
  tree->SetBranchAddress("jet_ihfrac", &jet_ihfrac, &b_jet_ihfrac);
  tree->SetBranchAddress("jet_ohfrac", &jet_ohfrac, &b_jet_ohfrac);
  tree->SetBranchAddress("jet_time", &jet_time, &b_jet_time);

  if (isMC) {   
    tree->SetBranchAddress("jet_pt_smear" , &jet_pt_smear , &b_jet_pt_smear );
    tree->SetBranchAddress("nTruthClusters", &nTruthClusters, &b_nTruthClusters);
    tree->SetBranchAddress("truth_cluster_pt" , &truth_cluster_pt , &b_truth_cluster_pt );
    tree->SetBranchAddress("truth_cluster_e"  , &truth_cluster_e  , &b_truth_cluster_e  );
    tree->SetBranchAddress("truth_cluster_eta", &truth_cluster_eta, &b_truth_cluster_eta);
    tree->SetBranchAddress("truth_cluster_phi", &truth_cluster_phi, &b_truth_cluster_phi);
    
    tree->SetBranchAddress("nTruthJets", nTruthJets, &b_nTruthJets);
    tree->SetBranchAddress("truth_jet_pt" , &truth_jet_pt , &b_truth_jet_pt );
    tree->SetBranchAddress("truth_jet_e"  , &truth_jet_e  , &b_truth_jet_e  );
    tree->SetBranchAddress("truth_jet_eta", &truth_jet_eta, &b_truth_jet_eta);
    tree->SetBranchAddress("truth_jet_phi", &truth_jet_phi, &b_truth_jet_phi);
  }
}
#endif
