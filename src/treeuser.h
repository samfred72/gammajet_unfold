#ifndef TREEUSER_H
#define TREEUSER_H

#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <TSystem.h>
#include <string>
#include "TH1D.h"
#include "TH2D.h"
#include "TBranch.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
using namespace std;

class treeuser {
  public:
    treeuser(string trigger, string sim = "pythia") {
      this->trigger = trigger;
      this->sim = sim;
      // Raw trees are large, pre-existing input data - read from the original gammajet
      // project rather than duplicating them here.
      // "Jet12_full" is not its own trigger - it's a 10x-higher-statistics copy of the
      // same Jet12 towerntup production (1.2M vs 120k entries, identical 39-branch
      // schema), kept in the purity_check project rather than gammajet/trees for its
      // own high-stat ABCD validation. Threshold selection below reuses "Jet12"'s
      // threshmap/threshmap_high/reco_threshmap_high entries since it's physically the
      // same trigger.
      f = (trigger == "Data" ?
          TFile::Open(Form("/home/samson72/sphnx/gammajet/trees/gammajet_%s.root",trigger.c_str()),"read") :
          trigger == "Jet12_full" ?
          TFile::Open("/home/samson72/sphnx/purity_check/gammajet_pythia_Jet12_full.root","read") :
          TFile::Open(Form("/home/samson72/sphnx/gammajet/trees/gammajet_%s_%s.root",sim.c_str(),trigger.c_str()),"read")
          );
      t = (TTree*)f->Get("towerntup");
      isMC = (trigger != "Data");
      treesetup();
      // -1: cluster
      //  0: jet R=0.2
      //  1: jet R=0.3
      //  2: jet R=0.4
      //  3: jet R=0.5
      //  4: jet R=0.6
      //  5: jet R=0.7
      //  6: jet R=0.8
      threshmap = (sim == "pythia" ? 
      std::map<int, std::map<std::string, int>>{
        {-1,{{"Jet5", 0},{"Jet8", 0},{"Jet12", 0},{"Jet20", 0},{"Jet30", 0},{"Jet50", 0},{"Jet70",  0},{"Photon5", 0},{"Photon10",12},{"Photon20", 24}}},
        { 0,{{"Jet5", 0},{"Jet8", 0},{"Jet12",12},{"Jet12_full",0},{"Jet20",20},{"Jet30",30},{"Jet50",50},{"Jet70", 70},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 1,{{"Jet5", 0},{"Jet8", 6},{"Jet12",13},{"Jet12_full",0},{"Jet20",21},{"Jet30",31},{"Jet50",51},{"Jet70", 71},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 2,{{"Jet5", 0},{"Jet8", 0},{"Jet12",14},{"Jet12_full",0},{"Jet20",21},{"Jet30",32},{"Jet50",52},{"Jet70", 71},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 3,{{"Jet5", 0},{"Jet8",10},{"Jet12",19},{"Jet12_full",0},{"Jet20",27},{"Jet30",38},{"Jet50",58},{"Jet70", 79},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 4,{{"Jet5", 0},{"Jet8",12},{"Jet12",22},{"Jet12_full",0},{"Jet20",29},{"Jet30",41},{"Jet50",63},{"Jet70", 79},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 5,{{"Jet5", 0},{"Jet8",14},{"Jet12",24},{"Jet12_full",0},{"Jet20",32},{"Jet30",45},{"Jet50",66},{"Jet70", 80},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 6,{{"Jet5", 0},{"Jet8",15},{"Jet12",25},{"Jet12_full",0},{"Jet20",34},{"Jet30",47},{"Jet50",68},{"Jet70", 80},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}}
      } : 
      std::map<int, std::map<std::string, int>>{
        {-1,{{"Jet5", 0},{"Jet8", 0},{"Jet12", 0},{"Jet20", 0},{"Jet30", 0},{"Jet50", 0},{"Jet70",  0},{"Photon5", 0},{"Photon10", 0},{"Photon20", 20}}}, // because no Photon5
        { 0,{{"Jet5", 0},{"Jet8", 5},{"Jet12",12},{"Jet12_full",0},{"Jet20",20},{"Jet30",30},{"Jet50",50},{"Jet70", 70},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 1,{{"Jet5", 0},{"Jet8", 6},{"Jet12",13},{"Jet12_full",0},{"Jet20",21},{"Jet30",31},{"Jet50",51},{"Jet70", 71},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 2,{{"Jet5", 0},{"Jet8", 7},{"Jet12",14},{"Jet12_full",0},{"Jet20",21},{"Jet30",32},{"Jet50",52},{"Jet70", 71},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 3,{{"Jet5", 0},{"Jet8",10},{"Jet12",19},{"Jet12_full",0},{"Jet20",27},{"Jet30",38},{"Jet50",58},{"Jet70", 79},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 4,{{"Jet5", 0},{"Jet8",12},{"Jet12",22},{"Jet12_full",0},{"Jet20",29},{"Jet30",41},{"Jet50",63},{"Jet70", 79},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 5,{{"Jet5", 0},{"Jet8",14},{"Jet12",24},{"Jet12_full",0},{"Jet20",32},{"Jet30",45},{"Jet50",66},{"Jet70", 80},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 6,{{"Jet5", 0},{"Jet8",15},{"Jet12",25},{"Jet12_full",0},{"Jet20",34},{"Jet30",47},{"Jet50",68},{"Jet70", 80},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}}
      }
      ); 
      threshmap_high = (sim == "pythia" ? 
      std::map<int, std::map<std::string, int>>{
        {-1,{{"Jet5", 0},{"Jet8", 0},{"Jet12", 0},{"Jet20", 0},{"Jet30", 0},{"Jet50", 0},{"Jet70",  0},{"Photon5", 12},{"Photon10",24},{"Photon20",100}}},
        { 0,{{"Jet5", 5},{"Jet8",12},{"Jet12",20},{"Jet12_full",100},{"Jet20",30},{"Jet30",50},{"Jet50",70},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 1,{{"Jet5", 6},{"Jet8",13},{"Jet12",21},{"Jet12_full",100},{"Jet20",31},{"Jet30",51},{"Jet50",71},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 2,{{"Jet5", 7},{"Jet8",14},{"Jet12",21},{"Jet12_full",100},{"Jet20",32},{"Jet30",52},{"Jet50",71},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 3,{{"Jet5",10},{"Jet8",19},{"Jet12",27},{"Jet12_full",100},{"Jet20",38},{"Jet30",58},{"Jet50",79},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 4,{{"Jet5",12},{"Jet8",22},{"Jet12",29},{"Jet12_full",100},{"Jet20",41},{"Jet30",63},{"Jet50",79},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 5,{{"Jet5",14},{"Jet8",24},{"Jet12",32},{"Jet12_full",100},{"Jet20",45},{"Jet30",66},{"Jet50",80},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 6,{{"Jet5",15},{"Jet8",25},{"Jet12",34},{"Jet12_full",100},{"Jet20",47},{"Jet30",68},{"Jet50",80},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}}
      } :
      std::map<int, std::map<std::string, int>>{
        {-1,{{"Jet5", 0},{"Jet8", 0},{"Jet12", 0},{"Jet20", 0},{"Jet30", 0},{"Jet50", 0},{"Jet70",  0},{"Photon5", 0},{"Photon10",20},{"Photon20",100}}}, // because no Photon5
        { 0,{{"Jet5", 5},{"Jet8",12},{"Jet12",20},{"Jet12_full",100},{"Jet20",30},{"Jet30",50},{"Jet50",70},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 1,{{"Jet5", 6},{"Jet8",13},{"Jet12",21},{"Jet12_full",100},{"Jet20",31},{"Jet30",51},{"Jet50",71},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 2,{{"Jet5", 7},{"Jet8",14},{"Jet12",21},{"Jet12_full",100},{"Jet20",32},{"Jet30",52},{"Jet50",71},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 3,{{"Jet5",10},{"Jet8",19},{"Jet12",27},{"Jet12_full",100},{"Jet20",38},{"Jet30",58},{"Jet50",79},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 4,{{"Jet5",12},{"Jet8",22},{"Jet12",29},{"Jet12_full",100},{"Jet20",41},{"Jet30",63},{"Jet50",79},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 5,{{"Jet5",14},{"Jet8",24},{"Jet12",32},{"Jet12_full",100},{"Jet20",45},{"Jet30",66},{"Jet50",80},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 6,{{"Jet5",15},{"Jet8",25},{"Jet12",34},{"Jet12_full",100},{"Jet20",47},{"Jet30",68},{"Jet50",80},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}}
      }
      );
      reco_threshmap_high = std::map<int, std::map<std::string, int>>{
        {-1,{{"Jet5",12},{"Jet8",20},{"Jet12",30},{"Jet12_full",100},{"Jet20",42},{"Jet30",60},{"Jet50",83},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20", 0}}},
        { 0,{{"Jet5",12},{"Jet8",20},{"Jet12",30},{"Jet12_full",100},{"Jet20",42},{"Jet30",60},{"Jet50",83},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20", 0}}},
        { 1,{{"Jet5",12},{"Jet8",20},{"Jet12",30},{"Jet12_full",100},{"Jet20",42},{"Jet30",60},{"Jet50",83},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20", 0}}},
        { 2,{{"Jet5",12},{"Jet8",20},{"Jet12",30},{"Jet12_full",100},{"Jet20",42},{"Jet30",60},{"Jet50",83},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20", 0}}},
        { 3,{{"Jet5",12},{"Jet8",20},{"Jet12",30},{"Jet12_full",100},{"Jet20",42},{"Jet30",60},{"Jet50",83},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20", 0}}},
        { 4,{{"Jet5",12},{"Jet8",20},{"Jet12",30},{"Jet12_full",100},{"Jet20",42},{"Jet30",60},{"Jet50",83},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20", 0}}},
        { 5,{{"Jet5",12},{"Jet8",20},{"Jet12",30},{"Jet12_full",100},{"Jet20",42},{"Jet30",60},{"Jet50",83},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20", 0}}},
        { 6,{{"Jet5",12},{"Jet8",20},{"Jet12",30},{"Jet12_full",100},{"Jet20",42},{"Jet30",60},{"Jet50",83},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20", 0}}}
      };
    }
    void treesetup();
    vector<bool> check_keep_MC(float pho_pt, float pt_reco_pho, float jet_pt[], float pt_reco[], string trigger);

    bool isMC;
    TFile * f;
    TTree * t;
    string trigger;
    string sim;
    const static int nJetR = ana::nJetR;
    map<int, map<string,int>> threshmap;
    map<int, map<string,int>> threshmap_high;
    map<int, map<string,int>> reco_threshmap_high;


    // TTree variables
    Int_t           RunNumber;
    Float_t         vz;
    Bool_t          ScaledTriggerBit[64];
    Bool_t          LiveTriggerBit[64];
    Int_t           Scaledowns[64];

    Float_t         mbd_time;

    Float_t         cluster_pt;
    Float_t         cluster_e;
    Float_t         cluster_eta;
    Float_t         cluster_phi;
    Float_t         cluster_time;
    Float_t         cluster_bdt_scores[11];
    Float_t         cluster_showershape[12]; // [10]=iso_topo_03, [11]=iso_topo_04

    Float_t         truth_cluster_pt;
    Float_t         truth_cluster_e;
    Float_t         truth_cluster_eta;
    Float_t         truth_cluster_phi;
    Float_t         truth_cluster_iso3;
    Float_t         truth_cluster_iso4;

    Float_t         jet_pt[nJetR];
    Float_t         jet_pt_calib[nJetR];
    Float_t         jet_pt_recalib[nJetR];
    Float_t         jet_pt_smear[nJetR];
    Float_t         jet_pt_smear_reco[nJetR];
    Float_t         jet_pt_smear_high[nJetR];
    Float_t         jet_pt_smear_low[nJetR];
    Float_t         jet_e[nJetR];
    Float_t         jet_eta[nJetR];
    Float_t         jet_phi[nJetR];
    Float_t         jet_emfrac[nJetR];
    Float_t         jet_time[nJetR]; 

    Bool_t          hasthirdjet[nJetR];
    Float_t         thirdjet_pt[nJetR];
    Float_t         thirdjet_dr[nJetR];

    Float_t         hadron_p[nJetR];

    Float_t         truth_jet_pt[nJetR];
    Float_t         truth_jet_e[nJetR];
    Float_t         truth_jet_eta[nJetR];
    Float_t         truth_jet_phi[nJetR];

    // List of branches
    TBranch        *b_RunNumber;   //!
    TBranch        *b_vz;   //!
    TBranch        *b_ScaledTriggerBit;   //!
    TBranch        *b_LiveTriggerBit;   //!
    TBranch        *b_Scaledowns;   //!
    TBranch        *b_mbd_time;   //!

    TBranch        *b_cluster_pt;
    TBranch        *b_cluster_e;
    TBranch        *b_cluster_eta;
    TBranch        *b_cluster_phi;
    TBranch        *b_cluster_bdt_scores;   //!
    TBranch        *b_cluster_showershape;   //!
    TBranch        *b_cluster_time;   //!

    TBranch        *b_truth_cluster_pt;
    TBranch        *b_truth_cluster_e;
    TBranch        *b_truth_cluster_eta;
    TBranch        *b_truth_cluster_phi;
    TBranch        *b_truth_cluster_iso3;
    TBranch        *b_truth_cluster_iso4;

    TBranch        *b_jet_pt;                           //!
    TBranch        *b_jet_pt_calib;                           //!
    TBranch        *b_jet_pt_recalib;                           //!
    TBranch        *b_jet_pt_smear;                           //!
    TBranch        *b_jet_pt_smear_reco;                           //!
    TBranch        *b_jet_pt_smear_high;                           //!
    TBranch        *b_jet_pt_smear_low;                           //!
    TBranch        *b_jet_e;                           //!
    TBranch        *b_jet_eta;                           //!
    TBranch        *b_jet_phi;                           //!
    TBranch        *b_jet_emfrac;   //!
    TBranch        *b_jet_time;   //!
    
    TBranch        *b_hasthirdjet;   //!
    TBranch        *b_thirdjet_pt;   //!
    TBranch        *b_thirdjet_dr;   //!

    TBranch        *b_hadron_p;   //!

    TBranch        *b_truth_jet_pt;                           //!
    TBranch        *b_truth_jet_e;                           //!
    TBranch        *b_truth_jet_eta;                           //!
    TBranch        *b_truth_jet_phi;                           //!
};
#endif // TREEUSER_H
