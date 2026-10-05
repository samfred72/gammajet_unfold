#include <iostream>
#include "TFile.h"
#include "TRandom.h"
#include "TF1.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "TLegend.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"    
#include <vector> 
#include <cmath>
#include "style.h"





void smearcheck() {

  // Local multiJet trees (see analysis.cc's header); was SDCC multiJet_skimmed/TREE_MULTIJET_SKIM_pythia_Jet12.root
  const char * filepath = "../trees";
  TFile f12(Form("%s/multijet_%s_Jet%i.root",filepath, "pythia", 12), "READ"); //sim 

  //==============================================
  // TREE SETUP
  //==============================================
  const int ntrees = 1;
  TTree *t12 = (TTree*) f12.Get("ttree");
  TTree * intree[ntrees] = {t12};
  TH2D * hct = new TH2D("hct", ";calib pt;smear pt", 100,0,50,100,0,50);
  TH2D * hcr = new TH2D("hcr", ";calib pt;smear pt", 100,0,50,100,0,50);
  TH1D * hcalib = new TH1D("hcalib",";calib pt;counts",100,0,100);
  TH1D * hreco = new TH1D("hreco",";calib pt;counts",100,0,100);
  TH1D * htruth = new TH1D("htruth",";calib pt;counts",100,0,100);

  std::cout << "Running over trees..." << std::endl;
  for (int i = 0; i < ntrees; i++) {
    std::cout << i << "..." << std::endl;
    bool isMC = true;

    std::vector<float>* TRUTH_pt         = nullptr;
    std::vector<float>* TRUTH_eta        = nullptr;
    std::vector<float>* TRUTH_phi        = nullptr;

    std::vector<float>* jet_pt_smear     = nullptr;
    std::vector<float>* jet_pt_smearRECO = nullptr;
    std::vector<float>* jet_pt_smearHIGH = nullptr;
    std::vector<float>* jet_pt_smearLOW  = nullptr;

    std::vector<float>* jet_pt_calib     = nullptr;
    std::vector<float>* jet_eta_det      = nullptr;
    std::vector<float>* jet_phi          = nullptr;
    std::vector<float>* jet_emfrac       = nullptr;
    std::vector<float>* jet_e            = nullptr;
    
    ULong64_t trigger = 0;
    float zvtx = 0;
    int mbd_hit = 0;
    double calib_lead_time = 0;
    double calib_delta_time = 0;
    int radius = 4;

    if (isMC) {
      intree[i]->SetBranchAddress(Form("truth_jet_pt_%d", radius), &TRUTH_pt);
      intree[i]->SetBranchAddress(Form("truth_jet_eta_%d", radius), &TRUTH_eta);
      intree[i]->SetBranchAddress(Form("truth_jet_phi_%d", radius), &TRUTH_phi);
      // DijetTreeMaker names (old skimmed trees: jet_pt_smear / jet_pt_smearRECO/HIGH/LOW)
      intree[i]->SetBranchAddress(Form("jet_pt_smear_truth_%d", radius), &jet_pt_smear);
      intree[i]->SetBranchAddress(Form("jet_pt_smear_reco_%d", radius), &jet_pt_smearRECO);
      intree[i]->SetBranchAddress(Form("jet_pt_smear_high_reco_%d", radius), &jet_pt_smearHIGH);
      intree[i]->SetBranchAddress(Form("jet_pt_smear_low_reco_%d", radius), &jet_pt_smearLOW);
    }
    else {
      intree[i]->SetBranchAddress("calib_lead_time", &calib_lead_time);
      intree[i]->SetBranchAddress("calib_delta_time", &calib_delta_time);
    }


    intree[i]->SetBranchAddress(Form("jet_pt_calib_%d", radius), &jet_pt_calib);
    intree[i]->SetBranchAddress(Form("jet_eta_det_%d", radius), &jet_eta_det);
    intree[i]->SetBranchAddress(Form("jet_phi_%d", radius), &jet_phi);
    intree[i]->SetBranchAddress(Form("jet_emcal_%d", radius), &jet_emfrac);
    intree[i]->SetBranchAddress(Form("jet_e_%d", radius), &jet_e);
    intree[i]->SetBranchAddress("gl1_scaled", &trigger);
    intree[i]->SetBranchAddress("mbd_vertex_z", &zvtx);
    intree[i]->SetBranchAddress("mbd_hit", &mbd_hit);

    int unpaircount = 0;
    int count = 0;
    for (int e = 0; e < intree[i]->GetEntries(); e++) {
      intree[i]->GetEntry(e);

      float maxcalibpt = 0;
      float maxrecopt = 0;
      float maxtruthpt = 0;
      for (int i = 0; i < jet_pt_calib->size(); i++) {
        if (jet_pt_calib->at(i) > maxcalibpt) {
          maxcalibpt = jet_pt_calib->at(i);
        }
        if (jet_pt_smearRECO->at(i) > maxrecopt) {
          maxrecopt = jet_pt_smearRECO->at(i);
        }
        if (jet_pt_smear->at(i) > maxtruthpt) {
          maxtruthpt = jet_pt_smear->at(i);
        }

        hct->Fill(jet_pt_calib->at(i), jet_pt_smear->at(i));
        hcr->Fill(jet_pt_calib->at(i), jet_pt_smearRECO->at(i));
        if (jet_pt_calib->at(i) > 7 && fabs(jet_pt_calib->at(i) - jet_pt_smear->at(i)) < 0.00001) {
          //cout << "Unpaired jet with pt: " << jet_pt_calib->at(i) << endl;
          unpaircount++;
        }
        count++;
      }
      hcalib->Fill(maxcalibpt);
      hreco->Fill(maxrecopt);
      htruth->Fill(maxtruthpt);
    }
    cout << "Fraction of unmatched jets: " << (float)unpaircount / (float)count << endl;
  }

  //////////////////////// WRITING ////////////////////////

  TFile * wf = TFile::Open(Form("checkhists.root"),"RECREATE");
  hcr->Write();
  hct->Write();
  hcalib->Write();
  hreco->Write();
  htruth->Write();

  std::cout << "done" << std::endl;

}


//g++ analysis.cc style.cc -std=c++11 \
  `root-config --cflags --libs` \
  -o analysis







