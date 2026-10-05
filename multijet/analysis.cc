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
#include "TMarker.h"
#include <vector>
#include <cmath>
#include <iomanip>
#include "style.h"

// Multijet balance analysis (from SDCC multiJet_legacy/analysis.cc, Oct 2026).
// Reads the multiJet analysis trees (multiJet/src/DijetTreeMaker.cc) - the same local
// copies the in-situ combined mode uses (../trees/multijet_*.root, see insitu/README.md):
//   Data: multijet_Data.root, MC: multijet_<sim>_Jet{8,12,20,30}.root
// Build: ./make.sh    Run (from this directory): ./analysis <0=pythia|1=herwig> [tree dir]
// Writes multijet_analysis_<sim>.root and pdfs/unmatched_event_display_<sim>.pdf; reads the
// MC reweighting fits aux/ratio<R>_<sim>.root (written by makeratio.C from a previous pass).

float deltaPhi(float phi1, float phi2) {
  float dphi = std::fabs(phi1 - phi2);
  if (dphi > M_PI) dphi = 2 * M_PI - dphi;
  return dphi;
}

std::vector<int> IDXGrab(int nJets, const std::vector<float>& pT) {
  std::vector<int> idx;
  idx.reserve(3);

  int max_idx  = -1; float max_pT  = -1;
  int max2_idx = -1; float max2_pT = -1;
  int max3_idx = -1; float max3_pT = -1;

  for (int i = 0; i < nJets; i++) {
    float current_pT = pT[i];

    if (current_pT > max_pT) {
      max3_pT  = max2_pT;
      max3_idx = max2_idx;

      max2_pT  = max_pT;
      max2_idx = max_idx;

      max_pT   = current_pT;
      max_idx  = i;

    } else if (current_pT > max2_pT) {
      max3_pT  = max2_pT;
      max3_idx = max2_idx;

      max2_pT  = current_pT;
      max2_idx = i;

    } else if (current_pT > max3_pT) {
      max3_pT  = current_pT;
      max3_idx = i;
    }
  }

  idx.push_back(max_idx);
  idx.push_back(max2_idx);
  idx.push_back(max3_idx);

  return idx;
}

float newCoordinates(float pT , float phi,  float pT2 , float phi2 , bool phiCoordinate = false){
  float pTNew_x = pT * cos(phi) + pT2 * cos(phi2);
  float pTNew_y = pT * sin(phi) + pT2 * sin(phi2);
  float phiNew = atan2(pTNew_y, pTNew_x);

  if(phiCoordinate) return phiNew;

  return sqrt(pTNew_x*pTNew_x + pTNew_y*pTNew_y);
}

int main(int argc, char* argv[]) {

  gROOT->SetBatch(true);
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <0=pythia|1=herwig> [tree dir, default ../trees]" << std::endl;
    return 1;
  }
  bool isherwig = (strcmp(argv[1], "1") == 0);
  const char * treedir = (argc > 2 ? argv[2] : "../trees");
  const char * sim = (isherwig ? "herwig" : "pythia");
  std::cout << "sim is: " << argv[1] << " " << sim << std::endl;

  // All 7 jet radii processed in one pass per sim - each input tree already
  // carries every radius's branches for a given event, so one GetEntry()
  // gives us all of them at once.
  const int nRadii = 7;
  const int radii[nRadii] = {2,3,4,5,6,7,8};
  const int reprRadiusIdx = 2; // radius 4 - representative sample for the debug event display only

  // Only the reco-based JER smear variants: RECO is the central value
  // (smearing applied directly to the calibrated reco pt), HIGH/LOW are its
  // systematic envelope. The truth-anchored "smear" variant is intentionally
  // not produced here.
  const int nsystypes = 3;
  const char * sysNames[nsystypes] = {"RECO","HIGH","LOW"};
  const char * fitFuncNames[nsystypes] = {"ratio_func_JERreco","ratio_func_JERhigh","ratio_func_JERlow"};

  //=====================================
  // file setup
  //=====================================
  TFile infile(Form("%s/multijet_Data.root", treedir), "READ"); // Data
  TFile f08(Form("%s/multijet_%s_Jet8.root",  treedir, sim), "READ"); //sim
  TFile f12(Form("%s/multijet_%s_Jet12.root", treedir, sim), "READ"); //sim
  TFile f20(Form("%s/multijet_%s_Jet20.root", treedir, sim), "READ"); //sim
  TFile f30(Form("%s/multijet_%s_Jet30.root", treedir, sim), "READ"); //sim

  if (infile.IsZombie()) return 1;
  if (f08.IsZombie()) return 1;
  if (f12.IsZombie()) return 1;
  if (f20.IsZombie()) return 1;
  if (f30.IsZombie()) return 1;

  //====================================
  // WEIGHTING
  //====================================

  float leading_pT_Cutoff = 20;
  float subleadingPTCutoff = 7;
  float SSLCutOff = 7;
  float SSLDPHI = M_PI/2.0;
  float SLDPHI = 3*M_PI/4.0;

  double weight8, weight12, weight20, weight30;
  //pythia weights
  if (!isherwig) {
    weight8  = 1.3013e+07 / 2.5298e+03;// * 10000000.0 / 9998000.0;
    weight12 = 1.4903e+06 / 2.5298e+03;
    weight20 = 6.2623e+04 / 2.5298e+03;
    weight30 = 1;
  }
  //herwig weights
  else {
    weight8  = 1.8437e+08 / 2.0694e+03;
    weight12 = 1.1324e+06 / 2.0694e+03 * 10000000.0 / 10001000.0;
    weight20 = 5.2613e+04 / 2.0694e+03 * 10000000.0 / 10913000.0;
    weight30 = 1 * 10000000.0 / 9999000.0;
  }
  const int nsimtrees = 4;
  const int ntrees = nsimtrees + 1;
  double weights[ntrees] = {1, weight8, weight12, weight20, weight30};

  // truth-pT-hat stitching windows, per radius (row) per tree (column: data,Jet8,Jet12,Jet20,Jet30)
  double cutLow[nRadii][4] = {
    {0,12,20,30}, {0,13,21,31}, {0,14,21,32}, {0,19,27,38},
    {0,22,29,41}, {0,24,32,45}, {0,25,34,47}
  };
  double cutHigh[nRadii][4] = {
    {12,20,30,100}, {13,21,31,100}, {14,21,32,100}, {19,27,38,100},
    {22,29,41,100}, {24,32,45,100}, {25,34,47,100}
  };
  double lowcuts[nRadii][ntrees];
  double highcuts[nRadii][ntrees];
  for (int ir = 0; ir < nRadii; ir++) {
    lowcuts[ir][0] = 0; highcuts[ir][0] = 100;
    for (int k = 0; k < 4; k++) {
      lowcuts[ir][k+1]  = cutLow[ir][k];
      highcuts[ir][k+1] = cutHigh[ir][k];
    }
  }

  //==============================================================
  // per-radius ratio/fit inputs (produced by makeratio.C from the previous
  // iteration's leadingJetPT/leadingJetPT_Pyth spectra below)
  //==============================================================
  TFile* fRatio[nRadii];
  TF1* myFit[nRadii][nsystypes];
  TH1D* zvtxRatio[nRadii];
  for (int ir = 0; ir < nRadii; ir++) {
    fRatio[ir] = TFile::Open(Form("aux/ratio%d_%s.root", radii[ir], sim), "READ");
    if (!fRatio[ir] || fRatio[ir]->IsZombie()) {
      std::cerr << "Missing ratio file for radius " << radii[ir] << std::endl;
      return 1;
    }
    for (int j = 0; j < nsystypes; j++) {
      myFit[ir][j] = (TF1*)fRatio[ir]->Get(fitFuncNames[j]);
    }
    zvtxRatio[ir] = (TH1D*)fRatio[ir]->Get(Form("hratio_zvtx%d_%s", radii[ir], sim));
  }

  //==============================================================
  // OUTPUT FILE - one file for this sim: all radii, all systematics
  //==============================================================
  TFile * wf = TFile::Open(Form("multijet_analysis_%s.root", sim), "RECREATE");
  wf->cd();
  gStyle->SetOptStat(0);

  //==============================================================
  // OUTPUT HISTOGRAM SETUP
  //==============================================================
  std::vector<float> pTBins = {20,25,30,35,40,50,60,70};
  const int nPtBins = 7;

  TH1D* leadingJetPT[nRadii];
  TH1D* leadingJetPT_Pyth[nRadii][nsystypes];
  TH1D* leadingJetPT_truth[nRadii];
  TH1D* dphi12[nRadii];
  TH1D* dphi13[nRadii];
  TH1D* dphi1_23[nRadii];
  TH1D* dphi1_23_data[nRadii];
  TH1D* heta0[nRadii];
  TH1D* heta1[nRadii];
  TH1D* heta2[nRadii];
  TH1D* zvtx_data[nRadii];
  TH1D* zvtx_MC[nRadii];
  const char * hxjnames[2] = {"data","sim"};
  TH1D* hxj[nRadii][nPtBins][nsystypes][2]; // for data ([...][0]) only the j=0 slot is used

  for (int ir = 0; ir < nRadii; ir++) {
    int radius = radii[ir];
    leadingJetPT[ir] = new TH1D(Form("leadingJet_r%d", radius), ";Calibrated Leading Jet p_{T} (GeV);Counts", 40, 20, 60);
    for (int j = 0; j < nsystypes; j++) {
      leadingJetPT_Pyth[ir][j] = new TH1D(Form("leadingJetPT_Pyth_r%d_%s", radius, sysNames[j]), ";Uncalibrated Leading Jet p_{T} (GeV);Counts", 40, 20, 60);
    }
    leadingJetPT_truth[ir] = new TH1D(Form("leadingJetPT_truth_r%d", radius), ";Uncalibrated Leading Jet p_{T} (GeV);Counts", 40, 20, 60);
    dphi12[ir]         = new TH1D(Form("dphi12_r%d", radius), ";dphi;counts", 100, 0, M_PI);
    dphi13[ir]         = new TH1D(Form("dphi13_r%d", radius), ";dphi;counts", 100, 0, M_PI);
    dphi1_23[ir]       = new TH1D(Form("dphi1_23_r%d", radius), ";dphi;counts", 100, 0, M_PI);
    dphi1_23_data[ir]  = new TH1D(Form("dphi1_23_data_r%d", radius), ";dphi;counts", 100, 0, M_PI);
    heta0[ir]          = new TH1D(Form("heta0_r%d", radius), ";eta;counts", 100, -1.2, 1.2);
    heta1[ir]          = new TH1D(Form("heta1_r%d", radius), ";eta;counts", 100, -1.2, 1.2);
    heta2[ir]          = new TH1D(Form("heta2_r%d", radius), ";eta;counts", 100, -1.2, 1.2);
    zvtx_data[ir]      = new TH1D(Form("zvtx_data_r%d", radius), ";z (mm);counts", 400, -1000, 1000);
    zvtx_MC[ir]        = new TH1D(Form("zvtx_MC_r%d", radius), ";z (mm);counts", 400, -1000, 1000);

    for (int k = 0; k < nPtBins; k++) {
      hxj[ir][k][0][0] = new TH1D(Form("hxj_r%d_%1.0f_data", radius, pTBins[k]), ";x_{j};#frac{1}{N}#frac{dN}{dx_j}", 45, 0.4, 2.65);
      hxj[ir][k][0][0]->Sumw2();
      for (int j = 0; j < nsystypes; j++) {
        hxj[ir][k][j][1] = new TH1D(Form("hxj_r%d_%1.0f_%s_sim", radius, pTBins[k], sysNames[j]), ";x_{j};#frac{1}{N}#frac{dN}{dx_j}", 45, 0.4, 2.65);
        hxj[ir][k][j][1]->Sumw2();
      }
    }
  }

  //==============================================================
  // OUTPUT TTREE SETUP
  //==============================================================
  const char * treenames[ntrees] = {"data", "Jet8", "Jet12", "Jet20", "Jet30"};

  TTree* outtree[nRadii][ntrees][nsystypes]; // for data (i=0) only the j=0 slot is used

  float LeadingPT[nRadii][ntrees][nsystypes];
  float SLPT     [nRadii][ntrees][nsystypes];
  float SSLPT    [nRadii][ntrees][nsystypes];
  float SLeta    [nRadii][ntrees][nsystypes];
  float SLphi    [nRadii][ntrees][nsystypes];
  float SSLeta   [nRadii][ntrees][nsystypes];
  float SSLphi   [nRadii][ntrees][nsystypes];
  float PT23     [nRadii][ntrees][nsystypes];
  float weight   [nRadii][ntrees][nsystypes];

  for (int ir = 0; ir < nRadii; ir++) {
    int radius = radii[ir];
    for (int i = 0; i < ntrees; i++) {
      bool isMC = (i != 0);
      int nJ = (isMC ? nsystypes : 1);
      for (int j = 0; j < nJ; j++) {
        std::string treename = (!isMC ? Form("ttree_data_r%d", radius)
                                       : Form("ttree_%s_r%d_%s", treenames[i], radius, sysNames[j]));
        outtree[ir][i][j] = new TTree(treename.c_str(), treename.c_str());
        outtree[ir][i][j]->Branch("leadingPT", &LeadingPT[ir][i][j], "leadingPT/F");
        outtree[ir][i][j]->Branch("SLPT",      &SLPT[ir][i][j],      "SLPT/F");
        outtree[ir][i][j]->Branch("SSLPT",     &SSLPT[ir][i][j],     "SSLPT/F");
        outtree[ir][i][j]->Branch("SLeta",     &SLeta[ir][i][j],     "SLeta/F");
        outtree[ir][i][j]->Branch("SLphi",     &SLphi[ir][i][j],     "SLphi/F");
        outtree[ir][i][j]->Branch("SSLeta",    &SSLeta[ir][i][j],    "SSLeta/F");
        outtree[ir][i][j]->Branch("SSLphi",    &SSLphi[ir][i][j],    "SSLphi/F");
        outtree[ir][i][j]->Branch("PT23",      &PT23[ir][i][j],      "PT23/F");
        outtree[ir][i][j]->Branch("weight",    &weight[ir][i][j],    "weight/F");
      }
    }
  }

  //==============================================
  // TREE SETUP
  //==============================================
  TTree* TreeRead = (TTree*) infile.Get("ttree");
  TTree *t08 = (TTree*) f08.Get("ttree");
  TTree *t12 = (TTree*) f12.Get("ttree");
  TTree *t20 = (TTree*) f20.Get("ttree");
  TTree *t30 = (TTree*) f30.Get("ttree");

  TTree * intree[ntrees] = {TreeRead, t08, t12, t20, t30};

  TCanvas * c = new TCanvas("c_unmatched","",500,1000);
  c->SaveAs(Form("pdfs/unmatched_event_display_%s.pdf[", sim));
  int icount = 0;

  std::cout << "Running over trees..." << std::endl;
  for (int i = 0; i < ntrees; i++) {
    std::cout << i << "..." << std::endl;
    bool isMC = (i != 0);
    int nJloop = (isMC ? nsystypes : 1);

    std::vector<float>* TRUTH_pt[nRadii]  = {nullptr};
    std::vector<float>* TRUTH_eta[nRadii] = {nullptr};
    std::vector<float>* TRUTH_phi[nRadii] = {nullptr};

    std::vector<float>* jet_pt_smearRECO[nRadii] = {nullptr};
    std::vector<float>* jet_pt_smearHIGH[nRadii] = {nullptr};
    std::vector<float>* jet_pt_smearLOW[nRadii]  = {nullptr};

    std::vector<float>* jet_pt_calib[nRadii] = {nullptr};
    std::vector<float>* jet_eta_det[nRadii]  = {nullptr};
    std::vector<float>* jet_phi[nRadii]      = {nullptr};

    ULong64_t trigger = 0;
    float zvtx = 0;
    int mbd_hit = 0;
    double calib_lead_time = 0;
    double calib_delta_time = 0;

    for (int ir = 0; ir < nRadii; ir++) {
      int radius = radii[ir];
      if (isMC) {
        intree[i]->SetBranchAddress(Form("truth_jet_pt_%d", radius), &TRUTH_pt[ir]);
        intree[i]->SetBranchAddress(Form("truth_jet_eta_%d", radius), &TRUTH_eta[ir]);
        intree[i]->SetBranchAddress(Form("truth_jet_phi_%d", radius), &TRUTH_phi[ir]);
        // DijetTreeMaker names (the old separately-skimmed trees used jet_pt_smearRECO/HIGH/LOW_<R>)
        intree[i]->SetBranchAddress(Form("jet_pt_smear_reco_%d", radius), &jet_pt_smearRECO[ir]);
        intree[i]->SetBranchAddress(Form("jet_pt_smear_high_reco_%d", radius), &jet_pt_smearHIGH[ir]);
        intree[i]->SetBranchAddress(Form("jet_pt_smear_low_reco_%d", radius), &jet_pt_smearLOW[ir]);
      }
      intree[i]->SetBranchAddress(Form("jet_pt_calib_%d", radius), &jet_pt_calib[ir]);
      intree[i]->SetBranchAddress(Form("jet_eta_det_%d", radius), &jet_eta_det[ir]);
      intree[i]->SetBranchAddress(Form("jet_phi_%d", radius), &jet_phi[ir]);
    }
    if (!isMC) {
      intree[i]->SetBranchAddress("calib_lead_time", &calib_lead_time);
      intree[i]->SetBranchAddress("calib_delta_time", &calib_delta_time);
    }
    intree[i]->SetBranchAddress("gl1_scaled", &trigger);
    intree[i]->SetBranchAddress("mbd_vertex_z", &zvtx);
    intree[i]->SetBranchAddress("mbd_hit", &mbd_hit);

    //////////////////////// TREE ANALYSIS ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int matched_events[nRadii] = {0};
    int passed_events[nRadii] = {0};

    Long64_t nentries = intree[i]->GetEntries();
    for (Long64_t e = 0; e < nentries; e++) {
      intree[i]->GetEntry(e);
      if (e % 1000 == 0) std::cout << "entry " << e << "/" << nentries << " (" << (float)e/nentries*100. << "%)" << "\t\r" << std::flush;

      // =========================================
      // Event Level Cuts (radius-independent)
      // =========================================

      if (std::fabs(zvtx) > 60) continue;
      //if (!mbd_hit) continue;

      if (!isMC) {
        // Timing cut
        double x = calib_lead_time;
        double y = x - calib_delta_time;
        if (!(fabs(x) < 6 && fabs(y) < 3)) continue; //requiring timing cut
      }

      // Trigger selection
      bool bit[64];
      for (int b = 0; b < 64; b++) {
        bit[b] = (((trigger >> b) & 0x1) == 0x1);
      }
      if (!isMC && !bit[22]) continue;

      for (int ir = 0; ir < nRadii; ir++) {
        int radius = radii[ir];
        float R = radius / 10.0f;

        if (isMC) zvtx_MC[ir]->Fill(zvtx, weights[i]);
        if (!isMC) zvtx_data[ir]->Fill(zvtx, weights[i]);

        if (jet_pt_calib[ir]->size() <= 2) continue; // requiring three jets

        if (isMC) {
          // Truth selection cuts
          if (TRUTH_pt[ir]->size() == 0) continue; // requiring at least one truth jet

          std::vector<int> IDXTRUTH = IDXGrab(TRUTH_pt[ir]->size(), *TRUTH_pt[ir]);
          if (lowcuts[ir][i] > TRUTH_pt[ir]->at(IDXTRUTH[0]) ||
              TRUTH_pt[ir]->at(IDXTRUTH[0]) >= highcuts[ir][i]) continue;
          leadingJetPT_truth[ir]->Fill(TRUTH_pt[ir]->at(IDXTRUTH[0]), weights[i]);
        }

        //==========================================================
        // LOOP OVER RECO / HIGH / LOW JER-SMEARED PT VECTORS
        //==========================================================

        std::vector<std::vector<float>*> ptVariations = (!isMC ?
            std::vector<std::vector<float>*>{ jet_pt_calib[ir] } :
            std::vector<std::vector<float>*>{
              jet_pt_smearRECO[ir],
              jet_pt_smearHIGH[ir],
              jet_pt_smearLOW[ir]
            });

        for (int j = 0; j < nJloop; j++) {

          std::vector<float>* currentPt = ptVariations[j];
          std::vector<int> Idx_Jets = IDXGrab(currentPt->size(), *currentPt);

          if (Idx_Jets[0] == -1 || Idx_Jets[1] == -1 || Idx_Jets[2] == -1) continue;

          float dPhi13 = deltaPhi(
              jet_phi[ir]->at(Idx_Jets[0]),
              jet_phi[ir]->at(Idx_Jets[2])
              );

          float dPhi12 = deltaPhi(
              jet_phi[ir]->at(Idx_Jets[0]),
              jet_phi[ir]->at(Idx_Jets[1])
              );
          float phi23 = newCoordinates(currentPt->at(Idx_Jets[1]), jet_phi[ir]->at(Idx_Jets[1]),
                                       currentPt->at(Idx_Jets[2]), jet_phi[ir]->at(Idx_Jets[2]), true);
          float dPhi123 = deltaPhi(
              jet_phi[ir]->at(Idx_Jets[0]),
              phi23
              );

          if (i == 1 && j == 0) {
            dphi12[ir]->Fill(dPhi12);
            dphi13[ir]->Fill(dPhi13);
            heta0[ir]->Fill(jet_eta_det[ir]->at(Idx_Jets[0]));
            heta1[ir]->Fill(jet_eta_det[ir]->at(Idx_Jets[1]));
            heta2[ir]->Fill(jet_eta_det[ir]->at(Idx_Jets[2]));
          }

          if (leading_pT_Cutoff > currentPt->at(Idx_Jets[0]) ||
              subleadingPTCutoff > currentPt->at(Idx_Jets[1]) ||
              SSLCutOff > currentPt->at(Idx_Jets[2]) ||
              30 <= currentPt->at(Idx_Jets[1]) ||
              30 <= currentPt->at(Idx_Jets[2]) ||
              std::fabs(jet_eta_det[ir]->at(Idx_Jets[0])) > (1.1-R) ||
              std::fabs(jet_eta_det[ir]->at(Idx_Jets[1])) > (1.1-R) ||
              std::fabs(jet_eta_det[ir]->at(Idx_Jets[2])) > (1.1-R) ||
              std::fabs(dPhi13) < SSLDPHI ||
              std::fabs(dPhi12) < SLDPHI
             ) continue;

          if (isMC && j == 0) {
            bool match0 = (std::fabs(currentPt->at(Idx_Jets[0]) - jet_pt_calib[ir]->at(Idx_Jets[0])) > 0.00001);
            bool match1 = (std::fabs(currentPt->at(Idx_Jets[1]) - jet_pt_calib[ir]->at(Idx_Jets[1])) > 0.00001);
            bool match2 = (std::fabs(currentPt->at(Idx_Jets[2]) - jet_pt_calib[ir]->at(Idx_Jets[2])) > 0.00001);

            if (match0 && match1 && match2) { matched_events[ir]++; }
            else if (ir == reprRadiusIdx && i == 1 && icount < 100) { // Draw a bad event (radius 4, Jet8 MC, representative sample)
              c->cd();

              TH2D *h = new TH2D(Form("heventdisplay_%d",icount), ";#eta;#phi", 100,-1.5,1.5, 100,-M_PI,M_PI);
              h->Draw();

              // Legend
              TLegend *leg = new TLegend(0.60,0.90,0.90,1.00);
              leg->SetBorderSize(0);
              leg->SetFillStyle(0);

              TMarker *recoLeg = new TMarker(0,0,107);
              recoLeg->SetMarkerColor(kRed);
              recoLeg->SetMarkerSize(2.0);

              TMarker *truthLeg = new TMarker(0,0,107);
              truthLeg->SetMarkerColor(kGreen+2);
              truthLeg->SetMarkerSize(2.0);

              leg->AddEntry(recoLeg,  "Reco jets",  "p");
              leg->AddEntry(truthLeg, "Truth jets", "p");

              // Reco jets
              for (int ijet = 0; ijet < jet_pt_calib[ir]->size(); ijet++) {
                double eta = jet_eta_det[ir]->at(ijet);
                double phi = jet_phi[ir]->at(ijet);
                double pt  = jet_pt_calib[ir]->at(ijet);

                TMarker *mjet = new TMarker(eta, phi, 107);
                mjet->SetMarkerColor(kRed);
                mjet->SetMarkerSize(15);
                mjet->Draw();
              }

              // Truth jets
              for (int ijet = 0; ijet < TRUTH_pt[ir]->size(); ijet++) {
                double eta = TRUTH_eta[ir]->at(ijet);
                double phi = TRUTH_phi[ir]->at(ijet);

                TMarker *mjet = new TMarker(eta, phi, 107);
                mjet->SetMarkerColor(kGreen+2);
                mjet->SetMarkerSize(15);
                mjet->Draw();
              }
              // Reco jets
              for (int ijet = 0; ijet < jet_pt_calib[ir]->size(); ijet++) {
                double eta = jet_eta_det[ir]->at(ijet);
                double phi = jet_phi[ir]->at(ijet);
                double pt  = jet_pt_calib[ir]->at(ijet);

                TLatex *lab = new TLatex(eta + 0.03, phi + 0.03,
                    Form("%.0f GeV", pt));
                lab->SetTextColor(kRed+2);
                lab->SetTextSize(0.025);
                lab->Draw();
              }

              // Truth jets
              for (int ijet = 0; ijet < TRUTH_pt[ir]->size(); ijet++) {
                double eta = TRUTH_eta[ir]->at(ijet);
                double phi = TRUTH_phi[ir]->at(ijet);
                double pt  = TRUTH_pt[ir]->at(ijet);

                TLatex *lab = new TLatex(eta + 0.03, phi + 0.03,
                    Form("%.0f GeV", pt));
                lab->SetTextColor(kGreen+3);
                lab->SetTextSize(0.025);
                lab->Draw();
              }

              leg->Draw();

              c->Modified();
              c->Update();
              c->SaveAs(Form("pdfs/unmatched_event_display_%s.pdf", sim));

              icount++;
              delete h;
            }
            passed_events[ir]++;
          }

          double fitValue = myFit[ir][j]->Eval(currentPt->at(Idx_Jets[0]));
          double zvtxValue = zvtxRatio[ir]->GetBinContent(zvtxRatio[ir]->FindBin(zvtx));
          double w_ratio = (isMC ? zvtxValue * fitValue * weights[i] : weights[i]);

          LeadingPT[ir][i][j] = currentPt->at(Idx_Jets[0]);
          SLPT[ir][i][j]      = currentPt->at(Idx_Jets[1]);
          SSLPT[ir][i][j]     = currentPt->at(Idx_Jets[2]);

          SLeta[ir][i][j]  = jet_eta_det[ir]->at(Idx_Jets[1]);
          SLphi[ir][i][j]  = jet_phi[ir]->at(Idx_Jets[1]);
          SSLeta[ir][i][j] = jet_eta_det[ir]->at(Idx_Jets[2]);
          SSLphi[ir][i][j] = jet_phi[ir]->at(Idx_Jets[2]);

          PT23[ir][i][j] = newCoordinates(SLPT[ir][i][j], SLphi[ir][i][j], SSLPT[ir][i][j], SSLphi[ir][i][j], false);
          if (i > 0 && j == 0) {
            dphi1_23[ir]->Fill(dPhi123, w_ratio);
          }
          if (i == 0 && j == 0) {
            dphi1_23_data[ir]->Fill(dPhi123);
          }
          weight[ir][i][j] = w_ratio;

          for (int k = 0; k < (int)pTBins.size()-1; k++){
            if (LeadingPT[ir][i][j] > pTBins[k] && LeadingPT[ir][i][j] < pTBins[k+1]){
              if (!isMC) hxj[ir][k][0][0]->Fill(LeadingPT[ir][i][j] / PT23[ir][i][j], w_ratio);
              else       hxj[ir][k][j][1]->Fill(LeadingPT[ir][i][j] / PT23[ir][i][j], w_ratio);
            }
          }
          if (!isMC) leadingJetPT[ir]->Fill(          LeadingPT[ir][i][j], weights[i]); // This is used to calculate w_ratio, so we don't want to use it here
          if (isMC)  leadingJetPT_Pyth[ir][j]->Fill(  LeadingPT[ir][i][j], weights[i]); // This is used to calculate w_ratio, so we don't want to use it here
          outtree[ir][i][j]->Fill();
        }
      }
    }
    std::cout << std::endl;
    if (i > 0) {
      for (int ir = 0; ir < nRadii; ir++) {
        std::cout << "  r0" << radii[ir] << " fraction of matched events: " << matched_events[ir] << " / " << passed_events[ir]
                   << " = " << (passed_events[ir] ? (float)matched_events[ir]/(float)passed_events[ir] : 0.f) << std::endl;
      }
    }
    if (i == 1) c->SaveAs(Form("pdfs/unmatched_event_display_%s.pdf]", sim));
  }

  //////////////////////// WRITING ////////////////////////

  std::cout << "Writing " << wf->GetName() << std::endl;
  wf->Write();
  wf->Close();

  std::cout << "done" << std::endl;
  return 0;

}

