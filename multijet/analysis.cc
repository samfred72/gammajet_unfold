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
#include "TSystem.h"
#include "TRandom3.h"
#include <vector>
#include <cmath>
#include <iomanip>
#include "style.h"

// Multijet balance analysis (from SDCC multiJet_legacy/analysis.cc). Reads the multiJet trees
// (../trees/multijet_*.root): Data multijet_Data.root, MC multijet_<sim>_Jet{8,12,20,30}.root.
// Build: ./make.sh    Run: ./analysis <0=pythia|1=herwig> [tree dir] [--no-jet8] [--tight] [--truth-smear]
//   --no-jet5, --no-jet8  leave out that sample (the lowest sample used starts at truth pT 0)
//   --truth-smear  MC jet pT from the truth-anchored JER smears (jet_pt_smear_{,high_,low_}truth)
//   --no-smear     MC jet pT without the extra JER smearing (jet_pt_calib, in all three slots; a test)
//   --recoil-smear truth2|all|hybrid|none   where jets 2+3 are added unsmeared and the sum smeared once:
//                  truth2 (default) MC events whose jets 2 and 3 belong to one truth jet (see above), width
//                  at that truth jet's pT; all: every MC event, width at the unsmeared recoil pT; hybrid:
//                  every MC event, width at the truth jet's pT when jets 2 and 3 belong to one, else at the
//                  unsmeared recoil pT; none: never
//   --truth-jet-min X  truth-jet pT threshold for that matching (default 7 GeV)
//   --no-recoil-cut  drop the recoil >= 14 GeV cut
//   --tight    test cut: leading > 25 GeV; writes multijet_analysis_<sim>_tight.root
// Writes multijet_analysis_<sim>.root; reads the MC reweighting fits aux/ratio<R>_<sim>.root
// (from makeratio.C on a previous pass) and the JER smearing template aux/jer_smear_templates.root.
//
// Recoil = |pT,2 + pT,3| (vector sum of the subleading and subsubleading jets). Selection: leading jet
// >= 20 GeV, jets 2 and 3 >= 7 GeV each (before any recoil smearing, so the selection stays inside both
// skims), recoil >= 14 GeV, |eta| < 1.1-R for the three jets, dphi12 > 3pi/4, dphi13 > pi/2. In MC the
// leading reco jet is also capped per sample (Jet12 <= 35 GeV, Jet20 <= 50 GeV).
// MC jets are JER-smeared. When jets 2 and 3 (ranked by unsmeared pT) belong to one truth jet, they are
// added unsmeared and the sum is smeared once at that truth jet's pT, with the same template; otherwise
// every jet is smeared on its own. "Belong to one truth jet" (truth jets >= 7 GeV, match radius 0.75R):
//   the leading jet matches a truth jet TL; jet 2 matches a different truth jet TR; jet 3 matches TR or no
//   truth jet (never TL); every other truth jet is more than 1.5R from jets 2 and 3; and the summed
//   recoil points within 0.75R of TR.
// Trees from the current treemaker carry the reco-truth match (jet_truth_pt) and a per-event recoil
// deviate (recoil_smear_z), used here when present; older trees fall back to a dR match and a drawn deviate.

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
    std::cerr << "Usage: " << argv[0] << " <0=pythia|1=herwig> [tree dir, default ../trees] [--no-jet5] [--no-jet8] [--tight] [--truth-smear] [--no-smear] [--recoil-smear truth2|all|hybrid|none] [--no-recoil-cut] [--truth-jet-min X]" << std::endl;
    return 1;
  }
  bool isherwig = (strcmp(argv[1], "1") == 0);
  const char * treedir = "../trees";
  bool useJet5 = true, useJet8 = true, tight = false, truthSmear = false, noSmear = false, recoilCut = true;
  std::string recoilMode = "truth2";
  float truthJetMin = 7;   // truth jets used in the recoil matching
  for (int a = 2; a < argc; a++) {
    if (!strcmp(argv[a], "--no-jet5")) useJet5 = false;
    else if (!strcmp(argv[a], "--no-jet8")) useJet8 = false;
    else if (!strcmp(argv[a], "--tight")) tight = true;
    else if (!strcmp(argv[a], "--truth-smear")) truthSmear = true;
    else if (!strcmp(argv[a], "--no-smear")) noSmear = true;
    else if (!strcmp(argv[a], "--no-recoil-cut")) recoilCut = false;
    else if (!strcmp(argv[a], "--recoil-smear") && a+1 < argc) recoilMode = argv[++a];
    else if (!strcmp(argv[a], "--truth-jet-min") && a+1 < argc) truthJetMin = atof(argv[++a]);
    else if (argv[a][0] != '-') treedir = argv[a];
    else { std::cerr << "Unknown option " << argv[a] << std::endl; return 1; }
  }
  const char * sim = (isherwig ? "herwig" : "pythia");
  const char * tag = (tight ? "_tight" : "");
  std::cout << "sim is: " << argv[1] << " " << sim << (useJet8 ? "" : ", no Jet8") << (tight ? ", tight cuts" : "") << (truthSmear ? ", truth-anchored smearing" : "") << (noSmear ? ", no extra MC smearing" : "")
            << ", recoil smearing " << recoilMode << (recoilCut ? "" : ", no recoil cut")
            << ", truth jets >= " << truthJetMin << " GeV" << std::endl;
  if (recoilMode != "truth2" && recoilMode != "all" && recoilMode != "hybrid" && recoilMode != "none") {
    std::cerr << "--recoil-smear must be truth2, all, hybrid or none" << std::endl;
    return 1;
  }
  const char * smearKind = truthSmear ? "truth" : "reco";

  // All 7 radii in one pass.
  const int nRadii = 7;
  const int radii[nRadii] = {2,3,4,5,6,7,8};
  const int reprRadiusIdx = 2; // radius 4, for the debug event display

  // Reco-based JER variants only: RECO central, HIGH/LOW its envelope.
  const int nsystypes = 3;
  const char * sysNames[nsystypes] = {"RECO","HIGH","LOW"};
  const char * fitFuncNames[nsystypes] = {"ratio_func_JERreco","ratio_func_JERhigh","ratio_func_JERlow"};

  // Selection. Data trees (not the histograms) keep events down to the lowest in-situ JES the
  // scans try, jesFloor: a jet passes there if pT/jesFloor passes the nominal cut.
  const float leading_pT_Cutoff = tight ? 25 : 20;
  const float recoilJetMin = 7;  // jets 2 and 3, each
  const float recoilPtMin = recoilCut ? 14 : 0;  // |pT,2 + pT,3|
  const float SSLDPHI = M_PI/2.0;
  const float SLDPHI = 3*M_PI/4.0;
  const float jesFloor = 0.9;

  // MC samples (index i-1; i = 0 is Data). Cross sections (pb) from the MDC2 table; herwig also
  // corrects for the generated event counts. Weights are normalized to Jet30.
  const int nsimtrees = 5;
  const int ntrees = nsimtrees + 1;
  const char * treenames[ntrees] = {"data", "Jet5", "Jet8", "Jet12", "Jet20", "Jet30"};
  const double xsec[2][nsimtrees] = {
    {1.3878e+08, 1.3013e+07, 1.4903e+06, 6.2623e+04, 2.5298e+03}, // pythia
    {1.8437e+08, 0, 1.132355e+06, 5.2613e+04, 2.0694e+03},        // herwig (no Jet8 sample)
  };
  const double neventCorr[2][nsimtrees] = {
    {1, 1, 1, 1, 1},
    {1, 1, 10000000.0/10001000.0, 10000000.0/10913000.0, 10000000.0/9999000.0},
  };
  // Highest leading reco jet pT (GeV) accepted from each sample (index = tree index, 0 = Data): keeps the
  // heavily weighted low-pT-hat samples out of the bins the higher samples cover.
  const float leadRecoMax[ntrees] = {1e9, 1e9, 1e9, 35, 50, 1e9};

  // Leading-truth-jet pT (GeV) above which each sample is fully efficient, per radius (MDC2 table).
  // Jet8 is not in that table: it can only be combined with Jet5 once these are filled in.
  const double truthThreshold[nsimtrees][nRadii] = {
    { 5,  6,  7, 10, 12, 14, 15}, // Jet5
    {-1, -1, -1, -1, -1, -1, -1}, // Jet8 (unknown)
    {12, 13, 14, 19, 22, 24, 25}, // Jet12
    {20, 21, 21, 27, 29, 32, 34}, // Jet20
    {30, 31, 32, 38, 41, 45, 47}, // Jet30
  };

  TFile * infile = TFile::Open(Form("%s/multijet_Data.root", treedir), "READ");
  if (!infile || infile->IsZombie()) return 1;
  TFile * simfile[nsimtrees] = {nullptr};
  for (int k = 0; k < nsimtrees; k++) {
    const char * name = Form("%s/multijet_%s_%s.root", treedir, sim, treenames[k+1]);
    bool required = (k >= 2); // Jet12, Jet20, Jet30
    const bool excluded = (k == 0 && !useJet5) || (k == 1 && !useJet8);
    if (excluded || gSystem->AccessPathName(name)) {
      if (required) { std::cerr << "Missing " << name << std::endl; return 1; }
      std::cout << "Not using " << treenames[k+1] << (excluded ? " (excluded)" : " (no file)") << std::endl;
      continue;
    }
    simfile[k] = TFile::Open(name, "READ");
    if (!simfile[k] || simfile[k]->IsZombie()) return 1;
  }
  if (simfile[0] && simfile[1]) {
    std::cerr << "Jet5 and Jet8 together need Jet8's truth-pT thresholds (truthThreshold); run with --no-jet8" << std::endl;
    return 1;
  }

  const int isim = isherwig ? 1 : 0;
  double weights[ntrees] = {1};
  for (int k = 0; k < nsimtrees; k++) weights[k+1] = xsec[isim][k] / xsec[isim][nsimtrees-1] * neventCorr[isim][k];

  // Truth-pT stitching windows: each used sample from its threshold up to the next used sample's;
  // the lowest used sample starts at 0 (nothing below it).
  double lowcuts[nRadii][ntrees];
  double highcuts[nRadii][ntrees];
  int firstMC = -1;
  for (int ir = 0; ir < nRadii; ir++) {
    lowcuts[ir][0] = 0; highcuts[ir][0] = 1e9;
    int prev = -1;
    for (int k = 0; k < nsimtrees; k++) {
      if (!simfile[k]) continue;
      lowcuts[ir][k+1] = (prev < 0) ? 0 : truthThreshold[k][ir];
      highcuts[ir][k+1] = 1e9;
      if (prev >= 0) highcuts[ir][prev+1] = truthThreshold[k][ir];
      if (prev < 0) firstMC = k+1;
      prev = k;
    }
  }
  for (int k = 0; k < nsimtrees; k++) {
    if (!simfile[k]) continue;
    std::cout << treenames[k+1] << ": weight " << weights[k+1] << ", R=0.4 truth window ["
              << lowcuts[2][k+1] << ", " << highcuts[2][k+1] << ")" << std::endl;
  }

  // per-radius reweighting fits (makeratio.C)
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

  // JER smearing template (nominal, up, down - the RECO/HIGH/LOW variations), as in the treemaker
  TFile * fJER = TFile::Open("aux/jer_smear_templates.root", "READ");
  if (!fJER || fJER->IsZombie()) { std::cerr << "Missing aux/jer_smear_templates.root" << std::endl; return 1; }
  TH1D * jerWidth[nsystypes];
  const char * jerNames[nsystypes] = {"nominal", "sysup", "sysdown"};
  for (int j = 0; j < nsystypes; j++) {
    jerWidth[j] = (TH1D*)fJER->Get(Form("h_jer_smear_r04_pileup_EMfracJES_%s", jerNames[j]));
    jerWidth[j]->SetDirectory(0);
  }
  fJER->Close();
  TRandom3 rnd(20261006);

  // output file: all radii, all systematics
  TFile * wf = TFile::Open(Form("multijet_analysis_%s%s.root", sim, tag), "RECREATE");
  wf->cd();
  gStyle->SetOptStat(0);

  // output histograms
  std::vector<float> pTBins = {20,25,30,35,50};
  const int nPtBins = 4;

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
    zvtx_data[ir]      = new TH1D(Form("zvtx_data_r%d", radius), ";z (cm);counts", 24, -60, 60);
    zvtx_MC[ir]        = new TH1D(Form("zvtx_MC_r%d", radius), ";z (cm);counts", 24, -60, 60);

    for (int k = 0; k < nPtBins; k++) {
      hxj[ir][k][0][0] = new TH1D(Form("hxj_r%d_%1.0f_data", radius, pTBins[k]), ";x_{j};#frac{1}{N}#frac{dN}{dx_j}", 45, 0.4, 2.65);
      hxj[ir][k][0][0]->Sumw2();
      for (int j = 0; j < nsystypes; j++) {
        hxj[ir][k][j][1] = new TH1D(Form("hxj_r%d_%1.0f_%s_sim", radius, pTBins[k], sysNames[j]), ";x_{j};#frac{1}{N}#frac{dN}{dx_j}", 45, 0.4, 2.65);
        hxj[ir][k][j][1]->Sumw2();
      }
    }
  }

  // output trees (one per MC sample slot, empty if the sample isn't used)
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

  // tree setup
  TTree * intree[ntrees] = {(TTree*) infile->Get("ttree")};
  for (int k = 0; k < nsimtrees; k++) intree[k+1] = simfile[k] ? (TTree*) simfile[k]->Get("ttree") : nullptr;

  TCanvas * c = new TCanvas("c_unmatched","",500,1000);
  c->SaveAs(Form("pdfs/unmatched_event_display_%s%s.pdf[", sim, tag));
  int icount = 0;

  std::cout << "Running over trees..." << std::endl;
  for (int i = 0; i < ntrees; i++) {
    std::cout << i << "..." << std::endl;
    if (!intree[i]) continue;
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
    std::vector<float>* jet_eta[nRadii]      = {nullptr};
    std::vector<float>* jet_phi[nRadii]      = {nullptr};
    std::vector<float>* jet_truth_pt[nRadii] = {nullptr};
    float recoil_z[nRadii] = {0};
    bool hasMatch = false, hasRecoilZ = false;

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
        intree[i]->SetBranchAddress(Form("jet_pt_smear_%s_%d", smearKind, radius), &jet_pt_smearRECO[ir]);
        intree[i]->SetBranchAddress(Form("jet_pt_smear_high_%s_%d", smearKind, radius), &jet_pt_smearHIGH[ir]);
        intree[i]->SetBranchAddress(Form("jet_pt_smear_low_%s_%d", smearKind, radius), &jet_pt_smearLOW[ir]);
        intree[i]->SetBranchAddress(Form("jet_eta_%d", radius), &jet_eta[ir]);
        hasMatch = intree[i]->GetBranch(Form("jet_truth_pt_%d", radius));
        hasRecoilZ = intree[i]->GetBranch(Form("recoil_smear_z_%d", radius));
        if (hasMatch) intree[i]->SetBranchAddress(Form("jet_truth_pt_%d", radius), &jet_truth_pt[ir]);
        if (hasRecoilZ) intree[i]->SetBranchAddress(Form("recoil_smear_z_%d", radius), &recoil_z[ir]);
      }
      intree[i]->SetBranchAddress(Form("jet_pt_calib_%d", radius), &jet_pt_calib[ir]);
      intree[i]->SetBranchAddress(Form("jet_eta_det_%d", radius), &jet_eta_det[ir]);
      intree[i]->SetBranchAddress(Form("jet_phi_%d", radius), &jet_phi[ir]);
    }
    if (!isMC) {
      intree[i]->SetBranchAddress("calib_lead_time", &calib_lead_time);
      intree[i]->SetBranchAddress("calib_delta_time", &calib_delta_time);
    }
    intree[i]->SetBranchAddress("mbd_vertex_z", &zvtx);
    intree[i]->SetBranchAddress("mbd_hit", &mbd_hit);

    // tree analysis
    int matched_events[nRadii] = {0};
    double nSelected[nRadii] = {0}, nRecoilSmeared[nRadii] = {0}; // selected RECO events, unweighted
    int passed_events[nRadii] = {0};

    Long64_t nentries = intree[i]->GetEntries();
    for (Long64_t e = 0; e < nentries; e++) {
      intree[i]->GetEntry(e);
      if (e % 1000 == 0) std::cout << "entry " << e << "/" << nentries << " (" << (float)e/nentries*100. << "%)" << "\t\r" << std::flush;

      // event-level cuts (radius-independent)

      if (std::fabs(zvtx) > 60) continue;

      if (!isMC) {
        // timing cut
        double x = calib_lead_time;
        double y = x - calib_delta_time;
        if (!(fabs(x) < 6 && fabs(y) < 3)) continue; //requiring timing cut
      }

      for (int ir = 0; ir < nRadii; ir++) {
        int radius = radii[ir];
        float R = radius / 10.0f;

        if (jet_pt_calib[ir]->size() <= 2) continue; // requiring three jets

        if (isMC) {
          // truth selection
          if (TRUTH_pt[ir]->size() == 0) continue; // requiring at least one truth jet

          std::vector<int> IDXTRUTH = IDXGrab(TRUTH_pt[ir]->size(), *TRUTH_pt[ir]);
          if (lowcuts[ir][i] > TRUTH_pt[ir]->at(IDXTRUTH[0]) ||
              TRUTH_pt[ir]->at(IDXTRUTH[0]) >= highcuts[ir][i]) continue;
          leadingJetPT_truth[ir]->Fill(TRUTH_pt[ir]->at(IDXTRUTH[0]), weights[i]);
        }
        std::vector<int> truthJets;
        if (isMC) for (size_t t = 0; t < TRUTH_pt[ir]->size(); t++) if (TRUTH_pt[ir]->at(t) >= truthJetMin) truthJets.push_back(t);
        // Do jets 2 and 3 (ranked by unsmeared pT) belong to one truth jet? (see the header comment)
        int recoilTruth = -1;
        if (isMC && !noSmear && truthJets.size() >= 2 && (recoilMode == "truth2" || recoilMode == "hybrid")) {
          std::vector<int> c = IDXGrab(jet_pt_calib[ir]->size(), *jet_pt_calib[ir]);
          if (c[2] != -1) {
            const float rMatch = 0.75*R, rIso = 1.5*R;
            auto dRj = [&](int jet, int t) {
              float de = jet_eta[ir]->at(jet) - TRUTH_eta[ir]->at(t), dp = deltaPhi(jet_phi[ir]->at(jet), TRUTH_phi[ir]->at(t));
              return std::sqrt(de*de + dp*dp);
            };
            auto nearest = [&](int jet, int exclude) {
              int best = -1; float bestDR = rMatch;
              for (int t : truthJets) if (t != exclude && dRj(jet, t) < bestDR) { bestDR = dRj(jet, t); best = t; }
              return best;
            };
            // leading jet: the treemaker's match when it is one of the truth jets, else the nearest within rMatch
            int tL = -1;
            const float leadMatch = hasMatch ? jet_truth_pt[ir]->at(c[0]) : -1;
            for (int t : truthJets) if (leadMatch > 0 && std::fabs(TRUTH_pt[ir]->at(t) - leadMatch) < 1e-3) tL = t;
            if (tL < 0) tL = nearest(c[0], -1);
            const int tR = (tL >= 0) ? nearest(c[1], tL) : -1;
            bool ok = (tL >= 0 && tR >= 0);
            if (ok) {
              const int t3 = nearest(c[2], -1); // jet 3: on the recoil truth jet, or on no truth jet
              if (dRj(c[2], tL) < rMatch || (t3 >= 0 && t3 != tR && dRj(c[2], tR) >= rMatch)) ok = false;
              for (int t : truthJets)               // other truth jets well away from jets 2 and 3
                if (t != tL && t != tR && (dRj(c[1], t) < rIso || dRj(c[2], t) < rIso)) ok = false;
            }
            if (ok) {                               // the summed recoil points at the recoil truth jet
              double px = 0, py = 0, pz = 0;
              for (int k : {c[1], c[2]}) {
                const double pt = jet_pt_calib[ir]->at(k), eta = jet_eta[ir]->at(k), phi = jet_phi[ir]->at(k);
                px += pt*std::cos(phi); py += pt*std::sin(phi); pz += pt*std::sinh(eta);
              }
              const double ptSum = std::hypot(px, py);
              if (ptSum > 0) {
                const double de = std::asinh(pz/ptSum) - TRUTH_eta[ir]->at(tR);
                const double dp = deltaPhi(std::atan2(py, px), TRUTH_phi[ir]->at(tR));
                if (std::sqrt(de*de + dp*dp) < rMatch) recoilTruth = tR;
              }
            }
          }
        }
        const bool recoilSmear = isMC && !noSmear &&
                                 (recoilMode == "all" || recoilMode == "hybrid" || recoilTruth >= 0);
        // one recoil deviate per event and radius, shared by the RECO/HIGH/LOW variations
        const float recoilZ = !recoilSmear ? 0.f : hasRecoilZ ? recoil_z[ir] : (float)rnd.Gaus(0, 1);

        // loop over RECO / HIGH / LOW smeared pT

        std::vector<std::vector<float>*> ptVariations = (!isMC ?
            std::vector<std::vector<float>*>{ jet_pt_calib[ir] } : noSmear ?
            std::vector<std::vector<float>*>{ jet_pt_calib[ir], jet_pt_calib[ir], jet_pt_calib[ir] } :
            std::vector<std::vector<float>*>{
              jet_pt_smearRECO[ir],
              jet_pt_smearHIGH[ir],
              jet_pt_smearLOW[ir]
            });

        for (int j = 0; j < nJloop; j++) {

          std::vector<float>* currentPt = ptVariations[j];
          std::vector<float>* rankPt = recoilSmear ? jet_pt_calib[ir] : currentPt;
          std::vector<int> Idx_Jets = IDXGrab(rankPt->size(), *rankPt);

          if (Idx_Jets[0] == -1 || Idx_Jets[1] == -1 || Idx_Jets[2] == -1) continue;

          const float lead = currentPt->at(Idx_Jets[0]);
          float pt2 = currentPt->at(Idx_Jets[1]), pt3 = currentPt->at(Idx_Jets[2]);
          if (recoilSmear) {
            pt2 = jet_pt_calib[ir]->at(Idx_Jets[1]);
            pt3 = jet_pt_calib[ir]->at(Idx_Jets[2]);
          }
          const float cut2 = pt2, cut3 = pt3; // the per-jet cuts use the values before any recoil smearing
          if (recoilSmear) {
            // add jets 2 and 3 unsmeared, smear the sum once at the recoil truth jet's pT, and scale both
            // jets by the same factor so their vector sum is the smeared recoil
            const float sum = newCoordinates(pt2, jet_phi[ir]->at(Idx_Jets[1]), pt3, jet_phi[ir]->at(Idx_Jets[2]), false);
            // width at the matched recoil truth jet's pT, or at the unsmeared recoil pT ("all", unmatched "hybrid")
            const float ptTruth = recoilTruth >= 0 ? TRUTH_pt[ir]->at(recoilTruth) : sum;
            const float width = jerWidth[j]->Interpolate(std::min(std::max(ptTruth, 5.01f), 79.9f));
            const float smeared = sum + recoilZ * ptTruth * width;
            if (sum <= 0 || smeared <= 0) continue;
            pt2 *= smeared / sum;
            pt3 *= smeared / sum;
          }

          float dPhi13 = deltaPhi(
              jet_phi[ir]->at(Idx_Jets[0]),
              jet_phi[ir]->at(Idx_Jets[2])
              );

          float dPhi12 = deltaPhi(
              jet_phi[ir]->at(Idx_Jets[0]),
              jet_phi[ir]->at(Idx_Jets[1])
              );
          float phi23 = newCoordinates(pt2, jet_phi[ir]->at(Idx_Jets[1]), pt3, jet_phi[ir]->at(Idx_Jets[2]), true);
          const float recoil = newCoordinates(pt2, jet_phi[ir]->at(Idx_Jets[1]), pt3, jet_phi[ir]->at(Idx_Jets[2]), false);
          float dPhi123 = deltaPhi(
              jet_phi[ir]->at(Idx_Jets[0]),
              phi23
              );

          const float floorScale = isMC ? 1.f : jesFloor;
          if (leading_pT_Cutoff*floorScale > lead ||
              lead > leadRecoMax[i] ||
              recoilJetMin*floorScale > cut2 ||
              recoilJetMin*floorScale > cut3 ||
              recoilPtMin*floorScale > recoil ||
              std::fabs(jet_eta_det[ir]->at(Idx_Jets[0])) > (1.1-R) ||
              std::fabs(jet_eta_det[ir]->at(Idx_Jets[1])) > (1.1-R) ||
              std::fabs(jet_eta_det[ir]->at(Idx_Jets[2])) > (1.1-R) ||
              std::fabs(dPhi13) < SSLDPHI ||
              std::fabs(dPhi12) < SLDPHI
             ) continue;
          // nominal selection (the loose Data events above only go to the output tree)
          const bool passNominal = lead >= leading_pT_Cutoff && cut2 >= recoilJetMin && cut3 >= recoilJetMin && recoil >= recoilPtMin;

          if (passNominal && i == firstMC && j == 0) {
            dphi12[ir]->Fill(dPhi12);
            dphi13[ir]->Fill(dPhi13);
            heta0[ir]->Fill(jet_eta_det[ir]->at(Idx_Jets[0]));
            heta1[ir]->Fill(jet_eta_det[ir]->at(Idx_Jets[1]));
            heta2[ir]->Fill(jet_eta_det[ir]->at(Idx_Jets[2]));
          }

          if (isMC && j == 0 && !recoilSmear) {
            bool match0 = (std::fabs(currentPt->at(Idx_Jets[0]) - jet_pt_calib[ir]->at(Idx_Jets[0])) > 0.00001);
            bool match1 = (std::fabs(currentPt->at(Idx_Jets[1]) - jet_pt_calib[ir]->at(Idx_Jets[1])) > 0.00001);
            bool match2 = (std::fabs(currentPt->at(Idx_Jets[2]) - jet_pt_calib[ir]->at(Idx_Jets[2])) > 0.00001);

            if (match0 && match1 && match2) { matched_events[ir]++; }
            else if (ir == reprRadiusIdx && i == firstMC && icount < 100) { // bad-event display (radius 4, lowest MC sample)
              c->cd();

              TH2D *h = new TH2D(Form("heventdisplay_%d",icount), ";#eta;#phi", 100,-1.5,1.5, 100,-M_PI,M_PI);
              h->Draw();

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

              for (int ijet = 0; ijet < jet_pt_calib[ir]->size(); ijet++) {
                double eta = jet_eta_det[ir]->at(ijet);
                double phi = jet_phi[ir]->at(ijet);
                double pt  = jet_pt_calib[ir]->at(ijet);

                TMarker *mjet = new TMarker(eta, phi, 107);
                mjet->SetMarkerColor(kRed);
                mjet->SetMarkerSize(15);
                mjet->Draw();
              }

              for (int ijet = 0; ijet < TRUTH_pt[ir]->size(); ijet++) {
                double eta = TRUTH_eta[ir]->at(ijet);
                double phi = TRUTH_phi[ir]->at(ijet);

                TMarker *mjet = new TMarker(eta, phi, 107);
                mjet->SetMarkerColor(kGreen+2);
                mjet->SetMarkerSize(15);
                mjet->Draw();
              }
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
              c->SaveAs(Form("pdfs/unmatched_event_display_%s%s.pdf", sim, tag));

              icount++;
              delete h;
            }
            passed_events[ir]++;
          }

          double fitValue = myFit[ir][j]->Eval(lead);
          double zvtxValue = zvtxRatio[ir]->GetBinContent(zvtxRatio[ir]->FindBin(zvtx));
          double w_ratio = (isMC ? zvtxValue * fitValue * weights[i] : weights[i]);

          LeadingPT[ir][i][j] = lead;
          SLPT[ir][i][j]      = pt2;
          SSLPT[ir][i][j]     = pt3;

          SLeta[ir][i][j]  = jet_eta_det[ir]->at(Idx_Jets[1]);
          SLphi[ir][i][j]  = jet_phi[ir]->at(Idx_Jets[1]);
          SSLeta[ir][i][j] = jet_eta_det[ir]->at(Idx_Jets[2]);
          SSLphi[ir][i][j] = jet_phi[ir]->at(Idx_Jets[2]);

          PT23[ir][i][j] = newCoordinates(SLPT[ir][i][j], SLphi[ir][i][j], SSLPT[ir][i][j], SSLphi[ir][i][j], false);
          weight[ir][i][j] = w_ratio;
          outtree[ir][i][j]->Fill();
          if (!passNominal) continue;
          if (isMC && j == 0) { nSelected[ir]++; if (recoilSmear) nRecoilSmeared[ir]++; }

          if (i > 0 && j == 0) {
            dphi1_23[ir]->Fill(dPhi123, w_ratio);
          }
          if (i == 0 && j == 0) {
            dphi1_23_data[ir]->Fill(dPhi123);
          }

          for (int k = 0; k < (int)pTBins.size()-1; k++){
            if (LeadingPT[ir][i][j] > pTBins[k] && LeadingPT[ir][i][j] < pTBins[k+1]){
              if (!isMC) hxj[ir][k][0][0]->Fill(LeadingPT[ir][i][j] / PT23[ir][i][j], w_ratio);
              else       hxj[ir][k][j][1]->Fill(LeadingPT[ir][i][j] / PT23[ir][i][j], w_ratio);
            }
          }
          // Inputs to makeratio.C's reweighting, after the full selection and with the
          // cross-section weights only: each correction is derived from the unreweighted MC.
          if (!isMC) leadingJetPT[ir]->Fill(LeadingPT[ir][i][j]);
          if (isMC)  leadingJetPT_Pyth[ir][j]->Fill(LeadingPT[ir][i][j], weights[i]);
          if (j == 0) (isMC ? zvtx_MC[ir] : zvtx_data[ir])->Fill(zvtx, weights[i]);
        }
      }
    }
    std::cout << std::endl;
    if (i > 0) {
      for (int ir = 0; ir < nRadii; ir++) {
        std::cout << "  r0" << radii[ir] << " fraction of matched events: " << matched_events[ir] << " / " << passed_events[ir]
                   << " = " << (passed_events[ir] ? (float)matched_events[ir]/(float)passed_events[ir] : 0.f)
                   << "; recoil-smeared: " << nRecoilSmeared[ir] << " / " << nSelected[ir] << " selected" << std::endl;
      }
    }
    if (i == firstMC) c->SaveAs(Form("pdfs/unmatched_event_display_%s%s.pdf]", sim, tag));
  }

  // writing

  std::cout << "Writing " << wf->GetName() << std::endl;
  wf->Write();
  wf->Close();

  std::cout << "done" << std::endl;
  return 0;

}

