// Draws the reco-level xJ shape (Data vs. MC, one page per reported photon-pT bin) using
// only the "Everything new" reconstruction - i.e. compare_old_new.C's case 5
// (cluster_pt/eta/phi, cluster_bdt_scores[9], cluster_showershape[11], and the "new" jet
// pt: jet_pt_smear_truth[ir] for MC / jet_pt_calib[ir] for Data) - NOT the purity-corrected
// case 6 there. Same event selection (vz cut, truth-pt trigger-stitching window, xJ floor,
// |eta| cuts, dphi cut, ABCD region-A cut), same representative jet radius (ir=2, R=0.4),
// and same cross-section/vz/cluster-pT reweighting as compare_old_new.C - see that file's
// header comment for the full rationale; this macro only pulls out the one case needed for
// a plain Data/MC shape comparison, so it doesn't carry the other five cases' bookkeeping.
//
// Both histograms are unit-normalized (independently) before drawing, since Data and MC
// have unrelated absolute normalizations (weighted cross-section-scaled MC vs. raw Data
// event counts) - this is a shape comparison, not an absolute-rate one.

#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/reweight_utility.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"
R__LOAD_LIBRARY(libgammajet_unfold.so)

#include <string>
#include <map>
#include <utility>
#include <cmath>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"

const int ir = 2; // R=0.4, representative radius - same convention as compare_old_new.C

// truth-level cluster-pT trigger-stitching window, pythia sim, from
// treeuser.cc's threshmap[-1]/threshmap_high[-1] (Photon5/10/20 rows). MC only.
std::map<std::string, std::pair<float, float>> truthPtWindow = {
  {"Photon5",  {0.f,  12.f}},
  {"Photon10", {12.f, 24.f}},
  {"Photon20", {24.f, 100.f}}
};

// Per-sample cross-section/luminosity normalization - literal copy of drawer.h's
// scalemap[/*isphoton=*/1][sample] for sim=="pythia" (src/drawer.h's `drawer` ctor).
std::map<std::string, double> sampleWeight = {
  {"Photon5",  146359.3},
  {"Photon10", 6944.675},
  {"Photon20", 130.4461}
};

// Exact unbinned mean/error accumulator - same as compare_old_new.C's Accum. Needed
// because the displayed histograms get width-scaled for the density plot (see the
// normalization comment below), and TH1::GetMean() on a width-scaled histogram is NOT
// the event-weighted mean of x_J (it implicitly re-weights bins by 1/width) - so the
// legend's <x_J> must come from this, not from the displayed histogram.
struct Accum {
  double sumw = 0, sumw2 = 0, sumwx = 0, sumwx2 = 0;
  Long64_t n = 0;
  void fill(double x, double w) {
    sumw += w; sumw2 += w * w; sumwx += w * x; sumwx2 += w * x * x; n++;
  }
  double mean() const { return sumw > 0 ? sumwx / sumw : 0; }
};

// Mirrors unfolder::check_pair (src/unfolder.cc) at fixed ir=2, floorScale=1, testPt=-1.
bool checkPair(float jetPt, float jetEta, float phoPt, float phoEta, float dphi) {
  int ptbin = ana::findPtBin(phoPt);
  if (ptbin == -1) return false;
  float val = jetPt / phoPt;
  float lowval = ana::jet_calib_pt_cut[ir] / ana::ptBins[ptbin];
  float lowbin = ana::unfoldXjBins[ana::findUnfoldXjBin(lowval) + 1];
  if (val < lowbin) return false;
  if (fabs(phoEta) > ana::etacut) return false;
  if (fabs(jetEta) > ana::etacut - ana::JetRs[ir]) return false;
  if (dphi < ana::oppcut) return false;
  return true;
}

float deltaPhi(float p1, float p2) {
  float dphi = fabs(p1 - p2);
  if (dphi > M_PI) dphi = 2 * M_PI - dphi;
  return dphi;
}

// Fills hxjPt[ip] (ana::nPtBinsUsed histograms, ana::unfoldXjBins binning) and accPt[ip]
// (exact unbinned mean) with the "Everything new" case's xJ - see file header comment.
// trigger is one of "Photon5"/"Photon10"/"Photon20" when isMC, or "Data" when not.
void processFile(const std::string & trigger, bool isMC, std::vector<TH1D*> & hxjPt,
    std::vector<Accum> & accPt, Reweighter & rw, Long64_t maxEntries = -1) {

  std::string fname = isMC ?
      "/home/samson72/sphnx/gammajet_unfold/trees/gammajet_pythia_" + trigger + ".root" :
      "/home/samson72/sphnx/gammajet_unfold/trees/gammajet_Data.root";
  TFile * f = TFile::Open(fname.c_str(), "read");
  if (!f || f->IsZombie()) { std::cout << "Could not open " << fname << std::endl; return; }
  TTree * t = (TTree*)f->Get("towerntup");
  Long64_t nentries = t->GetEntries();
  if (maxEntries >= 0 && maxEntries < nentries) nentries = maxEntries;
  std::cout << "Processing " << fname << " (" << nentries << " entries, "
            << (isMC ? "MC" : "Data") << ")" << std::endl;

  t->SetBranchStatus("*", 0);
  std::vector<std::string> branches = {"vz", "cluster_pt", "cluster_eta", "cluster_phi",
       "cluster_showershape", "cluster_bdt_scores", "jet_eta", "jet_phi"};
  if (isMC) { branches.push_back("truth_cluster_pt"); branches.push_back("jet_pt_smear_truth"); }
  else      { branches.push_back("jet_pt_calib"); }
  for (const auto & bn : branches) t->SetBranchStatus(bn.c_str(), 1);

  Float_t vz, cluster_pt, cluster_eta, cluster_phi;
  Float_t cluster_showershape[12], cluster_bdt_scores[11];
  Float_t truth_cluster_pt = 0;
  Float_t jet_pt_smear_truth[ana::nJetR] = {0}, jet_pt_calib[ana::nJetR] = {0};
  Float_t jet_eta[ana::nJetR], jet_phi[ana::nJetR];

  t->SetBranchAddress("vz", &vz);
  t->SetBranchAddress("cluster_pt", &cluster_pt);
  t->SetBranchAddress("cluster_eta", &cluster_eta);
  t->SetBranchAddress("cluster_phi", &cluster_phi);
  t->SetBranchAddress("cluster_showershape", cluster_showershape);
  t->SetBranchAddress("cluster_bdt_scores", cluster_bdt_scores);
  t->SetBranchAddress("jet_eta", jet_eta);
  t->SetBranchAddress("jet_phi", jet_phi);
  if (isMC) {
    t->SetBranchAddress("truth_cluster_pt", &truth_cluster_pt);
    t->SetBranchAddress("jet_pt_smear_truth", jet_pt_smear_truth);
  } else {
    t->SetBranchAddress("jet_pt_calib", jet_pt_calib);
  }

  auto windowIt = truthPtWindow.find(trigger);
  bool haveWindow = isMC && windowIt != truthPtWindow.end();
  double xsecWeight = isMC ? sampleWeight.at(trigger) : 1.0;

  for (Long64_t e = 0; e < nentries; e++) {
    t->GetEntry(e);
    if (e % 2000000 == 0)
      std::cout << "  entry " << e << "/" << nentries
                << " (" << (float)e / nentries * 100. << "%)" << std::endl;

    if (fabs(vz) > ana::vzcut) continue;
    if (haveWindow && !(truth_cluster_pt > windowIt->second.first && truth_cluster_pt < windowIt->second.second))
      continue;

    float mcWeight = (isMC ? rw.GetWeight(vz, cluster_pt) : 1.0f) * xsecWeight;

    float jetEta = jet_eta[ir];
    float jetPhi = jet_phi[ir];
    // "Everything new": cluster_pt/eta/phi (with-saturation), cluster_bdt_scores[9]
    // (new bdt), cluster_showershape[11] (new iso), and the "new" jet pt - see file
    // header comment for the isMC/Data split.
    float phoPt = cluster_pt, phoEta = cluster_eta, phoPhi = cluster_phi;
    float jetPt = isMC ? jet_pt_smear_truth[ir] : jet_pt_calib[ir];
    float iso = cluster_showershape[11];
    float bdt = cluster_bdt_scores[9];

    if (!std::isfinite(phoPt) || phoPt <= 0 || !std::isfinite(jetPt)) continue;
    if (!(phoPt >= ana::ptBins[0] && phoPt < ana::ptBins[ana::nPtBins] && jetPt > ana::jet_calib_pt_cut[ir]))
      continue;
    float dphi = deltaPhi(phoPhi, jetPhi);
    if (!checkPair(jetPt, jetEta, phoPt, phoEta, dphi)) continue;

    int iabcd = ana::findabcdBin(iso, bdt, 0);
    if (iabcd != 0) continue;

    float xj = jetPt / phoPt;
    if (!std::isfinite(xj)) continue;

    int ptbin = ana::findPtBin(phoPt);
    int iused = ptbin - ana::firstUsedPtBin;
    if (iused >= 0 && iused < ana::nPtBinsUsed) {
      hxjPt[iused]->Fill(xj, mcWeight);
      accPt[iused].fill(xj, mcWeight);
    }
  }
  f->Close();
}

void draw_xj_data_mc(Long64_t maxEntriesPerFile = -1) {
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();

  Reweighter rw;

  // drawer::drawAll (used below for the sPHENIX label) only needs its member function,
  // not the per-sample unfolding-output files its constructor opens - see
  // compare_old_new.C's identical comment for why the resulting TFile::Open error spam
  // is silenced here.
  int oldErrLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  drawer d("pythia", "nominal");
  gErrorIgnoreLevel = oldErrLevel;

  int nb = ana::nUnfoldXjBins;
  const double * edges = ana::unfoldXjBins;
  std::string jetFeature = Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV", ana::JetRs[ir], ana::jet_calib_pt_cut[ir]);

  std::vector<TH1D*> hxjMC(ana::nPtBinsUsed), hxjData(ana::nPtBinsUsed);
  std::vector<Accum> accMC(ana::nPtBinsUsed), accData(ana::nPtBinsUsed);
  for (int ip = 0; ip < ana::nPtBinsUsed; ip++) {
    hxjMC[ip] = new TH1D(Form("hxjMC_pt%d", ip), ";x_{J#gamma};normalized counts", nb, edges);
    hxjMC[ip]->SetStatOverflows(TH1::kConsider);
    hxjMC[ip]->SetLineColor(kMagenta + 1);
    hxjMC[ip]->SetLineWidth(2);

    hxjData[ip] = new TH1D(Form("hxjData_pt%d", ip), ";x_{J#gamma};normalized counts", nb, edges);
    hxjData[ip]->SetStatOverflows(TH1::kConsider);
    hxjData[ip]->SetLineColor(kBlue);
    hxjData[ip]->SetLineWidth(2);
  }

  for (const std::string & trigger : {"Photon5", "Photon10", "Photon20"})
    processFile(trigger, /*isMC=*/true, hxjMC, accMC, rw, maxEntriesPerFile);

  std::string dataFname = "/home/samson72/sphnx/gammajet_unfold/trees/gammajet_Data.root";
  bool haveData = !gSystem->AccessPathName(dataFname.c_str()); // AccessPathName returns 0 (false) iff the file exists
  if (haveData) {
    processFile("Data", /*isMC=*/false, hxjData, accData, rw, maxEntriesPerFile);
  } else {
    std::cout << "\nData tree not found at " << dataFname << " - Data histograms will be empty." << std::endl;
  }

  std::cout << "\n" << Form("%-22s %10s %10s %10s", "pT bin", "<xJ> MC", "<xJ> Data", "Data/MC") << std::endl;
  for (int ip = 0; ip < ana::nPtBinsUsed; ip++) {
    double mMC = accMC[ip].mean(), mData = accData[ip].mean();
    std::cout << Form("%-22s %10.4f %10.4f %10.4f",
        Form("%.0f-%.0f GeV", ana::ptBinsUsed[ip], ana::ptBinsUsed[ip + 1]), mMC, mData,
        (mMC > 0 && mData > 0) ? mData / mMC : 0) << std::endl;
  }

  std::string pdfPath = "/home/samson72/sphnx/gammajet_unfold/temporary_study/draw_xj_data_mc.pdf";
  TCanvas * c = new TCanvas("c_xj_data_mc", "", 700, 600);
  c->SetLeftMargin(.13);
  c->SetBottomMargin(.13);
  c->SaveAs((pdfPath + "[").c_str());

  for (int ip = 0; ip < ana::nPtBinsUsed; ip++) {
    c->Clear();
    // ana::unfoldXjBins is variable-width (0.1 wide up to 1.3, then 0.2, then 0.3 near
    // the tail) - Integral() here deliberately omits "width" (it's just the total
    // event/weight count, the correct normalization denominator), but Scale's "width"
    // option is required so the displayed bin content is a density (content/binwidth),
    // not a bare bin fraction - otherwise the wider tail bins would be shown too tall
    // relative to the narrower ones for the same underlying number of events.
    double integMC = hxjMC[ip]->Integral(0, hxjMC[ip]->GetNbinsX() + 1);
    double integData = hxjData[ip]->Integral(0, hxjData[ip]->GetNbinsX() + 1);
    if (integMC > 0) hxjMC[ip]->Scale(1.0 / integMC, "width");
    if (integData > 0) hxjData[ip]->Scale(1.0 / integData, "width");

    double ymax = std::max(hxjMC[ip]->GetMaximum(), hxjData[ip]->GetMaximum());
    hxjMC[ip]->SetTitle(";x_{J#gamma};(1/N) dN/dx_{J#gamma}");
    hxjMC[ip]->GetYaxis()->SetRangeUser(0, ymax * 1.5);
    hxjMC[ip]->Draw("hist");
    hxjData[ip]->Draw("hist same");

    TLegend * leg = new TLegend(0.47, 0.65, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    // <x_J> from the exact unbinned Accum, NOT TH1::GetMean() on the width-scaled
    // histogram above (that would implicitly re-weight bins by 1/width - see Accum's
    // comment).
    leg->AddEntry(hxjMC[ip], Form("Pythia8 #gamma+jet MC  #LTx_{J}#GT=%.3f", accMC[ip].mean()), "l");
    leg->AddEntry(hxjData[ip], Form("p+p Run24 Data  #LTx_{J}#GT=%.3f", accData[ip].mean()), "l");
    leg->Draw();

    d.drawAll({"Pythia8 MC vs. Data"}, {Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",
              ana::ptBinsUsed[ip], ana::ptBinsUsed[ip + 1]), jetFeature}, .17, .87, 14, gPad->GetWh() * 0.8);
    c->SaveAs(pdfPath.c_str());
  }

  c->SaveAs((pdfPath + "]").c_str());
  std::cout << "Wrote " << pdfPath << std::endl;
}
