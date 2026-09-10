// Draws the reco-level xJ shape for Data only, raw vs. purity-corrected, one page per
// reported photon-pT bin - i.e. compare_old_new.C's case 5 ("Everything new": cluster_pt/
// eta/phi, cluster_bdt_scores[9], cluster_showershape[11], jet_pt_calib[ir]) vs. its case 6
// (the same selection, purity-corrected). No MC here - see draw_xj_data_mc.C for the
// Data-vs-MC comparison (raw, uncorrected) this macro's "raw" curve matches exactly.
//
// Purity correction: same two-purity background subtraction as insitu/grid_insitu.C and
// compare_old_new.C's case 6 - case 5's Region A (findabcdBin==0) events combined with its
// Region C (findabcdBin==2: good iso, bad bdt) sibling events per used pT bin, via
// unfold_utility::purityCorrectCoeffs/purityCorrect, with purity P_A/P_C (+ErrorLow/
// ErrorHigh) from ana::getPurity/getPurityC at systag="nominal", ir=2 - i.e. the
// purity_nominal.root inputs, not re-derived here. Same event selection otherwise (vz cut,
// xJ floor, |eta| cuts, dphi cut), same representative jet radius (ir=2, R=0.4) as
// compare_old_new.C/draw_xj_data_mc.C.
//
// Both histograms are unit-normalized as a density (Scale(..., "width")) before drawing,
// since ana::unfoldXjBins is variable-width - see compare_old_new.C/draw_xj_data_mc.C's
// identical normalization comment. The legend's <x_J> is the exact unbinned mean (Accum),
// not TH1::GetMean() on the density-scaled histogram, for the same reason.

#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/unfold_utility.h"
R__LOAD_LIBRARY(libgammajet_unfold.so)

#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

const int ir = 2; // R=0.4, representative radius - same convention as compare_old_new.C

// Exact unbinned mean/error accumulator - same as compare_old_new.C/draw_xj_data_mc.C's
// Accum. Needed because the displayed histograms get width-scaled for the density plot,
// and TH1::GetMean() on a width-scaled histogram is NOT the event-weighted mean of x_J.
struct Accum {
  double sumw = 0, sumw2 = 0, sumwx = 0, sumwx2 = 0;
  Long64_t n = 0;
  void fill(double x, double w) {
    sumw += w; sumw2 += w * w; sumwx += w * x; sumwx2 += w * x * x; n++;
  }
  double mean() const { return sumw > 0 ? sumwx / sumw : 0; }
  double meanErr() const {
    if (sumw <= 0 || n < 2) return 0;
    double var = sumwx2 / sumw - mean() * mean();
    if (var < 0) var = 0;
    double neff = sumw2 > 0 ? sumw * sumw / sumw2 : n;
    return std::sqrt(var / neff);
  }
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

void draw_xj_data_purity(Long64_t maxEntries = -1) {
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();

  // drawer::drawAll (used below for the sPHENIX label) only needs its member function -
  // see compare_old_new.C's identical comment for why the resulting TFile::Open error
  // spam is silenced here.
  int oldErrLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  drawer d("pythia", "nominal");
  gErrorIgnoreLevel = oldErrLevel;

  int nb = ana::nUnfoldXjBins;
  const double * edges = ana::unfoldXjBins;
  std::string jetFeature = Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV", ana::JetRs[ir], ana::jet_calib_pt_cut[ir]);

  std::vector<TH1D*> hRaw(ana::nPtBinsUsed), hCorr(ana::nPtBinsUsed);
  std::vector<TH1D*> hA(ana::nPtBinsUsed), hC(ana::nPtBinsUsed); // Region A/C temps, not drawn
  std::vector<Accum> accRaw(ana::nPtBinsUsed), accA(ana::nPtBinsUsed), accC(ana::nPtBinsUsed), accCorr(ana::nPtBinsUsed);
  for (int ip = 0; ip < ana::nPtBinsUsed; ip++) {
    hRaw[ip] = new TH1D(Form("hxjRaw_pt%d", ip), ";x_{J#gamma};normalized counts", nb, edges);
    hRaw[ip]->SetStatOverflows(TH1::kConsider);
    hRaw[ip]->SetLineColor(kBlue);
    hRaw[ip]->SetLineWidth(2);

    hCorr[ip] = new TH1D(Form("hxjCorr_pt%d", ip), ";x_{J#gamma};normalized counts", nb, edges);
    hCorr[ip]->SetStatOverflows(TH1::kConsider);
    hCorr[ip]->SetLineColor(kOrange + 7);
    hCorr[ip]->SetLineWidth(2);

    hA[ip] = new TH1D(Form("hxjRegionA_pt%d", ip), "", nb, edges);
    hC[ip] = new TH1D(Form("hxjRegionC_pt%d", ip), "", nb, edges);
  }

  std::string fname = "/home/samson72/sphnx/gammajet/trees/gammajet_Data.root";
  TFile * f = TFile::Open(fname.c_str(), "read");
  if (!f || f->IsZombie()) { std::cout << "Could not open " << fname << std::endl; return; }
  TTree * t = (TTree*)f->Get("towerntup");
  Long64_t nentries = t->GetEntries();
  if (maxEntries >= 0 && maxEntries < nentries) nentries = maxEntries;
  std::cout << "Processing " << fname << " (" << nentries << " entries, Data)" << std::endl;

  t->SetBranchStatus("*", 0);
  for (const char * bn : {"vz", "cluster_pt", "cluster_eta", "cluster_phi",
       "cluster_showershape", "cluster_bdt_scores", "jet_eta", "jet_phi", "jet_pt_calib"})
    t->SetBranchStatus(bn, 1);

  Float_t vz, cluster_pt, cluster_eta, cluster_phi;
  Float_t cluster_showershape[12], cluster_bdt_scores[11];
  Float_t jet_pt_calib[ana::nJetR] = {0}, jet_eta[ana::nJetR], jet_phi[ana::nJetR];

  t->SetBranchAddress("vz", &vz);
  t->SetBranchAddress("cluster_pt", &cluster_pt);
  t->SetBranchAddress("cluster_eta", &cluster_eta);
  t->SetBranchAddress("cluster_phi", &cluster_phi);
  t->SetBranchAddress("cluster_showershape", cluster_showershape);
  t->SetBranchAddress("cluster_bdt_scores", cluster_bdt_scores);
  t->SetBranchAddress("jet_eta", jet_eta);
  t->SetBranchAddress("jet_phi", jet_phi);
  t->SetBranchAddress("jet_pt_calib", jet_pt_calib);

  for (Long64_t e = 0; e < nentries; e++) {
    t->GetEntry(e);
    if (e % 2000000 == 0)
      std::cout << "  entry " << e << "/" << nentries
                << " (" << (float)e / nentries * 100. << "%)" << std::endl;

    if (fabs(vz) > ana::vzcut) continue;

    // "Everything new": cluster_pt/eta/phi (with-saturation), cluster_bdt_scores[9] (new
    // bdt), cluster_showershape[11] (new iso), jet_pt_calib[ir] (new jet pt) - same as
    // compare_old_new.C's case 5/draw_xj_data_mc.C, Data side.
    float phoPt = cluster_pt, phoEta = cluster_eta, phoPhi = cluster_phi;
    float jetEta = jet_eta[ir], jetPhi = jet_phi[ir];
    float jetPt = jet_pt_calib[ir];
    float iso = cluster_showershape[11];
    float bdt = cluster_bdt_scores[9];

    if (!std::isfinite(phoPt) || phoPt <= 0 || !std::isfinite(jetPt)) continue;
    if (!(phoPt >= ana::ptBins[0] && phoPt < ana::ptBins[ana::nPtBins] && jetPt > ana::jet_calib_pt_cut[ir]))
      continue;
    float dphi = deltaPhi(phoPhi, jetPhi);
    if (!checkPair(jetPt, jetEta, phoPt, phoEta, dphi)) continue;

    int iabcd = ana::findabcdBin(iso, bdt, 0);
    if (iabcd != 0 && iabcd != 2) continue;

    float xj = jetPt / phoPt;
    if (!std::isfinite(xj)) continue;

    int ptbin = ana::findPtBin(phoPt);
    int iused = ptbin - ana::firstUsedPtBin;
    if (iused < 0 || iused >= ana::nPtBinsUsed) continue;

    if (iabcd == 0) {
      hRaw[iused]->Fill(xj);
      accRaw[iused].fill(xj, 1);
      hA[iused]->Fill(xj);
      accA[iused].fill(xj, 1);
    } else { // iabcd == 2
      hC[iused]->Fill(xj);
      accC[iused].fill(xj, 1);
    }
  }
  f->Close();

  // Combine Region A/Region C into the purity-corrected result, per used pT bin - see
  // file header comment. Mirrors compare_old_new.C's case 6 combination exactly (same
  // unfold_utility::purityCorrectCoeffs/purityCorrect calls, same ana::getPurity/
  // getPurityC "nominal"/ir=2 inputs).
  for (int ip = 0; ip < ana::nPtBinsUsed; ip++) {
    double ptlo = ana::ptBinsUsed[ip], pthi = ana::ptBinsUsed[ip + 1];
    float pA        = ana::getPurity(ptlo, pthi, "nominal", ir);
    float pAErrLow  = ana::getPurityErrorLow(ptlo, pthi, "nominal", ir);
    float pAErrHigh = ana::getPurityErrorHigh(ptlo, pthi, "nominal", ir);
    float pC        = ana::getPurityC(ptlo, pthi, "nominal", ir);
    float pCErrLow  = ana::getPurityCErrorLow(ptlo, pthi, "nominal", ir);
    float pCErrHigh = ana::getPurityCErrorHigh(ptlo, pthi, "nominal", ir);

    double NA = accA[ip].sumw, NC = accC[ip].sumw;
    float coeffA, coeffC;
    unfold_utility::purityCorrectCoeffs(pA, pC, NA, NC, coeffA, coeffC);
    double Ncorr = coeffA * NA - coeffC * NC;
    if (Ncorr > 0) {
      accCorr[ip].sumw   = Ncorr;
      accCorr[ip].sumw2  = Ncorr; // effective-N == Ncorr, matching insitu_utility::computeCorrectedMeans' err = sqrt(var/Ncorr)
      accCorr[ip].sumwx  = coeffA * accA[ip].sumwx  - coeffC * accC[ip].sumwx;
      accCorr[ip].sumwx2 = coeffA * accA[ip].sumwx2 - coeffC * accC[ip].sumwx2;
      accCorr[ip].n      = accA[ip].n + accC[ip].n;
    }

    TH1D * hcorr = unfold_utility::purityCorrect(hA[ip], hC[ip], pA, pAErrLow, pAErrHigh,
        pC, pCErrLow, pCErrHigh, Form("hPurityCorrTmp_pt%d", ip));
    if (hcorr) {
      hCorr[ip]->Add(hcorr);
      delete hcorr;
    }
    delete hA[ip];
    delete hC[ip];
  }

  std::cout << "\n" << Form("%-14s %10s %10s %10s", "pT bin", "<xJ> raw", "<xJ> corr", "corr/raw") << std::endl;
  for (int ip = 0; ip < ana::nPtBinsUsed; ip++) {
    double mRaw = accRaw[ip].mean(), mCorr = accCorr[ip].mean();
    std::cout << Form("%-14s %10.4f %10.4f %10.4f",
        Form("%.0f-%.0f GeV", ana::ptBinsUsed[ip], ana::ptBinsUsed[ip + 1]), mRaw, mCorr,
        (mRaw > 0 && mCorr != 0) ? mCorr / mRaw : 0) << std::endl;
  }

  std::string pdfPath = "/home/samson72/sphnx/gammajet_unfold/temporary_study/draw_xj_data_purity.pdf";
  TCanvas * c = new TCanvas("c_xj_data_purity", "", 700, 600);
  c->SetLeftMargin(.13);
  c->SetBottomMargin(.13);
  c->SaveAs((pdfPath + "[").c_str());

  for (int ip = 0; ip < ana::nPtBinsUsed; ip++) {
    c->Clear();
    // ana::unfoldXjBins is variable-width - see file header comment; Integral() omits
    // "width" (total count, the correct normalization denominator) but Scale needs it
    // so displayed content is a density, not a bare bin fraction.
    double integRaw = hRaw[ip]->Integral(0, hRaw[ip]->GetNbinsX() + 1);
    double integCorr = hCorr[ip]->Integral(0, hCorr[ip]->GetNbinsX() + 1);
    if (integRaw > 0) hRaw[ip]->Scale(1.0 / integRaw, "width");
    if (integCorr != 0) hCorr[ip]->Scale(1.0 / integCorr, "width");

    double ymax = std::max(hRaw[ip]->GetMaximum(), hCorr[ip]->GetMaximum());
    hRaw[ip]->SetTitle(";x_{J#gamma};(1/N) dN/dx_{J#gamma}");
    hRaw[ip]->GetYaxis()->SetRangeUser(0, ymax * 1.5);
    hRaw[ip]->Draw("hist");
    hCorr[ip]->Draw("hist same");

    TLegend * leg = new TLegend(0.47, 0.60, 0.88, 0.83);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    // <x_J> from the exact unbinned Accum, NOT TH1::GetMean() on the density-scaled
    // histogram above - see Accum's comment.
    leg->AddEntry(hRaw[ip], Form("Data, raw  #LTx_{J}#GT=%.3f#pm%.3f", accRaw[ip].mean(), accRaw[ip].meanErr()), "l");
    leg->AddEntry(hCorr[ip], Form("Data, purity-corr.  #LTx_{J}#GT=%.3f#pm%.3f", accCorr[ip].mean(), accCorr[ip].meanErr()), "l");
    leg->Draw();

    d.drawAll({"p+p Run24 Data: raw vs. purity-corrected"}, {Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",
              ana::ptBinsUsed[ip], ana::ptBinsUsed[ip + 1]), jetFeature}, .17, .87, 14, gPad->GetWh() * 0.8);
    c->SaveAs(pdfPath.c_str());
  }

  c->SaveAs((pdfPath + "]").c_str());
  std::cout << "Wrote " << pdfPath << std::endl;
}
