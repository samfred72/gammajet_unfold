#include "../src/ana.h"
#include "../src/drawer.h"
#include "../src/unfold_utility.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
// The original gammajet project builds its OWN, differently-laid-out drawer/ana classes
// into /home/samson72/root/lib/libgammajet.so, sitting on the same library search path
// as this project's libgammajet_unfold.so. Without forcing which one loads first, ROOT's
// implicit symbol autoload can bind drawer/ana calls to the wrong (mismatched-layout)
// library and segfault - explicit load removes the ambiguity.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Covariance/correlation matrix of the unfolded Data result, at every scanned iteration
// count. Every chi2/NDF metric built elsewhere in this directory (draw_purity_corrected.C,
// draw_iteration_halfclosure.C, toy_iterations.C,
// draw_refolding.C, draw_nonclosure.C, draw_prior_sensitivity.C) sums (residual/sigma_i)^2
// bin by bin - i.e. assumes a DIAGONAL covariance matrix. Bayesian unfolding doesn't
// produce one: every iteration redistributes weight across bins through the same
// migration matrix, so a fluctuation in one bin systematically pulls its neighbors too.
// The ATLAS dijet-xJ paper's own toy-based covariance shows exactly this structure
// (Phys. Lett. B 774 (2017) 379, referenced earlier this session): "Nearby xJ bins show a
// strong positive correlation that diminishes for bins separated in xJ... Bins well
// separated in xJ show an anti-correlation attributable to the normalisation of
// (1/N)dN/dxJ." This macro checks whether that structure is present here, and how it
// changes with niter (more iterations couple bins more - the covariance's off-diagonal
// growth is the other side of the same regularization tradeoff niter already balances
// against statistical noise and prior bias throughout this directory).
//
// The covariance itself isn't something toy-computed here - RooUnfoldBayes already
// computes it analytically on every single unfoldOnce() call (it's what the
// "Calculating covariances..." console line, silenced everywhere else via SetVerbose(-1),
// refers to). unfold_utility::unfoldOnce() only ever keeps Hreco()'s diagonal (via
// TH1D::Clone(), which carries just GetBinError()), so a RooUnfoldBayes object is built
// directly here (not through that helper) specifically to also call
// Eunfold(RooUnfold::kCovariance), which returns the full TMatrixD.
//
// Displayed as a CORRELATION matrix (cov(i,j)/sqrt(cov(i,i)*cov(j,j)), range [-1,1]), not
// raw covariance - correlation is scale-free, so it's comparable across pT bins with very
// different absolute yields. Restricted to one pT bin's own xJ x xJ block per page (the
// full flattened matrix also carries cross-pT-bin correlations, from the response matrix's
// migration crossing pT-bin boundaries, but the ATLAS quote above - and the physical
// intuition for it - is specifically about correlations WITHIN a measured xJ spectrum).
//
// A self-check is printed for every bin: this macro's own extracted diagonal
// (sqrt(cov(i,i))) is compared to Hreco()'s own GetBinError() for the identical bin at
// the identical niter - they MUST match (both come from the same RooUnfoldBayes error
// propagation, just read out two different ways), so a mismatch means the flat-index
// bookkeeping below is wrong, not that anything about the unfolding itself is.

const int ir = 2; // nominal jet radius index (R=0.4)
const int nPtBinsUsed = ana::nPtBinsUsed; // physics analysis only uses ana::ptBins[ana::firstUsedPtBin..]
const int niterate = 2; // matches draw_purity_corrected.C / draw_final_result.C's chosen nominal iteration count
const vector<int> iterationsToScan = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}; // same list as every other niter scan in this directory

// densityForDisplay now lives in unfold_utility - see src/unfold_utility.h.
// buildFullyCorrected now lives in unfold_utility (purity-corrects all ana::nPtBins
// slices via unfold_utility::purityCorrect and reflattens for RooUnfold) - see
// src/unfold_utility.h.

// Blue (corr=-1) - white (corr=0) - red (corr=+1) diverging palette, the standard
// convention for correlation-matrix heatmaps (a sequential palette like ROOT's default
// kBird would make zero-correlation and strong-correlation regions hard to tell apart).
void setDivergingPalette() {
  const int nStops = 3;
  double stops[nStops] = {0.0, 0.5, 1.0};
  double red[nStops]   = {0.13, 1.00, 0.70};
  double green[nStops] = {0.30, 1.00, 0.05};
  double blue[nStops]  = {0.75, 1.00, 0.10};
  TColor::CreateGradientColorTable(nStops, stops, red, green, blue, 255);
  gStyle->SetNumberContours(255);
}

// flatbin(ipt,ixj), 0-indexed, matching TMatrixD/TVectorD element ordering: identical to
// unfold_utility::unflattenXj's own "flatbin" (ipt*(nUnfoldXjBins+2)+ixj+1), since
// RooUnfoldResponse defaults to _overflow=false (never enabled anywhere in this project),
// so TVectorD element i <-> ROOT histogram bin i+1, with no extra offset from
// under/overflow bins being folded in.
int matrixIndex(int ipt, int ixj) {
  return ipt*(ana::nUnfoldXjBins+2) + ixj + 1;
}

void draw_covariance(string systag = "nominal") {
  gStyle->SetOptStat(0);
  setDivergingPalette();

  drawer d("pythia", systag);
  string pdfPath  = Form("%s/pdfs/draw_covariance_%s.pdf", ana::dir(), systag.c_str());
  string rootPath = Form("%s/hists/covariance_%s.root", ana::dir(), systag.c_str());

  // Response matrix + purity-corrected Data - same construction as draw_refolding.C/
  // draw_prior_sensitivity.C.
  TH1D * respRecoTemplate  = d.get(Form("hrecoxj%i",ir), 1);
  TH1D * respTruthTemplate = d.get(Form("htruthxj%i",ir), 1);
  TH2D * respMatrix2D      = d.get2d(Form("hxjresponse%i",ir), 1);
  RooUnfoldResponse * response = new RooUnfoldResponse(respRecoTemplate, respTruthTemplate, respMatrix2D);

  TH1D * flatA = d.get(Form("hrecoxj%i_0",ir), 0);
  TH1D * flatC = d.get(Form("hrecoxj%i_2",ir), 0);
  TH1D * flatMeasured = unfold_utility::buildFullyCorrected(flatA, flatC, "data", systag);

  TFile * fout = TFile::Open(rootPath.c_str(), "RECREATE");
  TCanvas * c = new TCanvas("c","",800,800);
  c->SaveAs(Form("%s[", pdfPath.c_str()));

  cout << "Covariance/correlation scan: pT bin, niter, mean |off-diagonal correlation|" << endl;
  TGraph * gMeanCorr = new TGraph((int)iterationsToScan.size());

  for (unsigned k = 0; k < iterationsToScan.size(); k++) {
    int iter = iterationsToScan[k];

    // Built directly (not via unfold_utility::unfoldOnce) so the RooUnfoldBayes object -
    // and its Eunfold() - stays alive long enough to read the covariance out of, not just
    // the diagonal-only clone that helper returns.
    RooUnfoldBayes unfold(response, flatMeasured, iter, 0, 1);
    unfold.SetVerbose(-1);
    TH1D * hUnfolded = (TH1D*)((TH1D*)unfold.Hreco())->Clone(Form("hUnfoldedFull_iter%d", iter));
    TMatrixD cov = unfold.Eunfold(RooUnfold::kCovariance);

    double meanAbsCorrSum = 0;
    int nPairsSum = 0;

    for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
      TH1D * hUnfoldPt = unfold_utility::unflattenXj(hUnfolded, ipt, Form("hUnfoldPt_iter%d_pt%d", iter, ipt));

      TH2D * hCorr = new TH2D(Form("hCorr_iter%d_pt%d", iter, ipt), ";x_{J#gamma};x_{J#gamma}",
          ana::nUnfoldXjBins, ana::unfoldXjBins, ana::nUnfoldXjBins, ana::unfoldXjBins);

      double meanAbsCorr = 0;
      int nPairs = 0;
      for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) {
        int idxI = matrixIndex(ipt, ixj);
        double covII = cov(idxI, idxI);

        // Self-check: this bin's diagonal covariance element must reproduce Hreco()'s own
        // bin error for the SAME bin at the SAME niter (both come from the same
        // RooUnfoldBayes error propagation) - a mismatch means matrixIndex() is wrong,
        // not that anything about the unfolding is.
        double diagErr = sqrt(std::max(covII, 0.0));
        double histErr = hUnfoldPt->GetBinError(ixj+1);
        if (fabs(diagErr-histErr) > 1e-3*std::max(histErr,1.0))
          cout << "WARNING: covariance/Hreco() mismatch at pt" << ipt << " ixj" << ixj
               << " (cov sqrt=" << diagErr << ", hist=" << histErr << ") - check matrixIndex()." << endl;

        for (int jxj = 0; jxj < ana::nUnfoldXjBins; jxj++) {
          int idxJ = matrixIndex(ipt, jxj);
          double covJJ = cov(idxJ, idxJ);
          double corr = (covII > 0 && covJJ > 0) ? cov(idxI, idxJ)/sqrt(covII*covJJ) : 0;
          hCorr->SetBinContent(ixj+1, jxj+1, corr);
          if (ixj != jxj) { meanAbsCorr += fabs(corr); nPairs++; }
        }
      }
      meanAbsCorrSum += meanAbsCorr;
      nPairsSum += nPairs;
      double meanAbsCorrPt = nPairs > 0 ? meanAbsCorr/nPairs : 0;
      cout << "  pt" << ipt << " (" << ana::ptBins[ipt] << "-" << ana::ptBins[ipt+1] << " GeV), niter=" << iter
           << ": mean |off-diag corr| = " << meanAbsCorrPt << endl;

      c->Clear();
      c->cd();
      gPad->SetRightMargin(0.18);
      gPad->SetLeftMargin(0.15);
      gPad->SetTicks(1,1);
      hCorr->GetZaxis()->SetRangeUser(-1,1);
      hCorr->GetZaxis()->SetTitle("Correlation");
      hCorr->GetYaxis()->SetTitleOffset(1.3);
      hCorr->Draw("colz");
      d.drawAll({"p+p Run24 Data"},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
          Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV, %d iterations", ana::JetRs[ir], ana::jet_calib_pt_cut[ir], iter)}, .18, .95, 14, 800);
      c->SaveAs(pdfPath.c_str());

      fout->cd();
      hCorr->Write();
      delete hUnfoldPt;
    }
    double meanAbsCorrIter = nPairsSum > 0 ? meanAbsCorrSum/nPairsSum : 0;
    gMeanCorr->SetPoint(k, iter, meanAbsCorrIter);
    delete hUnfolded;
  }

  // Bonus page: mean |off-diagonal correlation| (averaged over used pT bins) vs iteration
  // count - the other side of the niter tradeoff already explored throughout this
  // directory (statistical noise and prior bias both generally SHRINK with more
  // iterations; this should generally GROW, since more iterations couple bins more).
  c->Clear();
  c->cd();
  gPad->SetTicks(1,1);
  gPad->SetLeftMargin(.15);
  gMeanCorr->SetMarkerStyle(20);
  gMeanCorr->SetMarkerColor(kBlack);
  gMeanCorr->SetLineColor(kBlack);
  gMeanCorr->SetLineWidth(2);
  gMeanCorr->GetXaxis()->SetTitle("Bayesian unfolding iterations");
  gMeanCorr->GetYaxis()->SetTitle("Mean |off-diagonal correlation|");
  double ymax = 0;
  for (int k = 0; k < gMeanCorr->GetN(); k++) { double x,y; gMeanCorr->GetPoint(k,x,y); ymax = std::max(ymax,y); }
  gMeanCorr->SetMinimum(0);
  gMeanCorr->SetMaximum(std::max(ymax*1.3, 0.1));
  gMeanCorr->Draw("APL");
  TLine * lnom = new TLine(niterate, 0, niterate, std::max(ymax*1.3, 0.1));
  lnom->SetLineStyle(9);
  lnom->SetLineColor(kRed);
  lnom->Draw("same");
  d.drawAll({"p+p Run24 Data"},{Form("Jet R=%.1f",ana::JetRs[ir]),
      Form("Nominal: %d iterations (dashed line)",niterate)}, .5, .85, 16, 700);
  c->SaveAs(pdfPath.c_str());
  fout->cd();
  gMeanCorr->Write("gMeanAbsCorrelation");

  c->SaveAs(Form("%s]", pdfPath.c_str()));
  fout->Close();
  cout << "Done. Wrote " << pdfPath << " and " << rootPath << endl;
}
