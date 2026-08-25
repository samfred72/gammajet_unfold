#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/unfold_utility.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "TDecompSVD.h"
// The original gammajet project builds its OWN, differently-laid-out drawer/ana classes
// into /home/samson72/root/lib/libgammajet.so, sitting on the same library search path
// as this project's libgammajet_unfold.so. Without forcing which one loads first, ROOT's
// implicit symbol autoload can bind drawer/ana calls to the wrong (mismatched-layout)
// library and segfault - explicit load removes the ambiguity.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Follow-up to draw_covariance.C: does the correlation found there (mean |off-diagonal|
// ~0.2-0.25, decaying from niter=1 then levelling off) actually change the chi2/NDF
// conclusions this whole directory's niter selection rests on, or is the diagonal
// approximation (Sum (residual/sigma_i)^2, used everywhere else in this directory) good
// enough in practice? This computes BOTH versions of the SAME comparison - unfolded Data
// vs the Pythia8 truth prior, the exact reference draw_purity_corrected.C's own
// best-iteration scan is built on - at every niter in the standard scan, so the two can be
// read off side by side.
//
// Diagonal: Sum_i (residual_i)^2 / sigma_i^2, sigma_i^2 = Eunfold(kCovariance)'s own
// diagonal element for bin i (equivalent to what every computeChi2NDF-style function
// elsewhere in this directory already computes, just sourced from the covariance matrix's
// diagonal directly rather than TH1::GetBinError()).
// Full covariance: res^T C^-1 res, using the FULL (not diagonal-truncated) covariance
// submatrix restricted to the same used bins, inverted directly (RooUnfoldT::Chi2()
// itself does exactly this with TMatrixD::Invert() - not reused here only because it
// operates over the response's ENTIRE flattened space, including the low-pT/high-pT
// buffer bins and the low-count high-xJ tail this directory already excludes everywhere
// else via nXjBinsForChi2 - so the restricted submatrix is built and inverted directly).
//
// Only the UNFOLDED result's own covariance is used for both versions - the Pythia truth
// reference's own (much smaller, weighted-MC-sample) statistical error is not folded in,
// matching the simplification already used by every chi2 metric in this directory besides
// draw_purity_corrected.C's original computeChi2NDF (which combines both sides) - keeping
// this file to ONE side's uncertainty isolates the effect of the correlation itself,
// rather than mixing in a different error-combination convention too.
//
// Shape-normalized per pT bin before comparing, exactly like draw_purity_corrected.C's own
// unfolded-vs-truth chi2 (each side scaled to unit area over its own pT bin's full
// nUnfoldXjBins range, THEN restricted to the first nXjBinsForChi2 of those for the actual
// sum) - Data (raw counts) and the cross-section-weighted Pythia truth are on wildly
// different absolute scales, so subtracting them directly (as an earlier version of this
// file did) makes every residual enormous regardless of how correct the covariance
// treatment is. The unfolded side's normalization constant (1/Integral()) is treated as a
// fixed scalar here - so the covariance submatrix for a pT bin is scaled by that constant
// squared - the same level of rigor (not propagating the normalization's own uncertainty)
// every other shape-normalized comparison in this directory already uses.

const int ir = 2; // nominal jet radius index (R=0.4)
const int nPtBinsUsed = ana::nPtBinsUsed; // physics analysis only uses ana::ptBins[ana::firstUsedPtBin..]
const int niterate = 2; // matches draw_purity_corrected.C / draw_final_result.C's chosen nominal iteration count
const vector<int> iterationsToScan = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}; // same list as every other niter scan in this directory

// The last 3 xJ bins in each pT bin have very low counts, so chi2/NDF here would be
// dominated by their noise rather than genuine convergence behavior - excluded from both
// chi2 metrics below, same exclusion as every other macro in this directory.
const int nXjBinsForChi2 = ana::nUnfoldXjBins - 3;

// buildFullyCorrected now lives in unfold_utility (purity-corrects all ana::nPtBins
// slices via unfold_utility::purityCorrect and reflattens for RooUnfold) - see
// src/unfold_utility.h.

// flatbin(ipt,ixj), 0-indexed, matching TMatrixD/TVectorD element ordering - identical to
// draw_covariance.C's own matrixIndex() (and to unfold_utility::unflattenXj's "flatbin"),
// since RooUnfoldResponse defaults to _overflow=false everywhere in this project.
int matrixIndex(int ipt, int ixj) {
  return ipt*(ana::nUnfoldXjBins+2) + ixj + 1;
}

void draw_covariance_chi2(string systag = "nominal") {
  gStyle->SetOptStat(0);

  drawer d("pythia", systag);
  string pdfPath  = Form("/home/samson72/sphnx/gammajet_unfold/pdfs/draw_covariance_chi2_%s.pdf", systag.c_str());
  string rootPath = Form("/home/samson72/sphnx/gammajet_unfold/hists/covariance_chi2_%s.root", systag.c_str());

  TH1D * respRecoTemplate  = d.get(Form("hrecoxj%i",ir), 1);
  TH1D * respTruthTemplate = d.get(Form("htruthxj%i",ir), 1);
  TH2D * respMatrix2D      = d.get2d(Form("hxjresponse%i",ir), 1);
  RooUnfoldResponse * response = new RooUnfoldResponse(respRecoTemplate, respTruthTemplate, respMatrix2D);

  TH1D * flatA = d.get(Form("hrecoxj%i_0",ir), 0);
  TH1D * flatC = d.get(Form("hrecoxj%i_2",ir), 0);
  TH1D * flatMeasured = unfold_utility::buildFullyCorrected(flatA, flatC, "data", systag);

  TFile * fout = TFile::Open(rootPath.c_str(), "RECREATE");
  TCanvas * c = new TCanvas("c","",700,900);
  c->SaveAs(Form("%s[", pdfPath.c_str()));

  cout << "Diagonal vs full-covariance chi2/NDF (unfolded Data vs Pythia8 truth, shape-"
       << "normalized per pT bin), " << nXjBinsForChi2 << " used xJ bins x " << nPtBinsUsed
       << " pT bins:" << endl;
  TGraph * gDiag = new TGraph((int)iterationsToScan.size());
  TGraph * gFull = new TGraph((int)iterationsToScan.size());

  for (unsigned k = 0; k < iterationsToScan.size(); k++) {
    int iter = iterationsToScan[k];

    RooUnfoldBayes unfold(response, flatMeasured, iter, 0, 1);
    unfold.SetVerbose(-1);
    TH1D * hUnfolded = (TH1D*)((TH1D*)unfold.Hreco())->Clone(Form("hUnfoldedFull_iter%d", iter));
    TMatrixD cov = unfold.Eunfold(RooUnfold::kCovariance);

    // Summed across the nPtBinsUsed pT bins, each normalized (and inverted) separately -
    // matching draw_purity_corrected.C's own per-pT-bin "combined +=" convention, not one
    // cross-pT-bin block (shape-normalization is itself a per-pT-bin operation, so a
    // single combined block would need cross-pT-bin normalization factors that don't
    // mean anything physical here).
    double chi2DiagSum = 0, chi2FullSum = 0;
    int ndfDiagSum = 0, nKeptSum = 0, nSum = 0;

    for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
      TH1D * hUnfoldPt = unfold_utility::unflattenXj(hUnfolded, ipt, Form("hUnfoldPt_iter%d_pt%d", iter, ipt));
      TH1D * hTruthPt  = unfold_utility::unflattenXj(respTruthTemplate, ipt, Form("hTruthPt_iter%d_pt%d", iter, ipt));
      double cU = hUnfoldPt->Integral() > 0 ? 1.0/hUnfoldPt->Integral() : 0; // over the FULL nUnfoldXjBins range, matching densityForDisplay+Scale(1/Integral()) convention
      double cT = hTruthPt->Integral()  > 0 ? 1.0/hTruthPt->Integral()  : 0;

      int n = nXjBinsForChi2;
      TMatrixD subCov(n, n);
      TVectorD res(n);
      for (int a = 0; a < n; a++) {
        int ia = matrixIndex(ipt, a);
        res[a] = hUnfoldPt->GetBinContent(a+1)*cU - hTruthPt->GetBinContent(a+1)*cT;
        for (int b = 0; b < n; b++) {
          int ib = matrixIndex(ipt, b);
          // cU*cU: the normalization constant is treated as a fixed scalar (not itself a
          // random variable correlated with the bins it's built from) - the same level of
          // rigor every shape-normalized comparison elsewhere in this directory already
          // uses. Truth's own covariance isn't folded in - see file header.
          subCov(a,b) = cov(ia, ib) * cU * cU;
        }
      }

      // Diagonal chi2/NDF: same formula as every other computeChi2NDF-style function in
      // this directory, just reading sigma_i^2 straight off subCov's own diagonal.
      double chi2Diag = 0;
      int ndfDiag = 0;
      for (int a = 0; a < n; a++) {
        double sigma2 = subCov(a,a);
        if (sigma2 <= 0) continue;
        chi2Diag += res[a]*res[a]/sigma2;
        ndfDiag++;
      }

      // Full-covariance chi2: invert the restricted (shape-normalized) submatrix and form
      // the proper quadratic form res^T C^-1 res - mirrors what RooUnfoldT::Chi2() does
      // internally, but via a Moore-Penrose PSEUDO-inverse (SVD, relative cutoff on small
      // singular values) rather than TMatrixD::Invert() - D'Agostini unfolding ties the
      // total unfolded yield to the actual measured total at every iteration, which bakes
      // at least one exact (or near-exact) linear dependency into the covariance. A plain
      // inverse doesn't degrade gracefully against that: near-zero eigenvalues get
      // inverted into astronomically large ones, blowing up the quadratic form. This is a
      // known, documented property of RooUnfold's own covariance matrices -
      // RooUnfoldT::Chi2()'s own header comment says it "removes rows/cols with all their
      // elements equal to 0" before inverting for exactly this reason; the SVD cutoff
      // below is the general version of that same defense.
      TDecompSVD svd(subCov);
      const TVectorD & sigma = svd.GetSig(); // singular values, descending
      const TMatrixD & U = svd.GetU();
      const TMatrixD & V = svd.GetV();
      double sigmaMax = sigma[0];
      const double relCut = 1e-10; // standard relative pseudo-inverse cutoff
      int nKept = 0;
      TMatrixD sigmaPInv(n, n);
      for (int i = 0; i < n; i++) {
        if (sigma[i] > relCut*sigmaMax) { sigmaPInv(i,i) = 1.0/sigma[i]; nKept++; }
        // else: near-null direction, dropped (0) rather than blown up.
      }
      TMatrixD Ut(TMatrixD::kTransposed, U);
      TMatrixD subCovPInv = V * sigmaPInv * Ut;
      double chi2Full = 0;
      for (int a = 0; a < n; a++) {
        double rowSum = 0;
        for (int b = 0; b < n; b++) rowSum += subCovPInv(a,b)*res[b];
        chi2Full += res[a]*rowSum;
      }

      cout << "  pt" << ipt << " niter=" << iter << ": diagonal chi2=" << chi2Diag << "/" << ndfDiag
           << ", full-cov chi2=" << chi2Full << "/" << nKept << " (" << nKept << "/" << n << " directions kept)" << endl;

      chi2DiagSum += chi2Diag; ndfDiagSum += ndfDiag;
      chi2FullSum += chi2Full; nKeptSum += nKept; nSum += n;
      delete hUnfoldPt; delete hTruthPt;
    }

    // NDF is the EFFECTIVE number of independent directions actually used (nKeptSum for
    // the covariance version), not the raw bin count - directions zeroed out above carry
    // no constraint, so they shouldn't be counted as degrees of freedom either.
    double chi2NdfDiag = ndfDiagSum > 0 ? chi2DiagSum/ndfDiagSum : 0;
    double chi2NdfFull = nKeptSum   > 0 ? chi2FullSum/nKeptSum   : 0;
    gDiag->SetPoint(k, iter, chi2NdfDiag);
    gFull->SetPoint(k, iter, chi2NdfFull);
    cout << "  niter=" << iter << " TOTAL: diagonal chi2/NDF=" << chi2NdfDiag
         << ", full-covariance chi2/NDF=" << chi2NdfFull << " (" << nKeptSum << "/" << nSum << " directions kept)"
         << ", ratio (full/diag)=" << (chi2NdfDiag > 0 ? chi2NdfFull/chi2NdfDiag : 0) << endl;

    delete hUnfolded;
  }

  c->Clear();
  c->cd();
  gPad->SetTicks(1,1);
  gPad->SetLeftMargin(.15);
  double ymax = 0;
  for (int k = 0; k < gDiag->GetN(); k++) { double x,y; gDiag->GetPoint(k,x,y); ymax = std::max(ymax,y); }
  for (int k = 0; k < gFull->GetN(); k++) { double x,y; gFull->GetPoint(k,x,y); ymax = std::max(ymax,y); }
  gDiag->SetMarkerStyle(20);
  gDiag->SetMarkerColor(kBlack);
  gDiag->SetLineColor(kBlack);
  gDiag->SetLineWidth(2);
  gDiag->GetXaxis()->SetTitle("Bayesian unfolding iterations");
  gDiag->GetYaxis()->SetTitle("#chi^{2}/NDF (Unfolded Data vs Pythia8 truth)");
  gDiag->SetMinimum(0);
  gDiag->SetMaximum(ymax*1.3);
  gDiag->Draw("APL");
  gFull->SetMarkerStyle(21);
  gFull->SetMarkerColor(kRed);
  gFull->SetLineColor(kRed);
  gFull->SetLineWidth(2);
  gFull->Draw("PL same");
  TLine * lnom = new TLine(niterate, 0, niterate, ymax*1.3);
  lnom->SetLineStyle(9);
  lnom->SetLineColor(kGray+2);
  lnom->Draw("same");
  TLegend * l = new TLegend(.5,.72,.85,.88);
  l->SetLineWidth(0);
  l->SetTextSize(0.03);
  l->AddEntry(gDiag, "Diagonal (Sum (res/#sigma_{i})^{2})", "lp");
  l->AddEntry(gFull, "Full covariance (res^{T} C^{-1} res)", "lp");
  l->Draw();
  d.drawAll({"p+p Run24 Data"},{Form("Jet R=%.1f",ana::JetRs[ir]),
      Form("Nominal: %d iterations (dashed line)",niterate)}, .18, .85, 14, 700);
  c->SaveAs(pdfPath.c_str());

  fout->cd();
  gDiag->Write("gChi2Diag");
  gFull->Write("gChi2FullCov");

  c->SaveAs(Form("%s]", pdfPath.c_str()));
  fout->Close();
  cout << "Done. Wrote " << pdfPath << " and " << rootPath << endl;
}
