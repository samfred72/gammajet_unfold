#include "../src/ana.h"
#include "../src/drawer.h"
#include "../src/unfold_utility.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "TDecompSVD.h"
// Load explicitly: the sibling gammajet project's libgammajet.so has same-named classes.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Does the bin-to-bin correlation (draw_covariance.C) change the niter conclusions? Compares the
// diagonal chi2 with the full-covariance chi2 (res^T C^-1 res) of unfolded Data vs Pythia8 truth
// at every niter, per pT bin on the used bins only. Only the unfolded covariance enters.
// Both sides are shape-normalized over the full pT slice; the normalization is treated as a
// fixed scalar (covariance scaled by its square).

const int ir = 2; // nominal jet radius index (R=0.4)
const int nPtBinsUsed = ana::nPtBinsUsed; // reported bins start at ana::firstUsedPtBin
const int niterate = 2; // nominal iteration count
const vector<int> iterationsToScan = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}; // same list as every other niter scan in this directory

// The last 3 xJ bins per pT bin have too few counts: excluded from both chi2 metrics.
const int nXjBinsForChi2 = ana::nUnfoldXjBins - 3;

// flatbin(ipt,ixj), 0-indexed (TMatrixD order; RooUnfoldResponse has no overflow here).
int matrixIndex(int ipt, int ixj) {
  return ipt*(ana::nUnfoldXjBins+2) + ixj + 1;
}

void draw_covariance_chi2(string systag = "nominal") {
  gStyle->SetOptStat(0);

  drawer d("pythia", systag);
  string pdfPath  = Form("%s/pdfs/draw_covariance_chi2_%s.pdf", ana::dir(), systag.c_str());
  string rootPath = Form("%s/hists/covariance_chi2_%s.root", ana::dir(), systag.c_str());

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

    // Summed over pT bins, each normalized and inverted separately.
    double chi2DiagSum = 0, chi2FullSum = 0;
    int ndfDiagSum = 0, nKeptSum = 0, nSum = 0;

    for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
      TH1D * hUnfoldPt = unfold_utility::unflattenXj(hUnfolded, ipt, Form("hUnfoldPt_iter%d_pt%d", iter, ipt));
      TH1D * hTruthPt  = unfold_utility::unflattenXj(respTruthTemplate, ipt, Form("hTruthPt_iter%d_pt%d", iter, ipt));
      double cU = hUnfoldPt->Integral() > 0 ? 1.0/hUnfoldPt->Integral() : 0; // over the full nUnfoldXjBins range
      double cT = hTruthPt->Integral()  > 0 ? 1.0/hTruthPt->Integral()  : 0;

      int n = nXjBinsForChi2;
      TMatrixD subCov(n, n);
      TVectorD res(n);
      for (int a = 0; a < n; a++) {
        int ia = matrixIndex(ipt, a);
        res[a] = hUnfoldPt->GetBinContent(a+1)*cU - hTruthPt->GetBinContent(a+1)*cT;
        for (int b = 0; b < n; b++) {
          int ib = matrixIndex(ipt, b);
          subCov(a,b) = cov(ia, ib) * cU * cU;
        }
      }

      // Diagonal chi2/NDF from subCov's diagonal.
      double chi2Diag = 0;
      int ndfDiag = 0;
      for (int a = 0; a < n; a++) {
        double sigma2 = subCov(a,a);
        if (sigma2 <= 0) continue;
        chi2Diag += res[a]*res[a]/sigma2;
        ndfDiag++;
      }

      // Full-covariance chi2 via an SVD pseudo-inverse: D'Agostini fixes the total yield, so the
      // covariance is (near-)singular and a plain inverse blows up.
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
        // else: near-null direction, dropped.
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

    // NDF = number of kept directions.
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
