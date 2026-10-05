#include "../src/ana.h"
#include "../src/drawer.h"
#include "../src/unfold_utility.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
// Load explicitly: the sibling gammajet project's libgammajet.so has same-named classes.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Correlation matrix of the unfolded Data at every scanned niter, per pT bin's xJ block. The
// chi2 metrics elsewhere assume a diagonal covariance; Bayesian unfolding correlates neighbors
// (ATLAS PLB 774 (2017) 379). RooUnfoldBayes is built directly to read Eunfold(kCovariance).
// Self-check: sqrt(cov(i,i)) must equal Hreco()'s bin error.

const int ir = 2; // nominal jet radius index (R=0.4)
const int nPtBinsUsed = ana::nPtBinsUsed; // reported bins start at ana::firstUsedPtBin
const int niterate = 2; // nominal iteration count
const vector<int> iterationsToScan = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}; // same list as every other niter scan in this directory

// Blue-white-red diverging palette for correlations.
void setDivergingPalette() {
  const int nStops = 3;
  double stops[nStops] = {0.0, 0.5, 1.0};
  double red[nStops]   = {0.13, 1.00, 0.70};
  double green[nStops] = {0.30, 1.00, 0.05};
  double blue[nStops]  = {0.75, 1.00, 0.10};
  TColor::CreateGradientColorTable(nStops, stops, red, green, blue, 255);
  gStyle->SetNumberContours(255);
}

// flatbin(ipt,ixj), 0-indexed: TVectorD element i is histogram bin i+1 (no overflow).
int matrixIndex(int ipt, int ixj) {
  return ipt*(ana::nUnfoldXjBins+2) + ixj + 1;
}

void draw_covariance(string systag = "nominal") {
  gStyle->SetOptStat(0);
  setDivergingPalette();

  drawer d("pythia", systag);
  string pdfPath  = Form("%s/pdfs/draw_covariance_%s.pdf", ana::dir(), systag.c_str());
  string rootPath = Form("%s/hists/covariance_%s.root", ana::dir(), systag.c_str());

  // Response and purity-corrected Data, as in draw_refolding.C.
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

    // Built directly so Eunfold() is available.
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

        // Self-check against Hreco()'s bin error; a mismatch means the index bookkeeping is wrong.
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

  // Bonus page: mean |off-diagonal correlation| vs iteration count.
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
