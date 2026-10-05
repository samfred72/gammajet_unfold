#include "../src/ana.h"
#include "../src/drawer.h"
#include "../src/unfold_utility.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
// Load explicitly: the sibling gammajet project's libgammajet.so has same-named classes.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Iteration scan on Data: the purity-corrected xJ unfolded through the combined Photon5/10/20
// response, compared between consecutive iterations (and with the MC truth for reference).
// Not a strict closure test (Data has no truth). The half-MC closure version preferred 1
// iteration, since two halves of one sample share the same expectation.

const int ir = 2; // nominal jet radius index (R=0.4)
const int nPtBinsUsed = ana::nPtBinsUsed; // reported bins start at ana::firstUsedPtBin
const vector<int> iterationsToTest = {1,2,3,4,5,6,7,8,9,10};

// The last 3 xJ bins per pT bin are low-count noise: excluded from chi2 only (still drawn).
const int nXjBinsForChi2 = ana::nUnfoldXjBins - 3;

// Red (i=0) to purple (i=n-1): HSV hue 0-270.
int rainbowColor(int i, int n) {
  float hue = (n > 1) ? 270.0 * i / (n - 1) : 0;
  float r, g, b;
  TColor::HSV2RGB(hue, 1.0, 1.0, r, g, b);
  return TColor::GetColor(r, g, b);
}

void draw_iteration_halfclosure(string systag = "nominal") {
  gStyle->SetOptStat(0);

  drawer d("pythia", systag);
  string pdfPath = Form("%s/pdfs/iteration_halfclosure_%s.pdf", ana::dir(), systag.c_str());

  // Response: Photon5/10/20 combined (Data's stored response is a diagonal placeholder).
  TH1D * respRecoTemplate  = d.get(Form("hrecoxj%i",ir), 1);
  TH1D * respTruthTemplate = d.get(Form("htruthxj%i",ir), 1);
  TH2D * respMatrix2D      = d.get2d(Form("hxjresponse%i",ir), 1);
  RooUnfoldResponse * response = new RooUnfoldResponse(respRecoTemplate, respTruthTemplate, respMatrix2D);

  // Measured: Data's purity-corrected xJ; truth reference: the MC truth.
  TH1D * flatA = d.get(Form("hrecoxj%i_0",ir), 0);
  TH1D * flatC = d.get(Form("hrecoxj%i_2",ir), 0);
  TH1D * flatMeasured = unfold_utility::buildFullyCorrected(flatA, flatC, "data", systag);
  TH1D * flatTruth = respTruthTemplate;

  // Per iteration: mean fractional uncertainty of the unfolded result (noise amplification).
  vector<int> iters;
  vector<double> errorScore;
  vector<vector<TH1D*>> unfoldedByIter; // [iterIdx][ipt]

  for (unsigned k = 0; k < iterationsToTest.size(); k++) {
    int iter = iterationsToTest[k];
    TH1D * hUnfoldFull = unfold_utility::unfoldOnce(response, flatMeasured, iter, Form("hUnfoldData_iter%d", iter));

    double fracErrSum = 0;
    int nBinsCounted = 0;
    vector<TH1D*> perPt(ana::nPtBins);
    for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
      TH1D * hUnfold = unfold_utility::unflattenXj(hUnfoldFull, ipt, Form("hxjunfold_iter%d_pt%d", iter, ipt));
      for (int b = 1; b <= nXjBinsForChi2; b++) {
        double u  = hUnfold->GetBinContent(b);
        double ue = hUnfold->GetBinError(b);
        if (u > 0) { fracErrSum += ue/u; nBinsCounted++; }
      }
      perPt[ipt] = hUnfold;
    }
    iters.push_back(iter);
    errorScore.push_back(nBinsCounted > 0 ? fracErrSum/nBinsCounted : 0);
    unfoldedByIter.push_back(perPt);
  }

  // Convergence metric: chi2/ndf of the relative change between consecutive iterations, summed
  // over used pT bins (defined from the second tested iteration on).
  vector<int> pairIters;
  vector<double> pairBiasScore;
  for (unsigned k = 1; k < iters.size(); k++) {
    double chi2 = 0;
    int ndf = 0;
    for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
      TH1D * uPrev = unfoldedByIter[k-1][ipt];
      TH1D * uNext = unfoldedByIter[k][ipt];
      for (int b = 1; b <= nXjBinsForChi2; b++) {
        double up = uPrev->GetBinContent(b);
        double un = uNext->GetBinContent(b);
        if (up <= 0) continue;
        double relDiff = (un - up)/up;
        chi2 += relDiff*relDiff;
        ndf++;
      }
    }
    pairIters.push_back(iters[k]);
    pairBiasScore.push_back(ndf > 0 ? chi2/ndf : 0);
  }

  // Combined score: both metrics min-max normalized and summed; minimum is best.
  double biasMin = *min_element(pairBiasScore.begin(), pairBiasScore.end());
  double biasMax = *max_element(pairBiasScore.begin(), pairBiasScore.end());
  double errMin  = *min_element(errorScore.begin()+1, errorScore.end());
  double errMax  = *max_element(errorScore.begin()+1, errorScore.end());
  int bestPairIdx = 0;
  double bestCombined = 1e18;
  for (unsigned k = 0; k < pairIters.size(); k++) {
    double normBias = (biasMax > biasMin) ? (pairBiasScore[k]-biasMin)/(biasMax-biasMin) : 0;
    double normErr  = (errMax  > errMin ) ? (errorScore[k+1]-errMin)/(errMax-errMin)     : 0;
    double combined = normBias + normErr;
    if (combined < bestCombined) { bestCombined = combined; bestPairIdx = k; }
  }
  int bestIter = pairIters[bestPairIdx];
  int bestIdx = bestPairIdx + 1; // offset by the dropped first iteration
  cout << "Data iteration scan: best iteration = " << bestIter
       << " (relative-change chi2=" << pairBiasScore[bestPairIdx] << ", mean frac. error=" << errorScore[bestIdx] << ")" << endl;

  TCanvas * c = new TCanvas("c","",700,900);
  c->SaveAs(Form("%s[", pdfPath.c_str()));

  // Pages 1..N: the chosen iteration's unfolded Data vs MC truth per used pT bin.
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hUnfold = unfoldedByIter[bestIdx][ipt];
    TH1D * hTruth  = unfold_utility::unflattenXj(flatTruth, ipt, Form("hxjtruth_final_pt%d", ipt));
    // Shape-normalize: Data and weighted MC truth have very different scales.
    TH1D * hUnfoldDisp = unfold_utility::densityForDisplay(hUnfold, Form("hxjunfold_final_pt%d_disp", ipt));
    hUnfoldDisp->Scale(1./hUnfoldDisp->Integral());
    TH1D * hTruthDisp  = unfold_utility::densityForDisplay(hTruth,  Form("hxjtruth_final_pt%d_disp", ipt));
    hTruthDisp->Scale(1./hTruthDisp->Integral());

    c->Clear();
    c->cd();
    TPad * pp1 = new TPad(Form("pp1_%d",ipt),"",0,.35,1,1);
    TPad * pp2 = new TPad(Form("pp2_%d",ipt),"",0,0,1,.35);
    pp1->Draw();
    pp2->Draw();

    pp1->cd();
    pp1->SetBottomMargin(0.02);
    pp1->SetLeftMargin(.15);
    gPad->SetTicks(1,1);
    hUnfoldDisp->SetLineColor(kRed);
    hUnfoldDisp->SetMarkerColor(kRed);
    hUnfoldDisp->SetMarkerStyle(21);
    hUnfoldDisp->SetLineWidth(2);
    hUnfoldDisp->GetXaxis()->SetLabelSize(0);
    hUnfoldDisp->GetXaxis()->SetTitle("");
    hUnfoldDisp->GetYaxis()->SetRangeUser(0, std::max(hUnfoldDisp->GetMaximum(), hTruthDisp->GetMaximum())*1.4);
    hUnfoldDisp->Draw("p e");
    hTruthDisp->SetLineColor(kBlack);
    hTruthDisp->SetMarkerColor(kBlack);
    hTruthDisp->SetMarkerStyle(20);
    hTruthDisp->SetLineWidth(2);
    hTruthDisp->Draw("p e same");
    TLegend * l = new TLegend(.55,.65,.85,.80);
    l->SetLineWidth(0);
    l->SetTextSize(0.035);
    l->AddEntry(hUnfoldDisp, Form("Data unfolded (%d iter.)", bestIter));
    l->AddEntry(hTruthDisp,  "Truth (#gamma+jet MC)");
    l->Draw();
    d.drawAll({"p+p Run24 Data"},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
        Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV", ana::JetRs[ir], ana::jet_calib_pt_cut[ir])}, .5, .85, 16, gPad->GetWh()*0.8);

    pp2->cd();
    pp2->SetTopMargin(0.02);
    pp2->SetBottomMargin(0.3);
    pp2->SetLeftMargin(.15);
    gPad->SetTicks(1,1);
    TH1D * hratio = (TH1D*)hUnfoldDisp->Clone(Form("hxjratio_final_pt%d", ipt));
    hratio->Divide(hTruthDisp);
    hratio->SetLineColor(kBlack);
    hratio->SetMarkerColor(kBlack);
    hratio->SetMarkerStyle(20);
    hratio->GetYaxis()->SetRangeUser(0.5,1.5);
    hratio->GetYaxis()->SetTitle("Unfolded / Truth");
    hratio->GetYaxis()->SetTitleSize(0.09);
    hratio->GetYaxis()->SetTitleOffset(0.7);
    hratio->GetYaxis()->SetLabelSize(0.08);
    hratio->GetXaxis()->SetTitle("x_{J#gamma}");
    hratio->GetXaxis()->SetTitleSize(0.09);
    hratio->GetXaxis()->SetLabelSize(0.08);
    hratio->Draw("p e");
    TLine * line = new TLine(ana::unfoldXjBins[0],1,ana::unfoldXjBins[ana::nUnfoldXjBins],1);
    line->SetLineStyle(9);
    line->Draw("same");
    c->SaveAs(pdfPath.c_str());
  }

  // Next pages: all iterations overlaid, red (fewest) to purple (most), one page per pT bin.
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hTruth0     = unfold_utility::unflattenXj(flatTruth, ipt, Form("hxjtruth_perIter_pt%d", ipt));
    TH1D * hTruthDisp0 = unfold_utility::densityForDisplay(hTruth0, Form("hxjtruth_perIter_pt%d_disp", ipt));
    hTruthDisp0->Scale(1./hTruthDisp0->Integral());

    vector<TH1D*> hUnfoldDispAll(iters.size());
    double ymax0 = hTruthDisp0->GetMaximum();
    for (unsigned k = 0; k < iters.size(); k++) {
      TH1D * hUnfold = unfoldedByIter[k][ipt];
      hUnfoldDispAll[k] = unfold_utility::densityForDisplay(hUnfold, Form("hxjunfold_perIter_pt%d_iter%d_disp", ipt, iters[k]));
      hUnfoldDispAll[k]->Scale(1./hUnfoldDispAll[k]->Integral());
      ymax0 = std::max(ymax0, hUnfoldDispAll[k]->GetMaximum());
    }

    c->Clear();
    c->cd();
    TPad * ppi1 = new TPad(Form("ppi1_pt%d",ipt),"",0,.35,1,1);
    TPad * ppi2 = new TPad(Form("ppi2_pt%d",ipt),"",0,0,1,.35);
    ppi1->Draw();
    ppi2->Draw();

    ppi1->cd();
    ppi1->SetBottomMargin(0.02);
    ppi1->SetLeftMargin(.15);
    gPad->SetTicks(1,1);
    hTruthDisp0->SetLineColor(kBlack);
    hTruthDisp0->SetMarkerColor(kBlack);
    hTruthDisp0->SetMarkerStyle(20);
    hTruthDisp0->SetLineWidth(2);
    hTruthDisp0->GetXaxis()->SetLabelSize(0);
    hTruthDisp0->GetXaxis()->SetTitle("");
    hTruthDisp0->GetYaxis()->SetRangeUser(0, ymax0*1.4);
    hTruthDisp0->Draw("p e");
    TLegend * li = new TLegend(.63,.33,.88,.80);
    li->SetLineWidth(0);
    li->SetTextSize(0.026);
    li->AddEntry(hTruthDisp0, "Truth (#gamma+jet MC)");
    for (unsigned k = 0; k < iters.size(); k++) {
      int col = rainbowColor(k, iters.size());
      hUnfoldDispAll[k]->SetLineColor(col);
      hUnfoldDispAll[k]->SetMarkerColor(col);
      hUnfoldDispAll[k]->SetMarkerStyle(20);
      hUnfoldDispAll[k]->SetLineWidth(2);
      hUnfoldDispAll[k]->Draw("p e same");
      li->AddEntry(hUnfoldDispAll[k], Form("%d iter.%s", iters[k], iters[k] == bestIter ? " (best)" : ""));
    }
    li->Draw();
    d.drawAll({"p+p Run24 Data"},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
        Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV", ana::JetRs[ir], ana::jet_calib_pt_cut[ir])}, .5, .88, 14, gPad->GetWh()*0.8);

    ppi2->cd();
    ppi2->SetTopMargin(0.02);
    ppi2->SetBottomMargin(0.3);
    ppi2->SetLeftMargin(.15);
    gPad->SetTicks(1,1);
    TH1D * hratioFrame = (TH1D*)hTruthDisp0->Clone(Form("hratio_perIter_frame_pt%d", ipt));
    hratioFrame->Reset("ICES");
    hratioFrame->SetLineColor(kWhite);
    hratioFrame->SetMarkerColor(kWhite);
    hratioFrame->GetYaxis()->SetRangeUser(0.5,1.5);
    hratioFrame->GetYaxis()->SetTitle("Unfolded / Truth");
    hratioFrame->GetYaxis()->SetTitleSize(0.09);
    hratioFrame->GetYaxis()->SetTitleOffset(0.7);
    hratioFrame->GetYaxis()->SetLabelSize(0.08);
    hratioFrame->GetXaxis()->SetTitle("x_{J#gamma}");
    hratioFrame->GetXaxis()->SetTitleSize(0.09);
    hratioFrame->GetXaxis()->SetLabelSize(0.08);
    hratioFrame->Draw("p");
    for (unsigned k = 0; k < iters.size(); k++) {
      TH1D * hratioi = (TH1D*)hUnfoldDispAll[k]->Clone(Form("hxjratio_perIter_pt%d_iter%d", ipt, iters[k]));
      hratioi->Divide(hTruthDisp0);
      hratioi->Draw("p e same");
    }
    TLine * linei = new TLine(ana::unfoldXjBins[0],1,ana::unfoldXjBins[ana::nUnfoldXjBins],1);
    linei->SetLineStyle(9);
    linei->Draw("same");
    c->SaveAs(pdfPath.c_str());
  }

  // Last plot: relative-change chi2/ndf (top) and uncertainty (bottom) vs iteration count.
  c->Clear();
  c->cd();
  TPad * p1 = new TPad("p1","",0,.5,1,1);
  TPad * p2 = new TPad("p2","",0,0,1,.5);
  p1->Draw();
  p2->Draw();

  TGraph * gBias = new TGraph(pairIters.size());
  TGraph * gErr  = new TGraph(pairIters.size());
  for (unsigned k = 0; k < pairIters.size(); k++) {
    gBias->SetPoint(k, pairIters[k], pairBiasScore[k]);
    gErr->SetPoint(k, pairIters[k], errorScore[k+1]);
  }

  p1->cd();
  p1->SetBottomMargin(0.02);
  p1->SetLeftMargin(.15);
  gPad->SetTicks(1,1);
  gBias->SetMarkerStyle(20);
  gBias->SetMarkerColor(kBlack);
  gBias->SetLineColor(kBlack);
  gBias->SetLineWidth(2);
  gBias->GetXaxis()->SetLabelSize(0);
  gBias->GetYaxis()->SetTitle("#chi^{2}/NDF of (iter n - iter n-1)/(iter n-1)");
  gBias->SetMinimum(0);
  gBias->SetMaximum(biasMax*1.3);
  gBias->Draw("APL");
  TLine * lbias = new TLine(bestIter, 0, bestIter, biasMax*1.3);
  lbias->SetLineStyle(9);
  lbias->SetLineColor(kRed);
  lbias->Draw("same");
  d.drawAll({"p+p Run24 Data"},{Form("Jet R=%.1f",ana::JetRs[ir]),
      Form("Best (min combined score): %d iterations",bestIter)}, .5, .85, 16, 700);

  p2->cd();
  p2->SetTopMargin(0.02);
  p2->SetBottomMargin(0.3);
  p2->SetLeftMargin(.15);
  gPad->SetTicks(1,1);
  gErr->SetMarkerStyle(20);
  gErr->SetMarkerColor(kBlue+1);
  gErr->SetLineColor(kBlue+1);
  gErr->SetLineWidth(2);
  gErr->GetXaxis()->SetTitle("Bayesian unfolding iterations");
  gErr->GetYaxis()->SetTitle("Mean fractional uncertainty (variance)");
  gErr->GetYaxis()->SetTitleSize(0.05);
  gErr->GetXaxis()->SetTitleSize(0.05);
  gErr->SetMinimum(0);
  gErr->SetMaximum(errMax*1.3);
  gErr->Draw("APL");
  TLine * lerr = new TLine(bestIter, 0, bestIter, errMax*1.3);
  lerr->SetLineStyle(9);
  lerr->SetLineColor(kRed);
  lerr->Draw("same");
  c->SaveAs(pdfPath.c_str());

  c->SaveAs(Form("%s]", pdfPath.c_str()));
}
