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

// Iteration scan on real Data, unfolded through the full (non-half) combined Photon MC
// response matrix.
//
// This started as a half-closure test (see git history / prior discussion): unfolder.cc's
// per-event use_half coin flip splits each MC sample into two independent halves, letting
// you unfold one half's reco spectrum through a response trained on the other half and
// compare to that half's own truth - a genuine closure test, immune to the circularity of
// testing a response against the same sample's own prior. That version answered "does more
// iterations help or just add noise" cleanly: 1 iteration won, monotonically, because two
// halves of the same MC sample share the same expectation value - there's no real bias for
// extra iterations to remove, only variance for them to add.
//
// This version instead unfolds Data's purity-corrected xJ spectrum
// (unfold_utility::purityCorrect/buildFullyCorrected below) through the full, properly
// cross-section-weighted Photon5/10/20 response (same response construction as
// drawing/draw_purity_corrected.C too), comparing to the fixed #gamma+jet MC truth as a
// reference. Since Data has no real truth level, this is no longer a strict closure
// test - Data can genuinely differ from the MC prior for real physics reasons, not just
// noise - but it shows how Data's unfolded result actually depends on iteration count,
// using the same purity-corrected input the real analysis result is built from.

const int ir = 2; // nominal jet radius index (R=0.4)
const int nPtBinsUsed = ana::nPtBinsUsed; // physics analysis only uses ana::ptBins[ana::firstUsedPtBin..]
const vector<int> iterationsToTest = {1,2,3,4,5,6,7,8,9,10};

// The last 3 xJ bins in each pT bin have very low counts, so any chi2/NDF computed
// against iteration count is dominated by their noise rather than genuine convergence
// behavior - excluded from errorScore/pairBiasScore below (still drawn everywhere else).
const int nXjBinsForChi2 = ana::nUnfoldXjBins - 3;

// densityForDisplay now lives in unfold_utility - see src/unfold_utility.h.
// buildFullyCorrected now lives in unfold_utility (purity-corrects all ana::nPtBins
// slices via unfold_utility::purityCorrect and reflattens for RooUnfold) - see
// src/unfold_utility.h.

// Rainbow gradient from red (i=0) to purple (i=n-1) - HSV hue 0 is red, 270 is
// violet/purple; sweeping only that range (not the full 360, which would wrap back to
// red) gives the ROYGBIV ordering rather than a color wheel.
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

  // Response matrix: the full, cross-section-weighted combination of Photon5/10/20
  // (type=1, isample=-1 default), same construction as drawing/draw_purity_corrected.C -
  // this is the response that would actually be used to unfold the real result, unlike
  // Data's own trivially-diagonal stored response (unfolder.cc sets truth=reco as a
  // placeholder for non-MC).
  TH1D * respRecoTemplate  = d.get(Form("hrecoxj%i",ir), 1);
  TH1D * respTruthTemplate = d.get(Form("htruthxj%i",ir), 1);
  TH2D * respMatrix2D      = d.get2d(Form("hxjresponse%i",ir), 1);
  RooUnfoldResponse * response = new RooUnfoldResponse(respRecoTemplate, respTruthTemplate, respMatrix2D);

  // Measured: Data's purity-corrected xJ spectrum (region A minus the two-purity
  // background estimate, per pT bin - see buildFullyCorrected above / unfold_utility::purityCorrect),
  // matching what drawing/draw_purity_corrected.C actually feeds into unfolding for the
  // real result. Unfolding raw region-A reco (background and all) would answer a
  // different question than the one this scan is meant to inform. Truth reference: the
  // same fixed #gamma+jet MC truth the response was built from - Data has no truth of
  // its own.
  TH1D * flatA = d.get(Form("hrecoxj%i_0",ir), 0);
  TH1D * flatC = d.get(Form("hrecoxj%i_2",ir), 0);
  TH1D * flatMeasured = unfold_utility::buildFullyCorrected(flatA, flatC, "data", systag);
  TH1D * flatTruth = respTruthTemplate;

  // Per iteration: error = mean fractional uncertainty of the unfolded result itself,
  // tracking noise amplification directly. (The convergence/"bias" metric is computed
  // separately below, once all iterations are unfolded - see pairBiasScore.)
  vector<int> iters;
  vector<double> errorScore;
  vector<vector<TH1D*>> unfoldedByIter; // [iterIdx][ipt]

  for (unsigned k = 0; k < iterationsToTest.size(); k++) {
    int iter = iterationsToTest[k];
    // unfold_utility::unfoldOnce clones Hreco() before returning - RooUnfoldBayes owns
    // that histogram internally and goes out of scope at the end of this loop body, so a
    // raw Hreco() pointer stored below would otherwise dangle.
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

  // Convergence metric: chi2/ndf of the RELATIVE CHANGE between consecutive iterations,
  // (u_{n+1} - u_n)/u_n per bin, summed over all used pT bins - this asks "has the
  // result stabilized" using only Data's own output at successive iteration counts,
  // rather than comparing to the MC truth reference (which, for real Data, isn't a
  // trustworthy enough "ground truth" to score bias against - see prior discussion).
  // Only defined for a pair of iterations, so the first tested iteration has no partner
  // to compare against and this starts at the SECOND tested iteration.
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

  // Combined score: min-max normalize both metrics to [0,1] and sum, so neither metric's
  // arbitrary absolute scale dominates the other - the minimum is "best" by both at once.
  // Only over pairIters' range (iteration 2 onward), since pairBiasScore isn't defined
  // for the first tested iteration.
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
  int bestIdx = bestPairIdx + 1; // index into iters/unfoldedByIter/errorScore (offset by the dropped first iteration)
  cout << "Data iteration scan: best iteration = " << bestIter
       << " (relative-change chi2=" << pairBiasScore[bestPairIdx] << ", mean frac. error=" << errorScore[bestIdx] << ")" << endl;

  TCanvas * c = new TCanvas("c","",700,900);
  c->SaveAs(Form("%s[", pdfPath.c_str()));

  // Pages 1..N: sanity check - the winning iteration's unfolded Data result vs the
  // #gamma+jet MC truth reference, per used pT bin.
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hUnfold = unfoldedByIter[bestIdx][ipt];
    TH1D * hTruth  = unfold_utility::unflattenXj(flatTruth, ipt, Form("hxjtruth_final_pt%d", ipt));
    // Shape-normalize (bin-width density, then unit area) - Data (raw counts) and the
    // weighted MC truth (cross-section scaled, ~10^9) are on wildly different absolute
    // scales; without this Data's real, non-zero curve is invisible next to truth's.
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

  // Next pages: all iterations overlaid on one plot instead of one page per iteration -
  // colored red (fewest iterations) to purple (most), following the rainbow, so the
  // shape change with more iterations is visible directly rather than only in the
  // chi2/error summary. One page per used pT bin. Shape-normalized (bin-width density,
  // then unit area) for the same reason as the per-pT-bin pages above - Data and the
  // weighted MC truth are on wildly different absolute scales.
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

  // Last plot: chi2/ndf of the relative change between consecutive iterations (top) and
  // unfolded uncertainty (bottom) vs iteration count, both marking the combined-best
  // iteration. Both panels share the same x-range (iteration 2 onward) since the
  // relative-change metric has no value at the first tested iteration.
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
