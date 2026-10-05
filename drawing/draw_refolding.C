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

// Refolding closure test: take the nominal unfolded result and fold it back FORWARD
// through the same response matrix (RooUnfoldResponse::ApplyToTruth), then compare that
// prediction to the purity-corrected measured spectrum that was actually unfolded.
// Complementary to the other cross-checks in this directory: draw_iteration_halfclosure.C
// asks whether Data's unfolded shape converges as niter increases, and
// toy_resp_iterations.C/toy_data_iterations.C ask how much statistical noise the response
// matrix / Data itself contribute to that result - none of those ever compare back to the
// reco-level measurement, so none would catch a broken round-trip (a bad response matrix,
// an inconsistent purity correction, or a RooUnfold setup mistake). This does.
//
// ApplyToTruth's output carries no error of its own ("Errors not set, since we assume
// original truth has no errors" - RooUnfoldResponse.cxx) - it's a deterministic forward
// projection, not an independent measurement - so the only well-defined chi2 here uses the
// measured spectrum's own (purity-corrected) statistical error as the yardstick, rather
// than combining two uncertain errors.
//
// Not shape-normalized like most other comparisons in this directory: refolded and
// measured are already on the same absolute (reco-level, raw-count) scale by construction,
// so an absolute-count comparison IS the actual test here, not just a shape comparison.
//
// A perfect round-trip isn't expected even in a fully self-consistent chain: finite
// Bayesian iterations regularize (deliberately) rather than exactly inverting the
// response, and any truth-prior mismatch shows up here as residual bias. A residual
// comparable to or larger than the measured statistical uncertainty is the actual red flag.
//
// Fakes: RooUnfoldResponse computes Hfakes() once, when the response matrix is built, as
// (training reco template) - (projection of the response matrix onto reco) - i.e. reco-level
// content with no truth-level partner inside the binned truth range (out-of-acceptance
// migration, unmatched jets). ApplyToTruth() only applies the truth->reco migration matrix
// and never adds Hfakes() back, so it has to be added in by hand below.
//
// Hfakes() itself is on the MC TRAINING sample's own absolute (cross-section-weighted)
// scale, NOT Data's - it must NOT be added to the refolded prediction as-is. Checking how
// RooUnfoldBayes itself actually uses fakes (RooUnfoldBayes.cxx setup()/unfold()) confirms
// why: Hfakes() only ever enters as one more competing Bayesian "cause", contributing a
// TRAINING PROBABILITY (Nji/nCi, dimensionless) that then gets applied multiplicatively to
// the REAL measured input being unfolded - Data's own observed counts get reallocated
// between genuine truth causes and "fakes" in proportion to the training fakes FRACTION,
// never by adding/subtracting the MC's own absolute fakes count. Reproduced the same way
// below: fakes fraction = Hfakes()/Hmeasured() (both on the training's own scale, so the
// ratio is dimensionless), applied to Data's own measured content bin-by-bin to get a
// Data-scale-consistent predicted fakes contribution, added into flatRefolded so every
// "refolded" curve/metric from here on already includes it.
//
// A second, pure-MC check is also done below (flatRefoldedTruth): fold the ACTUAL Pythia
// truth-level spectrum (respTruthTemplate) forward through the response and compare to
// the ACTUAL Pythia reco-level spectrum (respRecoTemplate) - both from the exact same
// training sample, no Data, no purity correction, no unfolding, no niter involved at all.
// If this doesn't close, the response matrix itself has a structural problem, independent
// of anything Data- or unfolding-specific the rest of this file might otherwise be blamed
// on. Since truth/reco/fakes are all on the SAME MC training scale here, the raw MC-scale
// Hfakes() (flatFakesTraining) is added directly - no Data-rescaling needed this time.

const int ir = 2; // nominal jet radius index (R=0.4)
const int nPtBinsUsed = ana::nPtBinsUsed; // physics analysis only uses ana::ptBins[ana::firstUsedPtBin..]
const int niterate = 2; // matches draw_purity_corrected.C / draw_final_result.C's chosen nominal iteration count
const vector<int> iterationsToScan = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}; // for the bonus niter-dependence page

// The last 3 xJ bins in each pT bin have very low counts, so chi2/NDF here would be
// dominated by their noise rather than genuine refolding residual - excluded from the
// chi2 metric only (still drawn on the comparison pages). Same exclusion as
// draw_purity_corrected.C/draw_iteration_halfclosure.C/toy_resp_iterations.C/
// toy_data_iterations.C.
const int nXjBinsForChi2 = ana::nUnfoldXjBins - 3;

// densityForDisplay now lives in unfold_utility - see src/unfold_utility.h.
// buildFullyCorrected now lives in unfold_utility (purity-corrects all ana::nPtBins
// slices via unfold_utility::purityCorrect and reflattens for RooUnfold) - see
// src/unfold_utility.h.

// chi2/NDF of refolded vs measured, using ONLY the measured spectrum's own statistical
// error - see file header for why (ApplyToTruth's output has no meaningful error).
double computeChi2NDF(TH1D * hRefolded, TH1D * hMeasured) {
  double chi2 = 0;
  int ndf = 0;
  for (int b = 1; b <= nXjBinsForChi2; b++) {
    double vR = hRefolded->GetBinContent(b);
    double vM = hMeasured->GetBinContent(b);
    double eM = hMeasured->GetBinError(b);
    if (eM <= 0) continue;
    chi2 += pow(vR-vM,2)/(eM*eM);
    ndf++;
  }
  return (ndf > 0) ? chi2/ndf : 0;
}

void draw_refolding(string systag = "nominal") {
  gStyle->SetOptStat(0);

  drawer d("pythia", systag);
  string pdfPath  = Form("%s/pdfs/draw_refolding_%s.pdf", ana::dir(), systag.c_str());
  string rootPath = Form("%s/hists/refolding_%s.root", ana::dir(), systag.c_str());

  // Response matrix: full, cross-section-weighted combination of Photon5/10/20 - same
  // construction as draw_purity_corrected.C/draw_final_result.C. Built explicitly here
  // (rather than via unfold_utility::unfoldOnce's templates+matrix overload) so the SAME
  // RooUnfoldResponse object can also be used for ApplyToTruth below.
  TH1D * respRecoTemplate  = d.get(Form("hrecoxj%i",ir), 1);
  TH1D * respTruthTemplate = d.get(Form("htruthxj%i",ir), 1);
  TH2D * respMatrix2D      = d.get2d(Form("hxjresponse%i",ir), 1);
  RooUnfoldResponse * response = new RooUnfoldResponse(respRecoTemplate, respTruthTemplate, respMatrix2D);

  // Raw MC-training-scale fakes (see file header) - used both to derive the dimensionless
  // training fakes fraction below (for the Data-scale check) and directly (for the
  // pure-MC truth-refolding check, since that one never leaves the MC's own scale).
  TH1D * flatFakesTraining = (TH1D*)response->Hfakes()->Clone("hFakesTrainingFlat");

  // Pure-MC self-consistency check (see file header): fold the actual Pythia truth
  // forward and add back the (unweighted - already on the same scale) training fakes.
  TH1D * flatRefoldedTruth = (TH1D*)response->ApplyToTruth(respTruthTemplate, "hRefoldedTruth");
  flatRefoldedTruth->Add(flatFakesTraining);

  TH1D * flatA = d.get(Form("hrecoxj%i_0",ir), 0);
  TH1D * flatC = d.get(Form("hrecoxj%i_2",ir), 0);
  TH1D * flatMeasured = unfold_utility::buildFullyCorrected(flatA, flatC, "data", systag);

  // Data-scale predicted fakes: (training fakes fraction, dimensionless) x (Data's own
  // measured content), bin by bin - see file header for why this, not the raw MC-scale
  // Hfakes(), is the quantity that belongs on the same footing as flatMeasured/flatRefolded.
  TH1D * flatFakes = (TH1D*)flatMeasured->Clone("hFakesInDataFlat");
  flatFakes->Reset("ICES");
  for (int b = 1; b <= flatFakes->GetNbinsX(); b++) {
    double mesTraining = respRecoTemplate->GetBinContent(b);
    if (mesTraining <= 0) continue;
    double fakeFrac = flatFakesTraining->GetBinContent(b) / mesTraining;
    flatFakes->SetBinContent(b, fakeFrac * flatMeasured->GetBinContent(b));
  }

  // The nominal refolded prediction, fakes included - every "refolded" curve/metric below
  // uses this, not the bare ApplyToTruth() output.
  TH1D * flatUnfolded = unfold_utility::unfoldOnce(response, flatMeasured, niterate, "hUnfoldedNominal");
  TH1D * flatRefolded = (TH1D*)response->ApplyToTruth(flatUnfolded, "hRefoldedNominal");
  flatRefolded->Add(flatFakes);

  TFile * fout = TFile::Open(rootPath.c_str(), "RECREATE");
  TCanvas * c = new TCanvas("c","",700,900);
  c->SaveAs(Form("%s[", pdfPath.c_str()));

  cout << "Refolding closure (niter=" << niterate << "): pT bin, chi2/NDF (first "
       << nXjBinsForChi2 << " of " << ana::nUnfoldXjBins << " xJ bins), fakes integral" << endl;

  // Pages 1..nPtBinsUsed: measured (purity-corrected reco) vs refolded (prediction, fakes
  // included), absolute counts/bin-width (NOT shape-normalized - see file header), with a
  // ratio panel below.
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hMeasured = unfold_utility::unflattenXj(flatMeasured, ipt, Form("hMeasured_pt%d", ipt));
    TH1D * hRefolded = unfold_utility::unflattenXj(flatRefolded, ipt, Form("hRefolded_pt%d", ipt));
    TH1D * hFakes    = unfold_utility::unflattenXj(flatFakes, ipt, Form("hFakes_pt%d", ipt));
    double chi2ndf = computeChi2NDF(hRefolded, hMeasured);
    cout << "  pt" << ipt << " (" << ana::ptBins[ipt] << "-" << ana::ptBins[ipt+1] << " GeV): chi2/NDF = " << chi2ndf
         << ", fakes integral = " << hFakes->Integral() << endl;

    TH1D * hMeasuredDisp = unfold_utility::densityForDisplay(hMeasured, Form("hMeasuredDisp_pt%d", ipt));
    TH1D * hRefoldedDisp = unfold_utility::densityForDisplay(hRefolded, Form("hRefoldedDisp_pt%d", ipt));

    c->Clear();
    c->cd();
    TPad * p1 = new TPad(Form("p1_%d",ipt),"",0,.35,1,1);
    TPad * p2 = new TPad(Form("p2_%d",ipt),"",0,0,1,.35);
    p1->Draw();
    p2->Draw();

    p1->cd();
    p1->SetBottomMargin(0.02);
    p1->SetLeftMargin(.15);
    gPad->SetTicks(1,1);
    hMeasuredDisp->SetLineColor(kBlack);
    hMeasuredDisp->SetMarkerColor(kBlack);
    hMeasuredDisp->SetMarkerStyle(20);
    hMeasuredDisp->SetLineWidth(2);
    hMeasuredDisp->GetXaxis()->SetLabelSize(0);
    hMeasuredDisp->GetXaxis()->SetTitle("");
    hMeasuredDisp->GetYaxis()->SetRangeUser(0, std::max(hMeasuredDisp->GetMaximum(), hRefoldedDisp->GetMaximum())*1.4);
    hMeasuredDisp->Draw("p e");
    // Drawn as a line, not points-with-errors: ApplyToTruth's output (and the fakes
    // contribution added to it) has no meaningful per-bin error (see file header), so error
    // bars here would misrepresent it as an independent measurement rather than a
    // deterministic forward projection.
    hRefoldedDisp->SetLineColor(kRed);
    hRefoldedDisp->SetLineWidth(3);
    hRefoldedDisp->Draw("hist same");
    TLegend * l = new TLegend(.5,.65,.85,.80);
    l->SetLineWidth(0);
    l->SetTextSize(0.032);
    l->AddEntry(hMeasuredDisp, "Measured (purity-corrected)", "lp");
    l->AddEntry(hRefoldedDisp, Form("Refolded + fakes (%d iter.)", niterate), "l");
    l->Draw();
    d.drawAll({"p+p Run24 Data"},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
        Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV, #chi^{2}/NDF = %.2f", ana::JetRs[ir], ana::jet_calib_pt_cut[ir], chi2ndf)}, .18, .85, 14, gPad->GetWh()*0.8);

    p2->cd();
    p2->SetTopMargin(0.02);
    p2->SetBottomMargin(0.3);
    p2->SetLeftMargin(.15);
    gPad->SetTicks(1,1);
    TH1D * hratio = (TH1D*)hRefoldedDisp->Clone(Form("hratio_pt%d", ipt));
    hratio->Divide(hMeasuredDisp);
    hratio->SetLineColor(kBlack);
    hratio->SetMarkerColor(kBlack);
    hratio->SetMarkerStyle(20);
    hratio->GetYaxis()->SetRangeUser(0.5,1.5);
    hratio->GetYaxis()->SetTitle("Refolded / Measured");
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

    fout->cd();
    hMeasured->Write();
    hRefolded->Write();
    hFakes->Write();
    hratio->Write();
    delete hMeasured; delete hRefolded; delete hFakes;
  }

  // Pages (nPtBinsUsed+1)..(2 nPtBinsUsed): pure-MC check - actual Pythia reco vs actual
  // Pythia truth refolded forward (+ training fakes), absolute counts/bin-width (same
  // scale, same reasoning as above). No niter dependence here at all - ApplyToTruth is a
  // one-shot forward projection, not an iterative unfold, so unlike everything else in
  // this file there's nothing to scan over.
  cout << "Pure-MC truth-refolding closure: pT bin, chi2/NDF (first "
       << nXjBinsForChi2 << " of " << ana::nUnfoldXjBins << " xJ bins)" << endl;
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hPythiaReco      = unfold_utility::unflattenXj(respRecoTemplate, ipt, Form("hPythiaReco_pt%d", ipt));
    TH1D * hRefoldedTruthPt = unfold_utility::unflattenXj(flatRefoldedTruth, ipt, Form("hRefoldedTruth_pt%d", ipt));
    double chi2ndf = computeChi2NDF(hRefoldedTruthPt, hPythiaReco);
    cout << "  pt" << ipt << " (" << ana::ptBins[ipt] << "-" << ana::ptBins[ipt+1] << " GeV): chi2/NDF = " << chi2ndf << endl;

    TH1D * hPythiaRecoDisp      = unfold_utility::densityForDisplay(hPythiaReco, Form("hPythiaRecoDisp_pt%d", ipt));
    TH1D * hRefoldedTruthDisp   = unfold_utility::densityForDisplay(hRefoldedTruthPt, Form("hRefoldedTruthDisp_pt%d", ipt));

    c->Clear();
    c->cd();
    TPad * p1 = new TPad(Form("pt1_%d",ipt),"",0,.35,1,1);
    TPad * p2 = new TPad(Form("pt2_%d",ipt),"",0,0,1,.35);
    p1->Draw();
    p2->Draw();

    p1->cd();
    p1->SetBottomMargin(0.02);
    p1->SetLeftMargin(.15);
    gPad->SetTicks(1,1);
    hPythiaRecoDisp->SetLineColor(kBlack);
    hPythiaRecoDisp->SetMarkerColor(kBlack);
    hPythiaRecoDisp->SetMarkerStyle(20);
    hPythiaRecoDisp->SetLineWidth(2);
    hPythiaRecoDisp->GetXaxis()->SetLabelSize(0);
    hPythiaRecoDisp->GetXaxis()->SetTitle("");
    hPythiaRecoDisp->GetYaxis()->SetRangeUser(0, std::max(hPythiaRecoDisp->GetMaximum(), hRefoldedTruthDisp->GetMaximum())*1.4);
    hPythiaRecoDisp->Draw("p e");
    hRefoldedTruthDisp->SetLineColor(kRed);
    hRefoldedTruthDisp->SetLineWidth(3);
    hRefoldedTruthDisp->Draw("hist same");
    TLegend * l = new TLegend(.5,.65,.85,.80);
    l->SetLineWidth(0);
    l->SetTextSize(0.032);
    l->AddEntry(hPythiaRecoDisp, "Pythia8 reco (training)", "lp");
    l->AddEntry(hRefoldedTruthDisp, "Refolded truth + fakes", "l");
    l->Draw();
    d.drawAll({"Pythia8 #gamma+jet MC"},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
        Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV, #chi^{2}/NDF = %.2f", ana::JetRs[ir], ana::jet_calib_pt_cut[ir], chi2ndf)}, .18, .85, 14, gPad->GetWh()*0.8);

    p2->cd();
    p2->SetTopMargin(0.02);
    p2->SetBottomMargin(0.3);
    p2->SetLeftMargin(.15);
    gPad->SetTicks(1,1);
    TH1D * hratio = (TH1D*)hRefoldedTruthDisp->Clone(Form("hratio_truth_pt%d", ipt));
    hratio->Divide(hPythiaRecoDisp);
    hratio->SetLineColor(kBlack);
    hratio->SetMarkerColor(kBlack);
    hratio->SetMarkerStyle(20);
    hratio->GetYaxis()->SetRangeUser(0.5,1.5);
    hratio->GetYaxis()->SetTitle("Refolded truth / Pythia reco");
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

    fout->cd();
    hPythiaReco->Write();
    hRefoldedTruthPt->Write();
    hratio->Write();
    delete hPythiaReco; delete hRefoldedTruthPt;
  }

  // Last page: refolding chi2/NDF (mean over used pT bins) vs iteration count - does the
  // round trip get better or worse with more/fewer iterations? Mirrors the niter-dependence
  // scans in draw_iteration_halfclosure.C/toy_resp_iterations.C/toy_data_iterations.C,
  // using this file's own refolding residual (refolded+fakes vs measured) as the metric
  // instead. Fakes is niter-independent (a fixed property of the response matrix), so it's
  // added once per iteration point here, not recomputed.
  cout << "Refolding niter-dependence scan..." << endl;
  TGraph * gChi2 = new TGraph((int)iterationsToScan.size());
  for (unsigned k = 0; k < iterationsToScan.size(); k++) {
    int iter = iterationsToScan[k];
    TH1D * hUnfoldedIter = unfold_utility::unfoldOnce(response, flatMeasured, iter, Form("hUnfoldedIter_%d", iter));
    TH1D * hRefoldedIter = (TH1D*)response->ApplyToTruth(hUnfoldedIter, Form("hRefoldedIter_%d", iter));
    hRefoldedIter->Add(flatFakes);
    double chi2Sum = 0;
    int nCounted = 0;
    for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
      TH1D * hM = unfold_utility::unflattenXj(flatMeasured, ipt, Form("hMeasuredScan_%d_%d", iter, ipt));
      TH1D * hR = unfold_utility::unflattenXj(hRefoldedIter, ipt, Form("hRefoldedScan_%d_%d", iter, ipt));
      for (int b = 1; b <= nXjBinsForChi2; b++) {
        double vR = hR->GetBinContent(b);
        double vM = hM->GetBinContent(b);
        double eM = hM->GetBinError(b);
        if (eM <= 0) continue;
        chi2Sum += pow(vR-vM,2)/(eM*eM);
        nCounted++;
      }
      delete hM; delete hR;
    }
    double chi2ndf = nCounted > 0 ? chi2Sum/nCounted : 0;
    gChi2->SetPoint(k, iter, chi2ndf);
    cout << "  niter=" << iter << ": refolding chi2/NDF = " << chi2ndf << endl;
    delete hUnfoldedIter;
    delete hRefoldedIter;
  }

  c->Clear();
  c->cd();
  gPad->SetTicks(1,1);
  gPad->SetLeftMargin(.15);
  gChi2->SetMarkerStyle(20);
  gChi2->SetMarkerColor(kBlack);
  gChi2->SetLineColor(kBlack);
  gChi2->SetLineWidth(2);
  gChi2->GetXaxis()->SetTitle("Bayesian unfolding iterations");
  gChi2->GetYaxis()->SetTitle("#chi^{2}/NDF (Refolded+Fakes vs Measured)");
  double ymax = 0;
  for (int k = 0; k < gChi2->GetN(); k++) { double x,y; gChi2->GetPoint(k,x,y); ymax = std::max(ymax,y); }
  gChi2->SetMinimum(0);
  gChi2->SetMaximum(ymax*1.3);
  gChi2->Draw("APL");
  TLine * lnom = new TLine(niterate, 0, niterate, ymax*1.3);
  lnom->SetLineStyle(9);
  lnom->SetLineColor(kRed);
  lnom->Draw("same");
  d.drawAll({"p+p Run24 Data"},{Form("Jet R=%.1f",ana::JetRs[ir]),
      Form("Nominal: %d iterations (dashed line)",niterate)}, .5, .85, 16, 700);
  c->SaveAs(pdfPath.c_str());
  fout->cd();
  gChi2->Write("gRefoldingChi2");

  c->SaveAs(Form("%s]", pdfPath.c_str()));
  fout->Close();
  cout << "Done. Wrote " << pdfPath << " and " << rootPath << endl;
}
