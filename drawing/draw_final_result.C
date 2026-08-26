#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/unfold_utility.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
// The original gammajet project builds its OWN, differently-laid-out drawer/ana classes
// into /home/samson72/root/lib/libgammajet.so, sitting on the same library search path
// as this project's libgammajet_unfold.so. Without forcing which one loads first, ROOT's
// implicit symbol autoload can bind drawer/ana calls to the wrong (mismatched-layout)
// library and segfault - explicit load removes the ambiguity.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Final result: nominal unfolded (1/N)dN/dxJ per pT bin, with statistical and systematic
// uncertainty shown separately (standard practice - stat is uncorrelated bin-to-bin,
// syst is often correlated, so combining them into one band would obscure that), compared
// to the #gamma+jet MC truth prediction.
//
// Statistical uncertainty: the bin errors already carried on the unfolded histogram
// itself, from RooUnfoldBayes's own analytic covariance calculation (computed for every
// unfoldOnce() call - see the "Calculating covariances..." log line). toy_resp_iterations.C
// / toy_data_iterations.C exist specifically to cross-check this analytic covariance
// against a toy-based estimate; they're a validation of this number, not a replacement
// for it, so no toy loops are re-run here.
//
// Systematic uncertainty: the quadrature-sum fractional uncertainty (hquadsum_up/down_pt<N>)
// already computed and saved by draw_systematics.C - JERhigh/JERlow/emscale_high/
// emscale_low/EMRhigh/EMRlow/jes_high/jes_low/threejet/Purity (narrowBDT/narrowISO/
// narrowBDTbkg/narrowISObkg/wideISObkg's own combined ABCD sideband-boundary systematic,
// following PPG12's treatment - see draw_systematics.C's purityMembers)/herwig
// (generator-modeling)/priorSensitivity (data-informed prior reweighting - see
// drawing/draw_prior_sensitivity.C) combined in quadrature. Asymmetric: JER/JES/emscale/
// EMR are true two-point (high/low) systematics there, combined via their per-bin
// envelope rather than symmetrized (Purity is symmetric, like threejet/herwig/
// priorSensitivity - see draw_systematics.C's purityMembers), so the box below is a
// TGraphAsymmErrors, not a plain TH1 (which can only carry one symmetric error per bin).
//
// Statistical uncertainty stays on the data points' own error bars; systematic
// uncertainty is drawn as a separate azure box behind them. MC truth is overlaid for
// reference, with a ratio panel below.

// Jet radius index - mutable (not const) so draw_final_result(int) can set it at the top
// of the function, before any of the code below (all written against this global) runs.
// Defaults to the nominal R=0.4 working point used throughout the note.
int ir = 2;
const int nPtBinsUsed = ana::nPtBinsUsed; // physics analysis only uses ana::ptBins[ana::firstUsedPtBin..]
const int niterate = 2; // matches draw_purity_corrected.C's best-iteration scan result
// Set inside draw_final_result(int) from ir - the nominal R=0.4 default reproduces the
// unsuffixed filenames every other macro/main.tex reads; every other radius gets its own
// _<rname>-suffixed triple instead of clobbering the nominal file. systRootPath must match
// whatever suffix draw_systematics(int) was run with for the same radius.
string pdfPathStr, rootPathStr, systRootPathStr;
const char * pdfPath;
const char * rootPath;
const char * systRootPath;

// densityForDisplay now lives in unfold_utility - see src/unfold_utility.h.
// buildFullyCorrected now lives in unfold_utility (purity-corrects all ana::nPtBins
// slices via unfold_utility::purityCorrect and reflattens for RooUnfold) - see
// src/unfold_utility.h.

void draw_final_result(int jetRadiusIndex = 2) {
  gStyle->SetOptStat(0);

  ir = jetRadiusIndex;
  pdfPathStr      = (ir == 2) ? "/home/samson72/sphnx/gammajet_unfold/pdfs/final_result.pdf"
                              : Form("/home/samson72/sphnx/gammajet_unfold/pdfs/final_result_%s.pdf", ana::rnames[ir]);
  rootPathStr     = (ir == 2) ? "/home/samson72/sphnx/gammajet_unfold/hists/final_result.root"
                              : Form("/home/samson72/sphnx/gammajet_unfold/hists/final_result_%s.root", ana::rnames[ir]);
  systRootPathStr = (ir == 2) ? "/home/samson72/sphnx/gammajet_unfold/hists/systematics.root"
                              : Form("/home/samson72/sphnx/gammajet_unfold/hists/systematics_%s.root", ana::rnames[ir]);
  pdfPath      = pdfPathStr.c_str();
  rootPath     = rootPathStr.c_str();
  systRootPath = systRootPathStr.c_str();

  TFile * fsyst = TFile::Open(systRootPath);
  if (!fsyst || fsyst->IsZombie()) {
    cout << "ERROR: couldn't open " << systRootPath << " - run drawing/draw_systematics.C first." << endl;
    return;
  }

  TFile * fout = TFile::Open(rootPath, "RECREATE");
  TCanvas * c = new TCanvas("c","",700,900);
  c->SaveAs(Form("%s[", pdfPath));

  drawer d("pythia", "nominal");
  TH1D * respRecoTemplate  = d.get(Form("hrecoxj%i",ir), 1);
  TH1D * respTruthTemplate = d.get(Form("htruthxj%i",ir), 1);
  TH2D * respMatrix2D      = d.get2d(Form("hxjresponse%i",ir), 1);
  TH1D * flatA = d.get(Form("hrecoxj%i_0",ir), 0);
  TH1D * flatC = d.get(Form("hrecoxj%i_2",ir), 0);
  TH1D * flatCorrected = unfold_utility::buildFullyCorrected(flatA, flatC, "data", "nominal");
  TH1D * flatUnfolded  = unfold_utility::unfoldOnce(respRecoTemplate, respTruthTemplate, respMatrix2D, flatCorrected, niterate, "hUnfoldedFinal");

  cout << "pT bin, xJ bin: nominal, stat frac, syst frac, total frac" << endl;

  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hUnfold     = unfold_utility::unflattenXj(flatUnfolded, ipt, Form("hFinal_pt%d", ipt));
    TH1D * hUnfoldDisp = unfold_utility::densityForDisplay(hUnfold, Form("hFinalDisp_pt%d", ipt));
    hUnfoldDisp->Scale(1./hUnfoldDisp->Integral());
    hUnfoldDisp->GetYaxis()->SetTitle("(1/N) dN/dx_{J#gamma}");

    TH1D * hTruth     = unfold_utility::unflattenXj(respTruthTemplate, ipt, Form("hFinalTruth_pt%d", ipt));
    TH1D * hTruthDisp = unfold_utility::densityForDisplay(hTruth, Form("hFinalTruthDisp_pt%d", ipt));
    hTruthDisp->Scale(1./hTruthDisp->Integral());

    // Systematic-only uncertainty box - kept separate from the data points' own
    // statistical error bars (standard practice: stat and syst behave differently -
    // stat is uncorrelated bin-to-bin, syst is often correlated - so they're shown
    // separately rather than combined into one band). hquadsum_up/down_pt<N> are the
    // fractional systematic (from draw_systematics.C's shape-normalized comparison, so
    // already in the same "fraction of shape-normalized content" units used here) -
    // asymmetric because JER/JES/emscale are combined as true two-point (high/low)
    // systematics there, not symmetrized. A TGraphAsymmErrors carries the box instead of
    // a TH1 (whose SetBinError is inherently one symmetric value per bin).
    TH1D * hquadUp   = (TH1D*)fsyst->Get(Form("hquadsum_up_pt%d", ipt));
    TH1D * hquadDown = (TH1D*)fsyst->Get(Form("hquadsum_down_pt%d", ipt));
    int nb = hUnfoldDisp->GetNbinsX();
    TGraphAsymmErrors * gSystBox = new TGraphAsymmErrors(nb);
    gSystBox->SetName(Form("gFinalSystBox_pt%d", ipt));
    for (int b = 1; b <= nb; b++) {
      double content = hUnfoldDisp->GetBinContent(b);
      double x = hUnfoldDisp->GetBinCenter(b);
      double halfw = hUnfoldDisp->GetBinWidth(b) / 2.;
      double statFrac    = content > 0 ? hUnfoldDisp->GetBinError(b) / content : 0;
      double systFracUp   = (content > 0 && hquadUp)   ? hquadUp->GetBinContent(b)   : 0;
      double systFracDown = (content > 0 && hquadDown) ? hquadDown->GetBinContent(b) : 0;
      gSystBox->SetPoint(b-1, x, content);
      gSystBox->SetPointError(b-1, halfw, halfw, systFracDown*content, systFracUp*content);
      if (content > 0)
        cout << "  pt" << ipt << " bin" << b << ": " << content
             << ", stat=" << statFrac << ", syst up=" << systFracUp << ", syst down=" << systFracDown
             << ", total up=" << sqrt(statFrac*statFrac + systFracUp*systFracUp)
             << ", total down=" << sqrt(statFrac*statFrac + systFracDown*systFracDown) << endl;
    }

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
    // Frame-only histogram for the axes - TGraphAsymmErrors carries no axis-range/label
    // machinery of its own, so a zero-content clone sets up the pad exactly as the old
    // hSystBox draw used to before the box itself is drawn on top.
    TH1D * hFrame1 = (TH1D*)hUnfoldDisp->Clone(Form("hFinalFrame1_pt%d", ipt));
    hFrame1->Reset("ICES");
    hFrame1->SetLineColor(kWhite);
    hFrame1->GetXaxis()->SetLabelSize(0);
    hFrame1->GetXaxis()->SetTitle("");
    hFrame1->GetYaxis()->SetRangeUser(0, std::max(hUnfoldDisp->GetMaximum(), hTruthDisp->GetMaximum())*1.5);
    hFrame1->Draw("p");
    gSystBox->SetFillColorAlpha(kAzure+1, 0.35);
    gSystBox->SetLineColor(kWhite);
    gSystBox->Draw("2 same");
    hUnfoldDisp->SetLineColor(kBlack);
    hUnfoldDisp->SetMarkerColor(kBlack);
    hUnfoldDisp->SetMarkerStyle(20);
    hUnfoldDisp->SetLineWidth(2);
    hUnfoldDisp->Draw("p e same");
    hTruthDisp->SetLineColor(kRed);
    hTruthDisp->SetMarkerColor(kRed);
    hTruthDisp->SetMarkerStyle(24);
    hTruthDisp->SetLineWidth(2);
    hTruthDisp->Draw("p e same");
    TLegend * l = new TLegend(.55,.62,.85,.85);
    l->SetLineWidth(0);
    l->SetTextSize(0.032);
    l->AddEntry(hUnfoldDisp, "p+p Run24 Data", "lp");
    l->AddEntry(gSystBox, "Systematic uncertainty", "f");
    l->AddEntry(hTruthDisp, "Pythia8 #gamma+jet (truth)", "lp");
    l->Draw();
    d.drawAll({"p+p Run24 Data"},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
        Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV", ana::JetRs[ir], ana::jet_calib_pt_cut[ir])}, .18, .85, 14, gPad->GetWh()*0.8);

    p2->cd();
    p2->SetTopMargin(0.02);
    p2->SetBottomMargin(0.3);
    p2->SetLeftMargin(.15);
    gPad->SetTicks(1,1);
    // Ratio box: same x/y/errors as gSystBox, each scaled by that bin's own truth content
    // (not truth's own uncertainty - the box is meant to carry only the analysis's
    // systematic uncertainty, same quantity as the top panel, just rescaled to sit on the
    // Data/Truth axis).
    TGraphAsymmErrors * gRatioSystBox = new TGraphAsymmErrors(nb);
    gRatioSystBox->SetName(Form("gFinalRatioSystBox_pt%d", ipt));
    for (int b = 1; b <= nb; b++) {
      double truth = hTruthDisp->GetBinContent(b);
      double x, y;
      gSystBox->GetPoint(b-1, x, y);
      double exl = gSystBox->GetErrorXlow(b-1);
      double exh = gSystBox->GetErrorXhigh(b-1);
      double eyl = gSystBox->GetErrorYlow(b-1);
      double eyh = gSystBox->GetErrorYhigh(b-1);
      double ratioY = truth > 0 ? y/truth : 0;
      gRatioSystBox->SetPoint(b-1, x, ratioY);
      gRatioSystBox->SetPointError(b-1, exl, exh, truth > 0 ? eyl/truth : 0, truth > 0 ? eyh/truth : 0);
    }
    TH1D * hRatio = (TH1D*)hUnfoldDisp->Clone(Form("hFinalRatio_pt%d", ipt));
    hRatio->Divide(hTruthDisp);

    TH1D * hFrame2 = (TH1D*)hRatio->Clone(Form("hFinalFrame2_pt%d", ipt));
    hFrame2->Reset("ICES");
    hFrame2->SetLineColor(kWhite);
    hFrame2->GetYaxis()->SetRangeUser(0.5,1.5);
    hFrame2->GetYaxis()->SetTitle("Data / Pythia Truth");
    hFrame2->GetYaxis()->SetTitleSize(0.09);
    hFrame2->GetYaxis()->SetTitleOffset(0.7);
    hFrame2->GetYaxis()->SetLabelSize(0.08);
    hFrame2->GetXaxis()->SetTitle("x_{J#gamma}");
    hFrame2->GetXaxis()->SetTitleSize(0.09);
    hFrame2->GetXaxis()->SetLabelSize(0.08);
    hFrame2->Draw("p");
    gRatioSystBox->SetFillColorAlpha(kAzure+1, 0.35);
    gRatioSystBox->SetLineColor(kWhite);
    gRatioSystBox->Draw("2 same");
    hRatio->SetLineColor(kBlack);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetMarkerStyle(20);
    hRatio->Draw("p e same");
    TLine * line = new TLine(ana::unfoldXjBins[0],1,ana::unfoldXjBins[ana::nUnfoldXjBins],1);
    line->SetLineStyle(9);
    line->Draw("same");
    c->SaveAs(pdfPath);

    fout->cd();
    hUnfoldDisp->Write();
    gSystBox->Write();
    hTruthDisp->Write();
  }

  c->SaveAs(Form("%s]", pdfPath));
  fout->Close();
  cout << "Done. Wrote " << pdfPath << " and " << rootPath << endl;
}
