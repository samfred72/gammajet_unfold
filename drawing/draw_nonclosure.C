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

// Non-closure / prior-dependence systematic: does the unfolding correctly recover a truth
// shape that's DIFFERENT from the one the response matrix's own training sample used as its
// prior? RooUnfoldBayes (Adye, "Unfolding algorithms and tests using RooUnfold", section 3.1)
// starts its Bayes iteration from the training truth as the initial prior rather than a flat
// distribution - iterating is supposed to wash that out, but at a finite, deliberately-chosen
// niterate there can be residual "pull" toward the training shape. Neither draw_refolding.C
// (checks the forward+inverse round-trip is self-consistent, using whatever truth the
// unfolding actually produced) nor the niter-dependence scans elsewhere in this directory
// (check convergence/noise using the SAME Data each time) would catch this - both would look
// fine even if unfolding is quietly biased toward the MC prior's shape.
//
// Method, following the ATLAS dijet-xJ paper's own prior-variation study (Phys. Lett. B 774
// (2017) 379, page 6: reweight the truth to remove its peak, build the correspondingly
// reweighted reconstructed distribution from the SAME reweighted events, unfold that using
// the ORIGINAL (un-reweighted) response, and check whether the unfolded result correctly
// comes back looking like the reweighted truth rather than snapping back to the original
// peaked shape). Adapted here to not need event-level access: the response matrix
// (hxjresponse%i) already holds the joint (reco,truth) counts, so reweighting is done by
// scaling that matrix's truth (Y) axis by a per-truth-bin factor and re-projecting onto the
// reco (X) axis - mathematically identical to reweighting the underlying MC events and
// rebuilding both marginal histograms from scratch.
//
// The specific reweighting used here: flatten the training truth to be CONSTANT across the
// nUnfoldXjBins real xJ bins within each pT slice (preserving that slice's own total, so only
// shape - not overall normalization - changes). This is the other canonical prior choice
// discussed in the unfolding literature (Adye again: "RooUnfoldBayes takes the training truth
// as its initial prior, rather than a flat distribution") - i.e. this test asks what happens
// if the unfolding had instead been handed the OTHER standard prior choice, without needing
// to hand-tune a distortion shape.
//
// Fakes (see draw_refolding.C for the full derivation) have no truth-level partner, so they
// aren't reweighted by the truth-shape change - they're added back to the alternate reco
// unweighted, at their nominal (MC training) absolute scale. Unlike draw_refolding.C, no
// Data-scale rescaling is needed anywhere in this file: every histogram here (nominal and
// alternate truth/reco/unfolded) stays entirely on the MC training sample's own absolute
// scale - real Data is never touched.
//
// These alternate histograms are deterministic reweightings/projections, not independent
// measurements, so - matching draw_refolding.C's treatment of ApplyToTruth's output - their
// bin errors are only ever an approximate sqrt(content) stand-in, not rigorously propagated.

const int ir = 2; // nominal jet radius index (R=0.4)
const int nPtBinsUsed = ana::nPtBinsUsed; // physics analysis only uses ana::ptBins[ana::firstUsedPtBin..]
const int niterate = 2; // matches draw_purity_corrected.C / draw_final_result.C's chosen nominal iteration count
const vector<int> iterationsToScan = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}; // for the bonus niter-dependence page

// The last 3 xJ bins in each pT bin have very low counts, so chi2/NDF here would be
// dominated by their noise rather than genuine bias - excluded from the chi2 metric only
// (still drawn on the comparison pages). Same exclusion as draw_purity_corrected.C/
// draw_iteration_halfclosure.C/toy_iterations.C/draw_refolding.C.
const int nXjBinsForChi2 = ana::nUnfoldXjBins - 3;
// densityForDisplay now lives in unfold_utility - see src/unfold_utility.h.

// Sets bin errors to a simple Poisson-like sqrt(content) stand-in - see file header for why
// (these are deterministic reweightings/projections, not real independent measurements).
void setApproxErrors(TH1D * h) {
  for (int b = 0; b <= h->GetNbinsX()+1; b++) h->SetBinError(b, sqrt(std::max(h->GetBinContent(b), 0.)));
}

// Per flattened bin, the factor that flattens the training truth to a constant value across
// the nUnfoldXjBins real xJ bins within each pT slice, preserving that slice's own total -
// see file header. The 2 reserved under/overflow-xJ slots per pT bin (see unfold_utility.h)
// are left at weight=1 (untouched) - they aren't among the "real" xJ bins being reshaped.
TH1D * buildFlatteningWeights(TH1D * truthFlat) {
  TH1D * w = (TH1D*)truthFlat->Clone("hFlatteningWeights");
  for (int b = 0; b <= w->GetNbinsX()+1; b++) w->SetBinContent(b, 1.0);
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
    TH1D * truthPt = unfold_utility::unflattenXj(truthFlat, ipt, Form("hTruthPt_tmp_%d", ipt));
    double target = truthPt->Integral() / ana::nUnfoldXjBins;
    for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) {
      int flatbin = ipt*(ana::nUnfoldXjBins+2) + ixj + 1;
      double content = truthPt->GetBinContent(ixj+1);
      w->SetBinContent(flatbin+1, (content > 0) ? target/content : 0.);
    }
    delete truthPt;
  }
  return w;
}

// chi2/NDF of the unfolded result vs its target truth, using the unfolded result's OWN
// (RooUnfoldBayes-propagated) statistical error as the yardstick - i.e. is the prior-driven
// bias big or small compared to how uncertain the unfolded result already is, the same
// question ATLAS's delta_prior/delta_stat combination asks.
double computeChi2NDF(TH1D * hUnfolded, TH1D * hTruth) {
  double chi2 = 0;
  int ndf = 0;
  for (int b = 1; b <= nXjBinsForChi2; b++) {
    double vU = hUnfolded->GetBinContent(b);
    double vT = hTruth->GetBinContent(b);
    double eU = hUnfolded->GetBinError(b);
    if (eU <= 0) continue;
    chi2 += pow(vU-vT,2)/(eU*eU);
    ndf++;
  }
  return (ndf > 0) ? chi2/ndf : 0;
}

// Shared two-panel (overlay + ratio) page, used for all three comparison page sets below.
void drawComparisonPage(TCanvas * c, const char * pdfPath, drawer & d, int ipt,
    TH1D * hA, const char * labelA, TH1D * hB, const char * labelB, const char * ratioTitle,
    vector<string> extraLines) {
  TH1D * hAdisp = unfold_utility::densityForDisplay(hA, Form("%s_disp", hA->GetName()));
  TH1D * hBdisp = unfold_utility::densityForDisplay(hB, Form("%s_disp", hB->GetName()));

  c->Clear();
  c->cd();
  TPad * p1 = new TPad(Form("p1_%s",hA->GetName()),"",0,.35,1,1);
  TPad * p2 = new TPad(Form("p2_%s",hA->GetName()),"",0,0,1,.35);
  p1->Draw();
  p2->Draw();

  p1->cd();
  p1->SetBottomMargin(0.02);
  p1->SetLeftMargin(.15);
  gPad->SetTicks(1,1);
  hAdisp->SetLineColor(kBlack);
  hAdisp->SetMarkerColor(kBlack);
  hAdisp->SetMarkerStyle(20);
  hAdisp->SetLineWidth(2);
  hAdisp->GetXaxis()->SetLabelSize(0);
  hAdisp->GetXaxis()->SetTitle("");
  hAdisp->GetYaxis()->SetRangeUser(0, std::max(hAdisp->GetMaximum(), hBdisp->GetMaximum())*1.4);
  hAdisp->Draw("p e");
  hBdisp->SetLineColor(kRed);
  hBdisp->SetMarkerColor(kRed);
  hBdisp->SetMarkerStyle(21);
  hBdisp->SetLineWidth(2);
  hBdisp->Draw("p e same");
  TLegend * l = new TLegend(.5,.65,.85,.80);
  l->SetLineWidth(0);
  l->SetTextSize(0.032);
  l->AddEntry(hAdisp, labelA, "lp");
  l->AddEntry(hBdisp, labelB, "lp");
  l->Draw();
  vector<string> lines = {Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
      Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV", ana::JetRs[ir], ana::jet_calib_pt_cut[ir])};
  for (auto & s : extraLines) lines.push_back(s);
  d.drawAll({"Pythia8 #gamma+jet MC"}, lines, .18, .85, 14, gPad->GetWh()*0.8);

  p2->cd();
  p2->SetTopMargin(0.02);
  p2->SetBottomMargin(0.3);
  p2->SetLeftMargin(.15);
  gPad->SetTicks(1,1);
  TH1D * hratio = (TH1D*)hBdisp->Clone(Form("%s_ratio", hB->GetName()));
  hratio->Divide(hAdisp);
  hratio->SetLineColor(kBlack);
  hratio->SetMarkerColor(kBlack);
  hratio->SetMarkerStyle(20);
  hratio->GetYaxis()->SetRangeUser(0.5,1.5);
  hratio->GetYaxis()->SetTitle(ratioTitle);
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
  c->SaveAs(pdfPath);
  delete hAdisp; delete hBdisp; delete hratio;
}

void draw_nonclosure() {
  gStyle->SetOptStat(0);

  drawer d("pythia", "nominal");
  string pdfPath  = ana::path("pdfs/draw_nonclosure.pdf");
  string rootPath = ana::path("hists/nonclosure.root");

  // Nominal response matrix - full, cross-section-weighted combination of Photon5/10/20,
  // same construction as draw_refolding.C/draw_purity_corrected.C.
  TH1D * respRecoTemplate  = d.get(Form("hrecoxj%i",ir), 1);
  TH1D * respTruthTemplate = d.get(Form("htruthxj%i",ir), 1);
  TH2D * respMatrix2D      = d.get2d(Form("hxjresponse%i",ir), 1);
  RooUnfoldResponse * response = new RooUnfoldResponse(respRecoTemplate, respTruthTemplate, respMatrix2D);
  TH1D * flatFakes = (TH1D*)response->Hfakes()->Clone("hFakesFlat"); // MC-training scale - no rescaling needed here, see file header.

  // ---- Build the alternate (flattened-prior) truth, response, and reco ----
  TH1D * weights = buildFlatteningWeights(respTruthTemplate);

  TH1D * altTruthTemplate = (TH1D*)respTruthTemplate->Clone("hAltTruth");
  for (int b = 1; b <= altTruthTemplate->GetNbinsX(); b++)
    altTruthTemplate->SetBinContent(b, respTruthTemplate->GetBinContent(b) * weights->GetBinContent(b));
  setApproxErrors(altTruthTemplate);

  TH2D * altRespMatrix2D = (TH2D*)respMatrix2D->Clone("hAltResponse2D");
  for (int by = 1; by <= altRespMatrix2D->GetNbinsY(); by++) {
    double wgt = weights->GetBinContent(by);
    for (int bx = 1; bx <= altRespMatrix2D->GetNbinsX(); bx++)
      altRespMatrix2D->SetBinContent(bx, by, respMatrix2D->GetBinContent(bx,by) * wgt);
  }

  TH1D * altRecoMatched = (TH1D*)respRecoTemplate->Clone("hAltRecoMatched");
  for (int bx = 1; bx <= altRespMatrix2D->GetNbinsX(); bx++) {
    double sum = 0;
    for (int by = 1; by <= altRespMatrix2D->GetNbinsY(); by++) sum += altRespMatrix2D->GetBinContent(bx,by);
    altRecoMatched->SetBinContent(bx, sum);
  }
  TH1D * altReco = (TH1D*)altRecoMatched->Clone("hAltReco");
  altReco->Add(flatFakes); // unweighted - fakes have no truth-level partner to reweight by
  setApproxErrors(altReco);

  // ---- Unfold both the nominal and alternate reco through the SAME (un-reweighted) response ----
  TH1D * nomUnfolded = unfold_utility::unfoldOnce(response, respRecoTemplate, niterate, "hNomUnfolded");
  TH1D * altUnfolded = unfold_utility::unfoldOnce(response, altReco, niterate, "hAltUnfolded");

  TFile * fout = TFile::Open(rootPath.c_str(), "RECREATE");
  TCanvas * c = new TCanvas("c","",700,900);
  c->SaveAs(Form("%s[", pdfPath.c_str()));

  cout << "Non-closure / prior-dependence test (niter=" << niterate << "):" << endl;
  cout << "pT bin: bias chi2/NDF (unfolded alternate vs true alternate truth, first "
       << nXjBinsForChi2 << " of " << ana::nUnfoldXjBins << " xJ bins)" << endl;

  // ---- Pages 1..nPtBinsUsed: the injected distortion itself - nominal vs flattened truth ----
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hNomTruth = unfold_utility::unflattenXj(respTruthTemplate, ipt, Form("hNomTruth_pt%d", ipt));
    TH1D * hAltTruth = unfold_utility::unflattenXj(altTruthTemplate, ipt, Form("hAltTruthPage1_pt%d", ipt));
    drawComparisonPage(c, pdfPath.c_str(), d, ipt, hNomTruth, "Nominal truth (Pythia8 prior)",
        hAltTruth, "Alternate truth (flattened)", "Alternate / Nominal",
        {"Injected prior distortion"});
    delete hNomTruth; delete hAltTruth;
  }

  // ---- Pages (nPtBinsUsed+1)..(2 nPtBinsUsed): alternate reco vs nominal reco ----
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hNomReco = unfold_utility::unflattenXj(respRecoTemplate, ipt, Form("hNomReco_pt%d", ipt));
    TH1D * hAltRecoPt = unfold_utility::unflattenXj(altReco, ipt, Form("hAltReco_pt%d", ipt));
    drawComparisonPage(c, pdfPath.c_str(), d, ipt, hNomReco, "Nominal reco",
        hAltRecoPt, "Alternate reco (flattened-truth MC)", "Alternate / Nominal",
        {"Reconstructed-level effect of the distortion"});
    delete hNomReco; delete hAltRecoPt;
  }

  // ---- Pages (2 nPtBinsUsed+1)..(3 nPtBinsUsed): unfolded alternate vs unfolded nominal ----
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hNomUnfPt = unfold_utility::unflattenXj(nomUnfolded, ipt, Form("hNomUnf_pt%d", ipt));
    TH1D * hAltUnfPt = unfold_utility::unflattenXj(altUnfolded, ipt, Form("hAltUnf_pt%d", ipt));
    drawComparisonPage(c, pdfPath.c_str(), d, ipt, hNomUnfPt, Form("Unfolded nominal (%d iter.)", niterate),
        hAltUnfPt, Form("Unfolded alternate (%d iter.)", niterate), "Alternate / Nominal",
        {"Both unfolded through the SAME (nominal) response"});
    delete hNomUnfPt; delete hAltUnfPt;
  }

  // ---- Pages (3 nPtBinsUsed+1)..(4 nPtBinsUsed): THE bias check - unfolded alternate vs the
  // KNOWN alternate truth (with nominal truth overlaid as a "what bias would pull toward"
  // reference). This is the actual non-closure/bias systematic. ----
  vector<double> chi2ByPt(ana::nPtBins);
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hAltTruthPt = unfold_utility::unflattenXj(altTruthTemplate, ipt, Form("hAltTruth_pt%d", ipt));
    TH1D * hAltUnfPt   = unfold_utility::unflattenXj(altUnfolded, ipt, Form("hAltUnfBias_pt%d", ipt));
    TH1D * hNomTruthPt = unfold_utility::unflattenXj(respTruthTemplate, ipt, Form("hNomTruthRef_pt%d", ipt));
    double chi2ndf = computeChi2NDF(hAltUnfPt, hAltTruthPt);
    chi2ByPt[ipt] = chi2ndf;
    cout << "  pt" << ipt << " (" << ana::ptBins[ipt] << "-" << ana::ptBins[ipt+1] << " GeV): bias chi2/NDF = " << chi2ndf << endl;

    TH1D * hAltTruthDisp = unfold_utility::densityForDisplay(hAltTruthPt, Form("hAltTruthDisp_pt%d", ipt));
    TH1D * hAltUnfDisp   = unfold_utility::densityForDisplay(hAltUnfPt, Form("hAltUnfDisp_pt%d", ipt));
    TH1D * hNomTruthDisp = unfold_utility::densityForDisplay(hNomTruthPt, Form("hNomTruthRefDisp_pt%d", ipt));

    c->Clear();
    c->cd();
    TPad * p1 = new TPad(Form("pb1_%d",ipt),"",0,.35,1,1);
    TPad * p2 = new TPad(Form("pb2_%d",ipt),"",0,0,1,.35);
    p1->Draw();
    p2->Draw();

    p1->cd();
    p1->SetBottomMargin(0.02);
    p1->SetLeftMargin(.15);
    gPad->SetTicks(1,1);
    hAltTruthDisp->SetLineColor(kBlack);
    hAltTruthDisp->SetMarkerColor(kBlack);
    hAltTruthDisp->SetMarkerStyle(20);
    hAltTruthDisp->SetLineWidth(2);
    hAltTruthDisp->GetXaxis()->SetLabelSize(0);
    hAltTruthDisp->GetXaxis()->SetTitle("");
    hAltTruthDisp->GetYaxis()->SetRangeUser(0, std::max({hAltTruthDisp->GetMaximum(), hAltUnfDisp->GetMaximum(), hNomTruthDisp->GetMaximum()})*1.4);
    hAltTruthDisp->Draw("p e");
    hAltUnfDisp->SetLineColor(kRed);
    hAltUnfDisp->SetMarkerColor(kRed);
    hAltUnfDisp->SetMarkerStyle(21);
    hAltUnfDisp->SetLineWidth(2);
    hAltUnfDisp->Draw("p e same");
    hNomTruthDisp->SetLineColor(kGray+2);
    hNomTruthDisp->SetLineStyle(2);
    hNomTruthDisp->SetLineWidth(2);
    hNomTruthDisp->Draw("hist same");
    TLegend * l = new TLegend(.45,.58,.85,.80);
    l->SetLineWidth(0);
    l->SetTextSize(0.03);
    l->AddEntry(hAltTruthDisp, "True alternate truth (target)", "lp");
    l->AddEntry(hAltUnfDisp, Form("Unfolded alternate (%d iter.)", niterate), "lp");
    l->AddEntry(hNomTruthDisp, "Nominal truth (bias would pull here)", "l");
    l->Draw();
    d.drawAll({"Pythia8 #gamma+jet MC"},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
        Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV, bias #chi^{2}/NDF = %.2f", ana::JetRs[ir], ana::jet_calib_pt_cut[ir], chi2ndf)}, .18, .85, 14, gPad->GetWh()*0.8);

    p2->cd();
    p2->SetTopMargin(0.02);
    p2->SetBottomMargin(0.3);
    p2->SetLeftMargin(.15);
    gPad->SetTicks(1,1);
    TH1D * hratio = (TH1D*)hAltUnfDisp->Clone(Form("hBiasRatio_pt%d", ipt));
    hratio->Divide(hAltTruthDisp);
    hratio->SetLineColor(kBlack);
    hratio->SetMarkerColor(kBlack);
    hratio->SetMarkerStyle(20);
    hratio->GetYaxis()->SetRangeUser(0.5,1.5);
    hratio->GetYaxis()->SetTitle("Unfolded / True alt. truth");
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
    hAltTruthPt->Write(); hAltUnfPt->Write(); hNomTruthPt->Write(); hratio->Write();
    delete hAltTruthPt; delete hAltUnfPt; delete hNomTruthPt;
    delete hAltTruthDisp; delete hAltUnfDisp; delete hNomTruthDisp; delete hratio;
  }

  // ---- Summary page: fractional bias (unfolded alt - true alt truth)/true alt truth vs xJ,
  // all used pT bins overlaid - the size of the systematic at a glance, same style as
  // draw_systematics.C's final overlay page. ----
  c->Clear();
  c->cd();
  vector<TPad*> pads(nPtBinsUsed);
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    double y1 = 1.0 - (ipt+1)*(1.0/nPtBinsUsed);
    double y2 = 1.0 - ipt*(1.0/nPtBinsUsed);
    pads[ipt] = new TPad(Form("psum_%d",ipt),"",0,y1,1,y2);
    pads[ipt]->Draw();
  }
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    // idisplay: 0-based position among the nPtBinsUsed stacked pads (pads[]/"last panel"
    // checks need this), separate from ipt, the real ana::ptBins index (altTruthTemplate/
    // altUnfolded/chi2ByPt/ana::ptBins[] all still need the real index).
    int idisplay = ipt - ana::firstUsedPtBin;
    pads[idisplay]->cd();
    pads[idisplay]->SetLeftMargin(.15);
    pads[idisplay]->SetBottomMargin(idisplay == nPtBinsUsed-1 ? 0.2 : 0.02);
    pads[idisplay]->SetTopMargin(0.05);
    gPad->SetTicks(1,1);

    TH1D * hAltTruthPt = unfold_utility::unflattenXj(altTruthTemplate, ipt, Form("hAltTruthSum_pt%d", ipt));
    TH1D * hAltUnfPt   = unfold_utility::unflattenXj(altUnfolded, ipt, Form("hAltUnfSum_pt%d", ipt));
    TH1D * hfrac = (TH1D*)hAltUnfPt->Clone(Form("hFracBias_pt%d", ipt));
    for (int b = 1; b <= hfrac->GetNbinsX(); b++) {
      double t = hAltTruthPt->GetBinContent(b);
      hfrac->SetBinContent(b, t > 0 ? (hAltUnfPt->GetBinContent(b)-t)/t : 0);
      hfrac->SetBinError(b, 0); // deterministic curve, not an independent measurement - see file header
    }
    double ymax = std::max(0.2, std::max(hfrac->GetMaximum(), -hfrac->GetMinimum()));
    hfrac->GetYaxis()->SetRangeUser(-ymax*1.3, ymax*1.3);
    hfrac->GetYaxis()->SetTitle("(Unfolded-True)/True");
    hfrac->GetYaxis()->SetTitleSize(0.08);
    hfrac->GetYaxis()->SetTitleOffset(0.8);
    hfrac->GetYaxis()->SetLabelSize(0.07);
    hfrac->GetXaxis()->SetTitle(idisplay == nPtBinsUsed-1 ? "x_{J#gamma}" : "");
    hfrac->GetXaxis()->SetLabelSize(idisplay == nPtBinsUsed-1 ? 0.07 : 0);
    hfrac->GetXaxis()->SetTitleSize(0.08);
    hfrac->SetLineColor(kBlue+1);
    hfrac->SetLineWidth(2);
    hfrac->Draw("hist");
    TLine * zero = new TLine(ana::unfoldXjBins[0],0,ana::unfoldXjBins[ana::nUnfoldXjBins],0);
    zero->SetLineStyle(9);
    zero->Draw("same");
    TLatex * t = new TLatex(.18,.85,Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV, #chi^{2}/NDF=%.2f",ana::ptBins[ipt],ana::ptBins[ipt+1],chi2ByPt[ipt]));
    t->SetNDC();
    t->SetTextFont(43);
    t->SetTextSize(16);
    t->Draw();
    fout->cd();
    hfrac->Write();
    delete hAltTruthPt; delete hAltUnfPt;
  }
  c->SaveAs(pdfPath.c_str());

  // ---- Bonus page: bias chi2/NDF vs iteration count - does more regularization reduce the
  // prior-dependence bias (as expected) at the cost of amplifying noise (see the OTHER
  // niter-dependence scans in this directory)? Directly parallel to the ATLAS paper's
  // delta_prior(niter) contribution to its own niter selection. ----
  cout << "Non-closure bias niter-dependence scan..." << endl;
  TGraph * gChi2 = new TGraph((int)iterationsToScan.size());
  for (unsigned k = 0; k < iterationsToScan.size(); k++) {
    int iter = iterationsToScan[k];
    TH1D * hAltUnfIter = unfold_utility::unfoldOnce(response, altReco, iter, Form("hAltUnfIter_%d", iter));
    double chi2Sum = 0;
    int nCounted = 0;
    for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
      TH1D * hU = unfold_utility::unflattenXj(hAltUnfIter, ipt, Form("hAltUnfScan_%d_%d", iter, ipt));
      TH1D * hT = unfold_utility::unflattenXj(altTruthTemplate, ipt, Form("hAltTruthScan_%d_%d", iter, ipt));
      for (int b = 1; b <= nXjBinsForChi2; b++) {
        double eU = hU->GetBinError(b);
        if (eU <= 0) continue;
        chi2Sum += pow(hU->GetBinContent(b)-hT->GetBinContent(b),2)/(eU*eU);
        nCounted++;
      }
      delete hU; delete hT;
    }
    double chi2ndf = nCounted > 0 ? chi2Sum/nCounted : 0;
    gChi2->SetPoint(k, iter, chi2ndf);
    cout << "  niter=" << iter << ": bias chi2/NDF = " << chi2ndf << endl;
    delete hAltUnfIter;
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
  gChi2->GetYaxis()->SetTitle("Bias #chi^{2}/NDF (Unfolded alternate vs true alternate truth)");
  double ymax = 0;
  for (int k = 0; k < gChi2->GetN(); k++) { double x,y; gChi2->GetPoint(k,x,y); ymax = std::max(ymax,y); }
  gChi2->SetMinimum(0);
  gChi2->SetMaximum(ymax*1.3);
  gChi2->Draw("APL");
  TLine * lnom = new TLine(niterate, 0, niterate, ymax*1.3);
  lnom->SetLineStyle(9);
  lnom->SetLineColor(kRed);
  lnom->Draw("same");
  d.drawAll({"Pythia8 #gamma+jet MC"},{Form("Jet R=%.1f",ana::JetRs[ir]),
      Form("Nominal: %d iterations (dashed line)",niterate)}, .5, .85, 16, 700);
  c->SaveAs(pdfPath.c_str());
  fout->cd();
  gChi2->Write("gNonClosureChi2");

  c->SaveAs(Form("%s]", pdfPath.c_str()));
  fout->Close();
  cout << "Done. Wrote " << pdfPath << " and " << rootPath << endl;
}
