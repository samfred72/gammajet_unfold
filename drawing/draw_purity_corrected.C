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

// Physics-level analysis uses only ana::ptBins[ana::firstUsedPtBin..] (15-20,20-25,25-35
// GeV) - ana::ptBins also carries a low-pT migration-only buffer bin (13-15 GeV) below that
// and a high-pT overflow bin (35-100 GeV) above it, neither of which is reported. The
// response matrix and the flattened (pT,xJ) histograms span all ana::nPtBins bins though -
// see buildFullyCorrected() below, which corrects all of them to keep the unfold input
// dimensionally consistent with the response matrix, and only the ana::nPtBinsUsed used
// ones get displayed.
const int nPtBinsUsed = ana::nPtBinsUsed;
const int ir = 2; // nominal jet radius index (R=0.4); change here if a different radius is wanted.
const int niterate = 2; // best-iteration scan result (see draw_iteration_halfclosure.C / toy_iterations_chi2.pdf)

// The last 3 xJ bins in each pT bin have very low counts, so any chi2/NDF computed
// against iteration count is dominated by their noise rather than genuine convergence
// behavior - excluded from computeChi2NDF below, per-bin only (still drawn everywhere else).
const int nXjBinsForChi2 = ana::nUnfoldXjBins - 3;

// signal(xJ) = A(xJ) - (1-P)*(N_A/N_C)*C(xJ), with N_A, N_C the pT-bin-integrated
// region A/C totals (from the same xJ-differential histograms, so no cross-binning
// mismatch), and P = ana::getPurity(ptBins[ipt],ptBins[ipt+1]) - the actual
// puritymaker.C bootstrap point for this bin, not a smooth fit. pErrLow/pErrHigh are
// puritymaker.C's asymmetric (16th/84th percentile) bootstrap errors for that same
// point, propagated through d(signal)/dP = (N_A/N_C)*C; P (and its errors) is one
// shared value for the whole pT bin, so this term is fully correlated across xJ bins,
// not independent per bin - same simplification the rest of this codebase uses for
// such correction uncertainties.
//
// TH1D bins can only hold one symmetric error, so h's bin errors use the larger of
// errLow/errHigh (a conservative choice) - that's what feeds RooUnfold, Integral(),
// Chi2 calcs, etc. downstream, none of which support asymmetric errors anyway. If
// graphOut is non-null, *graphOut receives a TGraphAsymmErrors with the true
// asymmetric errors, for display where the asymmetry should actually be visible.
// Returns h=nullptr (and *graphOut=nullptr if requested) if region C has no statistics
// in this pT bin (N_A/N_C undefined) - callers must check.
TH1D * purityCorrectP(TH1D * A, TH1D * C, float p, float pErrLow, float pErrHigh, const char * name, TGraphAsymmErrors ** graphOut = nullptr) {
  float NA = A->Integral();
  float NC = C->Integral();
  if (NC <= 0) {
    cout << "WARNING: " << name << " has zero region-C statistics - cannot cross-normalize, skipping." << endl;
    if (graphOut) *graphOut = nullptr;
    return nullptr;
  }
  float scale = (1-p)*(NA/NC);
  TH1D * h = (TH1D*)A->Clone(name);
  TGraphAsymmErrors * g = graphOut ? new TGraphAsymmErrors(A->GetNbinsX()) : nullptr;
  if (g) g->SetName(Form("%s_graph", name));
  for (int i = 1; i <= A->GetNbinsX(); i++) {
    float a  = A->GetBinContent(i);
    float ae = A->GetBinError(i);
    float c  = C->GetBinContent(i);
    float ce = C->GetBinError(i);
    float dPurityLow  = (NA/NC)*c*pErrLow;
    float dPurityHigh = (NA/NC)*c*pErrHigh;
    float content = a - scale*c;
    float errLow  = sqrt(ae*ae + pow(scale*ce,2) + pow(dPurityLow,2));
    float errHigh = sqrt(ae*ae + pow(scale*ce,2) + pow(dPurityHigh,2));
    h->SetBinContent(i, content);
    h->SetBinError(i, std::max(errLow, errHigh));
    if (g) {
      double xc  = A->GetXaxis()->GetBinCenter(i);
      double xlo = xc - A->GetXaxis()->GetBinLowEdge(i);
      double xhi = A->GetXaxis()->GetBinUpEdge(i) - xc;
      g->SetPoint(i-1, xc, content);
      g->SetPointError(i-1, xlo, xhi, errLow, errHigh);
    }
  }
  if (graphOut) *graphOut = g;
  return h;
}

// Graph analogue of densityForDisplay() below: divide y-values and asymmetric y-errors
// by bin width (from refBinning's axis) for display, since ana::unfoldXjBins is
// non-uniform. X position/width are left as-is (already set from real bin edges).
TGraphAsymmErrors * densityForDisplayGraph(TGraphAsymmErrors * g, TH1D * refBinning, const char * name) {
  TGraphAsymmErrors * gd = (TGraphAsymmErrors*)g->Clone(name);
  for (int i = 0; i < gd->GetN(); i++) {
    double x, y;
    gd->GetPoint(i, x, y);
    double w = refBinning->GetXaxis()->GetBinWidth(i+1);
    gd->SetPoint(i, x, y/w);
    gd->SetPointEYlow(i, gd->GetErrorYlow(i)/w);
    gd->SetPointEYhigh(i, gd->GetErrorYhigh(i)/w);
  }
  return gd;
}

// (1-P)*(N_A/N_C)*C(x): the piece subtracted off region A, drawn separately.
TH1D * bkgFromRegionC(TH1D * C, float p, float NA, float NC, const char * name) {
  float scale = (1-p)*(NA/NC);
  TH1D * h = (TH1D*)C->Clone(name);
  h->Scale(scale);
  return h;
}

// ana::unfoldXjBins is non-uniform (0.1-wide up to xJ=1.3, then 0.2,0.2,0.3) - plotting
// raw bin content directly makes an otherwise-smooth density look like it has a
// shelf/cliff right where the bin width changes. Divide by bin width for display only;
// the histograms used for Integral()/N_A/N_C and the subtraction itself stay as raw counts.
TH1D * densityForDisplay(TH1D * h, const char * name) {
  TH1D * hd = (TH1D*)h->Clone(name);
  hd->Scale(1., "width");
  hd->GetYaxis()->SetTitle("Counts / bin width");
  return hd;
}

// Purity-correct all ana::nPtBins slices (not just the 3 used for physics results)
// of flatA/flatC and write the result into one full flattened histogram matching the
// response matrix's dimensionality, so RooUnfold sees a complete, consistently-binned
// "measured" vector - unfolding a 2D (pT,xJ) measurement needs the whole flattened
// vector at once because migration crosses pT-bin boundaries, not just xJ ones.
// If a pT slice's region C is empty (can't purity-correct), falls back to raw region A
// for that slice only - it isn't part of the displayed/physics-used bins anyway.
TH1D * buildFullyCorrected(TH1D * flatA, TH1D * flatC, const char * tag, string systag) {
  TH1D * flatCorrected = (TH1D*)flatA->Clone(Form("hxjcorrected_flat_%s", tag));
  flatCorrected->Reset("ICES");
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
    float ptlow  = ana::ptBins[ipt];
    float pthigh = ana::ptBins[ipt+1];
    float p        = ana::getPurity(ptlow, pthigh, systag);
    float pErrLow  = ana::getPurityErrorLow(ptlow, pthigh, systag);
    float pErrHigh = ana::getPurityErrorHigh(ptlow, pthigh, systag);
    TH1D * A = unfold_utility::unflattenXj(flatA, ipt, Form("htmpA_%s_pt%d", tag, ipt));
    TH1D * C = unfold_utility::unflattenXj(flatC, ipt, Form("htmpC_%s_pt%d", tag, ipt));
    TH1D * hcorr = purityCorrectP(A, C, p, pErrLow, pErrHigh, Form("htmpcorr_%s_pt%d", tag, ipt));
    unfold_utility::reflattenXj(hcorr ? hcorr : A, ipt, flatCorrected);
    delete A; delete C; if (hcorr) delete hcorr;
  }
  return flatCorrected;
}

// Bayesian-unfolding iteration scan: candidate iteration counts to test the chosen
// `niterate` above against.
const vector<int> iterationsToTest = {1,2,3,4,5,6,8,10,12,15};

double computeChi2NDF(TH1D * h1, TH1D * h2) {
  double chi2 = 0;
  int ndf = 0;
  for (int b = 1; b <= nXjBinsForChi2; b++) {
    double v1 = h1->GetBinContent(b), e1 = h1->GetBinError(b);
    double v2 = h2->GetBinContent(b), e2 = h2->GetBinError(b);
    double err2 = e1*e1 + e2*e2;
    if (err2 <= 0) continue;
    chi2 += pow(v1-v2,2)/err2;
    ndf++;
  }
  return (ndf > 0) ? chi2/ndf : 0;
}

// Scans iterationsToTest, scores each by chi2/ndf against the (fixed) photon+jet truth
// summed over the used pT bins, and picks the iteration count with the lowest combined
// score. Only that winning iteration's unfolded result gets drawn (not every candidate) -
// alongside one summary scan page showing where it converged, so the choice is visible.
// Written to its own pdf (iterationsPdfPath), separate from purity_corrected.pdf.
void draw_iteration_test(TH1D * flatCorrected, TH1D * flatTruth, RooUnfoldResponse * response,
    const char * label, const char * tag, TCanvas * c, TFile * fout,
    drawer & d, const char * iterationsPdfPath) {
  vector<TH1D*> truthDisp(ana::nPtBins);
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hTruth = unfold_utility::unflattenXj(flatTruth, ipt, Form("hxjtruth_iter_%s_pt%d", tag, ipt));
    truthDisp[ipt] = densityForDisplay(hTruth, Form("hxjtruth_iter_%s_pt%d_disp", tag, ipt));
    truthDisp[ipt]->Scale(1./truthDisp[ipt]->Integral());
  }

  TGraph * gcombined = new TGraph();
  gcombined->SetName(Form("giterscan_%s", tag));
  int bestIter = iterationsToTest[0];
  double bestScore = 1e18;
  map<int, vector<TH1D*>> unfoldedByIter; // iteration -> per-pT-bin shape-normalized unfolded

  for (unsigned k = 0; k < iterationsToTest.size(); k++) {
    int iter = iterationsToTest[k];
    TH1D * hUnfoldFull = unfold_utility::unfoldOnce(response, flatCorrected, iter, Form("hUnfoldIterScan_%s_%d", tag, iter));
    double combined = 0;
    vector<TH1D*> perPt(ana::nPtBins);
    for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
      TH1D * hUnfold = unfold_utility::unflattenXj(hUnfoldFull, ipt, Form("hxjunfold_iter%d_%s_pt%d", iter, tag, ipt));
      TH1D * hUnfolddisp = densityForDisplay(hUnfold, Form("hxjunfold_iter%d_%s_pt%d_disp", iter, tag, ipt));
      double integral = hUnfolddisp->Integral();
      if (integral > 0) hUnfolddisp->Scale(1./integral);
      combined += computeChi2NDF(hUnfolddisp, truthDisp[ipt]);
      perPt[ipt] = hUnfolddisp;
      delete hUnfold;
    }
    gcombined->SetPoint(k, iter, combined);
    if (combined < bestScore) {
      bestScore = combined;
      bestIter = iter;
    }
    unfoldedByIter[iter] = perPt;
  }

  // Page 1: chi2/ndf (summed over used pT bins) vs iteration count, marking the winner.
  c->Clear();
  c->cd();
  gPad->SetTicks();
  gPad->SetLeftMargin(.15);
  gcombined->SetMarkerStyle(20);
  gcombined->SetMarkerColor(kBlack);
  gcombined->SetLineColor(kBlack);
  gcombined->SetLineWidth(2);
  double ymax = 0;
  for (int k = 0; k < gcombined->GetN(); k++) {
    double x,y;
    gcombined->GetPoint(k,x,y);
    ymax = std::max(ymax,y);
  }
  gcombined->SetMinimum(0);
  gcombined->SetMaximum(ymax*1.3);
  gcombined->GetXaxis()->SetTitle("Bayesian unfolding iterations");
  gcombined->GetYaxis()->SetTitle("#Sigma #chi^{2}/NDF vs #gamma+jet truth (summed over p_{T} bins)");
  gcombined->Draw("APL");
  TLine * lbest = new TLine(bestIter, 0, bestIter, ymax*1.3);
  lbest->SetLineStyle(9);
  lbest->SetLineColor(kRed);
  lbest->Draw("same");
  d.drawAll({label},{Form("Jet R=%.1f", ana::JetRs[ir]),
      Form("Best: %d iterations (dashed line)", bestIter)}, .5, .85, 16, 700);
  c->SaveAs(iterationsPdfPath);
  fout->cd();
  gcombined->Write();

  // Pages 2..N: the winning iteration's unfolded result vs truth, per used pT bin.
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hUnfolddisp = unfoldedByIter[bestIter][ipt];
    c->Clear();
    c->cd();
    TPad * p1 = new TPad(Form("piter1_%s_%d",tag,ipt),"",0,.35,1,1);
    TPad * p2 = new TPad(Form("piter2_%s_%d",tag,ipt),"",0,0,1,.35);
    p1->Draw();
    p2->Draw();

    p1->cd();
    p1->SetBottomMargin(0.02);
    p1->SetLeftMargin(.15);
    gPad->SetTicks(1,1);
    hUnfolddisp->SetLineColor(kRed);
    hUnfolddisp->SetMarkerColor(kRed);
    hUnfolddisp->SetMarkerStyle(21);
    hUnfolddisp->SetLineWidth(2);
    hUnfolddisp->GetXaxis()->SetLabelSize(0);
    hUnfolddisp->GetXaxis()->SetTitle("");
    hUnfolddisp->GetYaxis()->SetTitle("Shape-normalized counts / bin width");
    hUnfolddisp->GetYaxis()->SetRangeUser(0, std::max(hUnfolddisp->GetMaximum(), truthDisp[ipt]->GetMaximum())*1.4);
    hUnfolddisp->Draw("p e");
    truthDisp[ipt]->SetLineColor(kBlack);
    truthDisp[ipt]->SetMarkerColor(kBlack);
    truthDisp[ipt]->SetMarkerStyle(20);
    truthDisp[ipt]->SetLineWidth(2);
    truthDisp[ipt]->Draw("p e same");
    TLegend * l = new TLegend(.55,.65,.85,.80);
    l->SetLineWidth(0);
    l->SetTextSize(0.035);
    l->AddEntry(hUnfolddisp,   Form("Unfolded (%d iter.)", bestIter));
    l->AddEntry(truthDisp[ipt],"Truth (#gamma+jet MC)");
    l->Draw();
    d.drawAll({label},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
        Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV", ana::JetRs[ir], ana::jet_calib_pt_cut[ir])}, .5, .85, 16, gPad->GetWh()*0.8);

    p2->cd();
    p2->SetTopMargin(0.02);
    p2->SetBottomMargin(0.3);
    p2->SetLeftMargin(.15);
    gPad->SetTicks(1,1);
    TH1D * hratio = (TH1D*)hUnfolddisp->Clone(Form("hxjratio_bestiter_%s_pt%d", tag, ipt));
    hratio->Divide(truthDisp[ipt]);
    hratio->SetLineColor(kBlack);
    hratio->SetMarkerColor(kBlack);
    hratio->SetMarkerStyle(20);
    hratio->GetYaxis()->SetRangeUser(0.5,1.5);
    hratio->GetYaxis()->SetTitle("Unfolded / Truth");
    hratio->GetYaxis()->SetTitleSize(0.09);
    hratio->GetYaxis()->SetLabelSize(0.08);
    hratio->GetXaxis()->SetTitle("x_{J#gamma}");
    hratio->GetXaxis()->SetTitleSize(0.09);
    hratio->GetXaxis()->SetLabelSize(0.08);
    hratio->Draw("p e");
    TLine * line = new TLine(ana::unfoldXjBins[0],1,ana::unfoldXjBins[ana::nUnfoldXjBins],1);
    line->SetLineStyle(9);
    line->Draw("same");
    c->SaveAs(iterationsPdfPath);
    fout->cd();
    hratio->Write();
    hUnfolddisp->Write();
  }
}

// flatTruth is always the photon+jet signal MC's truth-level flattened xJ (same sample
// the response matrix comes from - see draw_purity_corrected()), regardless of which
// sample's "measured" spectrum this call is unfolding. Data has no real truth level
// (unfolder.cc sets truth=reco as a placeholder for non-MC, which is why Data's own
// stored response matrix is trivially diagonal), and for a background-only sample like
// Jet MC, that sample's own "truth" would mean truth-level region-A-selected dijet
// events, not a meaningful physics reference either. The photon+jet truth is the one
// fixed, physically meaningful expectation to compare every unfolded result against.
void draw_one_sample(int type, const char * label, const char * tag, TCanvas * c, TFile * fout, RooUnfoldResponse * response, TH1D * flatTruth,
    drawer & d, const char * iterationsPdfPath, const char * purityPdfPath, string systag) {
  TH1D * flatA = d.get(Form("hrecoxj%i_0",ir), type);
  TH1D * flatC = d.get(Form("hrecoxj%i_2",ir), type);

  // Unfold the purity-corrected (not raw) measured spectrum - the whole point of
  // correcting before unfolding is that RooUnfold should never see the background.
  TH1D * flatCorrected = buildFullyCorrected(flatA, flatC, tag, systag);
  TH1D * flatUnfolded = unfold_utility::unfoldOnce(response, flatCorrected, niterate, Form("hUnfolded_flat_%s", tag));

  draw_iteration_test(flatCorrected, flatTruth, response, label, tag, c, fout, d, iterationsPdfPath);

  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    float ptlow  = ana::ptBins[ipt];
    float pthigh = ana::ptBins[ipt+1];
    float p        = ana::getPurity(ptlow, pthigh, systag);
    float pErrLow  = ana::getPurityErrorLow(ptlow, pthigh, systag);
    float pErrHigh = ana::getPurityErrorHigh(ptlow, pthigh, systag);

    TH1D * A = unfold_utility::unflattenXj(flatA, ipt, Form("hxjA_%s_pt%d", tag, ipt));
    TH1D * C = unfold_utility::unflattenXj(flatC, ipt, Form("hxjC_%s_pt%d", tag, ipt));

    float NA = A->Integral();
    float NC = C->Integral();
    TGraphAsymmErrors * hcorrGraph = nullptr;
    TH1D * hcorr = purityCorrectP(A, C, p, pErrLow, pErrHigh, Form("hxjcorrected_%s_pt%d", tag, ipt), &hcorrGraph);

    // Page 1: raw A, background subtracted from C (cross-normalized), and the corrected result.
    // Draw bin-width-normalized clones only - see densityForDisplay().
    c->Clear();
    c->cd();
    gPad->SetTicks();
    TH1D * Adisp = densityForDisplay(A, Form("hxjA_%s_pt%d_disp", tag, ipt));
    Adisp->SetLineColor(kBlack);
    Adisp->SetLineWidth(2);
    Adisp->GetXaxis()->SetTitle("x_{J#gamma}");
    Adisp->Draw("hist");
    if (hcorr) {
      TH1D * hbkg = bkgFromRegionC(C, p, NA, NC, Form("hxjbkg_%s_pt%d", tag, ipt));
      TH1D * hbkgdisp  = densityForDisplay(hbkg,  Form("hxjbkg_%s_pt%d_disp", tag, ipt));
      // hcorrdisp (TH1D, conservative symmetric error) is kept only for the Mean corr
      // text below - the actual plotted "purity-corrected signal" curve uses the
      // asymmetric-error graph so the true bootstrap errors are visible, not collapsed.
      TH1D * hcorrdisp = densityForDisplay(hcorr, Form("hxjcorrected_%s_pt%d_disp", tag, ipt));
      TGraphAsymmErrors * hcorrGraphDisp = densityForDisplayGraph(hcorrGraph, hcorr, Form("hxjcorrected_%s_pt%d_graphdisp", tag, ipt));
      hbkgdisp->SetLineColor(kAzure+2);
      hbkgdisp->SetLineWidth(2);
      hbkgdisp->Draw("hist same");
      hcorrGraphDisp->SetLineColor(kOrange+7);
      hcorrGraphDisp->SetMarkerColor(kOrange+7);
      hcorrGraphDisp->SetMarkerStyle(20);
      hcorrGraphDisp->SetLineWidth(2);
      hcorrGraphDisp->Draw("p same");
      TLegend * l = new TLegend(.55,.54,.85,.67);
      l->SetLineWidth(0);
      l->SetTextSize(0.028);
      l->AddEntry(Adisp,     "Region A (raw)");
      l->AddEntry(hbkgdisp,  Form("(1-P)#times#frac{N_{A}}{N_{C}}#times Region C"));
      l->AddEntry(hcorrGraphDisp, "Purity-corrected signal");
      l->Draw();
      fout->cd();
      hbkg->Write();
      hcorr->Write();
      hcorrGraph->Write();
      d.drawText(Form("Mean raw: %.2f #pm %.2f", Adisp->GetMean(), Adisp->GetMeanError()), 0.6,0.5);
      d.drawText(Form("Mean corr: %.2f #pm %.2f", hcorrdisp->GetMean(), hcorrdisp->GetMeanError()), 0.6,0.45);
    }
    else {
      d.drawText("Region C empty in this p_{T} bin - correction unavailable", .15, .5, kRed, 18);
    }
    // Info text block sits above the legend with a clear gap (drawAll's 5 lines here
    // span drawy down to roughly drawy-0.15).
    d.drawAll({label},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ptlow,pthigh),
        Form("P = %.3f +%.3f/-%.3f",p,pErrHigh,pErrLow),
        Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV", ana::JetRs[ir], ana::jet_calib_pt_cut[ir])}, .5, .85, 18, 700);
    c->SaveAs(purityPdfPath);

    fout->cd();
    A->Write(Form("hxjA_%s_pt%d", tag, ipt));
    C->Write(Form("hxjC_%s_pt%d", tag, ipt));

    // Page 2: purity-corrected (reco) vs unfolded vs truth, with a ratio panel below.
    // All three shape-normalized to unit area (Integral("width")) - reco/unfolded/truth
    // sit at different absolute scales due to reconstruction efficiency, so an absolute
    // comparison isn't meaningful; this matches the normalization convention already
    // used for this same comparison in drawing/draw_xj_unfold.C.
    if (hcorr) {
      TH1D * hUnfold = unfold_utility::unflattenXj(flatUnfolded, ipt, Form("hxjunfolded_%s_pt%d", tag, ipt));
      TH1D * hUnfolddisp = densityForDisplay(hUnfold, Form("hxjunfolded_%s_pt%d_disp", tag, ipt));
      hUnfolddisp->Scale(1./hUnfolddisp->Integral());
      TH1D * hCorrShape = (TH1D*)densityForDisplay(hcorr, Form("hxjcorrected_%s_pt%d_shape", tag, ipt));
      hCorrShape->Scale(1./hCorrShape->Integral());

      TH1D * hTruth = unfold_utility::unflattenXj(flatTruth, ipt, Form("hxjtruth_%s_pt%d", tag, ipt));
      TH1D * hTruthdisp = densityForDisplay(hTruth, Form("hxjtruth_%s_pt%d_disp", tag, ipt));
      hTruthdisp->Scale(1./hTruthdisp->Integral());
      fout->cd();
      hTruth->Write();
      hUnfold->Write();

      c->Clear();
      c->cd();
      TPad * p1 = new TPad(Form("p1_%s_%d",tag,ipt),"",0,.35,1,1);
      TPad * p2 = new TPad(Form("p2_%s_%d",tag,ipt),"",0,0,1,.35);
      p1->Draw();
      p2->Draw();

      p1->cd();
      p1->SetBottomMargin(0.02);
      p1->SetLeftMargin(.15);
      gPad->SetTicks(1,1);
      hCorrShape->SetLineColor(kAzure+2);
      hCorrShape->SetMarkerColor(kAzure+2);
      hCorrShape->SetMarkerStyle(24);
      hCorrShape->SetLineWidth(2);
      hCorrShape->GetXaxis()->SetLabelSize(0);
      hCorrShape->GetXaxis()->SetTitle("");
      hCorrShape->GetYaxis()->SetTitle("Shape-normalized counts / bin width");
      hCorrShape->GetYaxis()->SetRangeUser(0, std::max({hCorrShape->GetMaximum(),
            hUnfolddisp->GetMaximum(), hTruthdisp->GetMaximum()})*1.4);
      hCorrShape->Draw("p e");
      hUnfolddisp->SetLineColor(kRed);
      hUnfolddisp->SetMarkerColor(kRed);
      hUnfolddisp->SetMarkerStyle(21);
      hUnfolddisp->SetLineWidth(2);
      hUnfolddisp->Draw("p e same");
      hTruthdisp->SetLineColor(kBlack);
      hTruthdisp->SetMarkerColor(kBlack);
      hTruthdisp->SetMarkerStyle(20);
      hTruthdisp->SetLineWidth(2);
      hTruthdisp->Draw("p e same");
      TLegend * l2 = new TLegend(.55,.50,.85,.70);
      l2->SetLineWidth(0);
      l2->SetTextSize(0.035);
      l2->AddEntry(hCorrShape,  "Purity-corrected (reco)");
      l2->AddEntry(hUnfolddisp, "Unfolded (reco)");
      l2->AddEntry(hTruthdisp,  "Truth (#gamma+jet MC)");
      l2->Draw();
      d.drawAll({label},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ptlow,pthigh),
          Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV", ana::JetRs[ir], ana::jet_calib_pt_cut[ir])}, .5, .85, 16, gPad->GetWh()*0.8);

      p2->cd();
      p2->SetTopMargin(0.02);
      p2->SetBottomMargin(0.3);
      p2->SetLeftMargin(.15);
      gPad->SetTicks(1,1);
      TH1D * hratio = (TH1D*)hUnfolddisp->Clone(Form("hxjratio_%s_pt%d", tag, ipt));
      hratio->Divide(hTruthdisp);
      hratio->SetLineColor(kBlack);
      hratio->SetMarkerColor(kBlack);
      hratio->SetMarkerStyle(20);
      hratio->GetYaxis()->SetRangeUser(0.5,1.5);
      hratio->GetYaxis()->SetTitle("Unfolded / Pythia Truth");
      hratio->GetYaxis()->SetTitleSize(0.09);
      hratio->GetYaxis()->SetLabelSize(0.08);
      hratio->GetXaxis()->SetTitle("x_{J#gamma}");
      hratio->GetXaxis()->SetTitleSize(0.09);
      hratio->GetXaxis()->SetLabelSize(0.08);
      hratio->Draw("p e");
      TLine * line = new TLine(ana::unfoldXjBins[0],1,ana::unfoldXjBins[ana::nUnfoldXjBins],1);
      line->SetLineStyle(9);
      line->Draw("same");
      fout->cd();
      hratio->Write();
      c->SaveAs(purityPdfPath);
    }
  }
}

void draw_purity_corrected(string systag = "nominal") {
  gStyle->SetOptStat(0);

  drawer d("pythia", systag);
  string purityRootPath     = Form("/home/samson72/sphnx/gammajet_unfold/hists/purity_corrected_%s.root", systag.c_str());
  string purityPdfPath      = Form("/home/samson72/sphnx/gammajet_unfold/pdfs/purity_corrected_%s.pdf", systag.c_str());
  string iterationsPdfPath  = Form("/home/samson72/sphnx/gammajet_unfold/pdfs/purity_corrected_iterations_%s.pdf", systag.c_str());

  TFile * fout = TFile::Open(purityRootPath.c_str(),"RECREATE");
  TCanvas * c = new TCanvas("c","",700,700);
  gPad->SetLeftMargin(.15);
  c->SaveAs(Form("%s[", purityPdfPath.c_str()));
  c->SaveAs(Form("%s[", iterationsPdfPath.c_str()));

  // Response matrix: Data's own stored response (response_full_jetR2 etc.) is trivially
  // diagonal - unfolder.cc sets truth=reco as a placeholder for non-MC samples, so it
  // carries no real migration information. The physically meaningful response comes from
  // the signal (photon+jet) MC; type=1, isample=-1 gives the cross-section-weighted
  // combination of Photon5/10/20, matching the exact convention already established in
  // drawing/draw_xj_unfold.C. This one response is reused for every sample below, since
  // it characterizes the detector, independent of which "measured" spectrum is unfolded.
  TH1D * respRecoTemplate  = d.get(Form("hrecoxj%i",ir), 1);
  TH1D * respTruthTemplate = d.get(Form("htruthxj%i",ir), 1);
  TH2D * respMatrix2D      = d.get2d(Form("hxjresponse%i",ir), 1);
  RooUnfoldResponse * response = new RooUnfoldResponse(respRecoTemplate, respTruthTemplate, respMatrix2D);

  draw_one_sample(0, "p+p Run24 Data",     "data",     c, fout, response, respTruthTemplate, d, iterationsPdfPath.c_str(), purityPdfPath.c_str(), systag);
  //draw_one_sample(1, "Pythia8 #gamma+jet", "photonmc", c, fout, response, respTruthTemplate, d, iterationsPdfPath.c_str(), purityPdfPath.c_str(), systag);
  //draw_one_sample(2, "Pythia8 Jet",      "jetmc",    c, fout, response, respTruthTemplate, d, iterationsPdfPath.c_str(), purityPdfPath.c_str(), systag);

  c->SaveAs(Form("%s]", purityPdfPath.c_str()));
  c->SaveAs(Form("%s]", iterationsPdfPath.c_str()));
  fout->Close();
}
