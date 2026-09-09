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

// Toy/pseudo-experiment test of the response-matrix statistical uncertainty, following
// the ATLAS dijet-xJ paper (Phys. Lett. B 774 (2017) 379, Sec. 6): "An additional
// covariance is obtained from applying the pseudo-experiment procedure to the response
// matrix." Concretely: take the nominal (pT,xJ) response matrix, Poisson-fluctuate every
// bin many times (each toy assumes the nominal bin content is the mean of a Poisson
// distribution), rebuild a RooUnfoldResponse from each fluctuated matrix, unfold the same
// (fixed) purity-corrected Data spectrum through it, and take the spread of the resulting
// unfolded distributions across toys as the response-matrix's contribution to the
// statistical uncertainty on the unfolded result - independent of and complementary to
// RooUnfold's analytic covariance (the "mean fractional uncertainty" used in
// draw_iteration_halfclosure.C), to toy_data_iterations.C (which instead toys the measured
// Data spectrum and holds the response fixed), and to the JES/JER response-matrix
// systematic variations (a deterministic shift/smear, not a statistical resampling).
//
// Only the response matrix is toyed here; the measured (purity-corrected Data) spectrum
// is held fixed, isolating the response-matrix-statistics component specifically.

const int ir = 2; // nominal jet radius index (R=0.4)
const int nPtBinsUsed = ana::nPtBinsUsed; // physics analysis only uses ana::ptBins[ana::firstUsedPtBin..]
const int nToys = 1000; // toy count for the primary (fixed-iteration) ensemble
const int nToysScan = 200; // toy count per iteration point for the niter-dependence scan
const int niterPrimary = 2; // fixed iteration count for the main toy ensemble (see draw_iteration_halfclosure.C's best-iteration scan)
const int nSampleOverlay = 50; // number of individual toy curves drawn on the overlay pages
const vector<int> iterationsToScan = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}; // consecutive, so "iter n vs iter n-1" is well-defined at every point

// The last 3 xJ bins in each pT bin have very low counts, so any chi2/NDF computed
// against iteration count is dominated by their noise rather than genuine convergence
// behavior - excluded from chi2ToyByIter/pairChi2ByIter below (still drawn everywhere else).
const int nXjBinsForChi2 = ana::nUnfoldXjBins - 3;

// densityForDisplay now lives in unfold_utility - see src/unfold_utility.h.
// buildFullyCorrected now lives in unfold_utility (purity-corrects all ana::nPtBins
// slices via unfold_utility::purityCorrect and reflattens for RooUnfold) - see
// src/unfold_utility.h.

// Poisson-fluctuate every bin of a response matrix: each bin's nominal content is taken
// as the mean of a Poisson distribution and a new value is drawn from it. This treats the
// (cross-section-weighted) bin content as if it were a raw count, the standard
// simplification used for "Poisson toy" response-matrix uncertainty studies (matching the
// language "assuming Poissonian statistics" - a literal per-entry weighted-Poisson
// treatment would require the underlying unweighted entry counts, which aren't carried
// through the histogram alone).
TH2D * poissonToyMatrix(TH2D * nominal, const char * name) {
  TH2D * toy = (TH2D*)nominal->Clone(name);
  for (int bx = 1; bx <= nominal->GetNbinsX(); bx++) {
    for (int by = 1; by <= nominal->GetNbinsY(); by++) {
      double mean = nominal->GetBinContent(bx, by);
      double val = (mean > 0) ? gRandom->PoissonD(mean) : 0;
      toy->SetBinContent(bx, by, val);
      toy->SetBinError(bx, by, sqrt(std::max(val,0.0)));
    }
  }
  return toy;
}

void toy_resp_iterations(string systag = "nominal") {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gRandom->SetSeed(12345);

  drawer d("pythia", systag);
  string pdfPath = Form("/home/samson72/sphnx/gammajet_unfold/pdfs/toy_resp_iterations_%s.pdf", systag.c_str());

  // Response matrix: full, cross-section-weighted combination of Photon5/10/20 - same
  // construction as draw_purity_corrected.C / draw_iteration_halfclosure.C.
  TH1D * respRecoTemplate  = d.get(Form("hrecoxj%i",ir), 1);
  TH1D * respTruthTemplate = d.get(Form("htruthxj%i",ir), 1);
  TH2D * respMatrix2D      = d.get2d(Form("hxjresponse%i",ir), 1);

  TH1D * flatA = d.get(Form("hrecoxj%i_0",ir), 0);
  TH1D * flatC = d.get(Form("hrecoxj%i_2",ir), 0);
  TH1D * flatMeasured = unfold_utility::buildFullyCorrected(flatA, flatC, "data", systag);
  TH1D * flatTruth = respTruthTemplate;

  TH1D * hNominal = unfold_utility::unfoldOnce(respRecoTemplate, respTruthTemplate, respMatrix2D, flatMeasured, niterPrimary, "hUnfoldNominal");

  TCanvas * c = new TCanvas("c","",700,900);
  c->SaveAs(Form("%s[", pdfPath.c_str()));

  // ---- Page 1: validation - are the Poisson toys behaving like Poisson toys? ----
  // Pull = (toy - nominal)/sqrt(nominal) for every response-matrix bin with enough
  // content to be meaningful, pooled across a handful of toys. If the toy generator is
  // correct this should be a unit Gaussian centred at 0; a bias or wrong width here would
  // mean the whole test downstream is not to be trusted.
  cout << "Generating pull-distribution validation sample..." << endl;
  TH1D * hPull = new TH1D("hPull", ";(toy - nominal)/#sqrt{nominal};Response-matrix bins", 60, -5, 5);
  const int nPullToys = 100;
  for (int itoy = 0; itoy < nPullToys; itoy++) {
    TH2D * toyMat = poissonToyMatrix(respMatrix2D, Form("hToyPullCheck_%d", itoy));
    for (int bx = 1; bx <= respMatrix2D->GetNbinsX(); bx++) {
      for (int by = 1; by <= respMatrix2D->GetNbinsY(); by++) {
        double nominal = respMatrix2D->GetBinContent(bx, by);
        if (nominal < 5) continue; // low-stat bins have a skewed (non-Gaussian) Poisson - excluded from this Gaussian check only, not from the toy unfolding itself
        double toyval = toyMat->GetBinContent(bx, by);
        hPull->Fill((toyval - nominal)/sqrt(nominal));
      }
    }
    delete toyMat;
  }
  TF1 * fGaus = new TF1("fGaus", "gaus", -5, 5);
  hPull->Scale(1./hPull->Integral(), "width");
  hPull->Fit(fGaus, "Q0");
  c->Clear();
  c->cd();
  gPad->SetTicks(1,1);
  gPad->SetLeftMargin(.15);
  hPull->SetLineColor(kBlack);
  hPull->SetMarkerColor(kBlack);
  hPull->SetMarkerStyle(20);
  hPull->GetYaxis()->SetTitle("Probability density");
  hPull->Draw("p e");
  TF1 * fUnit = new TF1("fUnit", "TMath::Gaus(x,0,1,1)", -5, 5);
  fUnit->SetLineColor(kBlue);
  fUnit->SetLineStyle(2);
  fUnit->Draw("same");
  fGaus->SetLineColor(kRed);
  fGaus->Draw("same");
  TLegend * lp = new TLegend(.55,.65,.88,.85);
  lp->SetLineWidth(0);
  lp->SetTextSize(0.03);
  lp->AddEntry(hPull, "Toy pulls (100 toys pooled)");
  lp->AddEntry(fUnit, "Unit Gaussian (expectation)");
  lp->AddEntry(fGaus, Form("Fit: #mu=%.2f, #sigma=%.2f", fGaus->GetParameter(1), fGaus->GetParameter(2)));
  lp->Draw();
  d.drawAll({"Response matrix toy validation"},{Form("Jet R=%.1f",ana::JetRs[ir]), "Bins with nominal content > 5"}, .18, .85, 14, gPad->GetWh()*0.8);
  c->SaveAs(pdfPath.c_str());

  // ---- Primary ensemble: nToys toys at fixed niterPrimary ----
  cout << "Running primary toy ensemble: " << nToys << " toys at niter=" << niterPrimary << "..." << endl;
  // Running sums per flattened bin for mean/RMS (Welford not needed - single pass sum/sumsq is fine here).
  int nFlatBins = flatMeasured->GetNbinsX();
  vector<double> sum(nFlatBins+2, 0), sumsq(nFlatBins+2, 0);
  vector<TH1D*> sampleToys; // a subset kept for the overlay plot
  int keepEvery = std::max(1, nToys/nSampleOverlay);
  for (int itoy = 0; itoy < nToys; itoy++) {
    TH2D * toyMat = poissonToyMatrix(respMatrix2D, Form("hToyMat_%d", itoy));
    // includeSystematics=false: only GetBinContent() is read below, so skip the
    // (expensive - see unfold_utility.h's includeSystematics comment) response-matrix-
    // statistics covariance calculation whose result would just be discarded per toy.
    TH1D * hToy = unfold_utility::unfoldOnce(respRecoTemplate, respTruthTemplate, toyMat, flatMeasured, niterPrimary, Form("hUnfoldToy_%d", itoy), false);
    for (int b = 0; b <= nFlatBins+1; b++) {
      double v = hToy->GetBinContent(b);
      sum[b] += v;
      sumsq[b] += v*v;
    }
    if (itoy % keepEvery == 0 && (int)sampleToys.size() < nSampleOverlay) {
      sampleToys.push_back(hToy);
    } else {
      delete hToy;
    }
    delete toyMat;
    if (itoy % 200 == 0) cout << "  toy " << itoy << "/" << nToys << endl;
  }
  TH1D * hToyMean = (TH1D*)flatMeasured->Clone("hToyMean");
  TH1D * hToyRMS  = (TH1D*)flatMeasured->Clone("hToyRMS");
  for (int b = 0; b <= nFlatBins+1; b++) {
    double mean = sum[b]/nToys;
    double var = sumsq[b]/nToys - mean*mean;
    hToyMean->SetBinContent(b, mean);
    hToyRMS->SetBinContent(b, sqrt(std::max(var,0.0)));
  }
  // Nominal unfolded result, but with bin errors replaced by the toy-ensemble RMS - this
  // is the actual deliverable: the response-matrix-statistics contribution to the
  // uncertainty on the unfolded xJ distribution.
  TH1D * hNominalToyErr = (TH1D*)hNominal->Clone("hNominalToyErr");
  for (int b = 0; b <= nFlatBins+1; b++) hNominalToyErr->SetBinError(b, hToyRMS->GetBinContent(b));

  // ---- Pages 2..(2+nPtBinsUsed-1): sample toy curves overlaid per pT bin ----
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hNomPt = unfold_utility::unflattenXj(hNominal, ipt, Form("hNomPt_%d", ipt));
    TH1D * hNomDisp = unfold_utility::densityForDisplay(hNomPt, Form("hNomDisp_%d", ipt));
    double ymax = hNomDisp->GetMaximum();
    vector<TH1D*> sampleDisp;
    for (unsigned k = 0; k < sampleToys.size(); k++) {
      TH1D * hp = unfold_utility::unflattenXj(sampleToys[k], ipt, Form("hSamplePt_%d_%d", ipt, k));
      TH1D * hd = unfold_utility::densityForDisplay(hp, Form("hSampleDisp_%d_%d", ipt, k));
      sampleDisp.push_back(hd);
      ymax = std::max(ymax, hd->GetMaximum());
      delete hp;
    }
    c->Clear();
    c->cd();
    gPad->SetTicks(1,1);
    gPad->SetLeftMargin(.15);
    hNomDisp->SetLineColor(kRed);
    hNomDisp->SetLineWidth(3);
    hNomDisp->GetYaxis()->SetRangeUser(0, ymax*1.3);
    hNomDisp->GetXaxis()->SetTitle("x_{J#gamma}");
    hNomDisp->Draw("hist");
    for (unsigned k = 0; k < sampleDisp.size(); k++) {
      sampleDisp[k]->SetLineColorAlpha(kGray+2, 0.35);
      sampleDisp[k]->SetLineWidth(1);
      sampleDisp[k]->Draw("hist same");
    }
    hNomDisp->Draw("hist same");
    TLegend * lo = new TLegend(.5,.6,.85,.75);
    lo->SetLineWidth(0);
    lo->SetTextSize(0.03);
    lo->AddEntry(hNomDisp, "Nominal unfolded", "l");
    lo->AddEntry(sampleDisp[0], Form("%d of %d response-matrix toys", (int)sampleDisp.size(), nToys), "l");
    lo->Draw();
    d.drawAll({"p+p Run24 Data"},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
        Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV", ana::JetRs[ir], ana::jet_calib_pt_cut[ir])}, .5, .85, 14, gPad->GetWh()*0.8);
    c->SaveAs(pdfPath.c_str());
    delete hNomPt; delete hNomDisp;
    for (auto h : sampleDisp) delete h;
  }

  // ---- Pages (2+nPtBinsUsed)..: nominal result with toy-derived error band, per pT bin ----
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hNomPt = unfold_utility::unflattenXj(hNominalToyErr, ipt, Form("hNomErrPt_%d", ipt));
    TH1D * hNomDisp = unfold_utility::densityForDisplay(hNomPt, Form("hNomErrDisp_%d", ipt));
    c->Clear();
    c->cd();
    gPad->SetTicks(1,1);
    gPad->SetLeftMargin(.15);
    hNomDisp->SetLineColor(kBlack);
    hNomDisp->SetMarkerColor(kBlack);
    hNomDisp->SetMarkerStyle(20);
    hNomDisp->SetFillColorAlpha(kAzure+7, 0.35);
    hNomDisp->GetYaxis()->SetRangeUser(0, hNomDisp->GetMaximum()*1.4);
    hNomDisp->GetXaxis()->SetTitle("x_{J#gamma}");
    hNomDisp->Draw("e2");
    hNomDisp->Draw("p e same");
    TLegend * le = new TLegend(.5,.62,.85,.75);
    le->SetLineWidth(0);
    le->SetTextSize(0.03);
    le->AddEntry(hNomDisp, "Nominal unfolded #pm toy RMS", "lf");
    le->Draw();
    d.drawAll({"p+p Run24 Data"},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
        Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV, %d toys", ana::JetRs[ir], ana::jet_calib_pt_cut[ir], nToys)}, .5, .85, 14, gPad->GetWh()*0.8);
    c->SaveAs(pdfPath.c_str());
    delete hNomPt; delete hNomDisp;
  }

  // ---- niter-dependence scan: mean fractional toy uncertainty vs iteration count ----
  // Mirrors the motivation behind ATLAS Fig. 4: more iterations amplifies sensitivity to
  // statistical noise already present in the response matrix, so this should increase
  // (or at least not decrease) with niter.
  // Metric: chi2/NDF = mean over used bins of (toy RMS / nominal)^2 - the SAME units as
  // draw_iteration_halfclosure.C's pairwise "iter n vs iter n-1" convergence metric
  // (mean of squared relative differences), so the two are directly comparable: is the
  // toy-implied statistical noise at a given iteration count bigger or smaller than the
  // actually observed change between consecutive iterations? Also directly comparable to
  // toy_data_iterations.C's identically-defined metric, since both report the same
  // dimensionless chi2/NDF quantity rather than one being a fraction and the other an
  // absolute count sum.
  cout << "Running niter-dependence scan (" << nToysScan << " toys per iteration point)..." << endl;
  vector<double> chi2ToyByIter;
  vector<TH1D*> nomByIter; // nominal (non-toy) unfolded result at each scanned iteration - kept for the pairwise n-vs-(n-1) comparison below
  for (unsigned ii = 0; ii < iterationsToScan.size(); ii++) {
    int iter = iterationsToScan[ii];
    nomByIter.push_back(unfold_utility::unfoldOnce(respRecoTemplate, respTruthTemplate, respMatrix2D, flatMeasured, iter, Form("hNomIter_%d", iter)));
    vector<double> s(nFlatBins+2, 0), s2(nFlatBins+2, 0);
    for (int itoy = 0; itoy < nToysScan; itoy++) {
      TH2D * toyMat = poissonToyMatrix(respMatrix2D, Form("hToyScan_%d_%d", iter, itoy));
      // includeSystematics=false - see the primary toy loop's comment above.
      TH1D * hToy = unfold_utility::unfoldOnce(respRecoTemplate, respTruthTemplate, toyMat, flatMeasured, iter, Form("hUnfoldToyScan_%d_%d", iter, itoy), false);
      for (int b = 0; b <= nFlatBins+1; b++) {
        double v = hToy->GetBinContent(b);
        s[b] += v; s2[b] += v*v;
      }
      delete hToy;
      delete toyMat;
    }
    double chi2Sum = 0;
    int nCounted = 0;
    for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
      for (int ixj = 0; ixj < nXjBinsForChi2; ixj++) {
        int flatbin = ipt*(ana::nUnfoldXjBins+2) + ixj + 1 + 1;
        double mean = s[flatbin]/nToysScan;
        double var  = s2[flatbin]/nToysScan - mean*mean;
        if (mean > 0) { chi2Sum += var/(mean*mean); nCounted++; }
      }
    }
    chi2ToyByIter.push_back(nCounted > 0 ? chi2Sum/nCounted : 0);
    cout << "  niter=" << iter << ": response-matrix-toy chi2/NDF = " << chi2ToyByIter.back() << endl;
  }

  // Pairwise "real" convergence metric, computed exactly as in draw_iteration_halfclosure.C:
  // chi2/NDF of (iter n - iter n-1)/(iter n-1), from the single nominal (non-toy) unfolding
  // chain - no toys involved. Plotted alongside the toy curve above so it's possible to see
  // directly whether the observed change from one iteration to the next is bigger or
  // smaller than the noise implied by the response matrix's own finite MC statistics.
  vector<int> pairIters;
  vector<double> pairChi2ByIter;
  for (unsigned k = 0; k < iterationsToScan.size(); k++) {
    double chi2 = 0;
    int ndf = 0;
    for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
      TH1D * uPrev = (k == 0 ?  unfold_utility::unflattenXj(flatMeasured  , ipt, Form("hPairPrev_%d_%d", k, ipt)) :
                                unfold_utility::unflattenXj(nomByIter[k-1], ipt, Form("hPairPrev_%d_%d", k, ipt)) );
      TH1D * uNext = unfold_utility::unflattenXj(nomByIter[k],   ipt, Form("hPairNext_%d_%d", k, ipt));
      for (int b = 1; b <= nXjBinsForChi2; b++) {
        double up = uPrev->GetBinContent(b);
        double un = uNext->GetBinContent(b);
        if (up <= 0) continue;
        double relDiff = (un-up)/up;
        chi2 += relDiff*relDiff;
        ndf++;
      }
      delete uPrev; delete uNext;
    }
    pairIters.push_back(iterationsToScan[k]);
    pairChi2ByIter.push_back(ndf > 0 ? chi2/ndf : 0);
    cout << "  niter=" << iterationsToScan[k] << ": observed chi2/NDF (vs iter " << (k == 0 ? 0 : iterationsToScan[k-1]) << ") = " << pairChi2ByIter.back() << endl;
  }
  for (auto h : nomByIter) delete h;

  TGraph * gToyChi2 = new TGraph(iterationsToScan.size());
  for (unsigned k = 0; k < iterationsToScan.size(); k++) gToyChi2->SetPoint(k, iterationsToScan[k], chi2ToyByIter[k]);
  TGraph * gPairChi2 = new TGraph(pairIters.size());
  for (unsigned k = 0; k < pairIters.size(); k++) gPairChi2->SetPoint(k, pairIters[k], pairChi2ByIter[k]);

  // Close the multi-page PDF here, without this final comparison page - see below.
  c->SaveAs(Form("%s]", pdfPath.c_str()));

  // Save this macro's toy graphs to a small file, then render the combined chi2/NDF
  // comparison page (response-matrix toy + Data toy + their quadrature sum + the real
  // observed convergence curve, all on one PDF - see plot_toy_chi2_combined.C) in a
  // SEPARATE, freshly-started ROOT process: TGraph's log-scale point/line painting was
  // found to silently drop points once they land close enough to the frame's top edge -
  // reproducible regardless of draw option, frame construction, or canvas reuse, and
  // mitigated by generous y-axis headroom in the helper; running it as a fresh process is
  // a cheap extra precaution. toy_data_iterations.C writes/reads the analogous
  // ".toy_data_chi2_data.root" file - whichever of the two macros runs second picks up
  // both and produces the full four-curve comparison; running only one still produces a
  // valid (partial) plot.
  string respChi2DataFile = Form("/home/samson72/sphnx/gammajet_unfold/pdfs/.toy_resp_chi2_data_%s.root", systag.c_str());
  string dataChi2DataFile = Form("/home/samson72/sphnx/gammajet_unfold/pdfs/.toy_data_chi2_data_%s.root", systag.c_str());
  string chi2PdfPath      = Form("/home/samson72/sphnx/gammajet_unfold/pdfs/toy_iterations_chi2_%s.pdf", systag.c_str());
  TFile * fOut = TFile::Open(respChi2DataFile.c_str(), "RECREATE");
  gToyChi2->Write("gToy");
  gPairChi2->Write("gPair");
  fOut->Close();
  gSystem->Exec(Form(
    "cd /home/samson72/sphnx/gammajet_unfold/drawing && root -b -l -q 'plot_toy_chi2_combined.C(\"%s\",\"%s\",\"%s\",\"Jet R=%.1f\")'",
    respChi2DataFile.c_str(), dataChi2DataFile.c_str(), chi2PdfPath.c_str(), ana::JetRs[ir]));

  cout << "Done. Wrote " << pdfPath << " and " << chi2PdfPath << endl;
}
