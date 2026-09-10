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

// Complementary toy test to toy_resp_iterations.C, following the OTHER half of the ATLAS
// dijet-xJ paper's pseudo-experiment procedure (Phys. Lett. B 774 (2017) 379, Sec. 6):
// "Stochastic variations of the data are generated based on its statistical uncertainty
// and each variation is unfolded and projected into xJ. The statistical covariance of the
// set is taken as the statistical uncertainty."
//
// Here the response matrix is held fixed at nominal; instead, the measured (purity-
// corrected Data) spectrum is Poisson-fluctuated bin-by-bin, many times, unfolded through
// the same nominal response, and the spread of the resulting unfolded distributions is the
// Data statistics' own contribution to the uncertainty on the unfolded result - the
// counterpart to toy_resp_iterations.C's response-matrix-statistics contribution. The
// deliverable requested is the sum (not mean) of the per-bin RMS across the toy ensemble,
// evaluated at each iteration count and plotted vs niter.

const int ir = 2; // nominal jet radius index (R=0.4)
const int nPtBinsUsed = ana::nPtBinsUsed; // physics analysis only uses ana::ptBins[ana::firstUsedPtBin..]
const int nToys = 1000; // toy count for the primary (fixed-iteration) ensemble
const int nToysScan = 500; // toy count per iteration point for the sum-of-RMS-vs-niter scan
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

// Poisson-fluctuate every bin of a 1D histogram: each bin's nominal content is taken as
// the mean of a Poisson distribution and a new value is drawn from it - same simplification
// as toy_resp_iterations.C's poissonToyMatrix, applied here to the measured Data spectrum
// instead of the response matrix.
TH1D * poissonToyHist(TH1D * nominal, const char * name) {
  TH1D * toy = (TH1D*)nominal->Clone(name);
  for (int b = 0; b <= nominal->GetNbinsX()+1; b++) {
    double mean = nominal->GetBinContent(b);
    double val = (mean > 0) ? gRandom->PoissonD(mean) : 0;
    toy->SetBinContent(b, val);
    toy->SetBinError(b, sqrt(std::max(val,0.0)));
  }
  return toy;
}

void toy_data_iterations(string systag = "nominal") {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gRandom->SetSeed(12345);

  drawer d("pythia", systag);
  string pdfPath = Form("/home/samson72/sphnx/gammajet_unfold/pdfs/toy_data_iterations_%s.pdf", systag.c_str());

  // Response matrix: full, cross-section-weighted combination of Photon5/10/20 - same
  // construction as draw_purity_corrected.C / draw_iteration_halfclosure.C. Held fixed at
  // nominal throughout this macro (contrast with toy_resp_iterations.C).
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
  // Pull = (toy - nominal)/sqrt(nominal) for every measured-spectrum bin with enough
  // content to be meaningful, pooled across many toys. Should be a unit Gaussian if the
  // toy generator is correct.
  cout << "Generating pull-distribution validation sample..." << endl;
  TH1D * hPull = new TH1D("hPull", ";(toy - nominal)/#sqrt{nominal};Measured-spectrum bins", 60, -5, 5);
  const int nPullToys = 200;
  for (int itoy = 0; itoy < nPullToys; itoy++) {
    TH1D * toyData = poissonToyHist(flatMeasured, Form("hToyPullCheck_%d", itoy));
    for (int b = 1; b <= flatMeasured->GetNbinsX(); b++) {
      double nominal = flatMeasured->GetBinContent(b);
      if (nominal < 5) continue; // low-stat bins have a skewed (non-Gaussian) Poisson - excluded from this Gaussian check only
      double toyval = toyData->GetBinContent(b);
      hPull->Fill((toyval - nominal)/sqrt(nominal));
    }
    delete toyData;
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
  lp->AddEntry(hPull, "Toy pulls (200 toys pooled)");
  lp->AddEntry(fUnit, "Unit Gaussian (expectation)");
  lp->AddEntry(fGaus, Form("Fit: #mu=%.2f, #sigma=%.2f", fGaus->GetParameter(1), fGaus->GetParameter(2)));
  lp->Draw();
  d.drawAll({"Data toy validation"},{Form("Jet R=%.1f",ana::JetRs[ir]), "Bins with nominal content > 5"}, .18, .85, 14, gPad->GetWh()*0.8);
  c->SaveAs(pdfPath.c_str());

  // ---- Primary ensemble: nToys toys at fixed niterPrimary ----
  cout << "Running primary toy ensemble: " << nToys << " toys at niter=" << niterPrimary << "..." << endl;
  int nFlatBins = flatMeasured->GetNbinsX();
  vector<double> sum(nFlatBins+2, 0), sumsq(nFlatBins+2, 0);
  vector<TH1D*> sampleToys; // a subset kept for the overlay plot
  int keepEvery = std::max(1, nToys/nSampleOverlay);
  for (int itoy = 0; itoy < nToys; itoy++) {
    TH1D * toyData = poissonToyHist(flatMeasured, Form("hToyData_%d", itoy));
    // includeSystematics=false: only GetBinContent() is read below, so skip the
    // (expensive - see unfold_utility.h's includeSystematics comment) response-matrix-
    // statistics covariance calculation whose result would just be discarded per toy.
    TH1D * hToy = unfold_utility::unfoldOnce(respRecoTemplate, respTruthTemplate, respMatrix2D, toyData, niterPrimary, Form("hUnfoldToy_%d", itoy), false);
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
    delete toyData;
    if (itoy % 200 == 0) cout << "  toy " << itoy << "/" << nToys << endl;
  }
  TH1D * hToyRMS = (TH1D*)flatMeasured->Clone("hToyRMS");
  for (int b = 0; b <= nFlatBins+1; b++) {
    double mean = sum[b]/nToys;
    double var = sumsq[b]/nToys - mean*mean;
    hToyRMS->SetBinContent(b, sqrt(std::max(var,0.0)));
  }
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
    lo->AddEntry(sampleDisp[0], Form("%d of %d Data toys", (int)sampleDisp.size(), nToys), "l");
    lo->Draw();
    d.drawAll({"p+p Run24 Data"},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
        Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV", ana::JetRs[ir], ana::jet_calib_pt_cut[ir])}, .5, .85, 14, gPad->GetWh()*0.8);
    c->SaveAs(pdfPath.c_str());
    delete hNomPt; delete hNomDisp;
    for (auto h : sampleDisp) delete h;
  }
  // sampleToys itself (the primary ensemble's kept-aside toy results, read from via
  // unflattenXj() above but never owned by anything past this point) is never freed
  // otherwise - a real leak of up to nSampleOverlay histograms per run.
  for (auto h : sampleToys) delete h;

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
    hNomDisp->SetFillColorAlpha(kGreen+1, 0.35);
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

  // Metric: chi2/NDF = mean over used bins of (toy RMS / nominal)^2 - the SAME units as
  // draw_iteration_halfclosure.C's pairwise "iter n vs iter n-1" convergence metric, and
  // the SAME units as toy_resp_iterations.C's toy metric (mean RMS/mean before was a plain
  // fraction; the raw sum-of-RMS before was in absolute count units - neither was
  // comparable to the other test or to the pairwise metric; this is).
  cout << "Running niter-dependence scan (" << nToysScan << " toys per iteration point)..." << endl;
  vector<double> chi2ToyByIter;
  vector<TH1D*> nomByIter; // nominal (non-toy) unfolded result at each scanned iteration - kept for the pairwise n-vs-(n-1) comparison below
  for (unsigned ii = 0; ii < iterationsToScan.size(); ii++) {
    int iter = iterationsToScan[ii];
    nomByIter.push_back(unfold_utility::unfoldOnce(respRecoTemplate, respTruthTemplate, respMatrix2D, flatMeasured, iter, Form("hNomIter_%d", iter)));
    vector<double> s(nFlatBins+2, 0), s2(nFlatBins+2, 0);
    for (int itoy = 0; itoy < nToysScan; itoy++) {
      TH1D * toyData = poissonToyHist(flatMeasured, Form("hToyDataScan_%d_%d", iter, itoy));
      // includeSystematics=false - see the primary toy loop's comment above.
      TH1D * hToy = unfold_utility::unfoldOnce(respRecoTemplate, respTruthTemplate, respMatrix2D, toyData, iter, Form("hUnfoldToyScan_%d_%d", iter, itoy), false);
      for (int b = 0; b <= nFlatBins+1; b++) {
        double v = hToy->GetBinContent(b);
        s[b] += v; s2[b] += v*v;
      }
      delete hToy;
      delete toyData;
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
    cout << "  niter=" << iter << ": Data-toy chi2/NDF = " << chi2ToyByIter.back() << endl;
  }

  // Pairwise "real" convergence metric, computed exactly as in draw_iteration_halfclosure.C
  // (and identically in toy_resp_iterations.C, since both use the same nominal response +
  // nominal Data - this curve should come out the same in both macros' plots, providing a
  // cross-check between them): chi2/NDF of (iter n - iter n-1)/(iter n-1), from the single
  // nominal (non-toy) unfolding chain.
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
  // a cheap extra precaution. toy_resp_iterations.C writes/reads the analogous
  // ".toy_resp_chi2_data.root" file - whichever of the two macros runs second picks up
  // both and produces the full four-curve comparison; running only one still produces a
  // valid (partial) plot.
  string respChi2DataFile = Form("/home/samson72/sphnx/gammajet_unfold/pdfs/.toy_resp_chi2_data_%s.root", systag.c_str());
  string dataChi2DataFile = Form("/home/samson72/sphnx/gammajet_unfold/pdfs/.toy_data_chi2_data_%s.root", systag.c_str());
  string chi2PdfPath      = Form("/home/samson72/sphnx/gammajet_unfold/pdfs/toy_iterations_chi2_%s.pdf", systag.c_str());
  TFile * fOut = TFile::Open(dataChi2DataFile.c_str(), "RECREATE");
  gToyChi2->Write("gToy");
  gPairChi2->Write("gPair");
  fOut->Close();
  gSystem->Exec(Form(
    "cd /home/samson72/sphnx/gammajet_unfold/drawing && root -b -l -q 'plot_toy_chi2_combined.C(\"%s\",\"%s\",\"%s\",\"Jet R=%.1f\")'",
    respChi2DataFile.c_str(), dataChi2DataFile.c_str(), chi2PdfPath.c_str(), ana::JetRs[ir]));

  cout << "Done. Wrote " << pdfPath << " and " << chi2PdfPath << endl;
}
