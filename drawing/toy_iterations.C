#include "../src/ana.h"
#include "../src/drawer.h"
#include "../src/unfold_utility.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
// Explicit load: the sibling gammajet project's libgammajet.so has same-named classes.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Poisson-toy estimate of a statistical uncertainty on the unfolded x_J (ATLAS dijet x_J,
// PLB 774 (2017) 379, Sec. 6). source = "resp": toy the response matrix, Data fixed;
// source = "data": toy the purity-corrected Data spectrum, response fixed. Also scans the
// toy spread vs iteration count against the observed iteration-to-iteration change, and
// renders the combined chi2/NDF page (plot_toy_chi2_combined.C) once both sources exist.

const int ir = 2; // R = 0.4
const int nPtBinsUsed = ana::nPtBinsUsed;
const int nToys = 1000;
const int niterPrimary = 2;
const int nSampleOverlay = 50;
const vector<int> iterationsToScan = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
// Last 3 x_J bins per pT bin are too sparse for the chi2 (see draw_covariance_chi2.C).
const int nXjBinsForChi2 = ana::nUnfoldXjBins - 3;

// Bin contents treated as Poisson means (weighted contents as raw counts - the usual
// simplification; true weighted Poisson would need the unweighted entries).
void poissonBin(TH1 * nominal, TH1 * toy, int b) {
  double mean = nominal->GetBinContent(b);
  double val = (mean > 0) ? gRandom->PoissonD(mean) : 0;
  toy->SetBinContent(b, val);
  toy->SetBinError(b, sqrt(std::max(val, 0.0)));
}
TH1D * poissonToy(TH1D * nominal, const char * name) { // includes under/overflow
  TH1D * toy = (TH1D*)nominal->Clone(name);
  for (int b = 0; b <= nominal->GetNbinsX()+1; b++) poissonBin(nominal, toy, b);
  return toy;
}
TH2D * poissonToy(TH2D * nominal, const char * name) { // in-range bins, x outer (fixed random order)
  TH2D * toy = (TH2D*)nominal->Clone(name);
  for (int bx = 1; bx <= nominal->GetNbinsX(); bx++)
    for (int by = 1; by <= nominal->GetNbinsY(); by++) poissonBin(nominal, toy, nominal->GetBin(bx, by));
  return toy;
}

void toy_iterations(string systag = "nominal", string source = "resp") {
  if (source != "resp" && source != "data") { cout << "source must be \"resp\" or \"data\"" << endl; return; }
  const bool resp = (source == "resp");
  const int nToysScan = resp ? 200 : 500;
  const int nPullToys = resp ? 100 : 200;
  const char * what = resp ? "response-matrix" : "Data";

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gRandom->SetSeed(12345);

  drawer d("pythia", systag);
  string pdfPath = Form("%s/pdfs/toy_%s_iterations_%s.pdf", ana::dir(), source.c_str(), systag.c_str());

  TH1D * respReco  = d.get(Form("hrecoxj%i",ir), 1);
  TH1D * respTruth = d.get(Form("htruthxj%i",ir), 1);
  TH2D * respMatrix = d.get2d(Form("hxjresponse%i",ir), 1);
  TH1D * flatMeasured = unfold_utility::buildFullyCorrected(d.get(Form("hrecoxj%i_0",ir), 0), d.get(Form("hrecoxj%i_2",ir), 0), "data", systag);
  TH1D * hNominal = unfold_utility::unfoldOnce(respReco, respTruth, respMatrix, flatMeasured, niterPrimary, "hUnfoldNominal");

  // One toy unfolding: fluctuate the chosen input, unfold, free the toy input.
  // includeSystematics=false: only bin contents are used.
  auto unfoldToy = [&](int niter, const char * name) {
    TH1D * h;
    if (resp) {
      TH2D * m = poissonToy(respMatrix, Form("%s_in", name));
      h = unfold_utility::unfoldOnce(respReco, respTruth, m, flatMeasured, niter, name, false);
      delete m;
    } else {
      TH1D * m = poissonToy(flatMeasured, Form("%s_in", name));
      h = unfold_utility::unfoldOnce(respReco, respTruth, respMatrix, m, niter, name, false);
      delete m;
    }
    return h;
  };

  TCanvas * c = new TCanvas("c","",700,900);
  c->SaveAs(Form("%s[", pdfPath.c_str()));

  // Page 1: toy pulls (toy - nominal)/sqrt(nominal) should be a unit Gaussian.
  TH1D * hPull = new TH1D("hPull", Form(";(toy - nominal)/#sqrt{nominal};%s bins", resp ? "Response-matrix" : "Measured-spectrum"), 60, -5, 5);
  TH1 * pullRef = resp ? (TH1*)respMatrix : (TH1*)flatMeasured;
  for (int itoy = 0; itoy < nPullToys; itoy++) {
    TH1 * toy = resp ? (TH1*)poissonToy(respMatrix, Form("hToyPull_%d", itoy)) : (TH1*)poissonToy(flatMeasured, Form("hToyPull_%d", itoy));
    for (int b = 0; b < pullRef->GetNcells(); b++) {
      if (pullRef->IsBinUnderflow(b) || pullRef->IsBinOverflow(b)) continue;
      double nominal = pullRef->GetBinContent(b);
      if (nominal < 5) continue; // non-Gaussian Poisson regime; excluded from this check only
      hPull->Fill((toy->GetBinContent(b) - nominal)/sqrt(nominal));
    }
    delete toy;
  }
  TF1 * fGaus = new TF1("fGaus", "gaus", -5, 5);
  hPull->Scale(1./hPull->Integral(), "width");
  hPull->Fit(fGaus, "Q0");
  c->Clear(); c->cd(); gPad->SetTicks(1,1); gPad->SetLeftMargin(.15);
  hPull->SetMarkerStyle(20);
  hPull->GetYaxis()->SetTitle("Probability density");
  hPull->Draw("p e");
  TF1 * fUnit = new TF1("fUnit", "TMath::Gaus(x,0,1,1)", -5, 5);
  fUnit->SetLineColor(kBlue); fUnit->SetLineStyle(2); fUnit->Draw("same");
  fGaus->SetLineColor(kRed); fGaus->Draw("same");
  TLegend * lp = new TLegend(.55,.65,.88,.85);
  lp->SetLineWidth(0); lp->SetTextSize(0.03);
  lp->AddEntry(hPull, Form("Toy pulls (%d toys pooled)", nPullToys));
  lp->AddEntry(fUnit, "Unit Gaussian (expectation)");
  lp->AddEntry(fGaus, Form("Fit: #mu=%.2f, #sigma=%.2f", fGaus->GetParameter(1), fGaus->GetParameter(2)));
  lp->Draw();
  d.drawAll({resp ? "Response matrix toy validation" : "Data toy validation"},{Form("Jet R=%.1f",ana::JetRs[ir]), "Bins with nominal content > 5"}, .18, .85, 14, gPad->GetWh()*0.8);
  c->SaveAs(pdfPath.c_str());

  // Primary ensemble at niterPrimary: per-bin toy RMS = this source's uncertainty.
  int nFlatBins = flatMeasured->GetNbinsX();
  vector<double> sum(nFlatBins+2, 0), sumsq(nFlatBins+2, 0);
  vector<TH1D*> sampleToys;
  int keepEvery = std::max(1, nToys/nSampleOverlay);
  for (int itoy = 0; itoy < nToys; itoy++) {
    TH1D * hToy = unfoldToy(niterPrimary, Form("hUnfoldToy_%d", itoy));
    for (int b = 0; b <= nFlatBins+1; b++) { double v = hToy->GetBinContent(b); sum[b] += v; sumsq[b] += v*v; }
    if (itoy % keepEvery == 0 && (int)sampleToys.size() < nSampleOverlay) sampleToys.push_back(hToy);
    else delete hToy;
    if (itoy % 200 == 0) cout << "  toy " << itoy << "/" << nToys << endl;
  }
  TH1D * hNominalToyErr = (TH1D*)hNominal->Clone("hNominalToyErr");
  for (int b = 0; b <= nFlatBins+1; b++) {
    double mean = sum[b]/nToys;
    hNominalToyErr->SetBinError(b, sqrt(std::max(sumsq[b]/nToys - mean*mean, 0.0)));
  }

  // Sample toys overlaid on the nominal result, per pT bin.
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hNomPt = unfold_utility::unflattenXj(hNominal, ipt, Form("hNomPt_%d", ipt));
    TH1D * hNomDisp = unfold_utility::densityForDisplay(hNomPt, Form("hNomDisp_%d", ipt));
    double ymax = hNomDisp->GetMaximum();
    vector<TH1D*> sampleDisp;
    for (unsigned k = 0; k < sampleToys.size(); k++) {
      TH1D * hp = unfold_utility::unflattenXj(sampleToys[k], ipt, Form("hSamplePt_%d_%d", ipt, k));
      sampleDisp.push_back(unfold_utility::densityForDisplay(hp, Form("hSampleDisp_%d_%d", ipt, k)));
      ymax = std::max(ymax, sampleDisp.back()->GetMaximum());
      delete hp;
    }
    c->Clear(); c->cd(); gPad->SetTicks(1,1); gPad->SetLeftMargin(.15);
    hNomDisp->SetLineColor(kRed); hNomDisp->SetLineWidth(3);
    hNomDisp->GetYaxis()->SetRangeUser(0, ymax*1.3);
    hNomDisp->GetXaxis()->SetTitle("x_{J#gamma}");
    hNomDisp->Draw("hist");
    for (auto h : sampleDisp) { h->SetLineColorAlpha(kGray+2, 0.35); h->SetLineWidth(1); h->Draw("hist same"); }
    hNomDisp->Draw("hist same");
    TLegend * lo = new TLegend(.5,.6,.85,.75);
    lo->SetLineWidth(0); lo->SetTextSize(0.03);
    lo->AddEntry(hNomDisp, "Nominal unfolded", "l");
    lo->AddEntry(sampleDisp[0], Form("%d of %d %s toys", (int)sampleDisp.size(), nToys, what), "l");
    lo->Draw();
    d.drawAll({"p+p Run24 Data"},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
        Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV", ana::JetRs[ir], ana::jet_calib_pt_cut[ir])}, .5, .85, 14, gPad->GetWh()*0.8);
    c->SaveAs(pdfPath.c_str());
    delete hNomPt; delete hNomDisp;
    for (auto h : sampleDisp) delete h;
  }
  for (auto h : sampleToys) delete h;

  // Nominal result with the toy-RMS band, per pT bin.
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hNomPt = unfold_utility::unflattenXj(hNominalToyErr, ipt, Form("hNomErrPt_%d", ipt));
    TH1D * hNomDisp = unfold_utility::densityForDisplay(hNomPt, Form("hNomErrDisp_%d", ipt));
    c->Clear(); c->cd(); gPad->SetTicks(1,1); gPad->SetLeftMargin(.15);
    hNomDisp->SetMarkerStyle(20);
    hNomDisp->SetFillColorAlpha(resp ? kAzure+7 : kGreen+1, 0.35);
    hNomDisp->GetYaxis()->SetRangeUser(0, hNomDisp->GetMaximum()*1.4);
    hNomDisp->GetXaxis()->SetTitle("x_{J#gamma}");
    hNomDisp->Draw("e2");
    hNomDisp->Draw("p e same");
    TLegend * le = new TLegend(.5,.62,.85,.75);
    le->SetLineWidth(0); le->SetTextSize(0.03);
    le->AddEntry(hNomDisp, "Nominal unfolded #pm toy RMS", "lf");
    le->Draw();
    d.drawAll({"p+p Run24 Data"},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
        Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV, %d toys", ana::JetRs[ir], ana::jet_calib_pt_cut[ir], nToys)}, .5, .85, 14, gPad->GetWh()*0.8);
    c->SaveAs(pdfPath.c_str());
    delete hNomPt; delete hNomDisp;
  }

  // Toy chi2/NDF vs niter: mean over used bins of (toy RMS/mean)^2 - same units as the
  // observed iteration-to-iteration chi2/NDF below (draw_iteration_halfclosure.C's metric).
  vector<double> chi2ToyByIter, pairChi2ByIter;
  vector<TH1D*> nomByIter;
  for (int iter : iterationsToScan) {
    nomByIter.push_back(unfold_utility::unfoldOnce(respReco, respTruth, respMatrix, flatMeasured, iter, Form("hNomIter_%d", iter)));
    vector<double> s(nFlatBins+2, 0), s2(nFlatBins+2, 0);
    for (int itoy = 0; itoy < nToysScan; itoy++) {
      TH1D * hToy = unfoldToy(iter, Form("hUnfoldToyScan_%d_%d", iter, itoy));
      for (int b = 0; b <= nFlatBins+1; b++) { double v = hToy->GetBinContent(b); s[b] += v; s2[b] += v*v; }
      delete hToy;
    }
    double chi2Sum = 0;
    int n = 0;
    for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
      for (int ixj = 0; ixj < nXjBinsForChi2; ixj++) {
        int flatbin = ipt*(ana::nUnfoldXjBins+2) + ixj + 2;
        double mean = s[flatbin]/nToysScan;
        if (mean > 0) { chi2Sum += (s2[flatbin]/nToysScan - mean*mean)/(mean*mean); n++; }
      }
    }
    chi2ToyByIter.push_back(n > 0 ? chi2Sum/n : 0);
    cout << "  niter=" << iter << ": " << what << "-toy chi2/NDF = " << chi2ToyByIter.back() << endl;
  }
  // Observed chi2/NDF of (iter n - iter n-1)/(iter n-1) from the nominal chain (iter 0 = measured).
  for (unsigned k = 0; k < iterationsToScan.size(); k++) {
    double chi2 = 0;
    int ndf = 0;
    for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
      TH1D * uPrev = unfold_utility::unflattenXj(k == 0 ? flatMeasured : nomByIter[k-1], ipt, Form("hPairPrev_%d_%d", k, ipt));
      TH1D * uNext = unfold_utility::unflattenXj(nomByIter[k], ipt, Form("hPairNext_%d_%d", k, ipt));
      for (int b = 1; b <= nXjBinsForChi2; b++) {
        double up = uPrev->GetBinContent(b);
        if (up <= 0) continue;
        double rel = (uNext->GetBinContent(b) - up)/up;
        chi2 += rel*rel;
        ndf++;
      }
      delete uPrev; delete uNext;
    }
    pairChi2ByIter.push_back(ndf > 0 ? chi2/ndf : 0);
    cout << "  niter=" << iterationsToScan[k] << ": observed chi2/NDF (vs iter " << (k == 0 ? 0 : iterationsToScan[k-1]) << ") = " << pairChi2ByIter.back() << endl;
  }
  for (auto h : nomByIter) delete h;
  c->SaveAs(Form("%s]", pdfPath.c_str()));

  TGraph * gToy = new TGraph(iterationsToScan.size());
  TGraph * gPair = new TGraph(iterationsToScan.size());
  for (unsigned k = 0; k < iterationsToScan.size(); k++) {
    gToy->SetPoint(k, iterationsToScan[k], chi2ToyByIter[k]);
    gPair->SetPoint(k, iterationsToScan[k], pairChi2ByIter[k]);
  }
  // The combined page is drawn in a fresh ROOT process (TGraph log-scale painting dropped
  // points near the frame edge when drawn here); it uses whichever source files exist.
  string respFile = Form("%s/pdfs/.toy_resp_chi2_data_%s.root", ana::dir(), systag.c_str());
  string dataFile = Form("%s/pdfs/.toy_data_chi2_data_%s.root", ana::dir(), systag.c_str());
  string chi2Pdf  = Form("%s/pdfs/toy_iterations_chi2_%s.pdf", ana::dir(), systag.c_str());
  TFile * fOut = TFile::Open((resp ? respFile : dataFile).c_str(), "RECREATE");
  gToy->Write("gToy");
  gPair->Write("gPair");
  fOut->Close();
  gSystem->Exec(Form("cd $GAMMAJET_UNFOLD/drawing && root -b -l -q 'plot_toy_chi2_combined.C(\"%s\",\"%s\",\"%s\",\"Jet R=%.1f\")'",
    respFile.c_str(), dataFile.c_str(), chi2Pdf.c_str(), ana::JetRs[ir]));
  cout << "Done. Wrote " << pdfPath << " and " << chi2Pdf << endl;
}
