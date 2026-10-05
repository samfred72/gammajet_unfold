#include "../src/ana.h"
#include "../src/drawer.h"
#include "../src/unfold_utility.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include <map>
#include <iomanip>
#include <algorithm>
// Load explicitly: the sibling gammajet project's libgammajet.so has same-named classes.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// RooUnfoldBayes's analytic error vs an ATLAS-style Poisson-toy bootstrap (PLB 774 (2017) 379,
// Sec. 6) at draw_final_result.C's working point (nominal, R=0.4, niter=2): Data toys and
// response toys, combined in quadrature.
//
// Response toys are weighted Poisson: Poisson(N_eff), N_eff = content^2/error^2 (Kish), rescaled
// by content/N_eff, so mean and variance match (content, error^2). The response is cross-section
// weighted and some bins have only a few effective entries, so Poisson(content) is far too narrow.
//
// Data toys are plain Poisson on the raw region-A/C counts. Kish does not model the signed
// combination coeffA*A - coeffC*C: N_eff < 1 where A and C nearly cancel.
double poissonToyValue(double content) {
  return gRandom->PoissonD(std::max(content, 0.0));
}

TH1D * poissonToyHist(TH1D * nominal, const char * name) {
  TH1D * toy = (TH1D*)nominal->Clone(name);
  for (int b = 0; b <= nominal->GetNbinsX()+1; b++) {
    double val = poissonToyValue(nominal->GetBinContent(b));
    toy->SetBinContent(b, val);
    toy->SetBinError(b, sqrt(fabs(val)));
  }
  return toy;
}

// Toy purity combination with coeffA/coeffC fixed per pT bin, as purityCorrect's error formula
// assumes (it treats them as constants). The purity-measurement term is added deterministically
// (purityErrSqByFlatBin), not by re-drawing pA/pC, which is wider than the analytic two-point
// envelope.
TH1D * toyMeasuredFixedCoeffs(TH1D * toyFlatA, TH1D * toyFlatC, const vector<float> & coeffAByPt,
    const vector<float> & coeffCByPt, const char * name) {
  TH1D * flatCorrected = (TH1D*)toyFlatA->Clone(name);
  flatCorrected->Reset("ICES");
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
    TH1D * A = unfold_utility::unflattenXj(toyFlatA, ipt, Form("htmpToyA_%s_pt%d", name, ipt));
    TH1D * C = unfold_utility::unflattenXj(toyFlatC, ipt, Form("htmpToyC_%s_pt%d", name, ipt));
    for (int i = 1; i <= A->GetNbinsX(); i++) {
      float content = coeffAByPt[ipt]*A->GetBinContent(i) - coeffCByPt[ipt]*C->GetBinContent(i);
      A->SetBinContent(i, content);
      A->SetBinError(i, sqrt(fabs(content))); // only the central value is read
    }
    unfold_utility::reflattenXj(A, ipt, flatCorrected);
    delete A; delete C;
  }
  return flatCorrected;
}

const int ir = 2; // nominal R=0.4
const int nPtBinsUsed = ana::nPtBinsUsed;
const int nToys = 10000; // toy count per side (Data, response)
const int niterate = 2; // matches draw_final_result.C's actual working point
const int nXjBinsForChi2 = ana::nUnfoldXjBins - 3; // exclude the low-stat tail, same convention as elsewhere

// false: skip the response-side toys.
const bool includeResponseUncertainty = true;

// Kish weighted-Poisson draw for one bin.
double weightedPoissonToyValue(double content, double error) {
  if (error <= 0) return content; // never-filled bin: leave unchanged
  double neff = content*content/(error*error);
  if (neff <= 0) return content; // guard only - shouldn't occur given error>0
  double meanWeight = content/neff;
  double k = gRandom->PoissonD(neff);
  return k*meanWeight;
}

TH1D * weightedPoissonToyHist(TH1D * nominal, const char * name) {
  TH1D * toy = (TH1D*)nominal->Clone(name);
  for (int b = 0; b <= nominal->GetNbinsX()+1; b++) {
    double val = weightedPoissonToyValue(nominal->GetBinContent(b), nominal->GetBinError(b));
    toy->SetBinContent(b, val);
    toy->SetBinError(b, sqrt(fabs(val)));
  }
  return toy;
}

TH2D * weightedPoissonToyMatrix(TH2D * nominal, const char * name) {
  TH2D * toy = (TH2D*)nominal->Clone(name);
  for (int bx = 1; bx <= nominal->GetNbinsX(); bx++) {
    for (int by = 1; by <= nominal->GetNbinsY(); by++) {
      double val = weightedPoissonToyValue(nominal->GetBinContent(bx,by), nominal->GetBinError(bx,by));
      toy->SetBinContent(bx, by, val);
      toy->SetBinError(bx, by, sqrt(fabs(val)));
    }
  }
  return toy;
}

void draw_toy_vs_analytic(string systag = "nominal") {
  // Unique names: thousands of same-named objects in gDirectory slow every New()/Clone().
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gRandom->SetSeed(12345); // same seed toy_iterations.C use

  drawer d("pythia", systag);
  string pdfPath = Form("%s/pdfs/toy_vs_analytic_%s.pdf", ana::dir(), systag.c_str());
  string outfilename = Form("%s/hists/toy_vs_analytic_%s.root", ana::dir(), systag.c_str());

  TH1D * respRecoTemplate  = d.get(Form("hrecoxj%i",ir), 1);
  TH1D * respTruthTemplate = d.get(Form("htruthxj%i",ir), 1);
  TH2D * respMatrix2D      = d.get2d(Form("hxjresponse%i",ir), 1);

  TH1D * flatA = d.get(Form("hrecoxj%i_0",ir), 0);
  TH1D * flatC = d.get(Form("hrecoxj%i_2",ir), 0);
  TH1D * flatMeasured = unfold_utility::buildFullyCorrected(flatA, flatC, "data", systag);

  TH1D * hNominal = unfold_utility::unfoldOnce(respRecoTemplate, respTruthTemplate, respMatrix2D, flatMeasured, niterate, "hUnfoldNominal");
  int nFlatBins = flatMeasured->GetNbinsX();

  // Coefficients from the nominal (non-toy) counts.
  vector<float> coeffAByPt(ana::nPtBins), coeffCByPt(ana::nPtBins);
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
    float ptlow = ana::ptBins[ipt], pthigh = ana::ptBins[ipt+1];
    float pA = ana::getPurity(ptlow, pthigh, systag, ir);
    float pC = ana::getPurityC(ptlow, pthigh, systag, ir);
    TH1D * Anom = unfold_utility::unflattenXj(flatA, ipt, Form("hNomA_pt%d", ipt));
    TH1D * Cnom = unfold_utility::unflattenXj(flatC, ipt, Form("hNomC_pt%d", ipt));
    unfold_utility::purityCorrectCoeffs(pA, pC, Anom->Integral(), Cnom->Integral(), coeffAByPt[ipt], coeffCByPt[ipt]);
    delete Anom; delete Cnom;
  }

  // Purity-measurement variance per bin: flatMeasured's error^2 minus the counting term, which
  // reproduces the analytic purity piece exactly.
  vector<double> purityErrSqByFlatBin(nFlatBins+2, 0.0);
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
    for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) {
      int flatbin = ipt*(ana::nUnfoldXjBins+2) + ixj + 1 + 1;
      double ae = flatA->GetBinError(flatbin);
      double ce = flatC->GetBinError(flatbin);
      double statErr = sqrt(pow(coeffAByPt[ipt]*ae,2) + pow(coeffCByPt[ipt]*ce,2));
      double measuredErr = flatMeasured->GetBinError(flatbin);
      purityErrSqByFlatBin[flatbin] = std::max(measuredErr*measuredErr - statErr*statErr, 0.0);
    }
  }

  // Every toy value is kept for the percentile intervals.
  vector<vector<double>> toyValD(nFlatBins+2, vector<double>(nToys));
  vector<vector<double>> toyValR(nFlatBins+2, vector<double>(nToys));

  // ---- Data-side toys (raw region A/C, response fixed) ----
  vector<double> sumD(nFlatBins+2, 0), sumsqD(nFlatBins+2, 0);
  cout << "Running " << nToys << " Data-side plain-Poisson toys at niter=" << niterate << "..." << endl;
  for (int itoy = 0; itoy < nToys; itoy++) {
    TH1D * toyFlatA = poissonToyHist(flatA, Form("hToyA_%d", itoy));
    TH1D * toyFlatC = poissonToyHist(flatC, Form("hToyC_%d", itoy));
    TH1D * toyData = toyMeasuredFixedCoeffs(toyFlatA, toyFlatC, coeffAByPt, coeffCByPt, Form("data_toy_%d", itoy));
    for (int b = 0; b <= nFlatBins+1; b++) {
      if (purityErrSqByFlatBin[b] <= 0) continue;
      toyData->SetBinContent(b, toyData->GetBinContent(b) + gRandom->Gaus(0, sqrt(purityErrSqByFlatBin[b])));
    }
    // includeSystematics=false: only the content is read.
    TH1D * hToy = unfold_utility::unfoldOnce(respRecoTemplate, respTruthTemplate, respMatrix2D, toyData, niterate, Form("hUnfoldToyD_%d", itoy), false);
    for (int b = 0; b <= nFlatBins+1; b++) { double v = hToy->GetBinContent(b); sumD[b] += v; sumsqD[b] += v*v; toyValD[b][itoy] = v; }
    delete hToy; delete toyData; delete toyFlatA; delete toyFlatC;
    if (itoy % 200 == 0) cout << "  data toy " << itoy << "/" << nToys << endl;
  }

  // ---- Response-side toys (respMatrix2D, Data fixed) ----
  // When skipped, toyValR stays zero; guarded at the percentile call.
  vector<double> sumR(nFlatBins+2, 0), sumsqR(nFlatBins+2, 0);
  if (includeResponseUncertainty) {
    cout << "Running " << nToys << " response-matrix weighted-Poisson toys at niter=" << niterate << "..." << endl;
    for (int itoy = 0; itoy < nToys; itoy++) {
      TH2D * toyMat = weightedPoissonToyMatrix(respMatrix2D, Form("hToyMat_%d", itoy));
      TH1D * hToy = unfold_utility::unfoldOnce(respRecoTemplate, respTruthTemplate, toyMat, flatMeasured, niterate, Form("hUnfoldToyR_%d", itoy), false);
      for (int b = 0; b <= nFlatBins+1; b++) { double v = hToy->GetBinContent(b); sumR[b] += v; sumsqR[b] += v*v; toyValR[b][itoy] = v; }
      delete hToy; delete toyMat;
      if (itoy % 200 == 0) cout << "  resp toy " << itoy << "/" << nToys << endl;
    }
  } else {
    cout << "Skipping response-matrix toys (includeResponseUncertainty=false) - data-side uncertainty only." << endl;
  }

  // 16th/84th-percentile interval around the nominal content (robust to outlier toys in
  // low-count bins, as in puritymaker.C). Returns positive magnitudes.
  auto percentileInterval = [](vector<double> vals, double center, double & lo, double & hi) {
    std::sort(vals.begin(), vals.end());
    int n = (int)vals.size();
    int i16 = std::max(0, (int)std::lround(0.16*(n-1)));
    int i84 = std::min(n-1, (int)std::lround(0.84*(n-1)));
    lo = std::max(center - vals[i16], 0.0);
    hi = std::max(vals[i84] - center, 0.0);
  };

  // Toy-error copy of hNominal plus per-bin ratios; RMS and percentile numbers side by side.
  TH1D * hNominalToyErr = (TH1D*)hNominal->Clone("hNominalToyErr");
  vector<double> percLoByFlatBin(nFlatBins+2, 0.0), percHiByFlatBin(nFlatBins+2, 0.0);
  vector<double> ratioVals, ratioValsPerc, ratioX;
  vector<int> ratioIpt;
  int idx = 0;
  cout << "\n" << std::left << std::setw(6) << "ptbin" << std::setw(6) << "xjbin"
       << std::setw(12) << "content" << std::setw(12) << "analytic" << std::setw(12) << "dataToy"
       << std::setw(12) << "respToy" << std::setw(14) << "combinedToy" << std::setw(10) << "ratio"
       << std::setw(12) << "percLo" << std::setw(12) << "percHi" << std::setw(10) << "percRatio" << endl;
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    // Every xJ bin, so the plotted tail does not keep the analytic error.
    for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) {
      int flatbin = ipt*(ana::nUnfoldXjBins+2) + ixj + 1 + 1;
      double meanD = sumD[flatbin]/nToys, varD = sumsqD[flatbin]/nToys - meanD*meanD;
      double meanR = sumR[flatbin]/nToys, varR = sumsqR[flatbin]/nToys - meanR*meanR;
      double dataRMS = sqrt(std::max(varD,0.0));
      double respRMS = sqrt(std::max(varR,0.0));
      double combinedToy = sqrt(dataRMS*dataRMS + respRMS*respRMS);
      hNominalToyErr->SetBinError(flatbin, combinedToy);

      double content = hNominal->GetBinContent(flatbin);
      // Data and response combined per side in quadrature.
      double dataLo, dataHi, respLo = 0, respHi = 0;
      percentileInterval(toyValD[flatbin], content, dataLo, dataHi);
      if (includeResponseUncertainty) percentileInterval(toyValR[flatbin], content, respLo, respHi);
      double percLo = sqrt(dataLo*dataLo + respLo*respLo);
      double percHi = sqrt(dataHi*dataHi + respHi*respHi);
      percLoByFlatBin[flatbin] = percLo;
      percHiByFlatBin[flatbin] = percHi;

      // The summary excludes the low-stat tail (as draw_covariance_chi2.C).
      if (ixj >= nXjBinsForChi2) continue;
      double analyticErr = hNominal->GetBinError(flatbin);
      if (content <= 0 || analyticErr <= 0) continue;
      double ratio = combinedToy/analyticErr;
      double percRatio = std::max(percLo,percHi)/analyticErr;
      ratioVals.push_back(ratio);
      ratioValsPerc.push_back(percRatio);
      ratioX.push_back(idx++);
      ratioIpt.push_back(ipt);
      printf("%-6d%-6d%-12.4f%-12.4f%-12.4f%-12.4f%-14.4f%-10.3f%-12.4f%-12.4f%-10.3f\n",
        ipt, ixj, content, analyticErr, dataRMS, respRMS, combinedToy, ratio, percLo, percHi, percRatio);
    }
  }
  double sum_r = 0; for (double r : ratioVals) sum_r += r;
  double sum_rp = 0; for (double r : ratioValsPerc) sum_rp += r;
  vector<double> sortedRatios = ratioVals;
  vector<double> sortedRatiosPerc = ratioValsPerc;
  std::sort(sortedRatios.begin(), sortedRatios.end());
  std::sort(sortedRatiosPerc.begin(), sortedRatiosPerc.end());
  cout << "\nSummary over " << ratioVals.size() << " used bins (combined_toy / analytic):" << endl;
  cout << "  mean   = " << sum_r/ratioVals.size() << endl;
  cout << "  median = " << sortedRatios[sortedRatios.size()/2] << endl;
  cout << "  min    = " << sortedRatios.front() << endl;
  cout << "  max    = " << sortedRatios.back() << endl;
  cout << "\nSummary over " << ratioValsPerc.size() << " used bins (percentile combined / analytic):" << endl;
  cout << "  mean   = " << sum_rp/ratioValsPerc.size() << endl;
  cout << "  median = " << sortedRatiosPerc[sortedRatiosPerc.size()/2] << endl;
  cout << "  min    = " << sortedRatiosPerc.front() << endl;
  cout << "  max    = " << sortedRatiosPerc.back() << endl;

  TCanvas * c = new TCanvas("c","",900,700);
  c->SaveAs(Form("%s[", pdfPath.c_str()));

  // ---- Page 1: toy/analytic ratio for every used bin, colored by pT bin ----
  c->Clear();
  c->cd();
  gPad->SetTicks(1,1);
  gPad->SetLeftMargin(.12);
  gPad->SetBottomMargin(.12);
  TH1F * frame = gPad->DrawFrame(-1, 0.5, (double)ratioVals.size(), 2.5);
  frame->GetYaxis()->SetTitle("Combined toy / RooUnfoldBayes analytic error");
  frame->GetXaxis()->SetTitle("Bin index (grouped by p_{T}^{#gamma} bin, low-to-high x_{J#gamma} within each)");
  int colors[7] = {kRed+1, kOrange+1, kSpring+2, kGreen+2, kCyan+2, kAzure+1, kViolet+1};
  map<int,TGraph*> byPt;
  map<int,TGraph*> byPtPerc;
  for (size_t k = 0; k < ratioVals.size(); k++) {
    int ipt = ratioIpt[k];
    if (!byPt.count(ipt)) byPt[ipt] = new TGraph();
    if (!byPtPerc.count(ipt)) byPtPerc[ipt] = new TGraph();
    byPt[ipt]->SetPoint(byPt[ipt]->GetN(), ratioX[k], ratioVals[k]);
    // Open circles: the percentile series, same color as the RMS point.
    byPtPerc[ipt]->SetPoint(byPtPerc[ipt]->GetN(), ratioX[k], ratioValsPerc[k]);
  }
  TLegend * legR = new TLegend(.65,.5,.88,.88);
  legR->SetLineWidth(0);
  legR->SetTextSize(0.025);
  int ic = 0;
  for (auto & pr : byPt) {
    pr.second->SetMarkerStyle(20);
    pr.second->SetMarkerColor(colors[ic % 7]);
    pr.second->SetLineColor(colors[ic % 7]);
    pr.second->SetMarkerSize(1.1);
    pr.second->Draw("p same");
    legR->AddEntry(pr.second, Form("%.0f-%.0f GeV", ana::ptBins[pr.first], ana::ptBins[pr.first+1]), "p");
    ic++;
  }
  ic = 0;
  for (auto & pr : byPtPerc) {
    pr.second->SetMarkerStyle(24);
    pr.second->SetMarkerColor(colors[ic % 7]);
    pr.second->SetMarkerSize(1.1);
    pr.second->Draw("p same");
    ic++;
  }
  TLine * lone = new TLine(-1, 1.0, (double)ratioVals.size(), 1.0);
  lone->SetLineStyle(9);
  lone->SetLineColor(kBlack);
  lone->Draw("same");
  TGraph * legPercMarker = new TGraph();
  legPercMarker->SetMarkerStyle(24);
  legPercMarker->SetMarkerColor(kBlack);
  legR->AddEntry(legPercMarker, "16th/84th percentile ratio (open)", "p");
  legR->Draw();
  d.drawAll({"p+p Run24 Data - stat. uncertainty check"},
      {Form("Jet R=%.1f, niter=%d, %d toys%s", ana::JetRs[ir], niterate, nToys,
            includeResponseUncertainty ? "/side (weighted-Poisson)" : " (DATA-SIDE ONLY - TEMPORARY)"),
       "Dashed line: perfect agreement (ratio=1)"}, .15, .85, 14, gPad->GetWh()*0.8);
  c->SaveAs(pdfPath.c_str());

  // ---- Pages 2..: per-pT-bin spectrum with analytic vs toy errors ----
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hAnaPt = unfold_utility::unflattenXj(hNominal, ipt, Form("hAnaPt_%d", ipt));
    TH1D * hToyPt = unfold_utility::unflattenXj(hNominalToyErr, ipt, Form("hToyPt_%d", ipt));
    TH1D * hAnaDisp = unfold_utility::densityForDisplay(hAnaPt, Form("hAnaDisp_%d", ipt));
    TH1D * hToyDisp = unfold_utility::densityForDisplay(hToyPt, Form("hToyDisp_%d", ipt));

    // Percentile interval as density boxes (asymmetric, as draw_final_result.C's gSystBox).
    TGraphAsymmErrors * gPercBox = new TGraphAsymmErrors(ana::nUnfoldXjBins);
    TGraph * gPercPts = new TGraph(ana::nUnfoldXjBins);
    for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) {
      double xlo = ana::unfoldXjBins[ixj], xhi = ana::unfoldXjBins[ixj+1];
      double width = xhi - xlo, xc = 0.5*(xlo+xhi), halfw = 0.5*width;
      int flatbin = ipt*(ana::nUnfoldXjBins+2) + ixj + 1 + 1;
      double content = hNominal->GetBinContent(flatbin) / width;
      double lo = percLoByFlatBin[flatbin] / width, hi = percHiByFlatBin[flatbin] / width;
      gPercBox->SetPoint(ixj, xc, content);
      gPercBox->SetPointError(ixj, halfw, halfw, lo, hi);
      gPercPts->SetPoint(ixj, xc, content);
    }

    c->Clear();
    c->cd();
    gPad->SetTicks(1,1);
    gPad->SetLeftMargin(.13);
    hAnaDisp->SetLineColor(kAzure+2);
    hAnaDisp->SetMarkerColor(kAzure+2);
    hAnaDisp->SetFillColorAlpha(kAzure+2, 0.30);
    hAnaDisp->SetMarkerStyle(0);
    hAnaDisp->GetXaxis()->SetTitle("x_{J#gamma}");
    hAnaDisp->GetYaxis()->SetTitle("Shape-normalized counts / bin width");
    hAnaDisp->GetYaxis()->SetRangeUser(0, hAnaDisp->GetMaximum()*1.5);
    hAnaDisp->Draw("e2");

    hToyDisp->SetLineColor(kRed+1);
    hToyDisp->SetMarkerColor(kRed+1);
    hToyDisp->SetFillColorAlpha(kRed+1, 0.45);
    hToyDisp->SetMarkerStyle(0);
    hToyDisp->Draw("e2 same");

    gPercBox->SetFillColorAlpha(kGreen+2, 0.35);
    gPercBox->SetLineColor(kWhite);
    gPercBox->Draw("2 same");

    hAnaDisp->SetMarkerStyle(20);
    hAnaDisp->SetMarkerColor(kBlack);
    hAnaDisp->Draw("p same");

    gPercPts->SetMarkerStyle(24);
    gPercPts->SetMarkerColor(kGreen+3);
    gPercPts->Draw("p same");

    TLegend * le = new TLegend(.5,.65,.85,.85);
    le->SetLineWidth(0);
    le->SetTextSize(0.03);
    le->AddEntry(hAnaDisp, "Nominal #pm RooUnfold analytic error", "lf");
    le->AddEntry(hToyDisp, Form("Nominal #pm toy RMS (%d toys%s)", nToys,
        includeResponseUncertainty ? "/side, weighted-Poisson" : ", data-side only - TEMPORARY"), "lf");
    le->AddEntry(gPercBox, "Nominal #pm toy 16th/84th percentile", "fp");
    le->Draw();
    d.drawAll({"p+p Run24 Data"},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
        Form("Jet R=%.1f, niter=%d", ana::JetRs[ir], niterate)}, .5, .55, 14, gPad->GetWh()*0.8);
    c->SaveAs(pdfPath.c_str());
    delete hAnaPt; delete hToyPt; delete hAnaDisp; delete hToyDisp; delete gPercBox; delete gPercPts;
  }

  c->SaveAs(Form("%s]", pdfPath.c_str()));

  // Percentile interval for every xJ bin.
  TGraphAsymmErrors * gPercentileToyErr = new TGraphAsymmErrors();
  gPercentileToyErr->SetName("gNominalToyErrPercentile");
  gPercentileToyErr->SetTitle(";flat bin;content");
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) {
      int flatbin = ipt*(ana::nUnfoldXjBins+2) + ixj + 1 + 1;
      int n = gPercentileToyErr->GetN();
      gPercentileToyErr->SetPoint(n, flatbin, hNominal->GetBinContent(flatbin));
      gPercentileToyErr->SetPointError(n, 0, 0, percLoByFlatBin[flatbin], percHiByFlatBin[flatbin]);
    }
  }

  TFile * fout = TFile::Open(outfilename.c_str(), "RECREATE");
  hNominal->Write();
  hNominalToyErr->Write();
  gPercentileToyErr->Write();
  fout->Close();

  cout << "\nWrote " << pdfPath << endl;
  cout << "Wrote " << outfilename << endl;
}
