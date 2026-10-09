#include "../src/ana.h"
R__LOAD_LIBRARY(libgammajet_unfold.so);
// Overlays toy_iterations.C's response-toy and Data-toy chi2/NDF, their sum with the observed
// "iter n vs n-1" curve (chi2/NDF is already a variance, so quadrature is a plain sum; from iter 2),
// and the observed curve itself. Missing inputs are skipped. Linear y over [yAxisMin, yAxisMax].
// Defaults are the nominal paths: root -b -l -q plot_toy_chi2_combined.C
const double yAxisMin = 0, yAxisMax = 1; // iteration-1 points (vs the measured spectrum) run off the top

void plot_toy_chi2_combined(
    const char * respFile = ana::path("pdfs/.toy_resp_chi2_data_nominal.root"),
    const char * dataFile = ana::path("pdfs/.toy_data_chi2_data_nominal.root"),
    const char * outPdf   = ana::path("pdfs/toy_iterations_chi2_nominal.pdf"),
    const char * subtitle = "Jet R=0.4") {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TGraph * gResp = nullptr;
  TGraph * gData = nullptr;
  TGraph * gPair = nullptr;

  if (!gSystem->AccessPathName(respFile)) {
    TFile * f = TFile::Open(respFile);
    gResp = (TGraph*)f->Get("gToy");
    if (!gPair) gPair = (TGraph*)f->Get("gPair");
  }
  if (!gSystem->AccessPathName(dataFile)) {
    TFile * f = TFile::Open(dataFile);
    gData = (TGraph*)f->Get("gToy");
    if (!gPair) gPair = (TGraph*)f->Get("gPair");
  }

  // gResp/gData point i is iter i+1; gPair point k is iter k+2.
  TGraph * gCombined = nullptr;
  if (gResp && gData && gPair && gResp->GetN() == gData->GetN() && gResp->GetN() == gPair->GetN()) {
    gCombined = new TGraph(gPair->GetN());
    for (int k = 0; k < gPair->GetN(); k++) {
      double xr, yr, xd, yd, xp, yp;
      gResp->GetPoint(k, xr, yr);
      gData->GetPoint(k, xd, yd);
      gPair->GetPoint(k, xp, yp);
      gCombined->SetPoint(k, xp, yr + yd + yp);
    }
  }

  struct Curve { TGraph * g; int color; int marker; int style; const char * label; };
  vector<Curve> curves;
  if (gResp)     curves.push_back({gResp, kAzure+2, 20, 1, "Response-matrix toy"});
  if (gData)     curves.push_back({gData, kGreen+2, 22, 1, "Data toy"});
  if (gPair)     curves.push_back({gPair, kBlack, 21, 1, "Observed: iter n vs iter n-1"});
  if (gCombined) curves.push_back({gCombined, kMagenta+1, 23, 9, "Quadrature sum (resp #oplus data #oplus iter)"});

  double yMin = 1e300, yMax = -1e300;
  for (auto & cv : curves)
    for (int i = 0; i < cv.g->GetN(); i++) {
      double x,y; cv.g->GetPoint(i,x,y);
      if (y > 0) { yMin = std::min(yMin,y); yMax = std::max(yMax,y); }
    }

  TCanvas * c = new TCanvas("c","",700,900);
  gPad->SetTicks(1,1);
  gPad->SetLeftMargin(.15);

  TLegend * lc = new TLegend(.17,.68,.6,.85);
  lc->SetLineWidth(0);
  lc->SetTextSize(0.025);

  for (unsigned i = 0; i < curves.size(); i++) {
    Curve & cv = curves[i];
    cv.g->SetMarkerStyle(cv.marker);
    cv.g->SetMarkerColor(cv.color);
    cv.g->SetLineColor(cv.color);
    cv.g->SetLineWidth(2);
    cv.g->SetLineStyle(cv.style);
    if (i == 0) {
      cv.g->SetMinimum(yAxisMin);
      cv.g->SetMaximum(yAxisMax);
      cv.g->Draw("APL");
      cv.g->GetXaxis()->SetTitle("Bayesian unfolding iterations");
      cv.g->GetYaxis()->SetTitle("#chi^{2}/NDF");
    } else {
      cv.g->Draw("PL same");
    }
    lc->AddEntry(cv.g, cv.label, "lp");
  }
  lc->Draw();

  TLatex * t1 = new TLatex(.55,.85,"#bf{#it{sPHENIX}} #kern[0.5]{Internal}");
  t1->SetNDC(); t1->SetTextFont(43); t1->SetTextSize(17.5); t1->Draw();
  TLatex * t2 = new TLatex(.55,.826,"p+p Run24 Data");
  t2->SetNDC(); t2->SetTextFont(43); t2->SetTextSize(17.5); t2->Draw();
  TLatex * t3 = new TLatex(.55,.802,subtitle);
  t3->SetNDC(); t3->SetTextFont(43); t3->SetTextSize(14); t3->Draw();

  c->SaveAs(outPdf);
  cout << "plot_toy_chi2_combined: wrote " << outPdf
       << " (resp=" << (gResp ? gResp->GetN() : 0)
       << ", data=" << (gData ? gData->GetN() : 0)
       << ", combined=" << (gCombined ? gCombined->GetN() : 0)
       << ", pair=" << (gPair ? gPair->GetN() : 0) << ")" << endl;
}
