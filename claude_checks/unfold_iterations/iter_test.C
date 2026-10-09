#include "../../src/ana.h"
#include "../../src/drawer.h"
#include "../../src/unfold_utility.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Test: are the iteration >10 spikes caused by negative bins in the purity-corrected measured input?
// Same metrics as draw_refolding.C (refold chi2/NDF, reported pT bins, first 13 xJ bins, measured
// error) and toy_iterations.C (response-toy mean (RMS/mean)^2), with the measured input as is
// ("asis") and with negative bins set to 0 ("clamped").
const int ir = 2;
const int G = ana::nUnfoldXjBins + 2;
const int nXj = ana::nUnfoldXjBins - 3;
const int nToys = 60;

double refoldChi2(RooUnfoldResponse * R, TH1D * meas, TH1D * fakes, int it) {
  TH1D * u = unfold_utility::unfoldOnce(R, meas, it, "uR", false);
  TH1D * f = (TH1D*)R->ApplyToTruth(u, "fR"); f->Add(fakes);
  double c = 0; int n = 0;
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin + ana::nPtBinsUsed; ipt++)
    for (int ixj = 0; ixj < nXj; ixj++) {
      int b = ipt*G + ixj + 2; double e = meas->GetBinError(b); if (e <= 0) continue;
      c += pow(f->GetBinContent(b) - meas->GetBinContent(b), 2)/(e*e); n++;
    }
  delete u; delete f;
  return n ? c/n : 0;
}
double respToyChi2(TH1D * reco, TH1D * truth, TH2D * M, TH1D * meas, int it) {
  int nb = meas->GetNbinsX();
  vector<double> s(nb+2, 0), s2(nb+2, 0);
  for (int t = 0; t < nToys; t++) {
    TH2D * m = (TH2D*)M->Clone("mtoy");
    for (int bx = 1; bx <= M->GetNbinsX(); bx++) for (int by = 1; by <= M->GetNbinsY(); by++) {
      double v = M->GetBinContent(bx, by); m->SetBinContent(bx, by, v > 0 ? gRandom->PoissonD(v) : 0);
    }
    TH1D * u = unfold_utility::unfoldOnce(reco, truth, m, meas, it, "uT", false);
    for (int b = 0; b <= nb+1; b++) { double v = u->GetBinContent(b); s[b] += v; s2[b] += v*v; }
    delete u; delete m;
  }
  double c = 0; int n = 0;
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin + ana::nPtBinsUsed; ipt++)
    for (int ixj = 0; ixj < nXj; ixj++) {
      int b = ipt*G + ixj + 2; double mean = s[b]/nToys;
      if (mean > 0) { c += (s2[b]/nToys - mean*mean)/(mean*mean); n++; }
    }
  return n ? c/n : 0;
}
void iter_test() {
  drawer d("pythia", "nominal");
  TH1D * reco = d.get(Form("hrecoxj%i",ir), 1);
  TH1D * truth = d.get(Form("htruthxj%i",ir), 1);
  TH2D * M = d.get2d(Form("hxjresponse%i",ir), 1);
  TH1D * asis = unfold_utility::buildFullyCorrected(d.get(Form("hrecoxj%i_0",ir), 0), d.get(Form("hrecoxj%i_2",ir), 0), "data", "nominal", ir, true);
  TH1D * clamped = (TH1D*)asis->Clone("clamped");
  for (int b = 0; b <= clamped->GetNbinsX()+1; b++) if (clamped->GetBinContent(b) < 0) clamped->SetBinContent(b, 0);
  RooUnfoldResponse * R = new RooUnfoldResponse(reco, truth, M);
  TH1D * fakesTr = (TH1D*)R->Hfakes()->Clone("fakesTr");
  auto fakesFor = [&](TH1D * meas) {
    TH1D * f = (TH1D*)meas->Clone(Form("fk_%s", meas->GetName())); f->Reset("ICES");
    for (int b = 1; b <= f->GetNbinsX(); b++) { double r = reco->GetBinContent(b); if (r > 0) f->SetBinContent(b, fakesTr->GetBinContent(b)/r*meas->GetBinContent(b)); }
    return f;
  };
  TH1D * fkA = fakesFor(asis), * fkC = fakesFor(clamped);
  printf("iter | refold chi2 asis  clamped | resp-toy chi2 asis  clamped (%d toys)\n", nToys);
  for (int it = 1; it <= 15; it++) {
    gRandom->SetSeed(1000 + it); double ta = respToyChi2(reco, truth, M, asis, it);
    gRandom->SetSeed(1000 + it); double tc = respToyChi2(reco, truth, M, clamped, it);
    printf("%4d | %10.4f %10.4f | %12.4g %12.4g\n", it, refoldChi2(R, asis, fkA, it), refoldChi2(R, clamped, fkC, it), ta, tc);
    fflush(stdout);
  }
}
