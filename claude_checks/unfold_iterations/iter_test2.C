#include "../../src/ana.h"
#include "../../src/drawer.h"
#include "../../src/unfold_utility.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Follow-up to iter_test.C: (1) clamping only the buffer pT bin (35-100 GeV) vs only the reported
// bins, refold chi2 vs iteration; (2) effect of clamping all negative bins on the reported unfolded
// result at the nominal 2 iterations, relative to its statistical error.
const int ir = 2;
const int G = ana::nUnfoldXjBins + 2;
const int nXj = ana::nUnfoldXjBins - 3;
double refoldChi2(RooUnfoldResponse * R, TH1D * reco, TH1D * fakesTr, TH1D * meas, int it) {
  TH1D * fk = (TH1D*)meas->Clone("fk"); fk->Reset("ICES");
  for (int b = 1; b <= fk->GetNbinsX(); b++) { double r = reco->GetBinContent(b); if (r > 0) fk->SetBinContent(b, fakesTr->GetBinContent(b)/r*meas->GetBinContent(b)); }
  TH1D * u = unfold_utility::unfoldOnce(R, meas, it, "uR", false);
  TH1D * f = (TH1D*)R->ApplyToTruth(u, "fR"); f->Add(fk);
  double c = 0; int n = 0;
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin + ana::nPtBinsUsed; ipt++)
    for (int ixj = 0; ixj < nXj; ixj++) {
      int b = ipt*G + ixj + 2; double e = meas->GetBinError(b); if (e <= 0) continue;
      c += pow(f->GetBinContent(b) - meas->GetBinContent(b), 2)/(e*e); n++;
    }
  delete u; delete f; delete fk;
  return n ? c/n : 0;
}
void iter_test2() {
  drawer d("pythia", "nominal");
  TH1D * reco = d.get(Form("hrecoxj%i",ir), 1);
  TH1D * truth = d.get(Form("htruthxj%i",ir), 1);
  TH2D * M = d.get2d(Form("hxjresponse%i",ir), 1);
  TH1D * asis = unfold_utility::buildFullyCorrected(d.get(Form("hrecoxj%i_0",ir), 0), d.get(Form("hrecoxj%i_2",ir), 0), "data", "nominal", ir, true);
  auto clampIf = [&](const char * name, auto sel) {
    TH1D * h = (TH1D*)asis->Clone(name);
    for (int b = 0; b <= h->GetNbinsX()+1; b++) if (h->GetBinContent(b) < 0 && sel((b-1)/G)) h->SetBinContent(b, 0);
    return h;
  };
  int buf = ana::nPtBins - 1;
  TH1D * cBuf = clampIf("cBuf", [&](int ipt){ return ipt == buf; });
  TH1D * cRep = clampIf("cRep", [&](int ipt){ return ipt != buf; });
  TH1D * cAll = clampIf("cAll", [&](int){ return true; });
  RooUnfoldResponse * R = new RooUnfoldResponse(reco, truth, M);
  TH1D * fakesTr = (TH1D*)R->Hfakes()->Clone("fakesTr");
  printf("iter | refold chi2: asis  buffer-only-clamped  reported-only-clamped\n");
  for (int it = 1; it <= 15; it++)
    printf("%4d | %8.4f %8.4f %8.4f\n", it, refoldChi2(R, reco, fakesTr, asis, it), refoldChi2(R, reco, fakesTr, cBuf, it), refoldChi2(R, reco, fakesTr, cRep, it));
  // reported bins at niter=2
  TH1D * uA = unfold_utility::unfoldOnce(R, asis, 2, "uA", false);
  TH1D * uC = unfold_utility::unfoldOnce(R, cAll, 2, "uC", false);
  printf("\nniter=2, reported bins: (clamped - asis)/stat.err of asis, and relative shift\n");
  double worstPull = 0, worstRel = 0;
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin + ana::nPtBinsUsed; ipt++) {
    double nA = 0, nC = 0;
    for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) { int b = ipt*G + ixj + 2; nA += uA->GetBinContent(b); nC += uC->GetBinContent(b); }
    for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) {
      int b = ipt*G + ixj + 2; double a = uA->GetBinContent(b), c = uC->GetBinContent(b), e = uA->GetBinError(b);
      // shape-normalized, as the result is reported
      double rel = (a > 0) ? (c/nC - a/nA)/(a/nA) : 0, pull = e > 0 ? (c/nC - a/nA)/(e/nA) : 0;
      if (fabs(pull) > fabs(worstPull)) worstPull = pull;
      if (fabs(rel) > fabs(worstRel)) worstRel = rel;
      if (fabs(rel) > 0.005) printf("  pt%d xj%.1f-%.1f: rel %+.4f  pull %+.3f\n", ipt, ana::unfoldXjBins[ixj], ana::unfoldXjBins[ixj+1], rel, pull);
    }
  }
  printf("worst |rel| %.4f, worst |pull| %.3f\n", fabs(worstRel), fabs(worstPull));
}
