#include "../../src/ana.h"
#include "../../src/drawer.h"
#include "../../src/unfold_utility.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Verifies unfold_utility::clampBufferNegatives: per radius, the refolding chi2/NDF vs iteration
// (draw_refolding.C's metric) with the new input, and the shift of the reported, shape-normalized
// unfolded result at niter=2 against the unclamped input (rebuilt here with purityCorrect).
const int G = ana::nUnfoldXjBins + 2;
const int nXj = ana::nUnfoldXjBins - 3;

TH1D * unclamped(TH1D * flatA, TH1D * flatC, TH1D * clampedIn, int ir) {
  TH1D * h = (TH1D*)clampedIn->Clone("unclamped");
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
    bool isBuffer = ipt < ana::firstUsedPtBin || ipt >= ana::firstUsedPtBin + ana::nPtBinsUsed;
    if (!isBuffer) continue;
    float lo = ana::ptBins[ipt], hi = ana::ptBins[ipt+1];
    TH1D * A = unfold_utility::unflattenXj(flatA, ipt, "uA"); TH1D * C = unfold_utility::unflattenXj(flatC, ipt, "uC");
    TH1D * hc = unfold_utility::purityCorrect(A, C, ana::getPurity(lo, hi, "nominal", ir), ana::getPurityErrorLow(lo, hi, "nominal", ir),
        ana::getPurityErrorHigh(lo, hi, "nominal", ir), ana::getPurityC(lo, hi, "nominal", ir), ana::getPurityCErrorLow(lo, hi, "nominal", ir),
        ana::getPurityCErrorHigh(lo, hi, "nominal", ir), "hc", nullptr, true);
    unfold_utility::reflattenXj(hc ? hc : A, ipt, h);
    delete A; delete C; delete hc;
  }
  return h;
}
double refold(RooUnfoldResponse * R, TH1D * reco, TH1D * fakesTr, TH1D * meas, int it) {
  TH1D * fk = (TH1D*)meas->Clone("fk"); fk->Reset("ICES");
  for (int b = 1; b <= fk->GetNbinsX(); b++) { double r = reco->GetBinContent(b); if (r > 0) fk->SetBinContent(b, fakesTr->GetBinContent(b)/r*meas->GetBinContent(b)); }
  TH1D * u = unfold_utility::unfoldOnce(R, meas, it, "uR", false);
  TH1D * f = (TH1D*)R->ApplyToTruth(u, "fR"); f->Add(fk);
  double c = 0; int n = 0;
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin + ana::nPtBinsUsed; ipt++)
    for (int ixj = 0; ixj < nXj; ixj++) { int b = ipt*G + ixj + 2; double e = meas->GetBinError(b); if (e <= 0) continue; c += pow(f->GetBinContent(b) - meas->GetBinContent(b), 2)/(e*e); n++; }
  delete u; delete f; delete fk;
  return n ? c/n : 0;
}
void verify_clamp() {
  drawer d("pythia", "nominal");
  for (int ir = 0; ir < ana::nJetR; ir++) {
    TH1D * reco = d.get(Form("hrecoxj%i",ir), 1), * truth = d.get(Form("htruthxj%i",ir), 1);
    TH2D * M = d.get2d(Form("hxjresponse%i",ir), 1);
    TH1D * fA = d.get(Form("hrecoxj%i_0",ir), 0), * fC = d.get(Form("hrecoxj%i_2",ir), 0);
    TH1D * now = unfold_utility::buildFullyCorrected(fA, fC, Form("v%d", ir), "nominal", ir, true);
    TH1D * old = unclamped(fA, fC, now, ir);
    int nclamped = 0; for (int b = 1; b <= now->GetNbinsX(); b++) if (now->GetBinContent(b) != old->GetBinContent(b)) nclamped++;
    RooUnfoldResponse * R = new RooUnfoldResponse(reco, truth, M);
    TH1D * fakesTr = (TH1D*)R->Hfakes()->Clone("fakesTr");
    printf("R=%.1f: %d buffer bins clamped\n  refold chi2/NDF iter 1..15 new:", ana::JetRs[ir], nclamped);
    for (int it = 1; it <= 15; it++) printf(" %.3f", refold(R, reco, fakesTr, now, it));
    printf("\n  refold chi2/NDF iter 1..15 old:");
    for (int it = 1; it <= 15; it++) printf(" %.3f", refold(R, reco, fakesTr, old, it));
    TH1D * uN = unfold_utility::unfoldOnce(R, now, 2, "uN", false), * uO = unfold_utility::unfoldOnce(R, old, 2, "uO", false);
    double worstRel = 0, worstPull = 0; string where;
    for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin + ana::nPtBinsUsed; ipt++) {
      double nN = 0, nO = 0;
      for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) { int b = ipt*G + ixj + 2; nN += uN->GetBinContent(b); nO += uO->GetBinContent(b); }
      for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) {
        int b = ipt*G + ixj + 2; double o = uO->GetBinContent(b)/nO, n = uN->GetBinContent(b)/nN, e = uO->GetBinError(b)/nO;
        double pull = e > 0 ? (n - o)/e : 0, rel = o > 0 ? (n - o)/o : 0;
        if (fabs(pull) > fabs(worstPull)) { worstPull = pull; worstRel = rel; where = Form("pT %.0f-%.0f xJ %.1f-%.1f", ana::ptBins[ipt], ana::ptBins[ipt+1], ana::unfoldXjBins[ixj], ana::unfoldXjBins[ixj+1]); }
      }
    }
    printf("\n  niter=2 reported shape: largest shift %+.3f sigma (%+.2f%%) at %s\n", worstPull, 100*worstRel, where.c_str());
    delete uN; delete uO; delete R;
  }
}
