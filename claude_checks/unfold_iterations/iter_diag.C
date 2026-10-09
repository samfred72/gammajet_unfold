#include "../../src/ana.h"
#include "../../src/drawer.h"
#include "../../src/unfold_utility.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Why do the refolding chi2 (note Fig. 17) and the response-toy chi2 (Fig. 19) jump after 10
// iterations? Nominal Data unfolding at R=0.4, iterations 1..15: negative / non-finite bins in the
// measured input and in each unfolded result, and the bins that move most between iterations.
const int ir = 2;
const int G = ana::nUnfoldXjBins + 2;
string binName(int b) { // TH1 bin -> (pT bin, xJ bin)
  int g = b - 1, ipt = g / G, ixj = g % G - 1;
  return Form("b%d[pt%d %s]", b, ipt, ixj < 0 ? "uf" : ixj >= ana::nUnfoldXjBins ? "of" : Form("xj%.1f-%.1f", ana::unfoldXjBins[ixj], ana::unfoldXjBins[ixj+1]));
}
void iter_diag() {
  drawer d("pythia", "nominal");
  TH1D * reco = d.get(Form("hrecoxj%i",ir), 1);
  TH1D * truth = d.get(Form("htruthxj%i",ir), 1);
  TH2D * M = d.get2d(Form("hxjresponse%i",ir), 1);
  TH1D * meas = unfold_utility::buildFullyCorrected(d.get(Form("hrecoxj%i_0",ir), 0), d.get(Form("hrecoxj%i_2",ir), 0), "data", "nominal", ir, true);
  printf("measured: %d bins, integral %.1f\n", meas->GetNbinsX(), meas->Integral());
  for (int b = 0; b <= meas->GetNbinsX()+1; b++) if (meas->GetBinContent(b) < 0) printf("  NEGATIVE measured %s = %.3f +- %.3f\n", binName(b).c_str(), meas->GetBinContent(b), meas->GetBinError(b));
  // truth bins with no response (prior 0)
  TH1D * prev = nullptr;
  for (int it = 1; it <= 15; it++) {
    TH1D * u = unfold_utility::unfoldOnce(reco, truth, M, meas, it, Form("u%d", it), false);
    int nneg = 0, nnan = 0; double vmin = 1e30; int bmin = -1;
    for (int b = 1; b <= u->GetNbinsX(); b++) { double v = u->GetBinContent(b); if (!std::isfinite(v)) nnan++; else { if (v < 0) nneg++; if (v < vmin) { vmin = v; bmin = b; } } }
    printf("iter %2d: integral %.1f  neg %d  nonfinite %d  min %s=%.4g", it, u->Integral(), nneg, nnan, binName(bmin).c_str(), vmin);
    if (prev) { // largest |change| relative to previous
      double worst = 0; int bw = -1;
      for (int b = 1; b <= u->GetNbinsX(); b++) { double p = prev->GetBinContent(b), v = u->GetBinContent(b); if (p == 0 && v == 0) continue; double r = fabs(v - p)/std::max(fabs(p), 1e-9); if (r > worst) { worst = r; bw = b; } }
      printf("  | max rel change %s: %.4g -> %.4g", binName(bw).c_str(), prev->GetBinContent(bw), u->GetBinContent(bw));
    }
    printf("\n");
    delete prev; prev = u;
  }
}
