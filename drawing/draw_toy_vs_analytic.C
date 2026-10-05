#include "../src/ana.h"
#include "../src/drawer.h"
#include "../src/unfold_utility.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include <map>
#include <iomanip>
#include <algorithm>
// The original gammajet project builds its OWN, differently-laid-out drawer/ana classes
// into /home/samson72/root/lib/libgammajet.so, sitting on the same library search path
// as this project's libgammajet_unfold.so. Without forcing which one loads first, ROOT's
// implicit symbol autoload can bind drawer/ana calls to the wrong (mismatched-layout)
// library and segfault - explicit load removes the ambiguity.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Direct bin-by-bin comparison of RooUnfoldBayes's analytic statistical error
// (RooUnfoldBayes::Hreco(), "kErrors" - the diagonal of the covariance the Bayes
// algorithm itself computes, RooUnfold/src/RooUnfoldBayes.cxx's getCovariance()) against
// an ATLAS-style Poisson-toy bootstrap (Phys. Lett. B 774 (2017) 379, Sec. 6):
// "Stochastic variations of the data are generated based on its statistical
// uncertainty and each variation is unfolded and projected into xJ. The statistical
// covariance of the set is taken as the statistical uncertainty." / "An additional
// covariance is obtained from applying the pseudo-experiment procedure to the response
// matrix."
//
// This repo already implements both halves of that toy procedure separately -
// toy_iterations.C ("data") (Data-side Poisson toys) and toy_iterations.C ("resp")
// (response-matrix-side Poisson toys) - but each only compares its OWN toy ensemble's
// niter-convergence behavior to the "iter n vs n-1" pairwise metric; neither puts its
// toy-derived spread side by side against the analytic error RooUnfoldBayes actually
// reports (which is what ends up in draw_final_result.C's error bars - see that file's
// header comment). That direct comparison is what this macro adds: rerun both toy
// ensembles at the exact same fixed working point draw_final_result.C uses (nominal
// systag, R=0.4, niter=2), combine their bin-by-bin spread in quadrature (matching how
// RooUnfoldBayes's own getCovariance() sums its Data-statistics and response-matrix-
// statistics terms), and plot that against the analytic error directly.
//
// ---------------------------------------------------------------------------------
// Weighted-Poisson toy generator for the response matrix
// ---------------------------------------------------------------------------------
// toy_iterations.C ("resp")'s own poissonToyMatrix draws straight from Poisson(bin content)
// - correct ONLY if bin content is a raw, unweighted count. The response matrix
// (hxjresponse<ir>) isn't: it's a cross-section-weighted combination of Photon5+10+20
// (unfolder.h's TH1::SetDefaultSumw2() means every bin's GetBinError() already correctly
// reflects that - sqrt(sum of w_i^2), not sqrt(content)).
// Measured directly (see the git history for this file / the session that added it):
// some response-matrix bins have GetBinContent() in the hundreds of thousands but
// GetBinError() implying as few as ~1-3 EFFECTIVE independent entries - a single
// high-cross-section Photon5 event can outweigh hundreds of Photon20 events. Drawing
// Poisson(content) for a bin like that (content ~300000) gives a ~0.2% relative
// fluctuation; the bin's true effective statistics give more like 50-100%. This is why
// toy_iterations.C ("resp")'s own response-matrix-toy contribution came out negligible
// (orders of magnitude below the Data-toy term) in the quick check that motivated this
// macro - it was underestimating the true response-matrix statistical noise, not
// correctly finding it small.
//
// Fix: draw each bin's toy from Poisson(N_eff), where N_eff = content^2/error^2 is the
// standard effective-entries (Kish) approximation for a weighted sum, then rescale by
// content/N_eff (the mean weight per effective entry) so the toy's mean and variance
// still match (content, error^2) exactly. This reduces to ordinary Poisson(content)
// whenever error == sqrt(content) (a genuinely unweighted bin, N_eff == content), but
// correctly reflects the true statistics otherwise.
//
// ---------------------------------------------------------------------------------
// Plain-Poisson toy generator for the Data side - NOT weighted-Poisson
// ---------------------------------------------------------------------------------
// An earlier version of this macro also weighted-Poisson-toyed the purity-CORRECTED
// measured spectrum (unfold_utility::buildFullyCorrected's output) directly, reasoning
// that it's a signed combination (coeffA*A - coeffC*C, see unfold_utility::purityCorrect)
// rather than a raw count, so its content/error^2 wouldn't reduce to a real entry count
// either. That reasoning doesn't hold: the Kish N_eff approximation models a
// POSITIVE-WEIGHTED SUM (many events of unequal weight, like the response matrix above)
// - it is not a valid model for a SIGNED DIFFERENCE of two raw, unweighted counts, which
// is what the purity correction actually is (region A and region C are both genuine,
// unweighted data counts; only the *combination* coeffA*A-coeffC*C is a non-trivial
// linear form). When A and C nearly cancel (small content, comparably-sized or larger
// error - measured to happen in the sparse, high-xJ tail of the highest used pT bin,
// e.g. content=2.6 +/- 3.06 there - see the session that added this comment), N_eff =
// content^2/error^2 collapses below 1: the toy then draws from a coarse Poisson(<1) times
// a huge per-quantum weight (content/N_eff) instead of a smooth fluctuation around
// (content, error^2), wildly overestimating that bin's spread relative to the correctly-
// propagated analytic error - and since unfolding mixes bins via migration, that
// contaminates the unfolded RESULT bins' toy spread too, well beyond the true
// uncertainty. This showed up as toy/analytic disagreement of 5-15x specifically (and
// only) in that pT bin, growing bin-by-bin with xJ, while the response-side toy
// contribution there stayed small and unremarkable - i.e. the DATA-side toy method
// itself was at fault, not a missing response-matrix-statistics term.
//
// Fix: toy the RAW region-A/region-C counts with plain, ordinary Poisson (no Kish
// approximation needed at all - they're genuinely unweighted real-data counts).
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

// ---------------------------------------------------------------------------------
// Fixed-coefficient purity combination for toys - NOT unfold_utility::buildFullyCorrected
// ---------------------------------------------------------------------------------
// An earlier version of this fix ran each toy's raw region-A/region-C draw through
// unfold_utility::buildFullyCorrected (i.e. the real, production purityCorrect path)
// directly. That's wrong for a like-for-like toy/analytic comparison: purityCorrect
// derives its coefficients (coeffA, coeffC - see unfold_utility::purityCorrectCoeffs)
// from NA=A->Integral(), NC=C->Integral() - the pT slice's TOTAL region A/C counts
// (summed over all xJ bins) - and its error formula propagates ae/ce (the per-bin
// counts) treating coeffA/coeffC/K=NA/NC as FIXED CONSTANTS (see the errLow/errHigh
// derivation in purityCorrect - no term for d(coeffC)/d(a) or d(coeffC)/d(c) through
// NA/NC). Re-deriving coeffA/coeffC fresh from each toy's OWN fluctuated NA/NC (as
// buildFullyCorrected does) adds a real extra variance term the analytic formula never
// accounts for; a second earlier version instead randomly re-drew pA/pC themselves per
// toy (jointly, independently, from a split-normal around their own asymmetric bootstrap
// error) to capture purity-measurement uncertainty - also wrong, but for a different
// reason: purityCorrect's error formula does NOT treat pAErrLow/pAErrHigh/pCErrLow/
// pCErrHigh as a distribution to sample from - it evaluates the deviation from shifting
// ONE parameter at a time to its edge (devPAHigh/devPALow/devPCHigh/devPCLow) and
// combines same-signed ones in quadrature (CLAUDE.md's documented per-source sign-split
// convention, same as JER/emscale/JES elsewhere). A joint, continuous, per-toy random
// draw of BOTH pA and pC is a strictly wider calculation than that fixed two-point
// envelope, and gets amplified further wherever content is already a near-cancellation
// (a-c) - it showed up as toy/analytic agreement worsening everywhere (~1.5x even in
// well-populated low-pT bins) and blowing up worse than before at high pT/high xJ.
//
// Fix: hold coeffA/coeffC fixed per pT bin (computed once from the nominal, non-toy
// flatA/flatC - correct for the raw-counting part, see above) and handle the SEPARATE
// purity-measurement-uncertainty contribution deterministically instead of by toying -
// see purityErrSqByFlatBin below, added onto each toy's content as independent Gaussian
// noise of the SAME (precomputed, analytic-matching) magnitude rather than by randomly
// varying pA/pC.
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
      A->SetBinError(i, sqrt(fabs(content))); // only the toy's central value is read downstream - see unfoldOnce
    }
    unfold_utility::reflattenXj(A, ipt, flatCorrected);
    delete A; delete C;
  }
  return flatCorrected;
}

const int ir = 2; // nominal jet radius index (R=0.4) - same fixed working point as toy_iterations.C
const int nPtBinsUsed = ana::nPtBinsUsed;
const int nToys = 10000; // toy count per side (Data, response)
const int niterate = 2; // matches draw_final_result.C's actual working point
const int nXjBinsForChi2 = ana::nUnfoldXjBins - 3; // exclude the low-stat tail, same convention as elsewhere

// Set false to skip the response-side toy loop entirely and report data-side-only toy
// uncertainty. Currently true (both toy and analytic - see unfold_utility.cc's
// unfoldOnce - include the response-matrix-statistics term, for a like-for-like
// comparison; CLAUDE.md's "Local RooUnfold Patch" section has the RooUnfold-side bug
// that used to block enabling this analytically at all).
const bool includeResponseUncertainty = true;

// Effective-entries (Kish) weighted-Poisson draw for one bin - see the file-header
// comment above for why this replaces plain Poisson(content).
double weightedPoissonToyValue(double content, double error) {
  if (error <= 0) return content; // no statistical info to fluctuate on (e.g. a never-filled bin) - leave unchanged
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
  // Every toy histogram below gets a unique name (hUnfoldToyD_<i>, hUnfoldToyR_<i>, ...)
  // - without this, each one auto-registers into gDirectory and accumulates there for
  // the life of the process (never reclaimed just by our own `delete`, since that only
  // frees the object, not gDirectory's own bookkeeping list of it), which starts to
  // visibly slow down every subsequent New()/Clone() call once the count reaches into
  // the thousands - the toy_*_iterations.C precedent this macro is based on has the
  // same gap but was only ever run at 1000 toys, where it doesn't yet show up.
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

  // Fixed purity-correction coefficients per pT bin, from the NOMINAL (non-toy) flatA/
  // flatC - see toyMeasuredFixedCoeffs's header comment for why these must be held fixed
  // across toys rather than re-derived from each toy's own fluctuated counts.
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

  // Deterministic per-bin purity-measurement-uncertainty variance, to be injected as
  // independent Gaussian noise onto each toy's measured content (see
  // toyMeasuredFixedCoeffs's header comment for why NOT to instead randomly toy pA/pC).
  // flatMeasured's own GetBinError() (from unfold_utility::purityCorrect, at the top of
  // this function) already equals sqrt(statErr^2 + purityPiece^2) - the same
  // devPA/devPC-based quadrature term this macro would otherwise have to re-derive by
  // hand is recoverable by simply subtracting off statErr^2 (the pure counting-error
  // term, computed independently here from ae/ce and the SAME fixed coeffA/coeffC above)
  // in quadrature. This exactly reproduces the analytic formula's own purity contribution
  // - by construction, not approximation - so that piece of the comparison can never
  // disagree; only the genuinely toy-simulated counting statistics can.
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

  // Full per-toy value storage (not just running sum/sumsq) - needed for the percentile-
  // based error below, which isn't derivable from sum/sumsq alone. nFlatBins*nToys
  // doubles is a few MB even at nToys=10000 - cheap.
  vector<vector<double>> toyValD(nFlatBins+2, vector<double>(nToys));
  vector<vector<double>> toyValR(nFlatBins+2, vector<double>(nToys));

  // ---- Data-side toys (plain-Poisson-fluctuate raw region A/C, response fixed) ----
  // See the file-header comment above: toying flatMeasured directly (the already-
  // purity-corrected combination) with the weighted-Poisson/Kish generator was wrong -
  // toy the raw counts feeding into it instead (toyMeasuredFixedCoeffs), then add the
  // purity-measurement-uncertainty contribution deterministically (purityErrSqByFlatBin
  // above), rather than randomly toying pA/pC.
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
    // includeSystematics=false: only GetBinContent() is ever read below, so skip paying
    // for the (expensive) response-matrix-statistics covariance on every toy - see
    // unfold_utility.h's includeSystematics comment.
    TH1D * hToy = unfold_utility::unfoldOnce(respRecoTemplate, respTruthTemplate, respMatrix2D, toyData, niterate, Form("hUnfoldToyD_%d", itoy), false);
    for (int b = 0; b <= nFlatBins+1; b++) { double v = hToy->GetBinContent(b); sumD[b] += v; sumsqD[b] += v*v; toyValD[b][itoy] = v; }
    delete hToy; delete toyData; delete toyFlatA; delete toyFlatC;
    if (itoy % 200 == 0) cout << "  data toy " << itoy << "/" << nToys << endl;
  }

  // ---- Response-side toys (weighted-Poisson-fluctuate respMatrix2D, Data fixed) ----
  // Skipped entirely when includeResponseUncertainty is false (TEMPORARY, see that
  // const's comment) - sumR/sumsqR stay zero (correctly giving respRMS=0 below), but
  // toyValR stays all-zero too, which would NOT correctly give respLo/respHi=0 in
  // percentileInterval (it'd read the ensemble's own [0,0] spread relative to a nonzero
  // `content` as a huge fake interval) - guarded explicitly at the percentile call below.
  vector<double> sumR(nFlatBins+2, 0), sumsqR(nFlatBins+2, 0);
  if (includeResponseUncertainty) {
    cout << "Running " << nToys << " response-matrix weighted-Poisson toys at niter=" << niterate << "..." << endl;
    for (int itoy = 0; itoy < nToys; itoy++) {
      TH2D * toyMat = weightedPoissonToyMatrix(respMatrix2D, Form("hToyMat_%d", itoy));
      // includeSystematics=false - see the data-side toy loop's comment above.
      TH1D * hToy = unfold_utility::unfoldOnce(respRecoTemplate, respTruthTemplate, toyMat, flatMeasured, niterate, Form("hUnfoldToyR_%d", itoy), false);
      for (int b = 0; b <= nFlatBins+1; b++) { double v = hToy->GetBinContent(b); sumR[b] += v; sumsqR[b] += v*v; toyValR[b][itoy] = v; }
      delete hToy; delete toyMat;
      if (itoy % 200 == 0) cout << "  resp toy " << itoy << "/" << nToys << endl;
    }
  } else {
    cout << "Skipping response-matrix toys (includeResponseUncertainty=false) - data-side uncertainty only." << endl;
  }

  // 16th/84th-percentile interval of a toy ensemble - a ROBUST alternative to mean+-RMS,
  // matching how puritymaker.C's own bootstrap error is already defined in this project
  // (ana.cc's getPurityErrorLow/High comment: "asymmetric (16th/84th percentile around
  // the median)"). Standard deviation isn't robust to the rare-but-real outlier toy
  // draws that show up once the underlying raw counts are this low (0-3 events, the
  // highest-pT bin's highest-xJ tail) - a handful of such draws dominate the sum-of-
  // squares RMS while barely reflecting where the bulk of the ensemble actually sits.
  // `center` is the NOMINAL (non-toy) content, not the toy ensemble's own median/mean -
  // so this is a deviation from the same reference point the analytic error uses,
  // directly comparable to it. lo/hi are returned as positive magnitudes.
  auto percentileInterval = [](vector<double> vals, double center, double & lo, double & hi) {
    std::sort(vals.begin(), vals.end());
    int n = (int)vals.size();
    int i16 = std::max(0, (int)std::lround(0.16*(n-1)));
    int i84 = std::min(n-1, (int)std::lround(0.84*(n-1)));
    lo = std::max(center - vals[i16], 0.0);
    hi = std::max(vals[i84] - center, 0.0);
  };

  // Build a "toy-error" version of hNominal (same central values, toy-combined errors)
  // and collect per-bin ratios for the summary page. Percentile-based low/high are kept
  // in parallel arrays (indexed by flatbin) for the output file and the ratio-vs-analytic
  // summary below - RMS-based and percentile-based numbers are reported side by side
  // rather than one replacing the other, so they can be compared directly.
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
    // Loop over every xJ bin (not just the nXjBinsForChi2 chi2-summary subset) so
    // hNominalToyErr - which the per-pT spectrum pages below plot in full, low-stat
    // tail included - gets a real toy-derived error in every bin. Without this, the
    // last 3 (excluded-from-chi2) bins never had SetBinError() called on them at all
    // and silently kept hNominal->Clone()'s original analytic error instead, showing up
    // as an unphysical cliff (huge toy error, then a sudden drop back to the much
    // smaller analytic value) right at the chi2-cutoff boundary on those plots.
    for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) {
      int flatbin = ipt*(ana::nUnfoldXjBins+2) + ixj + 1 + 1;
      double meanD = sumD[flatbin]/nToys, varD = sumsqD[flatbin]/nToys - meanD*meanD;
      double meanR = sumR[flatbin]/nToys, varR = sumsqR[flatbin]/nToys - meanR*meanR;
      double dataRMS = sqrt(std::max(varD,0.0));
      double respRMS = sqrt(std::max(varR,0.0));
      double combinedToy = sqrt(dataRMS*dataRMS + respRMS*respRMS);
      hNominalToyErr->SetBinError(flatbin, combinedToy);

      double content = hNominal->GetBinContent(flatbin);
      // Percentile-based interval, combining data+response in quadrature per side
      // (same per-source sign-split convention as CLAUDE.md's asymmetric-systematics
      // handling elsewhere in this project - up-with-up, down-with-down).
      double dataLo, dataHi, respLo = 0, respHi = 0;
      percentileInterval(toyValD[flatbin], content, dataLo, dataHi);
      if (includeResponseUncertainty) percentileInterval(toyValR[flatbin], content, respLo, respHi);
      double percLo = sqrt(dataLo*dataLo + respLo*respLo);
      double percHi = sqrt(dataHi*dataHi + respHi*respHi);
      percLoByFlatBin[flatbin] = percLo;
      percHiByFlatBin[flatbin] = percHi;

      // The ratio/chi2-style summary itself still excludes the low-stat tail, per the
      // same convention as draw_covariance_chi2.C:60.
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

  // ---- Page 1: summary - toy/analytic ratio across every used bin, colored by pT bin ----
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
    // Open-circle series (matching draw_final_result.C's hTruthDisp convention for a
    // secondary overlay) - same per-pT-bin color as the filled RMS-based point, so the
    // two markers for one bin are visually paired.
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

  // ---- Pages 2..: per-pT-bin unfolded spectrum with analytic vs toy error bands ----
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hAnaPt = unfold_utility::unflattenXj(hNominal, ipt, Form("hAnaPt_%d", ipt));
    TH1D * hToyPt = unfold_utility::unflattenXj(hNominalToyErr, ipt, Form("hToyPt_%d", ipt));
    TH1D * hAnaDisp = unfold_utility::densityForDisplay(hAnaPt, Form("hAnaDisp_%d", ipt));
    TH1D * hToyDisp = unfold_utility::densityForDisplay(hToyPt, Form("hToyDisp_%d", ipt));

    // Percentile-based interval for this pT bin, density-scaled (divided by bin width)
    // the same way densityForDisplay treats hAnaDisp/hToyDisp above, so all three sit on
    // the same y-axis. A TGraphAsymmErrors, not a TH1D, since the interval is genuinely
    // asymmetric - same box-drawing technique as draw_final_result.C's gSystBox (per-
    // point box spanning the bin width in x, the asymmetric error in y), just green
    // instead of azure, plus open-circle markers at the central points (matching how
    // draw_final_result.C itself distinguishes an overlay series - see hTruthDisp's
    // marker style 24 there) rather than reusing the filled circles already used for
    // hAnaDisp's own data points.
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

  // Percentile-based (16th/84th) toy interval, as a TGraphAsymmErrors alongside the
  // RMS-based hNominalToyErr histogram - matching how this project already stores other
  // asymmetric errors (e.g. purityCorrect's own graphOut) rather than collapsing to one
  // symmetric number. Covers every xJ bin (not just the chi2-used subset), same as
  // hNominalToyErr.
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
