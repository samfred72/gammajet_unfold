#ifndef UNFOLD_UTILITY_H
#define UNFOLD_UTILITY_H

#include "ana.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphAsymmErrors.h"
#include "RooUnfoldResponse.h"

// Unfolding helpers shared by the drawing/*.C macros.
class unfold_utility {
  public:
    // hrecoxj-style histograms flatten (pT bin, x_J bin) as bin = ipt*(nUnfoldXjBins+2) + ixj + 1
    // (ana::findUnfoldBin). Extract one pT slice's nUnfoldXjBins x_J bins.
    static TH1D * unflattenXj(TH1D * flat, int ipt, const char * name);

    // Inverse of unflattenXj.
    static void reflattenXj(TH1D * perPt, int ipt, TH1D * flatOut);

    // Build a response from the templates and a (possibly toyed) matrix, then unfold once. Fakes are
    // recomputed from the templates vs the matrix projection; truth/efficiency stay with the templates,
    // so toying only `matrix` isolates the migration fluctuation.
    // includeSystematics adds RooUnfold's response-statistics covariance: ~3 s per call vs a few ms.
    // Toy loops read only bin contents and should pass false.
    static TH1D * unfoldOnce(TH1D * respRecoTemplate, TH1D * respTruthTemplate, TH2D * matrix, TH1D * flatMeasured, int niter, const char * name, bool includeSystematics = true);

    // Unfold through an already-built response (includeSystematics as above).
    static TH1D * unfoldOnce(RooUnfoldResponse * response, TH1D * flatMeasured, int niter, const char * name, bool includeSystematics = true);

    // Below this |P_A - P_C| the two regions cannot be separated; region C is then treated as pure background.
    static constexpr float minPurityDiff = 0.03;

    // Two-purity coefficients: the corrected value of any linear quantity measured as qA, qC in regions
    // A, C (counts, Sum x_J, Sum x_J^2, NA/NC) is coeffA*qA - coeffC*qC. NA, NC are the region totals.
    // Returns false with the single-purity fallback (coeffA = 1, coeffC = (1-pA)*NA/NC) if NC <= 0 or
    // |P_A - P_C| < minPurityDiff.
    static bool purityCorrectCoeffs(float pA, float pC, float NA, float NC, float & coeffA, float & coeffC);

    // Two-purity corrected x_J spectrum: solve A = P_A N_A s + (1-P_A) N_A b and
    // C = P_C N_C s + (1-P_C) N_C b for the shared shapes s, b and return P_A N_A s.
    // P_A and P_C errors are propagated as two independent two-point sources, sign-split per bin
    // (the sign of d(signal)/dP_C can change bin to bin; see draw_systematics.C). Bin errors take the
    // larger side; *graphOut, if given, gets the asymmetric errors for display. Returns nullptr if
    // region C is empty.
    // quiet suppresses the minPurityDiff fallback warning (expected every call for the 35-100 GeV
    // buffer bin, which has no purity point).
    static TH1D * purityCorrect(TH1D * A, TH1D * C, float pA, float pAErrLow, float pAErrHigh,
        float pC, float pCErrLow, float pCErrHigh, const char * name, TGraphAsymmErrors ** graphOut = nullptr,
        bool quiet = false);

    // Background subtracted from region A by purityCorrect (A - signal), for display; stat errors only.
    static TH1D * purityCorrectBkg(TH1D * A, TH1D * C, float pA, float pC, const char * name);

    // ana::unfoldXjBins is non-uniform: divide by bin width for display only.
    static TH1D * densityForDisplay(TH1D * h, const char * name);

    // Purity-correct every ana::nPtBins slice and reflatten into one histogram matching the response
    // matrix (unfolding needs the full vector: migration crosses pT bins). A slice with empty region C
    // falls back to raw region A.
    // ir selects the radius's purity curve (default R = 0.4; per-radius callers must pass it).
    // quietAll silences the fallback warning for every slice (repeated calls such as toy loops).
    // Negative purity-corrected contents in the buffer pT bins are set to 0 (clampBufferNegatives; see
    // the .cc). Reported bins are not changed.
    static constexpr bool clampBufferNegatives = true;
    static TH1D * buildFullyCorrected(TH1D * flatA, TH1D * flatC, const char * tag, string systag, int ir = 2, bool quietAll = false);
};

#endif // UNFOLD_UTILITY_H
