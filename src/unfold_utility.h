#ifndef UNFOLD_UTILITY_H
#define UNFOLD_UTILITY_H

#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphAsymmErrors.h"
#include "RooUnfoldResponse.h"

// Shared helpers for the drawing/*.C unfolding macros - previously duplicated
// (byte-for-byte, for unflattenXj/reflattenXj) across draw_purity_corrected.C,
// draw_iteration_halfclosure.C, toy_resp_iterations.C and toy_data_iterations.C.
class unfold_utility {
  public:
    // hrecoxj_abcd[ir][region] (written by unfolder.cc) is a single TH1D whose bin index
    // flattens (unfoldPtBin, xJ-bin) via ana::findUnfoldBin: bin = ipt*(nUnfoldXjBins+2)+ixj+1,
    // landing in ROOT bin (that value)+1. Pull out the nUnfoldXjBins real xJ bins for one
    // pT slice.
    static TH1D * unflattenXj(TH1D * flat, int ipt, const char * name);

    // Inverse of unflattenXj: write a per-pT-bin xJ histogram back into its slice of a
    // full flattened (unfoldPtBin, xJ) histogram, matching the same index convention.
    static void reflattenXj(TH1D * perPt, int ipt, TH1D * flatOut);

    // Build a RooUnfoldResponse from templates + a (possibly toyed) matrix and unfold
    // flatMeasured through it once. RooUnfoldResponse recomputes the Fakes histogram from
    // respRecoTemplate/respTruthTemplate vs. matrix's projection when they're non-empty,
    // but leaves the Truth/inefficiency accounting fixed to the templates - so toying only
    // `matrix` isolates the migration-probability fluctuation from also re-randomizing the
    // separate fakes/efficiency estimate.
    static TH1D * unfoldOnce(TH1D * respRecoTemplate, TH1D * respTruthTemplate, TH2D * matrix, TH1D * flatMeasured, int niter, const char * name);

    // Unfold flatMeasured through an already-built response.
    static TH1D * unfoldOnce(RooUnfoldResponse * response, TH1D * flatMeasured, int niter, const char * name);

    // Below this |P_A-P_C|, region C's own purity can't be disentangled from region
    // A's (the two regions have too similar a composition), and
    // purityCorrect/purityCorrectBkg/purityCorrectCoeffs fall back to treating region C
    // as 100% background - see the "Purity-Corrected Background Subtraction" derivation
    // in the analysis note / project conversation history.
    static constexpr float minPurityDiff = 0.03;

    // coeffA, coeffC such that the purity-corrected value of ANY linear quantity q
    // measured as qA in region A and qC in region C (raw per-bin counts, Sum(x_J),
    // Sum(x_J^2), or the regions' own totals NA/NC themselves) is coeffA*qA - coeffC*qC.
    // NA, NC are always the plain region A/C total counts, regardless of which q this is
    // applied to. This is the single source of truth for the two-purity formula's
    // coefficients - purityCorrect/purityCorrectBkg below both call this rather than
    // re-deriving it. Returns false (and the single-purity fallback coefficients,
    // coeffA=1, coeffC=(1-pA)*NA/NC) if NC<=0 or |P_A-P_C| < minPurityDiff.
    static bool purityCorrectCoeffs(float pA, float pC, float NA, float NC, float & coeffA, float & coeffC);

    // Exact two-purity purity-corrected xJ spectrum: solves
    //   A(xJ) = P_A*N_A*s(xJ) + (1-P_A)*N_A*bkg(xJ)
    //   C(xJ) = P_C*N_C*s(xJ) + (1-P_C)*N_C*bkg(xJ)
    // for the shared signal/background shapes s(xJ), bkg(xJ), returning P_A*N_A*s(xJ) -
    // i.e. region A background-subtracted using both regions' purity, not just region
    // A's. Falls back to the single-purity (region C = 100% background) formula if
    // P_A, P_C are too close to disentangle (minPurityDiff).
    // pA/pAErrLow/pAErrHigh, pC/pCErrLow/pCErrHigh are ana::getPurity/getPurityC's (and
    // their ErrorLow/ErrorHigh counterparts') asymmetric bootstrap values for this pT
    // bin. P_A and P_C uncertainty is propagated as two independent two-point sources,
    // sign-split per bin (see drawing/draw_systematics.C's asymmetricSystematics
    // convention) since the sign of d(signal)/dP_C can flip bin-to-bin - a fixed
    // "high variant -> up" mapping, valid for a single monotonic source, would be wrong
    // here. P_A and P_C are treated as uncorrelated - a simplification, not yet
    // cross-checked against their true (likely partially shared) bootstrap covariance.
    // TH1D bins can only hold one symmetric error, so h's bin errors use the larger of
    // errLow/errHigh (a conservative choice) - that's what feeds RooUnfold, Integral(),
    // Chi2 calcs, etc. downstream, none of which support asymmetric errors anyway. If
    // graphOut is non-null, *graphOut receives a TGraphAsymmErrors with the true
    // asymmetric errors, for display where the asymmetry should actually be visible.
    // Returns h=nullptr (and *graphOut=nullptr if requested) if region C has no
    // statistics in this pT bin (N_A/N_C undefined) - callers must check.
    static TH1D * purityCorrect(TH1D * A, TH1D * C, float pA, float pAErrLow, float pAErrHigh,
        float pC, float pCErrLow, float pCErrHigh, const char * name, TGraphAsymmErrors ** graphOut = nullptr);

    // The background piece subtracted off region A by purityCorrect above (Bkg = A -
    // signal), for display. Stat-only error propagation from A/C bin errors at fixed
    // P_A,P_C - the P_A/P_C uncertainty itself is already shown via purityCorrect's
    // asymmetric graph, not duplicated here since this curve is a secondary/illustrative
    // display only.
    static TH1D * purityCorrectBkg(TH1D * A, TH1D * C, float pA, float pC, const char * name);

    // ana::unfoldXjBins is non-uniform (0.1-wide up to xJ=1.3, then 0.2,0.2,0.3) -
    // plotting raw bin content directly makes an otherwise-smooth density look like it
    // has a shelf/cliff right where the bin width changes. Divide by bin width for
    // display only; histograms used for Integral()/purityCorrect etc. should stay raw.
    static TH1D * densityForDisplay(TH1D * h, const char * name);

    // Purity-correct all ana::nPtBins slices (not just the ana::nPtBinsUsed used for
    // physics results) of flatA/flatC via purityCorrect, and write the result into one
    // full flattened histogram matching the response matrix's dimensionality, so
    // RooUnfold sees a complete, consistently-binned "measured" vector - unfolding a 2D
    // (pT,xJ) measurement needs the whole flattened vector at once because migration
    // crosses pT-bin boundaries, not just xJ ones. If a pT slice's region C is empty
    // (can't purity-correct), falls back to raw region A for that slice only.
    static TH1D * buildFullyCorrected(TH1D * flatA, TH1D * flatC, const char * tag, string systag);
};

#endif // UNFOLD_UTILITY_H
