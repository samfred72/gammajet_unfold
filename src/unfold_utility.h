#ifndef UNFOLD_UTILITY_H
#define UNFOLD_UTILITY_H

#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "TH1D.h"
#include "TH2D.h"
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
};

#endif // UNFOLD_UTILITY_H
