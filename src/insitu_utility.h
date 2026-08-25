#ifndef INSITU_UTILITY_H
#define INSITU_UTILITY_H

#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include <vector>
#include <string>
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"

// Shared helpers for the insitu/*.C in-situ JES calibration macros - previously
// duplicated byte-for-byte across grid_insitu.C, grid_insitu_jet12.C,
// grid_insitu_unfolded.C, and (for some) draw_insitu_xj.C.
struct DataEvent { float pho_pt, jet_pt; int ptbin; };

class insitu_utility {
  public:
    // In-situ JES grid scan window: every grid_insitu*.C in insitu/ scans a trial overall
    // jet-energy-scale factor pa over na=scanN steps of (scanHigh-scanLow)/scanN each,
    // covering [scanLow, scanHigh) - centralized here (was six copies of
    // "const float lowa = 0.95, higha = 1.05;"/"const int na = 1000;" duplicated across
    // grid_insitu.C, grid_insitu_jet12.C, grid_insitu_shapechi2.C, grid_insitu_unfolded.C,
    // and grid_insitu_unfolded_shapechi2.C) so the scan range only needs to change in one
    // place. scanN is chosen to hold the step size (scanHigh-scanLow)/scanN fixed at 1e-4
    // (the original [0.95,1.05]/1000 granularity) as the window moves - unfolder.cc's
    // ispairedInsitu floorScale is tied directly to scanLow (see fill_matrix()), so
    // events are always kept in insitutree down to whatever scanLow is set to here.
    static constexpr float scanLow = 0.90;
    static constexpr float scanHigh = 1.00;
    static constexpr int scanN = 1000;

    // Reads the insitutree (pho_pt, jet_pt, abcd, ir) written by unfolder.cc, keeping
    // only events in the requested ABCD region AND the requested jet radius (the tree
    // holds every ana::nJetR radius's rows together, one file per systag - see
    // unfolder.h's insitu_tree construction; a single underlying event can appear as up
    // to ana::nJetR separate rows, one per radius it paired at). If restrictToUsed
    // (default), only events in ana::ptBinsUsed are kept and ptbin is re-indexed to
    // 0..ana::nPtBinsUsed-1 (grid_insitu.C/grid_insitu_jet12.C's convention); if false,
    // every ana::ptBins bin is kept with ptbin = ana::findPtBin's raw index
    // (grid_insitu_unfolded.C's convention - needed so the full flattened (pT,xJ)
    // measured vector has a complete input for cross-pT-bin migration during unfolding).
    //
    // Deliberately does not read the tree's "weight" branch (the vz/cluster_pt mcWeight
    // - see unfolder.cc's in-situ test tree comment): every caller of this function
    // reads a Data_*_insitu.root file, where that weight is always 1.0 (mcWeight is
    // isMC-only), so it would be a no-op here. MC insitutree consumers must read and
    // apply "weight" themselves - see grid_insitu.C's referenceMeans/buildMCXjByPtBin.
    static vector<DataEvent> cacheDataEvents(const char * filename, int abcdSelect, int ir, bool restrictToUsed = true);

    // Low-xJ floor for a given jet radius `ir` and photon-pT bin's LOWER edge `ptLow` -
    // identical formula to unfolder::check_pair's `lowbin` (src/unfolder.cc): the
    // enclosing ana::unfoldXjBins bin's UPPER edge above jet_calib_pt_cut[ir]/ptLow (i.e.
    // rounded up to the next full xJ bin, same as check_pair, not down to the raw ratio
    // itself). Every insitu/grid_insitu*.C mean(x_J)/shape calculation should drop events
    // below this floor at whatever trial pa is being evaluated, so it excludes exactly
    // the same low-xJ events unfolder.cc's ispaired/check_pair would exclude from the
    // response matrix and hrecoxj/htruthxj at the same jet radius and JES scale - jets
    // below jet_calib_pt_cut[ir] aren't a reconstruction-trustworthy denominator for xJ,
    // regardless of which macro is computing a mean or a shape from them.
    static double lowXjFloor(int ir, double ptLow);

    // Canonical insitutree filename for a given trigger/sim/systag - one file holds
    // every jet radius (see cacheDataEvents' ir parameter to select one back out) -
    // mirrors unfolder.h's insitu_tree construction exactly, so every insitu/*.C macro
    // locates the same file unfolder.cc actually wrote instead of re-deriving the
    // naming convention independently in each macro. Pass sim="" for Data (3-part
    // name, no sim component, matching unfolder.cc's isMC branch); a non-empty sim for
    // MC (4-part name: trigger_sim_systag).
    static string insituFilename(const char * insitu_dir, const char * trigger,
        const char * sim, const string & systag);

    // Scans outward from the minimum on a chi2-vs-pa graph for the two points where
    // chi2 first crosses minchisq+1 (68% CL for one parameter).
    static void findError(TGraph * g, int ibest, float minchisq, float & errLow, float & errHigh);

    // Mean(x_J) vs. photon pT, one point per ana::ptBinsUsed bin (x error = half bin width).
    static TGraphErrors * meanGraph(const float mean[], const float err[], const char * name);

    // Ratio of two mean(x_J) arrays (e.g. Data/MC) vs. photon pT, errors combined
    // assuming the numerator and denominator are independent.
    static TGraphErrors * ratioGraph(const float meanNum[], const float errNum[],
        const float meanDen[], const float errDen[], const char * name);

    // sPHENIX label block: bold-italic "sPHENIX Internal" title, then one line per
    // sample, then one line per feature - same text/font convention as
    // drawer::drawAll() (src/drawer.cc), reimplemented here so the insitu/ macros don't
    // have to construct a full drawer (which opens a batch of unrelated unfolding-output
    // files they have no other use for).
    static void drawSPhenixLabel(vector<string> samples, vector<string> features,
        float drawx, float drawy, int fontsize, float csize);

    // ----- Computational helpers shared by the six grid_insitu*.C JES-scan macros -----
    // Consolidated here after the same logic was found copy-pasted (and once, actually
    // buggy in one of its copies but not the other) across grid_insitu.C,
    // grid_insitu_shapechi2.C, grid_insitu_unfolded.C,
    // grid_insitu_unfolded_shapechi2.C, grid_insitu_jet12.C, and
    // grid_insitu_jet12_shapechi2.C. See those macros for how each is used - none of
    // the logic below changed in the move, only its location.

    // Sums a fine-binned array's [startBin, nBinsForFit) range into groups of groupSize
    // fine bins each - used by the shape-chi2 macros' pT-bin-2 coarse rebinning (see
    // grid_insitu_shapechi2.C's coarseRebinPtBin comment).
    static vector<double> coarsenSum(const vector<double> & fine, int startBin, int nBinsForFit, int groupSize);
    // Same grouping, but combines per-bin errors in quadrature.
    static vector<double> coarsenQuadrature(const vector<double> & fineErr, int startBin, int nBinsForFit, int groupSize);

    // Weighted mean/error of x=jet_pt/pho_pt per used photon-pT bin, combining several
    // (filename,weight) MC samples read from their insitutree - the fixed reference a
    // grid_insitu*.C pa scan compares Data against. abcdSelect/ir/lowXj[] filter events
    // the same way cacheDataEvents does on the Data side. Sample-agnostic: pass
    // Photon5+10+20 for the primary in-situ study or a single Jet12(-family) sample for
    // the dijet-MC cross-check - see referenceShape() below for the shape-chi2 analogue.
    static void referenceMeans(const vector<pair<string,double>> & samples, int abcdSelect, int ir,
        float refMean[], float refMeanErr[], const float lowXj[]);

    // Fixed MC reference xJ SHAPE (bin fraction, not density) and its per-bin error, per
    // used photon-pT bin - the shape-chi2 analogue of referenceMeans() above (see
    // grid_insitu_shapechi2.C's referenceShape() comment for the bin-fraction-vs-density
    // and error-convention rationale).
    static void referenceShape(const vector<pair<string,double>> & samples, int abcdSelect, int ir,
        vector<vector<double>> & refFrac, vector<vector<double>> & refFracErr, const float lowXj[]);

    // Region-A-only mean(x_J)/error per used photon-pT bin, at a given trial
    // jet-energy-scale factor pa - display-only where the fit criterion is shape chi2 or
    // unfolded-vs-truth.
    static void computeRegionAMeans(const vector<DataEvent> & dataA, float pa, float mean[], float err[], const float lowXj[]);

    // Purity-corrected (two-purity method) mean(x_J)/error per used photon-pT bin, at a
    // given trial pa - region A and region C are scaled by the same pa, purity[]/
    // purityC[] held fixed (see unfold_utility::purityCorrectCoeffs).
    static void computeCorrectedMeans(const vector<DataEvent> & dataA, const vector<DataEvent> & dataC,
        float pa, const float purity[], const float purityC[], float mean[], float err[], const float lowXj[]);

    // x_J histogram per pT bin (ana::unfoldXjBins binning) from a cached Data sample, at
    // a given trial jet-energy-scale factor pa. nBins is ana::nPtBinsUsed for the
    // reco-level/Jet12-referenced macros, or ana::nPtBins for the unfolded macros (which
    // need every pT bin, including the migration-buffer/overflow bins, for cross-pT-bin
    // migration during unfolding - see grid_insitu_unfolded.C's cacheDataEvents comment).
    static vector<TH1D*> buildXjByPtBin(const vector<DataEvent> & data, float pa, int nBins,
        const char * prefix, const float lowXj[]);

    // x_J histogram per used photon-pT bin (ana::nPtBinsUsed) for the fixed (never
    // rescaled) cross-section-weighted MC reference - same samples/weights convention as
    // referenceMeans()/referenceShape() above. Always ana::nPtBinsUsed-sized: only the
    // reco-level/Jet12-referenced macros scan a raw MC ntuple like this at all: the
    // unfolded macros' truth reference comes from the response-matrix template instead.
    static vector<TH1D*> buildMCXjByPtBin(const vector<pair<string,double>> & samples, int abcdSelect, int ir,
        const char * prefix, const float lowXj[]);

    // Purity-correct region A/C histograms per pT bin via unfold_utility::purityCorrect
    // (src/unfold_utility.h) - the exact two-purity method. nBins is ana::nPtBinsUsed
    // (grid_insitu.C/grid_insitu_shapechi2.C, with real asymmetric purity-error arrays)
    // or ana::nPtBins (the unfolded macros, which pass zero-filled error arrays since
    // purity is held fixed there and the asymmetric term isn't needed).
    static vector<TH1D*> purityCorrectByPtBin(const vector<TH1D*> & hA, const vector<TH1D*> & hC, int nBins,
        const float purity[], const float purityErrLow[], const float purityErrHigh[],
        const float purityC[], const float purityCErrLow[], const float purityCErrHigh[],
        const char * prefix);
};

#endif // INSITU_UTILITY_H
