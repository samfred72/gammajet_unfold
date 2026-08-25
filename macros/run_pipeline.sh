#!/bin/bash
# run_pipeline.sh - runs the full gammajet_unfold analysis chain top to bottom.
#
# Usage:
#   bash run_pipeline.sh          # full pipeline: nominal + every systematic variation
#   bash run_pipeline.sh nominal  # fast path: nominal only (compile+unfold+purity+draw)
#                                 # - use this while iterating on a code change, then run
#                                 # the full version before trusting/reporting results.
#
# ---------------------------------------------------------------------------------
# WHICH STAGE TO RERUN, BY WHAT YOU CHANGED
# ---------------------------------------------------------------------------------
# This is the mistake that prompted this script: puritymaker.C is a SEPARATE macro
# from unfolder.cc, easy to forget, and getting it out of sync with unfolder.cc's
# output silently corrupts specific bins rather than failing loudly (findPtBin-indexed
# lookups into a purity graph built for a different binning just read the wrong point).
#
#   Changed ana.h binning/cuts (ptBins, nUnfoldXjBins, findabcdBin thresholds, etc.),
#   unfolder.cc/.h selection or response-matrix logic, jet_object/pho_object/treeuser/
#   object .cc/.h, or the underlying ROOT trees themselves:
#     -> rerun EVERYTHING below, in order (stage 1-4). Binning/cut changes affect the
#        response matrix AND the region-A/C histograms puritymaker.C reads, so both
#        unfolder.cc and puritymaker.C are stale, and every drawing macro downstream
#        of either is stale too.
#
#   Changed only a systematic variation's OWN definition (e.g. a JES point, one entry
#   in systSources in draw_systematics.C) and nothing shared:
#     -> rerun stage 2 for just that systag, then stage 3 for just that systag if it's
#        a Data-affecting variation (jes_high/jes_low) or purity-selection-affecting
#        variation (narrowBDT/narrowISO/threejet), then stage 4.
#
#   Changed unfold_utility.h/.cc or drawer.h/.cc (shared helpers the drawing macros
#   call, not the production fill_matrix() step):
#     -> stages 1-3 are unaffected (those don't touch these files). Rerun stage 1
#        (recompile - these are part of libgammajet_unfold.so) then stage 4.
#
#   Changed only a drawing macro itself (cosmetic, a new plot, a metric definition
#   like nXjBinsForChi2's exclusion count):
#     -> rerun just that one macro. No need for any of stages 1-3.
#
#   Changed only puritymaker.C itself (not ana.h, not unfolder.cc's output):
#     -> rerun stage 3, then stage 4.
# ---------------------------------------------------------------------------------

set -e
cd "$(dirname "$0")"

MODE=${1:-full}
if [[ "$MODE" == "nominal" ]]; then
  SYSTAGS=(nominal)
else
  # Keep in sync with ana::systags (src/ana.h) - the definitive systag list - bash can't
  # read a C++ static vector<string> directly, so this is a duplicated, explicit list
  # (same convention as insitu/run_grid.sh's SYSTAGS).
  SYSTAGS=(nominal JERhigh JERlow emscale_high emscale_low jes_high jes_low threejet narrowBDT narrowISO)
fi
echo "Mode: $MODE (systags: ${SYSTAGS[*]})"

# ---- Stage 1: compile src/*.cc into libgammajet_unfold.so ----
# Needed whenever any file under src/ changed (ana.h, unfolder.cc/.h, object/pho_object/
# jet_object/treeuser/drawer/unfold_utility). Drawing macros under drawing/ are NOT
# compiled here - they're plain ROOT macros, picked up fresh (Cling JIT) on every run,
# so editing them never requires this step.
echo "=== Stage 1: compiling libgammajet_unfold.so ==="
bash make.sh

# ---- Stage 2: unfolder.cc - builds the response matrix + region A/B/C/D histograms ----
# Reads the raw ntuples; writes hists/<trigger>_<sim>_<systag>_unfolding.root (MC) and
# hists/<trigger>_<systag>_unfolding.root (Data). Full mode reads each trigger's tree
# once and fills every systag's histograms in that one pass (runall_unfold_allsys.sh /
# unfold_allsys.C) instead of once per systag - validated bin-for-bin identical to the
# old per-systag-loop output. Nominal-only mode stays on the plain single-systag
# runall_unfold.sh - it's already just one tree read, nothing to consolidate there, and
# this is the fast-iteration path so it shouldn't pay for the other 9 systags.
echo "=== Stage 2: unfolder.cc (response matrix + region A/B/C/D) ==="
if [[ "$MODE" == "nominal" ]]; then
  bash runall_unfold.sh nominal
else
  bash runall_unfold_allsys.sh
fi

# ---- Stage 3: puritymaker.C - the step that was missed ----
# Reads Data's region A/B/C/D histograms from stage 2's output and writes
# hists/purity_<systag>.root (the bootstrap "combined" graph ana::getPurity/
# getPurityError* read by array INDEX - see ana.cc). Must be rerun any time stage 2
# reran for Data with a given systag, or any time ana::ptBins/nPtBins changed (the
# whole reason this was missed this time: unfolder.cc's response matrix picked up the
# new binning automatically on rebuild, but purity_nominal.root silently kept the OLD
# binning's bootstrap points until this stage actually ran again).
echo "=== Stage 3: puritymaker.C (purity bootstrap) ==="
if [[ "$MODE" == "nominal" ]]; then
  bash run_puritymaker.sh nominal
else
  bash runall_puritymaker.sh
fi

# ---- Stage 4: drawing macros - read stages 2+3's hists/*.root, write pdfs/ + summary root files ----
# Internal order matters here:
#   - draw_prior_sensitivity.C writes hists/prior_sensitivity_nominal.root, which
#     draw_systematics.C now reads (its "priorSensitivity" source) - must run first, same
#     kind of dependency as puritymaker.C on unfolder.cc's output in stage 3. Only needs
#     the nominal chain, so it runs in both modes (unlike draw_systematics.C itself).
#   - draw_systematics.C writes hists/systematics.root, which draw_final_result.C reads -
#     it errors out cleanly (with a message) if run first.
# toy_resp_iterations.C / toy_data_iterations.C each write their own
# .toy_{resp,data}_chi2_data_<systag>.root and read the OTHER's if present, then invoke
# plot_toy_chi2_combined.C as a subprocess - run both for the full four-curve comparison
# (running only one still produces a valid partial plot, just not this script's problem
# to special-case).
echo "=== Stage 4: drawing macros ==="
cd ../drawing
root -b -l -q 'draw_purity_corrected.C("nominal")'
root -b -l -q 'draw_iteration_halfclosure.C("nominal")'
root -b -l -q 'toy_resp_iterations.C("nominal")'
root -b -l -q 'toy_data_iterations.C("nominal")'
root -b -l -q 'draw_refolding.C("nominal")'
root -b -l -q 'draw_nonclosure.C()'
root -b -l -q 'draw_prior_sensitivity.C("nominal")'
if [[ "$MODE" != "nominal" ]]; then
  root -b -l -q 'draw_systematics.C()'
  root -b -l -q 'draw_final_result.C()'
else
  echo "Skipping draw_systematics.C/draw_final_result.C in nominal-only mode (they need every systag's stage 2+3 output - rerun in full mode before trusting the final result)."
fi
cd ../macros

echo "=== Done ==="
