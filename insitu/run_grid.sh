#!/bin/bash
# Runs the in-situ JES scan (all four methods - grid_insitu.C, grid_insitu.C (shape method),
# grid_insitu_unfolded.C, grid_insitu_unfolded.C (shape method) - plus the fifth,
# Jet12_long-referenced grid_insitu_jet12.C method, plus draw_grid_chi2.C's comparison
# page) for every systag in ana::systags (src/ana.h) - nominal + the 14 named
# systematic-variation reprocessings. Every macro below now loops internally over every
# jet radius (ana::nJetR=7, R=0.2-0.8) and writes all seven into one file/one multi-page
# PDF per systag (one ana::rnames[ir] subdirectory/page per radius - see grid_insitu.C's
# header comment), so this script itself no longer loops over radius - only over systag.
#
# Usage: run_grid.sh [--systag TAG] [--mode gammajet|combined]
#   --systag TAG   Restrict the sweep to one systag (e.g. "nominal") instead of the
#                   full SYSTAGS list below - e.g. for a quick nominal-only comparison
#                   right after running unfold once, without waiting on the full
#                   15-systag sweep.
#   --mode MODE    gammajet (default): the constant-JES gamma+jet scan, followed by
#                   draw_jes_summary.C (which rewrites src/ana.h) and draw_jes_variations.C.
#                   combined: grid_insitu.C's gamma+jet + multijet linear-JES fit
#                   (output/pdfs grid_insitu_combined_<systag>.*). Cross-check only - it
#                   never runs draw_jes_summary.C, so ana.h is not touched. Needs
#                   ../multijet/multijet_analysis_pythia.root from ../multijet/analysis
#                   (see README.md). Run the gammajet mode first if you want its constant
#                   scale overlaid on the combined plots.
#
# SYSTAGS below must be kept in sync with ana::systags (src/ana.h) - bash can't read a
# C++ static vector<string> directly, so this is a duplicated, explicit list (same
# convention as this project's other cross-language constant duplications, e.g.
# grid_insitu.C (shape method)'s nXjBinsForChi2). Adding a systag to ana::systags does NOT
# automatically add it here - it does automatically make unfolder.cc/unfold_allsys.C
# fill that systag's insitu_tree (every radius) once you re-run the production
# pipeline, and automatically makes drawing/draw_systematics.C include it in the final
# systematic combination - only this shell-side loop needs a manual one-line update.
#
# Prerequisite: unfold_allsys.C must have been re-run for Data and every Photon5/10/20
# pythia sample (macros/run_pipeline.sh or runall_unfold_allsys.sh) AFTER the
# per-(systag,radius) insitu_tree schema change in src/unfolder.h/.cc - the old
# <trigger>_<systag>_insitu.root files (no radius suffix) predate that change and won't
# be read by anything below.
SYSTAGS=(nominal JERhigh JERlow emscale_high emscale_low EMRhigh EMRlow jes_high jes_low threejet narrowBDT narrowISO narrowBDTbkg narrowISObkg wideISObkg)

MODE=gammajet
while [[ $# -gt 0 ]]; do
  case "$1" in
    --systag)
      SYSTAGS=("$2")
      shift 2
      ;;
    --mode)
      MODE="$2"
      shift 2
      ;;
    *)
      echo "Usage: $0 [--systag TAG] [--mode gammajet|combined]" >&2
      exit 1
      ;;
  esac
done
if [[ "$MODE" != gammajet && "$MODE" != combined ]]; then
  echo "Unknown --mode '$MODE' (gammajet or combined)" >&2
  exit 1
fi

if [[ "$MODE" == combined ]]; then
  for systag in "${SYSTAGS[@]}"; do
    echo "=== systag=$systag (gammajet + multijet) ==="
    root -b -l -q "grid_insitu.C(\"$systag\",\"combined\")"
  done
  exit 0
fi

for systag in "${SYSTAGS[@]}"; do
  echo "=== systag=$systag ==="
  root -b -l -q "grid_insitu.C(\"$systag\")"
  #root -b -l -q "grid_insitu.C(\"$systag\",\"gammajet\",-1,\"shape\")"
  #root -b -l -q "grid_insitu_unfolded.C(\"$systag\")"
  #root -b -l -q "grid_insitu_unfolded.C(\"$systag\",insitu_utility::scanN,\"shape\")"
  # Fifth method (mean chi2) and its shape-chi2 counterpart: Data Region A (raw) vs.
  # Pythia8 Jet12_long Region A (raw) - requires
  # insitu/inputs/Jet12_long_pythia_<systag>_insitu.root, produced by running Jet12_long
  # through the normal production pipeline (runall_unfold_allsys.sh already includes
  # it - see macros/runall_unfold_allsys.sh). Jet12_long has no truth-level jet-pT cut
  # (src/treeuser.h) and is looped over every radius/systag here so draw_grid_chi2.C's
  # pads 3 and 6 can show it alongside the other four methods. (The legacy "Jet12_full"
  # cross-check that used to run as a nominal-only one-off here has been retired - see
  # src/treeuser.h/grid_insitu_jet12.C.)
  #root -b -l -q "grid_insitu_jet12.C(\"$systag\")"
  #root -b -l -q "grid_insitu_jet12.C(\"$systag\",\"shape\")"
  #root -b -l -q "draw_grid_chi2.C(\"$systag\")"
done

# draw_jes_summary.C collects grid_insitu.C's purity-corrected p_a for every systag above
# and, when this was the full sweep (all of ana::systags, not a --systag-restricted
# subset), rewrites src/ana.h's jesNominal, jesStatErrLow/High (the nominal fit's
# statistical error, used by jes_high/jes_low) and the per-systag jesBySystag table that
# unfolder.cc uses to correct Data in each variation with that variation's own p_a (PPG18
# review issue 5). Its guard skips the ana.h update (with a warning) if any systag/radius
# combination above didn't produce a result. It prints "Updated .../ana.h" when it did
# rewrite the file - macros/run_full_pipeline.sh checks for that line. Rebuild
# (src/make.sh) afterwards for the new constants to take effect.
root -b -l -q "draw_jes_summary.C"
# p_a under every variation, per radius (reads the same outputs; no scan).
root -b -l -q "draw_jes_variations.C"
