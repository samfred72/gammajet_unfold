#!/bin/bash
# Runs the in-situ JES scan (all four methods - grid_insitu.C, grid_insitu_shapechi2.C,
# grid_insitu_unfolded.C, grid_insitu_unfolded_shapechi2.C - plus the fifth,
# Jet12_long-referenced grid_insitu_jet12.C method, plus draw_grid_chi2.C's comparison
# page) for every systag in ana::systags (src/ana.h) - nominal + the 9 named
# systematic-variation reprocessings. Every macro below now loops internally over every
# jet radius (ana::nJetR=7, R=0.2-0.8) and writes all seven into one file/one multi-page
# PDF per systag (one ana::rnames[ir] subdirectory/page per radius - see grid_insitu.C's
# header comment), so this script itself no longer loops over radius - only over systag.
#
# Usage: run_grid.sh [--systag TAG]
#   --systag TAG   Restrict the sweep to one systag (e.g. "nominal") instead of the
#                   full SYSTAGS list below - e.g. for a quick nominal-only comparison
#                   right after running unfold once, without waiting on the full
#                   10-systag sweep.
#
# SYSTAGS below must be kept in sync with ana::systags (src/ana.h) - bash can't read a
# C++ static vector<string> directly, so this is a duplicated, explicit list (same
# convention as this project's other cross-language constant duplications, e.g.
# grid_insitu_shapechi2.C's nXjBinsForChi2). Adding a systag to ana::systags does NOT
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
SYSTAGS=(nominal JERhigh JERlow emscale_high emscale_low jes_high jes_low threejet narrowBDT narrowISO)

while [[ $# -gt 0 ]]; do
  case "$1" in
    --systag)
      SYSTAGS=("$2")
      shift 2
      ;;
    *)
      echo "Usage: $0 [--systag TAG]" >&2
      exit 1
      ;;
  esac
done

for systag in "${SYSTAGS[@]}"; do
  echo "=== systag=$systag ==="
  root -b -l -q "grid_insitu.C(\"$systag\")"
  root -b -l -q "grid_insitu_shapechi2.C(\"$systag\")"
  root -b -l -q "grid_insitu_unfolded.C(\"$systag\")"
  root -b -l -q "grid_insitu_unfolded_shapechi2.C(\"$systag\")"
  # Fifth method (mean chi2) and its shape-chi2 counterpart: Data Region A (raw) vs.
  # Pythia8 Jet12_long Region A (raw) - requires
  # insitu/inputs/Jet12_long_pythia_<systag>_insitu.root, produced by running Jet12_long
  # through the normal production pipeline (runall_unfold_allsys.sh already includes
  # it - see macros/runall_unfold_allsys.sh). Jet12_long has no truth-level jet-pT cut
  # (src/treeuser.h) and is looped over every radius/systag here so draw_grid_chi2.C's
  # pads 3 and 6 can show it alongside the other four methods. (The legacy "Jet12_full"
  # cross-check that used to run as a nominal-only one-off here has been retired - see
  # src/treeuser.h/grid_insitu_jet12.C.)
  root -b -l -q "grid_insitu_jet12.C(\"$systag\")"
  root -b -l -q "grid_insitu_jet12_shapechi2.C(\"$systag\")"
  root -b -l -q "draw_grid_chi2.C(\"$systag\")"
done
