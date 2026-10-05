#!/bin/bash
# In-situ JES scans for every systag (all radii inside each macro).
#
# Usage: run_grid.sh [--systag TAG] [--mode gammajet|combined]
#   --systag TAG  one systag only (quick check; ana.h is not updated).
#   --mode MODE   gammajet (default): constant-JES scan, then draw_jes_summary.C (rewrites
#                 src/ana.h) and draw_jes_variations.C.
#                 combined: gamma+jet + multijet linear fit (cross-check; ana.h untouched).
#                 Needs ../multijet/multijet_analysis_pythia.root.
#
# Keep SYSTAGS in sync with ana::systags (src/ana.h).
# Prerequisite: unfold_allsys.C run for Data and Photon5/10/20.
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
  # Other methods (optional; Jet12 needs the Jet12_long insitu input):
  # root -b -l -q "grid_insitu.C(\"$systag\",\"gammajet\",-1,\"shape\")"
  # root -b -l -q "grid_insitu_unfolded.C(\"$systag\")"
  # root -b -l -q "grid_insitu_unfolded.C(\"$systag\",insitu_utility::scanN,\"shape\")"
  # root -b -l -q "grid_insitu_jet12.C(\"$systag\")"
  # root -b -l -q "grid_insitu_jet12.C(\"$systag\",\"shape\")"
  # root -b -l -q "draw_grid_chi2.C(\"$systag\")"
done

# draw_jes_summary.C rewrites src/ana.h's JES constants only after a full sweep (prints
# "Updated .../ana.h", which run_full_pipeline.sh checks). Rebuild afterwards.
root -b -l -q "draw_jes_summary.C"
# p_a under every variation, per radius (no scan).
root -b -l -q "draw_jes_variations.C"
