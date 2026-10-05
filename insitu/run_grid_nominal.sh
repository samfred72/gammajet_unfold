#!/bin/bash
# Quick nominal-only in-situ check: recompile, unfold Data, Photon5/10/20 and Jet12_long, rerun
# puritymaker (the scans read the purities), then run_grid.sh --systag nominal. For results,
# run the full pipeline and plain run_grid.sh instead.
#
# Usage: run_grid_nominal.sh [--mode gammajet|combined]  (passed to run_grid.sh)
set -e
cd "$(dirname "$0")"

echo "=== Stage 1: compiling libgammajet_unfold.so ==="
(cd ../macros && bash make.sh)

echo "=== Stage 2: unfolder.cc (nominal) for Data + Photon5/10/20 pythia + Jet12_long pythia ==="
(
  cd ../macros
  bash run_unfold.sh Data pythia nominal &
  bash run_unfold.sh Photon5 pythia nominal &
  bash run_unfold.sh Photon10 pythia nominal &
  bash run_unfold.sh Photon20 pythia nominal &
  bash run_unfold.sh Jet12_long pythia nominal &
  wait
)

echo "=== Stage 3: puritymaker.C (nominal) ==="
(cd ../macros && bash run_puritymaker.sh nominal)

echo "=== Stage 4: in-situ grid (nominal) ==="
bash run_grid.sh --systag nominal "$@"

echo "=== Done ==="
