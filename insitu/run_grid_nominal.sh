#!/bin/bash
# Runs the full prerequisite chain for a nominal-only in-situ JES check: recompile,
# then unfolder.cc for every trigger/sim run_grid.sh's five methods read (Data,
# Photon5/10/20 pythia for grid_insitu[.C/_shapechi2.C/_unfolded.C/
# _unfolded_shapechi2.C], plus Jet12_long pythia for the grid_insitu_jet12[.C/
# _shapechi2.C] method - see run_grid.sh's header comment), then puritymaker.C
# (grid_insitu.C/_shapechi2.C/_unfolded.C/_unfolded_shapechi2.C all call
# ana::getPurity/getPurityC per pT bin, which read hists/purity_nominal.root - stale
# after any Data unfolder.cc rerun with a binning/cut change, e.g. ana::ptBins - same
# dependency run_pipeline.sh's stage 3 comment documents), then run_grid.sh itself
# restricted to nominal. Equivalent in scope to run_pipeline.sh's stages 1-3 (nominal
# mode) followed by run_grid.sh --systag nominal, but skips stage 4 (drawing/) since
# the in-situ grid doesn't read its output.
#
# Use this for a quick nominal-only in-situ check after changing unfolder.cc/ana.h -
# for the full systag sweep (needed before trusting/reporting in-situ results), use
# run_pipeline.sh (full mode) followed by plain run_grid.sh instead.
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
bash run_grid.sh --systag nominal

echo "=== Done ==="
