#!/bin/bash
# Runs the in-situ JES closure test (insitu_closure/) end to end.
#
# Stage 1 (make_closure_trees.C) splits ONE MC sample's events into a "data" half (jet
# pt scaled by a single injected JES factor) and a "sim" half (nominal reco, the fixed
# reference). It's run here once per Photon5/Photon10/Photon20 pythia sample, as three
# parallel background ROOT processes (same &/wait pattern as macros/runall_unfold.sh) -
# looping the ~9.8M-entry Photon20 and ~8.7M-entry Photon10 raw ntuples sequentially in
# one process would otherwise dominate the wall time, and the three samples are
# otherwise fully independent (each only touches its own output files). The injected
# scale is fixed ONCE, here, and passed identically to all three - each process getting
# its own random injected value would break the whole premise of combining them.
#
# Stage 2 (fit_closure.C) reads all three samples' insitu trees back, stitches them by
# cross-section weight (same weights make_closure_trees.C's photon_scale map uses), and
# scans for the best-fit jet-energy-scale factor at R=0.4 (this project's nominal jet
# radius - see fit_closure.C's header comment) - reported against the known injected
# value in insitu_closure/pdfs/fit_closure.pdf and insitu_closure/output/fit_closure.root.
#
# Usage: run_closure.sh [scale]
#   scale: the exact injected "data"-half JES-scale factor to test, e.g. 0.95. If
#   omitted, one is drawn uniformly from [0.9,1.1) instead. Either way, fit_closure.C's
#   own pa scan window is a fixed [0.80,1.20] (see that file) - a scale outside that
#   range can't be recovered, so this script warns (but does not refuse) if you pass one.
set -e

DEFAULT_LOW=0.9
DEFAULT_HIGH=1.1
SCANLOW=0.80
SCANHIGH=1.20

if [[ -n "$1" ]]; then
  SCALE=$1
  echo "=== Injected closure JES scale: $SCALE (user-specified) ==="
  if (( $(awk -v s="$SCALE" -v lo="$SCANLOW" -v hi="$SCANHIGH" 'BEGIN{print (s<lo || s>hi)}') )); then
    echo "WARNING: $SCALE is outside fit_closure.C's [$SCANLOW,$SCANHIGH] pa scan window - it will not be recoverable. Widen paLow/paHigh in fit_closure.C first." >&2
  fi
else
  # awk instead of a throwaway ROOT process, since this is just a single uniform-random draw.
  SCALE=$(awk -v lo="$DEFAULT_LOW" -v hi="$DEFAULT_HIGH" -v seed="$RANDOM$RANDOM$$" \
    'BEGIN{srand(seed); printf "%.6f", lo+rand()*(hi-lo)}')
  echo "=== Injected closure JES scale: $SCALE (drawn from [$DEFAULT_LOW,$DEFAULT_HIGH)) ==="
fi

# Distinct fixed per-sample seed offsets (from a random base) for the per-event
# response-half coin flip - keeps the three parallel processes' random streams
# independent of each other without relying on TRandom3's seed=0 auto-seed to
# desynchronize processes launched in the same instant.
BASESEED=$RANDOM
echo "=== Stage 1: building closure trees for Photon5/10/20 pythia (parallel) ==="
root -b -l -q "make_closure_trees.C(\"Photon5\",\"pythia\",$SCALE,0.9,1.1,$((BASESEED+1)))"  > output/make_closure_trees_Photon5.log  2>&1 &
root -b -l -q "make_closure_trees.C(\"Photon10\",\"pythia\",$SCALE,0.9,1.1,$((BASESEED+2)))" > output/make_closure_trees_Photon10.log 2>&1 &
root -b -l -q "make_closure_trees.C(\"Photon20\",\"pythia\",$SCALE,0.9,1.1,$((BASESEED+3)))" > output/make_closure_trees_Photon20.log 2>&1 &
wait
tail -n2 output/make_closure_trees_Photon5.log output/make_closure_trees_Photon10.log output/make_closure_trees_Photon20.log

echo "=== Stage 2: fitting for the best-fit JES scale per jet radius ==="
root -b -l -q "fit_closure.C()"

echo "=== Done - see pdfs/fit_closure.pdf and output/fit_closure.root ==="
