#!/bin/bash
# run_pipeline.sh: the gammajet_unfold analysis chain.
#
# Usage:
#   bash run_pipeline.sh          # nominal + every systematic variation
#   bash run_pipeline.sh nominal  # nominal only (fast iteration)
#   bash run_pipeline.sh draw     # stage 4 only, on existing stage 2+3 output
#
# What to rerun after a change:
#   ana.h binning/cuts, unfolder, object/treeuser code, or the trees -> stages 1-4.
#   One systag's own definition -> stage 2 for that systag, stage 3 if it affects Data or the
#     purity selection (jes_*, narrow*/wide*, threejet), then stage 4.
#   unfold_utility or drawer -> stage 1, then stage 4.
#   puritymaker.C only -> stage 3, then stage 4.
#   A drawing macro only -> just that macro.
# A purity file out of sync with unfolder's binning corrupts bins silently.

set -e
cd "$(dirname "$0")"

MODE=${1:-full}
if [[ "$MODE" != "full" && "$MODE" != "nominal" && "$MODE" != "draw" ]]; then
  echo "Usage: $0 [full|nominal|draw]" >&2
  exit 1
fi
if [[ "$MODE" == "nominal" ]]; then
  SYSTAGS=(nominal)
else
  # Keep in sync with ana::systags (src/ana.h).
  SYSTAGS=(nominal JERhigh JERlow emscale_high emscale_low EMRhigh EMRlow jes_high jes_low threejet narrowBDT narrowISO narrowBDTbkg narrowISObkg wideISObkg timingwide)
fi
echo "Mode: $MODE (systags: ${SYSTAGS[*]})"

# ---- Stage 1: compile src/*.cc into libgammajet_unfold.so ----
if [[ "$MODE" != "draw" ]]; then
echo "=== Stage 1: compiling libgammajet_unfold.so ==="
bash make.sh

# ---- Stage 2: unfolder (response matrix + region A/B/C/D histograms) ----
# Full mode fills every systag in one read of each tree (unfold_allsys.C).
echo "=== Stage 2: unfolder.cc (response matrix + region A/B/C/D) ==="
if [[ "$MODE" == "nominal" ]]; then
  bash runall_unfold.sh nominal
else
  bash runall_unfold_allsys.sh
fi

# ---- Stage 3: puritymaker.C -> hists/purity_<systag>.root ----
echo "=== Stage 3: puritymaker.C (purity bootstrap) ==="
if [[ "$MODE" == "nominal" ]]; then
  bash run_puritymaker.sh nominal
else
  bash runall_puritymaker.sh
fi

# ---- Stage 4: drawing ----
# Order: draw_prior_sensitivity -> draw_systematics -> draw_final_result (each reads the
# previous one's output); only these three run per radius. The other macros are R=0.4
# nominal diagnostics. The two toy_iterations runs read each other's output for the
# combined plot.
fi # end of stages 1-3 (skipped in draw mode)

echo "=== Stage 4: drawing macros ==="
cd ../drawing
root -b -l -q 'draw_purity_corrected.C("nominal")'
root -b -l -q 'draw_iteration_halfclosure.C("nominal")'
root -b -l -q 'toy_iterations.C("nominal","resp")'
root -b -l -q 'toy_iterations.C("nominal","data")'
root -b -l -q 'draw_refolding.C("nominal")'
root -b -l -q 'draw_nonclosure.C()'
root -b -l -q 'draw_prior_sensitivity.C("nominal")'
if [[ "$MODE" != "nominal" ]]; then
  # Keep in sync with ana::nJetR/ana::rnames (src/ana.h).
  RNAMES=(R02 R03 R04 R05 R06 R07 R08)
  for ir in 0 1 2 3 4 5 6; do
    echo "--- radius index $ir (${RNAMES[$ir]}) ---"
    if [[ $ir -ne 2 ]]; then
      root -b -l -q "draw_prior_sensitivity.C(\"nominal\", $ir)"
    fi
    root -b -l -q "draw_systematics.C($ir)"
    root -b -l -q "draw_final_result.C($ir)"
  done
else
  echo "Skipping draw_systematics.C/draw_final_result.C in nominal-only mode (they need every systag's stage 2+3 output - rerun in full mode before trusting the final result)."
fi
cd ../macros

echo "=== Done ==="
