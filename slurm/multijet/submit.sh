#!/bin/bash
# Submits the multijet balance analysis to Slurm as a dependency chain:
#   build -> analysis (array over sims) -> makeratio -> analysis -> ... -> draw_xj (array over radii)
# Each analysis pass after the first reweights the MC with the fits makeratio.C made from the
# previous pass (see multijet/README.md).
#
# Usage: bash slurm/multijet/submit.sh [--passes N (2)] [--sims "pythia herwig"] [--no-draw]
# Logs and each pass's output go to logs/slurm/multijet_<timestamp>/.
set -eo pipefail
: "${GAMMAJET_UNFOLD:?set GAMMAJET_UNFOLD to the gammajet_unfold checkout (see README.md)}"

HERE="$GAMMAJET_UNFOLD/slurm/multijet"
PASSES=2
SIMS="pythia herwig"
DRAW=1
while [[ $# -gt 0 ]]; do
  case "$1" in
    --passes)  PASSES="$2"; shift 2 ;;
    --sims)    SIMS="$2"; shift 2 ;;
    --no-draw) DRAW=0; shift ;;
    *) echo "Usage: $0 [--passes N] [--sims \"pythia herwig\"] [--no-draw]" >&2; exit 1 ;;
  esac
done

# analysis.cc's sim index: 0 = pythia, 1 = herwig
IDX=()
for sim in $SIMS; do
  case "$sim" in
    pythia) IDX+=(0) ;;
    herwig) IDX+=(1) ;;
    *) echo "Unknown sim '$sim' (pythia or herwig)" >&2; exit 1 ;;
  esac
done
ARRAY=$(IFS=,; echo "${IDX[*]}")

export MJ_LOGDIR="$GAMMAJET_UNFOLD/logs/slurm/multijet_$(date +%Y%m%d_%H%M%S)"
export MJ_SIMS="$SIMS"
mkdir -p "$MJ_LOGDIR"
echo "Sims: $MJ_SIMS   passes: $PASSES   draw: $DRAW"
echo "Logs: $MJ_LOGDIR"

sub() { sbatch --parsable --export=ALL "$@"; }

jid=$(sub -o "$MJ_LOGDIR/build.log" "$HERE/build.sbatch")
echo "build         $jid"
for (( p = 1; p <= PASSES; p++ )); do
  jid=$(sub --dependency=afterok:$jid --array="$ARRAY" -J mj_analysis_p$p \
            -o "$MJ_LOGDIR/analysis_pass${p}_%a.log" "$HERE/analysis.sbatch" "$p")
  echo "analysis p$p   $jid (array $ARRAY)"
  if (( p < PASSES )); then
    jid=$(sub --dependency=afterok:$jid -J mj_makeratio_p$p \
              -o "$MJ_LOGDIR/makeratio_pass$p.log" "$HERE/makeratio.sbatch")
    echo "makeratio p$p  $jid"
  fi
done
if (( DRAW )); then
  jid=$(sub --dependency=afterok:$jid --array=0-6 -o "$MJ_LOGDIR/draw_%a.log" "$HERE/draw.sbatch")
  echo "draw_xj       $jid (array 0-6 = R0.2-R0.8)"
fi
echo "Check with: squeue -u $USER"
