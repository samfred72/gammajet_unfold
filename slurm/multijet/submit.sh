#!/bin/bash
# Submits the multijet balance analysis to Slurm as a dependency chain:
#   build -> analysis (array over sims) -> makeratio -> analysis -> ... -> draw (one job)
# Each analysis pass after the first reweights the MC with the fits makeratio.C made from the
# previous pass (see multijet/README.md).
#
# Usage: bash slurm/multijet/submit.sh [--passes N (2)] [--sims "pythia herwig"] [--no-draw]
#                                      [--with-jet5] [--with-jet8] [--truth-smear] [--no-smear] [--tight]
#   Jet5 and Jet8 are left out by default (Jet12 then starts at truth pT 0); --with-jet5/--with-jet8 use them.
#   --truth-smear, --no-smear  analysis.cc's MC smearing options
#   --tight    test cut (analysis --tight): one pass with the current aux/ fits, no makeratio,
#              output multijet_analysis_<sim>_tight.root and pdfs/*_tight.pdf
# The draw job writes pdfs/hxj_means_<sim>.pdf, hxj_dists_<sim>.pdf and xj_samples_<sim>.pdf, each with
# every radius and JER variation.
# Logs and each pass's output go to logs/slurm/multijet_<timestamp>/.
set -eo pipefail
: "${GAMMAJET_UNFOLD:?set GAMMAJET_UNFOLD to the gammajet_unfold checkout (see README.md)}"

HERE="$GAMMAJET_UNFOLD/slurm/multijet"
PASSES=2
SIMS="pythia herwig"
DRAW=1
export MJ_ANAOPTS=""
export MJ_TAG=""
JET5=0
JET8=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --passes)  PASSES="$2"; shift 2 ;;
    --sims)    SIMS="$2"; shift 2 ;;
    --no-draw) DRAW=0; shift ;;
    --with-jet5)   JET5=1; shift ;;
    --with-jet8)   JET8=1; shift ;;
    --truth-smear) MJ_ANAOPTS+=" --truth-smear"; shift ;;
    --no-smear)    MJ_ANAOPTS+=" --no-smear"; shift ;;
    --tight)       MJ_ANAOPTS+=" --tight"; MJ_TAG="_tight"; shift ;;
    *) echo "Usage: $0 [--passes N] [--sims \"pythia herwig\"] [--no-draw] [--with-jet5] [--with-jet8] [--truth-smear] [--no-smear] [--tight]" >&2; exit 1 ;;
  esac
done
(( JET5 )) || MJ_ANAOPTS+=" --no-jet5"
(( JET8 )) || MJ_ANAOPTS+=" --no-jet8"

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
[[ -n "$MJ_TAG" ]] && PASSES=1  # the tight test reuses the current reweighting fits

export MJ_LOGDIR="$GAMMAJET_UNFOLD/logs/slurm/multijet_$(date +%Y%m%d_%H%M%S)"
export MJ_SIMS="$SIMS"
mkdir -p "$MJ_LOGDIR"
echo "Sims: $MJ_SIMS   passes: $PASSES   draw: $DRAW   analysis options:${MJ_ANAOPTS:- none}"
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
  jid=$(sub --dependency=afterok:$jid -o "$MJ_LOGDIR/draw.log" "$HERE/draw.sbatch")
  echo "draw          $jid"
fi
echo "Check with: squeue -u $USER"
