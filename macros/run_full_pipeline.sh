#!/bin/bash
# run_full_pipeline.sh: the whole analysis on this machine, without Slurm. The same stages as
# slurm/gammajet/submit.sh (both use macros/pipeline_steps.sh); the per-sample, per-systag and
# per-radius steps of a stage run in parallel, at most --jobs at a time. See pipeline_steps.sh for
# the stage list. macros/run_analysis.sh picks Slurm or this script.
#
# Usage: bash macros/run_full_pipeline.sh [--jobs N (8)] [--from STAGE] [--logdir DIR] [--no-draw]
#   --from STAGE  resume at STAGE (pipeline_steps.sh PS_STAGES), reusing the outputs of earlier stages
#                 (e.g. to finish locally a Slurm chain that stopped)
#   --logdir DIR  log into an earlier run's directory (keeps its starting JES table for jes_update)
# Logs go to logs/full_pipeline_<timestamp>/ (one file per task). Long-running; launch in the
# background (run_analysis.sh --local does).
set -eo pipefail
MACRODIR="$(cd "$(dirname "$0")" && pwd)"
export GAMMAJET_UNFOLD="${GAMMAJET_UNFOLD:-$(cd "$MACRODIR/.." && pwd)}"
source "$MACRODIR/pipeline_steps.sh"

JOBS=8
DRAW=1
FROM=build
LOGDIR=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --jobs)    JOBS="$2"; shift 2 ;;
    --from)    FROM="$2"; shift 2 ;;
    --logdir)  LOGDIR="$2"; shift 2 ;;
    --no-draw) DRAW=0; shift ;;
    *) echo "Usage: $0 [--jobs N] [--from STAGE] [--logdir DIR] [--no-draw]" >&2; exit 1 ;;
  esac
done
FROM_I=$(ps_stage_index "$FROM")
if [[ -n $LOGDIR ]]; then LOGDIR="$(cd "$LOGDIR" && pwd)"; ps_init_logdir "$LOGDIR" 0
else LOGDIR="$GAMMAJET_UNFOLD/logs/full_pipeline_$(date +%Y%m%d_%H%M%S)"; ps_init_logdir "$LOGDIR" 1; fi
SYSTAGS=$(ps_systags)
echo "Logs: $LOGDIR   (jobs $JOBS, from $FROM, draw $DRAW)"
echo "Systags: $SYSTAGS"
want() { (( $(ps_stage_index "$1") >= FROM_I )); }

# stage NAME FUNC ARG...: FUNC once per ARG, at most $JOBS at a time, each logged to
# $LOGDIR/NAME_<ARG>.log; stops the pipeline if any task fails.
stage() {
  local name=$1 fn=$2; shift 2
  local t0=$(date +%s) fail=0 arg
  echo "[$(date +%T)] $name ($# tasks) ..."
  for arg in "$@"; do
    while (( $(jobs -rp | wc -l) >= JOBS )); do wait -n || fail=1; done
    ( "$fn" "$arg" ) > "$LOGDIR/${name}_${arg//:/_}.log" 2>&1 &
  done
  while (( $(jobs -rp | wc -l) > 0 )); do wait -n || fail=1; done
  if (( fail )); then echo "[$(date +%T)] FAILED: $name - see $LOGDIR/${name}_*.log"; exit 1; fi
  echo "[$(date +%T)] $name done ($(( ($(date +%s) - t0) / 60 )) min)"
}

want build        && stage 01_build        ps_build      build
want unfold_all   && stage 02_unfold_all   ps_unfold     $PS_UNFOLD_ALL
want purity       && stage 03_purity       ps_purity     $SYSTAGS
want reweight     && stage 04_reweight     ps_reweight   reweight
want unfold_mc    && stage 05_unfold_mc    ps_unfold     $PS_UNFOLD_MC
want insitu       && stage 06_insitu       ps_insitu     $SYSTAGS
jes_update_here() { ps_jes_update "$LOGDIR"; }
if want jes_update; then
  stage 07_jes_update jes_update_here update
  grep -h "max |delta p_a|" "$LOGDIR"/07_jes_update_*.log | sed 's/^/    /'
fi
want unfold_data  && stage 08_unfold_data  ps_unfold     $PS_UNFOLD_DATA
want purity_final && stage 09_purity       ps_purity     $SYSTAGS
if [[ $DRAW -eq 1 ]]; then
  ps_final_task() { case $1 in insitu_plots) ps_insitu_plots ;; diagnostics) ps_diagnostics ;; *) ps_draw "$1" ;; esac; }
  stage 10_final           ps_final_task insitu_plots diagnostics $PS_RADII
fi
echo "[$(date +%T)] Full pipeline done. Logs: $LOGDIR"
