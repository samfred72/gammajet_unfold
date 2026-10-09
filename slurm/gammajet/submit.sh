#!/bin/bash
# Submits the gamma+jet pipeline to Slurm as one dependency chain, the stages of
# macros/pipeline_steps.sh with each per-sample, per-systag and per-radius step as an array task.
# macros/run_analysis.sh calls this when Slurm is up, else runs macros/run_full_pipeline.sh locally.
#
# Usage: bash slurm/gammajet/submit.sh [--from STAGE] [--logdir DIR] [--no-draw]
#   --from STAGE  resume at STAGE (pipeline_steps.sh PS_STAGES), reusing the outputs of earlier stages
#   --logdir DIR  log into an earlier run's directory (keeps its starting JES table for jes_update)
# Logs: logs/slurm/gammajet_<timestamp>/ (jes_update.log reports the in-situ result and its change
# from the starting table).
set -eo pipefail
: "${GAMMAJET_UNFOLD:?set GAMMAJET_UNFOLD to the gammajet_unfold checkout (see README.md)}"
source "$GAMMAJET_UNFOLD/slurm/gammajet/lib.sh"

export GJ_DRAW=1
FROM=build
LOGDIR_ARG=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --from)    FROM="$2"; shift 2 ;;
    --logdir)  LOGDIR_ARG="$2"; shift 2 ;;
    --no-draw) GJ_DRAW=0; shift ;;
    *) echo "Usage: $0 [--from STAGE] [--logdir DIR] [--no-draw]" >&2; exit 1 ;;
  esac
done
FROM_I=$(ps_stage_index "$FROM")
if [[ -n $LOGDIR_ARG ]]; then export GJ_LOGDIR="$(cd "$LOGDIR_ARG" && pwd)"; ps_init_logdir "$GJ_LOGDIR" 0
else export GJ_LOGDIR="$GAMMAJET_UNFOLD/logs/slurm/gammajet_$(date +%Y%m%d_%H%M%S)"; ps_init_logdir "$GJ_LOGDIR" 1; fi
echo "Systags: $(ps_systags)"
echo "from $FROM, draw $GJ_DRAW. Logs: $GJ_LOGDIR"
want() { (( $(ps_stage_index "$1") >= FROM_I )); }

jid=""
if want build;        then jid=$(gj_sub -o "$GJ_LOGDIR/build.log" "$GJ_SLURM/build.sbatch"); echo "build               $jid"; fi
if want unfold_all;   then jid=$(gj_submit_unfold all "$jid" "$PS_UNFOLD_ALL"); fi
if want purity;       then jid=$(gj_submit_systag_array purity purity.sbatch "$jid"); fi
if want reweight;     then jid=$(gj_sub $(gj_dep "$jid") -o "$GJ_LOGDIR/reweight.log" "$GJ_SLURM/reweight.sbatch"); echo "reweight            $jid"; fi
if want unfold_mc;    then jid=$(gj_submit_unfold mc "$jid" "$PS_UNFOLD_MC"); fi
if want insitu;       then jid=$(gj_submit_systag_array insitu insitu.sbatch "$jid"); fi
if want jes_update;   then jid=$(gj_sub $(gj_dep "$jid") -o "$GJ_LOGDIR/jes_update.log" "$GJ_SLURM/jes_update.sbatch"); echo "jes_update          $jid"; fi
if want unfold_data;  then jid=$(gj_submit_unfold data "$jid" "$PS_UNFOLD_DATA"); fi
if want purity_final; then jid=$(gj_submit_systag_array purity_final purity.sbatch "$jid"); fi
gj_submit_final "$jid" > /dev/null
echo "Check with: squeue -u $USER"
