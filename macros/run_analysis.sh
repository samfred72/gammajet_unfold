#!/bin/bash
# Runs the full analysis: on Slurm (slurm/gammajet/submit.sh) if it responds, otherwise locally in the
# background (macros/run_full_pipeline.sh, slower: one machine). Both run the same stages
# (macros/pipeline_steps.sh).
#
# Usage: bash macros/run_analysis.sh [--slurm | --local] [--jobs N (local, 8)] [--from STAGE] [--logdir DIR] [--no-draw]
#   default: Slurm if `sinfo` and `squeue` answer within 30 s, else local.
#   --from/--logdir resume a chain (see run_full_pipeline.sh), e.g. a Slurm chain that stopped:
#   bash macros/run_analysis.sh --local --from <failed stage> --logdir logs/slurm/gammajet_<timestamp>
set -eo pipefail
MACRODIR="$(cd "$(dirname "$0")" && pwd)"
export GAMMAJET_UNFOLD="${GAMMAJET_UNFOLD:-$(cd "$MACRODIR/.." && pwd)}"

WHERE=auto
PASS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --slurm)   WHERE=slurm; shift ;;
    --local)   WHERE=local; shift ;;
    --jobs)    PASS+=(--jobs "$2"); shift 2 ;;
    --no-draw) PASS+=(--no-draw); shift ;;
    --from)    PASS+=(--from "$2"); shift 2 ;;
    --logdir)  PASS+=(--logdir "$2"); shift 2 ;;
    *) echo "Usage: $0 [--slurm | --local] [--jobs N] [--from STAGE] [--logdir DIR] [--no-draw]" >&2; exit 1 ;;
  esac
done

slurm_up() {
  command -v sbatch > /dev/null && timeout 30 sinfo -h > /dev/null 2>&1 && timeout 30 squeue -h -u "$USER" > /dev/null 2>&1
}
if [[ $WHERE == auto ]]; then
  if slurm_up; then WHERE=slurm; else WHERE=local; echo "Slurm is not responding: running locally."; fi
fi

if [[ $WHERE == slurm ]]; then
  # --jobs is local only.
  slurm_args=()
  for (( i = 0; i < ${#PASS[@]}; i++ )); do
    case ${PASS[$i]} in
      --jobs) i=$((i+1)) ;;
      --from|--logdir) slurm_args+=("${PASS[$i]}" "${PASS[$((i+1))]}"); i=$((i+1)) ;;
      *) slurm_args+=("${PASS[$i]}") ;;
    esac
  done
  bash "$GAMMAJET_UNFOLD/slurm/gammajet/submit.sh" "${slurm_args[@]}"
else
  mkdir -p "$GAMMAJET_UNFOLD/logs"
  out="$GAMMAJET_UNFOLD/logs/full_pipeline_local_$(date +%Y%m%d_%H%M%S).out"
  nohup bash "$MACRODIR/run_full_pipeline.sh" "${PASS[@]}" > "$out" 2>&1 &
  echo "Local pipeline started (pid $!). Progress: $out"
fi
