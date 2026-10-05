#!/bin/bash
#SBATCH --job-name=mj_draw
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
# draw_xj for one radius (array index 0-6 = R0.2-R0.8), every JER variation and sim in $MJ_SIMS.
# Submitted by submit.sh.
source ~/.bashrc   # ROOT and GAMMAJET_UNFOLD (before set -e: /etc/bashrc may return nonzero)
set -eo pipefail
cd "$GAMMAJET_UNFOLD/multijet"
RADIUS=$(( SLURM_ARRAY_TASK_ID + 2 ))
for sim in $MJ_SIMS; do
  for sys in RECO HIGH LOW; do
    root -b -l -q "draw_xj.C($RADIUS,\"$sys\",\"$sim\")"
  done
done
