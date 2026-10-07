#!/bin/bash
#SBATCH --job-name=mj_draw
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
# For every sim in $MJ_SIMS: the combined x_j PDFs (every radius and JER variation in one file each).
# Submitted by submit.sh.
source ~/.bashrc   # ROOT and GAMMAJET_UNFOLD (before set -e: /etc/bashrc may return nonzero)
set -eo pipefail
cd "$GAMMAJET_UNFOLD/multijet"
LEAD=$([[ "$MJ_TAG" == "_tight" ]] && echo 25 || echo 20)
for sim in $MJ_SIMS; do
  root -b -l -q -e ".L draw_xj.C" -e "draw_xj_all(\"$sim\",\"$MJ_TAG\")"
  root -b -l -q "draw_xj_samples.C(\"multijet_analysis_$sim$MJ_TAG.root\",\"pdfs/xj_samples_$sim$MJ_TAG.pdf\",$LEAD,14)"
done
