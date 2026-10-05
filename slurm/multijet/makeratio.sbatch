#!/bin/bash
#SBATCH --job-name=mj_makeratio
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
# Refits the MC reweighting (aux/ratio<R>_<sim>.root) from the last pass's output.
# Submitted by submit.sh.
source ~/.bashrc   # ROOT and GAMMAJET_UNFOLD (before set -e: /etc/bashrc may return nonzero)
set -eo pipefail
cd "$GAMMAJET_UNFOLD/multijet"
root -b -l -q makeratio.C
mkdir -p "$MJ_LOGDIR/aux_after_${SLURM_JOB_ID}"
cp aux/ratio*.root "$MJ_LOGDIR/aux_after_${SLURM_JOB_ID}/"
