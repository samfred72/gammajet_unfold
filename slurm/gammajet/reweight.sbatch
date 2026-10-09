#!/bin/bash
#SBATCH --job-name=gj_reweight
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
# v_z / cluster-pT reweight from the stage-3 purities.
source ~/.bashrc   # ROOT and GAMMAJET_UNFOLD (before set -e: /etc/bashrc may return nonzero)
set -eo pipefail
source "$GAMMAJET_UNFOLD/slurm/gammajet/lib.sh"
ps_reweight
