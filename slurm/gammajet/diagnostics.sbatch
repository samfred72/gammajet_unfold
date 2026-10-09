#!/bin/bash
#SBATCH --job-name=gj_diagnostics
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=03:00:00
# R=0.4 nominal unfolding diagnostics.
source ~/.bashrc   # ROOT and GAMMAJET_UNFOLD (before set -e: /etc/bashrc may return nonzero)
set -eo pipefail
source "$GAMMAJET_UNFOLD/slurm/gammajet/lib.sh"
ps_diagnostics
