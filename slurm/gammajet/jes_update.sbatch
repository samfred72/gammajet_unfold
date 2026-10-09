#!/bin/bash
#SBATCH --job-name=gj_jes_update
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
# Combine the scans into src/ana.h, report the change from GJ_LOGDIR/jes_start.txt, recompile.
source ~/.bashrc   # ROOT and GAMMAJET_UNFOLD (before set -e: /etc/bashrc may return nonzero)
set -eo pipefail
source "$GAMMAJET_UNFOLD/slurm/gammajet/lib.sh"
ps_jes_update "$GJ_LOGDIR"
