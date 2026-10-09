#!/bin/bash
#SBATCH --job-name=gj_draw
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
# Per radius (array index): prior sensitivity -> systematics -> final result.
source ~/.bashrc   # ROOT and GAMMAJET_UNFOLD (before set -e: /etc/bashrc may return nonzero)
set -eo pipefail
source "$GAMMAJET_UNFOLD/slurm/gammajet/lib.sh"
ps_draw "$SLURM_ARRAY_TASK_ID"
