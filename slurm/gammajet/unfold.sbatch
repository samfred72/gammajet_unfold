#!/bin/bash
#SBATCH --job-name=gj_unfold
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=02:00:00
# Unfold one task: array index into GJ_UNFOLD_TASKS (trigger:sim:mode).
source ~/.bashrc   # ROOT and GAMMAJET_UNFOLD (before set -e: /etc/bashrc may return nonzero)
set -eo pipefail
source "$GAMMAJET_UNFOLD/slurm/gammajet/lib.sh"
tasks=($GJ_UNFOLD_TASKS)
ps_unfold "${tasks[$SLURM_ARRAY_TASK_ID]}"
