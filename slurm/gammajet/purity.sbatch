#!/bin/bash
#SBATCH --job-name=gj_purity
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
# puritymaker.C for one systag: array index into GJ_SYSTAGS.
source ~/.bashrc   # ROOT and GAMMAJET_UNFOLD (before set -e: /etc/bashrc may return nonzero)
set -eo pipefail
source "$GAMMAJET_UNFOLD/slurm/gammajet/lib.sh"
systags=($GJ_SYSTAGS)
ps_purity "${systags[$SLURM_ARRAY_TASK_ID]}"
