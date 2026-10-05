#!/bin/bash
#SBATCH --job-name=mj_analysis
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
# One analysis pass for one sim: array index = analysis.cc's sim argument (0 pythia, 1 herwig).
# $1 = pass number. Submitted by submit.sh.
source ~/.bashrc   # ROOT and GAMMAJET_UNFOLD (before set -e: /etc/bashrc may return nonzero)
set -eo pipefail
SIM=$([[ $SLURM_ARRAY_TASK_ID == 1 ]] && echo herwig || echo pythia)
cd "$GAMMAJET_UNFOLD/multijet"
./analysis "$SLURM_ARRAY_TASK_ID" "$GAMMAJET_UNFOLD/trees"
# keep each pass's output so the reweighting passes can be compared
mkdir -p "$MJ_LOGDIR/pass$1"
cp "multijet_analysis_$SIM.root" "$MJ_LOGDIR/pass$1/"
