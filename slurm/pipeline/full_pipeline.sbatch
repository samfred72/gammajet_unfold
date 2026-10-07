#!/bin/bash
#SBATCH --job-name=gj_full_pipeline
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
# macros/run_full_pipeline.sh on one node: unfold_allsys and puritymaker run their triggers and
# systags in parallel (up to 15 processes). Usage: sbatch slurm/pipeline/full_pipeline.sbatch [args]
source ~/.bashrc   # ROOT and GAMMAJET_UNFOLD (before set -e: /etc/bashrc may return nonzero)
set -eo pipefail
bash "$GAMMAJET_UNFOLD/macros/run_full_pipeline.sh" "$@"
