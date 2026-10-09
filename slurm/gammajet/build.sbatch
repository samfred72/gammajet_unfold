#!/bin/bash
#SBATCH --job-name=gj_build
#SBATCH --partition=acompile
#SBATCH --qos=compile
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
# Stage 1: src/make.sh. Submitted by submit.sh.
source ~/.bashrc   # ROOT and GAMMAJET_UNFOLD (before set -e: /etc/bashrc may return nonzero)
set -eo pipefail
source "$GAMMAJET_UNFOLD/slurm/gammajet/lib.sh"
ps_build
