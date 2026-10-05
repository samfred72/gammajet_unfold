#!/bin/bash
#SBATCH --job-name=mj_build
#SBATCH --partition=acompile
#SBATCH --qos=compile
#SBATCH --ntasks=1
#SBATCH --time=00:15:00
# Builds multijet/analysis (submitted by submit.sh).
source ~/.bashrc   # ROOT and GAMMAJET_UNFOLD (before set -e: /etc/bashrc may return nonzero)
set -eo pipefail
bash "$GAMMAJET_UNFOLD/multijet/make.sh"
ls -l "$GAMMAJET_UNFOLD/multijet/analysis"
