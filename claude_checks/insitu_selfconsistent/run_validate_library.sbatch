#!/bin/bash
#SBATCH --job-name=sc_validate_lib
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/projects/sali6653/PPG18_analysis/claude_checks/insitu_selfconsistent/validate_library_%j.out
source ~/.bashrc
cd /projects/sali6653/PPG18_analysis/claude_checks/insitu_selfconsistent && root -b -l -q validate_library.C 2>&1 | grep -v "^Info\|Jet[0-9]*_pythia"
