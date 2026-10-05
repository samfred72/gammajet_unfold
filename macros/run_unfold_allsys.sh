#!/bin/bash
# Usage: run_unfold_allsys.sh <trigger> [sim]
# One tree read fills every systag (unfold_allsys.C).
TRIGGER=$1
SIM=${2:-pythia}
root -b -l -q "unfold_allsys.C(\"$TRIGGER\",\"$SIM\")"
