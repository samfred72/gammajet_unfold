#!/bin/bash
# Usage: run_unfold_allsys.sh <trigger> [sim]
# One-pass, all-systematics unfold (see unfold_allsys.C) - reads the raw tree once and
# fills every systag's histograms in that one pass, instead of run_unfold.sh's one
# process per systag (one tree read per systag). Systags processed are unfold_allsys.C's
# default list (nominal + the 9 named variations) - see its header comment.
TRIGGER=$1
SIM=${2:-pythia}
root -b -l -q "unfold_allsys.C(\"$TRIGGER\",\"$SIM\")"
