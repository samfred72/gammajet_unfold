#!/bin/bash
# Usage: runall_unfold_allsys.sh [jet|photon|data]
# One tree read per trigger fills every systag (run_unfold_allsys.sh). herwig is nominal-only.
DODATA=0
DOPHOTON=0
DOJET=0

if [[ -z $1 ]]; then
  DOPHOTON=1
  DOJET=1
  DODATA=1
elif [[ $1 == "jet" ]]; then
  DOJET=1
elif [[ $1 == "photon" ]]; then
  DOPHOTON=1
elif [[ $1 == "data" ]]; then
  DODATA=1
fi

if [[ $DODATA == 1 ]]; then
  bash run_unfold_allsys.sh Data pythia &
fi
if [[ $DOPHOTON == 1 ]]; then
  bash run_unfold_allsys.sh Photon5 pythia &
  bash run_unfold_allsys.sh Photon10 pythia &
  bash run_unfold_allsys.sh Photon20 pythia &
  bash run_unfold.sh Photon5 herwig nominal &
  bash run_unfold.sh Photon10 herwig nominal &
  bash run_unfold.sh Photon20 herwig nominal &
fi
wait
if [[ $DOJET == 1 ]]; then
  bash run_unfold_allsys.sh Jet12_long pythia &
fi
wait
echo "All Done!"
