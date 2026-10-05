#!/bin/bash
# Usage: runall_unfold.sh [jet|photon|data] [systag]
# JER*/emscale_*/EMR* are MC-only (skip data); jes_high/low are Data-only (skip photon).
DODATA=0
DOPHOTON=0
DOJET=0
SYSTAG=${1:-nominal}

if [[ -z $2 ]]; then
  DOPHOTON=1
  DOJET=1
  DODATA=1
elif [[ $2 == "jet" ]]; then
  DOJET=1
elif [[ $2 == "photon" ]]; then
  DOPHOTON=1
elif [[ $2 == "data" ]]; then
  DODATA=1
fi

if [[ $DODATA == 1 ]]; then
  bash run_unfold.sh Data pythia $SYSTAG &
fi
if [[ $DOPHOTON == 1 ]]; then
  bash run_unfold.sh Photon5 pythia $SYSTAG &
  bash run_unfold.sh Photon10 pythia $SYSTAG &
  bash run_unfold.sh Photon20 pythia $SYSTAG &
  # herwig is only the generator systematic (compared at nominal), so run it at nominal only.
  if [[ $SYSTAG == "nominal" ]]; then
    bash run_unfold.sh Photon5 herwig nominal &
    bash run_unfold.sh Photon10 herwig nominal &
    bash run_unfold.sh Photon20 herwig nominal &
  fi
fi
wait
wait
echo "All Done!"
