#/bin/bash
# Usage: runall_unfold.sh [jet|photon|data] [systag]
# systag: nominal (default), JERhigh, JERlow, emscale_high, emscale_low, EMRhigh,
# EMRlow, jes_high, jes_low, threejet, narrowBDT, narrowISO, narrowBDTbkg,
# narrowISObkg, wideISObkg - see the unfolder constructor comment in src/unfolder.h.
# JERhigh/JERlow/emscale_high/emscale_low/
# EMRhigh/EMRlow are MC-only (Data silently reprocesses as nominal under that
# filename) - skip "data" for those six. jes_high/jes_low are the reverse: Data-only
# (Photon MC silently reprocesses as nominal under that filename) - skip "photon" for
# those two.
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
  # herwig is only used as the generator-modeling systematic (draw_systematics.C's
  # "herwig" source always compares against dataSystag "nominal") - no need to reprocess
  # it under every JER/emscale/etc. systag too, so only run it alongside nominal.
  if [[ $SYSTAG == "nominal" ]]; then
    bash run_unfold.sh Photon5 herwig nominal &
    bash run_unfold.sh Photon10 herwig nominal &
    bash run_unfold.sh Photon20 herwig nominal &
  fi
fi
wait
#if [[ $DOJET == 1 ]]; then
#  bash run_unfold.sh Jet8 pythia $SYSTAG &
#  bash run_unfold.sh Jet12 pythia $SYSTAG &
#  bash run_unfold.sh Jet20 pythia $SYSTAG &
#  bash run_unfold.sh Jet30 pythia $SYSTAG &
#  wait
#  bash run_unfold.sh Jet50 pythia $SYSTAG &
#  bash run_unfold.sh Jet60 pythia $SYSTAG &
#  bash run_unfold.sh Jet80 pythia $SYSTAG &
#fi
wait
echo "All Done!"
