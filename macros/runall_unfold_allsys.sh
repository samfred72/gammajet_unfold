#/bin/bash
# Usage: runall_unfold_allsys.sh [jet|photon|data]
# One-pass, all-systematics version of runall_unfold.sh: for each trigger/sim
# combination, calls run_unfold_allsys.sh once - reads that trigger's tree once and
# fills every systag's histograms in that one pass - instead of runall_unfold.sh's old
# per-systag loop (one tree read per systag, ten reads total for the full systag list).
# Produces exactly the same hists/*_unfolding.root and insitu/*_insitu.root files as
# calling runall_unfold.sh once per systag would (validated bin-for-bin against that
# approach) - see unfold_allsys.C's header comment for the full systag list.
#
# herwig is still nominal-only (only used as the generator-modeling systematic - see
# runall_unfold.sh's comment), so it stays on the plain single-systag run_unfold.sh
# rather than the all-systag path - no benefit to batching a single systag.
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
#  bash run_unfold_allsys.sh Jet8 pythia &
  bash run_unfold_allsys.sh Jet12_long pythia &
#  bash run_unfold_allsys.sh Jet20 pythia &
#  bash run_unfold_allsys.sh Jet30 pythia &
#  wait
#  bash run_unfold_allsys.sh Jet50 pythia &
#  bash run_unfold_allsys.sh Jet60 pythia &
#  bash run_unfold_allsys.sh Jet80 pythia &
fi
wait
echo "All Done!"
