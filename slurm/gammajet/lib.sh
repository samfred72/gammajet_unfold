# Shared by submit.sh and the gamma+jet job scripts (slurm/gammajet). The steps themselves are
# macros/pipeline_steps.sh (shared with the local runner, macros/run_full_pipeline.sh). Expects
# GAMMAJET_UNFOLD and, for submissions, GJ_LOGDIR and GJ_DRAW (exported by submit.sh).
: "${GAMMAJET_UNFOLD:?set GAMMAJET_UNFOLD to the gammajet_unfold checkout (see README.md)}"
GJ_SLURM="$GAMMAJET_UNFOLD/slurm/gammajet"
source "$GAMMAJET_UNFOLD/macros/pipeline_steps.sh"

gj_sub() { sbatch --parsable --export=ALL "$@"; }
gj_n() { set -- $1; echo $#; }
gj_dep() { [[ -n $1 ]] && echo "--dependency=afterok:$1"; }

# gj_submit_unfold NAME DEP TASKS: unfold array over TASKS (pipeline_steps.sh PS_UNFOLD_*). Prints the job id.
gj_submit_unfold() {
  local name=$1 dep=$2 tasks=$3 jid
  export GJ_UNFOLD_TASKS="$tasks"
  jid=$(gj_sub $(gj_dep $dep) --array=0-$(( $(gj_n "$tasks") - 1 )) -J gj_unfold_$name \
          -o "$GJ_LOGDIR/unfold_${name}_%a.log" "$GJ_SLURM/unfold.sbatch")
  echo "unfold $name $(printf '%*s' $((12-${#name})) '')$jid ($(gj_n "$tasks") tasks)" >&2
  echo $jid
}

# gj_submit_systag_array NAME SCRIPT DEP: one task per ana::systag. Prints the job id.
gj_submit_systag_array() {
  local name=$1 script=$2 dep=$3 jid
  export GJ_SYSTAGS="$(ps_systags)"
  jid=$(gj_sub $(gj_dep $dep) --array=0-$(( $(gj_n "$GJ_SYSTAGS") - 1 )) -J gj_$name \
          -o "$GJ_LOGDIR/${name}_%a.log" "$GJ_SLURM/$script")
  echo "$name $(printf '%*s' $((19-${#name})) '')$jid ($(gj_n "$GJ_SYSTAGS") systags)" >&2
  echo $jid
}

# gj_submit_final DEP: in-situ figures, R=0.4 diagnostics and the per-radius draw chain (array).
gj_submit_final() {
  local dep=$1 jid
  if [[ ${GJ_DRAW:-1} -eq 0 ]]; then echo "final stage skipped (--no-draw)" >&2; echo $dep; return; fi
  jid=$(gj_sub $(gj_dep $dep) -o "$GJ_LOGDIR/insitu_plots.log" "$GJ_SLURM/insitu_plots.sbatch")
  echo "in-situ plots       $jid" >&2
  jid=$(gj_sub $(gj_dep $dep) -o "$GJ_LOGDIR/diagnostics.log" "$GJ_SLURM/diagnostics.sbatch")
  echo "R=0.4 diagnostics   $jid" >&2
  jid=$(gj_sub $(gj_dep $dep) --array=0-$(( $(gj_n "$PS_RADII") - 1 )) -o "$GJ_LOGDIR/draw_ir%a.log" "$GJ_SLURM/draw.sbatch")
  echo "draw per radius     $jid ($(gj_n "$PS_RADII") radii)" >&2
  echo $jid
}
