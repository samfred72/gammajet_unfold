# The analysis chain's steps, defined once. Sourced by macros/run_full_pipeline.sh (local runner) and by
# the slurm/gammajet job scripts, which run the same stages in the same order:
#   1 build
#   2 unfold all samples (ps_unfold over PS_UNFOLD_ALL) with the current ana.h table; writes the
#     in-situ trees
#   3 puritymaker per systag            4 v_z / cluster-pT reweight (needs the stage-3 purities)
#   5 re-unfold the MC (PS_UNFOLD_MC) with the new reweight
#   6 in-situ scan per systag: grid_insitu.C recomputes the Data purity from the in-situ trees at each
#     trial p_a until the two agree, so no unfold -> purity -> scan loop is needed
#   7 ps_jes_update: draw_jes_summary.C rewrites src/ana.h, draw_jes_variations.C, recompile
#   8 re-unfold Data (PS_UNFOLD_DATA) with the new p_a       9 puritymaker per systag
#  10 in-situ figures, R=0.4 diagnostics, per-radius prior sensitivity -> systematics -> final result
# Only Data depends on the in-situ p_a (unfolder.cc applies it to Data jets), only the MC on the
# reweight. Each step function runs in its own directory and fails on error.
: "${GAMMAJET_UNFOLD:?set GAMMAJET_UNFOLD to the gammajet_unfold checkout (see README.md)}"
export GAMMAJET_NO_PROGRESS=1
source "$GAMMAJET_UNFOLD/macros/jes_table.sh"

# Unfold tasks, trigger:sim:mode (allsys = unfold_allsys.C over every systag; nominal = unfold.C).
PS_UNFOLD_DATA="Data:pythia:allsys"
PS_UNFOLD_MC="Photon5:pythia:allsys Photon10:pythia:allsys Photon20:pythia:allsys Photon5:herwig:nominal Photon10:herwig:nominal Photon20:herwig:nominal Jet12_long:pythia:allsys"
PS_UNFOLD_ALL="$PS_UNFOLD_DATA $PS_UNFOLD_MC"
PS_RADII="0 1 2 3 4 5 6"
# Stage names, in order, for --from (resume a chain at that stage; earlier outputs are reused).
PS_STAGES="build unfold_all purity reweight unfold_mc insitu jes_update unfold_data purity_final final"

# ps_stage_index NAME: position in PS_STAGES (exits with an error for an unknown name).
ps_stage_index() {
  local i=0 st
  for st in $PS_STAGES; do [[ $st == "$1" ]] && { echo $i; return; }; i=$((i+1)); done
  echo "Unknown stage '$1' (one of: $PS_STAGES)" >&2; return 1
}

# ps_init_logdir DIR NEW: with NEW=1 create DIR and record the starting JES table (jes_start.txt, which
# ps_jes_update compares against); else DIR must be an existing run's log directory (--logdir), whose
# jes_start.txt is kept.
ps_init_logdir() {
  if [[ $2 -eq 1 ]]; then
    mkdir -p "$1"
    jes_snapshot "$GAMMAJET_UNFOLD/src/ana.h" "$1/jes_start.txt"
  elif [[ ! -f $1/jes_start.txt ]]; then
    echo "--logdir $1: no jes_start.txt there (not a pipeline log directory?)" >&2; return 1
  fi
}

# ana::systags, read from src/ana.cc so there is no second list to keep in sync.
ps_systags() {
  python3 - "$GAMMAJET_UNFOLD/src/ana.cc" <<'PYEOF'
import re, sys
src = open(sys.argv[1]).read()
block = re.search(r"ana::systags\s*=\s*\{(.*?)\};", src, re.S).group(1)
print(" ".join(re.findall(r'"([^"]+)"', block)))
PYEOF
}

ps_build()    { cd "$GAMMAJET_UNFOLD/src" && bash make.sh; }
ps_unfold() {   # trigger:sim:mode
  local trigger sim mode
  IFS=: read -r trigger sim mode <<< "$1"
  cd "$GAMMAJET_UNFOLD/macros"
  if [[ $mode == allsys ]]; then root -b -l -q "unfold_allsys.C(\"$trigger\",\"$sim\")"
  else root -b -l -q "unfold.C(\"$trigger\",\"$sim\",\"nominal\")"; fi
}
ps_purity()   { cd "$GAMMAJET_UNFOLD/macros" && root -b -l -q "puritymaker.C(\"$1\")"; }
ps_reweight() { cd "$GAMMAJET_UNFOLD/reweight" && root -b -l -q 'make_vz_pt_reweight.C(true)'; }
ps_insitu()   { cd "$GAMMAJET_UNFOLD/insitu" && root -b -l -q "grid_insitu.C(\"$1\")"; }

# ps_jes_update LOGDIR: combine the scans into src/ana.h, report the change from LOGDIR/jes_start.txt
# (the table stage 2 used), recompile.
ps_jes_update() {
  local logdir=$1 ana="$GAMMAJET_UNFOLD/src/ana.h"
  cd "$GAMMAJET_UNFOLD/insitu"
  root -b -l -q "draw_jes_summary.C" | tee "$logdir/jes_summary.log"
  if ! grep -q "Updated .*ana.h" "$logdir/jes_summary.log"; then
    echo "draw_jes_summary.C did not update ana.h (missing systag scan output?)"; return 1
  fi
  root -b -l -q "draw_jes_variations.C"
  jes_snapshot "$ana" "$logdir/jes_final.txt"
  echo "max |delta p_a| vs the starting table: $(jes_maxdiff "$logdir/jes_start.txt" "$logdir/jes_final.txt")"
  cd "$GAMMAJET_UNFOLD/src" && bash make.sh
}

# In-situ figures for the note on the final inputs. The unfolded-level scans use 1000 points.
ps_insitu_plots() {
  cd "$GAMMAJET_UNFOLD/insitu"
  root -b -l -q 'grid_insitu.C("nominal")'
  root -b -l -q 'grid_insitu_unfolded.C("nominal", 1000)'
  root -b -l -q 'grid_insitu_unfolded.C("nominal", 1000, "shape")'
  root -b -l -q 'grid_insitu_jet12.C("nominal")'
  root -b -l -q 'draw_insitu_xj.C'
  root -b -l -q 'draw_insitu_xj_allR.C("nominal")'
  root -b -l -q 'draw_insitu_xj_jet12.C'
  root -b -l -q 'draw_jes_summary.C'
  root -b -l -q 'draw_jes_variations.C'
}
# R=0.4 nominal unfolding diagnostics. The two toy_iterations runs read each other's output.
ps_diagnostics() {
  cd "$GAMMAJET_UNFOLD/drawing"
  root -b -l -q 'draw_purity_corrected.C("nominal")'
  root -b -l -q 'draw_iteration_halfclosure.C("nominal")'
  root -b -l -q 'toy_iterations.C("nominal","resp")'
  root -b -l -q 'toy_iterations.C("nominal","data")'
  root -b -l -q 'draw_refolding.C("nominal")'
  root -b -l -q 'draw_nonclosure.C()'
}
# Per radius: each reads the previous one's output; outputs are per radius.
ps_draw() {
  cd "$GAMMAJET_UNFOLD/drawing"
  root -b -l -q "draw_prior_sensitivity.C(\"nominal\", $1)"
  root -b -l -q "draw_systematics.C($1)"
  root -b -l -q "draw_final_result.C($1)"
}
