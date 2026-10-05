#!/bin/bash
# run_full_pipeline.sh: the complete chain with the in-situ JES feedback.
#   1. compile
#   2. unfold_allsys + puritymaker (current ana.h constants)
#   2b. rederive the v_z / cluster-pT reweight (needs the pass-0 purities)
#   3. in-situ scan (insitu/run_grid.sh), which rewrites ana.h's JES constants
#   4. recompile; 5. unfold_allsys + puritymaker again
#   6. repeat 3-5 until every p_a moves by less than --tol
#   7. in-situ plots for the note; 8. run_pipeline.sh draw
#
# The loop is needed because the scan reads the purities, which depend on the JES-corrected
# pairing.
#
# Usage: bash macros/run_full_pipeline.sh [--max-iter N (3)] [--tol X (0.001)] [--no-draw]
# Logs go to logs/full_pipeline_<timestamp>/. Long-running; launch in the background.
# A scan result at insitu_utility::scanLow is a bound, not a measurement (drawn open red in
# insitu/pdfs/draw_jes_variations.pdf).

set -eo pipefail
MACRODIR="$(cd "$(dirname "$0")" && pwd)"
TOPDIR="$(cd "$MACRODIR/.." && pwd)"
ANA="$TOPDIR/src/ana.h"

MAXITER=3
TOL=0.001
DRAW=1
while [[ $# -gt 0 ]]; do
  case "$1" in
    --max-iter) MAXITER="$2"; shift 2 ;;
    --tol)      TOL="$2"; shift 2 ;;
    --no-draw)  DRAW=0; shift ;;
    *) echo "Usage: $0 [--max-iter N] [--tol X] [--no-draw]" >&2; exit 1 ;;
  esac
done

LOGDIR="$TOPDIR/logs/full_pipeline_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$LOGDIR"
export GAMMAJET_NO_PROGRESS=1
echo "Logs: $LOGDIR   (max-iter $MAXITER, tol $TOL, draw $DRAW)"

# step NAME CMD...: run CMD, logging to $LOGDIR/NAME.log; stop on failure
step() {
  local name="$1"; shift
  local t0=$(date +%s)
  echo "[$(date +%T)] $name ..."
  if ! ( "$@" ) > "$LOGDIR/$name.log" 2>&1; then
    echo "[$(date +%T)] FAILED: $name - see $LOGDIR/$name.log"
    exit 1
  fi
  echo "[$(date +%T)] $name done ($(( ($(date +%s) - t0) / 60 )) min)"
}

# The in-situ constants block of ana.h: jesNominal through "END jesBySystag".
jes_snapshot() { sed -n '/static constexpr double jesNominal\[nJetR\]/,/END jesBySystag/p' "$ANA" > "$1"; }

# Largest |difference| between the p_a values of two snapshots.
jes_maxdiff() {
  python3 - "$1" "$2" <<'PYEOF'
import re, sys
def table(path):
    rows = []
    for line in open(path):
        if "jesStatErr" in line or "{" not in line or "}" not in line:
            continue
        if "jesNominal" in line or line.strip().startswith("{"):
            rows.append([float(x) for x in re.findall(r"-?\d+\.\d+", line.split("//")[0])])
    return rows
a, b = table(sys.argv[1]), table(sys.argv[2])
if len(a) != len(b) or any(len(x) != len(y) for x, y in zip(a, b)):
    print("inf")
else:
    print(max((abs(x - y) for ra, rb in zip(a, b) for x, y in zip(ra, rb)), default=0.0))
PYEOF
}

compile()          { cd "$TOPDIR/src" && bash make.sh; }
unfold_all()       { cd "$MACRODIR" && bash runall_unfold_allsys.sh; }
purity_all()       { cd "$MACRODIR" && bash runall_puritymaker.sh; }
insitu_scan()      { cd "$TOPDIR/insitu" && bash run_grid.sh; }
# In-situ figures for the note, redrawn on the converged inputs. The unfolded-level scans use
# 1000 points to bound their runtime.
insitu_plots() {
  cd "$TOPDIR/insitu"
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

# ---- initial pass ----
step 01_compile compile
step 02_unfold_pass0 unfold_all
step 03_purity_pass0 purity_all
step 04_reweight bash -c "cd '$TOPDIR/reweight' && root -b -l -q 'make_vz_pt_reweight.C(true)'"

# ---- in-situ feedback loop ----
converged=0
for (( iter = 1; iter <= MAXITER; iter++ )); do
  jes_snapshot "$LOGDIR/jes_used_by_pass$((iter-1)).txt"
  step "1${iter}a_insitu_scan" insitu_scan
  if ! grep -q "Updated .*ana.h" "$LOGDIR/1${iter}a_insitu_scan.log"; then
    echo "draw_jes_summary.C did not update ana.h (incomplete scan sweep?) - see $LOGDIR/1${iter}a_insitu_scan.log"
    exit 1
  fi
  jes_snapshot "$LOGDIR/jes_scan$iter.txt"
  diff=$(jes_maxdiff "$LOGDIR/jes_used_by_pass$((iter-1)).txt" "$LOGDIR/jes_scan$iter.txt")
  echo "    scan $iter: max |delta p_a| vs the table used by pass $((iter-1)) = $diff"
  step "1${iter}b_compile" compile
  if [[ $iter -gt 1 ]] && python3 -c "import sys; sys.exit(0 if float('$diff') < float('$TOL') else 1)"; then
    echo "    converged (< $TOL): pass $((iter-1)) outputs are consistent with the in-situ result"
    converged=1
    break
  fi
  step "1${iter}c_unfold_pass$iter" unfold_all
  step "1${iter}d_purity_pass$iter" purity_all
done
if [[ $converged -eq 0 ]]; then
  echo "WARNING: in-situ table still moving after $MAXITER scans (last max |delta p_a| = $diff)."
  echo "         The unfold outputs use the table from the last scan; rerun with a larger --max-iter to converge."
fi

# ---- in-situ plots and drawing ----
if [[ $DRAW -eq 1 ]]; then
  step 19_insitu_plots insitu_plots
  step 20_draw bash -c "cd '$MACRODIR' && bash run_pipeline.sh draw"
fi
echo "[$(date +%T)] Full pipeline done. Logs: $LOGDIR"
