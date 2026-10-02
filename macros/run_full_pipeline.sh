#!/bin/bash
# run_full_pipeline.sh - the complete chain, including the in-situ JES feedback:
#
#   1. compile libgammajet_unfold.so
#   2. unfold_allsys (every trigger/sample, every systag) + puritymaker (every systag),
#      using whatever in-situ constants src/ana.h currently holds
#   2b. rederive the Data/MC v_z + cluster-pT reweight (reweight/make_vz_pt_reweight.C,
#      refilled from the trees) with the current selection and the pass-0 purities - the
#      weight uses the purity-corrected Data spectrum, so it needs a purity pass first;
#      every later unfold pass uses the new weight
#   3. in-situ scan (insitu/run_grid.sh: grid_insitu.C for every systag, then
#      draw_jes_summary.C, which rewrites src/ana.h's jesNominal/jesStatErrLow/High and
#      the per-systag jesBySystag table, then draw_jes_variations.C)
#   4. recompile, so unfolder.cc sees the new constants
#   5. unfold_allsys + puritymaker again, now with the updated in-situ correction
#   6. repeat 3-5 until the p_a table stops moving (see below), then
#   7. in-situ plots used in the note, on the final inputs (see insitu_plots below):
#      grid_insitu (nominal), the unfolded-level mean and shape-chi2 scans, the Jet12
#      reference scan, the reco-level xJ comparisons, draw_jes_summary and
#      draw_jes_variations
#   8. drawing stage (run_pipeline.sh draw: systematics, final result, etc. per radius)
#
# Why a loop: the in-situ scan itself reads the jet pT UNcorrected (insitu_tree's
# rawJetPt), but it also reads the purities, and the purities come from the Data region
# A/B/C/D histograms, whose photon-jet pairing uses the JES-CORRECTED jet pT. So a new
# p_a changes the purities slightly, which changes the next scan slightly. Each
# iteration compares the freshly scanned table with the one the last unfold pass used;
# once every entry (jesNominal and every jesBySystag row) agrees to within --tol, the
# last unfold pass is consistent with the in-situ result and the loop stops. The first
# scan always runs, so the minimum is unfold -> scan -> unfold.
#
# Usage (from anywhere):
#   bash macros/run_full_pipeline.sh [--max-iter N] [--tol X] [--no-draw]
#     --max-iter N   maximum number of in-situ scans (default 3)
#     --tol X        convergence tolerance on |delta p_a| (default 0.001)
#     --no-draw      stop after the unfold/in-situ loop, skip the in-situ plots and drawing stage
#
# Each step writes its own log under logs/full_pipeline_<timestamp>/ (the terminal only
# shows one line per step), the script stops at the first failing step, and the
# per-event progress line of the unfolder is turned off (GAMMAJET_NO_PROGRESS=1).
# Long-running (the unfold passes read every tree): launch it in the background, e.g.
#   nohup bash macros/run_full_pipeline.sh > logs/full_pipeline.out 2>&1 &
#
# Selection defaults (src/ana.h, src/insitu_utility.h, Sep 28 2026): photons |eta^gamma| <
# ana::photonEtaCut = 0.7 (reco and truth), in-situ scan range [0.80, 1.00]. Any scan result
# still at the lower edge (insitu_utility::scanLow) is a bound, not a measurement, and is
# still written into the table - check insitu/pdfs/draw_jes_variations.pdf (edge values
# are drawn open red) before trusting the systematics.

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

# step NAME CMD... : run CMD with its output in $LOGDIR/NAME.log; stop the pipeline if it fails
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

# Largest |difference| between the p_a values (jesNominal + every jesBySystag row) of two
# snapshots; the statistical-error lines are not compared.
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
# The in-situ figures the note uses, on the final (converged) inputs. run_grid.sh already
# wrote grid_insitu_<systag> for every systag plus draw_jes_summary/draw_jes_variations in
# the last loop iteration; they are redrawn here so the plots and ana.h agree even when the
# loop stopped on convergence. The two unfolded-level scans unfold the Data at every scan
# point, so they are given 1000 points (step 2e-4) instead of scanN = 2000 (step 1e-4) to
# keep their runtime where it was before the scan range was widened to 0.80-1.00.
insitu_plots() {
  cd "$TOPDIR/insitu"
  root -b -l -q 'grid_insitu.C("nominal")'
  root -b -l -q 'grid_insitu_unfolded.C("nominal", 1000)'
  root -b -l -q 'grid_insitu_unfolded_shapechi2.C("nominal", 1000)'
  root -b -l -q 'grid_insitu_jet12.C("nominal")'
  root -b -l -q 'draw_insitu_xj.C'
  root -b -l -q 'draw_insitu_xj_allR.C("nominal")'
  root -b -l -q 'draw_insitu_xj_jet12.C'
  root -b -l -q 'draw_jes_summary.C'
  root -b -l -q 'draw_jes_variations.C'
}

# ---- initial pass with the current constants ----
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
