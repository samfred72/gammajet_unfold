#!/bin/bash
# JER iteration: nominal in-situ p_a for one MC JER setting, using main's macros and library unchanged.
# GAMMAJET_JER = none (jet_pt_calib) or a JER template file re-smeared per jet (see ana::mcJerMode).
# Each pass: unfold Data + Photon5/10/20 -> puritymaker -> v_z/cluster-pT reweight -> re-unfold the MC
# -> self-consistent (MC-folded) grid_insitu scan -> new nominal p_a into the work dir's JES table
# (GAMMAJET_JES_TABLE). Stops when the table moves by less than TOL, or after MAXPASS passes.
# The hand-off file (per-radius nominal p_a) is written to WORK/insitu_jes_<tag>.root.
#
# usage: JER_iter/run_iter.sh TAG JER [START_TABLE]
#   TAG          name of the run; work dir JER_iter/work/TAG (made by setup_work.sh if missing)
#   JER          none, or a JER template file (export GAMMAJET_JER_HIST for a non-default TH1 name)
#   START_TABLE  optional JES table file to start from (default: the work dir's, else ana.h's nominal)
# Nominal systag only. Long-running: launch in the background or under Slurm.
set -eo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MAIN=$(dirname "$HERE")
TAG=${1:?tag}
JER=${2:?JER: none or a template file}
START=${3:-}
TOL=0.0005
MAXPASS=4
WORK=$HERE/work/$TAG
[ -d "$WORK" ] || bash "$HERE/setup_work.sh" "$TAG"
TABLE=$WORK/jes_table.txt
if [ -n "$START" ]; then cp "$START" "$TABLE"; fi
if [ ! -f "$TABLE" ]; then   # ana.h's nominal row, written by the library itself
  (cd "$HERE" && root -b -l -q 'write_table.C("'"$TABLE"'")') > /dev/null
fi
export GAMMAJET_UNFOLD=$WORK GAMMAJET_JER=$JER GAMMAJET_JES_TABLE=$TABLE GAMMAJET_NO_PROGRESS=1
export LD_LIBRARY_PATH=$MAIN/src:$LD_LIBRARY_PATH
LOGDIR=$WORK/logs
mkdir -p "$LOGDIR"

unfold() { (cd "$MAIN/macros" && root -b -l -q "unfold.C(\"$1\",\"pythia\",\"nominal\")") > "$LOGDIR/pass${pass}_unfold$2_$1.log" 2>&1; }
unfold_set() {   # unfold the given triggers in parallel; fail if any did
  local pids=() tr
  for tr in "${@:2}"; do unfold "$tr" "$1" & pids+=($!); done
  for p in "${pids[@]}"; do wait "$p"; done
}
nominal_row() { awk '$1=="nominal"{$1=""; print substr($0,2)}' "$TABLE"; }

echo "JER=$JER  start table: $(nominal_row)" | tee "$LOGDIR/table.txt"
for pass in $(seq 1 $MAXPASS); do
  echo "pass $pass $(date)"
  unfold_set A Data Photon5 Photon10 Photon20
  (cd "$MAIN/macros" && root -b -l -q 'puritymaker.C("nominal")') > "$LOGDIR/pass${pass}_purity.log" 2>&1
  (cd "$MAIN/reweight" && root -b -l -q 'make_vz_pt_reweight.C(true)') > "$LOGDIR/pass${pass}_reweight.log" 2>&1
  unfold_set B Photon5 Photon10 Photon20
  (cd "$MAIN/insitu" && root -b -l -q 'grid_insitu.C("nominal")') > "$LOGDIR/pass${pass}_insitu.log" 2>&1
  cp "$WORK/insitu/output/grid_insitu_nominal.root" "$LOGDIR/pass${pass}_grid_insitu_nominal.root"
  (cd "$HERE" && root -b -l -q read_pa.C) | grep '^PA' > "$LOGDIR/pass${pass}_pa.txt"
  cp "$TABLE" "$LOGDIR/pass${pass}_table_in.txt"
  maxd=$(python3 - "$TABLE" "$LOGDIR/pass${pass}_pa.txt" <<'PYEOF'
import sys
table, pafile = sys.argv[1], sys.argv[2]
pa = [float(l.split()[1]) for l in open(pafile)]
assert len(pa) == 7
lines = open(table).read().splitlines()
out, old = [], None
for l in lines:
    f = l.split('#')[0].split()
    if f and f[0] == "nominal":
        old = [float(x) for x in f[1:8]]
        l = "nominal " + " ".join("%.4f" % x for x in pa)
    out.append(l)
open(table, "w").write("\n".join(out) + "\n")
print("%.4f" % max(abs(a - b) for a, b in zip(pa, old)))
PYEOF
)
  echo "pass $pass: $(nominal_row)  max|dp_a| $maxd" | tee -a "$LOGDIR/table.txt"
  if python3 -c "import sys; sys.exit(0 if $maxd < $TOL else 1)"; then echo "converged" | tee -a "$LOGDIR/table.txt"; break; fi
done
(cd "$HERE" && root -b -l -q "write_jes_file.C(\"$WORK/insitu_jes_$TAG.root\", \"$JER ${GAMMAJET_JER_HIST:-}; radii without a TH1 use the nearest radius\")") | tee -a "$LOGDIR/table.txt"
echo "done $(date)"
