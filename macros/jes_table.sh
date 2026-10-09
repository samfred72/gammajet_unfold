# In-situ JES table helpers, sourced by run_full_pipeline.sh and slurm/gammajet/controller.sbatch.

# jes_snapshot ANA_H OUT: the in-situ constants block of ana.h, jesNominal through "END jesBySystag".
jes_snapshot() { sed -n '/static constexpr double jesNominal\[nJetR\]/,/END jesBySystag/p' "$1" > "$2"; }

# jes_maxdiff A B: largest |difference| between the p_a values of two snapshots ("inf" if the
# tables differ in shape).
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
